#!/usr/bin/env python3
"""Validate and exercise the local jet-radius physics study on Delphes ROOT.

The command is intentionally local-only.  It reads one ROOT file and writes a
JSON report plus diagnostic plots; it has no batch or remote-I/O path.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import sys
import traceback
from typing import Any, Sequence

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import FancyBboxPatch
import numpy as np


HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parent
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))


RADIUS_COLORS = {
    0.2: "#6B7280",  # gray
    0.4: "#7B2CBF",  # purple
    0.6: "#2563EB",  # blue
    0.8: "#111111",  # black: analysis default
    1.0: "#E76F51",  # orange
    1.2: "#D81B60",  # magenta
    1.4: "#0086A8",  # cyan-blue
    1.6: "#C62828",  # red
}

# Deliberately varied, non-green palette for individual dark-hadron ancestry
# groups. Descendants inherit the color of their initial dark-hadron parent.
DARK_HADRON_COLORS = (
    "#2563EB", "#D81B60", "#E76F51", "#7B2CBF", "#0086A8", "#C62828",
    "#F59E0B", "#6B7280", "#EC4899", "#4F46E5", "#92400E", "#0E7490",
)


def radius_color(radius: float) -> str:
    for configured, color in RADIUS_COLORS.items():
        if math.isclose(float(radius), configured, abs_tol=1.0e-9):
            return color
    fallback = ("#5B21B6", "#1D4ED8", "#B45309", "#BE123C", "#0E7490")
    return fallback[int(round(float(radius) * 10.0)) % len(fallback)]

from coffea.nanoevents import NanoEventsFactory  # noqa: E402
from common import DelphesSchema2  # noqa: E402
from core import (  # noqa: E402
    COLLECTIONS,
    DEFAULT_R_GRID,
    ClusteredJet,
    bootstrap_mean,
    build_fixed_truth_event,
    card_roundtrip_check,
    delta_r2,
    dark_hadron_group_metrics,
    fixed_axis_containment,
    inside_indices,
    json_ready,
    load_raw_ancestry,
    match_by_shared_visible_pt,
    matched_jet_metrics,
    particle_table,
    radial_percentiles,
    recluster_event,
    shape_observables,
    softdrop_jet,
    stored_jets,
    study_metadata,
    visible_storage_sets,
    wrap_delta_phi,
    write_json,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Local events.root")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--max-events", type=int, default=200)
    parser.add_argument(
        "--card", type=Path, default=REPOSITORY / "cards" / "delphes_card_CMS.tcl"
    )
    parser.add_argument(
        "--delphes", type=Path, default=REPOSITORY / "install" / "delphes"
    )
    parser.add_argument(
        "--collections",
        nargs="+",
        default=list(COLLECTIONS),
        choices=sorted(COLLECTIONS),
    )
    parser.add_argument(
        "--radii", nargs="+", type=float, default=list(DEFAULT_R_GRID)
    )
    parser.add_argument(
        "--diagnostic-radii",
        nargs="+",
        type=float,
        default=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6],
        help="Readable subset used in overlaid one-dimensional histograms",
    )
    parser.add_argument(
        "--event-displays",
        type=int,
        default=8,
        help="Number of different events to render with all diagnostic radii overlaid",
    )
    parser.add_argument("--bootstrap-resamples", type=int, default=500)
    args = parser.parse_args()
    if args.event_displays < 1:
        parser.error("--event-displays must be at least 1")
    return args


def load_events(path: Path, max_events: int) -> Any:
    if not path.is_file():
        raise FileNotFoundError(path)
    return NanoEventsFactory.from_root(
        {str(path.resolve()): "Delphes"},
        schemaclass=DelphesSchema2,
        entry_stop=max_events,
    ).events()


def check_card(card: Path, collections: Sequence[str]) -> dict[str, Any]:
    text = card.read_text()
    modules = sorted({COLLECTIONS[name].module for name in collections})
    results = []
    for module in modules:
        baseline = COLLECTIONS[
            next(name for name in collections if COLLECTIONS[name].module == module)
        ].default_r
        result = card_roundtrip_check(text, module, baseline + 0.137)
        result.pop("text")
        results.append(result)
    return {
        "passed": all(item["passed"] for item in results),
        "card": str(card.resolve()),
        "modules": results,
    }


def _dr(first: ClusteredJet, second: ClusteredJet) -> float:
    return float(math.sqrt(delta_r2(first.eta, first.phi, second.eta, second.phi)))


FLOAT32_EPS = float(np.finfo(np.float32).eps)


def float32_mass_noise(observed: float, stored: float, energy: float) -> bool:
    """True if a mass difference is within Delphes' Float_t storage noise.

    Delphes stores four-vector components as Float_t, so m^2 = E^2 - p^2 carries
    an absolute error of about 2*eps*E^2. For very forward jets (E >> pT) that
    is ~0.1-1 GeV, larger than the fixed mass tolerances, even though
    membership and pT close exactly. Compared on signed m^2 because FastJet
    returns negative mass for slightly negative m^2. The factor 4 is a 2x
    margin on the 2*eps*E^2 estimate.
    """
    signed_m2 = lambda m: m * abs(m)
    return abs(signed_m2(observed) - signed_m2(stored)) <= 4.0 * FLOAT32_EPS * energy**2


def check_collection_closure(
    events: Any,
    collection: str,
    n_events: int,
) -> dict[str, Any]:
    spec = COLLECTIONS[collection]
    required = {spec.jet_branch, spec.candidate_branch}
    missing = sorted(required - set(ak.fields(events)))
    if missing:
        return {
            "passed": False,
            "collection": collection,
            "missing_branches": missing,
        }

    count_mismatch = 0
    membership_mismatch = 0
    unmatched_stored_jets = 0
    extra_reclustered_jets = 0
    mass_tolerance_violations = 0
    float32_mass_noise_jets = 0
    compared_jets = 0
    max_pt_relative = 0.0
    max_mass_relative = 0.0
    max_mass_absolute = 0.0
    max_delta_r = 0.0
    examples: list[dict[str, int]] = []

    for event_index in range(n_events):
        event = events[event_index]
        candidates = particle_table(getattr(event, spec.candidate_branch))
        reclustered = recluster_event(candidates, spec.default_r, spec.pt_min)
        stored = stored_jets(getattr(event, spec.jet_branch), candidates)
        if len(reclustered) != len(stored):
            count_mismatch += 1
            if len(examples) < 5:
                examples.append(
                    {
                        "event": event_index,
                        "stored_n": len(stored),
                        "reclustered_n": len(reclustered),
                    }
                )
        reclustered_by_members = {
            jet.constituent_indices: jet for jet in reclustered
        }
        stored_member_sets = {jet.constituent_indices for jet in stored}
        extra_reclustered_jets += sum(
            jet.constituent_indices not in stored_member_sets for jet in reclustered
        )
        matched_pairs: list[tuple[ClusteredJet, ClusteredJet]] = []
        for stored_jet in stored:
            reclustered_jet = reclustered_by_members.get(
                stored_jet.constituent_indices
            )
            if reclustered_jet is None:
                unmatched_stored_jets += 1
                membership_mismatch += 1
                continue
            matched_pairs.append((reclustered_jet, stored_jet))
        for reclustered_jet, stored_jet in matched_pairs:
            compared_jets += 1
            pt_scale = max(abs(stored_jet.pt), 1.0e-12)
            mass_scale = max(abs(stored_jet.mass), 1.0e-12)
            max_pt_relative = max(
                max_pt_relative,
                abs(reclustered_jet.pt - stored_jet.pt) / pt_scale,
            )
            mass_absolute = abs(reclustered_jet.mass - stored_jet.mass)
            max_mass_absolute = max(max_mass_absolute, mass_absolute)
            if abs(stored_jet.mass) > 1.0e-6:
                max_mass_relative = max(
                    max_mass_relative, mass_absolute / mass_scale
                )
            if mass_absolute > max(0.3, 0.005 * abs(stored_jet.mass)):
                if float32_mass_noise(reclustered_jet.mass, stored_jet.mass,
                                      stored_jet.energy):
                    float32_mass_noise_jets += 1
                else:
                    mass_tolerance_violations += 1
            max_delta_r = max(max_delta_r, _dr(reclustered_jet, stored_jet))

    tolerances = {
        "pt_relative": 2.0e-5,
        "mass_relative": 5.0e-3,
        "mass_absolute_floor_GeV": 0.3,
        "float32_mass_squared_bound": "4 * eps_float32 * E^2",
        "delta_r": 2.0e-5,
    }
    raw_four_vector_passed = (
        max_pt_relative <= tolerances["pt_relative"]
        and mass_tolerance_violations == 0
        and max_delta_r <= tolerances["delta_r"]
    )
    if spec.four_vector_policy == "post_energy_scale":
        # The stored small-R Jet branch is UniqueObjectFinder/JetEnergyScale,
        # whereas its membership comes from FastJetFinder.  JEC changes pT and
        # mass by construction; closure here is exact membership plus direction.
        four_vector_passed = max_delta_r <= tolerances["delta_r"]
    else:
        four_vector_passed = raw_four_vector_passed
    count_policy_passed = (
        count_mismatch == 0
        if spec.four_vector_policy == "raw_cluster"
        else unmatched_stored_jets == 0
    )
    passed = count_policy_passed and membership_mismatch == 0 and four_vector_passed
    return {
        "passed": passed,
        "collection": collection,
        "candidate_branch": spec.candidate_branch,
        "radius": spec.default_r,
        "pt_min": spec.pt_min,
        "four_vector_policy": spec.four_vector_policy,
        "raw_four_vector_passed": raw_four_vector_passed,
        "events": n_events,
        "compared_jets": compared_jets,
        "count_mismatch_events": count_mismatch,
        "count_policy_passed": count_policy_passed,
        "membership_mismatch_jets": membership_mismatch,
        "unmatched_stored_jets": unmatched_stored_jets,
        "extra_reclustered_jets": extra_reclustered_jets,
        "mass_tolerance_violations": mass_tolerance_violations,
        "float32_mass_noise_jets": float32_mass_noise_jets,
        "max_pt_relative": max_pt_relative,
        "max_mass_relative": max_mass_relative,
        "max_mass_absolute": max_mass_absolute,
        "max_delta_r": max_delta_r,
        "tolerances": tolerances,
        "mass_precision_note": (
            "ROOT stores Delphes four-vector components as Float_t; the absolute "
            "floor covers cancellation in low-mass, high-momentum jets"
        ),
        "examples": examples,
    }


def check_truth_ancestry(
    events: Any, raw_ancestry: Sequence[Any], n_events: int
) -> dict[str, Any]:
    required = {
        "GenParticle",
        "GenCandidate",
        "DarkHadronCandidate",
        "DarkHadronJet",
        "DarkHadronVisibleJet",
    }
    missing = sorted(required - set(ak.fields(events)))
    if missing:
        return {"passed": False, "missing_branches": missing}

    compared = 0
    mismatch = 0
    jet_count_mismatch = 0
    ordering_permutation_events = 0
    ambiguous_visible = 0
    cache_reassigned_visible = 0
    examples: list[dict[str, Any]] = []
    for event_index in range(n_events):
        event = events[event_index]
        candidates = particle_table(event.GenCandidate)
        truth = build_fixed_truth_event(event, raw_ancestry[event_index])
        stored_sets = visible_storage_sets(event, candidates)
        if len(stored_sets) != len(truth.jets):
            jet_count_mismatch += 1
        ambiguous_visible += len(truth.unresolved_visible)
        cache_reassigned_visible += len(truth.cache_reassigned_visible)
        unused_stored = set(range(len(stored_sets)))
        event_mapping: list[int] = []
        for truth_jet in truth.jets:
            compared += 1
            rebuilt = set(truth_jet.visible_indices)
            exact_matches = [
                stored_index
                for stored_index in sorted(unused_stored)
                if rebuilt == stored_sets[stored_index]
            ]
            if exact_matches:
                stored_index = exact_matches[0]
                unused_stored.remove(stored_index)
                event_mapping.append(stored_index)
                continue
            else:
                mismatch += 1
                if len(examples) < 5:
                    best_index = max(
                        unused_stored,
                        key=lambda index: len(rebuilt & stored_sets[index]),
                        default=-1,
                    )
                    stored_set = (
                        stored_sets[best_index] if best_index >= 0 else set()
                    )
                    examples.append(
                        {
                            "event": event_index,
                            "truth_id": truth_jet.truth_id,
                            "rebuilt_only": sorted(rebuilt - stored_set)[:10],
                            "stored_only": sorted(stored_set - rebuilt)[:10],
                        }
                    )
        if len(event_mapping) == len(truth.jets) and event_mapping != list(
            range(len(truth.jets))
        ):
            ordering_permutation_events += 1
    return {
        "passed": jet_count_mismatch == 0 and mismatch == 0 and ambiguous_visible == 0,
        "events": n_events,
        "compared_truth_jets": compared,
        "jet_count_mismatch_events": jet_count_mismatch,
        "root_order_permutation_events": ordering_permutation_events,
        "ordering_note": (
            "TreeWriter pT-sorts DarkHadronVisibleJet independently; closure is "
            "therefore exact one-to-one constituent-set matching, not array index"
        ),
        "constituent_set_mismatch_jets": mismatch,
        "ambiguous_visible_candidates": ambiguous_visible,
        "cache_reassigned_visible_candidates": cache_reassigned_visible,
        "cache_reassignment_is_physics_warning": cache_reassigned_visible > 0,
        "examples": examples,
    }


def check_softdrop_closure(events: Any, n_events: int) -> dict[str, Any]:
    """Validate the local Soft Drop implementation against stored R=0.8 jets."""

    if not {"GenCandidate", "GenFatJet"}.issubset(set(ak.fields(events))):
        return {"passed": False, "reason": "missing GenCandidate or GenFatJet"}
    compared = 0
    count_mismatch_events = 0
    tolerance_violations = 0
    float32_mass_noise_jets = 0
    pt_tolerance_violations = 0
    max_absolute = 0.0
    max_relative = 0.0
    max_pt_relative = 0.0
    for event_index in range(n_events):
        event = events[event_index]
        candidates = particle_table(event.GenCandidate)
        jets = recluster_event(candidates, 0.8, COLLECTIONS["GenFatJet"].pt_min)
        stored_mass = np.asarray(
            ak.to_numpy(event.GenFatJet.SoftDroppedJet.mass), dtype=float
        )
        stored_pt = np.asarray(
            ak.to_numpy(event.GenFatJet.SoftDroppedJet.pt), dtype=float
        )
        if len(jets) != len(stored_mass):
            count_mismatch_events += 1
        for jet, expected, expected_pt in zip(jets, stored_mass, stored_pt):
            result = softdrop_jet(
                candidates, jet, beta=0.0, zcut=0.1, r0=0.8
            )
            observed = result.mass
            absolute = abs(observed - float(expected))
            relative = absolute / max(abs(float(expected)), 1.0e-12)
            pt_relative = abs(result.pt - float(expected_pt)) / max(
                abs(float(expected_pt)), 1.0e-12
            )
            max_absolute = max(max_absolute, absolute)
            max_relative = max(max_relative, relative)
            max_pt_relative = max(max_pt_relative, pt_relative)
            if absolute > max(0.03, 5.0e-4 * abs(float(expected))):
                energy = math.hypot(result.pt * math.cosh(result.eta), result.mass)
                if float32_mass_noise(observed, float(expected), energy):
                    float32_mass_noise_jets += 1
                else:
                    tolerance_violations += 1
            if pt_relative > 2.0e-5:
                pt_tolerance_violations += 1
            compared += 1
    return {
        "passed": (
            count_mismatch_events == 0
            and tolerance_violations == 0
            and pt_tolerance_violations == 0
        ),
        "definition": "C/A reclustering; Soft Drop beta=0, zcut=0.1, R0=0.8",
        "compared_jets": compared,
        "count_mismatch_events": count_mismatch_events,
        "tolerance_violations": tolerance_violations,
        "float32_mass_noise_jets": float32_mass_noise_jets,
        "pt_tolerance_violations": pt_tolerance_violations,
        "max_mass_absolute_GeV": max_absolute,
        "max_mass_relative": max_relative,
        "max_pt_relative": max_pt_relative,
        "mass_tolerance": (
            "max(0.03 GeV, 5e-4 * stored mass), or |m^2 diff| <= "
            "4 * eps_float32 * E^2 (counted in float32_mass_noise_jets)"
        ),
        "pt_relative_tolerance": 2.0e-5,
        "mass_precision_note": (
            "Soft Drop pT closes at Float_t precision; the absolute mass floor "
            "covers E^2-p^2 cancellation for low-mass, high-pT jets"
        ),
    }


def _visible_mass(candidates: Any, indices: Sequence[int]) -> float:
    if not indices:
        return math.nan
    idx = np.asarray(indices, dtype=np.int64)
    energy = float(np.sum(candidates.energy[idx]))
    px = float(np.sum(candidates.px[idx]))
    py = float(np.sum(candidates.py[idx]))
    pz = float(np.sum(candidates.pz[idx]))
    return math.sqrt(max(energy * energy - px * px - py * py - pz * pz, 0.0))


def _split_merge_counts(truth: Any, jets: Sequence[ClusteredJet]) -> tuple[int, int]:
    truth_sets = [set(jet.visible_indices) for jet in truth.jets]
    clustered_sets = [set(jet.constituent_indices) for jet in jets]
    split = sum(
        sum(bool(truth_set & clustered_set) for clustered_set in clustered_sets) > 1
        for truth_set in truth_sets
    )
    merge = sum(
        sum(bool(truth_set & clustered_set) for truth_set in truth_sets) > 1
        for clustered_set in clustered_sets
    )
    return split, merge


def run_primary_scan(
    events: Any,
    raw_ancestry: Sequence[Any],
    n_events: int,
    radii: Sequence[float],
    bootstrap_resamples: int,
) -> tuple[dict[str, Any], dict[float, dict[tuple[int, int], ClusteredJet]]]:
    """Scan exact GenCandidate jets against the fixed truth partition."""

    fields = set(ak.fields(events))
    required = {"GenCandidate", "DarkHadronCandidate", "DarkHadronJet", "GenParticle"}
    missing = sorted(required - fields)
    if missing:
        return {"passed": False, "missing_branches": missing}, {}

    values: dict[float, dict[str, list[float]]] = {
        radius: defaultdict(list) for radius in radii
    }
    matched_by_radius: dict[float, dict[tuple[int, int], ClusteredJet]] = {
        radius: {} for radius in radii
    }
    split_counts = defaultdict(int)
    merge_counts = defaultdict(int)
    truth_counts = defaultdict(int)
    clustered_counts = defaultdict(int)
    truth_keys: set[tuple[int, int]] = set()
    monotonic_failures = 0
    same_id_failures = 0
    cache_reassigned_visible = 0

    for event_index in range(n_events):
        event = events[event_index]
        candidates = particle_table(event.GenCandidate)
        truth = build_fixed_truth_event(event, raw_ancestry[event_index])
        cache_reassigned_visible += len(truth.cache_reassigned_visible)
        ids_before = tuple(jet.truth_id for jet in truth.jets)
        for truth_jet in truth.jets:
            key = (event_index, truth_jet.truth_id)
            truth_keys.add(key)
            containment = fixed_axis_containment(
                candidates,
                truth_jet.visible_indices,
                truth_jet.eta,
                truth_jet.phi,
                radii,
            )
            finite = containment[np.isfinite(containment)]
            if len(finite) > 1 and np.any(np.diff(finite) < -1.0e-14):
                monotonic_failures += 1

        for radius in radii:
            jets = recluster_event(candidates, radius, COLLECTIONS["GenFatJet"].pt_min)
            matches = match_by_shared_visible_pt(truth, jets, candidates.pt)
            split, merge = _split_merge_counts(truth, jets)
            split_counts[radius] += split
            merge_counts[radius] += merge
            truth_counts[radius] += len(truth.jets)
            clustered_counts[radius] += len(jets)
            if tuple(jet.truth_id for jet in truth.jets) != ids_before:
                same_id_failures += 1
            for match in matches:
                truth_jet = truth.jets[match.truth_id]
                jet = jets[match.clustered_index]
                key = (event_index, match.truth_id)
                matched_by_radius[radius][key] = jet
                metric = matched_jet_metrics(
                    truth, truth_jet, jet, candidates
                )
                for name in (
                    "containment_pt",
                    "containment_count",
                    "contamination_pt",
                    "leading_dark_hadron_containment",
                    "n_dark_hadrons",
                    "n_dark_hadrons_with_visible_descendants",
                    "n_fully_contained_dark_hadrons",
                    "n_not_fully_contained_dark_hadrons",
                    "fraction_fully_contained_dark_hadrons",
                    "all_dark_hadrons_fully_contained",
                    "n_contaminating_dark_hadrons",
                    "n_fully_contained_contaminating_dark_hadrons",
                    "cross_truth_dark_hadron_contamination_pt",
                    "n_non_dark_hadron_constituents",
                    "fraction_non_dark_hadron_constituents",
                    "non_dark_hadron_constituent_pt_fraction",
                    "has_non_dark_hadron_constituents",
                    "jet_pt",
                    "jet_mass",
                ):
                    values[radius][name].append(float(metric[name]))
                truth_mass = _visible_mass(candidates, truth_jet.visible_indices)
                response = jet.mass / truth_mass if truth_mass > 0.0 else math.nan
                values[radius]["visible_mass_response"].append(response)
                values[radius]["association_shared_fraction"].append(
                    match.shared_fraction
                )
                values[radius]["axis_drift_from_truth"].append(
                    math.sqrt(delta_r2(jet.eta, jet.phi, truth_jet.eta, truth_jet.phi))
                )
                for name, value in radial_percentiles(
                    candidates,
                    truth_jet.visible_indices,
                    truth_jet.eta,
                    truth_jet.phi,
                ).items():
                    values[radius][name].append(value)
                for name, value in shape_observables(
                    candidates, jet.constituent_indices, jet.eta, jet.phi
                ).items():
                    values[radius][name].append(value)

    reference_radius = min(radii, key=lambda radius: abs(radius - 0.8))
    for radius in radii:
        common_keys = set(matched_by_radius[radius]) & set(
            matched_by_radius[reference_radius]
        )
        for key in common_keys:
            jet = matched_by_radius[radius][key]
            reference = matched_by_radius[reference_radius][key]
            values[radius]["axis_drift_from_r08"].append(
                math.sqrt(delta_r2(jet.eta, jet.phi, reference.eta, reference.phi))
            )

    summary: dict[str, Any] = {}
    for radius in radii:
        radius_summary = {
            name: bootstrap_mean(
                entries,
                n_resamples=bootstrap_resamples,
                seed=12345 + int(round(radius * 1000)),
            )
            for name, entries in sorted(values[radius].items())
        }
        radius_summary.update(
            {
                "matched_truth_jets": len(matched_by_radius[radius]),
                "truth_jets": truth_counts[radius],
                "clustered_jets": clustered_counts[radius],
                "match_efficiency": len(matched_by_radius[radius])
                / truth_counts[radius]
                if truth_counts[radius]
                else math.nan,
                "split_truth_jets": split_counts[radius],
                "merge_clustered_jets": merge_counts[radius],
            }
        )
        summary[f"{radius:.2f}"] = radius_summary

    return (
        {
            "passed": monotonic_failures == 0 and same_id_failures == 0,
            "policy": (
                "fixed shipped R=0.8 DarkHadronJet truth; anti-kT GenCandidate "
                "scan; one-to-one maximum shared visible-descendant pT"
            ),
            "events": n_events,
            "truth_keys": len(truth_keys),
            "fixed_axis_monotonicity_failures": monotonic_failures,
            "truth_identifier_drift_failures": same_id_failures,
            "cache_reassigned_visible_candidates": cache_reassigned_visible,
            "per_dark_hadron_metrics_exclude_cache_reassigned_candidates": True,
            "reference_radius_for_axis_drift": reference_radius,
            "radii": summary,
        },
        matched_by_radius,
    )


def _save_figure(figure: Any, path: Path) -> str:
    figure.tight_layout()
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)
    return str(path.resolve())


def make_workflow_plot(output: Path) -> str:
    """Draw the collection definitions and the fixed-truth association flow."""

    figure, axis = plt.subplots(figsize=(15, 7.5))
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, 1.0)
    axis.axis("off")

    def box(
        xy: tuple[float, float],
        size: tuple[float, float],
        title: str,
        body: str,
        color: str,
    ) -> None:
        x, y = xy
        width, height = size
        patch = FancyBboxPatch(
            (x, y), width, height,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            facecolor=color, edgecolor="0.25", linewidth=1.4,
        )
        axis.add_patch(patch)
        axis.text(
            x + width / 2, y + height * 0.68, title,
            ha="center", va="center", fontsize=13, fontweight="bold",
        )
        axis.text(
            x + width / 2, y + height * 0.30, body,
            ha="center", va="center", fontsize=9.2, linespacing=1.25,
        )

    def arrow(start: tuple[float, float], end: tuple[float, float], label: str = "") -> None:
        axis.annotate(
            "", xy=end, xytext=start,
            arrowprops={"arrowstyle": "-|>", "lw": 1.8, "color": "0.28"},
        )
        if label:
            axis.text(
                (start[0] + end[0]) / 2,
                (start[1] + end[1]) / 2 + 0.025,
                label, ha="center", va="bottom", fontsize=9.5,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
            )

    box(
        (0.02, 0.58), (0.19, 0.22), "Initial dark hadrons",
        "Immediately after hidden-valley\nfragmentation; before decays",
        "#d9d2e9",
    )
    box(
        (0.28, 0.58), (0.19, 0.22), "DarkHadronJets",
        "anti-$k_T$ clustering of initial dark hadrons\n"
        "Fixed truth partition: R=0.8",
        "#c9daf8",
    )
    box(
        (0.54, 0.72), (0.22, 0.20), "DarkHadronStableJets",
        "Manual ancestry combination of\nall stable descendants\n"
        "(visible + invisible)",
        "#e4d7f5",
    )
    box(
        (0.54, 0.43), (0.22, 0.20), "DarkHadronVisibleJets",
        "Manual ancestry combination of\nvisible stable descendants only",
        "#fff2cc",
    )
    box(
        (0.02, 0.12), (0.25, 0.20), "Stable visible SM particles",
        "GenCandidate pool after\nneutrino/invisible filtering",
        "#fce5cd",
    )
    box(
        (0.34, 0.12), (0.22, 0.20), r"GenFatJets(R)",
        "anti-$k_T$ reclustering at R=0.2...1.6\n"
        "R is the quantity varied",
        "#f4cccc",
    )
    box(
        (0.80, 0.36), (0.19, 0.28), "Compare at fixed\ntruth ID",
        "One-to-one match maximizing\nshared visible-descendant pT\n"
        "Then measure containment,\ncontamination and mass",
        "#ead1dc",
    )

    arrow((0.21, 0.69), (0.28, 0.69), r"anti-$k_T$ R=0.8")
    arrow((0.47, 0.73), (0.54, 0.81), "all descendants")
    arrow((0.47, 0.64), (0.54, 0.53), "visible only")
    arrow((0.27, 0.22), (0.34, 0.22), r"anti-$k_T$ R scan")
    arrow((0.76, 0.53), (0.80, 0.50))
    arrow((0.56, 0.22), (0.80, 0.43))
    axis.text(
        0.5, 0.98, "Jet-radius study workflow",
        ha="center", va="top", fontsize=19, fontweight="bold",
    )
    axis.text(
        0.5, 0.02,
        "Primary scan keeps the DarkHadronJet R=0.8 partition fixed. "
        "Only GenFatJet R changes, so the same physical truth objects are compared.",
        ha="center", va="bottom", fontsize=11.5,
    )
    return _save_figure(figure, output)


def _probability_histogram(
    axis: Any,
    values_by_radius: dict[float, list[float]],
    bins: np.ndarray,
    xlabel: str,
    title: str,
) -> None:
    for radius, entries in sorted(values_by_radius.items()):
        color = radius_color(radius)
        values = np.asarray(entries, dtype=float)
        values = values[np.isfinite(values)]
        if len(values) == 0:
            continue
        weights = np.full(len(values), 1.0 / len(values))
        axis.hist(
            values, bins=bins, weights=weights, histtype="step",
            linewidth=2.0, color=color,
            label=f"R={radius:.1f}  (median={np.median(values):.2f})",
        )
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Fraction of matched truth jets / bin")
    axis.set_yscale("log", nonpositive="clip")
    axis.set_ylim(bottom=1.0e-5)
    axis.set_title(title)
    axis.grid(alpha=0.22, which="both")
    axis.legend(fontsize=8, frameon=False, loc="upper right", ncol=2)


def make_diagnostic_plots(
    events: Any,
    raw_ancestry: Sequence[Any],
    n_events: int,
    diagnostic_radii: Sequence[float],
    matched: dict[float, dict[tuple[int, int], ClusteredJet]],
    scan: dict[str, Any],
    outdir: Path,
) -> dict[str, Any]:
    """Produce explanatory plots using a common fixed-truth jet cohort."""

    radii = sorted(set(float(radius) for radius in diagnostic_radii))
    missing = [radius for radius in radii if radius not in matched]
    if missing:
        return {"passed": False, "reason": f"unscanned diagnostic radii: {missing}"}
    common_keys = set.intersection(*(set(matched[radius]) for radius in radii))
    if not common_keys:
        return {"passed": False, "reason": "no common matched truth-jet cohort"}

    softdrop_mass: dict[float, list[float]] = {radius: [] for radius in radii}
    count_containment: dict[float, list[float]] = {radius: [] for radius in radii}
    pt_containment: dict[float, list[float]] = {radius: [] for radius in radii}
    full_dark_hadron_count: dict[float, list[float]] = {
        radius: [] for radius in radii
    }
    contaminating_dark_hadron_count: dict[float, list[float]] = {
        radius: [] for radius in radii
    }
    fully_contained_contaminating_count: dict[float, list[float]] = {
        radius: [] for radius in radii
    }
    not_fully_contained_dark_hadron_count: dict[float, list[float]] = {
        radius: [] for radius in radii
    }
    non_dark_hadron_constituent_count: dict[float, list[float]] = {
        radius: [] for radius in radii
    }
    non_dark_hadron_constituent_fraction: dict[float, list[float]] = {
        radius: [] for radius in radii
    }
    non_dark_hadron_constituent_pt_fraction: dict[float, list[float]] = {
        radius: [] for radius in radii
    }
    collection_mass: dict[str, list[float]] = {
        "DarkHadronJets": [],
        "DarkHadronStableJets": [],
        "DarkHadronVisibleJets": [],
        "matched GenFatJets (R=0.8)": [],
    }
    truth_by_key: dict[tuple[int, int], Any] = {}
    truth_event_by_event: dict[int, Any] = {}
    candidates_by_event: dict[int, Any] = {}

    for event_index in range(n_events):
        event = events[event_index]
        candidates = particle_table(event.GenCandidate)
        candidates_by_event[event_index] = candidates
        truth_event = build_fixed_truth_event(event, raw_ancestry[event_index])
        truth_event_by_event[event_index] = truth_event
        for truth_jet in truth_event.jets:
            truth_by_key[(event_index, truth_jet.truth_id)] = truth_jet
        for branch, label in (
            ("DarkHadronJet", "DarkHadronJets"),
            ("DarkHadronStableJet", "DarkHadronStableJets"),
            ("DarkHadronVisibleJet", "DarkHadronVisibleJets"),
        ):
            jets = getattr(event, branch)
            masses = np.asarray(ak.to_numpy(jets.Mass), dtype=float)
            pts = np.asarray(ak.to_numpy(jets.PT), dtype=float)
            collection_mass[label].extend(masses[pts > 0.0].tolist())

    reference_radius = min(radii, key=lambda radius: abs(radius - 0.8))
    for key in sorted(common_keys):
        event_index, _ = key
        candidates = candidates_by_event[event_index]
        truth_jet = truth_by_key[key]
        truth_set = set(truth_jet.visible_indices)
        truth_pt = float(np.sum(candidates.pt[list(truth_set)])) if truth_set else 0.0
        collection_mass["matched GenFatJets (R=0.8)"].append(
            matched[reference_radius][key].mass
        )
        for radius in radii:
            jet = matched[radius][key]
            dh_metrics = dark_hadron_group_metrics(
                truth_event_by_event[event_index], truth_jet.truth_id, jet, candidates
            )
            jet_set = set(jet.constituent_indices)
            captured = truth_set & jet_set
            count_containment[radius].append(
                len(captured) / len(truth_set) if truth_set else math.nan
            )
            captured_pt = (
                float(np.sum(candidates.pt[list(captured)])) if captured else 0.0
            )
            pt_containment[radius].append(
                captured_pt / truth_pt if truth_pt > 0.0 else math.nan
            )
            softdrop_mass[radius].append(
                softdrop_jet(
                    candidates, jet, beta=0.0, zcut=0.1, r0=radius
                ).mass
            )
            full_dark_hadron_count[radius].append(
                dh_metrics["n_fully_contained_target_dark_hadrons"]
            )
            contaminating_dark_hadron_count[radius].append(
                dh_metrics["n_contaminating_dark_hadrons"]
            )
            fully_contained_contaminating_count[radius].append(
                dh_metrics["n_fully_contained_contaminating_dark_hadrons"]
            )
            not_fully_contained_dark_hadron_count[radius].append(
                dh_metrics["n_not_fully_contained_target_dark_hadrons"]
            )
            non_dark_hadron_constituent_count[radius].append(
                dh_metrics["n_non_dark_hadron_constituents"]
            )
            non_dark_hadron_constituent_fraction[radius].append(
                dh_metrics["fraction_non_dark_hadron_constituents"]
            )
            non_dark_hadron_constituent_pt_fraction[radius].append(
                dh_metrics["non_dark_hadron_constituent_pt_fraction"]
            )

    outputs: dict[str, str] = {}
    outputs["workflow"] = make_workflow_plot(outdir / "workflow_overview.png")

    all_sd = np.concatenate(
        [np.asarray(values)[np.isfinite(values)] for values in softdrop_mass.values()]
    )
    sd_upper = max(100.0, math.ceil(float(np.percentile(all_sd, 99.5)) / 25.0) * 25.0)
    clipped_sd = {
        radius: np.minimum(values, np.nextafter(sd_upper, 0.0)).tolist()
        for radius, values in ((r, np.asarray(v, dtype=float)) for r, v in softdrop_mass.items())
    }
    figure, axis = plt.subplots(figsize=(10, 7.2))
    _probability_histogram(
        axis, clipped_sd, np.linspace(0.0, sd_upper, 46),
        r"GenFatJet Soft Drop mass $m_{SD}$ [GeV]",
        r"How GenFatJet radius changes $m_{SD}$",
    )
    axis.text(
        0.02, 0.97,
        r"anti-$k_T$ GenFatJets; Soft Drop $\beta=0$, $z_{cut}=0.1$"
        "\nSame fixed-truth jets at every R; overflow in final bin",
        transform=axis.transAxes, ha="left", va="top", fontsize=10,
    )
    outputs["softdrop_mass_overlay"] = _save_figure(
        figure, outdir / "softdrop_mass_by_radius.png"
    )

    figure, axis = plt.subplots(figsize=(10, 7.2))
    _probability_histogram(
        axis, count_containment, np.linspace(0.0, 1.0, 31),
        "Visible-descendant count containment",
        "Fraction of visible descendants captured by GenFatJet",
    )
    axis.axvline(1.0, color="black", linestyle="--", linewidth=1.2)
    axis.text(
        0.98, 0.04, "1.0 = every visible descendant is inside the matched jet",
        transform=axis.transAxes, ha="right", va="bottom", fontsize=10,
    )
    outputs["count_containment_overlay"] = _save_figure(
        figure, outdir / "visible_descendant_containment_by_radius.png"
    )

    figure, axis = plt.subplots(figsize=(10, 7.2))
    _probability_histogram(
        axis, pt_containment, np.linspace(0.0, 1.0, 31),
        r"Visible-descendant $p_T$ containment",
        r"Fraction of visible-descendant $p_T$ captured by GenFatJet",
    )
    axis.axvline(1.0, color="black", linestyle="--", linewidth=1.2)
    outputs["pt_containment_overlay"] = _save_figure(
        figure, outdir / "visible_descendant_pt_containment_by_radius.png"
    )

    max_full = max(
        1,
        int(max((max(values) for values in full_dark_hadron_count.values() if values), default=1)),
    )
    figure, axis = plt.subplots(figsize=(10, 7.2))
    _probability_histogram(
        axis,
        full_dark_hadron_count,
        np.arange(-0.5, max_full + 1.5, 1.0),
        "Number of fully contained target dark hadrons",
        "Fully contained dark-hadron multiplicity versus GenFatJet radius",
    )
    axis.text(
        0.02, 0.96,
        "Full = every resolved visible descendant is an exact jet constituent",
        transform=axis.transAxes, ha="left", va="top", fontsize=9.5,
    )
    outputs["fully_contained_dark_hadron_count"] = _save_figure(
        figure, outdir / "fully_contained_dark_hadron_count_by_radius.png"
    )

    max_not_full = max(
        1,
        int(max((max(values) for values in not_fully_contained_dark_hadron_count.values() if values), default=1)),
    )
    figure, axis = plt.subplots(figsize=(10, 7.2))
    _probability_histogram(
        axis,
        not_fully_contained_dark_hadron_count,
        np.arange(-0.5, max_not_full + 1.5, 1.0),
        "Number of target dark hadrons not fully contained",
        "Not-fully-contained dark-hadron multiplicity versus GenFatJet radius",
    )
    outputs["not_fully_contained_dark_hadron_count"] = _save_figure(
        figure, outdir / "not_fully_contained_dark_hadron_count_by_radius.png"
    )

    all_non_dark_counts = np.concatenate([
        np.asarray(values, dtype=float)
        for values in non_dark_hadron_constituent_count.values()
        if values
    ])
    non_dark_count_upper = max(
        10, int(math.ceil(float(np.percentile(all_non_dark_counts, 99.5)) / 5.0) * 5)
    )
    clipped_non_dark_count = {
        radius: np.minimum(values, non_dark_count_upper).tolist()
        for radius, values in (
            (r, np.asarray(v, dtype=float))
            for r, v in non_dark_hadron_constituent_count.items()
        )
    }
    figure, axis = plt.subplots(figsize=(10, 7.2))
    _probability_histogram(
        axis,
        clipped_non_dark_count,
        np.linspace(-0.5, non_dark_count_upper + 0.5, 42),
        "Number of GenFatJet constituents not descended from any initial dark hadron",
        "Non-dark-hadron particle multiplicity versus GenFatJet radius",
    )
    axis.text(
        0.98, 0.04, "Overflow is included in the final bin",
        transform=axis.transAxes, ha="right", va="bottom", fontsize=9.2,
    )
    outputs["non_dark_hadron_constituent_count"] = _save_figure(
        figure, outdir / "non_dark_hadron_constituent_count_by_radius.png"
    )

    figure, axes = plt.subplots(1, 2, figsize=(15, 6.4), sharey=True)
    _probability_histogram(
        axes[0], non_dark_hadron_constituent_fraction, np.linspace(0.0, 1.0, 41),
        "Fraction of constituents not from any initial dark hadron",
        "Particle-count contamination fraction",
    )
    _probability_histogram(
        axes[1], non_dark_hadron_constituent_pt_fraction,
        np.linspace(0.0, 1.0, 41),
        r"Fraction of constituent $p_T$ not from any initial dark hadron",
        r"Non-DH $p_T$ contamination fraction",
    )
    outputs["non_dark_hadron_constituent_fractions"] = _save_figure(
        figure, outdir / "non_dark_hadron_constituent_fractions_by_radius.png"
    )

    eligible_keys = {
        key for key, truth_jet in truth_by_key.items()
        if any(truth_jet.visible_by_dark_hadron)
    }
    multi_keys = {
        key for key, truth_jet in truth_by_key.items()
        if sum(bool(group) for group in truth_jet.visible_by_dark_hadron) > 1
    }
    acceptance_summary: dict[str, Any] = {}
    acceptance = []
    conditional_efficiency = []
    multi_acceptance = []
    multi_conditional_efficiency = []
    no_full_acceptance = []
    multi_no_full_acceptance = []
    has_not_full_fraction = []
    has_non_dark_fraction = []
    mean_non_dark_count_fraction = []
    mean_non_dark_pt_fraction = []
    for radius in radii:
        eligible_matched = eligible_keys & set(matched[radius])
        multi_matched = multi_keys & set(matched[radius])
        passing = 0
        multi_passing = 0
        has_not_full = 0
        has_non_dark = 0
        matched_count_fractions: list[float] = []
        matched_pt_fractions: list[float] = []
        for key in eligible_matched:
            event_index, truth_id = key
            metrics = dark_hadron_group_metrics(
                truth_event_by_event[event_index], truth_id,
                matched[radius][key], candidates_by_event[event_index],
            )
            passes = metrics["n_fully_contained_target_dark_hadrons"] >= 1
            passing += int(passes)
            multi_passing += int(passes and key in multi_keys)
            has_not_full += int(
                metrics["n_not_fully_contained_target_dark_hadrons"] >= 1
            )
            has_non_dark += int(metrics["has_non_dark_hadron_constituents"])
            matched_count_fractions.append(
                metrics["fraction_non_dark_hadron_constituents"]
            )
            matched_pt_fractions.append(
                metrics["non_dark_hadron_constituent_pt_fraction"]
            )
        acc = passing / len(eligible_keys) if eligible_keys else math.nan
        eff = passing / len(eligible_matched) if eligible_matched else math.nan
        multi_acc = multi_passing / len(multi_keys) if multi_keys else math.nan
        multi_eff = multi_passing / len(multi_matched) if multi_matched else math.nan
        acceptance.append(acc)
        conditional_efficiency.append(eff)
        multi_acceptance.append(multi_acc)
        multi_conditional_efficiency.append(multi_eff)
        no_full_acceptance.append(1.0 - acc)
        multi_no_full_acceptance.append(1.0 - multi_acc)
        has_not_full_value = (
            has_not_full / len(eligible_matched) if eligible_matched else math.nan
        )
        has_non_dark_value = (
            has_non_dark / len(eligible_matched) if eligible_matched else math.nan
        )
        mean_count_value = float(np.nanmean(matched_count_fractions))
        mean_pt_value = float(np.nanmean(matched_pt_fractions))
        has_not_full_fraction.append(has_not_full_value)
        has_non_dark_fraction.append(has_non_dark_value)
        mean_non_dark_count_fraction.append(mean_count_value)
        mean_non_dark_pt_fraction.append(mean_pt_value)
        acceptance_summary[f"{radius:.2f}"] = {
            "eligible_fixed_truth_jets": len(eligible_keys),
            "eligible_matched_truth_jets": len(eligible_matched),
            "pass_at_least_one_fully_contained_dark_hadron": passing,
            "acceptance": acc,
            "conditional_efficiency_given_match": eff,
            "eligible_multi_dark_hadron_truth_jets": len(multi_keys),
            "eligible_matched_multi_dark_hadron_truth_jets": len(multi_matched),
            "multi_dark_hadron_pass": multi_passing,
            "multi_dark_hadron_acceptance": multi_acc,
            "multi_dark_hadron_conditional_efficiency_given_match": multi_eff,
            "fraction_matched_with_not_fully_contained_target_dark_hadron": has_not_full_value,
            "fraction_matched_with_non_dark_hadron_constituent": has_non_dark_value,
            "mean_non_dark_hadron_constituent_count_fraction": mean_count_value,
            "mean_non_dark_hadron_constituent_pt_fraction": mean_pt_value,
        }

    figure, axes = plt.subplots(2, 1, figsize=(10.8, 10.0), sharex=True)
    axes[0].plot(radii, acceptance, marker="o", linewidth=2.2,
                 color="#2563EB", label="acceptance: ≥1 fully contained DH")
    axes[0].plot(radii, no_full_acceptance, marker="s", linewidth=2.0,
                 color="#D81B60", label="fails acceptance: no fully contained DH")
    axes[0].plot(radii, multi_acceptance, marker="^", linestyle="--", linewidth=2.0,
                 color="#7B2CBF", label="acceptance: >1 visible DH")
    axes[0].plot(radii, multi_no_full_acceptance, marker="D", linestyle="--",
                 linewidth=1.8, color="#E76F51",
                 label="fails acceptance: >1 visible DH")
    axes[0].set_ylabel("Fraction of eligible truth jets")
    axes[0].set_ylim(0.0, 1.02)
    axes[0].set_title("Dark-hadron containment acceptance")
    axes[0].grid(alpha=0.22)
    axes[0].legend(frameon=False, fontsize=8.8, ncol=2)
    axes[0].text(
        0.02, 0.04,
        "Acceptance denominator includes unmatched eligible truth jets",
        transform=axes[0].transAxes, ha="left", va="bottom", fontsize=9.2,
    )
    axes[1].plot(radii, has_not_full_fraction, marker="o", color="#7B2CBF",
                 label="≥1 target DH not fully contained")
    axes[1].plot(radii, has_non_dark_fraction, marker="s", color="#C62828",
                 label="≥1 constituent not from any initial dark hadron")
    axes[1].plot(radii, mean_non_dark_count_fraction, marker="^", color="#0086A8",
                 label="mean non-DH constituent-count fraction")
    axes[1].plot(radii, mean_non_dark_pt_fraction, marker="D", color="#F59E0B",
                 label=r"mean non-DH constituent-$p_T$ fraction")
    for axis in axes:
        axis.axvline(0.8, color="#111111", linestyle=":", linewidth=1.8,
                     label="default R=0.8" if axis is axes[1] else None)
        axis.set_ylim(0.0, 1.02)
        axis.grid(alpha=0.22)
    axes[1].set_xlabel("GenFatJet radius R")
    axes[1].set_ylabel("Fraction among matched eligible jets")
    axes[1].set_title("Not-contained target DH and non-DH particle contamination")
    axes[1].legend(frameon=False, fontsize=8.7, ncol=2)
    outputs["any_fully_contained_dark_hadron_acceptance"] = _save_figure(
        figure, outdir / "any_fully_contained_dark_hadron_acceptance_by_radius.png"
    )

    r08_not_full = np.asarray(
        not_fully_contained_dark_hadron_count[reference_radius], dtype=float
    )
    r08_non_dark = np.asarray(
        non_dark_hadron_constituent_count[reference_radius], dtype=float
    )
    r08_mass = np.asarray(softdrop_mass[reference_radius], dtype=float)
    r08_non_dark_pt = np.asarray(
        non_dark_hadron_constituent_pt_fraction[reference_radius], dtype=float
    )
    figure, axes = plt.subplots(1, 2, figsize=(15, 6.3))
    first = axes[0].hist2d(
        r08_not_full,
        np.minimum(r08_non_dark, non_dark_count_upper),
        bins=(np.arange(-0.5, max_not_full + 1.5, 1.0), 40),
        norm=LogNorm(), cmap="magma",
    )
    figure.colorbar(first[3], ax=axes[0], label="Matched truth jets / bin")
    axes[0].set_xlabel("Number of target DH not fully contained")
    axes[0].set_ylabel("Number of non-DH constituents (overflow clipped)")
    axes[0].set_title(f"Containment versus particle contamination, R={reference_radius:.1f}")
    second = axes[1].hist2d(
        r08_non_dark_pt, r08_mass,
        bins=(np.linspace(0.0, 1.0, 41), 45), norm=LogNorm(), cmap="plasma",
    )
    figure.colorbar(second[3], ax=axes[1], label="Matched truth jets / bin")
    axes[1].set_xlabel(r"Non-DH constituent $p_T$ fraction")
    axes[1].set_ylabel(r"GenFatJet $m_{SD}$ [GeV]")
    axes[1].set_title(f"Soft Drop mass versus non-DH contamination, R={reference_radius:.1f}")
    for axis in axes:
        axis.grid(alpha=0.15)
    outputs["containment_contamination_correlations"] = _save_figure(
        figure, outdir / "containment_contamination_correlations_R08.png"
    )

    max_contaminating = max(
        1,
        int(max((max(values) for values in contaminating_dark_hadron_count.values() if values), default=1)),
    )
    figure, axes = plt.subplots(1, 2, figsize=(15, 6.4), sharey=True)
    bins = np.arange(-0.5, max_contaminating + 1.5, 1.0)
    _probability_histogram(
        axes[0], contaminating_dark_hadron_count, bins,
        "Number of contaminating dark hadrons",
        "Any visible descendant enters target jet",
    )
    _probability_histogram(
        axes[1], fully_contained_contaminating_count, bins,
        "Number fully contained contaminating dark hadrons",
        "All visible descendants enter target jet",
    )
    outputs["dark_hadron_contamination_counts"] = _save_figure(
        figure, outdir / "dark_hadron_contamination_counts_by_radius.png"
    )

    figure, axes = plt.subplots(2, 1, figsize=(10.5, 9.0), sharex=True)
    scan_radii = np.asarray(sorted(float(key) for key in scan["radii"]), dtype=float)
    for metric, label, color in (
        ("containment_count", "visible count containment", "tab:blue"),
        ("containment_pt", r"visible $p_T$ containment", "#7B2CBF"),
        ("contamination_pt", r"non-truth $p_T$ contamination", "tab:red"),
    ):
        entries = [scan["radii"][f"{r:.2f}"][metric] for r in scan_radii]
        mean = np.asarray([entry["mean"] for entry in entries], dtype=float)
        lower = np.asarray([entry["lower"] for entry in entries], dtype=float)
        upper = np.asarray([entry["upper"] for entry in entries], dtype=float)
        axes[0].plot(scan_radii, mean, marker="o", label=label, color=color)
        axes[0].fill_between(scan_radii, lower, upper, color=color, alpha=0.15)
    axes[0].axvline(
        0.8, color="#111111", linestyle="--", linewidth=1.6,
        label="default R=0.8",
    )
    axes[0].set_ylabel("Mean fraction")
    axes[0].set_ylim(0.0, 1.02)
    axes[0].legend(frameon=False, ncol=2)
    axes[0].grid(alpha=0.22)
    efficiency = []
    split_fraction = []
    merge_fraction = []
    for radius in scan_radii:
        entry = scan["radii"][f"{radius:.2f}"]
        efficiency.append(entry["match_efficiency"])
        split_fraction.append(entry["split_truth_jets"] / entry["truth_jets"])
        merge_fraction.append(
            entry["merge_clustered_jets"] / entry["clustered_jets"]
            if entry["clustered_jets"] else math.nan
        )
    axes[1].plot(scan_radii, efficiency, marker="o", label="truth match efficiency")
    axes[1].plot(scan_radii, split_fraction, marker="s", label="truth split fraction")
    axes[1].plot(
        scan_radii, merge_fraction, marker="^", color="#D81B60",
        label="clustered-jet merge fraction",
    )
    axes[1].axvline(0.8, color="#111111", linestyle="--", linewidth=1.6)
    axes[1].set_xlabel("GenFatJet radius R")
    axes[1].set_ylabel("Fraction")
    axes[1].set_ylim(0.0, 1.02)
    axes[1].grid(alpha=0.22)
    axes[1].legend(frameon=False, ncol=2)
    figure.suptitle("Radius tradeoffs for fixed DarkHadronJet truth objects", fontsize=16)
    figure.tight_layout()
    outputs["radius_summary"] = _save_figure(
        figure, outdir / "radius_performance_summary.png"
    )

    all_collection_mass = np.concatenate(
        [np.asarray(values, dtype=float) for values in collection_mass.values() if values]
    )
    mass_upper = max(
        150.0,
        math.ceil(float(np.percentile(all_collection_mass, 99.0)) / 25.0) * 25.0,
    )
    figure, axis = plt.subplots(figsize=(10, 7.2))
    colors = ("#2563EB", "#7B2CBF", "#F59E0B", "#D62728")
    for color, (label, entries) in zip(colors, collection_mass.items()):
        values = np.asarray(entries, dtype=float)
        values = values[np.isfinite(values)]
        values = np.minimum(values, np.nextafter(mass_upper, 0.0))
        axis.hist(
            values, bins=np.linspace(0.0, mass_upper, 46),
            weights=np.full(len(values), 1.0 / len(values)), histtype="step",
            linewidth=2.0, color=color, label=f"{label} (N={len(values)})",
        )
    axis.set_xlabel("Jet mass [GeV]")
    axis.set_ylabel("Fraction of jets / bin")
    axis.set_yscale("log", nonpositive="clip")
    axis.set_ylim(bottom=1.0e-5)
    axis.set_title("Mass distributions of the four study collections")
    axis.grid(alpha=0.22, which="both")
    axis.legend(frameon=False)
    axis.text(
        0.98, 0.96,
        "Dark collections: all nonzero stored jets\nGenFatJets: common truth-matched R=0.8 cohort",
        transform=axis.transAxes, ha="right", va="top", fontsize=9.5,
    )
    outputs["collection_mass_comparison"] = _save_figure(
        figure, outdir / "collection_mass_comparison.png"
    )

    radius_summary: dict[str, Any] = {}
    for radius in radii:
        sd = np.asarray(softdrop_mass[radius], dtype=float)
        count = np.asarray(count_containment[radius], dtype=float)
        pt = np.asarray(pt_containment[radius], dtype=float)
        radius_summary[f"{radius:.2f}"] = {
            "softdrop_mass_median_GeV": float(np.nanmedian(sd)),
            "softdrop_mass_16_84_GeV": np.nanpercentile(sd, [16, 84]).tolist(),
            "mean_count_containment": float(np.nanmean(count)),
            "fully_count_contained_fraction": float(np.nanmean(count >= 0.999999)),
            "mean_pt_containment": float(np.nanmean(pt)),
        }
    return {
        "passed": True,
        "definition": (
            "Common fixed-truth cohort matched at every displayed R; GenFatJet "
            "pT>15 GeV; Soft Drop beta=0 and zcut=0.1"
        ),
        "common_truth_jets": len(common_keys),
        "diagnostic_radii": radii,
        "radii": radius_summary,
        "any_fully_contained_dark_hadron_acceptance": acceptance_summary,
        "plots": outputs,
    }


def make_event_displays(
    events: Any,
    raw_ancestry: Sequence[Any],
    radii: Sequence[float],
    matched: dict[float, dict[tuple[int, int], ClusteredJet]],
    outdir: Path,
    max_displays: int,
) -> dict[str, Any]:
    """Render several events with every requested radius overlaid per image."""

    useful_radii = sorted(set(float(radius) for radius in radii))
    common = set.intersection(*(set(matched[radius]) for radius in useful_radii))
    if not common:
        return {"passed": False, "reason": "no truth jet matched at display radii"}

    # Rank jets by how much visible-descendant containment changes across the
    # radius range. Keep at most one jet per event so the displays show genuinely
    # different events rather than two jets from the same collision.
    scored: list[tuple[float, tuple[int, int]]] = []
    for key in sorted(common):
        event_index, truth_id = key
        candidates = particle_table(events[event_index].GenCandidate)
        truth_jets = build_fixed_truth_event(
            events[event_index], raw_ancestry[event_index]
        ).jets
        truth_jet = truth_jets[truth_id]
        truth_set = set(truth_jet.visible_indices)
        if len(truth_set) < 8:
            continue
        low = len(
            truth_set & set(matched[useful_radii[0]][key].constituent_indices)
        )
        high = len(
            truth_set & set(matched[useful_radii[-1]][key].constituent_indices)
        )
        containment_change = (high - low) / len(truth_set)
        truth_pt = float(np.sum(candidates.pt[list(truth_set)]))
        n_visible_dark_hadrons = sum(
            bool(group) for group in truth_jet.visible_by_dark_hadron
        )
        multi_dark_hadron_bonus = 10.0 if n_visible_dark_hadrons > 1 else 0.0
        scored.append(
            (multi_dark_hadron_bonus + containment_change + 1.0e-6 * truth_pt, key)
        )

    selected: list[tuple[int, int]] = []
    used_events: set[int] = set()
    for _, key in sorted(scored, reverse=True):
        if key[0] in used_events:
            continue
        selected.append(key)
        used_events.add(key[0])
        if len(selected) >= max_displays:
            break
    if not selected:
        return {"passed": False, "reason": "no display jet has enough descendants"}

    reference_radius = min(useful_radii, key=lambda radius: abs(radius - 0.8))
    zoom = max(1.35, max(useful_radii) + 0.35)
    agreement = True
    displays: list[dict[str, Any]] = []

    for display_number, selected_key in enumerate(selected, start=1):
        event_index, truth_id = selected_key
        event = events[event_index]
        candidates = particle_table(event.GenCandidate)
        truth_event = build_fixed_truth_event(event, raw_ancestry[event_index])
        truth = truth_event.jets[truth_id]
        truth_set = set(truth.visible_indices)
        truth_mask = np.zeros(len(candidates), dtype=bool)
        truth_mask[list(truth_set)] = True
        truth_pt = float(np.sum(candidates.pt[list(truth_set)]))
        dphi = wrap_delta_phi(candidates.phi, truth.phi)
        deta = candidates.eta - truth.eta
        nearby = (np.abs(dphi) < zoom) & (np.abs(deta) < zoom)
        sizes = np.clip(np.sqrt(candidates.pt) * 8.0, 8.0, 140.0)

        reference_jet = matched[reference_radius][selected_key]
        reference_members = np.zeros(len(candidates), dtype=bool)
        reference_members[list(reference_jet.constituent_indices)] = True
        contamination_reference = reference_members & (~truth_mask)
        background = nearby & (~truth_mask) & (~reference_members)

        figure, axis = plt.subplots(figsize=(15.5, 10.0))
        axis.scatter(
            dphi[background], deta[background], s=sizes[background],
            c="0.82", alpha=0.32, linewidths=0, label="other visible particles",
        )
        axis.scatter(
            dphi[contamination_reference], deta[contamination_reference],
            s=sizes[contamination_reference] + 20.0, facecolors="none",
            edgecolors="#C62828", linewidths=1.4,
            label=f"non-truth constituents at R={reference_radius:.1f}",
        )
        axis.scatter(
            0.0, 0.0, marker="*", s=210, c="gold", edgecolors="black",
            zorder=8, label="fixed DarkHadronJet truth axis",
        )

        dark_hadrons = particle_table(event.DarkHadronCandidate)
        dark_hadron_uid_index = {
            int(uid): index for index, uid in enumerate(dark_hadrons.uid)
        }
        target_dark_hadron_uids = set(truth.dark_hadron_uids)
        for position, (uid, group) in enumerate(
            zip(truth.dark_hadron_uids, truth.visible_by_dark_hadron)
        ):
            color = DARK_HADRON_COLORS[position % len(DARK_HADRON_COLORS)]
            group_indices = np.asarray(group, dtype=int)
            if len(group_indices):
                captured = group_indices[reference_members[group_indices]]
                missed = group_indices[~reference_members[group_indices]]
                if len(captured):
                    axis.scatter(
                        dphi[captured], deta[captured], s=sizes[captured],
                        c=color, alpha=0.86, edgecolors="white", linewidths=0.35,
                        label=f"DH#{position} descendants captured ({len(captured)})",
                        zorder=4,
                    )
                if len(missed):
                    axis.scatter(
                        dphi[missed], deta[missed], s=sizes[missed],
                        c=color, marker="x", linewidths=1.6,
                        label=f"DH#{position} descendants missed ({len(missed)})",
                        zorder=5,
                    )
            dh_index = dark_hadron_uid_index.get(int(uid))
            if dh_index is None:
                continue
            dh_x = float(wrap_delta_phi(dark_hadrons.phi[dh_index], truth.phi))
            dh_y = float(dark_hadrons.eta[dh_index] - truth.eta)
            axis.scatter(
                dh_x, dh_y, marker="*", s=330, c=color,
                edgecolors="black", linewidths=1.1, zorder=9,
                label=f"initial DH#{position} (pT={dark_hadrons.pt[dh_index]:.0f} GeV)",
            )
            axis.annotate(
                f"DH#{position}", (dh_x, dh_y), xytext=(6, 6),
                textcoords="offset points", color=color, fontsize=8.5,
                fontweight="bold", zorder=10,
            )

        # Show other initial dark hadrons in the same event faintly. Those that
        # actually feed the R=0.8 target jet are emphasized with red outlines.
        reference_dh_metrics = dark_hadron_group_metrics(
            truth_event, truth_id, reference_jet, candidates
        )
        contaminating_uids = {
            record["dark_hadron_uid"]
            for record in reference_dh_metrics["contaminating_dark_hadrons"]
        }
        other_indices = [
            index for index, uid in enumerate(dark_hadrons.uid)
            if int(uid) not in target_dark_hadron_uids
            and abs(float(wrap_delta_phi(dark_hadrons.phi[index], truth.phi))) < zoom
            and abs(float(dark_hadrons.eta[index] - truth.eta)) < zoom
        ]
        if other_indices:
            other_x = np.asarray([
                float(wrap_delta_phi(dark_hadrons.phi[index], truth.phi))
                for index in other_indices
            ])
            other_y = np.asarray([
                float(dark_hadrons.eta[index] - truth.eta) for index in other_indices
            ])
            other_is_contaminating = np.asarray([
                int(dark_hadrons.uid[index]) in contaminating_uids
                for index in other_indices
            ])
            if np.any(~other_is_contaminating):
                axis.scatter(
                    other_x[~other_is_contaminating], other_y[~other_is_contaminating],
                    marker="D", s=85, facecolors="none", edgecolors="0.55",
                    linewidths=1.0, alpha=0.65, label="other initial dark hadrons",
                    zorder=6,
                )
            if np.any(other_is_contaminating):
                axis.scatter(
                    other_x[other_is_contaminating], other_y[other_is_contaminating],
                    marker="D", s=140, facecolors="#FEE2E2", edgecolors="#C62828",
                    linewidths=1.8,
                    label=f"contaminating DH parent at R={reference_radius:.1f}",
                    zorder=8,
                )

        rows: list[list[str]] = []
        per_radius: dict[str, Any] = {}
        for radius in useful_radii:
            jet = matched[radius][selected_key]
            rendered_members = np.zeros(len(candidates), dtype=bool)
            rendered_members[list(jet.constituent_indices)] = True
            calculation_members = np.zeros(len(candidates), dtype=bool)
            calculation_members[np.asarray(jet.constituent_indices, dtype=int)] = True
            exact = np.array_equal(rendered_members, calculation_members)
            agreement &= exact

            captured_truth = rendered_members & truth_mask
            contamination = rendered_members & (~truth_mask)
            count_fraction = float(np.sum(captured_truth)) / len(truth_set)
            captured_pt = float(np.sum(candidates.pt[captured_truth]))
            pt_fraction = captured_pt / truth_pt if truth_pt > 0.0 else math.nan
            jet_pt = float(np.sum(candidates.pt[rendered_members]))
            contamination_pt = float(np.sum(candidates.pt[contamination]))
            contamination_fraction = (
                contamination_pt / jet_pt if jet_pt > 0.0 else math.nan
            )
            dh_metrics = dark_hadron_group_metrics(
                truth_event, truth_id, jet, candidates
            )
            sd_mass = softdrop_jet(
                candidates, jet, beta=0.0, zcut=0.1, r0=radius
            ).mass
            jet_dphi = float(wrap_delta_phi(jet.phi, truth.phi))
            jet_deta = jet.eta - truth.eta
            color = radius_color(radius)
            axis.add_patch(
                plt.Circle(
                    (jet_dphi, jet_deta), radius, fill=False, color=color,
                    linestyle="-" if math.isclose(radius, 0.8) else "--",
                    linewidth=2.4 if math.isclose(radius, 0.8) else 1.45,
                    alpha=0.92 if math.isclose(radius, 0.8) else 0.72,
                    label=f"GenFatJet R={radius:.1f}",
                )
            )
            axis.scatter(
                jet_dphi, jet_deta, marker="+", s=90, c=color,
                linewidths=1.6, zorder=7,
            )
            rows.append(
                [
                    f"{radius:.1f}", f"{count_fraction:.1%}",
                    f"{pt_fraction:.1%}",
                    (
                        f"{dh_metrics['n_fully_contained_target_dark_hadrons']}/"
                        f"{dh_metrics['n_target_dark_hadrons_with_visible_descendants']}"
                    ),
                    str(dh_metrics["n_contaminating_dark_hadrons"]),
                    str(dh_metrics["n_fully_contained_contaminating_dark_hadrons"]),
                    f"{contamination_fraction:.1%}",
                    f"{sd_mass:.1f}",
                ]
            )
            per_radius[f"{radius:.2f}"] = {
                "rendered_constituent_indices": np.flatnonzero(
                    rendered_members
                ).tolist(),
                "calculation_constituent_indices": list(jet.constituent_indices),
                "bitwise_equal": exact,
                "geometric_circle_is_guide_only": True,
                "count_containment": count_fraction,
                "pt_containment": pt_fraction,
                "contamination_pt_fraction": contamination_fraction,
                "softdrop_mass_GeV": sd_mass,
                "dark_hadron_ancestry": dh_metrics,
            }

        axis.set_title(
            f"Event {event_index}, fixed DarkHadronJet truth ID {truth_id} "
            f"({sum(bool(group) for group in truth.visible_by_dark_hadron)} visible DH)\n"
            "Each dark hadron and its descendants share a color; all radii overlaid",
            fontsize=15,
        )
        axis.set_xlabel(r"$\Delta\phi$ from fixed DarkHadronJet axis")
        axis.set_ylabel(r"$\Delta\eta$ from fixed DarkHadronJet axis")
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlim(-zoom, zoom)
        axis.set_ylim(-zoom, zoom)
        axis.grid(alpha=0.18)
        axis.legend(
            loc="lower left", fontsize=6.8, framealpha=0.92, ncol=2,
        )
        table = axis.table(
            cellText=rows,
            colLabels=[
                "R", "desc.", r"desc. $p_T$", "target DH\nfull",
                "other DH", "other DH\nfull", r"all contam. $p_T$", r"$m_{SD}$",
            ],
            cellLoc="center", colLoc="center", bbox=[1.03, 0.18, 0.98, 0.64],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8.5)
        axis.text(
            1.04, 0.86,
            "Exact constituent and ancestry metrics\nfor the same fixed truth jet",
            transform=axis.transAxes, ha="left", va="bottom", fontsize=10.5,
            fontweight="bold",
        )
        axis.text(
            0.98, 0.02,
            "Circles are geometric guides; table values use exact constituents. "
            "R=0.8 is black.",
            transform=axis.transAxes, ha="right", va="bottom", fontsize=9.0,
        )
        output = outdir / (
            f"event_display_{display_number:02d}_event_{event_index}_truth_{truth_id}.png"
        )
        saved = _save_figure(figure, output)
        displays.append(
            {
                "event": event_index,
                "truth_id": truth_id,
                "output": saved,
                "radii": per_radius,
            }
        )

    return {
        "passed": agreement,
        "description": (
            "one fixed truth jet per event image; all diagnostic GenFatJet "
            "radii overlaid; each initial dark hadron and its visible "
            "descendants share a color"
        ),
        "requested_displays": max_displays,
        "produced_displays": len(displays),
        "diagnostic_radii": useful_radii,
        "displays": displays,
    }


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    report_path = args.outdir / "validation_report.json"
    report: dict[str, Any] = {
        "status": "running",
        "scope": {
            "local_only": True,
            "training_data": False,
            "campaign": False,
        },
    }
    try:
        events = load_events(args.input, args.max_events)
        raw_ancestry = load_raw_ancestry(args.input, args.max_events)
        n_events = len(events)
        if n_events == 0:
            raise RuntimeError("The ROOT file contains no events")
        report["metadata"] = study_metadata(
            args.input,
            REPOSITORY,
            args.delphes,
            extra={"events_read": n_events},
        )
        report["card_roundtrip"] = check_card(args.card, args.collections)
        report["default_radius_closure"] = {
            collection: check_collection_closure(events, collection, n_events)
            for collection in args.collections
        }
        if len(raw_ancestry) != n_events:
            raise RuntimeError(
                f"NanoEvents/raw genealogy event mismatch: {n_events}/"
                f"{len(raw_ancestry)}"
            )
        report["truth_ancestry_closure"] = check_truth_ancestry(
            events, raw_ancestry, n_events
        )
        report["softdrop_default_radius_closure"] = check_softdrop_closure(
            events, n_events
        )
        scan_radii = sorted(set(args.radii) | set(args.diagnostic_radii))
        scan, matched = run_primary_scan(
            events,
            raw_ancestry,
            n_events,
            scan_radii,
            args.bootstrap_resamples,
        )
        report["primary_radius_scan"] = scan
        report["diagnostics"] = make_diagnostic_plots(
            events,
            raw_ancestry,
            n_events,
            args.diagnostic_radii,
            matched,
            scan,
            args.outdir,
        )
        report["event_display_agreement"] = make_event_displays(
            events,
            raw_ancestry,
            args.diagnostic_radii,
            matched,
            args.outdir,
            args.event_displays,
        )
        gates = {
            "card_roundtrip": bool(report["card_roundtrip"]["passed"]),
            "default_radius_closure": all(
                bool(result["passed"])
                for result in report["default_radius_closure"].values()
            ),
            "truth_ancestry_closure": bool(report["truth_ancestry_closure"]["passed"]),
            "softdrop_default_radius_closure": bool(
                report["softdrop_default_radius_closure"]["passed"]
            ),
            "primary_scan_invariants": bool(report["primary_radius_scan"]["passed"]),
            "diagnostic_plots": bool(report["diagnostics"]["passed"]),
            "event_display_agreement": bool(report["event_display_agreement"]["passed"]),
        }
        report["gates"] = gates
        report["status"] = "passed" if all(gates.values()) else "failed"
    except Exception as error:
        report["status"] = "error"
        report["error"] = {
            "type": type(error).__name__,
            "message": str(error),
            "traceback": traceback.format_exc(),
        }
    write_json(report_path, report)
    print(json.dumps(json_ready({
        "status": report["status"],
        "report": str(report_path.resolve()),
        "gates": report.get("gates"),
        "error": report.get("error", {}).get("message"),
    }), indent=2))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
