"""Physics utilities for the local jet-radius study.

This module works directly from Delphes ROOT objects.  It deliberately contains
no dataset-production, sampling, tensor, NPZ, batch, or upload machinery.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import subprocess
from typing import Any, Iterable, Mapping, Sequence

import awkward as ak
import fastjet as fj
import numpy as np

from ancestry import (
    ancestor_masks,
    descends_from_any,
    targets_for_particle,
)
from scipy.optimize import linear_sum_assignment
import uproot


DEFAULT_R_GRID = tuple(round(float(r), 2) for r in np.arange(0.2, 1.61, 0.1))
DEFAULT_PERCENTILES = (50.0, 68.0, 90.0, 95.0, 99.0)


@dataclass(frozen=True)
class CollectionSpec:
    """A Delphes jet branch and the candidate pool from which it was built."""

    jet_branch: str
    candidate_branch: str
    module: str
    default_r: float
    pt_min: float = 15.0
    truth_link: str = "exact"
    four_vector_policy: str = "raw_cluster"


COLLECTIONS: dict[str, CollectionSpec] = {
    "GenJet": CollectionSpec("GenJet", "GenCandidate", "GenJetFinder", 0.4),
    "GenFatJet": CollectionSpec(
        "GenFatJet", "GenCandidate", "GenFatJetFinder", 0.8
    ),
    "DarkPartonJet": CollectionSpec(
        "DarkPartonJet", "DarkPartonCandidate", "DarkPartonJetFinder", 0.8
    ),
    "DarkHadronJet": CollectionSpec(
        "DarkHadronJet", "DarkHadronCandidate", "DarkHadronJetFinder", 0.8
    ),
    "Jet": CollectionSpec(
        "Jet", "ParticleFlowCandidate", "FastJetFinder", 0.4,
        truth_link="pf_refs", four_vector_policy="post_energy_scale"
    ),
    "FatJet": CollectionSpec(
        "FatJet", "ParticleFlowCandidate", "FatJetFinder", 0.8,
        truth_link="pf_refs",
    ),
}


@dataclass(frozen=True)
class ParticleTable:
    pt: np.ndarray
    eta: np.ndarray
    phi: np.ndarray
    energy: np.ndarray
    px: np.ndarray
    py: np.ndarray
    pz: np.ndarray
    uid: np.ndarray

    def __len__(self) -> int:
        return len(self.pt)


@dataclass(frozen=True)
class ClusteredJet:
    pt: float
    eta: float
    phi: float
    mass: float
    energy: float
    px: float
    py: float
    pz: float
    constituent_indices: tuple[int, ...]


@dataclass(frozen=True)
class SoftDropResult:
    pt: float
    eta: float
    phi: float
    mass: float
    n_dropped_branches: int
    passed_two_prong: bool


@dataclass
class FixedTruthJet:
    """One shipped R=0.8 dark-hadron partition and its visible descendants."""

    truth_id: int
    dark_hadron_uids: tuple[int, ...]
    dark_hadron_pts: tuple[float, ...]
    visible_indices: tuple[int, ...]
    visible_by_dark_hadron: tuple[tuple[int, ...], ...]
    pt: float
    eta: float
    phi: float
    mass: float


@dataclass
class FixedTruthEvent:
    jets: list[FixedTruthJet]
    visible_owner: np.ndarray
    visible_dark_hadron: np.ndarray
    visible_has_any_dark_hadron_ancestor: np.ndarray
    unresolved_visible: tuple[int, ...] = ()
    cache_reassigned_visible: tuple[int, ...] = ()


@dataclass(frozen=True)
class RawAncestryEvent:
    """Unmodified integer generator genealogy read directly with uproot."""

    uid: np.ndarray
    mother1: np.ndarray
    mother2: np.ndarray


@dataclass(frozen=True)
class Match:
    truth_id: int
    clustered_index: int
    shared_visible_pt: float
    truth_visible_pt: float
    clustered_pt: float

    @property
    def shared_fraction(self) -> float:
        if self.truth_visible_pt <= 0.0:
            return math.nan
        return self.shared_visible_pt / self.truth_visible_pt


def _numpy(values: Any, dtype: Any = float) -> np.ndarray:
    return np.asarray(ak.to_numpy(values), dtype=dtype)


def particle_table(branch: Any) -> ParticleTable:
    """Convert one event's Delphes candidate branch to plain four-vectors."""

    fields = set(ak.fields(branch))
    pt = _numpy(branch.PT)
    eta = _numpy(branch.Eta)
    phi = _numpy(branch.Phi)
    energy = _numpy(branch.E)
    if {"Px", "Py", "Pz"}.issubset(fields):
        px = _numpy(branch.Px)
        py = _numpy(branch.Py)
        pz = _numpy(branch.Pz)
    else:
        px = pt * np.cos(phi)
        py = pt * np.sin(phi)
        pz = pt * np.sinh(eta)
    uid = _numpy(branch.fUniqueID, np.int64)
    return ParticleTable(pt, eta, phi, energy, px, py, pz, uid)


def wrap_delta_phi(phi1: Any, phi2: Any) -> Any:
    return (np.asarray(phi1) - np.asarray(phi2) + np.pi) % (2.0 * np.pi) - np.pi


def delta_r2(eta1: Any, phi1: Any, eta2: Any, phi2: Any) -> Any:
    return (np.asarray(eta1) - np.asarray(eta2)) ** 2 + wrap_delta_phi(
        phi1, phi2
    ) ** 2


def inside_indices(
    particles: ParticleTable, axis_eta: float, axis_phi: float, radius: float
) -> np.ndarray:
    """Indices inside an eta-phi circle; shared by calculations and displays."""

    return np.flatnonzero(
        delta_r2(particles.eta, particles.phi, axis_eta, axis_phi)
        < float(radius) ** 2
    )


def recluster_event(
    particles: ParticleTable, radius: float, pt_min: float = 15.0
) -> list[ClusteredJet]:
    """Anti-kT cluster one event and immediately materialize constituents."""

    pseudojets: list[fj.PseudoJet] = []
    for idx, (px, py, pz, energy) in enumerate(
        zip(particles.px, particles.py, particles.pz, particles.energy)
    ):
        pseudojet = fj.PseudoJet(float(px), float(py), float(pz), float(energy))
        pseudojet.set_user_index(idx)
        pseudojets.append(pseudojet)
    if not pseudojets:
        return []

    definition = fj.JetDefinition(fj.antikt_algorithm, float(radius))
    sequence = fj.ClusterSequence(pseudojets, definition)
    fastjet_jets = fj.sorted_by_pt(sequence.inclusive_jets(float(pt_min)))
    output: list[ClusteredJet] = []
    for jet in fastjet_jets:
        indices = tuple(sorted(part.user_index() for part in jet.constituents()))
        output.append(
            ClusteredJet(
                pt=float(jet.pt()),
                eta=float(jet.eta()),
                phi=float(jet.phi_std()),
                mass=float(jet.m()),
                energy=float(jet.e()),
                px=float(jet.px()),
                py=float(jet.py()),
                pz=float(jet.pz()),
                constituent_indices=indices,
            )
        )
    return output


def cluster_input(
    events_px: Sequence[np.ndarray],
    events_py: Sequence[np.ndarray],
    events_pz: Sequence[np.ndarray],
    events_energy: Sequence[np.ndarray],
) -> Any:
    """Build the awkward record that :func:`recluster_events` clusters.

    Kept separate from the clustering call on purpose.  Assembling this array
    costs about as much as one clustering pass, so building it once and reusing
    it across the whole radius grid is what actually makes the batched path
    pay -- rebuilding it per radius gives back most of the gain.

    No dependency on the ``vector`` package: FastJet's awkward interface only
    needs the px/py/pz/E field names and a list-of-records layout.  The
    ``Momentum4D`` name is set for readability, and clustering gives identical
    results without it.
    """

    return ak.zip(
        {
            "px": ak.Array([np.asarray(x, dtype=np.float64) for x in events_px]),
            "py": ak.Array([np.asarray(x, dtype=np.float64) for x in events_py]),
            "pz": ak.Array([np.asarray(x, dtype=np.float64) for x in events_pz]),
            "E": ak.Array([np.asarray(x, dtype=np.float64) for x in events_energy]),
        },
        with_name="Momentum4D",
    )


def recluster_events(
    particles: Any, radius: float, pt_min: float = 15.0
) -> list[list[np.ndarray]]:
    """Anti-kT cluster MANY events in one FastJet call.

    Batched counterpart of :func:`recluster_event`, which builds one Python
    ``fj.PseudoJet`` per particle per event per radius.  ``particles`` is the
    record from :func:`cluster_input`; build it once and pass it to every
    radius.

    Measured on 60 events of the CMS validation sample over a 15-point radius
    grid, with the input built once: 2.82 -> 1.35 ms/event/radius, identical
    constituent partitions across all 900 (radius, event) combinations.

    Returns one list of constituent-index arrays per event.  Only indices come
    back, not four-vectors: every downstream quantity in this module is
    computed from exact constituent membership.

    :func:`recluster_event` is deliberately left in place -- the closure gates
    in ``validate.py`` check offline clustering against stored Delphes
    membership, and that comparison should keep exercising the simple path.
    """

    sequence = fj.ClusterSequence(
        particles, fj.JetDefinition(fj.antikt_algorithm, float(radius))
    )
    constituents = ak.to_list(sequence.constituent_index(min_pt=float(pt_min)))
    return [
        [np.asarray(group, dtype=np.int64) for group in event]
        for event in constituents
    ]


def softdrop_jet(
    particles: ParticleTable,
    jet: ClusteredJet,
    beta: float = 0.0,
    zcut: float = 0.1,
    r0: float = 0.8,
) -> SoftDropResult:
    """Apply the Delphes/FastJet-contrib Soft Drop definition to one jet.

    FastJet contrib reclusters the input constituents with C/A at its maximum
    allowable radius before walking backward through the hardest branch.  This
    implementation follows that behavior and uses the scalar-z condition.
    """

    pseudojets: list[fj.PseudoJet] = []
    for index in jet.constituent_indices:
        pseudojet = fj.PseudoJet(
            float(particles.px[index]),
            float(particles.py[index]),
            float(particles.pz[index]),
            float(particles.energy[index]),
        )
        pseudojet.set_user_index(int(index))
        pseudojets.append(pseudojet)
    if not pseudojets:
        return SoftDropResult(0.0, 0.0, 0.0, 0.0, 0, False)

    definition = fj.JetDefinition(
        fj.cambridge_algorithm, fj.JetDefinition.max_allowable_R
    )
    sequence = fj.ClusterSequence(pseudojets, definition)
    current = fj.sorted_by_pt(sequence.inclusive_jets(0.0))[0]
    dropped = 0
    passed = False
    while True:
        parent1 = fj.PseudoJet()
        parent2 = fj.PseudoJet()
        if not current.has_parents(parent1, parent2):
            break
        if parent2.pt2() > parent1.pt2():
            parent1, parent2 = parent2, parent1
        denominator = parent1.pt() + parent2.pt()
        scalar_z = min(parent1.pt(), parent2.pt()) / denominator if denominator else 0.0
        angular_factor = (
            parent1.squared_distance(parent2) / (float(r0) ** 2)
        ) ** (0.5 * float(beta))
        if scalar_z > float(zcut) * angular_factor:
            passed = True
            break
        current = parent1
        dropped += 1
    return SoftDropResult(
        pt=float(current.pt()),
        eta=float(current.eta()),
        phi=float(current.phi_std()),
        mass=float(current.m()),
        n_dropped_branches=dropped,
        passed_two_prong=passed,
    )


def stored_jets(branch: Any, candidates: ParticleTable) -> list[ClusteredJet]:
    """Materialize a stored Delphes jet branch and resolve TRef constituents."""

    candidate_index = {int(uid): idx for idx, uid in enumerate(candidates.uid)}
    refs = ak.to_list(branch.Constituents.refs)
    output: list[ClusteredJet] = []
    for idx, jet_refs in enumerate(refs):
        indices = tuple(
            sorted(candidate_index[int(uid)] for uid in jet_refs if int(uid) in candidate_index)
        )
        pt = float(branch.PT[idx])
        eta = float(branch.Eta[idx])
        phi = float(branch.Phi[idx])
        mass = float(branch.Mass[idx])
        px = pt * math.cos(phi)
        py = pt * math.sin(phi)
        pz = pt * math.sinh(eta)
        energy = math.sqrt(max(px * px + py * py + pz * pz + mass * mass, 0.0))
        output.append(
            ClusteredJet(
                pt, eta, phi, mass, energy, px, py, pz, indices
            )
        )
    return output


def load_raw_ancestry(
    path: str | Path,
    entry_stop: int | None = None,
    entry_start: int | None = None,
) -> list[RawAncestryEvent]:
    """Read genealogy outside NanoEvents to preserve signed integer M1/M2.

    Coffea's generic Delphes schema can reinterpret ``GenParticle.M2`` when it
    builds vector records.  Ancestry is discrete physics state, so the study
    reads these three branches directly and refuses lossy values.
    """

    branches = (
        "GenParticle.fUniqueID",
        "GenParticle.M1",
        "GenParticle.M2",
    )
    with uproot.open(path) as root_file:
        arrays = root_file["Delphes"].arrays(
            branches, entry_start=entry_start, entry_stop=entry_stop
        )
    output: list[RawAncestryEvent] = []
    for uid, mother1, mother2 in zip(
        arrays[branches[0]], arrays[branches[1]], arrays[branches[2]]
    ):
        output.append(
            RawAncestryEvent(
                _numpy(uid, np.int64),
                _numpy(mother1, np.int64),
                _numpy(mother2, np.int64),
            )
        )
    return output


def ancestor_targets(
    start_index: int,
    mother1: np.ndarray,
    mother2: np.ndarray,
    targets: set[int],
) -> set[int]:
    """Find target indices by Delphes' cycle-safe, zero-based mother walk."""

    found: set[int] = set()
    queue = [int(start_index)]
    visited: set[int] = set()
    while queue:
        current = queue.pop(0)
        if current in visited:
            continue
        visited.add(current)
        if current in targets:
            found.add(current)
            continue
        if current <= 1 or current >= len(mother1):
            continue
        for parent in (int(mother1[current]), int(mother2[current])):
            if 1 < parent < len(mother1) and parent not in visited and parent not in queue:
                queue.append(parent)
    return found


def has_ancestor_target(
    start_index: int,
    mother1: np.ndarray,
    mother2: np.ndarray,
    targets: set[int],
    cache: dict[int, bool] | None = None,
) -> bool:
    """Return whether a particle descends from any target generator index."""

    memo = cache if cache is not None else {}
    visiting: set[int] = set()

    def resolve(current: int) -> bool:
        if current in memo:
            return memo[current]
        if current in targets:
            memo[current] = True
            return True
        if current <= 1 or current >= len(mother1):
            memo[current] = False
            return False
        if current in visiting:
            return False
        visiting.add(current)
        for parent in (int(mother1[current]), int(mother2[current])):
            if 1 < parent < len(mother1) and resolve(parent):
                visiting.remove(current)
                memo[current] = True
                return True
        visiting.remove(current)
        memo[current] = False
        return False

    return resolve(int(start_index))


def build_fixed_truth_event(
    event: Any, raw_ancestry: RawAncestryEvent
) -> FixedTruthEvent:
    """Reproduce the algo-20 ancestry assignment from raw generator records.

    The dark-hadron jet partition is the stored R=0.8 partition.  It remains
    fixed throughout the primary radius scan, so truth identifiers cannot drift
    with the scan radius.
    """

    gen_candidates = particle_table(event.GenCandidate)
    dark_hadrons = particle_table(event.DarkHadronCandidate)
    dh_index = {int(uid): idx for idx, uid in enumerate(dark_hadrons.uid)}
    gen_index = {int(uid): idx for idx, uid in enumerate(raw_ancestry.uid)}
    jet_refs = ak.to_list(event.DarkHadronJet.Constituents.refs)

    dh_owner: dict[int, tuple[int, int]] = {}
    ordered_dh_uids: list[tuple[int, ...]] = []
    for truth_id, refs in enumerate(jet_refs):
        uids = tuple(int(uid) for uid in refs if int(uid) in dh_index)
        ordered_dh_uids.append(uids)
        for position, uid in enumerate(uids):
            if uid in gen_index:
                dh_owner[gen_index[uid]] = (truth_id, position)

    targets = set(dh_owner)
    all_dark_hadron_targets = {
        gen_index[uid] for uid in dh_index if uid in gen_index
    }
    per_jet: list[list[list[int]]] = [
        [[] for _ in uids] for uids in ordered_dh_uids
    ]
    per_jet_visible: list[list[int]] = [[] for _ in ordered_dh_uids]
    visible_owner = np.full(len(gen_candidates), -1, dtype=np.int64)
    visible_dh = np.full(len(gen_candidates), -1, dtype=np.int64)
    visible_has_any_dark_hadron_ancestor = np.zeros(
        len(gen_candidates), dtype=bool
    )
    unresolved: list[int] = []
    cache_reassigned: list[int] = []

    # Mirror the installed C++ algo-20 traversal, including its event-global
    # cache and first-match stopping rule, for bitwise stored-branch closure.
    # That walk stays as-is: its cache is already event-global, so it is cheap,
    # and its ordering semantics are what the closure gate checks.
    ancestry_cache: dict[int, int] = {}

    # The two physically-defined ancestry questions -- "which dark hadrons is
    # this particle descended from" and "does it descend from any dark hadron"
    # -- are resolved for the WHOLE event up front instead of once per
    # candidate.  core.ancestor_targets ran an uncached BFS per GenCandidate
    # (~600/event over a ~1700-node graph); this is one fixed-point iteration.
    # Measured 8.32 -> 1.01 ms/event on the CMS validation sample, with
    # identical results on 120,004 candidates.  See ancestry.py.
    _target_list = np.array(sorted(targets), dtype=np.int64)
    _target_mask = ancestor_masks(
        raw_ancestry.mother1, raw_ancestry.mother2, _target_list
    )
    _any_dark_hadron = descends_from_any(
        raw_ancestry.mother1,
        raw_ancestry.mother2,
        np.array(sorted(all_dark_hadron_targets), dtype=np.int64),
    )

    for visible_index, uid in enumerate(gen_candidates.uid):
        if int(uid) not in gen_index:
            continue
        start = gen_index[int(uid)]
        visible_has_any_dark_hadron_ancestor[visible_index] = bool(
            _any_dark_hadron[start]
        )
        matched_jet = -1
        found_match = False
        visited: list[int] = []
        queue = [start]
        in_queue = {start}
        while queue and not found_match:
            current = queue.pop(0)
            visited.append(current)
            if current in ancestry_cache:
                matched_jet = ancestry_cache[current]
                found_match = True
                break
            if current in dh_owner:
                matched_jet = dh_owner[current][0]
                found_match = True
                break
            if current <= 1 or current >= len(raw_ancestry.mother1):
                continue
            for parent in (
                int(raw_ancestry.mother1[current]),
                int(raw_ancestry.mother2[current]),
            ):
                if (
                    1 < parent < len(raw_ancestry.mother1)
                    and parent not in in_queue
                ):
                    queue.append(parent)
                    in_queue.add(parent)
        for index in visited:
            ancestry_cache[index] = matched_jet
        for index in queue:
            ancestry_cache[index] = matched_jet

        if matched_jet >= 0:
            visible_owner[visible_index] = matched_jet
            per_jet_visible[matched_jet].append(visible_index)

        # Independently resolve the individual dark hadron.  This is used only
        # when the physical ancestry and cached jet owner agree; otherwise the
        # per-DH observable would silently inherit the cache-order artifact.
        matches = targets_for_particle(_target_mask, start, _target_list)
        owners = {dh_owner[target] for target in matches}
        if len(owners) == 1:
            physical_truth_id, dh_position = next(iter(owners))
            if physical_truth_id != matched_jet:
                cache_reassigned.append(visible_index)
                continue
            visible_dh[visible_index] = dh_position
            per_jet[matched_jet][dh_position].append(visible_index)
        elif len(owners) > 1:
            unresolved.append(visible_index)

    jets: list[FixedTruthJet] = []
    for truth_id, uids in enumerate(ordered_dh_uids):
        groups = tuple(tuple(sorted(group)) for group in per_jet[truth_id])
        visible = tuple(sorted(per_jet_visible[truth_id]))
        dh_pts = tuple(float(dark_hadrons.pt[dh_index[uid]]) for uid in uids)
        stored = event.DarkHadronJet[truth_id]
        jets.append(
            FixedTruthJet(
                truth_id=truth_id,
                dark_hadron_uids=uids,
                dark_hadron_pts=dh_pts,
                visible_indices=visible,
                visible_by_dark_hadron=groups,
                pt=float(stored.PT),
                eta=float(stored.Eta),
                phi=float(stored.Phi),
                mass=float(stored.Mass),
            )
        )
    return FixedTruthEvent(
        jets,
        visible_owner,
        visible_dh,
        visible_has_any_dark_hadron_ancestor,
        tuple(unresolved),
        tuple(cache_reassigned),
    )


def visible_storage_sets(event: Any, candidates: ParticleTable) -> list[set[int]]:
    """Resolve stored algo-20 visible-jet TRefs to GenCandidate indices."""

    if "DarkHadronVisibleJet" not in ak.fields(event):
        return []
    uid_index = {int(uid): idx for idx, uid in enumerate(candidates.uid)}
    return [
        {uid_index[int(uid)] for uid in refs if int(uid) in uid_index}
        for refs in ak.to_list(event.DarkHadronVisibleJet.Constituents.refs)
    ]


def match_by_shared_visible_pt(
    truth: FixedTruthEvent,
    clustered: Sequence[ClusteredJet],
    visible_pt: np.ndarray,
) -> list[Match]:
    """One-to-one truth association maximizing shared visible-descendant pT."""

    if not truth.jets or not clustered:
        return []
    overlap = np.zeros((len(truth.jets), len(clustered)), dtype=float)
    clustered_sets = [set(jet.constituent_indices) for jet in clustered]
    for truth_id, truth_jet in enumerate(truth.jets):
        truth_set = set(truth_jet.visible_indices)
        for clustered_id, clustered_set in enumerate(clustered_sets):
            shared = list(truth_set & clustered_set)
            overlap[truth_id, clustered_id] = float(np.sum(visible_pt[shared]))
    truth_rows, clustered_cols = linear_sum_assignment(-overlap)
    matches: list[Match] = []
    for truth_id, clustered_id in zip(truth_rows, clustered_cols):
        shared_pt = float(overlap[truth_id, clustered_id])
        if shared_pt <= 0.0:
            continue
        truth_indices = list(truth.jets[truth_id].visible_indices)
        truth_pt = float(np.sum(visible_pt[truth_indices]))
        matches.append(
            Match(
                truth_id=int(truth_id),
                clustered_index=int(clustered_id),
                shared_visible_pt=shared_pt,
                truth_visible_pt=truth_pt,
                clustered_pt=float(clustered[clustered_id].pt),
            )
        )
    return sorted(matches, key=lambda match: match.truth_id)


def dark_hadron_group_metrics(
    truth: FixedTruthEvent,
    target_truth_id: int,
    clustered_jet: ClusteredJet,
    candidates: ParticleTable,
) -> dict[str, Any]:
    """Exact per-dark-hadron containment and cross-truth contamination.

    A dark hadron is fully contained only when every one of its resolved
    visible final-state descendants is an exact constituent of the clustered
    jet. Dark hadrons with no resolved visible descendants are reported but
    are not classified as fully contained. A contaminating dark hadron belongs
    to a different fixed DarkHadronJet truth partition and contributes at
    least one visible descendant to the target's matched GenFatJet.
    """

    target = truth.jets[target_truth_id]
    clustered_indices = set(clustered_jet.constituent_indices)

    def group_record(
        truth_jet: FixedTruthJet, position: int, group: Sequence[int]
    ) -> dict[str, Any]:
        group_set = set(group)
        captured = group_set & clustered_indices
        total_pt = (
            float(np.sum(candidates.pt[list(group_set)])) if group_set else 0.0
        )
        captured_pt = (
            float(np.sum(candidates.pt[list(captured)])) if captured else 0.0
        )
        return {
            "truth_id": int(truth_jet.truth_id),
            "dark_hadron_position": int(position),
            "dark_hadron_uid": int(truth_jet.dark_hadron_uids[position]),
            "dark_hadron_pt_GeV": float(truth_jet.dark_hadron_pts[position]),
            "n_visible_descendants": len(group_set),
            "n_captured_visible_descendants": len(captured),
            "visible_descendant_pt_fraction": (
                captured_pt / total_pt if total_pt > 0.0 else math.nan
            ),
            "captured_pt_GeV": captured_pt,
            "fully_contained": bool(group_set and group_set <= clustered_indices),
        }

    target_records = [
        group_record(target, position, group)
        for position, group in enumerate(target.visible_by_dark_hadron)
    ]
    contaminating_records: list[dict[str, Any]] = []
    for other in truth.jets:
        if other.truth_id == target_truth_id:
            continue
        for position, group in enumerate(other.visible_by_dark_hadron):
            record = group_record(other, position, group)
            if record["n_captured_visible_descendants"] > 0:
                contaminating_records.append(record)

    target_visible = [
        record for record in target_records if record["n_visible_descendants"] > 0
    ]
    n_target_full = sum(record["fully_contained"] for record in target_visible)
    n_target_not_full = len(target_visible) - n_target_full
    other_owner = (
        (truth.visible_owner >= 0)
        & (truth.visible_owner != target_truth_id)
    )
    other_indices = {
        index for index in clustered_indices if bool(other_owner[index])
    }
    cross_truth_pt = (
        float(np.sum(candidates.pt[list(other_indices)])) if other_indices else 0.0
    )
    non_dark_indices = {
        index
        for index in clustered_indices
        if not truth.visible_has_any_dark_hadron_ancestor[index]
    }
    non_dark_pt = (
        float(np.sum(candidates.pt[list(non_dark_indices)]))
        if non_dark_indices
        else 0.0
    )
    jet_constituent_pt = (
        float(np.sum(candidates.pt[list(clustered_indices)]))
        if clustered_indices
        else 0.0
    )

    # pT the jet picked up from dark hadrons it did NOT capture whole, from
    # ANY truth partition.  This is the quantity the radius study actually
    # needs, and it is not the same as cross-truth contamination: a dark hadron
    # of *this* truth jet that is only half captured is just as damaging as a
    # neighbour's leaking in, while a neighbour's dark hadron that arrives
    # fully contained still carries its own correct mass.
    #
    # A partially captured dark hadron is the worst of the three cases.  Fully
    # contained, its decay products sum to its mass.  Entirely absent, it costs
    # acceptance but tells no lies.  Half in, the jet gains a random fraction
    # of its momentum with no mass meaning attached -- heavier and less correct
    # at the same time.
    partial_records = [
        record
        for record in target_records + contaminating_records
        if record["n_captured_visible_descendants"] > 0
        and not record["fully_contained"]
    ]
    partial_pt = float(sum(record["captured_pt_GeV"] for record in partial_records))
    full_records = [
        record
        for record in target_records + contaminating_records
        if record["n_captured_visible_descendants"] > 0 and record["fully_contained"]
    ]
    full_pt = float(sum(record["captured_pt_GeV"] for record in full_records))

    return {
        "target_dark_hadrons": target_records,
        "contaminating_dark_hadrons": contaminating_records,
        "n_target_dark_hadrons": len(target_records),
        "n_target_dark_hadrons_with_visible_descendants": len(target_visible),
        "n_fully_contained_target_dark_hadrons": int(n_target_full),
        "n_not_fully_contained_target_dark_hadrons": int(n_target_not_full),
        "fraction_fully_contained_target_dark_hadrons": (
            n_target_full / len(target_visible) if target_visible else math.nan
        ),
        "all_target_dark_hadrons_fully_contained": bool(
            target_visible and n_target_full == len(target_visible)
        ),
        "n_contaminating_dark_hadrons": len(contaminating_records),
        "n_fully_contained_contaminating_dark_hadrons": int(
            sum(record["fully_contained"] for record in contaminating_records)
        ),
        "cross_truth_dark_hadron_contamination_pt": (
            cross_truth_pt / jet_constituent_pt
            if jet_constituent_pt > 0.0
            else math.nan
        ),
        "n_non_dark_hadron_constituents": len(non_dark_indices),
        "fraction_non_dark_hadron_constituents": (
            len(non_dark_indices) / len(clustered_indices)
            if clustered_indices
            else math.nan
        ),
        "non_dark_hadron_constituent_pt_fraction": (
            non_dark_pt / jet_constituent_pt
            if jet_constituent_pt > 0.0
            else math.nan
        ),
        "has_non_dark_hadron_constituents": bool(non_dark_indices),
        "n_partially_contained_dark_hadrons_in_jet": len(partial_records),
        "n_fully_contained_dark_hadrons_in_jet": len(full_records),
        "partially_contained_dark_hadron_pt_fraction": (
            partial_pt / jet_constituent_pt if jet_constituent_pt > 0.0 else math.nan
        ),
        "fully_contained_dark_hadron_pt_fraction": (
            full_pt / jet_constituent_pt if jet_constituent_pt > 0.0 else math.nan
        ),
    }


def matched_jet_metrics(
    truth: FixedTruthEvent,
    truth_jet: FixedTruthJet,
    clustered_jet: ClusteredJet,
    candidates: ParticleTable,
) -> dict[str, float | int]:
    """Exact generator-level containment and contamination for one match."""

    truth_indices = set(truth_jet.visible_indices)
    clustered_indices = set(clustered_jet.constituent_indices)
    captured = truth_indices & clustered_indices
    all_truth_visible = {
        index for jet in truth.jets for index in jet.visible_indices
    }
    contaminating = clustered_indices - all_truth_visible
    truth_pt = float(np.sum(candidates.pt[list(truth_indices)])) if truth_indices else 0.0
    captured_pt = float(np.sum(candidates.pt[list(captured)])) if captured else 0.0
    jet_constituent_pt = (
        float(np.sum(candidates.pt[list(clustered_indices)])) if clustered_indices else 0.0
    )
    contamination_pt = (
        float(np.sum(candidates.pt[list(contaminating)])) if contaminating else 0.0
    )

    dh_fractions: list[float] = []
    for group in truth_jet.visible_by_dark_hadron:
        group_set = set(group)
        denominator = float(np.sum(candidates.pt[list(group_set)])) if group_set else 0.0
        numerator = float(np.sum(candidates.pt[list(group_set & clustered_indices)]))
        dh_fractions.append(numerator / denominator if denominator > 0.0 else math.nan)
    leading = int(np.argmax(truth_jet.dark_hadron_pts)) if truth_jet.dark_hadron_pts else -1
    leading_fraction = dh_fractions[leading] if leading >= 0 else math.nan
    count_fraction = len(captured) / len(truth_indices) if truth_indices else math.nan
    per_dark_hadron = dark_hadron_group_metrics(
        truth, truth_jet.truth_id, clustered_jet, candidates
    )
    return {
        "truth_id": truth_jet.truth_id,
        "containment_pt": captured_pt / truth_pt if truth_pt > 0.0 else math.nan,
        "containment_count": count_fraction,
        "contamination_pt": contamination_pt / jet_constituent_pt
        if jet_constituent_pt > 0.0
        else math.nan,
        "leading_dark_hadron_containment": leading_fraction,
        "n_dark_hadrons": per_dark_hadron["n_target_dark_hadrons"],
        "n_dark_hadrons_with_visible_descendants": per_dark_hadron[
            "n_target_dark_hadrons_with_visible_descendants"
        ],
        "n_fully_contained_dark_hadrons": per_dark_hadron[
            "n_fully_contained_target_dark_hadrons"
        ],
        "n_not_fully_contained_dark_hadrons": per_dark_hadron[
            "n_not_fully_contained_target_dark_hadrons"
        ],
        "fraction_fully_contained_dark_hadrons": per_dark_hadron[
            "fraction_fully_contained_target_dark_hadrons"
        ],
        "all_dark_hadrons_fully_contained": per_dark_hadron[
            "all_target_dark_hadrons_fully_contained"
        ],
        "n_contaminating_dark_hadrons": per_dark_hadron[
            "n_contaminating_dark_hadrons"
        ],
        "n_fully_contained_contaminating_dark_hadrons": per_dark_hadron[
            "n_fully_contained_contaminating_dark_hadrons"
        ],
        "cross_truth_dark_hadron_contamination_pt": per_dark_hadron[
            "cross_truth_dark_hadron_contamination_pt"
        ],
        "n_partially_contained_dark_hadrons_in_jet": per_dark_hadron[
            "n_partially_contained_dark_hadrons_in_jet"
        ],
        "partially_contained_dark_hadron_pt_fraction": per_dark_hadron[
            "partially_contained_dark_hadron_pt_fraction"
        ],
        "fully_contained_dark_hadron_pt_fraction": per_dark_hadron[
            "fully_contained_dark_hadron_pt_fraction"
        ],
        "n_non_dark_hadron_constituents": per_dark_hadron[
            "n_non_dark_hadron_constituents"
        ],
        "fraction_non_dark_hadron_constituents": per_dark_hadron[
            "fraction_non_dark_hadron_constituents"
        ],
        "non_dark_hadron_constituent_pt_fraction": per_dark_hadron[
            "non_dark_hadron_constituent_pt_fraction"
        ],
        "has_non_dark_hadron_constituents": per_dark_hadron[
            "has_non_dark_hadron_constituents"
        ],
        "n_visible_truth": len(truth_indices),
        "n_visible_captured": len(captured),
        "jet_pt": clustered_jet.pt,
        "jet_mass": clustered_jet.mass,
    }


def radial_percentiles(
    particles: ParticleTable,
    indices: Sequence[int],
    axis_eta: float,
    axis_phi: float,
    percentiles: Sequence[float] = DEFAULT_PERCENTILES,
) -> dict[str, float]:
    """pT-weighted radii containing the requested fractions."""

    idx = np.asarray(indices, dtype=np.int64)
    if len(idx) == 0:
        return {f"r{int(p)}": math.nan for p in percentiles}
    radius = np.sqrt(delta_r2(particles.eta[idx], particles.phi[idx], axis_eta, axis_phi))
    order = np.argsort(radius)
    sorted_radius = radius[order]
    sorted_pt = particles.pt[idx][order]
    total = float(np.sum(sorted_pt))
    if total <= 0.0:
        return {f"r{int(p)}": math.nan for p in percentiles}
    cumulative = np.cumsum(sorted_pt) / total
    output: dict[str, float] = {}
    for percentile in percentiles:
        target = float(percentile) / 100.0
        position = min(int(np.searchsorted(cumulative, target, side="left")), len(idx) - 1)
        output[f"r{int(percentile)}"] = float(sorted_radius[position])
    return output


def shape_observables(
    particles: ParticleTable,
    indices: Sequence[int],
    axis_eta: float,
    axis_phi: float,
) -> dict[str, float]:
    """Girth, pTD, and pT-weighted major/minor axes with explicit axes."""

    idx = np.asarray(indices, dtype=np.int64)
    if len(idx) == 0:
        return {"girth": math.nan, "ptd": math.nan, "major": math.nan, "minor": math.nan}
    weights = particles.pt[idx]
    sum_pt = float(np.sum(weights))
    if sum_pt <= 0.0:
        return {"girth": math.nan, "ptd": math.nan, "major": math.nan, "minor": math.nan}
    deta = particles.eta[idx] - axis_eta
    dphi = wrap_delta_phi(particles.phi[idx], axis_phi)
    radius = np.hypot(deta, dphi)
    covariance = np.array(
        [
            [np.sum(weights * deta * deta), np.sum(weights * deta * dphi)],
            [np.sum(weights * deta * dphi), np.sum(weights * dphi * dphi)],
        ]
    ) / sum_pt
    eigenvalues = np.maximum(np.linalg.eigvalsh(covariance), 0.0)
    return {
        "girth": float(np.sum(weights * radius) / sum_pt),
        "ptd": float(np.sqrt(np.sum(weights * weights)) / sum_pt),
        "major": float(np.sqrt(eigenvalues[-1])),
        "minor": float(np.sqrt(eigenvalues[0])),
    }


def fixed_axis_containment(
    particles: ParticleTable,
    truth_indices: Sequence[int],
    axis_eta: float,
    axis_phi: float,
    radii: Sequence[float],
) -> np.ndarray:
    idx = np.asarray(truth_indices, dtype=np.int64)
    denominator = float(np.sum(particles.pt[idx])) if len(idx) else 0.0
    values: list[float] = []
    for radius in radii:
        mask = delta_r2(
            particles.eta[idx], particles.phi[idx], axis_eta, axis_phi
        ) < float(radius) ** 2
        numerator = float(np.sum(particles.pt[idx][mask]))
        values.append(numerator / denominator if denominator > 0.0 else math.nan)
    return np.asarray(values)


def bootstrap_mean(
    values: Sequence[float],
    n_resamples: int = 1000,
    seed: int = 12345,
) -> dict[str, float]:
    clean = np.asarray(values, dtype=float)
    clean = clean[np.isfinite(clean)]
    if len(clean) == 0:
        return {"mean": math.nan, "lower": math.nan, "upper": math.nan, "n": 0}
    rng = np.random.default_rng(seed)
    samples = rng.choice(clean, size=(int(n_resamples), len(clean)), replace=True)
    means = np.mean(samples, axis=1)
    return {
        "mean": float(np.mean(clean)),
        "lower": float(np.percentile(means, 2.5)),
        "upper": float(np.percentile(means, 97.5)),
        "n": int(len(clean)),
    }


_MODULE_START = re.compile(r"^\s*module\s+\S+\s+(\S+)\s*\{")
_PARAMETER_R = re.compile(
    r"^(?P<indent>\s*)set\s+ParameterR\s+(?P<value>[-+0-9.eE]+)(?P<tail>\s*(?:#.*)?)$"
)


def module_spans(card_text: str) -> dict[str, tuple[int, int]]:
    """Return line spans for top-level Delphes Tcl modules."""

    lines = card_text.splitlines(keepends=True)
    spans: dict[str, tuple[int, int]] = {}
    active: str | None = None
    start = -1
    depth = 0
    for index, line in enumerate(lines):
        code = line.split("#", 1)[0]
        if active is None:
            match = _MODULE_START.match(code)
            if match:
                active = match.group(1)
                start = index
                depth = code.count("{") - code.count("}")
                if depth == 0:
                    spans[active] = (start, index + 1)
                    active = None
        else:
            depth += code.count("{") - code.count("}")
            if depth == 0:
                spans[active] = (start, index + 1)
                active = None
    if active is not None:
        raise ValueError(f"Unclosed module {active!r}")
    return spans


def module_parameter_r(card_text: str) -> dict[str, float]:
    lines = card_text.splitlines()
    output: dict[str, float] = {}
    for module, (start, stop) in module_spans(card_text).items():
        for line in lines[start:stop]:
            match = _PARAMETER_R.match(line)
            if match:
                output[module] = float(match.group("value"))
    return output


def override_module_parameter_r(card_text: str, module: str, radius: float) -> str:
    """Change exactly one module-local ParameterR and nothing else."""

    lines = card_text.splitlines(keepends=True)
    spans = module_spans(card_text)
    if module not in spans:
        raise KeyError(f"Module {module!r} is absent from the card")
    start, stop = spans[module]
    changed = 0
    for index in range(start, stop):
        newline = "\n" if lines[index].endswith("\n") else ""
        body = lines[index][:-1] if newline else lines[index]
        match = _PARAMETER_R.match(body)
        if match:
            lines[index] = (
                f"{match.group('indent')}set ParameterR {float(radius):.12g}"
                f"{match.group('tail')}{newline}"
            )
            changed += 1
    if changed != 1:
        raise ValueError(
            f"Expected one ParameterR in module {module!r}, found {changed}"
        )
    return "".join(lines)


def card_roundtrip_check(card_text: str, module: str, radius: float) -> dict[str, Any]:
    before = module_parameter_r(card_text)
    changed_text = override_module_parameter_r(card_text, module, radius)
    after = module_parameter_r(changed_text)
    differences = {
        key: (before.get(key), after.get(key))
        for key in sorted(set(before) | set(after))
        if before.get(key) != after.get(key)
    }
    expected = {module: (before.get(module), float(radius))}
    return {
        "passed": differences == expected,
        "module": module,
        "requested_radius": float(radius),
        "differences": differences,
        "text": changed_text,
    }


def git_sha(path: str | Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def study_metadata(
    input_file: str | Path,
    repository: str | Path,
    delphes_path: str | Path | None = None,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    packages = {}
    for package in ("awkward", "coffea", "fastjet", "numpy", "scipy", "uproot"):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = None
    metadata: dict[str, Any] = {
        "purpose": "jet-radius physics study (not training-data production)",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_file": str(Path(input_file).resolve()),
        "repository_sha": git_sha(repository),
        "delphes_sha": git_sha(delphes_path) if delphes_path else None,
        "pythia8_minor": os.environ.get("PYTHIA8MINOR"),
        "packages": packages,
        "radius_grid": list(DEFAULT_R_GRID),
        "truth_policy": (
            "fixed shipped R=0.8 DarkHadronJet partition; one-to-one maximum "
            "shared visible-descendant pT association"
        ),
    }
    if extra:
        metadata.update(extra)
    return metadata


def json_ready(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if hasattr(value, "__dataclass_fields__"):
        return {key: json_ready(item) for key, item in asdict(value).items()}
    if isinstance(value, Mapping):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path: str | Path, payload: Any) -> None:
    Path(path).write_text(
        json.dumps(json_ready(payload), indent=2, sort_keys=True) + "\n"
    )
