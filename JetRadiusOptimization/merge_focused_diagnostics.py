#!/usr/bin/env python3
"""Merge independent radius jobs and make explanatory overlay plots."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


RADIUS_COLORS = {
    0.2: "#6B7280",
    0.4: "#7B2CBF",
    0.6: "#2563EB",
    0.8: "#111111",
    1.0: "#E76F51",
    1.2: "#D81B60",
    1.4: "#0086A8",
    1.6: "#C62828",
}
DEFAULT_DIAGNOSTIC_RADII = (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6)
ARRAY_METRICS = (
    "n_genfatjets_per_event",
    "n_fully_contained_target_dark_hadrons",
    "n_not_fully_contained_target_dark_hadrons",
    "n_contaminating_dark_hadrons",
    "n_fully_contained_contaminating_dark_hadrons",
    "n_non_dark_hadron_constituents",
    "fraction_non_dark_hadron_constituents",
    "non_dark_hadron_constituent_pt_fraction",
    "has_non_dark_hadron_constituents",
    "softdrop_mass_GeV",
)
COUNT_METRICS = (
    "eligible_matched_truth_jets",
    "eligible_matched_multi_dark_hadron_truth_jets",
    "pass_at_least_one_fully_contained_dark_hadron",
    "multi_dark_hadron_pass",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument(
        "--diagnostic-radii",
        nargs="+",
        type=float,
        default=list(DEFAULT_DIAGNOSTIC_RADII),
    )
    return parser.parse_args()


def radius_color(radius: float) -> str:
    for configured, color in RADIUS_COLORS.items():
        if math.isclose(radius, configured, abs_tol=1.0e-9):
            return color
    fallback = ("#5B21B6", "#1D4ED8", "#B45309", "#BE123C", "#0E7490")
    return fallback[int(round(radius * 10.0)) % len(fallback)]


def save_figure(figure: Any, path: Path) -> str:
    figure.tight_layout()
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)
    return str(path.resolve())


def probability_histogram(
    axis: Any,
    values_by_radius: dict[float, np.ndarray],
    bins: np.ndarray,
    xlabel: str,
    title: str,
) -> None:
    for radius, raw in sorted(values_by_radius.items()):
        values = np.asarray(raw, dtype=float)
        values = values[np.isfinite(values)]
        if not len(values):
            continue
        axis.hist(
            values,
            bins=bins,
            weights=np.full(len(values), 1.0 / len(values)),
            histtype="step",
            linewidth=2.0,
            color=radius_color(radius),
            label=f"R={radius:.1f} (median={np.median(values):.2f})",
        )
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Fraction of matched truth jets / bin")
    axis.set_yscale("log", nonpositive="clip")
    axis.set_ylim(bottom=1.0e-6)
    axis.set_title(title)
    axis.grid(alpha=0.22, which="both")
    axis.legend(fontsize=8, frameon=False, ncol=2)


def scalar(payload: Any, key: str) -> int:
    return int(np.asarray(payload[key]).reshape(-1)[0])


def merge_payloads(paths: Sequence[Path]) -> dict[float, dict[str, Any]]:
    merged: dict[float, dict[str, Any]] = {}
    coverage: dict[float, list[tuple[int, int, Path]]] = {}
    for path in paths:
        with np.load(path) as payload:
            radii = np.asarray(payload["radii"], dtype=float)
            start = scalar(payload, "entry_start")
            stop = scalar(payload, "entry_stop")
            for radius_index, radius_value in enumerate(radii):
                radius = float(radius_value)
                entry = merged.setdefault(
                    radius,
                    {
                        "arrays": {name: [] for name in ARRAY_METRICS},
                        "counts": {name: 0 for name in COUNT_METRICS},
                        "eligible_truth_jets": 0,
                        "eligible_multi_dark_hadron_truth_jets": 0,
                        "events": 0,
                        "input_files": [],
                    },
                )
                coverage.setdefault(radius, []).append((start, stop, path))
                entry["eligible_truth_jets"] += scalar(
                    payload, "eligible_truth_jets"
                )
                entry["eligible_multi_dark_hadron_truth_jets"] += scalar(
                    payload, "eligible_multi_dark_hadron_truth_jets"
                )
                entry["events"] += scalar(payload, "events_read")
                entry["input_files"].append(str(path.resolve()))
                prefix = f"r{radius_index}_"
                for name in ARRAY_METRICS:
                    entry["arrays"][name].append(
                        np.asarray(payload[prefix + name], dtype=float)
                    )
                for name in COUNT_METRICS:
                    entry["counts"][name] += scalar(payload, prefix + name)

    for radius, ranges in coverage.items():
        ordered = sorted(ranges)
        for previous, current in zip(ordered, ordered[1:]):
            if current[0] < previous[1]:
                raise RuntimeError(
                    f"overlapping event ranges for R={radius}: {previous}/{current}"
                )
        for name in ARRAY_METRICS:
            merged[radius]["arrays"][name] = np.concatenate(
                merged[radius]["arrays"][name]
            )
    coverage_signatures = {
        radius: [(start, stop) for start, stop, _ in sorted(ranges)]
        for radius, ranges in coverage.items()
    }
    expected_coverage = next(iter(coverage_signatures.values()))
    inconsistent = {
        radius: signature
        for radius, signature in coverage_signatures.items()
        if signature != expected_coverage
    }
    if inconsistent:
        raise RuntimeError(
            "radius jobs do not cover identical event ranges: "
            f"expected {expected_coverage}, found {inconsistent}"
        )
    denominators = {
        (
            entry["eligible_truth_jets"],
            entry["eligible_multi_dark_hadron_truth_jets"],
        )
        for entry in merged.values()
    }
    if len(denominators) != 1:
        raise RuntimeError(
            f"fixed-truth eligibility denominators differ by radius: {denominators}"
        )
    return dict(sorted(merged.items()))


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    merged = merge_payloads(args.inputs)
    radii = np.asarray(sorted(merged), dtype=float)
    requested = set(float(radius) for radius in args.diagnostic_radii)
    diagnostic_radii = [radius for radius in radii if float(radius) in requested]
    if not diagnostic_radii:
        raise RuntimeError("none of the requested diagnostic radii were produced")

    def arrays(name: str, subset: Sequence[float] = radii) -> dict[float, np.ndarray]:
        return {
            float(radius): merged[float(radius)]["arrays"][name]
            for radius in subset
        }

    acceptance = []
    no_full = []
    multi_acceptance = []
    multi_no_full = []
    conditional_efficiency = []
    has_not_full = []
    has_non_dark = []
    mean_non_dark_count_fraction = []
    mean_non_dark_pt_fraction = []
    mean_non_dark_count = []
    radius_summary: dict[str, Any] = {}
    for radius_value in radii:
        radius = float(radius_value)
        entry = merged[radius]
        denominator = entry["eligible_truth_jets"]
        multi_denominator = entry["eligible_multi_dark_hadron_truth_jets"]
        matched = entry["counts"]["eligible_matched_truth_jets"]
        passing = entry["counts"]["pass_at_least_one_fully_contained_dark_hadron"]
        multi_passing = entry["counts"]["multi_dark_hadron_pass"]
        acc = passing / denominator if denominator else math.nan
        multi_acc = multi_passing / multi_denominator if multi_denominator else math.nan
        matched_arrays = entry["arrays"]
        not_full_values = matched_arrays[
            "n_not_fully_contained_target_dark_hadrons"
        ]
        non_dark_values = matched_arrays["n_non_dark_hadron_constituents"]
        count_fraction_values = matched_arrays[
            "fraction_non_dark_hadron_constituents"
        ]
        pt_fraction_values = matched_arrays[
            "non_dark_hadron_constituent_pt_fraction"
        ]
        acceptance.append(acc)
        no_full.append(1.0 - acc)
        multi_acceptance.append(multi_acc)
        multi_no_full.append(1.0 - multi_acc)
        conditional_efficiency.append(passing / matched if matched else math.nan)
        has_not_full.append(float(np.mean(not_full_values >= 1)))
        has_non_dark.append(float(np.mean(non_dark_values >= 1)))
        mean_non_dark_count_fraction.append(float(np.nanmean(count_fraction_values)))
        mean_non_dark_pt_fraction.append(float(np.nanmean(pt_fraction_values)))
        mean_non_dark_count.append(float(np.nanmean(non_dark_values)))
        radius_summary[f"{radius:.2f}"] = {
            "events": entry["events"],
            "mean_genfatjets_per_event": float(
                np.mean(matched_arrays["n_genfatjets_per_event"])
            ),
            "median_genfatjets_per_event": float(
                np.median(matched_arrays["n_genfatjets_per_event"])
            ),
            "eligible_truth_jets": denominator,
            "eligible_matched_truth_jets": matched,
            "pass_at_least_one_fully_contained_dark_hadron": passing,
            "acceptance": acc,
            "conditional_efficiency_given_match": passing / matched
            if matched
            else math.nan,
            "eligible_multi_dark_hadron_truth_jets": multi_denominator,
            "multi_dark_hadron_pass": multi_passing,
            "multi_dark_hadron_acceptance": multi_acc,
            "fraction_matched_with_not_fully_contained_target_dark_hadron": has_not_full[-1],
            "fraction_matched_with_non_dark_hadron_constituent": has_non_dark[-1],
            "mean_non_dark_hadron_constituent_count": mean_non_dark_count[-1],
            "mean_non_dark_hadron_constituent_count_fraction": mean_non_dark_count_fraction[-1],
            "mean_non_dark_hadron_constituent_pt_fraction": mean_non_dark_pt_fraction[-1],
        }

    outputs: dict[str, str] = {}
    figure, axes = plt.subplots(2, 1, figsize=(11.0, 10.2), sharex=True)
    axes[0].plot(radii, acceptance, marker="o", color="#2563EB", linewidth=2.2,
                 label="acceptance: ≥1 fully contained DH")
    axes[0].plot(radii, no_full, marker="s", color="#D81B60", linewidth=2.0,
                 label="fails acceptance: no fully contained DH")
    axes[0].plot(radii, multi_acceptance, marker="^", color="#7B2CBF",
                 linestyle="--", linewidth=2.0, label="acceptance: >1 visible DH")
    axes[0].plot(radii, multi_no_full, marker="D", color="#E76F51",
                 linestyle="--", linewidth=1.8, label="fails acceptance: >1 visible DH")
    axes[0].set_ylabel("Fraction of eligible truth jets")
    axes[0].set_title("Dark-hadron containment acceptance")
    axes[0].legend(frameon=False, fontsize=8.8, ncol=2)
    axes[0].text(
        0.02, 0.04, "Acceptance denominator includes unmatched eligible truth jets",
        transform=axes[0].transAxes, fontsize=9.2,
    )
    axes[1].plot(radii, has_not_full, marker="o", color="#7B2CBF",
                 label="≥1 target DH not fully contained")
    axes[1].plot(radii, has_non_dark, marker="s", color="#C62828",
                 label="≥1 constituent not from any initial dark hadron")
    axes[1].plot(radii, mean_non_dark_count_fraction, marker="^", color="#0086A8",
                 label="mean non-DH constituent-count fraction")
    axes[1].plot(radii, mean_non_dark_pt_fraction, marker="D", color="#F59E0B",
                 label=r"mean non-DH constituent-$p_T$ fraction")
    count_axis = axes[1].twinx()
    count_axis.plot(radii, mean_non_dark_count, marker="P", color="#92400E",
                    linestyle="--", label="mean number of non-DH particles")
    count_axis.set_ylabel("Mean number of non-DH particles", color="#92400E")
    count_axis.tick_params(axis="y", colors="#92400E")
    axes[1].set_xlabel("GenFatJet radius R")
    axes[1].set_ylabel("Fraction among matched eligible jets")
    axes[1].set_title("Not-contained target DH and non-DH particle contamination")
    lines, labels = axes[1].get_legend_handles_labels()
    count_lines, count_labels = count_axis.get_legend_handles_labels()
    axes[1].legend(lines + count_lines, labels + count_labels,
                   frameon=False, fontsize=8.4, ncol=2)
    for axis in axes:
        axis.axvline(0.8, color="#111111", linestyle=":", linewidth=1.8)
        axis.set_ylim(0.0, 1.02)
        axis.grid(alpha=0.22)
    outputs["acceptance_and_non_dark_overlay"] = save_figure(
        figure, args.outdir / "any_fully_contained_dark_hadron_acceptance_by_radius.png"
    )

    jet_multiplicity = arrays("n_genfatjets_per_event", diagnostic_radii)
    all_jet_multiplicity = np.concatenate(list(jet_multiplicity.values()))
    jet_count_upper = max(
        5,
        int(math.ceil(float(np.percentile(all_jet_multiplicity, 99.9))))
    )
    clipped_jet_multiplicity = {
        radius: np.minimum(values, jet_count_upper)
        for radius, values in jet_multiplicity.items()
    }
    figure, axes = plt.subplots(1, 2, figsize=(15, 6.4))
    probability_histogram(
        axes[0],
        clipped_jet_multiplicity,
        np.arange(-0.5, jet_count_upper + 1.5, 1.0),
        r"Number of GenFatJets per event ($p_T>15$ GeV)",
        "Event-level GenFatJet multiplicity",
    )
    axes[0].text(
        0.98, 0.04, "Overflow in final bin",
        transform=axes[0].transAxes, ha="right", fontsize=9.2,
    )
    jet_means = []
    jet_medians = []
    jet_low = []
    jet_high = []
    for radius_value in radii:
        values = merged[float(radius_value)]["arrays"]["n_genfatjets_per_event"]
        jet_means.append(float(np.mean(values)))
        jet_medians.append(float(np.median(values)))
        low, high = np.percentile(values, [16.0, 84.0])
        jet_low.append(float(low))
        jet_high.append(float(high))
    axes[1].plot(
        radii, jet_means, marker="o", color="#2563EB", linewidth=2.1,
        label="mean jets / event",
    )
    axes[1].plot(
        radii, jet_medians, marker="s", color="#D81B60", linewidth=1.9,
        label="median jets / event",
    )
    axes[1].fill_between(
        radii, jet_low, jet_high, color="#7B2CBF", alpha=0.16,
        label="event 16–84% interval",
    )
    axes[1].axvline(
        0.8, color="#111111", linestyle=":", linewidth=1.8,
        label="default R=0.8",
    )
    axes[1].set_xlabel("GenFatJet radius R")
    axes[1].set_ylabel(r"Number of GenFatJets per event ($p_T>15$ GeV)")
    axes[1].set_title("How jet multiplicity changes with radius")
    axes[1].grid(alpha=0.22)
    axes[1].legend(frameon=False, fontsize=9)
    outputs["genfatjet_multiplicity"] = save_figure(
        figure, args.outdir / "genfatjet_multiplicity_by_radius.png"
    )

    diagnostic_count = arrays("n_non_dark_hadron_constituents", diagnostic_radii)
    all_counts = np.concatenate(list(diagnostic_count.values()))
    count_upper = max(
        10, int(math.ceil(float(np.percentile(all_counts, 99.5)) / 5.0) * 5)
    )
    clipped_count = {
        radius: np.minimum(values, count_upper)
        for radius, values in diagnostic_count.items()
    }
    figure, axis = plt.subplots(figsize=(10, 7.2))
    probability_histogram(
        axis, clipped_count, np.linspace(-0.5, count_upper + 0.5, 42),
        "Number of constituents not descended from any initial dark hadron",
        "Non-dark-hadron particle multiplicity",
    )
    axis.text(0.98, 0.04, "Overflow in final bin", transform=axis.transAxes,
              ha="right", fontsize=9.2)
    outputs["non_dark_count"] = save_figure(
        figure, args.outdir / "non_dark_hadron_constituent_count_by_radius.png"
    )

    figure, axes = plt.subplots(1, 2, figsize=(15, 6.4), sharey=True)
    probability_histogram(
        axes[0], arrays("fraction_non_dark_hadron_constituents", diagnostic_radii),
        np.linspace(0.0, 1.0, 41),
        "Fraction of constituents not from any initial dark hadron",
        "Particle-count contamination fraction",
    )
    probability_histogram(
        axes[1], arrays("non_dark_hadron_constituent_pt_fraction", diagnostic_radii),
        np.linspace(0.0, 1.0, 41),
        r"Fraction of constituent $p_T$ not from any initial dark hadron",
        r"Non-DH $p_T$ contamination fraction",
    )
    outputs["non_dark_fractions"] = save_figure(
        figure, args.outdir / "non_dark_hadron_constituent_fractions_by_radius.png"
    )

    figure, axes = plt.subplots(1, 3, figsize=(19, 6.2), sharey=True)
    count_specs = (
        ("n_fully_contained_target_dark_hadrons", "Number fully contained", "Target DH fully contained"),
        ("n_not_fully_contained_target_dark_hadrons", "Number not fully contained", "Target DH not fully contained"),
        ("n_contaminating_dark_hadrons", "Number from other truth jets", "Contaminating dark hadrons"),
    )
    for axis, (name, xlabel, title) in zip(axes, count_specs):
        data = arrays(name, diagnostic_radii)
        upper = max(1, int(max(np.max(value) for value in data.values())))
        probability_histogram(
            axis, data, np.arange(-0.5, upper + 1.5, 1.0), xlabel, title
        )
    outputs["dark_hadron_multiplicities"] = save_figure(
        figure, args.outdir / "dark_hadron_multiplicity_diagnostics_by_radius.png"
    )

    reference_radius = min(merged, key=lambda radius: abs(radius - 0.8))
    reference = merged[reference_radius]["arrays"]
    not_full_values = reference["n_not_fully_contained_target_dark_hadrons"]
    non_dark_values = reference["n_non_dark_hadron_constituents"]
    non_dark_pt = reference["non_dark_hadron_constituent_pt_fraction"]
    softdrop_mass = reference["softdrop_mass_GeV"]
    figure, axes = plt.subplots(1, 2, figsize=(15, 6.3))
    first = axes[0].hist2d(
        not_full_values, np.minimum(non_dark_values, count_upper),
        bins=(np.arange(-0.5, max(1, int(np.max(not_full_values))) + 1.5, 1.0), 40),
        norm=LogNorm(), cmap="magma",
    )
    figure.colorbar(first[3], ax=axes[0], label="Matched truth jets / bin")
    axes[0].set_xlabel("Number of target DH not fully contained")
    axes[0].set_ylabel("Number of non-DH constituents (overflow clipped)")
    axes[0].set_title(f"Containment versus contamination, R={reference_radius:.1f}")
    mass_upper = max(100.0, float(np.percentile(softdrop_mass, 99.5)))
    second = axes[1].hist2d(
        non_dark_pt, np.minimum(softdrop_mass, mass_upper),
        bins=(np.linspace(0.0, 1.0, 41), 45), norm=LogNorm(), cmap="plasma",
    )
    figure.colorbar(second[3], ax=axes[1], label="Matched truth jets / bin")
    axes[1].set_xlabel(r"Non-DH constituent $p_T$ fraction")
    axes[1].set_ylabel(r"GenFatJet $m_{SD}$ [GeV] (overflow clipped)")
    axes[1].set_title(f"Soft Drop mass versus contamination, R={reference_radius:.1f}")
    outputs["correlations"] = save_figure(
        figure, args.outdir / "containment_contamination_correlations_R08.png"
    )

    report = {
        "status": "complete",
        "execution": "independent radius jobs merged from exact event-level arrays",
        "radii": radius_summary,
        "diagnostic_radii": [float(radius) for radius in diagnostic_radii],
        "definition": (
            "non-DH constituent = exact GenFatJet constituent with no ancestor "
            "among any initial DarkHadronCandidate in the event"
        ),
        "plots": outputs,
    }
    with (args.outdir / "focused_diagnostics_report.json").open("w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
