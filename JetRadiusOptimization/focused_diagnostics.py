#!/usr/bin/env python3
"""Build mergeable containment/contamination metrics for one event shard.

This is the fast path for explanatory plots. It does not repeat closure tests,
event displays, or the full 0.1-spaced optimization scan performed by
``validate.py``. Independent event ranges can therefore run on separate batch
nodes and be merged without changing any event-level physics definition.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parent
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from coffea.nanoevents import NanoEventsFactory  # noqa: E402
from common import DelphesSchema2  # noqa: E402
from core import (  # noqa: E402
    COLLECTIONS,
    build_fixed_truth_event,
    dark_hadron_group_metrics,
    load_raw_ancestry,
    match_by_shared_visible_pt,
    particle_table,
    recluster_event,
    softdrop_jet,
)


DEFAULT_RADII = (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--entry-start", required=True, type=int)
    parser.add_argument("--max-events", required=True, type=int)
    parser.add_argument("--shard-id", required=True, type=int)
    parser.add_argument("--radii", nargs="+", type=float, default=list(DEFAULT_RADII))
    args = parser.parse_args()
    if args.entry_start < 0 or args.max_events < 1 or args.shard_id < 0:
        parser.error("entry start/shard ID must be nonnegative and max events positive")
    return args


def load_event_range(path: Path, entry_start: int, entry_stop: int) -> Any:
    if not path.is_file():
        raise FileNotFoundError(path)
    return NanoEventsFactory.from_root(
        {str(path.resolve()): "Delphes"},
        schemaclass=DelphesSchema2,
        entry_start=entry_start,
        entry_stop=entry_stop,
    ).events()


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    radii = sorted(set(float(radius) for radius in args.radii))
    entry_stop = args.entry_start + args.max_events
    events = load_event_range(args.input, args.entry_start, entry_stop)
    raw_ancestry = load_raw_ancestry(
        args.input, entry_start=args.entry_start, entry_stop=entry_stop
    )
    if len(events) != len(raw_ancestry):
        raise RuntimeError(
            f"NanoEvents/raw ancestry length mismatch: {len(events)}/{len(raw_ancestry)}"
        )

    metric_names = (
        "n_fully_contained_target_dark_hadrons",
        "n_not_fully_contained_target_dark_hadrons",
        "n_contaminating_dark_hadrons",
        "n_fully_contained_contaminating_dark_hadrons",
        "n_non_dark_hadron_constituents",
        "fraction_non_dark_hadron_constituents",
        "non_dark_hadron_constituent_pt_fraction",
        "has_non_dark_hadron_constituents",
    )
    values: dict[float, dict[str, list[float]]] = {
        radius: defaultdict(list) for radius in radii
    }
    counts: dict[float, dict[str, int]] = {
        radius: defaultdict(int) for radius in radii
    }
    eligible_truth_jets = 0
    eligible_multi_dark_hadron_truth_jets = 0

    for event_index in range(len(events)):
        event = events[event_index]
        candidates = particle_table(event.GenCandidate)
        truth = build_fixed_truth_event(event, raw_ancestry[event_index])
        eligibility: dict[int, tuple[bool, bool]] = {}
        for truth_jet in truth.jets:
            n_visible_groups = sum(
                bool(group) for group in truth_jet.visible_by_dark_hadron
            )
            eligible = n_visible_groups >= 1
            multi = n_visible_groups > 1
            eligibility[truth_jet.truth_id] = (eligible, multi)
            eligible_truth_jets += int(eligible)
            eligible_multi_dark_hadron_truth_jets += int(multi)

        for radius in radii:
            jets = recluster_event(
                candidates, radius, COLLECTIONS["GenFatJet"].pt_min
            )
            values[radius]["n_genfatjets_per_event"].append(float(len(jets)))
            matches = match_by_shared_visible_pt(truth, jets, candidates.pt)
            for match in matches:
                eligible, multi = eligibility[match.truth_id]
                if not eligible:
                    continue
                jet = jets[match.clustered_index]
                metrics = dark_hadron_group_metrics(
                    truth, match.truth_id, jet, candidates
                )
                counts[radius]["eligible_matched_truth_jets"] += 1
                counts[radius]["eligible_matched_multi_dark_hadron_truth_jets"] += int(
                    multi
                )
                passes = metrics["n_fully_contained_target_dark_hadrons"] >= 1
                counts[radius]["pass_at_least_one_fully_contained_dark_hadron"] += int(
                    passes
                )
                counts[radius]["multi_dark_hadron_pass"] += int(passes and multi)
                for name in metric_names:
                    values[radius][name].append(float(metrics[name]))
                values[radius]["softdrop_mass_GeV"].append(
                    softdrop_jet(
                        candidates, jet, beta=0.0, zcut=0.1, r0=radius
                    ).mass
                )

    payload: dict[str, np.ndarray] = {
        "radii": np.asarray(radii, dtype=np.float64),
        "entry_start": np.asarray([args.entry_start], dtype=np.int64),
        "entry_stop": np.asarray([args.entry_start + len(events)], dtype=np.int64),
        "events_read": np.asarray([len(events)], dtype=np.int64),
        "eligible_truth_jets": np.asarray([eligible_truth_jets], dtype=np.int64),
        "eligible_multi_dark_hadron_truth_jets": np.asarray(
            [eligible_multi_dark_hadron_truth_jets], dtype=np.int64
        ),
    }
    for radius_index, radius in enumerate(radii):
        prefix = f"r{radius_index}_"
        for name, entries in values[radius].items():
            payload[prefix + name] = np.asarray(entries, dtype=np.float64)
        for name, count in counts[radius].items():
            payload[prefix + name] = np.asarray([count], dtype=np.int64)

    metrics_path = args.outdir / f"focused_metrics_shard_{args.shard_id:02d}.npz"
    np.savez_compressed(metrics_path, **payload)
    summary = {
        "status": "complete",
        "shard_id": args.shard_id,
        "input": str(args.input.resolve()),
        "entry_start": args.entry_start,
        "entry_stop": args.entry_start + len(events),
        "events_read": len(events),
        "radii": radii,
        "eligible_truth_jets": eligible_truth_jets,
        "eligible_multi_dark_hadron_truth_jets": eligible_multi_dark_hadron_truth_jets,
        "metrics": str(metrics_path.resolve()),
        "definition": (
            "non-DH constituent = exact GenFatJet constituent with no ancestor "
            "among any initial DarkHadronCandidate in the event"
        ),
    }
    with (args.outdir / f"focused_summary_shard_{args.shard_id:02d}.json").open(
        "w"
    ) as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
