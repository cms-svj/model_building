#!/usr/bin/env python3
"""Stage 1: events.root -> radius-independent skim.

Everything in the radius scan that does not depend on R happens once, here:
ROOT I/O, dark-hadron ancestry, and the fixed R=0.8 truth partition.  Stage 2
(``scan.py``) then runs every radius on top of the skim.

Two reasons the split is worth it beyond raw CPU:

  * ``condor/focused_by_radius.dag`` queues one job per radius over the same
    events, so it recomputes the identical radius-independent ancestry once per
    radius.  Sharding over events instead keeps that at 1x.
  * changing a metric definition currently means re-reading ROOT and redoing
    ancestry across the whole campaign.  Against a skim it is seconds.

The skim stores exactly the fields needed to rebuild ``core.ParticleTable``
and ``core.FixedTruthEvent``, so stage 2 calls ``core.matched_jet_metrics``
and ``core.match_by_shared_visible_pt`` unchanged rather than reimplementing
them.  Nothing about the physics definitions lives in this file.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import awkward as ak
import numpy as np

HERE = Path(__file__).resolve().parent
for path in (str(HERE.parent), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

from coffea.nanoevents import NanoEventsFactory  # noqa: E402
from common import DelphesSchema2  # noqa: E402
import core  # noqa: E402

SCHEMA_VERSION = 1

CANDIDATE_FIELDS = (
    "pt", "eta", "phi", "energy", "px", "py", "pz", "uid",
    "owner", "dh_position", "has_dh_ancestor",
)
TRUTH_FIELDS = (
    "truth_id", "pt", "eta", "phi", "mass",
    "dark_hadron_uids", "dark_hadron_pts",
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input", required=True, type=str,
                   help="Local events.root path, or a root://... xrootd URL.")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--entry-start", type=int, default=0)
    p.add_argument("--max-events", type=int, default=None)
    p.add_argument("--label", default="unknown")
    p.add_argument("--point-id", type=int, default=-1)
    p.add_argument("--params", default=None,
                   help="JSON dict of dark-sector parameters to embed.")
    return p.parse_args()


def resolve_input(spec: str) -> str:
    """Absolute local path, or an xrootd/http URL passed through unchanged.

    ``Path(spec).resolve()`` corrupts a URL like ``root://host//eos/...``: it
    collapses the double slash after the host and, since the result no
    longer looks absolute, prepends the current working directory --
    ``root://host//eos/x`` silently becomes
    ``<cwd>/root:/host/eos/x``, which then fails as "file not found" with no
    hint that the URL was ever mangled. First hit running against a real
    remote sample (skim.py had only been exercised against
    local files until then).
    """
    return spec if "://" in spec else str(Path(spec).resolve())


def build_skim(input_spec: str, entry_start: int, max_events: int | None):
    entry_stop = None if max_events is None else entry_start + max_events
    resolved = resolve_input(input_spec)

    events = NanoEventsFactory.from_root(
        {resolved: "Delphes"},
        schemaclass=DelphesSchema2,
        entry_start=entry_start,
        entry_stop=entry_stop,
    ).events()
    raw_ancestry = core.load_raw_ancestry(
        resolved, entry_start=entry_start, entry_stop=entry_stop
    )
    n_events = len(events)
    if n_events != len(raw_ancestry):
        raise RuntimeError(
            f"NanoEvents/raw ancestry length mismatch: {n_events}/{len(raw_ancestry)}"
        )

    cand: dict[str, list] = {k: [] for k in CANDIDATE_FIELDS}
    truth: dict[str, list] = {k: [] for k in TRUTH_FIELDS}
    n_unresolved = n_reassigned = 0

    for i in range(n_events):
        event = events[i]
        table = core.particle_table(event.GenCandidate)
        fixed = core.build_fixed_truth_event(event, raw_ancestry[i])
        n_unresolved += len(fixed.unresolved_visible)
        n_reassigned += len(fixed.cache_reassigned_visible)

        cand["pt"].append(np.asarray(table.pt, dtype=np.float32))
        cand["eta"].append(np.asarray(table.eta, dtype=np.float32))
        cand["phi"].append(np.asarray(table.phi, dtype=np.float32))
        cand["energy"].append(np.asarray(table.energy, dtype=np.float32))
        cand["px"].append(np.asarray(table.px, dtype=np.float32))
        cand["py"].append(np.asarray(table.py, dtype=np.float32))
        cand["pz"].append(np.asarray(table.pz, dtype=np.float32))
        cand["uid"].append(np.asarray(table.uid, dtype=np.int64))
        cand["owner"].append(np.asarray(fixed.visible_owner, dtype=np.int16))
        cand["dh_position"].append(
            np.asarray(fixed.visible_dark_hadron, dtype=np.int16))
        cand["has_dh_ancestor"].append(
            np.asarray(fixed.visible_has_any_dark_hadron_ancestor, dtype=bool))

        truth["truth_id"].append(
            np.array([j.truth_id for j in fixed.jets], dtype=np.int16))
        truth["pt"].append(np.array([j.pt for j in fixed.jets], dtype=np.float32))
        truth["eta"].append(np.array([j.eta for j in fixed.jets], dtype=np.float32))
        truth["phi"].append(np.array([j.phi for j in fixed.jets], dtype=np.float32))
        truth["mass"].append(np.array([j.mass for j in fixed.jets], dtype=np.float32))
        truth["dark_hadron_uids"].append(
            [list(int(u) for u in j.dark_hadron_uids) for j in fixed.jets])
        truth["dark_hadron_pts"].append(
            [list(float(p) for p in j.dark_hadron_pts) for j in fixed.jets])

    skim = ak.Array({
        "candidate": ak.zip({k: ak.Array(v) for k, v in cand.items()}),
        "truth_id": ak.Array(truth["truth_id"]),
        "truth_pt": ak.Array(truth["pt"]),
        "truth_eta": ak.Array(truth["eta"]),
        "truth_phi": ak.Array(truth["phi"]),
        "truth_mass": ak.Array(truth["mass"]),
        "truth_dark_hadron_uids": ak.Array(truth["dark_hadron_uids"]),
        "truth_dark_hadron_pts": ak.Array(truth["dark_hadron_pts"]),
    })
    stats = {
        "events": n_events,
        "candidates": int(sum(len(x) for x in cand["pt"])),
        "truth_jets": int(sum(len(x) for x in truth["pt"])),
        "unresolved_visible": n_unresolved,
        "cache_reassigned_visible": n_reassigned,
    }
    return skim, stats


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    skim, stats = build_skim(args.input, args.entry_start, args.max_events)

    meta = {
        "schema_version": SCHEMA_VERSION,
        "stage": "skim",
        "label": args.label,
        "point_id": args.point_id,
        "input": resolve_input(args.input),
        "entry_start": args.entry_start,
        "entry_stop": args.entry_start + stats["events"],
        "params": json.loads(args.params) if args.params else {},
        "repository_sha": core.git_sha(HERE),
        **stats,
        "note": (
            "Radius-independent. Rebuild core.FixedTruthEvent with "
            "scan.rebuild_truth_event: truth jet t's visible descendants are "
            "candidate.owner == t; dark hadron p of that jet additionally has "
            "candidate.dh_position == p."
        ),
    }
    ak.to_parquet(skim, args.output, extensionarray=False)
    args.output.with_suffix(".meta.json").write_text(
        json.dumps(meta, indent=2, sort_keys=True) + "\n")

    size_mb = args.output.stat().st_size / 1e6
    print(json.dumps(meta, indent=2, sort_keys=True))
    print(f"[OK] {args.output}  {size_mb:.1f} MB "
          f"({1e3*size_mb/max(stats['events'],1):.1f} kB/event)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
