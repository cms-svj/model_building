#!/usr/bin/env python3
"""Stage 2: skim -> per-truth-jet metric table at every radius.

Reads only the Parquet skim from ``skim.py``.  No ROOT, no coffea, no
ancestry -- which is why this is the stage you can re-run whenever a metric
definition changes.

The physics is *not* reimplemented here.  This module rebuilds
``core.ParticleTable`` and ``core.FixedTruthEvent`` from the skim and then
calls ``core.match_by_shared_visible_pt`` and ``core.matched_jet_metrics``
unchanged, so containment and contamination keep exactly one definition in the
repository.  The only additions are the derived scalars in
``radius_metrics.py`` and the zero-fill for unmatched truth jets.

Clustering goes through ``core.recluster_events``, which hands FastJet the
whole chunk at once: 7.81 -> 1.02 ms/event/radius against the per-event
``PseudoJet`` loop, with identical constituent partitions.

Output is one flat row per (event, truth_id, radius), so the cross-model step
is a groupby rather than a bespoke merge.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import awkward as ak
import numpy as np

HERE = Path(__file__).resolve().parent
for path in (str(HERE.parent), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

import core  # noqa: E402
import radius_metrics as rm  # noqa: E402

DEFAULT_RADII = tuple(round(0.2 + 0.1 * i, 2) for i in range(15))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--skim", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--radii", default=",".join(f"{r:g}" for r in DEFAULT_RADII))
    p.add_argument("--pt-min", type=float,
                   default=core.COLLECTIONS["GenFatJet"].pt_min,
                   help="Jet pT threshold. This interacts with R -- larger R "
                        "crosses threshold more often -- so bin results in "
                        "truth visible pT before comparing radii.")
    return p.parse_args()


def rebuild_particle_table(candidate) -> core.ParticleTable:
    """core.ParticleTable for one event, straight from the skim columns."""
    return core.ParticleTable(
        pt=np.asarray(candidate["pt"], dtype=np.float64),
        eta=np.asarray(candidate["eta"], dtype=np.float64),
        phi=np.asarray(candidate["phi"], dtype=np.float64),
        energy=np.asarray(candidate["energy"], dtype=np.float64),
        px=np.asarray(candidate["px"], dtype=np.float64),
        py=np.asarray(candidate["py"], dtype=np.float64),
        pz=np.asarray(candidate["pz"], dtype=np.float64),
        uid=np.asarray(candidate["uid"], dtype=np.int64),
    )


def rebuild_truth_event(event) -> core.FixedTruthEvent:
    """core.FixedTruthEvent for one event, from the stored owner labels.

    ``visible_indices`` and ``visible_by_dark_hadron`` are recovered from the
    per-candidate ``owner`` / ``dh_position`` columns rather than stored as
    nested index lists, which keeps the skim two-level jagged.
    """
    owner = np.asarray(event["candidate"]["owner"], dtype=np.int64)
    dh_pos = np.asarray(event["candidate"]["dh_position"], dtype=np.int64)
    has_anc = np.asarray(event["candidate"]["has_dh_ancestor"], dtype=bool)

    uids_all = ak.to_list(event["truth_dark_hadron_uids"])
    pts_all = ak.to_list(event["truth_dark_hadron_pts"])

    jets: list[core.FixedTruthJet] = []
    for slot, truth_id in enumerate(np.asarray(event["truth_id"], dtype=np.int64)):
        truth_id = int(truth_id)
        mine = owner == truth_id
        uids = tuple(int(u) for u in uids_all[slot])
        groups = tuple(
            tuple(int(i) for i in np.flatnonzero(mine & (dh_pos == p)))
            for p in range(len(uids))
        )
        jets.append(core.FixedTruthJet(
            truth_id=truth_id,
            dark_hadron_uids=uids,
            dark_hadron_pts=tuple(float(x) for x in pts_all[slot]),
            visible_indices=tuple(int(i) for i in np.flatnonzero(mine)),
            visible_by_dark_hadron=groups,
            pt=float(event["truth_pt"][slot]),
            eta=float(event["truth_eta"][slot]),
            phi=float(event["truth_phi"][slot]),
            mass=float(event["truth_mass"][slot]),
        ))

    return core.FixedTruthEvent(jets, owner, dh_pos, has_anc)


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    radii = [float(r) for r in args.radii.split(",") if r.strip()]

    skim = ak.from_parquet(args.skim)
    meta_path = args.skim.with_suffix(".meta.json")
    meta = json.loads(meta_path.read_text()) if meta_path.is_file() else {}
    n_events = len(skim)

    tables = [rebuild_particle_table(skim[i]["candidate"]) for i in range(n_events)]
    truths = [rebuild_truth_event(skim[i]) for i in range(n_events)]

    rows: dict[str, list] = {}

    def emit(**kwargs) -> None:
        for key, value in kwargs.items():
            rows.setdefault(key, []).append(value)

    # Built ONCE and reused across the whole radius grid. Assembling it costs
    # about as much as one clustering pass, so rebuilding it per radius gives
    # back most of the batching gain.
    particles = core.cluster_input(
        [t.px for t in tables], [t.py for t in tables],
        [t.pz for t in tables], [t.energy for t in tables],
    )

    t0 = time.perf_counter()
    for radius in radii:
        clustered_all = core.recluster_events(particles, radius, args.pt_min)
        for ev in range(n_events):
            table, truth = tables[ev], truths[ev]
            jets = [
                core.ClusteredJet(
                    pt=float(np.hypot(table.px[idx].sum(), table.py[idx].sum())),
                    eta=0.0, phi=0.0, mass=0.0, energy=float(table.energy[idx].sum()),
                    px=float(table.px[idx].sum()), py=float(table.py[idx].sum()),
                    pz=float(table.pz[idx].sum()),
                    constituent_indices=tuple(sorted(int(i) for i in idx)),
                )
                for idx in clustered_all[ev]
            ]
            matches = core.match_by_shared_visible_pt(truth, jets, table.pt)
            matched = {m.truth_id: m for m in matches}

            for truth_jet in truth.jets:
                base = dict(event=ev, truth_id=truth_jet.truth_id, radius=radius,
                            truth_jet_pt=truth_jet.pt,
                            n_dark_hadrons=len(truth_jet.dark_hadron_uids))
                n_eligible = sum(
                    1 for g in truth_jet.visible_by_dark_hadron if len(g))
                if truth_jet.truth_id in matched:
                    match = matched[truth_jet.truth_id]
                    metrics = core.matched_jet_metrics(
                        truth, truth_jet, jets[match.clustered_index], table
                    )
                    emit(**base, **rm.scalars(metrics),
                         any_dark_hadron_full=float(
                             metrics["n_fully_contained_dark_hadrons"] >= 1),
                         frac_xdh_pt_legacy=float(
                             metrics["cross_truth_dark_hadron_contamination_pt"]),
                         contamination_pt_legacy=float(metrics["contamination_pt"]),
                         jet_pt=float(metrics["jet_pt"]),
                         jet_mass=float(metrics["jet_mass"]))
                else:
                    emit(**base, **rm.unmatched_scalars(n_eligible),
                         any_dark_hadron_full=0.0,
                         frac_xdh_pt_legacy=np.nan,
                         contamination_pt_legacy=np.nan,
                         jet_pt=0.0, jet_mass=0.0)

    elapsed = time.perf_counter() - t0
    table_out = ak.Array({k: np.asarray(v, dtype=np.float64)
                          for k, v in rows.items()})

    out_meta = {
        "stage": "scan",
        "skim": str(args.skim.resolve()),
        "label": meta.get("label", "unknown"),
        "point_id": meta.get("point_id", -1),
        "params": meta.get("params", {}),
        "radii": radii,
        "pt_min": args.pt_min,
        "events": n_events,
        "rows": len(table_out),
        "seconds": round(elapsed, 2),
        "ms_per_event_per_radius": round(
            1e3 * elapsed / max(n_events, 1) / max(len(radii), 1), 3),
        "note": (
            "quality = acceptance * purity, both from core.matched_jet_metrics "
            "via radius_metrics.py. acceptance = fully-contained dark hadrons / "
            "eligible dark hadrons (a count). purity = 1 - frac_partial_pt - "
            "frac_nodh_pt. The *_legacy columns are core's ownership-based "
            "split, kept for continuity but not used by the scores. Unmatched "
            "truth jets carry matched=0, acceptance=0, quality=0."
        ),
    }
    ak.to_parquet(table_out, args.output, extensionarray=False)
    args.output.with_suffix(".meta.json").write_text(
        json.dumps(out_meta, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out_meta, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
