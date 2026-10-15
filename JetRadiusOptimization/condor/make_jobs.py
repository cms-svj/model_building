#!/usr/bin/env python3
"""Build the (point_id, chunk) work list and the parameter table for a scan.

Two files come out:

  points.txt   one row per model point: the dark-sector parameters, plus the
               EOS path of its events.root.  Same column convention as
               FlatSignal_generationcampaign/*/mass_points_*.txt so generation
               and analysis stay indexed by a shared point_id.

  jobs.txt     one row per HTCondor job: (point_id, chunk).  Chunk size is
               chosen from a target wall time rather than a round number of
               events -- 5-minute jobs waste more in startup and transfer than
               they save in latency, and LPC will not thank you for 40,000 of
               them.

Sizing note.  Timed 2026-09-06 on the 200-event `validation_sample` (3 repeats
each, subtracting process-import overhead): skim.py measured ~61 ms/event,
scan.py ~2.4 ms/event/radius.  At 15 radii that is ~97 ms/event -- about 4.6x
the ~21 ms/event this file assumed before anyone had run skim.py.  This was
measured on local disk with no xrootd/EOS read; a real condor job pulling
events.root over the network will likely see a higher skim-stage number.
Remeasure at the real campaign scale before trusting the chunk sizing here.
"""

from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path

# Per-event cost model, seconds.  Measured 2026-09-06 on the 200-event
# validation_sample (local disk, no xrootd) -- see the module docstring.
# Update from a real timing run at campaign scale; these only set the chunk
# size, so being 30% off is harmless, but the local number above was off by
# 4.6x from the previous placeholder, so don't assume this one is final either.
SEC_PER_EVENT_ANCESTRY = 61.0e-3  # skim stage: coffea event access dominates, not ancestry
SEC_PER_EVENT_PER_RADIUS = 2.4e-3  # scan stage: clustering + matching + metrics


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", type=Path, default=Path("."))
    p.add_argument("--events-per-point", type=int, default=100_000,
                   help="Generated events available at each model point.")
    p.add_argument("--n-radii", type=int, default=15)
    p.add_argument("--target-minutes", type=float, default=45.0,
                   help="Wall-time target per job; sets the chunk size.")
    p.add_argument("--eos-base", default=None,
                   help="EOS directory holding <point_label>/events.root.")
    p.add_argument("--mmed", type=float, nargs="+", default=[1000.0])
    p.add_argument("--mdark", type=float, nargs="+", default=[20.0],
                   help="m_pi; m_rho is taken equal unless --mrho is given.")
    p.add_argument("--mrho", type=float, nargs="+", default=None)
    p.add_argument("--rinv", type=float, nargs="+",
                   default=[0.0, 0.1, 0.3, 0.5, 0.7, 0.9])
    p.add_argument("--nc", type=int, nargs="+", default=[2])
    p.add_argument("--nf", type=int, nargs="+", default=[2])
    p.add_argument("--pvector", type=float, nargs="+", default=[0.75])
    return p.parse_args()


def chunk_events(n_radii: int, target_minutes: float) -> int:
    per_event = SEC_PER_EVENT_ANCESTRY + n_radii * SEC_PER_EVENT_PER_RADIUS
    raw = target_minutes * 60.0 / per_event
    # round to a clean 5k so the last chunk is not a sliver
    return max(5_000, int(round(raw / 5_000.0)) * 5_000)


def label(mmed, nc, nf, mpi, mrho, pvector, rinv) -> str:
    def f(x):
        return f"{x:g}".replace(".", "p")
    return (f"mmed-{f(mmed)}_Nc-{nc}_Nf-{nf}_mpi-{f(mpi)}_mrho-{f(mrho)}"
            f"_pvector-{f(pvector)}_rinv-{f(rinv)}")


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    mrho_list = args.mrho if args.mrho is not None else None
    rows = []
    point_id = 0
    for mmed, nc, nf, mpi, pvector, rinv in itertools.product(
        args.mmed, args.nc, args.nf, args.mdark, args.pvector, args.rinv
    ):
        for mrho in (mrho_list if mrho_list is not None else [mpi]):
            name = label(mmed, nc, nf, mpi, mrho, pvector, rinv)
            path = f"{args.eos_base}/{name}/events.root" if args.eos_base else "AUTO"
            rows.append((point_id, mmed, nc, nf, mpi, mrho, pvector, rinv, name, path))
            point_id += 1

    points = args.outdir / "points.txt"
    with points.open("w") as fh:
        fh.write("# point_id mmed Nc Nf mpi mrho pvector rinv label events_root\n")
        for r in rows:
            fh.write("{} {:g} {} {} {:g} {:g} {:g} {:g} {} {}\n".format(*r))

    per_chunk = min(chunk_events(args.n_radii, args.target_minutes),
                    args.events_per_point)
    n_chunks = math.ceil(args.events_per_point / per_chunk)
    jobs = args.outdir / "jobs.txt"
    with jobs.open("w") as fh:
        for r in rows:
            for c in range(n_chunks):
                fh.write(f"{r[0]} {c}\n")

    n_jobs = len(rows) * n_chunks
    per_event = SEC_PER_EVENT_ANCESTRY + args.n_radii * SEC_PER_EVENT_PER_RADIUS
    core_hours = len(rows) * args.events_per_point * per_event / 3600.0
    print(f"model points        : {len(rows)}")
    print(f"events per point    : {args.events_per_point:,}")
    print(f"chunk size          : {per_chunk:,} events "
          f"(~{per_chunk*per_event/60:.0f} min/job)")
    print(f"chunks per point    : {n_chunks}")
    print(f"total jobs          : {n_jobs:,}")
    print(f"estimated core-hours: {core_hours:,.0f}")
    print(f"wrote {points} and {jobs}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
