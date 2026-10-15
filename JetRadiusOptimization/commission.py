#!/usr/bin/env python3
"""Commissioning: prove this checkout works on THIS machine before using it.

Run once on a fresh checkout, and again after changing core.py, skim.py or
scan.py. It checks the environment, the two benchmark correctness checks, the
validate.py gates, and a closure test of the skim/scan path against the
original per-event path.

    python3 JetRadiusOptimization/commission.py --input <events.root>

Every step prints PASS or FAIL with a diagnostic. Non-zero exit means do not
proceed to a campaign.

The important step is 6. It recomputes the metrics through the ORIGINAL
per-event path (build_fixed_truth_event -> recluster_event ->
matched_jet_metrics) and diffs them jet-by-jet, radius-by-radius, against what
skim.py + scan.py produced. If those agree the two-stage refactor is proven on
real data; if they do not, the scan output is wrong no matter how good the
plots look.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
for p in (str(HERE.parent), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

RADII = [0.4, 0.8, 1.2]
TOLERANCE = 1e-9

results: list[tuple[str, bool, str]] = []


def step(name: str, ok: bool, detail: str = "") -> bool:
    results.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""),
          flush=True)
    return ok


def run(cmd: list[str], cwd: Path | None = None):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def reference_metrics(input_path: Path, n_events: int, radii, pt_min):
    """Metrics via the ORIGINAL per-event path, with nothing from skim/scan."""
    from coffea.nanoevents import NanoEventsFactory
    from common import DelphesSchema2
    import core
    import radius_metrics as rm

    events = NanoEventsFactory.from_root(
        {str(input_path.resolve()): "Delphes"},
        schemaclass=DelphesSchema2, entry_start=0, entry_stop=n_events,
    ).events()
    raw = core.load_raw_ancestry(input_path, entry_start=0, entry_stop=n_events)

    out = {}
    for i in range(len(events)):
        table = core.particle_table(events[i].GenCandidate)
        truth = core.build_fixed_truth_event(events[i], raw[i])
        for radius in radii:
            jets = core.recluster_event(table, radius, pt_min)
            matched = {m.truth_id: m.clustered_index
                       for m in core.match_by_shared_visible_pt(truth, jets, table.pt)}
            for tj in truth.jets:
                key = (i, tj.truth_id, round(float(radius), 6))
                if tj.truth_id in matched:
                    out[key] = rm.scalars(
                        core.matched_jet_metrics(truth, tj, jets[matched[tj.truth_id]], table))
                else:
                    n_elig = sum(1 for g in tj.visible_by_dark_hadron if len(g))
                    out[key] = rm.unmatched_scalars(n_elig)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, required=True, help="An events.root.")
    ap.add_argument("--events", type=int, default=200)
    ap.add_argument("--workdir", type=Path, default=None)
    ap.add_argument("--skip-validate", action="store_true",
                    help="Skip the (slow) validate.py gate run.")
    args = ap.parse_args()

    work = args.workdir or Path(tempfile.mkdtemp(prefix="commission_"))
    work.mkdir(parents=True, exist_ok=True)
    print(f"input   : {args.input}\nevents  : {args.events}\nworkdir : {work}\n", flush=True)

    if not args.input.is_file():
        return step("input exists", False, f"not found: {args.input}") or 1

    # --- 1. environment ---------------------------------------------------
    missing = []
    for mod in ("numpy", "awkward", "uproot", "fastjet", "coffea", "pyarrow",
                "matplotlib", "scipy"):
        try:
            __import__(mod)
        except ImportError:
            missing.append(mod)
    if not step("environment", not missing,
                f"missing: {', '.join(missing)} -- pip install them" if missing
                else "all imports present"):
        return 1
    import pyarrow
    if not step("pyarrow >= 17", int(pyarrow.__version__.split(".")[0]) >= 17,
                f"found {pyarrow.__version__}" +
                ("" if int(pyarrow.__version__.split(".")[0]) >= 17 else
                 '; awkward\'s Parquet I/O needs >= 17 -- '
                 'python3 -m pip install "pyarrow==25.0.1"')):
        return 1

    try:
        import core
        import radius_metrics as rm  # noqa: F401
        import ancestry, optimize, plots  # noqa: F401
    except Exception as exc:
        return step("repo modules import", False, repr(exc)) or 1
    step("repo modules import", True)

    pt_min = core.COLLECTIONS["GenFatJet"].pt_min

    # --- 2 & 3. the two benchmarks that double as correctness checks ------
    for script, needle in (("bench_ancestry.py", "PHYSICS IDENTICAL"),
                           ("bench_cluster.py", "CLUSTERING IDENTICAL")):
        path = HERE / "benchmarks" / script
        if not path.is_file():
            step(f"benchmark {script}", False, "not found")
            continue
        proc = run([sys.executable, str(path), str(args.input)], cwd=HERE.parent)
        step(f"benchmark {script}", needle in proc.stdout,
             (proc.stdout.strip().splitlines() or ["no output"])[-1])

    # --- 4. validate.py gates --------------------------------------------
    if args.skip_validate:
        step("validate.py gates", True, "skipped by request")
    else:
        vout = work / "validate"
        proc = run([sys.executable, str(HERE / "validate.py"),
                    "--input", str(args.input), "--outdir", str(vout),
                    "--max-events", str(min(args.events, 200))])
        report = vout / "validation_report.json"
        if report.is_file():
            data = json.loads(report.read_text())
            gates = data.get("gates", {})
            failed = [k for k, v in gates.items() if not v]
            step("validate.py gates", data.get("status") == "passed" and not failed,
                 f"status={data.get('status')}" +
                 (f", failed: {', '.join(failed)}" if failed else ""))
        else:
            step("validate.py gates", False,
                 (proc.stderr.strip().splitlines() or ["no report written"])[-1])

    # --- 5. skim.py ----------------------------------------------------
    skim_path = work / "commission_skim.parquet"
    t0 = time.perf_counter()
    proc = run([sys.executable, str(HERE / "skim.py"),
                "--input", str(args.input), "--output", str(skim_path),
                "--max-events", str(args.events), "--label", "commission",
                "--params", json.dumps({"mmed": 0, "mpi": 0, "rinv": 0})])
    skim_secs = time.perf_counter() - t0
    if not step("skim.py runs", proc.returncode == 0 and skim_path.is_file(),
                f"{skim_secs:.1f}s, {skim_path.stat().st_size/1e6:.1f} MB"
                if skim_path.is_file()
                else (proc.stderr.strip().splitlines() or ["no stderr"])[-1]):
        print("\nskim.py failed. Read the traceback above.", flush=True)
        return 1

    # --- 6. scan.py, then the closure test that matters -------------------
    metrics_path = work / "commission_metrics.parquet"
    t0 = time.perf_counter()
    proc = run([sys.executable, str(HERE / "scan.py"),
                "--skim", str(skim_path), "--output", str(metrics_path),
                "--radii", ",".join(f"{r:g}" for r in RADII)])
    scan_secs = time.perf_counter() - t0
    if not step("scan.py runs", proc.returncode == 0 and metrics_path.is_file(),
                f"{scan_secs:.1f}s" if metrics_path.is_file()
                else (proc.stderr.strip().splitlines() or ["no stderr"])[-1]):
        return 1

    import awkward as ak
    import numpy as np

    table = ak.from_parquet(metrics_path)
    got = {}
    cols = [c for c in ak.fields(table)
            if c not in ("event", "truth_id", "radius")]
    ev = np.asarray(table["event"]).astype(int)
    ti = np.asarray(table["truth_id"]).astype(int)
    rr = np.asarray(table["radius"])
    arrays = {c: np.asarray(table[c]) for c in cols}
    for k in range(len(ev)):
        got[(ev[k], ti[k], round(float(rr[k]), 6))] = {c: arrays[c][k] for c in cols}

    print("\n  computing reference via the original per-event path ...", flush=True)
    t0 = time.perf_counter()
    ref = reference_metrics(args.input, args.events, RADII, pt_min)
    ref_secs = time.perf_counter() - t0

    missing_keys = set(ref) - set(got)
    extra_keys = set(got) - set(ref)
    diffs = []
    for key, refrow in ref.items():
        gotrow = got.get(key)
        if gotrow is None:
            continue
        for name, rv in refrow.items():
            gv = gotrow.get(name)
            if gv is None:
                continue
            rv, gv = float(rv), float(gv)
            if math.isnan(rv) and math.isnan(gv):
                continue
            if math.isnan(rv) != math.isnan(gv) or abs(rv - gv) > TOLERANCE:
                diffs.append((key, name, rv, gv))

    ok = not missing_keys and not extra_keys and not diffs
    step("CLOSURE: skim+scan == direct core path", ok,
         f"{len(ref)} (jet, radius) rows compared, {len(diffs)} value diffs, "
         f"{len(missing_keys)} missing, {len(extra_keys)} extra")
    if not ok:
        for key, name, rv, gv in diffs[:10]:
            print(f"      {key} {name}: direct={rv!r} scan={gv!r}", flush=True)
        if missing_keys:
            print(f"      missing from scan output: {sorted(missing_keys)[:5]}", flush=True)
        return 1

    speedup = ref_secs / max(scan_secs + skim_secs, 1e-9)
    print(f"      direct path {ref_secs:.1f}s vs skim+scan {skim_secs + scan_secs:.1f}s "
          f"for {len(RADII)} radii  ->  x{speedup:.1f}", flush=True)

    # --- 7. plots.py smoke ------------------------------------------------
    proc = run([sys.executable, str(HERE / "plots.py"),
                "--metrics", str(metrics_path), "--outdir", str(work / "plots"),
                "--x", "rinv", "--y", "mpi"])
    made = sorted(p.name for p in (work / "plots").glob("*.png")) \
        if (work / "plots").is_dir() else []
    step("plots.py runs", len(made) >= 3,
         f"{len(made)} figures" if made
         else (proc.stderr.strip().splitlines() or ["no figures"])[-1])

    # --- verdict ----------------------------------------------------------
    n_fail = sum(1 for _, ok_, _ in results if not ok_)
    print("\n" + "=" * 66)
    if n_fail:
        print(f"NOT READY -- {n_fail} step(s) failed. See above.")
    else:
        print("READY. Closure holds on real data; the two-stage path is proven.")
        print("\nRemaining before a campaign, in order:")
        print("  1. time stage 2 on a realistic chunk and put the real number")
        print("     into SEC_PER_EVENT_PER_RADIUS in condor/make_jobs.py")
        print("  2. submit ONE condor job and confirm it lands on EOS")
        print("  3. profile the ~40 ms/event remainder before scaling up")
        print(f"\n  measured here: skim {skim_secs*1000/args.events:.1f} ms/event, "
              f"scan {scan_secs*1000/args.events/len(RADII):.2f} ms/event/radius")
    print("=" * 66, flush=True)
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
