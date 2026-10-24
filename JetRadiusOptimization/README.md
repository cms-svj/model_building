# Jet Radius Optimization

Tools for the LHC Dark Showers Task Force study *"determine optimal jet radius
for each jet collection / stage of the shower"*.

Two workflows share one physics layer (`core.py`):

| | What it answers | Scale |
|---|---|---|
| **A. Validate** (`validate.py`) | Is the machinery correct, and what does R do on *this* sample? | one ROOT file, closure gates, event displays |
| **B. Scan** (`skim.py` → `scan.py` → `optimize.py`) | Which single radius works across *many* model points? | a parameter grid on HTCondor |

Workflow A is the original tool and is unchanged. Workflow B was added for
parameter scans; it calls `core.matched_jet_metrics` and
`core.match_by_shared_visible_pt` unchanged, so containment and contamination
have exactly one definition in this repository.

Neither workflow generates signal. Both read an existing `events.root` from
the patched `DarkHadronJets` Delphes configuration.

- **REFERENCE.md** — collection definitions, every plot, the full metric list,
  matching discipline, validation gates.
- **CHANGES.md** — what the scan patch changed in `core.py`, measured
  speedups, and corrections to earlier claims.

## Setup

Everything runs from a plain clone of this repository on an EL9 machine with
CVMFS (e.g. LPC). From the repository root:

```bash
git clone git@github.com:cms-svj/model_building
cd model_building
./install.sh          # once: Pythia, Delphes, and the python venv
source init.sh        # every new shell
python3 -m pip install "pyarrow==25.0.1"   # once, into the venv (see below)
```

`init.sh` gives you LCG_106 plus `coffea`, `fastjet` and `magiconfig` from
`install/python_packages.sh`. The Parquet stages also need pyarrow >= 17:
the awkward version that coffea pulls in refuses LCG_106's pyarrow 15 in
`ak.to_parquet`. 25.0.1 is the version this code was tested with.

## Quickstart: generate a sample and commission

Neither workflow generates signal. For a first run, make the 200-event CMS
benchmark sample with the repository's own generator (a few minutes):

```bash
./run_model helper -C configs/model_cms.py --steps all --events 200 \
  --dir JetRadiusOptimization/validation_sample
```

Then commission the checkout. One command runs every gate, including a closure
test of the skim/scan path against the original per-event path:

```bash
python3 JetRadiusOptimization/commission.py \
  --input JetRadiusOptimization/validation_sample/*/events.root
```

It ends with `READY` on success. A non-zero exit means don't proceed.
Generated samples and all outputs are git-ignored.

---

# Workflow A — validate one sample

Start here whenever you touch `core.py`, and start here on a new sample before
trusting any radius number from it.

## Smoke test, 10 events

```bash
# The Quickstart sample; for your own run_model output use models/<name>
MODEL_DIR="JetRadiusOptimization/validation_sample/s-channel_mmed-1000_Nc-2_Nf-2_scale-35.1539_mq-10_mpi-20_mrho-20_pvector-0.75_spectrum-cms_gq-0.25_gchi-0.5_rinv-0.3"

python3 JetRadiusOptimization/validate.py \
  --input   "$MODEL_DIR/events.root" \
  --outdir  "$MODEL_DIR/radius_diagnostics_10events" \
  --max-events 10 \
  --collections GenFatJet
```

Success ends with `{"status": "passed"}`. Ten events verifies the workflow and
gives you event displays to look at. It supports no physics conclusion.

## Full validation

Omit `--collections` to run closure on every clustered collection
(`GenJet GenFatJet DarkPartonJet DarkHadronJet Jet FatJet`):

```bash
python3 JetRadiusOptimization/validate.py \
  --input  "$MODEL_DIR/events.root" \
  --outdir "$MODEL_DIR/radius_validation" \
  --max-events 200
```

## Read the gates

```bash
python3 - "$MODEL_DIR/radius_validation/validation_report.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1]))
print("status:", report["status"])
for name, passed in report.get("gates", {}).items():
    print(f"{name:32s} {passed}")
PY
```

A failed gate is not a plotting warning. It means the associated radius result
should not be used. Gate meanings are in REFERENCE.md.

---

# Workflow B — scan a parameter grid

Three stages. The split exists because ROOT I/O and ancestry do not depend on
R, and re-running them per radius is the single biggest waste in the old
per-radius DAG.

```
skim.py     events.root      ──▶  skim.parquet       once per (point, chunk)
scan.py     skim.parquet     ──▶  metrics.parquet    all radii, no ROOT
plots.py    metrics × points ──▶  the decision figures + R*
```

The skim is also what makes the metric cheap to iterate: changing a definition
means re-running stage 2 only, which is seconds rather than a re-read of the
whole campaign.

## Stage 1 — skim

```bash
python3 JetRadiusOptimization/skim.py \
  --input   "$MODEL_DIR/events.root" \
  --output  "$MODEL_DIR/skim_test.parquet" \
  --max-events 200 \
  --label   "mmed-1000_mpi-20_rinv-0p3" \
  --params  '{"mmed":1000,"mpi":20,"mrho":20,"rinv":0.3,"Nc":2,"Nf":2}'
```

`--params` is free-form JSON carried through to the metrics table, so
`optimize.py` can regress R* against whatever you scanned. Use `--entry-start`
and `--max-events` to shard.

A sidecar `skim_test.meta.json` records the input path, entry range, repository
SHA, and counts.

## Stage 2 — scan

```bash
python3 JetRadiusOptimization/scan.py \
  --skim   "$MODEL_DIR/skim_test.parquet" \
  --output "$MODEL_DIR/metrics_test.parquet" \
  --radii  0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6
```

Default grid is 0.2 to 1.6 in steps of 0.1. Output is one flat row per
`(event, truth_id, radius)`, so everything downstream is a groupby.

> **Run this first.** On the 200-event `validation_sample`, check the R=0.8
> numbers reproduce what `validate.py` already gives you. That reproduction is
> the acceptance test for the whole patch. `commission.py` automates this
> check (step 6, the closure test).

## Understand the columns

Every constituent is classified by the **containment status of the dark hadron
it came from** — not by which truth jet owns that dark hadron.

| Class | Meaning | Effect on jet mass |
|---|---|---|
| FULL | from a dark hadron captured whole | decay products sum to its mass — real information |
| PARTIAL | from a dark hadron only fragments of which are in the jet | random fraction of its momentum, no mass meaning |
| NODH | no dark-hadron ancestor at all | pure background |
| ORPHAN | from a dark hadron that exists but was never assigned to any truth `DarkHadronJet` | currently uncounted — see caveat below |

PARTIAL is the worst case, and worse than missing a dark hadron entirely. Fully
contained, it tells you its mass. Absent, it costs acceptance but tells no
lies. Half in, the jet gets heavier and less correct at once.

| Column | Meaning |
|---|---|
| `acceptance` | fully-contained dark hadrons ÷ eligible dark hadrons. **A count** — a dark hadron is whole or it isn't |
| `purity` | `1 − frac_partial_pt − frac_nodh_pt` |
| `quality` | **the headline.** `acceptance × purity` |
| `frac_full_pt` | jet pT from dark hadrons captured whole |
| `frac_partial_pt` | jet pT from dark hadrons the jet shredded |
| `frac_nodh_pt` | jet pT with no dark-hadron ancestor |
| `frac_dh_pt` | jet pT from *any* dark hadron at all, whole/partial/orphan. `1 − frac_nodh_pt`, not `frac_full_pt + frac_partial_pt` — see caveat |
| `frac_orphan_dh_pt` | jet pT from an ORPHAN dark hadron (see caveat) |
| `matched` | 0 if no clustered jet matched this truth jet at this radius |
| `iou_pt` | pT-weighted Jaccard. Cross-check only, see below |
| `*_legacy` | `core`'s ownership-based split, kept for continuity |

A partially captured dark hadron is penalised **twice** — it fails to count in
`acceptance`, and its fragments drag down `purity`. One entirely outside the
jet costs only acceptance. Junk with no dark-hadron ancestor costs only purity.
So the score tracks the physics with no tuning constant.

**Caveat found 2026-09 while extending the mass scan.** FULL, PARTIAL, and
NODH do not actually partition jet pT to 1.0: `frac_full_pt` and
`frac_partial_pt` only sum over dark hadrons assigned to *some* truth jet's
`DarkHadronJet` partition (`core.dark_hadron_group_metrics` iterates
`truth.jets`), while `frac_nodh_pt`'s ancestry test
(`visible_has_any_dark_hadron_ancestor`) is resolved against *every*
`DarkHadronCandidate` in the event regardless of assignment. A dark hadron
that exists but was never assigned to a truth partition is therefore in
neither bucket — `frac_orphan_dh_pt` makes that residual explicit rather than
leaving it silently absorbed into `purity`. Measured mean ~0.014, up to 0.85
on individual jets, on both the DRAGON `mmed=2000` grid and the local
`mmed=1000` validation sample — so `purity`/`quality` are a slight
overestimate whenever it's non-zero. Not corrected in `core.py`; that's a
physics-definition call for whoever owns this metric, not a silent patch.

**Unmatched truth jets score `quality = 0`**, not NaN and not dropped.
Excluding them rewards radii that fail to reconstruct the object at all.

`iou_pt` is a cross-check, not the north star. It scores set overlap, which is
a different object from "how many complete dark hadrons survived", and it
weights constituents by pT linearly while jet mass responds to pT·ΔR² — so it
under-charges the soft wide junk a too-large radius admits.

## Stage 3 — the decision figures

```bash
python3 JetRadiusOptimization/plots.py \
  --metrics "$MODEL_DIR/metrics_*.parquet" \
  --outdir  "$MODEL_DIR/radius_decision_plots" \
  --x rinv --y mpi \
  --fixed mmed=1000 \
  --quantile 1.0
```

Five figures plus `decision_summary.json`:

**`optimal_radius_map.png` — the money plot.** One panel, colour = the radius
that wins at each grid point. Flat means one radius works everywhere and you
ship a constant. A gradient means the optimum tracks a physical scale and you
should be fitting `R(m_med, pT)` instead of quoting a number.

**`regret_map.png`** — cost of the compromise at the chosen global R, same
axes. Near zero everywhere means the recommendation holds. Wherever it isn't is
what you caveat in the writeup.

**`quality_vs_radius.png`** — one curve per model point. This is where the
**plateau** lives. Whether the peak is sharp or whether 0.8–1.0 are
indistinguishable is invisible in any heatmap, and it's the difference between
"R = 0.9" and "anything from 0.8 to 1.0, we suggest 0.8".

**`decomposition_maps.png`** — why a corner fails: acceptance,
`frac_partial_pt`, `frac_nodh_pt`, and matched fraction at the chosen radius.
Losing whole dark hadrons, shredding them, or drowning in UE are three
different problems with three different fixes.

**`regret_by_radius.png`** — the per-radius panel grid. Note it plots **regret,
not raw quality**: regret is normalised per model point so cells are comparable
across panels. Raw quality is not — its colour scale gets consumed by "which
model point is easy" rather than "which radius is good", which is why a grid of
raw-metric heatmaps hides the very thing you're looking for.

### Reading them

Three things will bite:

**`mmed` is a hidden third axis** and drives the boost harder than `rinv` or
`mpi`. Use `--fixed mmed=1000`, or the map is a projection averaging over it.
The value is printed in the title either way.

**High-`rinv` cells go sparse.** Fewer visible constituents means jets fall
below `pt_min` and matching drops, so the `rinv = 0.9` row goes noisy in a way
that looks like physics. Read the matched-fraction panel before trusting it.
Cells with under 50 jets are dropped automatically.

**4 × 6 is a table, not a heatmap.** To *fit* a gradient you want ~8 points per
axis. The scan is cheap; generation is the constraint.

`plots.py` warns if any model point's argmax sits on the first or last radius —
the optimum may lie outside the grid, and the scaling fit will be biased toward
the boundary. Widen `--radii` rather than quoting the edge.

### Scaling fit

If the map shows a gradient, test whether it collapses onto one curve:

```python
import sys; sys.path.insert(0, "JetRadiusOptimization")
from optimize import fit_scaling
print(fit_scaling(best_radius_per_model, predictor_per_model))
```

Candidates in rough order of expected power: `2·m_med/pT_visible`,
`m_dark/pT_visible`, `Λ/m_med`, `m_π/Λ`. A high R² turns the deliverable into a
formula valid outside the scanned grid — much more useful to experiments than a
constant. A low one justifies recommending a fixed value. Either is a result.

---

# Batch submission

The old `condor/focused_by_radius.dag` sharded over **radius**: 15 jobs, same
10,000 events, one radius each. That copied `events.root` 15 times, redid the
radius-independent ancestry 15 times, and pinned the job count at 15. It still
works; use it if you only want the focused diagnostic plots refreshed. Steps are
in [condor/README.md](condor/README.md).

`condor/scan.sub` shards over **(model point × event chunk)** with all radii
inside each job.

```bash
mkdir -p JetRadiusOptimization/condor/logs

python3 JetRadiusOptimization/condor/make_jobs.py \
  --outdir  JetRadiusOptimization/condor \
  --mmed    500 1000 2000 3000 \
  --mdark   5 10 20 40 \
  --rinv    0.0 0.1 0.3 0.5 0.7 0.9 \
  --events-per-point 100000 \
  --eos-base /store/user/$USER/DarkSectorStudies/MC/radius_grid
```

That writes `points.txt` (the parameter table, same column convention as
`FlatSignal_generationcampaign`, so generation and analysis share a `point_id`)
and `jobs.txt` (the work list), then prints the job count and a core-hour
estimate. Chunk size is set from a target wall time, not a round number of
events — 5-minute jobs lose more to startup and transfer than they gain in
latency.

> **Not yet commissioned.** No `scan.sub` job has been submitted yet, and
> three gaps are known: `scan_worker.sh` does not set up LCG or the venv;
> it expects the payload under `payload/`; and the chunk size that
> `make_jobs.py` computes is not passed to the worker as `EVENTS_PER_CHUNK`.
> Fix these and run one job before launching a campaign.

`scan.sub` uses paths relative to `condor/`, so submit from there:

```bash
cd JetRadiusOptimization/condor
condor_submit scan.sub
condor_q $USER -batch
```

Resource requests are 2000 MB memory and 4000 MB disk, roughly 2× the
934–1007 MB and 2.28 GB observed in the existing logs. The old submit file
asked for 8 GB of each, which cuts how many LPC slots you match for no benefit.

`SEC_PER_EVENT_PER_RADIUS` in `make_jobs.py` is a **placeholder estimate**. It
only sets chunk size, so being wrong is harmless, but replace it once you have
timed stage 2 for real.

---

# Three traps in the metric

Worth settling on the sample you already have, before spending core-hours.

**The `pt_min` threshold interacts with R.** Larger R sweeps up more pT and
crosses 15 GeV more often, which appears in the acceptance-vs-R curve as
physics when it is selection. Bin in truth visible pT.

**`rinv` partly acts through the boost.** Raising it removes visible pT, which
lowers the visible jet pT, which widens the optimal radius — the same mechanism
as lowering `mmed`, not an independent effect. Check whether `rinv` still moves
R* after regressing on `2·m_med/pT_visible`.

**Cross-check against mass response.** The radius minimising the width of
`m_SD / m_dark` is what the `StrategiesForTrainingData` mass regression
actually cares about. If it disagrees with the IoU optimum, the disagreement is
the interesting result and IoU needs reweighting.

---

# Troubleshooting

**`events.root` is missing.** Generate signal separately with `run_model` or
the `HTCondorSignalGeneration/` tutorial. Neither workflow silently generates an
input sample.

**A required branch is missing.** The input must come from the patched
`DarkHadronJets` Delphes configuration. `validate.py` reports the exact missing
branches and exits nonzero.

**`validation_report.json` says `failed`.** Read `gates`. Do not use the radius
results from that run.

**`ImportError: pyarrow 17.0.0 or later required`.** LCG_106's pyarrow 15 is
too old for awkward's Parquet I/O. Run `source init.sh`, then
`python3 -m pip install "pyarrow==25.0.1"`.

**`ancestry bitmask did not reach a fixed point`.** The mother graph has an
unexpected structure. Inspect the event rather than raising `max_passes` — the
fixed point is reached in 10–13 passes on normal input.

**`N targets exceeds the 64-bit mask`.** More than 64 dark hadrons assigned to
DarkHadronJets in one event. The CMS sample peaks at 19. If a model genuinely
exceeds it, split the target set rather than widening the guard silently.

**Plots look jagged or empty.** Not enough events. A 10-event run is a workflow
demonstration.

**`Predicted rinv = nan` during generation.** For the simplified CMS
configuration `rinv` is supplied directly rather than predicted from a flavour
count. It does not invalidate the configured `rinv` or the containment plots.

---

# Reproducibility

Reports record the input ROOT path, repository SHA, patched Delphes SHA,
`PYTHIA8MINOR`, package versions, radius grid, event count, and matching
policy. Resolved generation cards and logs stay beside `events.root`.

Do not combine results from different Pythia versions without labelling them.
The shipped 10-event diagnostic sample predates the current local install,
which reports Pythia 8.317.

`benchmarks/` reproduces the measured speedups. `bench_seeded.py` is kept
deliberately: clustering once at `R_max` and reclustering inside each seed jet
is **not** equivalent to anti-kT at smaller R (45 of 50 events disagreed at
R=0.8) and was far slower. Do not try it again.
