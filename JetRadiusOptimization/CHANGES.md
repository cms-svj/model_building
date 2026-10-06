# Changes to JetRadiusOptimization

A patch, not a parallel framework. Existing files keep their behaviour; the
additions import from `core.py` rather than reimplementing it.

Unzip over your existing `JetRadiusOptimization/`. The ROOT samples and
`condor/logs` are excluded, so your copies stay where they are.

## Corrections to what I told you earlier

Three of my earlier numbers were wrong. In order of how much they matter:

**1. The end-to-end speedup is not ~30x.** I extrapolated it from two
component benchmarks without checking they accounted for the observed runtime.
They do not. The measured per-event budget at one radius:

```
ancestry (both BFS queries)   15.43 ms   ->  0.94 ms   fixed, 16.5x
clustering                     2.82 ms   ->  1.35 ms   fixed,  2.1x
soft drop                      1.21 ms   ->     --     NOT fixed
observed total (your logs)    59.60 ms
unexplained remainder         40.14 ms   ->     --     NOT profiled
```

That 40 ms is the bulk of your runtime and I never looked at it. It is
matching, `matched_jet_metrics`, and — most likely dominant — per-event
`events[i]` indexing into the coffea/awkward record, plus `ak.to_list` on the
jet constituent refs inside `build_fixed_truth_event`. Per-event indexing of
an awkward array is slow, and you do it once per event per job.

The skim helps here not by making that work faster but by doing it **once**
instead of once per radius, and by leaving stage 2 in pure numpy. A defensible
estimate is roughly **5-8x** end to end for a 15-radius scan, most of it from
removing the 15x job duplication rather than from the vectorisation. I cannot
measure it properly without running `skim.py`, which needs coffea.

**2. Clustering is 2.1x, not 7.6x.** My earlier benchmark compared against a
`PseudoJet` loop fed *awkward* arrays. `core.recluster_event` is fed a
`ParticleTable` holding *numpy* arrays, which is much faster. Against the real
baseline the gain is 2.1x. Building the awkward record once and reusing it
across the grid matters — rebuilding per radius gave back most of the gain
(1.6x), which is why `cluster_input` is a separate function.

**3. `dark_hadron_group_metrics` did not have the NODH bug I claimed.** It
already uses `visible_has_any_dark_hadron_ancestor` and is correct. See below
for what the actual discrepancy is.

The ancestry number went **up**: 16.5x rather than 8.2x, because both
`ancestor_targets` and `has_ancestor_target` are now vectorized, and the
earlier figure only replaced the first.

## Two contamination definitions coexist in core.py

`matched_jet_metrics` returns both of these, and they are different quantities:

| Field | Numerator | Basis |
|---|---|---|
| `contamination_pt` | `clustered_indices - all_truth_visible` | owner (`visible_owner == -1`) |
| `non_dark_hadron_constituent_pt_fraction` | `not visible_has_any_dark_hadron_ancestor` | ancestry |

The first is a **superset** of the second. A candidate whose ancestry is
ambiguous across several truth jets gets `owner == -1` from the algo-20 walk
but still descends from a dark hadron, so it counts as contamination under the
first definition and not under the second. It should be charged as cross-truth
contamination, not as unrelated junk.

`README.md` line 251 describes `contamination_pt` as *"the fraction of
clustered-jet constituent pT that does not descend from any fixed truth
dark-hadron jet"* — which is what the **second** field computes, not the first.
Worth reconciling.

I did not change `contamination_pt`. Existing plots and closure gates consume
it, and silently redefining it would invalidate results you already have.
`radius_metrics.py` uses the ancestry-based pair and `scan.py` also writes the
owner-based value as `contamination_pt_legacy` so the two can be compared.

## Modified

### `core.py`

Four edits, all inside `build_fixed_truth_event`, plus two new functions.

- imports `ancestor_masks`, `descends_from_any`, `targets_for_particle` from
  the new `ancestry.py`;
- resolves both ancestry questions once per event before the candidate loop
  instead of once per candidate;
- the algo-20 cache-ordered owner walk is **untouched** — its cache is already
  event-global so it is cheap, and its ordering is what your closure gate
  checks;
- adds `partially_contained_dark_hadron_pt_fraction`,
  `fully_contained_dark_hadron_pt_fraction`, and the matching counts to
  `dark_hadron_group_metrics`, forwarded through `matched_jet_metrics`. Purely
  additive — no existing field changed value;
- adds `cluster_input()` and `recluster_events()` for batched clustering.
  `recluster_event()` is left in place; the closure gates should keep
  exercising the simple path.

Validated on all 200 events of `validation_sample`:

```
has_ancestor_target diffs : 0 / 120,004
ancestor_targets   diffs  : 0 / 120,004
partition mismatches      : 0 / 900   (15 radii x 60 events)
```

`ancestor_targets` and `has_ancestor_target` keep their signatures and still
work; nothing in the call path is now using them.

## Added

| File | What |
|---|---|
| `ancestry.py` | Vectorized bitmask ancestry. `ancestor_masks`, `descends_from_any`, `targets_for_particle`, `owner_labels`. |
| `radius_metrics.py` | Derived scalars only — acceptance/purity/IoU **from** `matched_jet_metrics` output. Recomputes nothing. |
| `skim.py` | Stage 1: `events.root` -> radius-independent Parquet skim. |
| `scan.py` | Stage 2: skim -> per-jet metrics at all radii. Rebuilds `core` dataclasses and calls `core.match_by_shared_visible_pt` and `core.matched_jet_metrics` unchanged. |
| `optimize.py` | Cross-model radius choice by minimax regret, plus scaling-law fit. |
| `plots.py` | The decision figures: optimal-radius map, regret map, quality-vs-R curves, decomposition panels, per-radius regret grid. |
| `condor/scan.sub`, `scan_worker.sh`, `make_jobs.py` | (point x chunk) campaign, replacing the per-radius DAG. |
| `benchmarks/` | The three benchmarks, including the seeding idea that fails. |

## The metric

`core` already had the pieces except one. The missing quantity was **the pT
fraction of the jet coming from dark hadrons it did not capture whole**, which
is now `partially_contained_dark_hadron_pt_fraction` in
`dark_hadron_group_metrics`.

Constituents are classified by the **containment status** of their source dark
hadron, not by ownership:

```
FULL      from a dark hadron captured whole by this jet
PARTIAL   from a dark hadron only fragments of which are in this jet
NODH      no dark-hadron ancestor at all
```

Ownership was the wrong axis. `cross_truth_dark_hadron_contamination_pt` lumps
a neighbour's *fully contained* dark hadron (harmless, carries its own correct
mass) together with a *half-caught* one (poison), and it misses a dark hadron
of this very truth jet that the jet shredded. Both `*_legacy` columns are still
written, but nothing scores on them.

PARTIAL is the worst of the three, and worse than missing a dark hadron
outright. Fully contained, its decay products sum to its mass. Absent, it costs
acceptance but tells no lies. Half in, the jet gains a random fraction of its
momentum with no mass meaning — heavier and less correct at once.

```
A = n_fully_contained_dark_hadrons / n_eligible_dark_hadrons
P = 1 - frac_partial_pt - frac_nodh_pt
quality = A * P
```

A is a **count**, deliberately: a dark hadron is whole or it isn't, and
pT-weighting a binary means nothing.

A partially captured dark hadron is penalised **twice** — it fails to count in
A's numerator and its fragments drag down P. One entirely outside the jet costs
only A. Junk with no dark-hadron ancestor costs only P. No tuning constant.

`iou_pt` is demoted to a cross-check. It scores set overlap, a different object
from "how many complete dark hadrons survived", and its linear pT weighting
under-charges the soft wide junk a large radius admits, since jet mass responds
to pT·ΔR².

Unmatched truth jets get `quality = 0`, not NaN and not omitted.

## Choosing one radius

`optimize.choose_radius` works in regret, `regret_m(R) = Q_m(R_m*) - Q_m(R)`,
and minimises its high quantile across model points. Averaging the metric and
taking the argmax lets easy models mask hard ones.

Report `quantile=1.0` and `quantile=0.9` both. If they disagree, one radius
genuinely cannot cover the grid — a finding, not a problem to average away.
The plateau and the outlier list matter more than the number.

`optimize.fit_scaling` tests whether per-model optima collapse against
`2*m_med/pT_visible` and friends. If they do, the deliverable is a formula
valid outside the scanned grid, which is what "directly provide generator
settings" wants. If they do not, the low R-squared justifies a fixed value.

## Three traps in the metric

- **`pt_min` interacts with R.** Larger R crosses the 15 GeV threshold more
  often, which looks like physics and is selection. Bin in truth visible pT.
- **`rinv` partly acts through the boost.** It removes visible pT, lowering jet
  pT, widening the optimum — same mechanism as lowering `mmed`. Check whether
  it still moves R* after regressing on `2*m_med/pT_visible`.
- **Cross-check against mass response.** The radius minimising the width of
  `m_SD / m_dark` is what your `StrategiesForTrainingData` regression cares
  about. Disagreement with the IoU optimum is the interesting result.

## Order of work

1. `skim.py` then `scan.py` on the 200-event `validation_sample`; check the
   R=0.8 numbers reproduce `validate.py`. **This is the acceptance test for
   the whole patch.** `skim.py` is the one file I could not execute — it needs
   coffea and the patched Delphes branches — so treat its first run as
   debugging.
2. Re-run your closure gates, especially gate 3 (truth ancestry) and gate 2
   (default-radius closure).
3. Time stage 2 properly and put the real number into
   `condor/make_jobs.py`, where `SEC_PER_EVENT_PER_RADIUS` is currently a
   placeholder estimate that only sets chunk size.
4. Profile the 40 ms remainder before assuming the scan is fast enough.
5. Fix the resource requests: you ask for 8 GB and 8 GB, you use ~1 GB and
   2.28 GB. `condor/scan.sub` uses 2000/4000 MB.
6. Only then generate the `rinv` grid.
