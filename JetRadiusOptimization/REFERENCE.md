# Jet Radius Optimization: reference

Detailed reference for the single-sample validation workflow: collection
definitions, every generated plot, the metrics in `validation_report.json`,
matching discipline, and the validation gates.

For how to *run* things, start at **README.md**. For what changed in the
two-stage scan patch, see **CHANGES.md**.

This directory studies how the anti-$k_T$ radius of generator-level fat jets
changes dark-hadron descendant containment, contamination, jet mass, and Soft
Drop mass.

It is a **local physics-analysis and validation tool**. It reads an existing
Delphes `events.root`; it does not generate signal, submit HTCondor jobs, create
training data, or write to EOS. General signal generation lives separately in
`HTCondorSignalGeneration/`.

## Current scope

The implemented primary scan is:

```text
fixed DarkHadronJet truth partition at R=0.8
                    ↓
visible stable descendants found through exact generator ancestry
                    ↓
anti-kT GenFatJets re-clustered at R=0.2 ... 1.6
                    ↓
one-to-one matching by maximum shared visible-descendant pT
                    ↓
containment + contamination + mass + jet-shape diagnostics
```

The framework currently provides a validated **GenFatJet radius study for one
local ROOT sample**. It does not yet produce a final global optimal-radius
recommendation across the full model-parameter grid.

## Collection definitions

- `DarkHadronJet`: anti-$k_T$ clustering of initial dark hadrons immediately
  after fragmentation and before decay. Its shipped $R=0.8$ partition defines
  the fixed truth objects used in this scan.
- `DarkHadronStableJet`: manual ancestry combination of all stable descendants,
  visible and invisible, belonging to a `DarkHadronJet`.
- `DarkHadronVisibleJet`: manual ancestry combination of visible stable
  descendants belonging to a `DarkHadronJet`.
- `GenFatJet`: anti-$k_T$ clustering of stable visible generator particles. This
  is the collection whose radius is varied by the current primary scan.

`DarkHadronVisibleJet` and `DarkHadronStableJet` use the patched Delphes ancestry
matching algorithm, not ordinary radius-based clustering. Their `ParameterR`
does not control membership and must not be interpreted as a containment radius.

Changing the `DarkHadronJet` clustering radius would redefine the truth
partition itself. That is a separate upstream-coupled study and is not mixed
into the current fixed-truth scan.

## Files

- `core.py`: reusable physics layer—candidate conversion, anti-$k_T$
  re-clustering, ancestry reconstruction, truth matching, containment,
  contamination, jet shapes, Soft Drop, card overrides, bootstrapping, and
  provenance.
- `validate.py`: command-line driver that runs closure tests, scans GenFatJet
  radius, writes the JSON report, and creates explanatory plots and event
  displays.

## Quick start with the 10-event CMS sample

From the repository root:

```bash
cd model_building
source init.sh
```

Define the generated model directory:

```bash
# The Quickstart sample; for your own run_model output use models/<name>
MODEL_DIR="JetRadiusOptimization/validation_sample/s-channel_mmed-1000_Nc-2_Nf-2_scale-35.1539_mq-10_mpi-20_mrho-20_pvector-0.75_spectrum-cms_gq-0.25_gchi-0.5_rinv-0.3"
```

Run the fast GenFatJet-only validation and plotting example:

```bash
python3 JetRadiusOptimization/validate.py \
  --input "$MODEL_DIR/events.root" \
  --outdir "$MODEL_DIR/radius_diagnostics_10events" \
  --max-events 10 \
  --collections GenFatJet
```

A successful run ends with JSON containing:

```json
{
  "status": "passed"
}
```

Ten events are enough to verify the workflow and inspect event displays. They
are not enough to support a physics conclusion about the preferred radius.

## Full local validation

Omit `--collections` to test every genuinely clustered collection registered in
`core.py`:

```bash
python3 JetRadiusOptimization/validate.py \
  --input "$MODEL_DIR/events.root" \
  --outdir "$MODEL_DIR/radius_validation" \
  --max-events 200
```

The closure collections are:

```text
GenJet  GenFatJet  DarkPartonJet  DarkHadronJet  Jet  FatJet
```

The primary physics scan remains GenFatJet-focused; `--collections` selects
which stored collection closures are required.

## Radius configuration

The default physics grid is:

```text
0.2, 0.3, 0.4, ..., 1.6
```

The less crowded plotting subset is:

```text
0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6
```

The analysis default $R=0.8$ is always black. Other radii use gray, purple,
blue, orange, magenta, cyan-blue, and red rather than a green color scale.

To choose a smaller test grid:

```bash
python3 JetRadiusOptimization/validate.py \
  --input "$MODEL_DIR/events.root" \
  --outdir "$MODEL_DIR/radius_diagnostics_small_grid" \
  --max-events 10 \
  --collections GenFatJet \
  --radii 0.4 0.6 0.8 1.0 1.2 \
  --diagnostic-radii 0.4 0.6 0.8 1.0 1.2
```

Every diagnostic radius must also be present in the scan radius list.

## Generated plots

### `workflow_overview.png`

Explains the four collection definitions and the fixed-truth matching flow. Use
this first when presenting the analysis to a new collaborator.

### `softdrop_mass_by_radius.png`

One-dimensional overlays of GenFatJet Soft Drop mass for each displayed radius.
The definition is Soft Drop $\beta=0$, $z_{cut}=0.1$, with $R_0$ equal to the
jet radius. All overlays use the same common matched truth-jet cohort.

### `visible_descendant_containment_by_radius.png`

Count-weighted containment:

```text
number of the truth jet's visible descendants captured by GenFatJet
-------------------------------------------------------------------
total number of that truth jet's visible descendants
```

A value of 1 means every visible descendant is captured.

### `visible_descendant_pt_containment_by_radius.png`

$p_T$-weighted containment:

```text
captured visible-descendant pT
------------------------------
total visible-descendant pT
```

This distinguishes missing many soft descendants from missing one hard
descendant.

### `any_fully_contained_dark_hadron_acceptance_by_radius.png`

The requested "at least one fully contained dark hadron" result. A dark hadron
passes only when every one of its resolved visible final-state descendants is
an exact constituent of the matched GenFatJet. Dark hadrons with no resolved
visible descendants are not called fully contained.

The plot separates two denominators:

- **acceptance:** passing jets divided by all eligible fixed truth jets;
- **conditional efficiency:** passing jets divided by eligible truth jets that
  were successfully matched at that radius.

Dashed curves repeat both definitions for fixed truth jets containing more than
one dark hadron with visible descendants. This is the analogue of the
`strat1_5_any_full` selection in `StrategiesForTrainingData.py`, but is evaluated
here only as a physics-study diagnostic.

The updated two-panel version also shows the fraction with no fully contained
dark hadron, the fraction of matched jets containing at least one particle not
descended from any fixed dark-hadron truth jet, the mean non-DH constituent
count and $p_T$ fractions, and the mean number of such particles.

### `fully_contained_dark_hadron_count_by_radius.png`

One-dimensional, log-y overlays of the number of target dark hadrons that are
fully contained in each matched GenFatJet.

### `dark_hadron_contamination_counts_by_radius.png`

The left panel counts dark hadrons belonging to another fixed truth jet that
place at least one visible descendant in the target GenFatJet. The right panel
counts the subset whose complete visible-descendant group enters the target
jet.

### Additional focused diagnostic plots

- `non_dark_hadron_constituent_count_by_radius.png`: log-y overlays of the
  number of exact GenFatJet constituents with no ancestor among any initial
  `DarkHadronCandidate` in the event.
- `non_dark_hadron_constituent_fractions_by_radius.png`: log-y particle-count
  and $p_T$ contamination-fraction overlays.
- `dark_hadron_multiplicity_diagnostics_by_radius.png`: fully contained,
  not-fully-contained, and cross-truth contaminating dark-hadron counts.
- `containment_contamination_correlations_R08.png`: R=0.8 heatmaps connecting
  not-fully-contained DH multiplicity to non-DH particle count, and Soft Drop
  mass to non-DH $p_T$ fraction.
- `genfatjet_multiplicity_by_radius.png`: log-y event-level GenFatJet
  multiplicity overlays plus the mean, median, and event 16--84% interval as a
  function of R. It uses the same anti-$k_T$, GenCandidate, $p_T>15$ GeV jet
  definition as the focused radius scan.

For faster updates after the full validation has already been performed, see
`condor/README.md`. The focused DAG runs one 10,000-event job per radius and
merges their exact event-level arrays into these plots.

### `radius_performance_summary.png`

Shows containment and contamination together as functions of radius, plus truth
matching efficiency and split/merge behavior. Bands are 95% bootstrap intervals
on the mean. A large radius is not automatically preferred: increasing
containment can also increase contamination and merging.

The contamination definition is the fraction of clustered-jet constituent
$p_T$ that does not descend from any fixed truth dark-hadron jet.

### `collection_mass_comparison.png`

Compares stored `DarkHadronJet`, `DarkHadronStableJet`, and
`DarkHadronVisibleJet` masses with the common truth-matched GenFatJet $R=0.8$
cohort.

### `event_display_*.png`

Each image displays one large, zoomed matched jet with all requested radii
overlaid:

- each initial dark hadron is a labeled star with its own color;
- its visible descendants use the same parent color;
- round colored points are descendants captured at the reference $R=0.8$;
- colored crosses are descendants missed at the reference $R=0.8$;
- red open circles: non-truth constituents clustered into the jet;
- gray open diamonds: initial dark hadrons from other truth jets;
- red diamonds: other dark-hadron parents contaminating the $R=0.8$ jet;
- gold star: fixed DarkHadronJet truth axis;
- plus markers: re-clustered GenFatJet axes;
- colored circles: geometric radius guides, with the default $R=0.8$ in black.

The plotted membership is taken from the exact FastJet constituent list, not
from whether a marker happens to lie inside a circle. The table reports visible
descendant containment, target-DH full containment, cross-truth dark-hadron
contamination, total contamination $p_T$, and Soft Drop mass at every radius.

The selected jet is intentionally one whose containment changes substantially
between small and large radius. The display is a diagnostic failure case, not a
gallery of unusually clean jets.

## Metrics written to `validation_report.json`

For every radius, the report includes:

- visible-descendant containment by count and $p_T$;
- non-truth $p_T$ contamination;
- leading-dark-hadron containment;
- number of fully contained dark hadrons;
- fraction of target dark hadrons fully contained and whether at least one passes;
- number of contaminating dark hadrons from other truth jets, including how
  many of those are themselves fully contained;
- cross-truth dark-hadron contamination $p_T$ fraction;
- GenFatJet $p_T$ and mass;
- visible mass response;
- truth-association shared fraction;
- axis drift from fixed truth and from the $R=0.8$ jet;
- jet-shape radii `r50`, `r68`, `r90`, `r95`, and `r99`;
- girth, $p_TD$, and major/minor axes;
- truth match efficiency and split/merge counts.

Bootstrap resampling uses deterministic seeds. Change the number of resamples
with:

```bash
--bootstrap-resamples 1000
```

## Matching and same-event discipline

At each radius, visible stable `GenCandidate` particles are re-clustered with
anti-$k_T$ and a 15 GeV jet threshold. Truth jets and clustered jets are matched
one-to-one by maximizing shared visible-descendant $p_T$; matches with zero
shared $p_T$ are discarded.

The histogram overlays use the intersection of truth IDs successfully matched
at every displayed radius. Therefore every radius is compared on the same
physical truth-jet cohort. Do not build each radius from an independently capped
or independently selected sample; that creates a sample-dilution artifact.

The fixed-axis containment calculation is required to be monotonic in radius.
The physically re-clustered jet can move, split, or merge, so its observed
containment does not have to be monotonic.

## Validation gates

`validate.py` exits nonzero unless all requested gates pass:

1. **Card round-trip:** changing one module's `ParameterR` changes only that
   module and parses back to the requested value.
2. **Default-radius closure:** offline anti-$k_T$ re-clustering reproduces stored
   Delphes membership and four-vectors within the documented tolerances.
3. **Truth ancestry closure:** reconstructed visible-descendant sets agree with
   stored `DarkHadronVisibleJet` sets.
4. **Soft Drop closure:** the local Soft Drop calculation agrees with stored
   GenFatJet $R=0.8$ values.
5. **Primary scan invariants:** fixed truth IDs do not drift, and fixed-axis
   containment remains monotonic.
6. **Diagnostic plots:** a common matched truth cohort exists and all plots are
   produced.
7. **Display/metric agreement:** the exact constituent indices drawn in the
   event display equal those used in the numerical calculation.

Inspect the gate summary with:

```bash
python3 - "$MODEL_DIR/radius_diagnostics_10events/validation_report.json" <<'PY'
import json
import sys

report = json.load(open(sys.argv[1]))
print("status:", report["status"])
for name, passed in report.get("gates", {}).items():
    print(f"{name:32s} {passed}")
PY
```

## Reproducibility notes

The report records the input ROOT path, repository SHA, patched Delphes SHA,
`PYTHIA8MINOR` visible in the analysis environment, package versions, radius
grid, event count, and matching policy. The resolved generation cards and logs
remain beside `events.root` in the model directory.

Do not combine results made with different Pythia versions without labeling
them. The existing 10-event diagnostic sample was generated with the earlier
Pythia setup, while the fresh local installation currently reports Pythia 8.317.

## Common problems

### `events.root` is missing

Generate the signal separately with `run_model` or the top-level
`HTCondorSignalGeneration/` tutorial. This analysis never silently generates an
input sample.

### A required branch is missing

The input must come from the patched `DarkHadronJets` Delphes configuration. The
validator reports the exact missing branches and exits nonzero.

### The report status is `failed`

Read `gates` in `validation_report.json`. A failed closure is not a plotting
warning; it means the associated radius result should not be used.

### The plots look jagged or mostly empty

Increase the number of generated events. A 10-event run is intended only as a
workflow demonstration.

### `Predicted rinv = nan` appeared during generation

For the simplified CMS configuration, `rinv` is supplied directly rather than
predicted from a complete-model flavor count. That message does not by itself
invalidate the configured `rinv` or the jet-level containment plots.
