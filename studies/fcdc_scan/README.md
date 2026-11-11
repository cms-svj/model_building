# FCDC Λ / probVector scan

This study tests how the FCDC dark-shower model depends on the confinement scale Λ and on probVector, and whether the simplified model reproduces the complete one.

- **Complete model** (`spectrum = 'fcdc'`): Nf = 3 dark flavours, one of them stable. Mesons carrying the stable flavour are invisible.
- **Simplified model** (`spectrum = 'fcdcSimp'`): a single dark quark species; each meson goes invisible with probability rinv.

Each scan is studied in two ways: in reco-level distributions, and along the dark jet chain to see at which step of the event generation a parameter acts.

## Files

| File | What it does |
|---|---|
| `run_scan.sh` | generates the samples with `run_model` (cards → Pythia → Delphes → histograms) |
| `plot_scan.py` | plots and summary tables for each part of the study |

The model points come from [configs/model_fcdc_lambdascan.py](../../configs/model_fcdc_lambdascan.py). There is one config object per mπ/Λ and model, e.g. `r1p7_full`, `r1p7_simp`. The probVector variations reuse these objects with the `--pvector` override.

## Samples

All samples share:
- s-channel Z′ with m_Z′ = 1000 GeV, g_q = 0.25, g_χ = 1/√(Nc·Nf);
- Nc = Nf = 3, mπ = 6 GeV, η′ excluded;
- 10k events each.

The complete model uses Ns = 1 (rinv = 0.5); the simplified model uses input rinv = 0.5. Λ = mπ/(mπ/Λ), with mq and mρ from the Snowmass mass formulas. ρ → ππ is open below mπ/Λ ≈ 1.52 (2-body regime); above, ρ decays to π q q̄ or q q̄ (3-body regime).

| Scan | Points | Models | Samples |
|---|---|---|---|
| mπ/Λ, 2-body | 0.3, **0.6**, 1.0, 1.4 at probVector = 0.333 | complete + simplified | 8 |
| mπ/Λ, 3-body | 1.6, **1.7**, 2.0, 2.5 at probVector = 0.333 | complete + simplified | 8 |
| probVector | 0.1, 0.2, **0.333**, 0.5, 0.75 at mπ/Λ = 0.6 and 1.7 | complete + simplified | 16 |

Reference points are in bold. probVector = 0.333 at mπ/Λ = 0.6 and 1.7 is part of the mπ/Λ scan, so the probVector scan adds 16 new samples.

## Running

Run from the repository root on a machine where the repository environment (`init.sh`) works. Generation needs the installed Pythia and Delphes.

```bash
source init.sh
export FCDC_SAMPLES=/path/to/samples   # where the samples are written and read

# 1. generate all 32 samples (12 in parallel)
OUT=$FCDC_SAMPLES ./studies/fcdc_scan/run_scan.sh

# 2. plots and tables (about 2 minutes), written to plots/fcdc_scan/
python studies/fcdc_scan/plot_scan.py
```

`run_scan.sh` also sources `init.sh` itself, so it can be run from a fresh shell.

`run_scan.sh` options, set as environment variables:

| Variable | Default | Use |
|---|---|---|
| `OUT` | `models/fcdc` | output directory for the samples |
| `STEPS` | `all` | `"cards delphes hist"` reruns Delphes on the existing Pythia output; `hist` remakes only the histograms |
| `EVENTS` | `10000` | events per sample |
| `NPROC` | `12` | samples run in parallel |
| `ONLY` | all | regex on the sample name, e.g. `ONLY=r1p7` or `ONLY=pv0.75` |
| `DRY` | unset | `DRY=1` prints the commands without running them |

Logs go to `$OUT/logs/<sample>.log`. At the end, the script lists the samples whose `run_model` call failed and then exits with an error.

`plot_scan.py` options:

| Option | Default | Use |
|---|---|---|
| `--samples` | `$FCDC_SAMPLES`, else `models/fcdc` | where the samples are read from |
| `--outdir` | `plots/fcdc_scan` | where the output goes |
| `--only` | all four | make only some parts, e.g. `--only lambda_scan probvector_scan` |

After changing only the plotting, step 2 is enough: it reads the existing `Hists.pkl` files. `plot_scan.py` stops with a list of the missing samples if the scan is incomplete.

## Output

There is one directory per part of the study in `plots/fcdc_scan/`. Each plot is written as `.png` and `.pdf`.

**Reading the plots:**
- **Distribution plots:** unit-normalized overlays with a ratio to the reference sample. Leading two jets for jet variables.
- **Chain plots:** the mean of jet pT, constituent count, girth and Lund multiplicity (kT > 1 GeV) at each level, as a relative shift against the reference sample. Error bars are statistical only. The levels are:
  - **DP:** dark partons after the dark shower, anti-kT 0.8;
  - **DH:** dark hadrons;
  - **DH stable:** stable descendants of the DH jet constituents, including invisible particles;
  - **DH visible:** the same without the invisible particles;
  - **Reco:** Delphes jets.

  DH stable and DH visible are not reclustered and can extend beyond R = 0.8.

### `lambda_scan/`: Λ at fixed probVector = 0.333, complete model

| File | Content |
|---|---|
| `distributions_2body`, `distributions_3body` | substructure (girth, major/minor axis, ECF₂, Lund multiplicity) and event-level variables (MET, mT, m_jj, jet pT, Δφ) for each mπ/Λ of the regime |
| `chain_2body`, `chain_3body` | mean shift vs the reference mπ/Λ along the dark jet chain |
| `means_vs_mpi_over_lambda` | mean Lund multiplicity, girth, mT and per-jet rinv vs mπ/Λ across both regimes, complete and simplified |
| `rinv_vs_mpi_over_lambda` | mean per-jet rinv vs mπ/Λ, generated vs prediction (0.5 in 2-body, effective-rinv formula `fcdc_rinv_3body[_simp]` in 3-body), both models |
| `chain_shifts.csv` | every point of the chain plots |

### `probvector_scan/`: probVector at mπ/Λ = 0.6 (2-body) and 1.7 (3-body), complete model

| File | Content |
|---|---|
| `distributions_2body`, `distributions_3body` | substructure and event-level variables for each probVector |
| `chain_2body`, `chain_3body` | mean shift vs probVector = 0.333 along the dark jet chain |
| `collections_2body`, `collections_3body` | distributions of constituent count, girth and Lund multiplicity for dark partons, dark hadrons, DH visible and reco jets |
| `rinv_vs_pvector` | mean per-jet rinv vs probVector, generated vs prediction (0.5 in 2-body, effective-rinv formula `fcdc_rinv_3body[_simp]` in 3-body), both models |
| `chain_shifts.csv` | every point of the chain plots |

### `complete_vs_simplified/`: simplified vs complete model at equal parameters

| File | Content |
|---|---|
| `distributions` | girth, Lund multiplicity, MET, mT and per-jet rinv at mπ/Λ = 0.6 and 1.7 |
| `chain` | simplified ÷ complete along the dark jet chain at every mπ/Λ (2-body blue, 3-body orange) |
| `rinv_per_jet` | per-jet rinv distribution at mπ/Λ = 0.6 and 1.7, with the fraction of jets without invisible energy |
| `pt_offset_3body` | left: fraction of ρ decays to π + q q̄ vs mπ/Λ in 3-body, generated and expected from the model inputs; right: simplified vs complete jet pT at DH-visible and reco level vs mπ/Λ, with the prediction (1 − rinv_eff,simp)/(1 − rinv_eff,complete) − 1 |
| `chain_shifts.csv` | every point of the chain plot |

### `lund_kt_threshold/`: Lund multiplicity at higher kT thresholds

| File | Content |
|---|---|
| `shift_vs_kt_threshold` | shift of the mean primary Lund multiplicity vs the kT threshold (1, 2, 5, 10 GeV), at dark-parton and reco level, for the Λ and probVector scans in both regimes |
| `shifts.csv` | every point of the plot |

### Tables

The shift tables (`chain_shifts.csv`, `shifts.csv`) have one row per point. The columns are:
- the sample: `model`, `mpi_over_lambda`, `pvector`;
- the reference it is compared to: `ref_model`, `ref_mpi_over_lambda`, `ref_pvector`;
- `variable` and `level`;
- `shift_percent` and `stat_error_percent`.

`sample_summary.csv` has one row per sample:
- the model parameters;
- the means of the reco-level substructure and event variables;
- the mean Lund multiplicity for each kT threshold at dark-parton and reco level;
- the per-jet rinv (generated and predicted);
- the dark-hadron counters (ρ → π q q̄ and ρ → q q̄ decays per event, vector fraction);
- the ratio of mean reco jet pT to DH-visible jet pT.

## Notes

- **Means of counting variables** (Lund multiplicity, hadron counts) use the integer value of each unit bin [n, n+1), not the bin center.
- **Vector fraction:** the `n_pion` histogram stops at 20 and overflows at low probVector, so the vector fraction in `sample_summary.csv` is only reliable for probVector ≥ 0.333.
- **Pythia version:** the samples were generated with Pythia 8.310.
- **Plot style:** comes from `common.set_plot_style`.
