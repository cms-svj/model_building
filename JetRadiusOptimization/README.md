# Jet radius optimization

Which jet radius best reconstructs the dark shower, at different stages of the
shower? The study looks at jets of **dark hadrons** at radius R_DH (the truth
jets) and jets of **visible particles** at radius R. For each R_DH it asks which
R best reconstructs the dark-hadron jets. It is built on the dark-hadron
ancestry in `darkancestry.py`.

## Definitions: two stages of the shower

The dark shower is probed at two stages, each clustered into jets with anti-kT
at its own radius (pT > 15 GeV for both, as in the Delphes card):

| Radius | Stage | Clustered from | Standard Delphes jets on this axis |
|---|---|---|---|
| **R_DH** | dark hadrons, after the dark shower has hadronized and before they decay | `DarkHadronCandidate` (the dark hadrons' four-momenta) | R_DH = 0.8: `DarkHadronJet` |
| **R** | visible particles, after the dark hadrons decay to SM particles and hadronize | `GenCandidate` (stable visible particles; neutrinos and stable dark hadrons removed) | R = 0.4: `GenJet`; R = 0.8: `GenFatJet` |

- **R_DH** is the radius of the jets **built from the dark hadrons**; it sets
  how the dark-hadron stage of the shower is grouped into objects (small R_DH:
  several jets per dark quark; large R_DH: about one jet per dark quark). These
  dark-hadron jets are the truth jets.
- **R** is the radius of the **generator-level visible jets**, i.e. a
  `GenJet`/`GenFatJet` with a variable radius. This is the stage an experiment
  reconstructs.
- The question at each R_DH: which visible-jet radius R best reconstructs the
  dark-hadron jets of radius R_DH? Scanning both gives the mapping between the
  two stages.

## What it does

For each (R_DH, R):

1. Cluster the dark hadrons at R_DH (truth jets) and the visible particles at R
   (visible jets), as defined above. At R_DH = R = 0.8 these are exactly the
   Delphes `DarkHadronJet` and `GenFatJet`.
2. Match truth jets one-to-one to visible jets by shared visible pT.
3. Split the matched jet's pT by where it came from, using `darkancestry`:

| | pT fraction from |
|---|---|
| `f_full` | dark hadrons **completely** inside the jet |
| `f_partial` | dark hadrons only **partly** inside the jet |
| `f_nodh` | particles with **no** dark-hadron ancestor (ISR, UE, ...) |

These add up to 1.

## Scores

All are fractions of pT. For a truth jet T matched to a visible jet J, with
pT(S) the scalar pT sum over a set of particles S, V_d the visible descendants
of dark hadron d, and V_T the union of V_d over the dark hadrons of T:

- d is **complete** in J if pT(V_d ∩ J) ≥ 0.95 · pT(V_d) (`--min-containment`,
  default 0.95; 1 requires every visible descendant)

| Score | Definition | Meaning |
|---|---|---|
| **acceptance** | Σ_{d ∈ T, complete} pT(V_d ∩ J) / pT(V_T) | fraction of the truth jet's visible pT that J holds through complete dark hadrons (penalizes R too small) |
| **purity** | Σ_{d complete in J} pT(V_d ∩ J) / pT(J) = `f_full` | fraction of J's pT from complete dark hadrons, of any truth jet (penalizes UE/ISR and partial dark hadrons) |
| **quality** | acceptance × purity | |
| **IoU** | pT(V_T ∩ J) / pT(V_T ∪ J) | overlap of J with the truth jet's visible particles; does not depend on the containment threshold |

An unmatched truth jet scores 0. Truth jets whose dark hadrons have no visible
descendants are skipped. The output also keeps the counts `n_dh` (dark hadrons
of T with visible descendants) and `n_full` (complete ones). Results are shown
both inclusively and in bins of the truth jet's visible pT, pT(V_T).

## Files and running

- `radius.py` (library): clustering, matching, metrics. Run directly, it checks
  our clustering against the stored `GenFatJet` and `DarkHadronJet`.
- `scan.py`: runs the scan for one `run_model` directory (default R = 0.2–2.0
  and R_DH = 0.8, 1.0; 0.8 and 1.0 are always included in both grids; set
  `--truth-radii` for the full 2D grid) and writes `JetRadius.pkl` next to
  `events.root`. Model parameters come from the directory's `config.py` via
  `svjHelper`.
- `plot.py`: scores vs R at one R_DH (inclusive and in bins of truth visible
  pT), 2D maps of quality and IoU over (R, R_DH) with the best R per R_DH
  marked, and the best-R points of all models overlaid on the same axes; PDF
  by default. Takes a model-directory pattern through `common.resolve_models`,
  so EOS (`/store/...`) works too.

From the repository root, after `source init.sh`:

```bash
# check our clustering against the stored GenFatJet and DarkHadronJet
python3 JetRadiusOptimization/radius.py -i models/<model>/events.root

# full 2D scan for one model directory made by run_model
python3 JetRadiusOptimization/scan.py -d models/<model> --truth-radii $(seq 0.2 0.2 2.0)

# plots for all scanned models, one curve per value of a model parameter
python3 JetRadiusOptimization/plot.py -m "models/*" -x rinv
```

Use `--help` on any script for all options.

## Output

- `JetRadius.pkl` (from `scan.py`): the model metadata, the radius grids, and
  one row per (event, truth jet, R_DH, R) with the scores, `f_full`,
  `f_partial`, `f_nodh`, `n_dh`, `n_full`, and the truth and matched-jet pT.
- `JetRadius_plots/` (from `plot.py`, `--dir`):

| File | Shows |
|---|---|
| `<score>.pdf` | each score vs R at one R_DH (`--truth-radius`, default 0.8), one curve per model |
| `<score>_ptbins_<model>.pdf` | each score vs R in bins of the truth jet's visible pT (`--pt-bins`) |
| `best_radius_vs_pt_<score>.pdf` | best R per truth-pT bin |
| `quality_2d_<model>.pdf`, `iou_2d_<model>.pdf` | the mean score over (R, R_DH), with the best R per R_DH marked |
| `best_radius_quality.pdf`, `best_radius_iou.pdf` | the best-R points of all models overlaid |

## Validation (CMS benchmark, 10k events)

- Our dark-hadron clustering at R_DH = 0.8 reproduces all 40,268 stored
  `DarkHadronJet`s exactly. `darkancestry` check: 40,268 jets, 0 mismatches.
- Our visible clustering at R = 0.8 reproduces the stored `GenFatJet`s except
  21 of 66,473 jets, all from numerical near-ties in anti-kT: Delphes clustered
  the unrounded momenta, while we cluster the float32 values stored in the
  file. Nudging the inputs by 1e-7 to 1e-6 (relative) changes our own
  clustering in all of these events, versus under 1% of other events. The
  check reports these as `unstable_jets: 21` and `mismatched_jets: 0`.

## Results (preliminary): CMS benchmark, 10k events

Model: s-channel, m_med = 1 TeV, N_c = N_f = 2, Λ = 35.15 GeV,
m_π = m_ρ = 20 GeV, p_vector = 0.75, r_inv = 0.3 (`configs/model_cms.py`).
Containment threshold 0.95.

Best visible-jet radius R for each truth radius R_DH (mean over all truth
jets; `best_radius_quality.pdf`, `best_radius_iou.pdf`):

| R_DH | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 | 1.2 | 1.4 | 1.6 | 1.8 | 2.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| best R by quality | 1.2 | 1.3 | 1.5 | 1.5 | 1.6 | 1.6 | 1.7 | 1.7 | 1.8 | 1.8 |
| best R by IoU | 0.3 | 0.5 | 0.6 | 0.7 | 0.9 | 1.0 | 1.1 | 1.2 | 1.3 | 1.3 |

In bins of the truth jet's visible pT, at R_DH = 0.8 (best R, with the score
there in parentheses; `<score>_ptbins_<model>.pdf`):

| truth visible pT | truth jets | quality | IoU | acceptance | purity |
|---|---|---|---|---|---|
| 0–50 GeV | 13,160 | 1.2 (0.03) | 0.6 (0.41) | 1.6 (0.06) | 2.0 (0.07) |
| 50–100 GeV | 5,943 | 1.2 (0.22) | 0.7 (0.59) | 1.7 (0.37) | 2.0 (0.34) |
| 100–200 GeV | 6,865 | 1.6 (0.54) | 0.8 (0.77) | 1.9 (0.80) | 1.5 (0.63) |
| > 200 GeV | 10,116 | 1.6 (0.78) | 0.9 (0.90) | 2.0 (0.96) | 1.2 (0.83) |

- **Quality** prefers R well above R_DH: to hold 95% of each dark hadron's
  visible pT, the visible jet has to be wider than the dark-hadron jet, since
  the decay products spread around each dark hadron's direction.
- **IoU** prefers R slightly below R_DH, approaching R_DH for hard jets,
  because it also charges for UE/ISR pT in the jet.
- **Purity** turns over only for hard jets (best R ≈ 1.2–1.5 above 100 GeV,
  where extra UE/ISR outweighs the gain in complete dark hadrons).
  **Acceptance** always rises with R.
- 36% of truth jets have less than 50 GeV of visible pT; their quality stays
  near 0.03 at every radius, so they dilute the inclusive means. The pT-binned
  results are the more informative ones.
- The two headline scores disagree systematically, so the answer depends on
  what the jet is meant to capture (whole dark hadrons vs a clean match to the
  dark-hadron jet). Both scores also rise with R_DH, so the result is a mapping
  from R_DH to R rather than a single optimal pair.
- Requiring every visible descendant (`--min-containment 1`) makes almost every
  dark hadron "partial" (one soft particle outside is enough), so quality only
  grows with R; hence the 0.95 default.

## Known limitations / follow-ups

- **Dark-parton stage:** the same scan one step earlier (dark partons → dark
  hadrons), using `darkancestry` with partons as targets, would complete the
  parton → hadron → visible chain.
- **Generator level only:** visible jets are clustered from `GenCandidate`.
  Reconstructed (particle-flow) jets need PF → generator-particle links.
- **One model point:** the results above are for the CMS benchmark
  (m_med = 1 TeV, m_dark = 20 GeV, r_inv = 0.3). The scan should be run over
  the parameter grid.
- **Choice of score:** quality and IoU point to different radii; which one to
  use is open.
