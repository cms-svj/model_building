# HTCondor jet-radius analysis

This directory submits one analysis job for the completed 10,000-event CMS
sample. It is separate from signal generation. The worker copies `events.root`
from EOS, runs the full R=0.2--1.6 scan, and copies plots and the validation JSON
back to EOS. It transfers a plain runtime directory and does not use a tarball.

Submit from the repository root. The first command builds
`condor/model_building_transfer/` (`init.sh`, `common.py` and the installed
python venv), which every worker receives. Rerun it whenever those change.

```bash
source init.sh
JetRadiusOptimization/condor/make_transfer_dir.sh
mkdir -p JetRadiusOptimization/condor/logs
condor_submit JetRadiusOptimization/condor/analysis.sub
```

By default the input and output are your own EOS area
(`/store/user/$USER/DarkSectorStudies/...`, see the macros at the top of each
`.sub` file). To point at another sample, override them at submit time:

```bash
condor_submit \
  -append 'INPUT_URL = root://cmseos.fnal.gov//store/user/<user>/<sample>/events.root' \
  -append 'EOS_OUT = /store/user/$ENV(USER)/<analysis output dir>' \
  JetRadiusOptimization/condor/analysis.sub
```

Monitor with the cluster ID printed by `condor_submit`:

```bash
condor_q CLUSTER.PROCESS
condor_tail -name lpcschedd4.fnal.gov CLUSTER.PROCESS
```

The output is complete only when `_SUCCESS.json` exists under:

```text
/store/user/$USER/DarkSectorStudies/Analysis/cms_10k_radius_study/
  point-00000_mdark-20_shard-00/jet_radius_validation/
```

## Faster focused diagnostics: one job per radius

Use this path when the full closure suite and event displays have already been
run and only containment, non-dark-hadron contamination, and Soft Drop
diagnostic plots need updating. It starts 15 independent jobs: one for each
radius from 0.2 through 1.6. Each job processes the same 10,000 events at only
its assigned radius and writes mergeable event-level arrays. A dependent merge
job runs automatically after all radii finish.

No signal is generated and no training data is made.

From the repository root (after `make_transfer_dir.sh`, as above):

```bash
source init.sh
mkdir -p JetRadiusOptimization/condor/logs

condor_submit_dag JetRadiusOptimization/condor/focused_by_radius.dag
```

The DAG contains:

```text
15 radius jobs (R=0.2, 0.3, ..., 1.6)
                    |
                    v
             1 merge/plot job
```

Monitor all jobs:

```bash
condor_q "$USER" -batch
```

The final output is ready only when this marker exists:

```bash
xrdfs root://cmseos.fnal.gov stat \
  /store/user/$USER/DarkSectorStudies/Analysis/cms_10k_radius_study/\
point-00000_mdark-20_shard-00/jet_radius_validation/\
_FOCUSED_DIAGNOSTICS_SUCCESS.json
```

The radius jobs use 15 single-CPU slots rather than attempting to stretch one
Python process across machines. This reduces wall time when LPC slots are
available, at the cost of copying the input ROOT file once per radius.
