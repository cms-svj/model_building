# Dark Sector Observable Automation

This folder is part of the [`model_building`](https://github.com/cms-svj/model_building) repository. It contains a workflow for adding new jet-level observables to semi-visible jet (SVJ) / dark sector samples. The goal is to automate the process of running dark sector model configurations, calculating additional observables, and saving those observables into a ROOT friend tree for later analysis.

Anyone using this code should first follow the installation instructions and environment setup described in the top-level README of this repository. The rest of `model_building` provides the baseline dark sector model-building framework, the model configuration files, and the tools needed to run the event generation workflow.

---

## Purpose

This folder extends the original dark-sector model-building workflow in two main ways.

First, it provides a more user-friendly way to automate the generation of dark-sector models across different parameter choices. Instead of manually editing configuration files for each model, the workflow allows important dark-sector parameters to be passed in a more systematic way, making it easier to scan parameter space and generate the corresponding model directories.

Second, it adds new observables to the generated ROOT output through a ROOT friend tree. The friend tree contains additional FatJet-level observables, such as jet-shape and substructure variables, while keeping them separate from the original `events.root` Delphes file. This allows the original detector-level output to remain unchanged while still making the new variables available for analysis.

Together, these additions make the workflow more practical for studying the phenomenology of a mirror dark-QCD model, from automated model generation to the analysis of new dark-sector observables.

The friend tree currently includes the following observables:

```python
jets_like = {
    "N2":         FatJet_N2,
    "N3":         FatJet_N3,
    "LundX":      FatJet_LundX,
    "LundY":      FatJet_LundY
}
```

These variables are stored as FatJet-level branches inside a friend ROOT tree called `DelphesFriend`.

The output friend file is saved with a name like:

```text
events_friend_<model_tag>_j<cluster>.<proc>.root
```

Here, `<model_tag>` refers to the full model directory name generated from the dark-sector parameters used in that run. For example, it can include the mediator mass, number of colors/flavors, dark pion mass, dark rho mass, spectrum choice, and invisible fraction. The `<cluster>` and `<proc>` values come from the HTCondor job identifiers, where `<cluster>` labels the submitted job cluster and `<proc>` labels the individual process within that cluster. In practice, `<proc>` is useful for distinguishing different jobs submitted together, such as jobs scanning over different values of a dark sector parameter. Together, these labels make the output file names unique for each submitted job.

---

## Folder layout

```text
model_building/
├── configs/
│   ├── master_snowmass.py            # scan-ready Snowmass-based config (added here)
│   ├── master_cms.py                 # scan-ready CMS-based config (added here)
│   └── ...                           # original model configs
└── automation_of_dark_sector_variables/
    ├── condor/
    │   ├── automating_jobs.jdl       # HTCondor submission file
    │   ├── pack.sh                   # packages the repo and uploads it to EOS
    │   ├── run_job.sh                # script that runs on the worker node
    │   └── logs/
    │       └── .gitkeep              # keeps the folder in git; Condor logs go here
    ├── examples/
    │   └── condor_3086704_1_rinv-0.5.out   # sample output of a successful job
    ├── scripts/
    │   ├── add_new_DS_observables.ipynb
    │   ├── add_new_DS_observables.py
    │   ├── merge_trees.C
    │   └── common.py                 # modified copy of the top-level common.py (see Notes)
    ├── .gitignore
    └── README.md
```

---

## How to run

Run these commands from the `condor/` folder of this directory.

```bash
# 1. Get a grid proxy (needed to read from and write to EOS)
voms-proxy-init -voms cms -valid 192:00

# 2. Go to the condor folder
cd automation_of_dark_sector_variables/condor

# 3. Package the repository and upload it to EOS
bash pack.sh

# 4. Submit the jobs
condor_submit automating_jobs.jdl

# 5. Monitor them
condor_q
```

Before submitting, open `automating_jobs.jdl` and check:

- **Model**: the second word of the `arguments` line. Options are `cms`, `snowmass`, `both`, `master_snowmass` and `master_cms`. The scan variable `rinv` is only used by the two `master_*` models.
- **Number of events**: the third word of the `arguments` line.
- **Scan values**: the `queue rinv in ...` line at the end of the file.
- **`PROJECT_NAME`** in the `environment` line must match the name of the folder you cloned the repository into (`model_building` by default).

The EOS user name is taken automatically from the user who submits the jobs. Note that `pack.sh` overwrites the bundle on EOS, so wait for any jobs that are still running with an older version of the code before packing again.

---

## Main notebook

The main analysis notebook is:

```text
scripts/add_new_DS_observables.ipynb
```

This notebook is used to develop and test the calculation of new dark sector observables. The purpose of using the notebook is to work interactively with multiple variables at the same time without having to rerun the entire full production workflow every time a change is made. Run it from inside the `scripts/` folder, so that it uses the `common.py` located there.

Once the notebook is ready to be used in batch mode, it is converted into a Python script. From the `scripts/` folder:

```bash
jupyter nbconvert --to script add_new_DS_observables.ipynb --output add_new_DS_observables
```

This creates:

```text
scripts/add_new_DS_observables.py
```

The Python script is then used by the batch workflow.


To run the notebook interactively, start Jupyter after `source init.sh` and select the kernel that uses the `mbenv` virtual environment created by `./install.sh`. Install this kernel once with:

```bash
python3 -m ipykernel install --user --name mbenv --display-name "Python (mbenv)"
```

Then, in the notebook, choose Kernel → Change kernel → Python (mbenv). Without this step, the notebook may fail with `ModuleNotFoundError: No module named 'magiconfig'`.

---

## Local packaging workflow

The file `condor/pack.sh` packages the local repository and uploads it to EOS so that HTCondor jobs can access the same code.

The script does not need to be edited. It finds the repository from its own location, so it works wherever the repository was cloned:

- The repository is the folder two levels above `pack.sh`, and its name is used for the bundle name.
- The bundle is created in the folder that contains the repository (not inside it, so it is never added to git).
- The EOS user defaults to your user name. To use a different one: `EOS_USER=name bash pack.sh`.

The script creates a compressed tarball:

```text
<repo_folder>_bundle.tgz
```

while excluding files and directories that should not be sent to the worker node:

```text
.git/
__pycache__/
.ipynb_checkpoints/
logs/
jobs/
models/
automation_of_dark_sector_variables/condor/logs/
*.root
*.pdf
*.png
```

The tarball is then uploaded to EOS under:

```text
/store/user/<user_name>/<repo_folder>/inputs
```

This makes the code bundle available to the Condor jobs.

---

## Current scan setup (automation)

The current workflow is designed to scan over values of:

```text
rinv
```

The scan is controlled through the master configuration files in the top-level `configs/` folder of the repository:

```text
configs/master_snowmass.py      (Snowmass model based)
configs/master_cms.py           (CMS-model based)
```

These configuration files are read by `run_job.sh` during the Condor workflow. Their purpose is to provide configurable model setups where selected dark sector parameters can be changed without manually editing multiple model configuration files.

The `master_snowmass.py` file is based on the Snowmass model, while `master_cms.py` is the CMS-like master configuration. Both allow the scan variable to be changed directly from the Condor workflow.

Both files are based respectively on existing configurations in the original repository:

```text
configs/model_snowmass_cmslike.py
configs/model_cms.py
```

In the JDL file, different values of `rinv` are passed as arguments to `run_job.sh`. The shell script then exports the corresponding value (as the environment variable `RINV`) so that the master configuration can use it when running the model.

For example, the JDL file can submit jobs for a single value such as:

```text
0.5
```

or an expanded scan such as:

```text
0.1, 0.3, 0.5, 0.7, 0.9
```

Each value of `rinv` is passed to the shell script, which runs the corresponding model and produces output files for that scan point.

---

## HTCondor workflow

The file `condor/run_job.sh` is the main shell script that runs on the worker node. It handles the full job workflow.

In broad terms, `run_job.sh` does the following:

1. Sets up the required environment.
2. Stages in the packaged code from EOS.
3. Runs the requested dark-sector model configuration.
4. Runs the observable calculation script `scripts/add_new_DS_observables.py`.
5. The `add_new_DS_observables.py` script creates the ROOT friend tree containing the new observables.
6. Using the `scripts/merge_trees.C` macro, it merges the main ROOT tree and the friend ROOT tree into one single tree.
7. Saves the output files back to EOS.

The file `condor/automating_jobs.jdl` is the HTCondor submission file. It controls how jobs are submitted to the cluster. Its arguments are, in order: the job tag, the model, the number of events, the random seed, and the value of `rinv`.

In the JDL file, the variable to scan is chosen. So far, the main scan variable is `rinv`. The JDL file sends the chosen `rinv` values to `run_job.sh` as arguments. The shell script then uses those values to run the model configurations described in the original model-building repository.

The JDL file also controls the creation of Condor log files. The output, error, and log files are written to the `condor/logs/` directory, with names that include the Condor cluster/process ID and the value of `rinv`.

Example log naming pattern:

```text
logs/condor_$(Cluster)_$(Process)_rinv-$(rinv).out
logs/condor_$(Cluster)_$(Process)_rinv-$(rinv).err
logs/condor_$(Cluster)_rinv-$(rinv).log
```

The `logs/` folder must exist when the jobs are submitted, otherwise Condor puts the jobs on hold. It is already included in the repository.

---

## Output files

At the end of each run, EOS will contain a list of job output directories. For the `master_snowmass` and `master_cms` models they look like:

```text
/store/user/<user_name>/model_building/outputs/j<cluster>.<proc>-<model>-rinv=<value>/<model_tag>
```

where `<model_tag>`, `<cluster>`, and `<proc>` are the same as described above in the "Purpose" section, `<model>` is `master_snowmass` or `master_cms`, and `<value>` is the value of `rinv` used for that specific run. For the other models the job directory is simply `j<cluster>.<proc>`.

In the `<model_tag>` directory, the following files can be found:

```text
events.root
pythia_card.txt
delphes_card.txt
config.py
events_friend_<model_tag>_j<cluster>.<proc>.root
merged_tree_<model_tag>_j<cluster>.<proc>.root
```

The first four files listed come directly from the model generation stage. The last two are, respectively, the friend ROOT file containing only the new observables and the merged ROOT file containing variables from both ROOT trees.

The `events.root` file contains the Delphes-level event information and the original reconstructed object kinematics. This is the main input used by the observable workflow to compute the additional FatJet variables.

The `config.py` file records the dark-sector model parameters used for that specific run, such as the mediator mass, dark pion mass, dark rho mass, number of colors/flavors, and invisible fraction.

The `pythia_card.txt` file stores the Pythia settings used during event generation and hadronization, including the dark-sector shower and decay settings.

The `delphes_card.txt` file stores the Delphes detector-simulation settings, including the detector response, object reconstruction, and jet definitions used to produce the final ROOT output.

The `merged_tree_<model_tag>_j<cluster>.<proc>.root` file contains the original Delphes tree information together with the additional FatJet observables from the friend tree. The added observables are written with an `F_` prefix to avoid name collisions with branches already present in the original Delphes tree.

In the job directory (one level above `<model_tag>`), the logs of the job are also saved: `run_j<cluster>.<proc>.log` (output of the observable calculation) and `gen_j<cluster>.<proc>_*.log` (output of the event generation).

All these output files are saved to EOS and are available for later analysis.

---

## Notes

Large generated files are intentionally excluded from version control through the `.gitignore` in this folder. This includes:

```text
*.root
*.hepmc
*.lhe
*.png
*.pdf
*.txt
*.pyc
condor/logs/*      (except the placeholder .gitkeep)
work_*/
.root_hist
.wget-hsts
models/
install/
jobs/
cards/
test_downloads/
__pycache__/
.ipynb_checkpoints/
```

This folder should mainly contain source code, notebooks, shell scripts, and HTCondor submission files. The master configuration files live in the top-level `configs/` folder.

**Modified `common.py`.** The `scripts/` folder contains a modified version of `common.py` from the original model-building repository. The `fix_delphes_mass_units` function was updated to check whether each particle collection and its `Mass` field exist before applying the unit conversion, since some ROOT files do not contain all expected collections, such as `GenStableCandidate`. Similarly, `init_constituents` now verifies that each jet collection, associated constituent collection, and `Constituents` field are available before processing them, allowing the code to safely skip missing collections instead of failing.

When Python runs a script it looks for imports in the script's own folder first, so `from common import load_sample` in `add_new_DS_observables.py` uses this modified copy and not the top-level `common.py`. Other top-level modules (for example `svjHelper`) are still found through the repository root, which `run_job.sh` adds to `PYTHONPATH`. This copy is based on the top-level `common.py` at the time of writing. If the top-level file changes in ways that this workflow needs, the changes have to be copied over.

**Sample output.** The folder `examples/` contains a sample file named `condor_3086704_1_rinv-0.5.out`. This is an output file from Condor for a job submission of only 20 events, so it should not be used for analysis. It should only be used as a sample of what a successful submission looks like.