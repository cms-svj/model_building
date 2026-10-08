#!/usr/bin/env bash

export PYTHONUNBUFFERED=1
export TERM=dumb   # silence tput in non-interactive jobs (use simplest terminal)

JOBTAG="${1:?need <jobtag>}"
MODEL="${2:?need <model>: cms|snowmass|both|master_snowmass|master_cms}"
NEV="${3:?need <events>}"
SEED="${4:?need <seed>}"



RINV="${5:-}"

if [[ "${MODEL}" == "master_snowmass" && -z "${RINV}" ]]; then
  echo "[FATAL] For MODEL=master_snowmass you must pass a 5th argument: RINV"
  echo "        Example args: j123.0 master_snowmass 20000 100001 0.3"
  exit 2
fi

if [[ "${MODEL}" == "master_cms" && -z "${RINV}" ]]; then
  echo "[FATAL] For MODEL=master_cms you must pass a 5th argument: RINV"
  echo "        Example args: j123.0 master_cms 20000 100001 0.3"
  exit 2
fi

if [[ -n "${RINV}" ]]; then
  export RINV
fi



# Project name (default = model_building, but can be overridden with env PROJECT_NAME)
PROJECT="${PROJECT_NAME:-model_building}"
EOS_USER="${EOS_USER:?need EOS_USER (set in the .jdl environment line)}"
AUTO_DIR="automation_of_dark_sector_variables"


EOS_AREA="${EOS_AREA:-$PROJECT}"





echo "===== NODE INFO ====="
date; hostname
echo "PWD: $(pwd)"
# echo "ARGS: JOBTAG=${JOBTAG} MODEL=${MODEL} NEV=${NEV} SEED=${SEED}"    
echo "ARGS: JOBTAG=${JOBTAG} MODEL=${MODEL} NEV=${NEV} SEED=${SEED} RINV=${RINV:-<unset>}"
echo "ENV : PROJECT=${PROJECT}"
echo "====================="

# Working area for this job
WORKDIR="$PWD/work_${JOBTAG}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"







# EOS_BASE_RSE="root://cmseos.fnal.gov//store/user/${EOS_USER}/${PROJECT}"
EOS_BASE_RSE="root://cmseos.fnal.gov//store/user/${EOS_USER}/${EOS_AREA}"












echo "ENV : PROJECT=${PROJECT} EOS_USER=${EOS_USER}"

echo "[stage-in] fetch code bundle ..."
xrdcp -f "${EOS_BASE_RSE}/inputs/${PROJECT}_bundle.tgz" bundle.tgz

echo "[untar] expanding ..."
tar xzf bundle.tgz

cd "${PROJECT}"

# [B] verify what’s in the repo + that configs exist
echo "[debug] repo tree (top level)"
ls -1
echo "[debug] configs available:"
ls -1 configs || true

SNOW_CFG="configs/model_snowmass_cmslike.py"
CMS_CFG="configs/model_cms.py"
MASTER_CFG="configs/master_snowmass.py"
MASTER_CMS_CFG="configs/master_cms.py"
[ -f "$SNOW_CFG" ] || echo "[ERROR] Missing $SNOW_CFG in bundle!"
[ -f "$CMS_CFG" ]   || echo "[ERROR] Missing $CMS_CFG in bundle!"
[ -f "$MASTER_CFG" ] || echo "[ERROR] Missing $MASTER_CFG in bundle!"
[ -f "$MASTER_CMS_CFG" ] || echo "[ERROR] Missing $MASTER_CMS_CFG in bundle!"



# Comands from original notebook to set up environment
./install.sh
source init.sh










echo "[debug] before fix: python3=$(which python3) VIRTUAL_ENV=${VIRTUAL_ENV:-<unset>}"
python3 -c "import coffea; print('[debug] coffea', coffea.__version__, coffea.__file__)"

# The bundled mbenv was built under a different absolute path; make sure the unpacked copy is the one used
MBENV="$PWD/install/python_packages/mbenv"
if [ -x "${MBENV}/bin/python3" ]; then
  export VIRTUAL_ENV="${MBENV}"
  export PATH="${MBENV}/bin:${PATH}"
  hash -r
fi

echo "[debug] after fix: python3=$(which python3)"
python3 -c "import coffea; print('[debug] coffea', coffea.__version__, coffea.__file__)"



















echo "[check] python and key packages"
python3 -V
python3 - <<'PY'
import sys
print("python ok:", sys.version.split()[0])
for m in ["uproot","awkward","numpy"]:
    __import__(m)
    print("import ok:", m)
PY

chmod -x ./run_model 2>/dev/null || true
chmod +x ./run_model || true




# Fail fast if the config needed for the chosen model is missing
case "$MODEL" in
  cms)
    [ -f "$CMS_CFG" ] || { echo "[FATAL] $CMS_CFG not found inside the tarball; rebuild bundle to include it."; exit 3; }
    ;;
  snowmass)
    [ -f "$SNOW_CFG" ] || { echo "[FATAL] $SNOW_CFG not found inside the tarball; rebuild bundle to include it."; exit 3; }
    ;;
    master_snowmass)
    [ -f "$MASTER_CFG" ] || { echo "[FATAL] $MASTER_CFG not found inside the tarball; rebuild bundle to include it."; exit 3; }
    ;;
    master_cms)
    [ -f "$MASTER_CMS_CFG" ] || { echo "[FATAL] $MASTER_CMS_CFG not found inside the tarball; rebuild bundle to include it."; exit 3; }
    ;;
  both)
    for cfg in "$CMS_CFG" "$SNOW_CFG"; do
      [ -f "$cfg" ] || { echo "[FATAL] $cfg not found inside the tarball; rebuild bundle to include it."; exit 3; }
    done
    ;;
  *)
    echo "[FATAL] MODEL must be one of: cms | snowmass | both | master_snowmass | master_cms"; exit 2
    ;;
esac


# Time total time og generating events in the chosen model
echo "[timing] Starting event generation..."
echo ""
GEN_START=$(date +%s)




echo "[debug] python import test of master_cms config:"
python3 - <<PY
import os, traceback, importlib.util
os.environ["RINV"] = "${RINV:-0.3}"
try:
    spec = importlib.util.spec_from_file_location("master_cms_cfg", "${MASTER_CMS_CFG}")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    print("loaded ok; has config =", hasattr(m, "config"))
    if hasattr(m, "config"):
        print("config.spectrum =", getattr(m.config, "spectrum", None))
        print("config.rinv =", getattr(m.config, "rinv", None))
except Exception:
    traceback.print_exc()
    raise
PY



# --- 1) Generate events ---
case "$MODEL" in
  cms)
    echo "[generate] CMS model"
    ./run_model helper -C "$CMS_CFG" \
      --steps all --events "${NEV}" --verbose \
      | tee "${WORKDIR}/gen_${JOBTAG}_cms.log"
    ;;
  snowmass)
    echo "[generate] Snowmass CMS-like model (mrho=83.666)"
    ./run_model helper -C "$SNOW_CFG" \
      --steps all --events "${NEV}" --verbose \
      --mrho=83.666 \
      | tee "${WORKDIR}/gen_${JOBTAG}_snowmass.log"
    ;;
  master_snowmass)
    echo "[generate] MASTER config (uses configs/master_snowmass.py)"
    if [[ -z "${RINV:-}" ]]; then
      echo "[warn] RINV not provided; master_snowmass.py will default to 0.3"
    else
      echo "[cfg] Using rinv=${RINV}"
    fi

    ./run_model helper -C "$MASTER_CFG" \
      --steps all --events "${NEV}" --verbose \
      | tee "${WORKDIR}/gen_${JOBTAG}_master.log"
    ;;
  master_cms)
    echo "[generate] MASTER CMS config (uses configs/master_cms.py)"
    if [[ -z "${RINV:-}" ]]; then
      echo "[warn] RINV not provided; master_cms.py will default to 0.3"
    else
      echo "[cfg] Using rinv=${RINV}"
    fi

    ./run_model helper -C "$MASTER_CMS_CFG" \
      --steps all --events "${NEV}" --verbose \
      | tee "${WORKDIR}/gen_${JOBTAG}_master_cms.log"
    ;;
  both)
    echo "[generate] CMS model"
    ./run_model helper -C "$CMS_CFG" \
      --steps all --events "${NEV}" --verbose \
      | tee "${WORKDIR}/gen_${JOBTAG}_cms.log"

    echo "[generate] Snowmass CMS-like model (mrho=83.666)"
    ./run_model helper -C "$SNOW_CFG" \
      --steps all --events "${NEV}" --verbose \
      --mrho=83.666 \
      | tee "${WORKDIR}/gen_${JOBTAG}_snowmass.log"
    ;;
esac



GEN_END=$(date +%s)
GEN_TIME=$((GEN_END - GEN_START))    # Finish timing to run event generation


echo "[timing] Event generation finished"
echo "[timing] Generation time: ${GEN_TIME} seconds (~$((GEN_TIME/60)) min)"
echo ""



echo "[models created listed RIGHT HERE!!!!!!!!!!!!!!!]"
ls -1 models || echo "[no models directory found]"

echo "[models with config.py]"
find models -maxdepth 2 -type f -name config.py -printf '%h\n' | sort || true




echo "[debug] show produced model config rinv"
python3 - <<'PY'
import glob, os
paths = sorted(glob.glob("models/s-channel_*/config.py"))
print("found", len(paths), "config.py")
if paths:
    p = paths[-1]
    print("latest:", p)
    txt = open(p).read()
    for line in txt.splitlines():
        if "rinv" in line:
            print(line)
PY





# Make models visible to analysis
export MB_MODELS_BASE="$(pwd)/models"  
export MB_JOBTAG="${JOBTAG}"
echo "[env] MB_MODELS_BASE=${MB_MODELS_BASE}"

















# Run the analysis file
export PYTHONPATH="$(pwd):${PYTHONPATH:-}"
echo "[analysis] python add_new_DS_observables.py"
python3 "${AUTO_DIR}/scripts/add_new_DS_observables.py" 2>&1 | tee "${WORKDIR}/run_${JOBTAG}.log"




# Add root trees merging to the automation
echo "[merge] Looking for friend trees to merge ..."

while IFS= read -r friend_file; do
  model_dir="$(dirname "$friend_file")"
  model_tag="$(basename "$model_dir")"

  main_file="${model_dir}/events.root"
  merged_file="${model_dir}/merged_tree_${model_tag}_${JOBTAG}.root"

  echo "[merge] model_dir   = ${model_dir}"
  echo "[merge] main_file   = ${main_file}"
  echo "[merge] friend_file = ${friend_file}"
  echo "[merge] merged_file = ${merged_file}"

  if [[ ! -f "${main_file}" ]]; then
    echo "[merge][warn] Missing ${main_file}; skipping this model."
    continue
  fi

  if [[ ! -f "${friend_file}" ]]; then
    echo "[merge][warn] Missing ${friend_file}; skipping this model."
    continue
  fi

  root -l -b -q "${AUTO_DIR}/scripts/merge_trees.C(\"${main_file}\",\"${friend_file}\",\"${merged_file}\")"

  MERGE_STATUS=$?
  if [[ ${MERGE_STATUS} -ne 0 ]]; then
    echo "[merge][fatal] merge_trees.C failed with status ${MERGE_STATUS}"
    exit ${MERGE_STATUS}
  fi

  echo "[merge] Successfully wrote ${merged_file}"

done < <(find "models" -type f -name 'events_friend*.root')




# Collect outputs and stage to EOS
cd "${WORKDIR}"
echo "[stage-out] collecting outputs ..."



# Condor job identifiers (prefer env; fallback to JOBTAG like j82103242.0)
CLUSTER="${ClusterId:-}"
PROC="${ProcId:-}"
if [[ -z "$CLUSTER" || -z "$PROC" ]]; then
  if [[ "$JOBTAG" =~ ^j([0-9]+)\.([0-9]+)$ ]]; then
    CLUSTER="${BASH_REMATCH[1]}"
    PROC="${BASH_REMATCH[2]}"
  else
    CLUSTER="unknown"
    PROC="0"
  fi
fi


# Include the job tag and the value of the variable in the name for EOS









# OUTROOT_BASE="/store/user/${EOS_USER}/${PROJECT}/outputs"
OUTROOT_BASE="/store/user/${EOS_USER}/${EOS_AREA}/outputs"















if [[ "${MODEL}" == "master_snowmass" || "${MODEL}" == "master_cms" ]]; then
  JOB_EOS_DIR="${OUTROOT_BASE}/j${CLUSTER}.${PROC}-${MODEL}-rinv=${RINV}"
else
  JOB_EOS_DIR="${OUTROOT_BASE}/j${CLUSTER}.${PROC}"
fi
JOB_EOS_RSE="root://cmseos.fnal.gov/${JOB_EOS_DIR}"


# Ensure job-level directory exists on EOS
xrdfs root://cmseos.fnal.gov mkdir -p "${JOB_EOS_DIR}" || true



# Stage friend ROOT files into per-model subfolders
echo "[stage-out] friend trees by model ..."
while IFS= read -r f; do
  # f is something like: ${PROJECT}/models/<model_tag>/events_friend_*.root
  model_dir="$(dirname "$f")"                                             
  model_dirname="$(basename "$model_dir")"                                
  model_dirname="${model_dirname//[^A-Za-z0-9._-]/-}"
  model_eos_dir="${JOB_EOS_DIR}/${model_dirname}"
  xrdfs root://cmseos.fnal.gov mkdir -p "${model_eos_dir}" || true

  echo "  -> $(basename "$f")  -->  ${model_eos_dir}/"
  xrdcp -f "$f" "root://cmseos.fnal.gov/${model_eos_dir}/"


  if [[ -f "${model_dir}/events.root" ]]; then
    echo "  -> events.root  -->  ${model_eos_dir}/"
    xrdcp -f "${model_dir}/events.root" "root://cmseos.fnal.gov/${model_eos_dir}/"
  else
    alt_events="$(find "${model_dir}" -maxdepth 1 -type f \( -name 'events_*.root' -o -name 'Events.root' \) | head -n1)"
    if [[ -n "$alt_events" ]]; then
      echo "  -> $(basename "$alt_events")  -->  ${model_eos_dir}/"
      xrdcp -f "$alt_events" "root://cmseos.fnal.gov/${model_eos_dir}/"
    fi
  fi
  
  
  
    # Stage merged ROOT files
  shopt -s nullglob
  for merged in "${model_dir}"/merged_tree*.root; do
    echo "  -> $(basename "$merged")  -->  ${model_eos_dir}/"
    xrdcp -f "$merged" "root://cmseos.fnal.gov/${model_eos_dir}/"
  done
  shopt -u nullglob
  
  
  

  # Stage model card/config files
  for card in config.py pythia_card.txt delphes_card.txt; do
    if [[ -f "${model_dir}/${card}" ]]; then
      echo "  -> ${card}  -->  ${model_eos_dir}/"
      xrdcp -f "${model_dir}/${card}" "root://cmseos.fnal.gov/${model_eos_dir}/"
    else
      echo "  [warn] ${card} not found in ${model_dir}"
    fi
  done
done < <(find "${PROJECT}/models" -type f -name 'events_friend*.root')



# # Stage plots (pdf/png)
# echo "[stage-out] plots ..."
# shopt -s nullglob
# for p in $(find "${PROJECT}" -type f \( -name '*.pdf' -o -name '*.png' \)); do
#   xrdcp -f "$p" "${JOB_EOS_RSE}/"
# done
# shopt -u nullglob


# Stage logs
echo "[stage-out] logs ..."
xrdcp -f "run_${JOBTAG}.log" "${JOB_EOS_RSE}/" 2>/dev/null || true
shopt -s nullglob
for g in gen_${JOBTAG}_*.log; do
  xrdcp -f "$g" "${JOB_EOS_RSE}/"
done
shopt -u nullglob

echo "Done."
date
