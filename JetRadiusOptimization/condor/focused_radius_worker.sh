#!/bin/bash

set -euo pipefail

INPUT_URL="${1:?missing EOS input URL}"
EOS_SHARD_BASE="${2:?missing EOS shard base}"
RADIUS="${3:?missing radius}"
RADIUS_INDEX="${4:?missing radius index}"
MAX_EVENTS="${5:-10000}"
JOB_SCRATCH="${_CONDOR_SCRATCH_DIR:-$(pwd)}"
EOS_HOST="root://cmseos.fnal.gov"
WORK_REPOSITORY="${JOB_SCRATCH}/model_building_transfer"
LOCAL_INPUT="${JOB_SCRATCH}/events.root"
LOCAL_OUTPUT="${JOB_SCRATCH}/focused_output"
RADIUS_TAG="R$(printf '%.1f' "${RADIUS}" | tr '.' 'p')"
EOS_OUTPUT_DIR="${EOS_SHARD_BASE}/${RADIUS_TAG}"

echo "============================================================"
echo "Focused jet-radius worker"
echo "radius         : ${RADIUS}"
echo "radius index   : ${RADIUS_INDEX}"
echo "events         : ${MAX_EVENTS}"
echo "host           : $(hostname)"
echo "EOS output     : ${EOS_HOST}/${EOS_OUTPUT_DIR}"
echo "============================================================"

for required in \
    "${WORK_REPOSITORY}/init.sh" \
    "${WORK_REPOSITORY}/common.py" \
    "${JOB_SCRATCH}/focused_diagnostics.py" \
    "${JOB_SCRATCH}/core.py" \
    "${JOB_SCRATCH}/ancestry.py"; do
    if [[ ! -e "${required}" ]]; then
        echo "[ERROR] Missing transferred input: ${required}" >&2
        exit 10
    fi
done

mkdir -p "${WORK_REPOSITORY}/JetRadiusOptimization" "${LOCAL_OUTPUT}"
cp "${JOB_SCRATCH}/focused_diagnostics.py" \
    "${WORK_REPOSITORY}/JetRadiusOptimization/focused_diagnostics.py"
cp "${JOB_SCRATCH}/core.py" "${WORK_REPOSITORY}/JetRadiusOptimization/core.py"
cp "${JOB_SCRATCH}/ancestry.py" "${WORK_REPOSITORY}/JetRadiusOptimization/ancestry.py"

cd "${WORK_REPOSITORY}"
set +u
source init.sh
set -u
TRANSFERRED_VENV="${WORK_REPOSITORY}/install/python_packages/mbenv"
unset VIRTUAL_ENV
export PATH="/cvmfs/sft.cern.ch/lcg/views/${LCG_VIEW}/${LCG_ARCH}/bin:${PATH}"
export PYTHONPATH="${TRANSFERRED_VENV}/lib/python3.11/site-packages:${PYTHONPATH:-}"
export PYTHONNOUSERSITE=1

echo "[INFO] Copying ROOT input to local worker scratch"
xrdcp -f --nopbar "${INPUT_URL}" "${LOCAL_INPUT}"

echo "[INFO] Processing R=${RADIUS}"
/usr/bin/time -v python3 JetRadiusOptimization/focused_diagnostics.py \
    --input "${LOCAL_INPUT}" \
    --outdir "${LOCAL_OUTPUT}" \
    --entry-start 0 \
    --max-events "${MAX_EVENTS}" \
    --shard-id "${RADIUS_INDEX}" \
    --radii "${RADIUS}"

echo "[INFO] Uploading mergeable metrics"
xrdfs "${EOS_HOST}" mkdir -p "${EOS_OUTPUT_DIR}"
for local_file in "${LOCAL_OUTPUT}"/*; do
    xrdcp -f --nopbar "${local_file}" \
        "${EOS_HOST}/${EOS_OUTPUT_DIR}/$(basename "${local_file}")"
done

echo "[DONE] R=${RADIUS}"
