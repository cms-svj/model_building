#!/bin/bash

set -euo pipefail

INPUT_URL="${1:?missing EOS input URL}"
EOS_OUTPUT_DIR="${2:?missing EOS output directory}"
MAX_EVENTS="${3:-10000}"
JOB_SCRATCH="${_CONDOR_SCRATCH_DIR:-$(pwd)}"
EOS_HOST="root://cmseos.fnal.gov"
WORK_REPOSITORY="${JOB_SCRATCH}/model_building_transfer"
LOCAL_INPUT="${JOB_SCRATCH}/events.root"
LOCAL_OUTPUT="${JOB_SCRATCH}/jet_radius_validation"

echo "============================================================"
echo "Jet-radius analysis worker"
echo "host           : $(hostname)"
echo "start time     : $(date --iso-8601=seconds)"
echo "scratch        : ${JOB_SCRATCH}"
echo "input          : ${INPUT_URL}"
echo "events         : ${MAX_EVENTS}"
echo "EOS output     : ${EOS_HOST}/${EOS_OUTPUT_DIR}"
echo "============================================================"

for required in \
    "${WORK_REPOSITORY}/init.sh" \
    "${WORK_REPOSITORY}/common.py" \
    "${JOB_SCRATCH}/validate.py" \
    "${JOB_SCRATCH}/core.py" \
    "${JOB_SCRATCH}/ancestry.py"; do
    if [[ ! -e "${required}" ]]; then
        echo "[ERROR] Transferred input is missing: ${required}" >&2
        exit 10
    fi
done

mkdir -p "${WORK_REPOSITORY}/JetRadiusOptimization" "${LOCAL_OUTPUT}"
cp "${JOB_SCRATCH}/validate.py" "${WORK_REPOSITORY}/JetRadiusOptimization/validate.py"
cp "${JOB_SCRATCH}/core.py" "${WORK_REPOSITORY}/JetRadiusOptimization/core.py"
cp "${JOB_SCRATCH}/ancestry.py" "${WORK_REPOSITORY}/JetRadiusOptimization/ancestry.py"

cd "${WORK_REPOSITORY}"
set +u
source init.sh
set -u

# Use the pinned LCG interpreter. The transferred repository overlay supplies
# only packages absent from LCG, notably Python fastjet and magiconfig.
TRANSFERRED_VENV="${WORK_REPOSITORY}/install/python_packages/mbenv"
unset VIRTUAL_ENV
export PATH="/cvmfs/sft.cern.ch/lcg/views/${LCG_VIEW}/${LCG_ARCH}/bin:${PATH}"
export PYTHONPATH="${TRANSFERRED_VENV}/lib/python3.11/site-packages:${PYTHONPATH:-}"
export PYTHONNOUSERSITE=1

for command in python3 xrdcp xrdfs; do
    if ! command -v "${command}" >/dev/null 2>&1; then
        echo "[ERROR] Runtime command is unavailable: ${command}" >&2
        exit 11
    fi
done

python3 - <<'PY'
import fastjet
import magiconfig
import sys

print(f"[INFO] Python: {sys.executable}")
print(f"[INFO] magiconfig: {magiconfig.__file__}")
print(f"[INFO] fastjet: {fastjet.__file__}")
PY

if xrdfs "${EOS_HOST}" stat "${EOS_OUTPUT_DIR}/_SUCCESS.json" >/dev/null 2>&1; then
    echo "[SKIP] Analysis success marker already exists."
    exit 0
fi

echo "[INFO] Copying the Delphes ROOT input to worker scratch"
xrdcp -f --nopbar "${INPUT_URL}" "${LOCAL_INPUT}"

echo "[INFO] Running the 10,000-event radius analysis"
set +e
/usr/bin/time -v python3 JetRadiusOptimization/validate.py \
    --input "${LOCAL_INPUT}" \
    --outdir "${LOCAL_OUTPUT}" \
    --max-events "${MAX_EVENTS}" \
    --collections GenFatJet \
    --radii \
        0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 \
        1.0 1.1 1.2 1.3 1.4 1.5 1.6 \
    --diagnostic-radii \
        0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 \
    --event-displays 12 \
    --bootstrap-resamples 500 \
    2>&1 | tee "${LOCAL_OUTPUT}/analysis.log"
analysis_status=${PIPESTATUS[0]}
set -e
if [[ ! -f "${LOCAL_OUTPUT}/validation_report.json" ]]; then
    echo "[ERROR] Radius analysis produced no validation_report.json" >&2
    if [[ "${analysis_status}" -eq 0 ]]; then
        analysis_status=1
    fi
    exit "${analysis_status}"
fi

VALIDATION_STATUS="$(python3 - "${LOCAL_OUTPUT}/validation_report.json" <<'PY'
import json
import sys

print(json.load(open(sys.argv[1])).get("status", "unknown"))
PY
)"
if [[ "${analysis_status}" -ne 0 ]]; then
    echo "[WARN] Validation status is ${VALIDATION_STATUS}; uploading diagnostic artifacts anyway."
fi

echo "[INFO] Copying plots, JSON report, and analysis log to EOS"
xrdfs "${EOS_HOST}" mkdir -p "${EOS_OUTPUT_DIR}"
while IFS= read -r -d '' local_file; do
    relative="${local_file#${LOCAL_OUTPUT}/}"
    remote_parent="${EOS_OUTPUT_DIR}/$(dirname "${relative}")"
    xrdfs "${EOS_HOST}" mkdir -p "${remote_parent}"
    xrdcp -f --nopbar "${local_file}" "${EOS_HOST}/${EOS_OUTPUT_DIR}/${relative}"
done < <(find "${LOCAL_OUTPUT}" -type f ! -name '_SUCCESS.json' -print0)

if [[ "${VALIDATION_STATUS}" == "passed" ]]; then
    MARKER_NAME="_SUCCESS.json"
else
    MARKER_NAME="_ARTIFACTS_READY.json"
fi
export INPUT_URL EOS_OUTPUT_DIR MAX_EVENTS VALIDATION_STATUS
python3 - "${LOCAL_OUTPUT}/${MARKER_NAME}" <<'PY'
from datetime import datetime, timezone
import json
import os
import platform
import sys

output = {
    "status": "complete" if os.environ["VALIDATION_STATUS"] == "passed" else "artifacts_ready",
    "validation_status": os.environ["VALIDATION_STATUS"],
    "completed_utc": datetime.now(timezone.utc).isoformat(),
    "host": platform.node(),
    "input": os.environ["INPUT_URL"],
    "eos_output_dir": os.environ["EOS_OUTPUT_DIR"],
    "events_requested": int(os.environ["MAX_EVENTS"]),
    "lcg_view": os.environ.get("LCG_VIEW"),
    "lcg_arch": os.environ.get("LCG_ARCH"),
}
with open(sys.argv[1], "w") as handle:
    json.dump(output, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
xrdcp -f --nopbar \
    "${LOCAL_OUTPUT}/${MARKER_NAME}" \
    "${EOS_HOST}/${EOS_OUTPUT_DIR}/${MARKER_NAME}"

echo "[DONE] Radius analysis artifacts copied to EOS (${MARKER_NAME})"
