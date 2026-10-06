#!/bin/bash

set -euo pipefail

EOS_SHARD_BASE="${1:?missing EOS shard base}"
EOS_FINAL_DIR="${2:?missing final EOS directory}"
JOB_SCRATCH="${_CONDOR_SCRATCH_DIR:-$(pwd)}"
EOS_HOST="root://cmseos.fnal.gov"
LOCAL_INPUTS="${JOB_SCRATCH}/radius_inputs"
LOCAL_OUTPUT="${JOB_SCRATCH}/focused_merged"

set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
set -u
export PYTHONNOUSERSITE=1
mkdir -p "${LOCAL_INPUTS}" "${LOCAL_OUTPUT}"

radii=(0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6)
for index in "${!radii[@]}"; do
    radius="${radii[$index]}"
    radius_tag="R$(printf '%.1f' "${radius}" | tr '.' 'p')"
    shard_id="$(printf '%02d' "${index}")"
    remote="${EOS_HOST}/${EOS_SHARD_BASE}/${radius_tag}/focused_metrics_shard_${shard_id}.npz"
    local_file="${LOCAL_INPUTS}/focused_metrics_R${radius}.npz"
    echo "[INFO] Downloading ${remote}"
    xrdcp -f --nopbar "${remote}" "${local_file}"
done

python3 "${JOB_SCRATCH}/merge_focused_diagnostics.py" \
    --inputs "${LOCAL_INPUTS}"/*.npz \
    --outdir "${LOCAL_OUTPUT}" \
    --diagnostic-radii 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6

xrdfs "${EOS_HOST}" mkdir -p "${EOS_FINAL_DIR}"
for local_file in "${LOCAL_OUTPUT}"/*; do
    xrdcp -f --nopbar "${local_file}" \
        "${EOS_HOST}/${EOS_FINAL_DIR}/$(basename "${local_file}")"
done

python3 - "${LOCAL_OUTPUT}/_FOCUSED_DIAGNOSTICS_SUCCESS.json" <<'PY'
from datetime import datetime, timezone
import json
import sys

with open(sys.argv[1], "w") as handle:
    json.dump(
        {
            "status": "complete",
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "execution": "15 independent radius jobs followed by one merge job",
        },
        handle,
        indent=2,
        sort_keys=True,
    )
    handle.write("\n")
PY
xrdcp -f --nopbar "${LOCAL_OUTPUT}/_FOCUSED_DIAGNOSTICS_SUCCESS.json" \
    "${EOS_HOST}/${EOS_FINAL_DIR}/_FOCUSED_DIAGNOSTICS_SUCCESS.json"

echo "[DONE] Merged focused diagnostics uploaded"
