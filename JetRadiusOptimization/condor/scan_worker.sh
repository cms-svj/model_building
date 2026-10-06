#!/bin/bash
# Worker for one (model point, event chunk) of the jet-radius scan.
#
# Runs the skim stage once, then every radius on top of it, and ships a single
# Parquet shard back to EOS.  Unlike the per-radius layout it replaces, the
# input ROOT file is copied once and the dark-hadron ancestry is resolved once,
# no matter how many radii are in the grid.
set -euo pipefail

POINT_ID="$1"
CHUNK="$2"
POINTS_FILE="$3"
RADII="$4"
EOS_OUT="$5"

EOS_HOST="${EOS_HOST:-root://cmseos.fnal.gov}"
EVENTS_PER_CHUNK="${EVENTS_PER_CHUNK:-100000}"

echo "============================================================"
echo "point_id : ${POINT_ID}"
echo "chunk    : ${CHUNK}"
echo "radii    : ${RADII}"
echo "host     : $(hostname)"
echo "started  : $(date -Is)"
echo "============================================================"

# --- resolve this point's parameters from the shared table -----------------
read -r _ MMED NC NF MPI MRHO PVECTOR RINV LABEL EVENTS_ROOT < <(
    awk -v id="${POINT_ID}" '!/^#/ && $1 == id' "${POINTS_FILE}"
)
if [[ -z "${LABEL:-}" ]]; then
    echo "[FATAL] point_id ${POINT_ID} not found in ${POINTS_FILE}" >&2
    exit 2
fi
echo "[INFO] ${LABEL}  (mmed=${MMED} mpi=${MPI} mrho=${MRHO} rinv=${RINV})"

tar xzf payload.tgz
export PYTHONPATH="${PWD}/payload:${PYTHONPATH:-}"

# --- stage the input -------------------------------------------------------
ENTRY_START=$(( CHUNK * EVENTS_PER_CHUNK ))
echo "[INFO] copying ${EVENTS_ROOT}"
for attempt in 1 2 3 4 5; do
    if xrdcp -f "${EVENTS_ROOT}" ./events.root; then break; fi
    echo "[WARN] xrdcp attempt ${attempt} failed" >&2
    sleep $(( attempt * 15 ))
    [[ ${attempt} -eq 5 ]] && { echo "[FATAL] input copy failed" >&2; exit 3; }
done

# --- stage 1: radius-independent skim (ancestry resolved exactly once) -----
echo "[INFO] skim"
python3 payload/JetRadiusOptimization/skim.py \
    --input events.root \
    --output "skim_${POINT_ID}_c${CHUNK}.parquet" \
    --entry-start "${ENTRY_START}" \
    --max-events "${EVENTS_PER_CHUNK}" \
    --label "${LABEL}" \
    --point-id "${POINT_ID}"

# --- stage 2: every radius on top of the skim ------------------------------
echo "[INFO] scan"
python3 payload/JetRadiusOptimization/scan.py \
    --skim "skim_${POINT_ID}_c${CHUNK}.parquet" \
    --output "metrics_${POINT_ID}_c${CHUNK}.parquet" \
    --radii "${RADII}"

# --- ship results ----------------------------------------------------------
DEST="${EOS_HOST}/${EOS_OUT}/${LABEL}"
xrdfs "${EOS_HOST}" mkdir -p "${EOS_OUT}/${LABEL}" || true
for f in "metrics_${POINT_ID}_c${CHUNK}.parquet"; do
    for attempt in 1 2 3; do
        if xrdcp -f "${f}" "${DEST}/${f}"; then break; fi
        echo "[WARN] upload attempt ${attempt} failed for ${f}" >&2
        sleep $(( attempt * 15 ))
        [[ ${attempt} -eq 3 ]] && { echo "[FATAL] upload failed" >&2; exit 4; }
    done
done

# Keep the skim only if asked -- it is ~10x the metrics table, but it is also
# what lets you re-run a changed metric definition without touching ROOT again.
if [[ "${KEEP_SKIM:-0}" == "1" ]]; then
    xrdcp -f "skim_${POINT_ID}_c${CHUNK}.parquet" \
        "${DEST}/skim_${POINT_ID}_c${CHUNK}.parquet"
fi

echo "[DONE] ${LABEL} chunk ${CHUNK} at $(date -Is)"
