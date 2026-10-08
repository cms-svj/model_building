#!/usr/bin/env bash
set -euo pipefail

# Find the repo from this script's location
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$(dirname "${REPO_DIR}")"
REPO="$(basename "${REPO_DIR}")"

# EOS user defaults to your own username
EOS_USER="${EOS_USER:-$USER}"

# EOS area where the bundle and outputs go
EOS_AREA="${EOS_AREA:-model_building_fork}"

TARBALL="${BASE}/${REPO}_bundle.tgz"
EOS_DIR="/store/user/${EOS_USER}/${EOS_AREA}/inputs"
EOS_URL="root://cmseos.fnal.gov/${EOS_DIR}"

echo "[pack] making ${TARBALL} from ${REPO_DIR}"
cd "${BASE}"

tar -czf "${TARBALL}" \
  --exclude="${REPO}/.git" \
  --exclude="${REPO}/logs" \
  --exclude="${REPO}/jobs" \
  --exclude="${REPO}/models" \
  --exclude="${REPO}/automation_of_dark_sector_variables/condor/logs" \
  --exclude="${REPO}/__pycache__" \
  --exclude="${REPO}/.ipynb_checkpoints" \
  --exclude="*.root" \
  --exclude="*.pdf" \
  --exclude="*.png" \
  "${REPO}"

echo "[pack] uploading to ${EOS_URL}"
xrdfs root://cmseos.fnal.gov/ mkdir -p "${EOS_DIR}" || true
xrdcp -f "${TARBALL}" "${EOS_URL}/${REPO}_bundle.tgz"

echo "[pack] done."