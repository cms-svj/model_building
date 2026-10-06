#!/bin/bash
# Build the runtime snapshot that analysis.sub and focused_by_radius.sub ship to
# each worker as model_building_transfer/: init.sh, common.py and the installed
# python_packages venv (fastjet, magiconfig, coffea overlay on LCG).
#
# Run once after ./install.sh, and again whenever init.sh, common.py or the
# venv changes:
#     JetRadiusOptimization/condor/make_transfer_dir.sh

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DEST="${REPO}/JetRadiusOptimization/condor/model_building_transfer"

if [[ ! -d "${REPO}/install/python_packages/mbenv" ]]; then
    echo "[ERROR] ${REPO}/install/python_packages/mbenv not found; run ./install.sh first" >&2
    exit 1
fi

rm -rf "${DEST}"
mkdir -p "${DEST}/install"
cp "${REPO}/init.sh" "${REPO}/common.py" "${DEST}/"
cp -a "${REPO}/install/python_packages" "${DEST}/install/"
echo "[DONE] ${DEST} ($(du -sh "${DEST}" | cut -f1))"
