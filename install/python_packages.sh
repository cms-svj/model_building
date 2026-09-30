#!/bin/bash

mkdir -p python_packages
cd python_packages

wget https://github.com/kpedro88/lcg-venv/raw/refs/heads/main/lcg-venv
chmod +x lcg-venv
./lcg-venv mbenv

cat << 'EOF' > mb_init.sh
source ${MODEL_BUILDING}/install/python_packages/mbenv/bin/activate
export LD_LIBRARY_PATH=${MODEL_BUILDING}/install/python_packages/mbenv/lib/python3.11/site-packages/fastjet/lib64:${LD_LIBRARY_PATH}
EOF

source mb_init.sh

# vector 1.9.0 validates that a record has only one 4th coordinate, but the
# coffea Delphes Particle mixin carries both E and Mass, so every particle
# collection fails to load. coffea only pins vector>=1.4.1, so pin it here.
PKGS_UPGRADE=(
vector==1.8.1 \
coffea==2025.12.0 \
mplhep \
)

for PKG in ${PKGS_UPGRADE[@]}; do
	pip install --upgrade $PKG
done

PKGS=(
magiconfig \
fastjet \
)

for PKG in ${PKGS[@]}; do
	pip install $PKG
done
