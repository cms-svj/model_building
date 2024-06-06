#!/bin/bash

mkdir -p python_packages
cd python_packages

PKGS=(
magiconfig \
)

for PKG in ${PKGS[@]}; do
	HOME=$PWD pip install --user --no-cache-dir $PKG
done

cat << 'EOF' > mb_init.sh
export PYTHONPATH=${MODEL_BUILDING}/install/python_packages/.local/lib/python3.9/site-packages/:${PYTHONPATH}
EOF
