#!/bin/bash

git clone ${ACCESS_GITHUB}cms-svj/delphes

if [ ! -d delphes ]; then
	echo "Failed to clone delphes"
	exit 1
fi

cd delphes
git fetch origin
git checkout DarkHadronJets_3.5.2pre02

make -j 8

cat << 'EOF' > mb_init.sh
export DELPHES=${MODEL_BUILDING}/install/delphes
export PATH=${DELPHES}/:${PATH}
export LD_LIBRARY_PATH=${DELPHES}/:${LD_LIBRARY_PATH}
EOF
