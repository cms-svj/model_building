#!/bin/bash

PYTHIA_VERSION=pythia8316

git clone git@github.com:kpedro88/pythia8 -b cms/316
cd pythia8

HEPMC_BASE=(/cvmfs/sft.cern.ch/lcg/releases/${LCG_VIEW}/HepMC/*/${LCG_ARCH})
LHAPDF_BASE=(/cvmfs/sft.cern.ch/lcg/releases/${LCG_VIEW}/MCGenerators/lhapdf/*/${LCG_ARCH})

./configure --with-hepmc2=${HEPMC_BASE} --with-lhapdf6=${LHAPDF_BASE} --with-python --with-gzip
make -j 8
make install

(cd examples
make main132
)

cat << 'EOF' > mb_init.sh
export PYTHIA8=${MODEL_BUILDING}/install/pythia8
export PATH=${PYTHIA8}/bin/:${PATH}
export LD_LIBRARY_PATH=${PYTHIA8}/lib/:${LD_LIBRARY_PATH}
export PYTHIA8DATA=$(pythia8-config --xmldoc)
export PYTHIA8RUNNER=${PYTHIA8}/examples/main132
EOF
