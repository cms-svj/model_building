#!/bin/bash

PYTHIA_VERSION=pythia83${PYTHIA8MINOR}
PYTHIA_EXE=main42
if [ "$PYTHIA8MINOR" -ge 12 ]; then
	PYTHIA_EXE=main132
fi

# remove an old unversioned install (it can't be moved: its library rpath points to install/pythia8/lib)
if [ -d pythia8 ] && [ ! -L pythia8 ]; then
	rm -rf pythia8
fi

if [ ! -d ${PYTHIA_VERSION} ]; then
	wget "https://pythia.org/releases/pythia83/${PYTHIA_VERSION}.tgz"
	if [ ! -e ${PYTHIA_VERSION}.tgz ]; then
		echo "Failed to download pythia"
		exit 1
	fi
	tar -xzf ${PYTHIA_VERSION}.tgz
	rm ${PYTHIA_VERSION}.tgz
fi
# make this version the default
ln -sfn ${PYTHIA_VERSION} pythia8
cd ${PYTHIA_VERSION}

HEPMC_BASE=(/cvmfs/sft.cern.ch/lcg/releases/${LCG_VIEW}/hepmc3/*/${LCG_ARCH})
LHAPDF_BASE=(/cvmfs/sft.cern.ch/lcg/releases/${LCG_VIEW}/MCGenerators/lhapdf/*/${LCG_ARCH})

./configure --with-hepmc3=${HEPMC_BASE} --with-hepmc3-lib=${HEPMC_BASE}/lib64 --with-lhapdf6=${LHAPDF_BASE} --with-python --with-gzip
make -j 8
make install

(cd examples
make ${PYTHIA_EXE}
)

cat << EOF > mb_init.sh
export PYTHIA8=\${MODEL_BUILDING}/install/${PYTHIA_VERSION}
export PATH=\${PYTHIA8}/bin/:\${PATH}
export LD_LIBRARY_PATH=\${PYTHIA8}/lib/:\${LD_LIBRARY_PATH}
export PYTHIA8DATA=\$(pythia8-config --xmldoc)
export PYTHIA8MINOR=${PYTHIA8MINOR}
export PYTHIA8RUNNER=\${PYTHIA8}/examples/${PYTHIA_EXE}
EOF
