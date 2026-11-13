#!/bin/bash

# optional argument: Pythia8.3 minor version to use (default = install/pythia8)
PYTHIA8SELECT=$1

export LCG_VIEW=LCG_106
export LCG_ARCH=x86_64-el9-gcc13-opt
source /cvmfs/sft.cern.ch/lcg/views/${LCG_VIEW}/${LCG_ARCH}/setup.sh

export MODEL_BUILDING=$PWD
export MODEL_BUILDING_EXTERNALS="pythia8 python_packages delphes"

for EXTERNAL in $MODEL_BUILDING_EXTERNALS; do
	EXTERNAL_DIR=${MODEL_BUILDING}/install/${EXTERNAL}
	if [ "$EXTERNAL" = pythia8 ] && [ -n "$PYTHIA8SELECT" ]; then
		EXTERNAL_DIR=${MODEL_BUILDING}/install/pythia83${PYTHIA8SELECT}
		if [ ! -d ${EXTERNAL_DIR} ]; then
			echo "Pythia 8.3${PYTHIA8SELECT} is not installed (use ./install.sh -e pythia8 -p ${PYTHIA8SELECT})"
			return 1
		fi
	fi
	if [ -d ${EXTERNAL_DIR} ]; then
		source ${EXTERNAL_DIR}/mb_init.sh
	fi
done
