#!/bin/bash

source init.sh

ECHO="echo -e"
export PYTHIA8MINOR=10
export ACCESS=ssh

usage() {
	$ECHO "install.sh [options]"
	$ECHO
	$ECHO "Options:"
	$ECHO "-p [version]    \tPythia8.3 minor version to install (default = ${PYTHIA8MINOR})"
	$ECHO "-c [protocol]   \tuse protocol for git clone (default = ${ACCESS}, alternative = https)"
	$ECHO "-h              \tprint this message and exit"
	exit $1
}

while getopts "p:s:h" opt; do
	case "$opt" in
		p) export PYTHIA8MINOR="$OPTARG"
		;;
		s) export ACCESS="$OPTARG"
		;;
		h) usage 0
		;;
	esac
done

if [ "$ACCESS" = "ssh" ]; then
	export ACCESS_GITHUB=git@github.com:
elif [ "$ACCESS" = "https" ]; then
	export ACCESS_GITHUB=https://github.com/
else
	usage 1
fi

cd install

for EXTERNAL in $MODEL_BUILDING_EXTERNALS; do
	echo $EXTERNAL
	./${EXTERNAL}.sh
done
