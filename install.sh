#!/bin/bash

# empty argument for init.sh
source init.sh ""

ECHO="echo -e"
export PYTHIA8MINOR=10
export ACCESS=ssh
export REINSTALL_EXTERNALS=

usage() {
	$ECHO "install.sh [options]"
	$ECHO
	$ECHO "Options:"
	$ECHO "-e [external(s)] \t(re)install specific external(s) (options = ${MODEL_BUILDING_EXTERNALS})"
	$ECHO "-p [version]     \tPythia8.3 minor version to install (default = ${PYTHIA8MINOR})"
	$ECHO "-s [protocol]    \tuse protocol for git clone (default = ${ACCESS}, alternative = https)"
	$ECHO "-h               \tprint this message and exit"
	exit $1
}

while getopts "e:p:s:h" opt; do
	case "$opt" in
		e) export REINSTALL_EXTERNALS="$OPTARG"
		;;
		p) export PYTHIA8MINOR="$OPTARG"
		;;
		s) export ACCESS="$OPTARG"
		;;
		h) usage 0
		;;
	esac
done

if [ -z "$REINSTALL_EXTERNALS" ]; then
	REINSTALL_EXTERNALS="$MODEL_BUILDING_EXTERNALS"
else
	FOUND_EXTERNALS=true
	for EXT in $REINSTALL_EXTERNALS; do
		if [[ ! " $MODEL_BUILDING_EXTERNALS " =~ " $EXT " ]]; then
			echo "External $EXT is not in allowed list"
			FOUND_EXTERNALS=false
		fi
	done
	if [ "$FOUND_EXTERNALS" = false ]; then
		echo "Allowed list: $MODEL_BUILDING_EXTERNALS"
		exit 1
	fi
fi

if [ "$ACCESS" = "ssh" ]; then
	export ACCESS_GITHUB=git@github.com:
elif [ "$ACCESS" = "https" ]; then
	export ACCESS_GITHUB=https://github.com/
else
	usage 1
fi

cd install

for EXTERNAL in $REINSTALL_EXTERNALS; do
	echo $EXTERNAL
	./${EXTERNAL}.sh
done
