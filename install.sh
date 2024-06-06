#!/bin/bash

source init.sh

cd install

for EXTERNAL in $MODEL_BUILDING_EXTERNALS; do
	echo $EXTERNAL
	./${EXTERNAL}.sh
done
