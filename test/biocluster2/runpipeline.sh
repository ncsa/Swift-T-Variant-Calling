#!/bin/bash

set -x

scriptdir="/home/a-m/azzaea/azza/swift_T_project/src/Swift-T-Variant-Calling"

swift-t -m slurm -l -V -s settings.sh -u -I $scriptdir/src -r $scriptdir/src/bioapps $scriptdir/src/VariantCalling.swift -runfile=$PWD/H3A_NextGen_assessment.Chr1_50X.set3.runfile | tee $PWD/log

