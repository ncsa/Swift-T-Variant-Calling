#!/bin/bash

set -x

scriptdir="/home1/04525/tg838247/swift_T_project/src/Swift-T-Variant-Calling"

swift-t -m slurm -l -V -s settings.sh -u -I $scriptdir/src -r $scriptdir/src/bioapps $scriptdir/src/VariantCalling.swift -runfile=$PWD/GIAB_Garvan_NA12878_HG001_HiSeq_Exome.runfile | tee $PWD/log
