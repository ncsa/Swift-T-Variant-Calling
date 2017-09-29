#!/bin/bash

module load swift-t

set -x

scriptdir="/home/a-m/azzaea/swift_T_project/src/Swift-T-Variant-Calling"

swift-t -m slurm -l -V -s settings.sh -u -I $scriptdir/src -r $scriptdir/src/bioapps $scriptdir/src/VariantCalling.swift -runfile=$PWD/synthetic_reads.runfile | tee $PWD/log.swift_t_run

