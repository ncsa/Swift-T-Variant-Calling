#!/bin/bash

#var=$(./autoincrement.sh)
#sed -i '/^OUTPUT/s/[0-9]*$/'"$var"'/' HgG0.lowcoverage.chr20.runfile

export TURBINE_LOG=1
#export ADLB_PERF_COUNTERS=1

### the variables below didn't really chip in much!
#export  ADLB_REPORT_LEAKS=1
#export  ADLB_TRACE=true
#export ADLB_DEBUG=true 
export  ADLB_DEBUG_RANKS=1
export TURBINE_OUTPUT=$PWD
export export SWIFT_TMP=/work/04525/tg838247/stampede2/swift_T_project/swift_tmp

scriptdir="/home1/04525/tg838247/swift_T_project/src/Swift-T-Variant-Calling"

swift-t -I $scriptdir/src -r $scriptdir/src/bioapps -u -n 3  $scriptdir/src/VariantCalling.swift -runfile=$PWD/GIAB_Garvan_NA12878_HG001_HiSeq_Exome.runfile | tee $PWD/log
