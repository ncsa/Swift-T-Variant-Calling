#!/bin/bash

module load swift-t

set -x
export PROCS=6 #total number of MPI processes
export TURBINE_LOG=1
export ADLB_DEBUG_RANKS=1
export SWIFT_TMP=/home/aeahmed/swift_T_project/results/WGS_chr1_30X_E0.005/tmp_swift

dstat -t -c --disk-util --top-cpu --top-io --top-mem -s --output $PWD/profile.complete_pipeline.log 1 18000 &
dsprocess=$(ps -aux | grep dstat | awk '{ print $2}' |head -n 1)

scriptdir="/home/aeahmed/swift_T_project/src/Swift-T-Variant-Calling"

swift-t -O3 -l -u -I $scriptdir/src -r $scriptdir/src/bioapps $scriptdir/src/VariantCalling.swift -runfile=$PWD/WGS_chr1_30X_E0.005.runfile | tee $PWD/log.swift_t_run

kill -9 $dsprocess
echo -e "Swift-T pipeline run on $HOSTNAME has concluded successfully!" | mail -s "swift_t_pipeline" "azzaea@gmail.com"
