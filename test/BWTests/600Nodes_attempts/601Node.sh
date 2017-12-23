#!/usr/bin/env bash
set -eu

export PPN=1
export NODES=601
export PROCS=$(($PPN * $NODES))
export WALLTIME=12:00:00
export PROJECT=baib
export QUEUE=high
# export SWIFT_TMP=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/tmp
# export TMPDIR=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/tmp
export SWIFT_TMP=/scratch/sciteam/$USER/variant-tmp
export TMPDIR=$SWIFT_TMP

# CRAY specific settings:
export CRAY_PPN=true

# (Optional variables to set)
export TURBINE_APP_RETRIES=6
export ADLB_SERVERS=1
export TURBINE_LOG=1    # This produces verbose logging info; great for debugging
export ADLB_DEBUG_RANKS=1	# Displays layout of ranks and nodes
# This specifies where the log info will be stored; defaults to one's home directory
# export TURBINE_OUTPUT=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/601Nodes/log_directory
export TURBINE_OUTPUT=$HOME/variant-output
export TURBINE_OUTPUT_FORMAT="X%Q"

PATH=/u/sciteam/wozniak/Public/sfw/compute/swift-t/stc/bin:$PATH

THIS=$( dirname $0 )
VARIANTS_HOME=$( cd $THIS/../../.. ; /bin/pwd )

swift-t -m cray -O3 -n $PROCS \
        -I $VARIANTS_HOME/src/ \
        -r $VARIANTS_HOME/src/bioapps \
        $VARIANTS_HOME/src/VariantCalling.swift \
        -runfile=$THIS/601Node.runfile
