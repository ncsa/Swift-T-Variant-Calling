#!/usr/bin/env bash

export PPN=1
export NODES=601
export PROCS=$(($PPN * $NODES))
export WALLTIME=00:10:00
# export PROJECT=baib
export PROJECT=bahk
export QUEUE=high
# export QUEUE=debug

SCRATCH=/scratch/sciteam/$USER
DATA=$SCRATCH/601Node_tests

echo SCRATCH=$SCRATCH

# export SWIFT_TMP=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/tmp
# export TMPDIR=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/tmp
export SWIFT_TMP=$DATA/tmp
export TMPDIR=$DATA/tmp

# CRAY specific settings:
export CRAY_PPN=true

# (Optional variables to set)
export TURBINE_APP_RETRIES=6
# export ADLB_SERVERS=0
export TURBINE_LOG=1    # This produces verbose logging info; great for debugging
export ADLB_DEBUG_RANKS=1	# Displays layout of ranks and nodes

export TURBINE_OUTPUT_ROOT=$SCRATCH/turbine-output
export TURBINE_OUTPUT_FORMAT=%Q

# export TURBINE_OUTPUT=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/Justins_Testing/601Node_tests/1samplePerNode_log_directory	# This specifies where the log info will be stored; defaults to one's home directory

# PATH=/u/sciteam/wozniak/Public/sfw/compute/swift-t/stc/bin:$PATH
PATH=/u/sciteam/wozniak/Public/sfw/compute/swift-t-2018-06-05/stc/bin:$PATH
which swift-t

# STVC = Swift/T Variant Calling main repo root
STVC=$( readlink -f $( dirname $0 )/../../../.. )

echo STVC=$STVC

swift-t -m cray -O3 -n $PROCS \
	-I $STVC/src/ \
	-r $STVC/src/bioapps \
        -e SWIFT_TMP \
	$STVC/src/VariantCalling.swift \
	-runfile=$STVC/test/BWTests/justin/601Node/601Node_1samplePerNode.runfile
