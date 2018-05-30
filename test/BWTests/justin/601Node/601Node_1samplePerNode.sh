#!/usr/bin/env bash

export PPN=1
export NODES=2
export PROCS=$(($PPN * $NODES))
export WALLTIME=00:05:00
export PROJECT=baib
# export QUEUE=high
export QUEUE=debug

TEST=/scratch/sciteam/wozniak/601Node_tests

# export SWIFT_TMP=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/tmp
# export TMPDIR=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/tmp
export SWIFT_TMP=$TEST/tmp
export TMPDIR=$TEST/tmp

# CRAY specific settings:
export CRAY_PPN=true

# (Optional variables to set)
export TURBINE_APP_RETRIES=6
# export ADLB_SERVERS=0
export TURBINE_LOG=0    # This produces verbose logging info; great for debugging
export ADLB_DEBUG_RANKS=1	# Displays layout of ranks and nodes

export TURBINE_OUTPUT_ROOT=/scratch/sciteam/wozniak/turbine-output
export TURBINE_OUTPUT_FORMAT=%Q

# export TURBINE_OUTPUT=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/Justins_Testing/601Node_tests/1samplePerNode_log_directory	# This specifies where the log info will be stored; defaults to one's home directory

PATH=/u/sciteam/wozniak/Public/sfw/compute/swift-t/stc/bin:$PATH

# STVC = Swift/T Variant Calling
STVC=$HOME/proj/Variants

# /projects/sciteam/baib/SwiftTesting/Swift-T-Variant-Calling

swift-t -m cray -O3 -n $PROCS \
	-I $STVC/src/ \
	-r $STVC/src/bioapps \
	$STVC/src/VariantCalling.swift \
	-runfile=$TEST/601Node_1samplePerNode.runfile
