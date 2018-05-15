#!/usr/bin/env bash

export PPN=2
export NODES=601
export PROCS=$(($PPN * $NODES))
export WALLTIME=7:00:00
export PROJECT=baib
export QUEUE=high
export SWIFT_TMP=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/tmp
export TMPDIR=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/tmp

# CRAY specific settings:
export CRAY_PPN=true

# (Optional variables to set)
export TURBINE_APP_RETRIES=6
export ADLB_SERVERS=1
export TURBINE_LOG=1    # This produces verbose logging info; great for debugging
export ADBL_DEBUG_RANKS=1	# Displays layout of ranks and nodes
export TURBINE_OUTPUT=/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/Justins_Testing/601Node_tests/2samplesPerNode_log_directory	# This specifies where the log info will be stored; defaults to one's home directory

PATH=/u/sciteam/wozniak/Public/sfw/compute/swift-t/stc/bin:$PATH

swift-t -m cray -O3 -n $PROCS -o /scratch/sciteam/jacobrh/purge_exempt/Swift_testing/compiled.tic \
	-I /projects/sciteam/baib/SwiftTesting/Swift-T-Variant-Calling/src/ \
	-r /projects/sciteam/baib/SwiftTesting/Swift-T-Variant-Calling/src/bioapps \
	/projects/sciteam/baib/SwiftTesting/Swift-T-Variant-Calling/src/VariantCalling.swift \
	-runfile=/projects/sciteam/baib/SwiftTesting/Swift-T-Variant-Calling/test/BWTests/Justins_Testing/601Node_tests/601Node_2samplesPerNode.runfile
