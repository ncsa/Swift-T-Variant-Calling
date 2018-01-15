#!/bin/bash

module load swift-t

set -x

swift-t -l -m slurm -O3 -l -V -s settings.sh tst.swift | tee $PWD/log.swift_t_run

