#!/bin/bash

module load swift-t

set -x

swift-t -m slurm -O3 -l -V -s settings.sh retries.swift | tee $PWD/log.swift_t_run

