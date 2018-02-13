#!/bin/bash

# This is a scirpt to run a swift/t script that is passed in as the first argument. It merely reduces the burden of calling the module and passing the settings file!

module load swift-t

set -x

swift-t -m slurm -O3 -l -V -s settings.sh $1 

