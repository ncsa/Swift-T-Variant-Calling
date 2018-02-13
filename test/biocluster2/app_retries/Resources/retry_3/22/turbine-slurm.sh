#!/bin/bash -l
# We changed the M4 comment to d-n-l, not hash
# We need 'bash -l' for the module system

# Copyright 2013 University of Chicago and Argonne National Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License

# TURBINE-SLURM.SH

# Created: Sun 14 Jan 05:10:22 CST 2018


# Define convenience macros
# This simply does environment variable substition when m4 runs



# Other key settings are on the sbatch command line
# See turbine-slurm-run.zsh
#SBATCH --time=48:00:00
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=3
#SBATCH --workdir=/home/a-m/azzaea/turbine-output/2018/01/14/05/10/22

# M4 conditional to optionally perform user email notifications
#SBATCH --mail-user=azzaea@gmail.com
#SBATCH --mail-type=ALL


# User directives:


echo TURBINE-SLURM.SH

export TURBINE_HOME=$( cd "$(dirname "$0")/../../.." ; /bin/pwd )

VERBOSE=0
if (( ${VERBOSE} ))
then
 set -x
fi

TURBINE_HOME=/home/apps/software/swift-t/1.3-IGB-gcc-4.9.4/turbine
source ${TURBINE_HOME}/scripts/turbine-config.sh

COMMAND="/home/apps/software/Tcl/8.6.6-IGB-gcc-4.9.4/bin/tclsh8.6 /home/a-m/azzaea/turbine-output/2018/01/14/05/10/22/swift-t-retries.lIx.tic "

# Use this on Midway:
# module load openmpi gcc/4.9

${TURBINE_LAUNCHER} ${COMMAND}
# Return exit code from mpirun
