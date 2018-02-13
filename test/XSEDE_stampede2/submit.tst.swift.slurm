#!/bin/bash

#SBATCH -J tstSwift 
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 2               # Total # of nodes (must be 1 for serial)
#SBATCH -n 5               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 02:30:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=azzaea@gmail.com
#SBATCH --mail-type=all    # Send email at begin and end of job

# Other commands must follow all #SBATCH directives...


export TURBINE_LOG=1
export  ADLB_DEBUG_RANKS=1
export ADLB_DEBUG_HOSTMAP=1

swift-t -l -n 2 tst.swift.2
