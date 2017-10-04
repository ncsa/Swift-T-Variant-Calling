#!/bin/bash
#SBATCH -J PrepRefGenome           # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 20               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 12:30:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=azzaea@gmail.com
#SBATCH --mail-type=all    # Send email at begin and end of job



module load novocraft
set -x
mkdir -p /home/groups/hpcbio_shared/azza/swift_T_project/data/genome

novoindex /home/groups/hpcbio_shared/azza/swift_T_project/data/genome/human.nix /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome/ucsc.hg19.fasta
