#!/bin/bash
#----------------------------------------------------
# Sample SLURM job script
#   for TACC Stampede2 KNL nodes
#
#   *** Serial Job on Normal Queue ***
#
# Last revised: 27 Jun 2017
#
# Notes:
#
#   -- Copy/edit this script as desired.  Launch by executing
#      "sbatch knl.serial.slurm" on a Stampede2 login node.
#
#   -- Serial codes run on a single node (upper case N = 1).
#        A serial code ignores the value of lower case n,
#        but slurm needs a plausible value to schedule the job.
#
#   -- For a good way to run multiple serial executables at the
#        same time, execute "module load launcher" followed
#        by "module help launcher".

#----------------------------------------------------

#SBATCH -J PrepRefGenome           # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 12:30:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=azzaea@gmail.com
#SBATCH --mail-type=all    # Send email at begin and end of job

# Other commands must follow all #SBATCH directives...

bwa=/home1/04525/tg838247/software/bwa/bwa-0.7.16/bwa
novocraft=/home1/04525/tg838247/software/novocraft/novocraft-3.08.00/
samtools=/home1/04525/tg838247/software/samtools/samtools-1.5/samtools
picardjar=/home1/04525/tg838247/software/picard/picard-2.10.7/picard.jar
java=/bin/java
gatkjar=/home1/04525/tg838247/software/gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
refdir=/work/04525/tg838247/stampede2/swift_T_project/data/genome/gatk-bundle/b37
# Launch serial code...

set -x

cd $refdir

#gunzip human_g1k_v37.fasta.gz
#gunzip human_g1k_v37.fasta.fai.gz
#gunzip human_g1k_v37.dict.gz

#$bwa index -a bwtsw human_g1k_v37.fasta 
#$novocraft/novoindex human_g1k_v37.nix human_g1k_v37.fasta 
#$samtools faidx human_g1k_v37.fasta

$java -jar $picardjar CreateSequenceDictionary \
	REFERENCE=human_g1k_v37.fasta \
	OUTPUT=human_g1k_v37.dict 

# ---------------------------------------------------


