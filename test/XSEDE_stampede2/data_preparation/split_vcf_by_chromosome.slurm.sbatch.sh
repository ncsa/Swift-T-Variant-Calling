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

#SBATCH -J SplitByChrom 
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 10:30:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=azzaea@gmail.com
#SBATCH --mail-type=all    # Send email at begin and end of job

# Other commands must follow all #SBATCH directives...

bwa=/home1/04525/tg838247/software/bwa/bwa-0.7.16/bwa
novocraft=/home1/04525/tg838247/software/novocraft/novocraft-3.08.00/
samtools=/home1/04525/tg838247/software/samtools/samtools-1.5/samtools
picardjar=/home1/04525/tg838247/software/picard/picard-2.10.5/picard.jar
java=/bin/java
gatkjar=/home1/04525/tg838247/software/gatk/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
refdir=/work/04525/tg838247/stampede2/swift_T_project/data/genome/gatk-bundle/b37
# Launch serial code...

set -x

cd $refdir
/home1/04525/tg838247/software/htslib/htslib-1.5/tabix -fp vcf Mills_and_1000G_gold_standard.indels.b37.vcf.gz
/home1/04525/tg838247/software/htslib/htslib-1.5/tabix -fp vcf 1000G_phase1.indels.b37.vcf.gz
for i in `seq 1 22` ; do
        $java -jar $gatkjar\
                -R ${refdir}/human_g1k_v37.fasta\
                -T SelectVariants \
                -V ${refdir}/1000G_phase1.indels.b37.vcf.gz\
                -L ${i} \
                -o ${refdir}/IndelsPerChrom/1000G.${i}.vcf
            
        echo "splitting the 1000G file;       chr: ${i};         exit with $?"

        $java -jar $gatkjar\
                -R ${refdir}/human_g1k_v37.fasta\
                -T SelectVariants \
                -V ${refdir}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz\
                -L ${i} \
                -o ${refdir}/IndelsPerChrom/Mills.${i}.vcf

        echo "splitting the Mills file;       chr: ${i};         exit with $?"
done

for i in X Y MT ; do
        $java -jar $gatkjar \
                -R ${refdir}/human_g1k_v37.fasta\
                -T SelectVariants \
                -V ${refdir}/1000G_phase1.indels.b37.vcf.gz\
                -L ${i} \
                -o ${refdir}/IndelsPerChrom/1000G.${i}.vcf
    
        echo "splitting the 1000G file;       chr: ${i};         exit with $?"

        $java -jar $gatkjar \
                -R ${refdir}/human_g1k_v37.fasta \
                -T SelectVariants \
                -V ${refdir}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz\
                -L ${i} \
                -o ${refdir}/IndelsPerChrom/Mills.${i}.vcf

        echo "splitting the Mills file;       chr: ${i};         exit with $?"
done


echo "Finishing the splitting job :)"

