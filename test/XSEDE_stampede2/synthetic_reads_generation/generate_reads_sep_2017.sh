#!/bin/bash
#SBATCH -p normal
#SBATCH -c 5
#SBATCH -n 30
#SBATCH --time=24:00:00
#SBATCH --mail-user=azzaea@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J synthetic_reads_generation 
#SBATCH -e /home1/04525/tg838247/swift_T_project/src/Swift-T-Variant-Calling/test/XSEDE_stampede2/synthetic_reads_generation/log.err
#SBATCH -o /home1/04525/tg838247/swift_T_project/src/Swift-T-Variant-Calling/test/XSEDE_stampede2/synthetic_reads_generation/log.out 

set -x

# This file takes care of all the steps necessary to generate a synthetic dataset of number of samples specified by (NumberOfSamples)
 
## set some paths
# Tools:
NEAT2_path="/home1/04525/tg838247/software/neat-genreads"
samtools_path="/home1/04525/tg838247/software/samtools/samtools-1.5"
reference="/work/04525/tg838247/stampede2/swift_T_project/data/genome/gatk-bundle/b37/human_g1k_v37.fasta"
targeted_region=/work/04525/tg838247/stampede2/swift_T_project/data/truseq-exome-targeted-regions-manifest-v1-2.bed

OutputDir=/work/04525/tg838247/stampede2/swift_T_project/data/reads/synthetic_reads_Sep_2017

rm -rf $OutputDir

DatasetName=synthetic_WES_Sep_2017

NumberOfSamples=3
NumberOfChunks=2

echo "sample_number read_depth" > $OutputDir/simulation_parameters.txt

echo "NEAT parameters are:" >> $OutputDir/neat_parameters.txt

for i in $(seq 1 ${NumberOfSamples})
do
        # Set random error rate each time you run the program. However, we need this to be predictable:
        read_depth=$(awk -v i=$i 'BEGIN{print (i*20+10) }') 
	echo ${i} ${read_depth} >> $OutputDir/simulation_parameters.txt
	NeatParams="-R 101 -c $read_depth -p 2 -t $targeted_region --pe 300 30 --rng 100 --vcf "
	echo ${NeatParams} >> $OutputDir/neat_parameters.txt

	#create output folders
	OutputFolder=$OutputDir/sample_${i}
	OutputLogsFolder=${OutputFolder}/Logs
	if [ ! -d ${OutputLogsFolder} ]
	then
	   mkdir -p ${OutputLogsFolder}
	fi
	OutputChunksFolder=${OutputFolder}/Chunks
	if [ ! -d ${OutputChunksFolder} ]
	then
	   mkdir -p ${OutputChunksFolder}
	fi
	OutputMergedFolder=${OutputFolder}/Merged
	if [ ! -d ${OutputMergedFolder} ]
	then
	   mkdir -p ${OutputMergedFolder}
	fi
	
	####
	AllChunksJobIds=""
	for (( jobnumber=1; jobnumber<=${NumberOfChunks}; jobnumber++ ))
	do
	   echo " " >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}

	   set +x		
	   module load python
	   set -x
 
	   python -u ${NEAT2_path}/genReads.py -r ${reference} -o ${OutputChunksFolder}/${DatasetName} ${NeatParams} --job ${jobnumber} ${NumberOfChunks} 
 
	done
 
	python ${NEAT2_path}/mergeJobs.py -i ${OutputChunksFolder}/${DatasetName} -o ${OutputMergedFolder}/${DatasetName}_merged -s  ${samtools_path}/samtools
 
done

