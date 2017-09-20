#!/bin/bash
set -x

# This file takes care of all the steps necessary to generate a synthetic dataset of number of samples specified by (NumberOfSamples)
 
## set some paths
# Tools:
NEAT2_path="/home1/04525/tg838247/software/neat-genreads"
samtools_path="/home1/04525/tg838247/software/samtools/samtools-1.5"
reference="/work/04525/tg838247/stampede2/swift_T_project/data/genome/gatk-bundle/b37/human_g1k_v37.fasta"

# Simulation-specific
targeted_region=/work/04525/tg838247/stampede2/swift_T_project/data/truseq-exome-targeted-regions-manifest-v1-2.bed

OutputDir=/work/04525/tg838247/stampede2/swift_T_project/data/reads/synthetic_reads

rm -rf $OutputDir

DatasetName=synthetic_WES_Sep_2017

NumberOfChunks=1

NumberOfSamples=3

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
	OutputFolder=$OutputDir/${i}
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
	 
	echo "#!/bin/bash" > ${OutputLogsFolder}/GenericSbatch.header
	echo "#SBATCH -p normal" >> ${OutputLogsFolder}/GenericSbatch.header
	echo "#SBATCH -c 5" >> ${OutputLogsFolder}/GenericSbatch.header
	echo "#SBATCH -n 30" >>  ${OutputLogsFolder}/GenericSbatch.header
	echo "#SBATCH --time=24:00:00" >> ${OutputLogsFolder}/GenericSbatch.header
	echo "#SBATCH --mail-user=azzaea@gmail.com" >> ${OutputLogsFolder}/GenericSbatch.header
	echo "#SBATCH --mail-type=ALL" >> ${OutputLogsFolder}/GenericSbatch.header
	 
	#### generate and submit the chunks jobs
	####
	AllChunksJobIds=""
	for (( jobnumber=1; jobnumber<=${NumberOfChunks}; jobnumber++ ))
	do
	   cat ${OutputLogsFolder}/GenericSbatch.header > ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.sbatch
	   echo "#SBATCH -J ${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}" >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.submit
	   echo "#SBATCH -e ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.er" >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.sbatch
	   echo "#SBATCH -o ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.ou" >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.sbatch

	   echo " " >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}
		
	   echo "module load python" >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.sbatch	  
 
	   echo "python -u ${NEAT2_path}/genReads.py -r ${reference} -o ${OutputChunksFolder}/${DatasetName} ${NeatParams} --job ${jobnumber} ${NumberOfChunks} >  ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.log " >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.sbatch
 
	   # ChunkJobId is a temporary variable, used to fill out the total list of all job ids
	   ChunkJobId=`sbatch ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.sbatch`
	   AllChunksJobIds=${AllChunksJobIds}":"${ChunkJobID}
	done
 
 
	#### generate and submit the merge job
	####
	cat ${OutputLogsFolder}/GenericSbatch.header > ${OutputLogsFolder}/${DatasetName}.MergeChunks.sbatch
	echo "#SBATCH -J ${DatasetName}.MergeChunks" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.sbatch
	echo "#SBATCH -e ${OutputLogsFolder}/${DatasetName}.MergeChunks.er" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.sbatch
	echo "#SBATCH -o ${OutputLogsFolder}/${DatasetName}.MergeChunks.ou" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.sbatch
	echo "#SBATCH -depend=afterok:${AllChunksJobIds}" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.sbatch
 
	echo " " >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.sbatch

	echo "module load python" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.sbatch
 
	echo "python ${NEAT2_path}/mergeJobs.py -i ${OutputChunksFolder}/${DatasetName} -o ${OutputMergedFolder}/${DatasetName}_merged -s  ${samtools_path}/samtools" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.sbatch
 
	sbatch ${OutputLogsFolder}/${DatasetName}.MergeChunks.sbatch


done

