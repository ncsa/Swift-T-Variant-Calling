#Data Used
##Reference
#Hg38 chromosomes 20,21,22
Homo_sapiens_assembly38_chr20_21_22.fasta

#dbSNP
dbsnp_146.hg38.chr20_21_22.vcf.gz

#Indel vcfs split by chromosome
Mills_and_1000G_gold_standard.indels.hg38.chr20.vcf, Mills_and_1000G_gold_standard.indels.hg38.chr21.vcf, Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf


#Samples
Samples were generated with NEAT read simulator based off of the 3 chromosome reference. The NEAT parameters used were:
	-r Homo_sapiens_assembly38_chr20_21_22.fasta
	-R 100 \
	--pe 300 30 \
	-M 0.001 \
	-E 0.001 \
	--bam \
	--vcf
We then used NEATs parralellizing simulation method to generate the output.
We copied the pairs of fastqs in order to run Swift with multiple samples.




