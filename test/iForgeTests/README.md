# Data Used
## Reference
We used the 3 smallest from Hg38: chromosomes 20,21,22.
Homo_sapiens_assembly38_chr20_21_22.fasta

## dbSNP
We parsed the hg38 dbSNP to just contain chromosomes 20,21,22.
dbsnp_146.hg38.chr20_21_22.vcf.gz

## Indel vcfs split by chromosome
Mills_and_1000G_gold_standard.indels.hg38.chr20.vcf, Mills_and_1000G_gold_standard.indels.hg38.chr21.vcf, Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf


## Samples
Samples were generated with NEAT read simulator. https://github.com/zstephens/neat-genreads#parallelizing-simulation
We simulated paired end reads based off of the off of the 3 chromosome reference described above. The NEAT parameters used were:

```
	-r Homo_sapiens_assembly38_chr20_21_22.fasta
	-R 100 \
	--pe 300 30 \
	-M 0.001 \
	-E 0.001 \
	--bam \
	--vcf
```
We used NEATs parallelizing simulation method to generate the output.
We then copied the pairs of fastqs in order to run Swift with multiple samples. For example, 4 samples ran were:
* hg38_chr_20_21_22_merged_copy1_read1.fq.gz
* hg38_chr_20_21_22_merged_copy1_read2.fq.gz
* hg38_chr_20_21_22_merged_copy2_read1.fq.gz
* hg38_chr_20_21_22_merged_copy2_read2.fq.gz
* hg38_chr_20_21_22_merged_copy3_read1.fq.gz
* hg38_chr_20_21_22_merged_copy3_read2.fq.gz
* hg38_chr_20_21_22_merged_copy4_read1.fq.gz
* hg38_chr_20_21_22_merged_copy4_read2.fq.gz





