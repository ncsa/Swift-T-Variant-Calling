# Data Used
## Reference
We used the 3 smallest chromosomes from Hg38: chromosomes 20,21,22.
Location of 3 chromsome reference on iForge:

`/projects/bioinformatics/HudsonSoybeanProject/swift_T_variant_calling_updated/Swift-T-Variant-Calling/test/iForgeTests/reference/Homo_sapiens_assembly38_chr20_21_22.fasta`

## dbSNP
We parsed the hg38 dbSNP to just contain chromosomes 20,21,22.
Location of dbSNP on iForge:

`/projects/bioinformatics/HudsonSoybeanProject/swift_T_variant_calling_updated/Swift-T-Variant-Calling/test/iForgeTests/reference/dbNSP/dbsnp_146.hg38.chr20_21_22.vcf.gz`

## Indel vcfs split by chromosome
We parsed the Mills and 1000G gold standard indels vcf for chromosomes 20,21,22. The workflow requires the indels vcf to be split by chromosome and in the same directory.
Location of indels directory on iForge:

`/projects/bioinformatics/HudsonSoybeanProject/swift_T_variant_calling_updated/Swift-T-Variant-Calling/test/iForgeTests/reference/indels`

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
We used NEATs parallelizing simulation method to generate the output. NEAT will output a 'golden' bam and 'golden' vcf file (refer to NEAT documentation). 

The output from the NEAT simulation is in:

`/projects/bioinformatics/HudsonSoybeanProject/swift_T_variant_calling_updated/Swift-T-Variant-Calling/test/iForgeTests/NEAT_simulated_data_anisimov_launcher/anisimov_scripts/Merged`

We then copied the pairs of fastqs in order to run Swift with multiple samples.  For example, 4 samples ran were:
```
hg38_chr_20_21_22_merged_copy1_read1.fq.gz
hg38_chr_20_21_22_merged_copy1_read2.fq.gz
hg38_chr_20_21_22_merged_copy2_read1.fq.gz
hg38_chr_20_21_22_merged_copy2_read2.fq.gz
hg38_chr_20_21_22_merged_copy3_read1.fq.gz
hg38_chr_20_21_22_merged_copy3_read2.fq.gz
hg38_chr_20_21_22_merged_copy4_read1.fq.gz
hg38_chr_20_21_22_merged_copy4_read2.fq.gz
```




