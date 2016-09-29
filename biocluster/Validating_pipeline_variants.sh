#!/bin/bash
#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o /home/groups/hpcbio_shared/azza/GIAB/src/log.hap.py.ou
#PBS -e /home/groups/hpcbio_shared/azza/GIAB/src/log.hap.py.e 
#PBS -M aeahmed@illinois.edu
#PBS -m abe

######### Paths defintions:
set -x
output=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/results/run2
goldenFile=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/reads/synthetic_rare_variants/sample.1/Chunks/synthetic_rare_variants_golden.vcf
workflowFile=$output/delivery/jointVCFs/jointVCFcalled.vcf

export HGREF=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/genome/HG19_GATKbundle2.8_noDecoys.fa

########################### Preparatory stages:
set +x
module load tabix

set -x
bgzip -c $goldenFile > $goldenFile.gz
tabix -p vcf $goldenFile.gz

bgzip -c $workflowFile > $workflowFile.gz
tabix -p vcf $workflowFile.gz

########################### Comparison stage:
set +x
module load hap.py/0.3.0

set -x
hap.py $goldenFile.gz $workflowFile.gz -o $output/variant_compare_hap.py \
 --logfile $output/hap.py.log \
 --no-internal-leftshift \
  --no-internal-preprocessing \
  --include-nonpass-truth -V \
  --roc Q_GQ \
  --roc-filter LowGQX
