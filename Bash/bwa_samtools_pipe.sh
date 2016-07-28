#!/bin/bash

# BWA - samtools pipe:
# This script does the following:
# $bwamemdir/bwa mem $bwamem_parms -t $thr -R "${rgheader}" $bwa_index $R1 $R2 | $samtoolsdir/samtools view -@ $thr -bSu -> $alignedbam 
# So, it needs the following input parameters:
# parameter name	parameter number
# $bwamem_parms 	1
# $thr 		2
# "${rgheader}" 	3
#$bwa_index 	4
#$R1 		5
#$R2 		6
$alignedbam 	7

bwa mem $1 -t $2 -R $3 $4 $5 $6 | samtools view -@ $2 -bSu -> $7

