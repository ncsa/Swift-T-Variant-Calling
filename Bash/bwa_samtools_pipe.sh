#!/bin/bash

############ BWA - samtools pipe:
############# This script does the following:
# $bwamemdir/bwa mem $bwamem_parms -t $thr -R "${rgheader}" $bwa_index $R1 $R2 | $samtoolsdir/samtools view -@ $thr -bSu -> $alignedbam 
# So, it needs the following input parameters (in this very order when called):
# parameter name	parameter order	parameter meaning
# $bwamem_parms 	1		bwa specific parameters (from run file)
# $thr 		2		number of threads (pbs cores)
# "${rgheader}" 	3		read group header
#$bwa_index 	4		bwa indexed reference genome file
#$R1 		5		read1 file	
#$R2 		6		read2 file	
#$alignedbam 	7		name of the output, aligned bam file
############################################################################################################

bwa mem $1 -t $2 -R $3 $4 $5 $6 | samtools view -@ $2 -bSu -> $7

