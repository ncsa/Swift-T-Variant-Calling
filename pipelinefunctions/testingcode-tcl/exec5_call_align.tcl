#!/usr/bin/env tclsh

lappend auto_path [pwd]

package require align 0.0

set index /home/azza/swift-project/Workshop_low_pass/ref/human_g1k_v37_chr20.fa
set R1 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_1.fastq
set R2 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_2.fastq
set rgheader  {@RG\tID:synthetic\tLB:synthetic\tPL:illumina\tPU:synthetic\tSM:synthetic\tCN:synthetic}

set bwadir {/usr/bin/bwa}
set samtoolsdir {/usr/local/bin/samtools}

set b [bwa $bwadir $index $R1 $R2 $rgheader]
set s [ samtools $samtoolsdir]

exec {*}$b | {*}$s > aligned.bam

#pipe {aligned.bam} $b $s


