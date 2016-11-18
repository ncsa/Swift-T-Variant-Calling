#!/usr/bin/env tclsh

set index /home/azza/swift-project/Workshop_low_pass/ref/human_g1k_v37_chr20.fa
set R1 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_1.fastq
set R2 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_2.fastq
set rgheader  {@RG\tID:synthetic\tLB:synthetic\tPL:illumina\tPU:synthetic\tSM:synthetic\tCN:synthetic}
set bwadir {/usr/bin/bwa}
set samtoolsdir {/usr/local/bin/samtools}

if { [catch {exec -keepnewline -- $bwadir mem $index $R1 $R2 -R $rgheader > tmp.sam} msg] } {
#	puts "bwa run details"
#	puts "Information about bwa run are: $::errorInfo"
}
 #exec -keepnewline -- $bwadir mem $index $R1 $R2 -R $rgheader > tmp.sam 

#exec $bwadir mem $index $R1 $R2 -R $rgheader > tmp.sam; #| $samtoolsdir view -b > tmp.bam 
