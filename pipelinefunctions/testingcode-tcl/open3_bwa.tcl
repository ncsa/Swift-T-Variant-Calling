#!/usr/bin/env tclsh
package require Tcl 8.6

set index /home/azza/swift-project/Workshop_low_pass/ref/human_g1k_v37_chr20.fa
set R1 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_1.fastq
set R2 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_2.fastq
set rgheader  {@RG\tID:synthetic\tLB:synthetic\tPL:illumina\tPU:synthetic\tSM:synthetic\tCN:synthetic}
set bwadir {/usr/bin/bwa}
set samtoolsdir {/usr/local/bin/samtools}

set readchan [open "|$bwadir mem $index $R1 $R2 -R $rgheader"  r+ ]
set writechan [open "tmp.sam" w+]
fconfigure $readchan -buffering none;  
fileevent $readchan readable  [puts $writechan [read $readchan]]

set readchan2 [open "|$samtoolsdir view -b tmp.sam " r+]
set writechan2 [open "tmp.bam" w+]
fconfigure $readchan2 -buffering none
fileevent $readchan2 readable  [puts $writechan2 [read $readchan2]]

#lassign [chan pipe] readchan writechan

close $writechan
close $writechan2

