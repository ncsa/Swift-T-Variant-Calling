#!/usr/bin/env tclsh
set index /home/azza/swift-project/Workshop_low_pass/ref/human_g1k_v37_chr20.fa
set R1 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_1.fastq
set R2 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_2.fastq
set rgheader  {@RG\tID:synthetic\tLB:synthetic\tPL:illumina\tPU:synthetic\tSM:synthetic\tCN:synthetic}
set bwadir {/usr/bin/bwa}
set samtoolsdir {/usr/local/bin/samtools}
set readchan [open "|$bwadir mem $index $R1 $R2 -R $rgheader"  r+ ]
set writechan [open "tmp.bam" w+]
fconfigure $readchan -buffering none;  
set readchan2 [open "|$samtoolsdir view -b tmp.sam " r+] ;# for this to work, tmp.sam must be created. This is really against the idea of creating pipes in the first place!!!!
fconfigure $readchan2 -buffering none

fileevent $readchan readable  [puts $writechan [
	fileevent $readchan2 readable  [puts $writechan2 [read $readchan2]]
						]]

#lassign [chan pipe] readchan writechan
close $readchan
close $writechan
close $readchan2 
close $writechan2

