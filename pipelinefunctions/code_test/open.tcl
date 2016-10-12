#!/usr/bin/env tclsh

set index /home/azza/swift-project/Workshop_low_pass/ref/human_g1k_v37_chr20.fa
set R1 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_1.fastq
set R2 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_2.fastq
set rgheader  {@RG\tID:synthetic\tLB:synthetic\tPL:illumina\tPU:synthetic\tSM:synthetic\tCN:synthetic}
set bwadir {/usr/bin/bwa}
set samtoolsdir {/usr/local/bin/samtools}


set alignedsam [open "|/usr/bin/bwa mem $index $R1 $R2 " w+  ]

fconfigure $alignedsam -buffering none -blocking 0;

fileevent $alignedsam readable [puts [open "tmp.sam" w+ ] $alignedsam]

#
#open  "|/usr/local/bin/samtools view -b" ]
#gets $alignedsam #

close $alignedsam

# > /home/azza/aligned.bam

##### Working with pipes:
#set pipe [open "|/usr/bin/bwa mem $index $R1 $R2 -R $rgheader tmp.sam" RDWR ]

