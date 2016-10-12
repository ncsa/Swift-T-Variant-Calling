#!/usr/bin/env tclsh

set index /home/azza/swift-project/Workshop_low_pass/ref/human_g1k_v37_chr20.fa
set R1 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_1.fastq
set R2 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_2.fastq
set rgheader  {@RG\tID:synthetic\tLB:synthetic\tPL:illumina\tPU:synthetic\tSM:synthetic\tCN:synthetic}
set bwadir {/usr/bin/bwa}
set samtoolsdir {/usr/local/bin/samtools}

#exec /usr/bin/bwa mem $index $R1 $R2 | /usr/local/bin/samtools view -b > /home/azza/aligned.bam

proc bwa {bwaDir Index r1 r2 rgheader args} {
	exec $bwaDir mem $Index $r1 $r2 -R $rgheader $args >tmp.sam 
}

proc samtools_view {samtoolsdir inputfile thr args} {
	exec $samtoolsdir view -@ $thr -bS $inputfile $args >tmp.bam
}

eval samtools_view $samtoolsdir [eval bwa $bwadir $index $R1 $R2 $rgheader -t 2 -M ] 2 

##### Working with pipes:
#set pipe [open "|/usr/bin/bwa mem $index $R1 $R2 -R $rgheader tmp.sam" RDWR ]

#fconfigure $pipe -blocking 0 -buffering none; #
#fileevent $pipe readable [puts "bwa is done!"]


