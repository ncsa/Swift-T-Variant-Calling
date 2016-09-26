#!/usr/bin/env tclsh

set index /home/azza/swift-project/Workshop_low_pass/ref/human_g1k_v37_chr20.fa
set R1 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_1.fastq
set R2 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_2.fastq

#exec /usr/bin/bwa mem $index $R1 $R2 | /usr/local/bin/samtools view -b > /home/azza/aligned.bam

#proc bwa {Index r1 r2 rgheader args} {
#	exec /usr/bin/bwa mem $Index $r1 $r2 -R $rgheader $args >tmp.sam 
#}

#proc samtools_view {thr file args} {
#	exec /usr/local/bin/samtools view -@ $thr -bS $args file >tmp.bam
#}


set rgheader  "@RG\tID:synthetic\tLB:synthetic\tPL:illumina\tPU:synthetic\tSM:synthetic\tCN:synthetic"
#eval bwa $index $R1 $R2 {@RG\tID:synthetic\tLB:synthetic\tPL:illumina\tPU:synthetic\tSM:synthetic\tCN:synthetic} -t 2 -M

set pipe [open "|/usr/bin/bwa mem $index $R1 $R2 -R $rgheader tmp.sam" RDWR ]

fconfigure $pipe -blocking 0 -buffering none; #
fileevent $pipe readable [puts "bwa is done!"]


