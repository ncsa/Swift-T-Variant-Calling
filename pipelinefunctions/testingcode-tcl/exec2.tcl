#!/usr/bin/env tclsh

set index /home/azza/swift-project/Workshop_low_pass/ref/human_g1k_v37_chr20.fa
set R1 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_1.fastq
set R2 /home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_2.fastq
set rgheader  {@RG\tID:synthetic\tLB:synthetic\tPL:illumina\tPU:synthetic\tSM:synthetic\tCN:synthetic}

set bwadir {/usr/bin/bwa}
set samtoolsdir {/usr/local/bin/samtools}

proc bwa {bwadir index R1 R2 rgheader args} {
	set result "$bwadir mem $index $R1 $R2 -R $rgheader"
	return $result
}

proc samtools {samtoolsdir args} {
	set result "$samtoolsdir view -b $args"
	return $result
}

set b [bwa $bwadir $index $R1 $R2 $rgheader]
set s [ samtools $samtoolsdir]
exec {*}$b | {*}$s   > tmp.bam 


