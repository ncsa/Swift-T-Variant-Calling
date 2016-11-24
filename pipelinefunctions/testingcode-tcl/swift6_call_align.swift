// call this using: swift-t -r $PWD swift6_call_align.swift

string bwadir = "/usr/bin/bwa";
string samtoolsdir = "/usr/local/bin/samtools";
string index = "/home/azza/swift-project/Workshop_low_pass/ref/human_g1k_v37_chr20.fa";
string R1 = "/home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_1.fastq";
string R2 = "/home/azza/swift-project/Dataset/HG00108.lowcoverage.chr20.smallregion_2.fastq";
string rgheader = "{@RG\tID:synthetic\tLB:synthetic\tPL:illumina\tPU:synthetic\tSM:synthetic\tCN:synthetic}";

(string out) bwa (string bwadir, string index, string R1, string R2, string rgheader) "align" "0.2" [
	"set <<out>> [ alignment::bwa <<bwadir>> <<index>> <<R1>> <<R2>> <<rgheader>>]"
];

(string out) samtools_view (string samtoolsdir, int|float|string|boolean... args) "align" "0.2" "samtools" [
	"set <<out>> [ alignment::samtools <<samtoolsdir>> <<args>> ]"
];

//proc pipe {out proc1 proc2 } {           
  //     exec {*}$proc1 | {*}$proc2  > $out
//}

(file output) pipe (string proc1, string proc2, string out) "align" "0.2" [
	" set <<output>> [alignment::pipe <<out>> <<proc1>> <<proc2>> ] " ];


string b = bwa(bwadir, index, R1, R2, rgheader);
string s =  samtools_view(samtoolsdir);

//trace("\n\n" + b + "\n\n"); //trace("\n\n" + s + "\n\n");

file out <"piped_ouput.bam">;
out = pipe(b,s,"piped.bam");

// comment: The piping part requires further tuning!

