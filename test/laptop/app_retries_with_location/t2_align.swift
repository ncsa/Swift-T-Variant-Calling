import location;
import string;

@dispatch=WORKER
app (file output) bwa_mem (string bwaexe, string read1, string read2, string INDEX, string bwamemparams[], int PBSCORES,  string rgheader){
        bwaexe "mem" "-M" bwamemparams "-t" PBSCORES "-R" rgheader INDEX read1 read2 @stdout=output ;
}




string read1= strcat("/home/azza/swift-project/data/reads/HG00108.lowcoverage.chr20.smallregion_1.fastq");

string read2= strcat("/home/azza/swift-project/data/reads/HG00108.lowcoverage.chr20.smallregion_2.fastq");

string bwaindex = strcat("/home/azza/swift-project/data/genome/human_g1k_v37_chr20.fa");

string rgheader = sprintf("@RG\tID:%s\tLB:%s\tPL:%s\tPU:%s\tSM:%s\tCN:%s", "synthetic", "illumina", "illumina", "illumina", "illumina","SAMPLECN");
string bwaparams[] = ["-k 12"];


file alignedsam < strcat("/home/azza/swift-project/tmp", "/align/", "test", ".nodups.sam") >;

alignedsam = @location=locationFromRank(3) bwa_mem("/usr/src/bwa/bwa-0.7.15/bwa", read1, read2, bwaindex, bwaparams, 3, rgheader                         );


