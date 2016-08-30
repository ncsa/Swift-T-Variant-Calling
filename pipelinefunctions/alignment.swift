/////// Alignment functions:
@dispatch=WORKER
app (file output) bwa (string bwadir, string read1, string read2, string INDEX, string bwamemparams, int PBSCORES,  string rgheader){
	bwadir "mem" "-M" bwamemparams "-t" PBSCORES "-R" rgheader  INDEX read1 read2 @stdout=output;
}

@dispatch=WORKER
app (file output) novoalign (string novoaligndir, string read1, string read2, string INDEX, string novoalignparams, int PBSCORES, string rgheader) {
	novoaligndir "-c" PBSCORES "-d" INDEX "-f" read1 read2 "-o" "SAM" rgheader @stdout=output; // if empty, novoalignparams causes problems,
}

@dispatch=WORKER
app (file output) samtools_view(string samtoolsdir, file inputFilename, int thr, string u){
	samtoolsdir "view" "-@" thr "-bS" u inputFilename @stdout=output;
}

@dispatch=WORKER
app (file output) samblaster(string samblasterdir, string inputFilename){
	samblasterdir "-M" inputFilename @stdout=output;
}

@dispatch=WORKER
app (file output) novosort (string novosortdir, file inputFilename, string tmpdir, int thr, string sortoptions){
	novosortdir "--index" "--tmpdir" tmpdir "--threads" thr inputFilename "-o" output  sortoptions;
}

@dispatch=WORKER
app (file outputfile, file metricsfile) picard (string javadir, string picarddir, string tmpdir, file inputFilename ){
        javadir "-jar" picarddir "MarkDuplicates"  "INPUT=" inputFilename "OUTPUT=" outputfile "METRICS_FILE=" metricsfile;

// This is the command line that is working:
// $javadir -Xmx8g -jar $picarddir MarkDuplicates INPUT=/home/azza/swift-project/Results/align/HG00108.lowcoverage.chr20.smallregion.nodups.bam
}

////////////////////////////////////////////////////////////// Testing and validation functions:

@dispatch=WORKER
app (file output) samtools_flagstat(string samtoolsdir, string inputFilename){
	samtoolsdir "flagstat" inputFilename @stdout=output;
}

