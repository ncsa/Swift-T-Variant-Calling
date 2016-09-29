/////// Alignment functions:
@dispatch=WORKER
app (file output) bwa_mem (string bwadir, string read1, string read2, string INDEX, string bwamemparams, int PBSCORES,  string rgheader){
	bwadir "mem" "-M" bwamemparams "-t" PBSCORES "-R" rgheader INDEX read1 read2 @stdout=output;
}

@dispatch=WORKER
app (file output) novoalign (string novoaligndir, string read1, string read2, string INDEX, string novoalignparams, int PBSCORES, string rgheader) {
	novoaligndir "-c" PBSCORES "-d" INDEX "-f" read1 read2 "-o" "SAM" rgheader @stdout=output; 
}

@dispatch=WORKER
app (file output) samtools_view(string samtoolsdir, file inputFile, int thr, string u){
	samtoolsdir "view" "-@" thr "-bS" inputFile u @stdout=output;
	/// samtools view has dual functions, as it also counts aligns!!
	/// samtools view -c bamfile
}

@dispatch=WORKER
app (file output) samblaster(string samblasterdir, file inputFile){
	samblasterdir "-M" "-i" inputFile @stdout=output;
}


@dispatch=WORKER
app (file output) novosort (string novosortdir, string inputFilename, string tmpdir, int thr, string sortoptions){
	novosortdir "--index" "--tmpdir" tmpdir "--threads" thr inputFilename "-o" output; 
	// novosort has dual function to also mark duplicates
}
app (file output) novosort (string novosortdir, string inputFilename[], string tmpdir, int thr, string sortoptions){
	novosortdir "--index" "--tmpdir" tmpdir "--threads" thr inputFilename "-o" output; 
	// novosort has dual function to also mark duplicates
}

@dispatch=WORKER
app (file outputfile, file metricsfile) picard (string javadir, string picarddir, string tmpdir, file inputFile ){
        javadir "-Xmx8g" "-jar" picarddir "MarkDuplicates" "INPUT=" inputFile "OUTPUT=" outputfile "METRICS_FILE=" metricsfile "TMP_DIR=" tmpdir "ASSUME_SORTED=true MAX_RECORDS_IN_RAM=null CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT";

// This is the command line that is working:
// $javadir -Xmx8g -jar $picarddir MarkDuplicates INPUT=/home/azza/swift-project/Results/align/HG00108.lowcoverage.chr20.smallregion.nodups.bam
}

//////////////////////////////////////// Testing and validation functions:

@dispatch=WORKER
app (file output) samtools_flagstat(string samtoolsdir, file inputFile){
	samtoolsdir "flagstat" inputFile @stdout=output;
}

