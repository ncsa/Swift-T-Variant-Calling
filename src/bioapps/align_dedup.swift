/////// Alignment functions:
@dispatch=WORKER
app (file output, file outLog) bwa_mem (string bwaexe, string read1, string read2, string INDEX, string bwamemparams[], int PBSCORES,  string rgheader){
	bwaexe "mem" bwamemparams "-t" PBSCORES "-R" rgheader INDEX read1 read2 @stdout=output @stderr=outLog;
}

@dispatch=WORKER
app (file output, file outLog) bwa_mem (string bwaexe, string read1, string INDEX, string bwamemparams[], int PBSCORES,  string rgheader){
	bwaexe "mem" bwamemparams "-t" PBSCORES "-R" rgheader INDEX read1 @stdout=output @stderr=outLog;
}

@dispatch=WORKER
app (file output, file outLog) novoalign (string novoalignexe, string read1, string read2, string INDEX, string novoalignparams[], int PBSCORES, string rgheader) {
	novoalignexe "-c" PBSCORES "-d" INDEX "-f" read1 read2 "-o" "SAM" rgheader @stdout=output @stderr=outLog; 
}

@dispatch=WORKER                                                                                                           
app (file output, file outLog) novoalign (string novoalignexe, string read1, string INDEX, string novoalignparams[], int PBSCORES, string rgheader) {
        novoalignexe "-c" PBSCORES "-d" INDEX "-f" read1 "-o" "SAM" rgheader @stdout=output @stderr=outLog;          
}

@dispatch=WORKER
app (file output) samtools_view(string samtoolsexe, file inputFile, int thr, string args[]){
	// Converting sam to bam
	samtoolsexe "view" "-@" thr "-bS" inputFile args @stdout=output;
}

@dispatch=WORKER                                                                                                        
app (file output) samtools_bam2sam(string samtoolsexe, file inputBam, int thr, string args[]){                            
        // Converting sam to bam                                                                                        
        samtoolsexe "view" "-h" "-@" thr inputBam args @stdout=output;                                  
}              

// Counting the number of alignments
@dispatch=WORKER
(int numAlignments) samtools_view2(string samtoolsexe, string inputFile)
	"align" "0.2" [	
	"set <<numAlignments>> [ alignment::samtools_view <<samtoolsexe>> <<inputFile>> ]"
];

@dispatch=WORKER
app (file output, file outLog) samblaster(string samblasterexe, file inputSam){
	samblasterexe "-M" "-i" inputSam @stdout=output @stderr=outLog;
}

@dispatch=WORKER
app (file output, file outLog) novosort (string novosortexe, file inputFile, string tmpdir, int thr, string sortoptions[], int memFlag){
	// processing a single file (sorting and indexing input)
	novosortexe "--index" "-m" memFlag "--tmpdir" tmpdir "--threads" thr inputFile "-o" output sortoptions @stderr=outLog; 
	// novosort has dual function to also mark duplicates
}

@dispatch=WORKER
app (file output, file outLog) novosort (string novosortexe, string inputFile[], string tmpdir, int thr, string sortoptions[], int memFlag){
	// processing multi-input files together (merging files)
	novosortexe "--index" "-m" memFlag "--tmpdir" tmpdir "--threads" thr inputFile "-o" output sortoptions @stderr=outLog; 
	// novosort has dual function to also mark duplicates
}
@dispatch=WORKER
app (file outputfile, file outLog, file metricsfile) picard (string javaexe, string java_heap, string picardjar, string tmpdir, file inputFile ){
	javaexe java_heap "-jar" picardjar "MarkDuplicates" "INPUT=" inputFile "OUTPUT=" outputfile "METRICS_FILE=" metricsfile "TMP_DIR=" tmpdir "ASSUME_SORTED=true" "MAX_RECORDS_IN_RAM=null" "CREATE_INDEX=true" "VALIDATION_STRINGENCY=SILENT"@stderr=outLog;

}

@dispatch=WORKER
app (file output) samtools_flagstat(string samtoolsexe, file inputFile){
	samtoolsexe "flagstat" inputFile @stdout=output;
}

