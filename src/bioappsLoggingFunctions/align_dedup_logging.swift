// This script creates wraper functions around the alignment and deduplication bioapps; with the objective of aiding logging 

import bioapps.align_dedup;

/*
 
//Ideally, this function would work, and save all the re-coding effort; but not at the moment :(

(boolean done) logme (string function_name, string function_call, file timingLog){
	string startmsg = strcat(function_name, " start", "\t", toString(clock_seconds()), "\n") =>
	function_call =>
	string endmsg = strcat(function_name, " end", "\t", toString(clock_seconds()), "\n") =>
	timeLog = echo(strcat(startmsg, endmsg)) =>
	done = true;
}
*/

/////// Alignment functions:

(file output, file outLog, file tmptimeLog) bwa_mem_logged (string bwaexe, string read1, string read2, string INDEX, string bwamemparams[], int PBSCORES,  string rgheader){
	string startmsg = strcat("BWAMEM start", "\t", toString(clock_seconds()), "\n") =>
	output, outLog = bwa_mem (bwaexe, read1, read2, INDEX, bwamemparams, PBSCORES, rgheader) =>
	string endmsg = strcat("BWAMEM end", "\t", toString(clock_seconds()), "\n") => 
	tmptimeLog = write(strcat(startmsg, endmsg));
}


(file output, file outLog, file tmptimeLog) bwa_mem_logged (string bwaexe, string read1, string INDEX, string bwamemparams[], int PBSCORES,  string rgheader){
	string startmsg = strcat("BWAMEM start", "\t", toString(clock_seconds()), "\n") =>
	output, outLog = bwa_mem (bwaexe, read1, INDEX, bwamemparams, PBSCORES, rgheader) =>
	string endmsg = strcat("BWAMEM end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file output, file outLog, file tmptimeLog) novoalign_logged (string novoalignexe, string read1, string read2, string INDEX, string novoalignparams[], int PBSCORES, string rgheader){
	string startmsg = strcat("NOVOALIGN start", "\t", toString(clock_seconds()), "\n");
	output, outLog =  novoalign (novoalignexe, read1, read2, INDEX, novoalignparams, PBSCORES, rgheader) => 
	string endmsg = strcat("NOVOALIGN end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg)); 
}

(file output, file outLog, file tmptimeLog) novoalign_logged (string novoalignexe, string read1, string INDEX, string novoalignparams[], int PBSCORES, string rgheader){
	string startmsg = strcat("NOVOALIGN start", "\t", toString(clock_seconds()), "\n") =>
	output, outLog =  novoalign (novoalignexe, read1, INDEX, novoalignparams, PBSCORES, rgheader) => 
	string endmsg = strcat("NOVOALIGN end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg)); 
}

(file output, file tmptimeLog) samtools_view_logged (string samtoolsexe, file inputFile, int thr, string args[]){
	string startmsg = strcat("SAMTOOLS_VIEW start", "\t", toString(clock_seconds()), "\n") => 
	output =  samtools_view(samtoolsexe, inputFile, thr, args) =>
	string endmsg =  strcat("SAMTOOLS_VIEW end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg)); 
}

(int numAlignments, file tmptimeLog) samtools_view2_logged(string samtoolsexe, string inputFile){
	string startmsg = strcat("SAMTOOLS_VIEW start", "\t", toString(clock_seconds()), "\n") =>
	numAlignments = samtools_view2(samtoolsexe, inputFile) =>
	string endmsg = strcat("SAMTOOLS_VIEW end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file output, file outLog, file tmptimeLog) samblaster_logged(string samblasterexe, file inputFile){
	string startmsg = strcat("SAMBLASTER start", "\t", toString(clock_seconds()), "\n") =>
	output, outLog = samblaster( samblasterexe, inputFile) =>
	string endmsg = strcat("SAMBLASTER end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file output, file outLog, file tmptimeLog) novosort_logged (string novosortexe, file inputFile, string tmpdir, int thr, string sortoptions[], int memFlag){
	string startmsg = strcat("NOVOSORT start", "\t", toString(clock_seconds()), "\n") =>
	output, outLog = novosort (novosortexe, inputFile, tmpdir, thr, sortoptions, memFlag) =>
	string endmsg = strcat("NOVOSORT end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file output, file outLog, file tmptimeLog) novosort_logged (string novosortexe, string inputFile[], string tmpdir, int thr, string sortoptions[], int memFlag){
	string startmsg = strcat("NOVOSORT start", "\t", toString(clock_seconds()), "\n") =>
	output, outLog= novosort (novosortexe, inputFile, tmpdir, thr, sortoptions, memFlag) =>
	string endmsg = strcat("NOVOSORT end", "\t", toString(clock_seconds()), "\n") => 
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file outputfile, file outLog, file metricsfile, file tmptimeLog) picard_logged (string javaexe, string picardjar, string tmpdir, file inputFile){
	string startmsg = strcat("PICARD start", "\t", toString(clock_seconds()), "\n") =>
	outputfile, outLog, metricsfile = picard (javaexe, picardjar, tmpdir, inputFile) =>
	string endmsg = strcat("PICARD end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file output, file tmptimeLog) samtools_flagstat_logged(string samtoolsexe, file inputFile){
	string startmsg = strcat("SAMTOOLS_FLAGSTATS start", "\t", toString(clock_seconds()), "\n") =>
	output = samtools_flagstat(samtoolsexe, inputFile) =>
	string endmsg = strcat("SAMTOOLS_FLAGSTATS end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}



