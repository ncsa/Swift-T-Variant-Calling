// This script creates wraper functions around the alignment and deduplication bioapps; with the objective of aiding logging 

import bioapps.align_dedup;

//////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////

(file output, file outLog, file tmptimeLog) bwa_mem_logged (string bwaexe, string read1, string read2, string INDEX, string bwamemparams[], int PBSCORES,  string rgheader, string sampleName){
	string startmsg = strcat(sampleName, "\t", "ALL", "\t", "BWAMEM start", "\t", toString(clock_seconds()), "Align\n") =>
	output, outLog = bwa_mem (bwaexe, read1, read2, INDEX, bwamemparams, PBSCORES, rgheader) =>
	string endmsg = strcat(sampleName, "\t", "ALL", "\t", "BWAMEM end", "\t", toString(clock_seconds()), "Align\n") => 
	tmptimeLog = write(strcat(startmsg, endmsg));
}


(file output, file outLog, file tmptimeLog) bwa_mem_logged (string bwaexe, string read1, string INDEX, string bwamemparams[], int PBSCORES,  string rgheader, string sampleName){
	string startmsg = strcat(sampleName, "\t", "ALL", "\t", "BWAMEM start", "\t", toString(clock_seconds()), "Align\n") =>
	output, outLog = bwa_mem (bwaexe, read1, INDEX, bwamemparams, PBSCORES, rgheader) =>
	string endmsg = strcat(sampleName, "\t", "ALL", "\t", "BWAMEM end", "\t", toString(clock_seconds()), "Align\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file output, file outLog, file tmptimeLog) novoalign_logged (string novoalignexe, string read1, string read2, string INDEX, string novoalignparams[], int PBSCORES, string rgheader, string sampleName){
	string startmsg = strcat(sampleName, "\t", "ALL", "\t", "NOVOALIGN start", "\t", toString(clock_seconds()), "Align\n");
	output, outLog =  novoalign (novoalignexe, read1, read2, INDEX, novoalignparams, PBSCORES, rgheader) => 
	string endmsg = strcat(sampleName, "\t", "ALL", "\t", "NOVOALIGN end", "\t", toString(clock_seconds()), "Align\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg)); 
}

(file output, file outLog, file tmptimeLog) novoalign_logged (string novoalignexe, string read1, string INDEX, string novoalignparams[], int PBSCORES, string rgheader, string sampleName){
	string startmsg = strcat(sampleName, "\t", "ALL", "\t", "NOVOALIGN start", "\t", toString(clock_seconds()), "Align\n") =>
	output, outLog =  novoalign (novoalignexe, read1, INDEX, novoalignparams, PBSCORES, rgheader) => 
	string endmsg = strcat(sampleName, "\t", "ALL", "\t", "NOVOALIGN end", "\t", toString(clock_seconds()), "Align\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg)); 
}

(file output, file tmptimeLog) samtools_view_logged (string samtoolsexe, file inputFile, int thr, string args[], string sampleName){
	string startmsg = strcat(sampleName, "\t", "ALL", "\t", "SAMTOOLS_VIEW start", "\t", toString(clock_seconds()), "Align\n") => 
	output =  samtools_view(samtoolsexe, inputFile, thr, args) =>
	string endmsg =  strcat(sampleName, "\t", "ALL", "\t", "SAMTOOLS_VIEW end", "\t", toString(clock_seconds()), "Align\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg)); 
}

(file output, file outLog, file tmptimeLog) samblaster_logged(string samblasterexe, file inputFile, string sampleName){
	string startmsg = strcat(sampleName, "\t", "ALL", "\t", "SAMBLASTER start", "\t", toString(clock_seconds()), "DedupSort\n") =>
	output, outLog = samblaster( samblasterexe, inputFile) =>
	string endmsg = strcat(sampleName, "\t", "ALL", "\t", "SAMBLASTER end", "\t", toString(clock_seconds()), "DedupSort\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file output, file outLog, file tmptimeLog) novosort_logged (string novosortexe, file inputFile, string tmpdir, int thr, string sortoptions[], int memFlag, string sampleName){
	string startmsg = strcat(sampleName, "\t", "ALL", "\t", "NOVOSORT start", "\t", toString(clock_seconds()), "DedupSort\n") =>
	output, outLog = novosort (novosortexe, inputFile, tmpdir, thr, sortoptions, memFlag) =>
	string endmsg = strcat(sampleName, "\t", "ALL", "\t", "NOVOSORT end", "\t", toString(clock_seconds()), "DedupSort\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file output, file outLog, file tmptimeLog) novosort_logged (string novosortexe, string inputFile[], string tmpdir, int thr, string sortoptions[], int memFlag){
	string startmsg = strcat("NOVOSORT start", "\t", toString(clock_seconds()), "DedupSort\n") =>
	output, outLog= novosort (novosortexe, inputFile, tmpdir, thr, sortoptions, memFlag) =>
	string endmsg = strcat("NOVOSORT end", "\t", toString(clock_seconds()), "DedupSort\n") => 
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file outputfile, file outLog, file metricsfile, file tmptimeLog) picard_logged (string javaexe, string java_heap, string picardjar, string tmpdir, file inputFile, string sampleName){
	string startmsg = strcat(sampleName, "\t", "ALL", "\t", "PICARD start", "\t", toString(clock_seconds()), "DedupSort\n") =>
	outputfile, outLog, metricsfile = picard (javaexe, strcat("-Xmx", java_heap), picardjar, tmpdir, inputFile) =>
	string endmsg = strcat(sampleName, "\t", "ALL", "\t", "PICARD end", "\t", toString(clock_seconds()), "DedupSort\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}

/* flagstats command doesn't take much time, so it is not really that meaningful to log it
(file output, file tmptimeLog) samtools_flagstat_logged(string samtoolsexe, file inputFile, string sampleName){
	string startmsg = strcat(sampleName, "\t", "ALL", "\t", "SAMTOOLS_FLAGSTATS start", "\t", toString(clock_seconds()), "DedupSort\n") =>
	output = samtools_flagstat(samtoolsexe, inputFile) =>
	string endmsg = strcat(sampleName, "\t", "ALL", "\t", "SAMTOOLS_FLAGSTATS end", "\t", toString(clock_seconds()), "DedupSort\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}
*/
