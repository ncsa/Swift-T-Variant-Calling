// This script creates wraper functions around the realignment and variant calling bioapps; with the objective of aiding logging 

import bioapps.realign_varcal;

////////////////////////////////////////////////////////////////
/* indexing doesn't take much time, so it is not really that meaningful to log it
(void signal, file tmptimeLog) samtools_index_logged (string samtoolsdir, file inputFilename) {
	string startmsg = strcat("SAMTOOLS INDEX start", "\t", toString(clock_seconds()), "Realign\n") =>
	signal = samtools_index (samtoolsdir, inputFilename) =>
        string endmsg = strcat("SAMTOOLS INDEX end", "\t", toString(clock_seconds()), "Realign\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}
*/
///////////////////////////////////////////////////////////////

(file outputfile, file outLog, file tmptimeLog) RealignerTargetCreator_logged (string javaexe, string java_heap, string gatkjar, string reference,
   file inputFile, int thr, string known[], string sampleName, string chr
  ) {
	string startmsg = strcat(sampleName, "\t", chr, "\t", "RealignerTargetCreator start", "\t", toString(clock_seconds()), "VCRun\n") =>
	outputfile, outLog = RealignerTargetCreator (javaexe, strcat("-Xmx", java_heap), gatkjar, reference, inputFile, thr, known) =>
	string endmsg = strcat(sampleName, "\t", chr, "\t", "RealignerTargetCreator end", "\t", toString(clock_seconds()), "VCRun\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file outputfile, file outLog, file tmptimeLog) IndelRealigner_logged (string javaexe, string java_heap, string gatkjar, string reference,
				   file inputFile, string known[], file intervals, string sampleName, string chr
				  ) {
	string startmsg = strcat(sampleName, "\t", chr, "\t", "IndelRealigner start", "\t", toString(clock_seconds()), "VCRun\n") =>
	outputfile, outLog = IndelRealigner (javaexe, strcat("-Xmx", java_heap), gatkjar, reference, inputFile, known, intervals) =>
	string endmsg = strcat(sampleName, "\t", chr, "\t", "IndelRealigner end", "\t", toString(clock_seconds()), "VCRun\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file outputfile, file outLog, file tmptimeLog) BaseRecalibrator_logged (string javaexe, string java_heap, string gatkjar, string reference,
				     file inputFile, int thr, string knownindels[], string dbsnp, string sampleName, string chr
				    ) {
	string startmsg = strcat(sampleName, "\t", chr, "\t", "BaseRecalibrator start", "\t", toString(clock_seconds()), "VCRun\n") =>
	outputfile, outLog = BaseRecalibrator (javaexe, strcat("-Xmx", java_heap), gatkjar, reference, inputFile, thr, knownindels, dbsnp) =>
	string endmsg = strcat(sampleName, "\t", chr, "\t", "BaseRecalibrator end", "\t", toString(clock_seconds()), "VCRun\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file outputfile, file outLog, file tmptimeLog) PrintReads_logged (string javaexe, string java_heap, string gatkjar, string reference,
			       file inputFile, int thr, file recalreport, string sampleName, string chr
			      ) {
	string startmsg = strcat(sampleName, "\t", chr, "\t", "PrintReads start", "\t", toString(clock_seconds()), "VCRun\n") =>
	outputfile, outLog = PrintReads (javaexe, strcat("-Xmx", java_heap), gatkjar, reference, inputFile, thr, recalreport) =>
	string endmsg = strcat(sampleName, "\t", chr, "\t", "PrintReads end", "\t", toString(clock_seconds()), "VCRun\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}


// The chromosome splitting version
(file outputfile, file outLog, file tmptimeLog) HaplotypeCaller_logged (string javaexe, string java_heap, string gatkjar, string reference,
						   file inputFile, string dbsnp, int thr, int ploidy, string chr, string sampleName
						  ) {
	string startmsg = strcat(sampleName, "\t", chr, "\t", "HaplotypeCaller start", "\t", toString(clock_seconds()), "VCRun\n") =>
	outputfile, outLog = HaplotypeCaller (javaexe, strcat("-Xmx", java_heap), gatkjar, reference, inputFile, dbsnp, thr, ploidy, chr) =>
	string endmsg = strcat(sampleName, "\t", chr, "\t", "HaplotypeCaller end", "\t", toString(clock_seconds()), "VCRun\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}


// The whole genome version
(file outputfile, file outLog, file tmptimeLog) HaplotypeCaller_logged (string javaexe, string java_heap, string gatkjar, string reference,      
						   file inputFile, string dbsnp, int thr, string sampleName 
						  ) {
	string startmsg = strcat(sampleName, "\t", "ALL", "\t", "HaplotypeCaller start", "\t", toString(clock_seconds()), "VCRun\n") =>
	outputfile, outLog = HaplotypeCaller (javaexe, strcat("-Xmx", java_heap), gatkjar, reference, inputFile, dbsnp, thr) =>
	string endmsg = strcat(sampleName, "\t", "ALL", "\t", "HaplotypeCaller end", "\t", toString(clock_seconds()), "VCRun\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}
