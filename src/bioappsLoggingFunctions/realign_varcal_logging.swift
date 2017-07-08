// This script creates wraper functions around the realignment and variant calling bioapps; with the objective of aiding logging 

import bioapps.realign_varcal;

/* indexing doesn't take much time, so it is not really that meaningful to log it
(void signal, file tmptimeLog) samtools_index_logged (string samtoolsdir, file inputFilename) {
	string startmsg = strcat("SAMTOOLS INDEX start", "\t", toString(clock_seconds()), "\n") =>
	signal = samtools_index (samtoolsdir, inputFilename) =>
        string endmsg = strcat("SAMTOOLS INDEX end", "\t", toString(clock_seconds()), "\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}
*/


(file outputfile, file outLog, file tmptimeLog) RealignerTargetCreator_logged (string javaexe, string gatkjar, string reference,
   file inputFile, int thr, string known[]
  ) {
	string startmsg = strcat("RealignerTargetCreator start", "\t", toString(clock_seconds()), "\n") =>
	outputfile, outLog = RealignerTargetCreator (javaexe, gatkjar, reference, inputFile, thr, known) =>
	string endmsg = strcat("RealignerTargetCreator end", "\t", toString(clock_seconds()), "\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file outputfile, file outLog, file tmptimeLog) IndelRealigner_logged (string javaexe, string gatkjar, string reference,
				   file inputFile, string known[], file intervals
				  ) {
	string startmsg = strcat("IndelRealigner start", "\t", toString(clock_seconds()), "\n") =>
	outputfile, outLog = IndelRealigner (javaexe, gatkjar, reference, inputFile, known, intervals) =>
	string endmsg = strcat("IndelRealigner end", "\t", toString(clock_seconds()), "\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file outputfile, file outLog, file tmptimeLog) BaseRecalibrator_logged (string javaexe, string gatkjar, string reference,
				     file inputFile, int thr, string knownindels[], string dbsnp
				    ) {
	string startmsg = strcat("BaseRecalibrator start", "\t", toString(clock_seconds()), "\n") =>
	outputfile, outLog = BaseRecalibrator (javaexe, gatkjar, reference, inputFile, thr, knownindels, dbsnp) =>
	string endmsg = strcat("BaseRecalibrator end", "\t", toString(clock_seconds()), "\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}

(file outputfile, file outLog, file tmptimeLog) PrintReads_logged (string javaexe, string gatkjar, string reference,
			       file inputFile, int thr, file recalreport
			      ) {
	string startmsg = strcat("PrintReads start", "\t", toString(clock_seconds()), "\n") =>
	outputfile, outLog = PrintReads (javaexe, gatkjar, reference, inputFile, thr, recalreport) =>
	string endmsg = strcat("PrintReads end", "\t", toString(clock_seconds()), "\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}


// The chromosome splitting version
(file outputfile, file outLog, file tmptimeLog) HaplotypeCaller_logged (string javaexe, string gatkjar, string reference,
						   file inputFile, string dbsnp, int thr, int ploidy, string chr
						  ) {
	string startmsg = strcat("HaplotypeCaller start", "\t", toString(clock_seconds()), "\n") =>
	outputfile, outLog = HaplotypeCaller (javaexe, gatkjar, reference, inputFile, dbsnp, thr, ploidy, chr) =>
	string endmsg = strcat("HaplotypeCaller end", "\t", toString(clock_seconds()), "\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}


// The whole genome version
(file outputfile, file outLog, file tmptimeLog) HaplotypeCaller_logged (string javaexe, string gatkjar, string reference,      
						   file inputFile, string dbsnp, int thr	  
						  ) {
	string startmsg = strcat("HaplotypeCaller start", "\t", toString(clock_seconds()), "\n") =>
	outputfile, outLog = HaplotypeCaller (javaexe, gatkjar, reference, inputFile, dbsnp, thr) =>
	string endmsg = strcat("HaplotypeCaller end", "\t", toString(clock_seconds()), "\n") =>  
	tmptimeLog = write(strcat(startmsg, endmsg));
}
