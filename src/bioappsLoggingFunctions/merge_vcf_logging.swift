// This script creates wraper functions around the merge_vcf  bioapps; with the objective of aiding logging

import bioapps.merge_vcf;

(file outputfile, file outLog, file tmptimeLog) CombineGVCFs_logged (string javaexe, string java_heap, string gatkjar, string reference, string dbsnp, string variants[], string sampleName ) {
	string startmsg = strcat(sampleName, "\t", "ALL", "\t", "CombineGVCFs start", "\t", toString(clock_seconds()), "\n") =>
	outputfile, outLog = CombineGVCFs (javaexe, strcat("-Xmx", java_heap), gatkjar, reference, dbsnp, variants) =>
	string endmsg = strcat(sampleName, "\t", "ALL", "\t", "CombineGVCFs end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}


