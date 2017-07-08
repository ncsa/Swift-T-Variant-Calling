// This script creates wraper functions around the merge_vcf  bioapps; with the objective of aiding logging

import bioapps.merge_vcf;

(file outputfile, file outLog, file tmptimeLog) CombineGVCFs_logged (string javaexe, string gatkjar, string reference, string dbsnp, string variants[] ) {
	string startmsg = strcat("CombineGVCFs start", "\t", toString(clock_seconds()), "\n") =>
	outputfile, outLog = CombineGVCFs (javaexe, gatkjar, reference, dbsnp, variants) =>
	string endmsg = strcat("CombineGVCFs end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}


