// This script creates wraper functions around the joint_vcf bioapps; with the objective of aiding logging 

import bioapps.joint_vcf;

(file outputfile, file outLog, file tmptimeLog) GenotypeGVCFs_logged (string javaexe, string java_heap, string gatkjar, string reference, string variants[], string threads) {
	string startmsg = strcat("ALL", "\t", "ALL", "\t", "GenotypeGVCFs start", "\t", toString(clock_seconds()), "\tjointGenotypeRun\n") =>
	outputfile, outLog = GenotypeGVCFs (javaexe, strcat("-Xmx", java_heap),  gatkjar, reference, variants, threads) =>
	string endmsg = strcat("ALL", "\t", "ALL", "\t", "GenotypeGVCFs end", "\t", toString(clock_seconds()), "\tjointGenotypeRun\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}


