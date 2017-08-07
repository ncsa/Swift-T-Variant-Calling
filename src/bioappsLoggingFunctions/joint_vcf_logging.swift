// This script creates wraper functions around the joint_vcf bioapps; with the objective of aiding logging 

import bioapps.joint_vcf;

(file outputfile, file outLog, file tmptimeLog) GenotypeGVCFs_logged (string javaexe, string gatkjar, string reference, string variants[], string threads) {
	string startmsg = strcat("ALL", "\t", "ALL", "\t", "GenotypeGVCFs start", "\t", toString(clock_seconds()), "\n") =>
	outputfile, outLog = GenotypeGVCFs (javaexe, gatkjar, reference, variants, threads) =>
	string endmsg = strcat("ALL", "\t", "ALL", "\t", "GenotypeGVCFs end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}

