// This script creates wraper functions around the joint_vcf bioapps; with the objective of aiding logging 

import bioapps.joint_vcf;

(file outputfile, file outLog, file tmptimeLog) GenotypeGVCFs_logged (string javaexe, string gatkjar, string reference, string variants[], string threads) {
	string startmsg = strcat("GenotypeGVCFs start", "\t", toString(clock_seconds()), "\n") =>
	outputfile, outLog = GenotypeGVCFs (javaexe, gatkjar, reference, variants, threads) =>
	string endmsg = strcat("GenotypeGVCFs end", "\t", toString(clock_seconds()), "\n") =>
	tmptimeLog = write(strcat(startmsg, endmsg));
}


