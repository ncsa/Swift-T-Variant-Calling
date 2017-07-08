/*

*****************************
 Pseudocode of Main Function
*****************************

() jointGenotypingMain(file inputVCFs[]) {
	if (JOINT_GENOTYPING_STAGE variable == "Y") {
		**************************
		*** EXECUTE THIS STAGE ***
		**************************

		- Create output and log file handles

		foreach sampleVCF and sampleVCFs {
			- Perform Joint Genotyping
		}
	}
}
*/

import files;
import string;
import sys;
import io;

import bioapps.joint_vcf;
import generalfunctions.general;
import bioappsLoggingFunctions.joint_vcf_logging;

/******************
 Main function for Joint Genotyping
*******************/

jointGenotypingMain(file inputVCFs[], string vars[string]) {
	// Since this is the last step, I only check to make sure this step is one of the executed stages.
	// If it is not, then nothing happens.
	if (vars["JOINT_GENOTYPING_STAGE"] == "Y") {

		// The joint genotype output file
		file jointVCF < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/jointVCFs/jointVCFcalled.vcf") >;
	
		// Log file for Joint Genotyping
		file jointLog < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/jointVCFs/jointVCF.log") >;
		mkdir(strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/jointVCFs"));

		string tmpLogDir = strcat(vars["TMPDIR"], "/timinglogs/" );
		file tmpjointLog < strcat(tmpLogDir, "jointVCF.log") >;	
				
		// This array holds the vcf file for each sample along with the "--variants" flags necessary for GenotypeGVCFs
		string variantSampleArray[];
		
		foreach sampleVCF, index in inputVCFs {
			/*
			For each sample.vcf added to the variantSampleArray needs "--variant" in front of it  
			Note: simply concatenating "--variant " with the sample name then adding to the array   
			produces an error in the command line call
			*/											 
			int varFlagPos = index * 2; 
			int namePos = varFlagPos + 1;
			variantSampleArray[varFlagPos] = "--variant";
			variantSampleArray[namePos] = filename(sampleVCF);						  
		}
		
		jointVCF, jointLog, tmpjointLog = GenotypeGVCFs_logged (vars["JAVAEXE"], vars["GATKJAR"], strcat(vars["REFGENOMEDIR"],
						   "/", vars["REFGENOME"]), variantSampleArray, vars["CORES"]
						  );
	}
}
