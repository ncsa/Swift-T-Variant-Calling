/*

*****************************
 Pseudocode of Run Function
*****************************

() jointGenotypingRun(file inputVCFs[]) {
	if (JOINT_GENOTYPING_STAGE variable == "Y") {
		**************************
		*** EXECUTE THIS STAGE ***
		**************************
		Check the GATKEXE is properly specified in the runfile

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
 Run function for Joint Genotyping
*******************/

jointGenotypingRun(file inputVCFs[], string vars[string], file timeLog) {
	// Since this is the last step, I only check to make sure this step is one of the executed stages.
	// If it is not, then nothing happens.
	if (vars["JOINT_GENOTYPING_STAGE"] == "Y" ||
	    vars["JOINT_GENOTYPING_STAGE"] == "Yes" ||
	    vars["JOINT_GENOTYPING_STAGE"] == "YES" ||
	    vars["JOINT_GENOTYPING_STAGE"] == "y" ||
	    vars["JOINT_GENOTYPING_STAGE"] == "yes" ||
	    vars["JOINT_GENOTYPING_STAGE"] == "End" ||
	    vars["JOINT_GENOTYPING_STAGE"] == "end" ||
	    vars["JOINT_GENOTYPING_STAGE"] == "E" ||
	    vars["JOINT_GENOTYPING_STAGE"] == "e"
	   ) {
		exec_check(vars["GATKJAR"], "GATKJAR");
		exec_check(vars["JAVAEXE"], "JAVAEXE");

		// The joint genotype output file
		file jointVCF < strcat(vars["OUTPUTDIR"], "/deliverables/jointVCFs/jointVCFcalled.vcf") >;
	
		// Log file for Joint Genotyping
		file jointLog < strcat(vars["OUTPUTDIR"], "/deliverables/jointVCFs/jointVCF.log") >;
		mkdir(strcat(vars["OUTPUTDIR"], "/deliverables/jointVCFs"));

		string tmpLogDir = strcat(vars["TMPDIR"], "/timinglogs/jointGenologs/"); // Path is stage specific
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
		
		jointVCF, jointLog, tmpjointLog = GenotypeGVCFs_logged (vars["JAVAEXE"], vars["JAVA_MAX_HEAP_SIZE"], vars["GATKJAR"], vars["REFGENOME"], variantSampleArray, vars["CORES_PER_NODE"]
					  ) =>
		// This should be the only location where the 'logging' function should be called
		//   outside of VariantCalling.swift
		logging(variables["TMPDIR"], timeLog, "jointGenologs");
	}
}
