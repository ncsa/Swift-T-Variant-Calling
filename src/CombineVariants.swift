/*

*****************************
 Pseudocode of Main Function
*****************************

(file vcfOutfiles[]) combineVariantsMain(file inputVCFs) {
	
	foreach sample in samples {
	
		- Get the output name

		if (COMBINE_VARIANT_STAGE variable == "Y") {
			**************************
			*** EXECUTE THIS STAGE ***
			**************************

			- CombineGVCFs
			- Add the VCF file to the output array

		} else {
			**********************************
			*** THIS STAGE IS NOT EXECUTED ***
			**********************************

			if (this samples output file can be found) {
				*** SUCCESS ***
				- Add the VCF file to the output array
			} else {
				*** FAILURE ***
				- Write an error message to the fail log
			}
		}
	}
}
*/

import files;
import string;
import sys;

import generalfunctions.general;
import bioapps.merge_vcf;

(file vcfOutfiles[]) combineVariantsMain(file inputVCFs[][], string vars[string], file failLog) {

	foreach sampleSet, sampleIndex in inputVCFs {

		// Get the sample name by looking at the first item in the samples chromosome set
		// It has the form: sampleName.wDedups.sorted.recalibrated.chrA.g.vcf
		string baseName = basename(sampleSet[0]);
		string pieces[] = split(baseName, ".");
		string sampleName = pieces[0]; // Note: this assumes that the sample name has no '.' character in it

		string outputName = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/", sampleName,
					   ".wDedups.sorted.recalibrated.g.vcf"
					  );
	
		if (vars["COMBINE_VARIANT_STAGE"] == "Y") {

			// This array holds all of the chromosome vcf files for this sample along with the --variants flags
			string chrSampleArray[];

			foreach chrBam, index in sampleSet {
				/*												      
 				 For each sample.vcf added to the chrSampleArray needs "--variant" added in			      
 				 front of it. Note: simply concatenating "--variant " with the sample name			       
 				 then adding to the array produces an error in the command line call				     
 				*/			
				int varFlagPos = index * 2;
				int namePos = varFlagPos + 1; 
	
				chrSampleArray[varFlagPos] = "--variant";
				chrSampleArray[namePos] = filename(chrBam);

			}
			// CombineGVCF output file
			file gvcfSample < outputName >;

			// Log file for CombineGVCFs
			file combineLog < strcat(vars["OUTPUTDIR"], "/", sampleName,
						 "/logs/", sampleName, "_CombineGVCFs.log"
						) >;	 
	
			gvcfSample, combineLog = CombineGVCFs(vars["JAVAEXE"],							     
							      vars["GATKDIR"],							     
							      strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),			
							      strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),			    
							      chrSampleArray							       
							     );
			vcfOutfiles[sampleIndex] = gvcfSample;
		}
		// If this section is to be skipped
		else {
			if (file_exists(outputName)) {
				vcfOutfiles[sampleIndex] = input(outputName);
			}
			else {
				string message = strcat("ERROR: ", outputName, " not found. Did you set ",
							"COMBINE_VARIANT_STAGE to 'N' by accident?\n"
						       );
				append(failLog, message) =>
				exitIfFlagGiven(vars, message);
			}
		}
	}
}
