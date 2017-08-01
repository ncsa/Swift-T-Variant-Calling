/*

*****************************
 Pseudocode of Run Function
*****************************

(file vcfOutfiles[]) combineVariantsRun(file inputVCFs) {
	
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
import bioappsLoggingFunctions.merge_vcf_logging;

(file vcfOutfiles[]) combineVariantsRun(file inputVCFs[][], string vars[string], file failLog ) {

	foreach sampleSet, sampleIndex in inputVCFs {

		// Get the sample name by looking at the first item in the samples chromosome set
		// It has the form: sampleName.wDedups.sorted.recalibrated.chrA.g.vcf
		string baseName = basename(sampleSet[0]);
		string trimmed = substring(baseName, 0, strlen(baseName) - 10);  // gets rid of '.g.vcf' extension
		string pieces[] = split(trimmed, ".");
		string chr = pieces[size(pieces) - 1];          // Grabs the last part, which is the chromosome part

		//string sampleName = pieces[0]; // Note: this assumes that the sample name has no '.' character in it- 
		//This is a strange assumption given that it was never asserted in previous scripts that are able to deal with such scenarios

                // Removes '.wDedups.sorted.recalibrated.chr' part of the sample's name
                string sampleName = substring(trimmed, 0, strlen(trimmed) - strlen(chr) - 29); //verified 

		string outputName = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/", sampleName,
					   ".wDedups.sorted.recalibrated.g.vcf"
					  );
	
		if (vars["COMBINE_VARIANT_STAGE"] == "Y" ||
		    vars["COMBINE_VARIANT_STAGE"] == "Yes" ||
		    vars["COMBINE_VARIANT_STAGE"] == "YES" ||
		    vars["COMBINE_VARIANT_STAGE"] == "y" ||
		    vars["COMBINE_VARIANT_STAGE"] == "yes" ||
		    vars["COMBINE_VARIANT_STAGE"] == "End" ||
		    vars["COMBINE_VARIANT_STAGE"] == "end" ||
		    vars["COMBINE_VARIANT_STAGE"] == "E" ||
		    vars["COMBINE_VARIANT_STAGE"] == "e"
		   ) {
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

			string tmpLogDir = strcat(vars["TMPDIR"], "/timinglogs/" );
			file tmpcombineLog < strcat(tmpLogDir, sampleName, ".", chr, "_CombineGVCFs.log") > ; 
	
			gvcfSample, combineLog, tmpcombineLog = CombineGVCFs_logged (vars["JAVAEXE"], 
							      vars["GATKJAR"],							     
							      strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),			
							      strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),			    
							      chrSampleArray, sampleName 
							     ) =>
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
