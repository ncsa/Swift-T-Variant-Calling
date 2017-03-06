/**********************************************************
 Summary
***********************************************************
Input:
	- Runfile
Output:
	- g.vcf files for each sample after their chromosome separated vcf files were combined 
*/

import files;
import string;
import sys;

import generalfunctions.general;
import bioapps.merge_vcf;


/************************	  
 Parse command-line args	   
*************************/	 
	   
argv_accept("runfile"); // return error if other flags are given	   
	   
string configFilename = argv("runfile");

/**************	    
 Parse runfile	     
***************/

file configFile = input_file(configFilename);	      
string configFileData[] = file_lines(configFile);

// Gather variables
string vars[string] = getConfigVariables(configFileData);
	   
// Get input sample list	   
file sampleInfoFile = input_file(vars["SAMPLEINFORMATION"]);	       
string sampleLines[] = file_lines(sampleInfoFile);	 
	   
string indices[] = split(vars["CHRNAMES"], ":");

foreach sample in sampleLines {
	string sampleName = split(sample, " ")[0];

	/*
	Get the names of each of the chromosome split sample files, add a "--variants" string
	before them, and feed the array to CombineGVCF
	*/

	// This array holds all of the chromosome vcf files for this sample along with the --variants flags
	string chrSampleArray[];

	foreach chr, index in indices {
		/*
		For each sample.vcf added to the chrSampleArray needs "--variant" added in
		front of it. Note: simply concatenating "--variant " with the sample name
		then adding to the array produces an error in the command line call
		*/

		// Generate path to the g.vcf file for this sample and chromosome
		string location = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/", 
					 sampleName, ".wDedups.sorted.recalibrated.", chr, ".g.vcf"
					);

		int varFlagPos = index * 2;
		int namePos = varFlagPos + 1;

		chrSampleArray[varFlagPos] = "--variant";
		chrSampleArray[namePos] = location;
		
	}

	// CombineGVCF output file
	file gvcfSample < strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/",
				 sampleName, ".wDedups.sorted.recalibrated.g.vcf"
				) >;
	
	// Log file for CombineGVCFs
	file combineLog < strcat(vars["OUTPUTDIR"], "/", sampleName, "/logs/", sampleName, "_CombineGVCFs.log") >;

	gvcfSample, combineLog = CombineGVCFs(vars["JAVADIR"],
					      vars["GATKDIR"], 
					      strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),
					      strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),
					      chrSampleArray
					     );
}


