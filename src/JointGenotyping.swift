/**********************************************************
 Summary
***********************************************************
Input:
	- Runfile

Output:
	- Joint Genotype vcf file
*/

import files;
import string;
import sys;

import bioapps.joint_vcf;
import generalfunctions.general;

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

/******************
 Joint Genotyping
*******************/

// The joint genotype output file
file jointVCF < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/jointVCFs/jointVCFcalled.vcf") >;

// Log file for Joint Genotyping
file jointLog < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/jointVCFs/jointVCF.log") >;

mkdir(strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/jointVCFs"));

// This array holds the vcf file for each sample along with the "--variants" flags necessary for GenotypeGVCFs
string variantSampleArray[];

foreach sample, index in sampleLines {

	string sampleName = split(sample, " ")[0];

	// Generate path to the g.vcf file for this sample
	string location = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/",
				 sampleName, ".wDedups.sorted.recalibrated.g.vcf"
				);

	/*
	For each sample.vcf added to the variantSampleArray needs "--variant" in front of it  
	Note: simply concatenating "--variant " with the sample name then adding to the array   
	produces an error in the command line call
	*/											 
	int varFlagPos = index * 2; 
	int namePos = varFlagPos + 1;
	variantSampleArray[varFlagPos] = "--variant";
	variantSampleArray[namePos] = location;						  
}

jointVCF, jointLog = GenotypeGVCFs(vars["JAVADIR"], vars["GATKDIR"], strcat(vars["REFGENOMEDIR"],
				   "/", vars["REFGENOME"]), variantSampleArray
				  );

