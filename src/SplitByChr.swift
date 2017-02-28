/**********************************************************
 Summary
***********************************************************
Input:
	- Runfile

Output:
	- Dedupped and sorted bam files for each of the input samples
*/

import files;
import string;
import sys;

import bioapps.align_dedup;
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

indices = split(vars["CHRNAMES"], ":");

/**********************
 Find the aligned bams
***********************/

// Create input bam array
string inputBams[];

// Fill in the array
foreach sample, index in sampleLines {
	string sampleName = split(sample, " ")[0];
	
	// Generate the path to the aligned file for this sample
	string location = strcat(vars["OUTPUTDIR"], "/", sampleName, "/align/", sampleName, ".wDedups.sorted.bam");
 
	// Add the location to the array
	inputBams[index] = location;
}

/******************************
 Begin Splitting by Chromosome
*******************************/

foreach chr in indices {

	foreach bam in inputBams {

		// Get everything in the file name except for the '.bam' extension
		string prefix = substring( bam, 0, strlen(bam) - 4 );

		string fileName = basename_string(prefix);
		// Grab the part that is not '.wDedups.sorted'
		string sampleName = substring( fileName, 0, strlen(fileName) - 15);
		
		/********************
 		 Split by chromosome
 		*********************/
		
		/*
		  Create chromosome split dedupped and sorted bam file
		*/
		file chrDedupSortedBam < strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/",
						fileName, ".", chr, ".bam"
					       ) >;

		// Subtract 1 because the main thread takes up is a thread, "threads" defines the number of additional
		//   threads
		int threads = ( string2int(vars["PBSCORES"]) %/ string2int(vars["PROCPERNODE"]) ) - 1;

		chrDedupSortedBam = samtools_view(vars["SAMTOOLSDIR"], input(bam),
						  threads, [strcat(chr)]
						 );
		
		// Check whether the splitting was successful
		if ( ! checkBam(vars, chrDedupSortedBam) ) {
			/*
			 FAILURE
			*/
			string message = strcat("FAILURE: ", filename(chrDedupSortedBam), " contains no alignments. ",
						"Chromosome splitting failed\n"
					       );
			// Update the failure log file
			append(input(strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/docs/Failures.log")),
				     message
			      );
		}
	}
}

