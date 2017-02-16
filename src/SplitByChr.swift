/**********************************************************
 Summary
***********************************************************
Input:
	- Runfile
	- File with list of bam files to be split by chromosome/contig 
	    - Specified with --input flag on the command line
Output:
	- Dedupped and sorted bam files for each of the input samples
	- File with list of the output bams for a particular chromosome
	    - Numbers of files = number of chromosomes
	    - Files for each chromosome located within the same directory

*/

import files;
import string;
import sys;

import bioapps.align_dedup;
import generalfunctions.general;

/************************
 Parse command-line args
*************************/

argv_accept("runfile", "input"); // return error if other flags are given

string inFile = argv("input");
string configFilename = argv("runfile");

// Create input bam array
string inputBams[] = file_lines(input_file(inFile));

/**************
 Parse runfile
***************/

file configFile = input_file(configFilename);									   
string configFileData[] = file_lines(configFile);

// Gather variables
string vars[string] = getConfigVariables(configFileData);

/******************************
 Begin Splitting by Chromosome
*******************************/

indices = split(vars["CHRNAMES"], ":");

mkdir( strcat(vars["OUTPUTDIR"], "/piping_files/chr_split_files") ) =>

foreach chr in indices {
	
	// Initialized with a blank string so this file can be appended to later
	file chrOutList < strcat(vars["OUTPUTDIR"], "/piping_files/chr_split_files/", chr, "_bam_list.txt"
				) > = write("") =>

	foreach bam in inputBams {

		// Get everything in the file name except for the '.bam' extension
		string prefix = substring( bam, 0, strlen(bam) - 4 );

		/*
		 Create chromosome split dedupped and sorted bam file
		*/
		string fileName = basename_string(prefix);
		// Grab the part that is not '.wDedups.sorted'
		string sampleName = substring( fileName, 0, strlen(fileName) - 15);
		
		/********************										   
 		 Split by chromosome										    
 		*********************/
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
		if ( checkBam(vars, chrDedupSortedBam) ) {
			/*
			 SUCCESS
			*/
			append( chrOutList, strcat(filename(chrDedupSortedBam), "\n") );
		}
		else {
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

