import files;
import string;
import sys;

import bioapps.align_dedup;
import generalfunctions.general;

/******************************
 Main function
*******************************/

(file outputBams[][]) splitByChrMain(file inputBams[], string vars[string], file failLog) {
	string indices[] = split(vars["CHRNAMES"], ":");

	foreach chr, chrIndex in indices {
		foreach bam, sampleIndex in inputBams {

			// Get everything in the file name except for the '.bam' extension
			string prefix = substring( filename(bam), 0, strlen(filename(bam)) - 4 );
	
			string fileName = basename_string(prefix);
			// Grab the part that is not '.wDedups.sorted'
			string sampleName = substring( fileName, 0, strlen(fileName) - 15);

			string outputName = strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/", fileName, ".",
						   chr, ".bam"
						  );			

			if (vars["CHR_SPLIT_STAGE"] == "Y") {

				/********************
	 			 Split by chromosome
	 			*********************/
			
				/*
				  Create chromosome split dedupped and sorted bam file
				*/
				file chrDedupSortedBam < outputName >;
	
				// Subtract 1 because the main thread takes up is a thread, "threads"
				//   defines the number of additional threads
				int threads = ( string2int(vars["PBSCORES"]) %/ string2int(vars["PROCPERNODE"]) ) - 1;
	
				chrDedupSortedBam = samtools_view(vars["SAMTOOLSDIR"], bam,
								  threads, [strcat(chr)]
								 );
			
				// Check whether the splitting was successful
				if ( ! checkBam(vars, chrDedupSortedBam) ) {
					/*
					 FAILURE
					*/
					string message = strcat("FAILURE: ", filename(chrDedupSortedBam), 
								" contains no alignments. Chromosome splitting failed\n"
							       );
					// Update the failure log file
					append(input(strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"],
						            "/docs/Failures.log")
							   ),
						     message
					      );
				}
				else {
					/*
					  SUCCESS
					*/
					outputBams[chrIndex][sampleIndex] = chrDedupSortedBam;
				}
			}
			// If this stage is to be skipped
			else {
				if (file_exists(outputName)) {
					outputBams[chrIndex][sampleIndex] = input(outputName);
				} 
				else {
					append(failLog, strcat("ERROR: ", outputName, 
							       " not found. Did you set CHR_SPLIT_STAGE to",   
							       " 'N' by accident?\n"					    
							      )
				       );
				}
			}
		}
	}
}	
