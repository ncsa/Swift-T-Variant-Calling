/*
*****************************
 Pseudocode of Run Function
*****************************
(file outBamArray[]) alignRun() {
	foreach sample in samples {
		- Parse sample specific information and construct RG header
		- Create the sample output directories
		if (ALIGN_DEDUP_STAGE variable == "Y") {
			**************************
			*** EXECUTE THIS STAGE ***
			**************************
			- Create output file handles
			- Perform alignment
			if (alignment was Successful) {
				- add it to the output array
			} else {
				*** FAILURE ***
				- Write an error message to the fail log
			}
		} else {
			**********************************
			*** THIS STAGE IS NOT EXECUTED ***
			**********************************
			if (this samples output file can be found) {
				*** SUCCESS ***
				- Add the raw aligned bam file to the output array
			} else {
				*** FAILURE ***
				- Write an error message to the fail log
			}
		}
	}
	- return the outBamArray once all samples are processed
}
*/

import unix;					       
import files;					      
import string;					     
import sys;						
import io;

import bioapps.align_dedup;
import bioappsLoggingFunctions.align_dedup_logging;

import generalfunctions.general;

/********************************************************
 Helper Functions
********************************************************/

/*********
 Alignment
**********/
(file outputSam) alignReads(string vars[string], string sampleName, string reads[], string rgheader) {
	/*
	This function returns a .sam file because samblaster requires it
	To minimize memory usage, delete the .sam file after a .bam file is made from it
	*/

	int threads = string2int(vars["CORES_PER_NODE"]) %/ string2int(vars["PROGRAMS_PER_NODE"]);

	// Log file
	string LogDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/logs/");
	file alignedLog < strcat(LogDir, sampleName, "_Alignment.log") >;
	string tmpLogDir = strcat(vars["TMPDIR"], "/timinglogs/" );
 	file tmpalignedLog < strcat(tmpLogDir, sampleName, "_Alignment.log")>;
	
	if (vars["PAIRED"] == "1" ||
	    vars["PAIRED"] == "YES" ||
	    vars["PAIRED"] == "Yes" ||
	    vars["PAIRED"] == "yes" ||
	    vars["PAIRED"] == "Y" ||
	    vars["PAIRED"] == "y" ||
	    vars["PAIRED"] == "TRUE" ||
	    vars["PAIRED"] == "True" ||
	    vars["PAIRED"] == "true" ||
	    vars["PAIRED"] == "T" ||
	    vars["PAIRED"] == "t"
	    ) {
		// Use the specified alignment tool
		if (vars["ALIGNERTOOL"] == "BWAMEM") {
			// Directly return the .sam file created from bwa_mem
			outputSam, alignedLog, tmpalignedLog = bwa_mem_logged(vars["BWAEXE"], reads[0], reads[1], vars["BWAINDEX"], 
					    [vars["BWAMEMPARAMS"]], threads, rgheader, sampleName
					   ) ;
		} 
		else { // Novoalign is the default aligner
			// Directly return the .sam file created from novoalign
			outputSam, alignedLog, tmpalignedLog = novoalign_logged(vars["NOVOALIGNEXE"], reads[0], reads[1],
					      vars["NOVOALIGNINDEX"], [vars["NOVOALIGNPARAMS"]], threads, rgheader, sampleName
					     ) ; 
		}
	}
	else {
		// Use the specified alignment tool
		if (vars["ALIGNERTOOL"] == "BWAMEM") {
			// Directly return the .sam file created from bwa_mem
			outputSam, alignedLog, tmpalignedLog = bwa_mem_logged(vars["BWAEXE"], reads[0], vars["BWAINDEX"],
					    [vars["BWAMEMPARAMS"]], threads, rgheader, sampleName
					   ) ; 
		}
		else { // Novoalign is the default aligner								 
			// Directly return the .sam file created from novoalign
			outputSam, alignedLog, tmpalignedLog = novoalign_logged(vars["NOVOALIGNEXE"], reads[0], 
					      vars["NOVOALIGNINDEX"], [vars["NOVOALIGNPARAMS"]], threads, rgheader, sampleName 
					     ) ; 
		}
	}
}

/*************************
 Run function
**************************/

// For now, we will just feed in the lines array as the starting point. If fastq quality control STAGES are
//   implemented, this would probably need be to be altered
(file outputBam[]) alignRun(string lines[], string vars[string], file failLog) {
	foreach sample, index in lines {
		/*****
		Parse sample specific information and construct RG header
		*****/
		string sampleInfo[] = split(sample, " ");
		string sampleName = sampleInfo[0];
		string rgheader = sprintf("@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s\\tCN:%s", sampleName,
					  vars["SAMPLELB"], vars["SAMPLEPL"], sampleName, sampleName, vars["SAMPLECN"] 
					 );

		/****************************************************************************
		Alignment and Deduplication (per sample)
		****************************************************************************/

		/*****
		Create the sample output directories
		*****/

		string LogDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/logs/");
		string AlignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/align/");
		string RealignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/"); 
		string VarcallDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/");
		string tmpLogDir = strcat(vars["TMPDIR"], "/timinglogs/" );

		if (vars["ALIGN_STAGE"] == "Y" ||
		    vars["ALIGN_STAGE"] == "Yes" ||
		    vars["ALIGN_STAGE"] == "YES" ||
		    vars["ALIGN_STAGE"] == "y" ||
		    vars["ALIGN_STAGE"] == "yes" ||
		    vars["ALIGN_STAGE"] == "End" ||
		    vars["ALIGN_STAGE"] == "end" ||
		    vars["ALIGN_STAGE"] == "E" ||
		    vars["ALIGN_STAGE"] == "e"
		   ) {

			mkdir(LogDir) =>
			mkdir(AlignDir) =>
			mkdir(RealignDir) => /* Explicit waiting added to try to fix mkdir collision: Issue #22 */
			mkdir(VarcallDir) =>
			void mkdirSignal = mkdir(tmpLogDir);

			/*****
			Create output file handles
			*****/
			file alignedbam < strcat(AlignDir, sampleName, ".noDedups.bam") >;
	
			// These are temporary files: If piping is implemented, they would not be needed.
			file alignedsam < strcat(vars["TMPDIR"], "/align/", sampleName, ".noDedups.sam") >;
			file tmpsamtoolsLog < strcat(tmpLogDir, sampleName, "_samtools.log")>;
	
			/*****
			Alignment
			*****/
			if (vars["PAIRED"] == "1" ||
			    vars["PAIRED"] == "YES" ||
			    vars["PAIRED"] == "Yes" ||
			    vars["PAIRED"] == "yes" ||
			    vars["PAIRED"] == "Y" ||
			    vars["PAIRED"] == "y" ||
			    vars["PAIRED"] == "TRUE" ||
			    vars["PAIRED"] == "True" ||
			    vars["PAIRED"] == "true" ||
			    vars["PAIRED"] == "T" ||
			    vars["PAIRED"] == "t"
			   ) {
				string read1 = sampleInfo[1];
				string read2 = sampleInfo[2];
				string reads[] = [read1, read2];
				wait (mkdirSignal) {
					alignedsam = alignReads(vars, sampleName, reads, rgheader);
				}
			}
			else {
				string read1 = sampleInfo[1];
				string reads[] = [read1];
				wait (mkdirSignal) {
					alignedsam = alignReads(vars, sampleName, reads, rgheader);
				}
			}

			int threads = string2int(vars["CORES_PER_NODE"]) %/ string2int(vars["PROGRAMS_PER_NODE"]);
			alignedbam, tmpsamtoolsLog = samtools_view_logged(vars["SAMTOOLSEXE"], alignedsam, threads, ["-u"], sampleName);
	
			// Verify alignment was successful
			if ( checkBam(vars, alignedbam) ) {	
				
				/*****
				Add the raw aligned bam to the output array
				******/
				outputBam[index] = alignedbam;
						}
			// If the alignment process fails, write a message to the failure log file
			else {
				string message = strcat("FAILURE: ", filename(alignedbam), " contains no alignments. ", 
							"Check ", LogDir, sampleName, "_Alignment.log for details.\n"
						       );
				append(failLog, message) =>
				exitIfFlagGiven(vars, message);
			}		
		}
		// If this STAGE is to be skipped
		else {
			string outputFile = strcat(AlignDir, sampleName, ".noDedups.bam");

			if (file_exists(outputFile)) {
				outputBam[index] = input(outputFile);
			}
			else {
				string message_string = strcat("ERROR: ", outputFile, " not found. Did you set ",
							       "ALIGN_STAGE to 'N' by accident?\n"
							      );
				append(failLog, message_string) =>
				exitIfFlagGiven(vars, message_string);
			}
		}
	}
}
