/*
*****************************
 Pseudocode of Run Function
*****************************
(file outBamArray[]) alignDedupRun() {
	foreach sample in samples {
		- deduce sample names from input bam file names
		if (ALIGN_DEDUP_STAGE variable == "Y") {
			**************************
			*** EXECUTE THIS STAGE ***
			**************************
			- Perform sorting and deduplication

			if (dedupsorting was successful) {
				- Perform quality control
			
				if (quality control passed) {
					*** SUCCESS ***
					- Add the dedupsorted bam file to the output array
				} else {
					*** FAILURE ***
					- Write an error message to the fail log
				}
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
				- Add the dedupsorted bam file to the output array
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
 Helper Function
********************************************************/

/***************
 Mark Duplicates
****************/
(file dedupSortedBam) markDuplicates(string vars[string], string sampleName, file alignedSam, file alignedBam) {
	/*
	The input file used depends on the dedup program being used:
		alignedSam => Samblaster
		alignedBam => Picard or Novosort
	*/

	int threads = string2int(vars["CORES_PER_NODE"]) %/ string2int(vars["PROGRAMS_PER_NODE"]);

	string LogDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/logs/");
	string AlignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/align/");
	string tmpLogDir = strcat(vars["TMPDIR"], "/timinglogs/" );

	if (vars["MARKDUPLICATESTOOL"] == "SAMBLASTER") {
		exec_check(vars["SAMBLASTEREXE"], "SAMBLASTEREXE");
		exec_check(vars["SAMTOOLSEXE"], "SAMTOOLSEXE");
		file dedupsam < strcat(vars["TMPDIR"], "/align/", sampleName, ".wDedups.sam") >;
		file dedupbam < strcat(AlignDir, sampleName, ".wDedups.bam") >;
		file samLog < strcat(LogDir, sampleName, "_SamblasterDedup.log") >; 
		file sortLog < strcat(LogDir, sampleName, "_Sort.log") >;
		file tmpsamblasterLog < strcat(tmpLogDir, sampleName, "_SamblasterDedup.log")>;	
		file tmpsamtoolsLog < strcat(tmpLogDir, sampleName, "_samtoolsDedup.log")>;	
		file tmpnovosortLog < strcat(tmpLogDir, sampleName, "_NovoSortDedup.log")>;	

		// Mark Duplicates
		dedupsam, samLog, tmpsamblasterLog = samblaster_logged(vars["SAMBLASTEREXE"], alignedSam, sampleName);
		dedupbam, tmpsamtoolsLog = samtools_view_logged(vars["SAMTOOLSEXE"],
								dedupsam,
								threads,
								["-u"],
								sampleName
							       ) =>
		// Delete the dedupsam file once dedupbam has been created (wait for dedupbam to be finished)
		rm(dedupsam);
		
		// Sort
		dedupSortedBam, sortLog, tmpnovosortLog = novosort_logged(vars["NOVOSORTEXE"],
									  dedupbam, 
									  vars["TMPDIR"],
						   			  threads, ["--compression", "1"], 
									  string2int(vars["NOVOSORT_MEMLIMIT"]),
									  sampleName
									 );
	}
	else if (vars["MARKDUPLICATESTOOL"] == "PICARD") {
		exec_check(vars["PICARDJAR"], "PICARDJAR");
		exec_check(vars["JAVAEXE"], "JAVAEXE");
		exec_check(vars["NOVOSORTEXE"], "NOVOSORTEXE");
		// Picard is unique in that it has a metrics file
		file metricsfile < strcat(AlignDir, sampleName, ".picard.metrics") >;
		file alignedsortedbam < strcat(AlignDir, sampleName, ".noDedups.sorted.bam") >;
		file picardLog < strcat(LogDir, sampleName, "_PicardDedup.log") >;
		file sortLog < strcat(LogDir, sampleName, "_Sort.log") >;
		file tmpnovosortLog < strcat(tmpLogDir, sampleName, "_NovoSortDedup.log")>;	
		file tmppicardLog < strcat(tmpLogDir, sampleName, "_PicardDedup.log")>;	

		// Sort
		alignedsortedbam, sortLog, tmpnovosortLog = novosort_logged(vars["NOVOSORTEXE"],
									    alignedBam,
									    vars["TMPDIR"],
									    threads,
									    [],
									    string2int(vars["NOVOSORT_MEMLIMIT"]),
									    sampleName
									   );
		// Mark Duplicates
		dedupSortedBam, picardLog, metricsfile, tmppicardLog = picard_logged(vars["JAVAEXE"],
										     vars["PICARDJAR"],
										     vars["TMPDIR"],
										     alignedsortedbam,
										     sampleName
										    );
	}
	else {	//Novosort is the default duplicate marker
		exec_check(vars["NOVOSORTEXE"], "NOVOSORTEXE");
		file novoLog < strcat(LogDir, sampleName, "_NovosortDedup.log") >;
		file tmpnovosortLog < strcat(tmpLogDir, sampleName, "_NovoSortDedup.log")>;

		// Sort and Mark Duplicates in one step
		dedupSortedBam, novoLog, tmpnovosortLog = novosort_logged(vars["NOVOSORTEXE"],
									  alignedBam,
									  vars["TMPDIR"],
									  threads,
									  ["--markDuplicates"],
									  string2int(vars["NOVOSORT_MEMLIMIT"]),
									  sampleName
									 );
	}
}

/*************************
 Run function
**************************/

(file outputBam[]) dedupSortRun(file inputBams[], string vars[string], file failLog) {
	foreach inputBam, index in inputBams {

		/*******
		Get the sample name
		********/
		string baseName = basename_string(filename(inputBam)); 
		// Remove ".noDedups.bam" extension
		string sampleName = substring(baseName, 0, strlen(baseName) - 13);

		string LogDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/logs/");
		string AlignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/align/");

		/****************************************************************************
		Deduplication (per sample)
		****************************************************************************/
		if (vars["DEDUP_SORT_STAGE"] == "Y" ||
		    vars["DEDUP_SORT_STAGE"] == "Yes" ||
		    vars["DEDUP_SORT_STAGE"] == "YES" ||
		    vars["DEDUP_SORT_STAGE"] == "y" ||
		    vars["DEDUP_SORT_STAGE"] == "yes" ||
		    vars["DEDUP_SORT_STAGE"] == "End" ||
		    vars["DEDUP_SORT_STAGE"] == "end" ||
		    vars["DEDUP_SORT_STAGE"] == "E" ||
		    vars["DEDUP_SORT_STAGE"] == "e"
		   ) {	
			/*****
			Deduplication
			*****/
			// Locate the alignedsam file, a leftover in temp from the previous stage
			file alignedsam = input(strcat(vars["TMPDIR"], "/align/", sampleName, ".noDedups.sam"));
			
			file dedupsortedbam < strcat(AlignDir, sampleName, ".wDedups.sorted.bam") >; 
			dedupsortedbam = markDuplicates(vars, sampleName, alignedsam, inputBam);
	
			// Verify deduplication was successful
			if ( checkBam(vars, dedupsortedbam) ) {
				// This will delete the raw aligned sam only after dedupsortedbam is written
				rm(alignedsam); 

				/*****
				Quality control of deduplicated file
				*****/
	
				file flagstats < strcat(AlignDir, sampleName, ".wDedups.sorted.bam", ".flagstats") >;
	
				flagstats = samtools_flagstat(vars["SAMTOOLSEXE"], dedupsortedbam);

				string stat[] = file_lines(flagstats);
				tot_mapped =  split(stat[4], " ")[0];
				tot_reads = split(stat[0], " ")[0];
				tot_dups = split(stat[3], " ")[0];
	
				perc_dup= string2float(tot_dups) * 100 / string2float(tot_reads);
				perc_mapped= string2float(tot_mapped) * 100 / string2float(tot_reads);
	
				// Message cutoff information
				string cutoff_info = strcat(sampleName, "\tPercentDuplication=", perc_dup,
							    ";DuplicationCutoff=", vars["DUP_CUTOFF"],
							    "\tPercentMapped=", perc_mapped,
							    ";MappingCutoff=", vars["MAP_CUTOFF"]
							   );
				// If both % duplicated and % mapped meet their cutoffs
				if ( perc_dup < string2float(vars["DUP_CUTOFF"]) && 
				     perc_mapped > string2float(vars["MAP_CUTOFF"])
				   ) {
					/*
					 SUCCESS
					*/
					printf(strcat("QC-Test SUCCESS: ", cutoff_info));
					//Add the aligned and dedupped file to the output list
					outputBam[index] = dedupsortedbam;
				}
				else {
					/*
					FAILURE
					*/
					string m = strcat("FAILURE: ", filename(dedupsortedbam),
							  " failed the QC test: ", cutoff_info
							 );
					append(failLog, m) =>
					exitIfFlagGiven(vars, m);
				}
			}
			// If the deduplication process fails, write a message to the failure log file
			else {
				string mssg = strcat("FAILURE: ",filename(dedupsortedbam), " contains no alignments. ",
						     "Check the log files within ", LogDir, 
						     sampleName, " for details.\n"
						    );
				append(failLog, mssg) =>
				exitIfFlagGiven(vars, mssg);
			}
		}
		// If this STAGE is to be skipped
		else {
			string outputFile = strcat(AlignDir, sampleName, ".wDedups.sorted.bam");

			if (file_exists(outputFile)) {
				outputBam[index] = input(outputFile);
			}
			else {
				string message_string = strcat("ERROR: ", outputFile, " not found. Did you set ",
							       "DEDUP_SORT_STAGE to 'N' by accident?\n"
							      );
				append(failLog, message_string) =>
				exitIfFlagGiven(vars, message_string);
			}
		}
	}
}
