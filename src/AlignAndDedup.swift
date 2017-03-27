/*

*****************************
 Pseudocode of Main Function
*****************************

(file outBamArray[]) alignDedupMain() {
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
import generalfunctions.general;

/********************************************************
 Helper Functions
********************************************************/

/*********
 Alignment
**********/
(file outputSam) alignReads(string vars[string], string sampleName, string read1, string read2, string rgheader) {
	/*
	This function returns a .sam file because samblaster requires it
	To minimize memory usage, delete the .sam file after a .bam file is made from it
	*/

	int threads = string2int(vars["PBSCORES"]) %/ string2int(vars["PROCPERNODE"]);

	// Log file
	string LogDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/logs/");
	file alignedLog < strcat(LogDir, sampleName, "_Alignment.log") >;
 	
	// Use the specified alignment tool
	if (vars["ALIGNERTOOL"] == "BWAMEM") {
		// Directly return the .sam file created from bwa_mem
		outputSam, alignedLog = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], 
				    [vars["BWAMEMPARAMS"]], threads, rgheader
				   );
	} 
	else { // Novoalign is the default aligner
		// Directly return the .sam file created from novoalign
		outputSam, alignedLog = novoalign(vars["NOVOALIGNDIR"], read1, read2, vars["NOVOALIGNINDEX"],
				      [vars["NOVOALIGNPARAMS"]], threads, rgheader
				     );
	}
}

/***************
 Mark Duplicates
****************/
(file dedupSortedBam) markDuplicates(string vars[string], string sampleName, file alignedSam, file alignedBam) {
	/*
	The input file used depends on the dedup program being used:
		alignedSam => Samblaster
		alignedBam => Picard or Novosort
	*/

	int threads = string2int(vars["PBSCORES"]) %/ string2int(vars["PROCPERNODE"]);

	string LogDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/logs/");
	string AlignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/align/");

	if (vars["MARKDUPLICATESTOOL"] == "SAMBLASTER") {
		file dedupsam < strcat(vars["TMPDIR"], "/align/", sampleName, ".wDedups.sam") >;
		file dedupbam < strcat(AlignDir, sampleName, ".wDedups.bam") >;
		file samLog < strcat(LogDir, sampleName, "_SamblasterDedup.log") >; 
		file sortLog < strcat(LogDir, sampleName, "_Sort.log") >;       	

		// Mark Duplicates
		dedupsam, samLog = samblaster(vars["SAMBLASTERDIR"], alignedSam);
		dedupbam = samtools_view(vars["SAMTOOLSDIR"], dedupsam, threads, ["-u"]);
		// Delete the dedupsam file once dedupbam has been created
		rm(dedupsam);
		
		// Sort
		dedupSortedBam, sortLog = novosort(vars["NOVOSORTDIR"], dedupbam, vars["TMPDIR"],
					  threads, ["--compression", "1"]
					 );
	}
	else if (vars["MARKDUPLICATESTOOL"] == "PICARD") {
		// Picard is unique in that it has a metrics file
		file metricsfile < strcat(AlignDir, sampleName, ".picard.metrics") >;
		file alignedsortedbam < strcat(AlignDir, sampleName, ".noDedups.sorted.bam") >;
		file picardLog < strcat(LogDir, sampleName, "_PicardDedup.log") >;
		file sortLog < strcat(LogDir, sampleName, "_Sort.log") >;
	
		// Sort
		alignedsortedbam, sortLog = novosort(vars["NOVOSORTDIR"], alignedBam, vars["TMPDIR"],
						     threads, []
						    );
		// Mark Duplicates
		dedupSortedBam, picardLog, metricsfile = picard(vars["JAVADIR"], vars["PICARDDIR"],
							 	vars["TMPDIR"], alignedsortedbam
							       ); 
	}
	else {	//Novosort is the default duplicate marker
		file novoLog < strcat(LogDir, sampleName, "_NovosortDedup.log") >;

		// Sort and Mark Duplicates in one step
		dedupSortedBam, novoLog = novosort(vars["NOVOSORTDIR"], alignedBam, vars["TMPDIR"],
					  threads, ["--markDuplicates"]
					 );
	}
}

/*************************
 Main function
**************************/

// For now, we will just feed in the lines array as the starting point. If fastq quality control STAGES are
//   implemented, this would probably need be to be altered
(file outputBam[]) alignDedupMain(string lines[], string vars[string], file failLog) {
	foreach sample, index in lines {
		/*****
		Parse sample specific information and construct RG header
		*****/
		string sampleInfo[] = split(sample, " ");
		string sampleName = sampleInfo[0];
		string read1 = sampleInfo[1];
		string read2 = sampleInfo[2];

		string rgheader = sprintf("@RG\tID:%s\tLB:%s\tPL:%s\tPU:%s\tSM:%s\tCN:%s", sampleName,
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

		if (vars["ALIGN_DEDUP_STAGE"] == "Y") {

			mkdir(LogDir);
			mkdir(AlignDir);
			mkdir(RealignDir);
			mkdir(VarcallDir);

			/*****
			Create output file handles
			*****/
			file alignedbam < strcat(AlignDir, sampleName, ".noDedups.bam") >;
	
			// These are temporary files: If piping is implemented, they would not be needed.
			file alignedsam < strcat(vars["TMPDIR"], "/align/", sampleName, ".noDedups.sam") >;
	
			/*****
			Alignment
			*****/
			alignedsam = alignReads(vars, sampleName, read1, read2, rgheader);
			int threads = string2int(vars["PBSCORES"]) %/ string2int(vars["PROCPERNODE"]);
			alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, threads, ["-u"]);
	
			// Verify alignment was successful
			if ( checkBam(vars, alignedbam) ) {
	
				/*****
				Deduplication
				*****/
			
				file dedupsortedbam < strcat(AlignDir, sampleName, ".wDedups.sorted.bam") >; 
				dedupsortedbam = markDuplicates(vars, sampleName, alignedsam, alignedbam);
	
				// Verify deduplication was successful
				if ( checkBam(vars, dedupsortedbam) ) {
	
					// This will delete the raw aligned sam only after dedupsortedbam is written
					rm(alignedsam); 
	
					/*****
					Quality control of deduplicated file
					*****/
	
					file flagstats < strcat(AlignDir, sampleName, ".wDedups.sorted.bam", ".flagstats") >;
	
					flagstats = samtools_flagstat(vars["SAMTOOLSDIR"], dedupsortedbam);
	
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
					     perc_mapped > string2float(vars["MAP_CUTOFF"]) ) {
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
						append(failLog, m);
					}
				}
				// If the deduplication process fails, write a message to the failure log file
				else {
					string mssg = strcat("FAILURE: ", filename(dedupsortedbam), " contains no alignments. ",
							     "Check the log files within ", LogDir, 
							     sampleName, " for details.\n"
							    );
					append(failLog, mssg);
				}
			}
			// If the alignment process fails, write a message to the failure log file
			else {
				string message = strcat("FAILURE: ", filename(alignedbam), " contains no alignments. ", 
							"Check ", LogDir, sampleName, "_Alignment.log for details.\n"
						       );
				append(failLog, message);
			}		
		}
		// If this STAGE is to be skipped
		else {
			string outputFile = strcat(AlignDir, sampleName, ".wDedups.sorted.bam");

			if (file_exists(outputFile)) {
				outputBam[index] = input(outputFile);
			}
			else {
				append(failLog, strcat("ERROR: ", outputFile, " not found. Did you set ",
							  "ALIGN_DEDUP_STAGE to 'N' by accident?\n"
						         )
				       );
			}
		}
	}
}
