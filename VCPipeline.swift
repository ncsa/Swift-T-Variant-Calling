// The way to run this script is:
// swift-t -r $PWD/pipelinefunctions VCcallingPipeline.swift --runfile=<Insert Runfile Location>

/* 
Note:
You may need additional logging, for example, how compilation works 
(add `-L <stc log file name>` to the command above); or how 
implementaion is mapped onto servers (you can enable turbine logging 
by `export TURBINE_LOG=1` before calling the command above)
*/

/**********************************************************************************************************************
********************************************   SUMMARY OF WORKFLOW   **************************************************
***********************************************************************************************************************

*****************
Helper functions
*****************

(file outputSam) alignReads(string read1, string read2, string rgheader);
	- Performs an alignment on the input files using the aligner specifed in the runfile

(file dedupSortedBam) markDuplicates(string sampleName, file alignedSam, file alignedBam);
	- Sorts and marks the duplicates of the input bam file using the deduplication to specified in the runfile 

(file realignedbam) realignBam(string sampleName, string chr, string realparms[], file inputBam);
	- Performs realignment utilizing GATK's RealignerTargetCreator and IndelRealigner

(file recalibratedbam) recalibrateBam(string sampleName, string chr, file inputBam, string recalparmsindels[]);
	- Performs recalibration utilizing GATK's BaseRecalibrator and PrintReads

() checkBam(file bamFile, string sampleName, string workflowStep);
	- Checks whether the given bam file is empty. If so, throw an error.

*****************************
Main loop through all samples
*****************************

// Because Swift/T syntax requires that control flow be specified with wait statements, waiting requires an output
     from the function to be waited on, and the workflow requires joint genotyping to be completed only when
     individual variant calling of all the samples has finished, it makes sense to nest the loop through all of the
     samples within a function with an output.

// The following is pseudocode that outlines the inner workings of the main loop through the samples:

(boolean finished) mainSampleLoop(boolean alignOnly) {
	foreach sample in Samples {
		
		- Create sample output directories
		- Create output file handles
		- alignReads (and verify alignment success)
		- markDuplicates (and verify dedup success)
		- QC on deduplicated file

		if (!alignOnly) {   // If the runfile did not specify ANALYSIS=ALIGN_ONLY

			foreach chromosome in sample {
				
				- Create output file handles
				- Separate dedupped and sorted bam by chromosome (and verify separation success)
				- Gather list of reference recalibration files and retrieve info from them
			
				if (Variant calling with realignment) {  // If runfile specifies ANALYSIS=VC_REALIGN
					- realignBam
					- recalibrateBam
				}
				else {		// If runfile specifies ANALYSIS=<nothing> or ANALYSIS=VC_NO_REALIGN 
					- recalibrateBam
				}

				- Call variants using HaplotypeCaller
			}
		}
	}
	finished = true;
}	

***************************
Implementation of Workflow
***************************

// Some parts of the workflow are not compartmentalized into functions, but are instead located within global scope
// I will refer to these as < SECTION X: What section X does > 

************
SECTION ONE: Parse Files
************
- Parse the variables from the runfile
- Parse the samples file
- Copy the runfile and samples file into <OutputFolder>/<DeliveryFolder>/docs for documentation purposes

***********
SECTION TWO: Immediate Error Checking 
***********
- Check for analysis variable errors
	- If the analysis variable is empty, warn that the default (variant calling without realignment) is being
	  used, and continue.
	- If the analysis variable is neither empty nor any of the permitted options (ALIGN_ONLY, VC_ALIGN, 
	  VC_NO_ALIGN), throw an error outlining what the acceptable inputs are and kill the workflow. This is to 
	  prevent confusing the end user. For example, one might think that ANALYSIS=ALIGN will start alignment only, 
	  but if this check was not here, it would default to full variant calling without alignment. This check also
	  prevents typos from causing the wrong analysis from being executed (Imagine setting ANALYSIS to
	  'align_only' or 'ALIIGN_ONLY'). 

***********
SECTION THREE: The Main Loop and Joint Genotyping 
***********
- Call the main loop
- (If not alignOnly) wait for main loop to end, then perform joint genotyping

*********************************************************************************************************************** 
***********************************************************************************************************************
**********************************************************************************************************************/

/*
Swift/T does not appear to provide an easy way to append to files that are already written, as the write method for
files requires a variable assignment and multiple assignment is not allowed.

A naive workaround might be to simply copy the contents of the original file to a string, add to it, and create 
a new updated file after deleting the original. However, this approach may lead to concurrency bugs (when multiple
processes try to add to the same file simultaneously, one may delete the original just as another tries to access
it) 

This is the precise problem encountered with the workflows QC file. As an easy solution, we simply test that each
process finished correctly with an assert method that prints an error message to stdout if the process is not
successful. Checking whether the workflow succeeded for each sample is as simple as checking stdout.
*/

import sys;
import io;
import files;
import string;
import unix;
import assert;
import pipelinefunctions.align_dedup;
import pipelinefunctions.realign_varcall_by_chr;
import pipelinefunctions.merge_vcf;
import pipelinefunctions.joint_vcf;
import pipelinefunctions.miscellaneous; 

/****************************************************************************
Helper functions (Easily handles alignment and duplicate marker choices)
****************************************************************************/

/*****
Alignment
*****/

(file outputSam) alignReads(string read1, string read2, string rgheader) {
	/*
	 This function returns a .sam file because samblaster requires it
	 To minimize memory usage, delete the .sam file after a .bam file is made from it
	*/

	// Use the specified alignment tool
	if (vars["ALIGNERTOOL"] == "BWAMEM") {
		// Directly return the .sam file created from bwa_mem
		outputSam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], 
				    [vars["BWAMEMPARAMS"]], string2int(vars["PBSCORES"]), rgheader
				   );
	} 
	else { // Novoalign is the default aligner
		// Directly return the .sam file created from novoalign
		outputSam =  novoalign(vars["NOVOALIGNDIR"], read1, read2, vars["NOVOALIGNINDEX"],
				       [vars["NOVOALIGNPARAMS"]], string2int(vars["PBSCORES"]), rgheader
				      );
	}
}

/*****
Mark Duplicates
*****/

(file dedupSortedBam) markDuplicates(string sampleName, file alignedSam, file alignedBam) {
	/*
	The input file used depends on the dedup program being used:
		alignedSam => Samblaster
		alignedBam => Picard or Novosort
	*/

	/*
	Note: this local directory location variable is necessary, as the directory handles
	defined in the main loop of the workflow are not in global scope
	*/

	string AlignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/align/");

	if (vars["MARKDUPLICATESTOOL"] == "SAMBLASTER") {
		file dedupsam < strcat(vars["TMPDIR"], "/align/", sampleName, ".wdups.sam") >;
		file dedupbam < strcat(AlignDir, sampleName, ".wdups.bam") >;
		
		// Mark Duplicates
		dedupsam = samblaster(vars["SAMBLASTERDIR"], alignedSam);
		dedupbam = samtools_view(vars["SAMTOOLSDIR"], dedupsam, string2int(vars["PBSCORES"]), ["-u"]);
		// Delete the dedupsam file once dedupbam has been created
		rm(dedupsam);
		
		// Sort
		dedupSortedBam = novosort(vars["NOVOSORTDIR"], dedupbam, vars["TMPDIR"],
					  string2int(vars["PBSCORES"]), ["--compression", "1"]
					 );
	}
	else if (vars["MARKDUPLICATESTOOL"] == "PICARD") {
		// Picard is unique in that it has a metrics file
		file metricsfile < strcat(AlignDir, sampleName, ".picard.metrics") >;
		file alignedsortedbam < strcat(AlignDir, sampleName, ".nodups.sorted.bam") >;
		
		// Sort
		alignedsortedbam = novosort(vars["NOVOSORTDIR"], alignedBam, vars["TMPDIR"],
					    string2int(vars["PBSCORES"]), []
					   );
		// Mark Duplicates
		dedupSortedBam, metricsfile = picard(vars["JAVADIR"], vars["PICARDDIR"],
						     vars["TMPDIR"], alignedsortedbam
						    ); 
	}
	else {	//Novosort is the default duplicate marker
		
		// Sort and Mark Duplicates in one step
		dedupSortedBam = novosort(vars["NOVOSORTDIR"], alignedBam, vars["TMPDIR"],
					  string2int(vars["PBSCORES"]), ["--markDuplicates"]
					 );
	}
}

/*****
Realignment
*****/

(file realignedbam) realignBam(string sampleName, string chr, string realparms[], file inputBam) { 
	
	printf("\n\n\n\n\nThis printed from realignBam\n\n\n\n\n\n\n");
	string RealignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/");

	file intervals < strcat(RealignDir, sampleName, ".", chr, ".realignTargetCreator.intervals") >;
		
	// The inputBam should be indexed before this function is called				
	intervals = RealignerTargetCreator(vars["JAVADIR"], vars["GATKDIR"],
					   strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),
					   inputBam, string2int(vars["PBSCORES"]), realparms		   
					  );
	realignedbam = IndelRealigner(vars["JAVADIR"], vars["GATKDIR"],
				      strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),
				      inputBam, realparms, intervals					   
				     );
	checkBam(realignedbam, sampleName, "realignment");
}

/*****
Recalibration
*****/

(file recalibratedbam) recalibrateBam(string sampleName, string chr, file inputBam, string recalparmsindels[]) {
	string RealignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/");				    
 
	file recalreport < strcat(RealignDir, sampleName, ".", chr, ".recal_report.grp") >;     	

	// The inputBam should be indexed before this function is called
	recalreport = BaseRecalibrator(vars["JAVADIR"], vars["GATKDIR"],
				       strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]), inputBam,
				       string2int(vars["PBSCORES"]), recalparmsindels,
				       strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"])
				      );
	recalibratedbam = PrintReads(vars["JAVADIR"], vars["GATKDIR"],
				     strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]), inputBam,
				     string2int(vars["PBSCORES"]), recalreport
				    );
	checkBam(recalibratedbam, sampleName, "recalibration");
}

/*****
BAM File Verification
*****/

() checkBam(file bamFile, string sampleName, string workflowStep) {
	// For some reason, Swift/T would try to use samtools_view2 even before bamFile is ready
	// Added explicit wait command
	wait (bamFile) {
		int alignNum = samtools_view2(vars["SAMTOOLSDIR"], filename(bamFile));
	
		string message = strcat("FAILURE: ", filename(bamFile), " contains no alignments. ",
					"Likely error during the ", workflowStep, " step when processing",
					sampleName, "."
			       	       );
		assert(alignNum > 0, message);
	}
}

(file samplesOut[]) mainSampleLoop(boolean alignOnly) {
	/*
	Output array is filled with the names of the final files for each sample in main loop
	Files are in gvcf format (list will be empty if alignment only) 
	*/

	/****************************************************************************
	Main loop begins
	*****************************************************************************/
	foreach sample, index in sampleLines {

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

		string AlignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/align/");
		mkdir(AlignDir);

		string VarcallDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/");
		mkdir(VarcallDir);
													   
		string RealignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/");
		mkdir(RealignDir);

		/*****
		Create output file handles
		*****/
		file alignedbam < strcat(AlignDir, sampleName, ".nodups.bam") >;
		file dedupsortedbam < strcat(AlignDir, sampleName, ".wdups.sorted.bam") >;		
											   
		// These are not specifically defined!
		file flagstats < strcat(AlignDir, sampleName, ".wdups.sorted.bam", ".flagstats") >;

		// These are temporary files: If piping is implemented, they would not be needed.
		file alignedsam < strcat(vars["TMPDIR"], "/align/", sampleName, ".nodups.sam") >;

		/*****
		Alignment
		*****/
		alignedsam = alignReads(read1, read2, rgheader);
		alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, string2int(vars["PBSCORES"]), ["-u"]);

		// Verify alignment was successful
		checkBam(alignedbam, sampleName, "alignment");
															
		/*****
		Deduplication
		*****/
		dedupsortedbam = markDuplicates(sampleName, alignedsam, alignedbam);
															
		// Wait until deduplication is done before deleting the raw aligned sam (in case samblaster was used)
		wait(dedupsortedbam) {
			rm(alignedsam);
		}
															
		// Verify deduplication was successful
		checkBam(dedupsortedbam, sampleName, "deduplication");
															
		/*****
		Quality control of deduplicated file
		*****/
															
		flagstats = samtools_flagstat(vars["SAMTOOLSDIR"], dedupsortedbam);
															
		string stat[] = file_lines(flagstats);
		tot_mapped =  split(stat[4], " ")[0];
		tot_reads = split(stat[0], " ")[0];
		tot_dups = split(stat[3], " ")[0];
															
		perc_dup= string2float(tot_dups) * 100 / string2float(tot_reads);
		perc_mapped= string2float(tot_mapped) * 100 / string2float(tot_reads);
															
		// Message cutoff information
		string cutoff_info = strcat(sampleName,
					    "\tPercentDuplication=", perc_dup,
					    ";DuplicationCutoff=", vars["DUP_CUTOFF"],
					    "\tPercentMapped=", perc_mapped, ";MappingCutoff=", vars["MAP_CUTOFF"]
					   );
		// If both % duplicated and % mapped meet their cutoffs
		if ( perc_dup < string2float(vars["DUP_CUTOFF"]) && perc_mapped > string2float(vars["MAP_CUTOFF"]) ) {
			/*
			 SUCCESS
			*/
			printf(strcat("QC-Test SUCCESS: ", cutoff_info));
		}
		else {
			/*
			FAILURE
			*/
			// Trigger failure and print error message
			assert(false, strcat("QC-Test FAILURE: ", cutoff_info));
		}
		
		/*******************************************************************
		IF ALIGN ONLY, THEN THIS IS THE END OF THE WORKFLOW FOR THIS SAMPLE
		*******************************************************************/

		if (!alignOnly) {

			/****************************************************************************
			Split by Chromosome
			****************************************************************************/

			indices = split(vars["CHRNAMES"], ":");							 
			
			file gvcfSample < strcat(VarcallDir, sampleName, ".raw.g.vcf") >;
			// List of chrgvcf files
			file gvcfChrCollection[];
															
			foreach chr, chrIndex in indices {
				/*****										  
				Create output file handles 
				*****/										  
				file chrdedupsortedbam < strcat(RealignDir, sampleName, ".", chr,		       
								".wdups.sorted.bam"				     
							       ) >;						     
				file recalibratedbam < strcat(RealignDir, sampleName, ".", chr, ".recalibrated.bam") >; 

				// Put gvcfvariant handle in the chrCompleted array
				file gvcfSampleChr < strcat(VarcallDir, sampleName, ".", chr, ".raw.g.vcf") >;
				// Temporary file								       
				file recalfiles < strcat(vars["TMPDIR"], "/", sampleName, ".", chr,		     
							 ".recal_foundfiles.txt"					
							) >;							    
															
				/*****										  
				Separate dedupped and sorted bam by chromosome 
				*****/										  
				chrdedupsortedbam = samtools_view(vars["SAMTOOLSDIR"], dedupsortedbam,		  
								  string2int(vars["PBSCORES"]), [strcat(chr)]	   
								 ) =>						   
															
				// Verify chromosome separation worked						  
				checkBam(chrdedupsortedbam, sampleName, "realignment");  

				string recalparmsindels[];							      
				string realparms[];								     
				wait (chrdedupsortedbam) {							      
					/*****									  
					Gather list of all reference recalibration files,			       
					then retrieve indels and parameters from them				   
					*****/									  
					recalfiles = find_files(							
						strcat(vars["REFGENOMEDIR"], "/", vars["INDELDIR"], "/"),	       
						strcat("*", chr, ".vcf" )					       
						       	       );

					recalparmsindels = split(						       
						trim(replace_all(read(sed(recalfiles, "s/^/--knownSites /g"	     
									 )					      
								     ), "\n", " ", 0				    
								)						       
					       	    ), " "
								);
															
					realparms = split(							      
						trim(replace_all(read(sed(recalfiles, "s/^/-known /g"		   
									 )					      
								     ), "\n", " ", 0				    
								)						       
						    ), " "							      
							  );						    
				}										       
															
				/****************************************************************************	   
				Realignment and/or Recalibration							
				****************************************************************************/	   
															
				if (vars["ANALYSIS"] == "VC_REALIGN") {						 
					/*****									  
					Realignment and Recalibration						   
					*****/
					file realignedbam < strcat(RealignDir, sampleName, ".", chr,
								   ".realigned.bam"
								  ) >;
					// Wait for index to get created before moving on
					samtools_index(vars["SAMTOOLSDIR"], chrdedupsortedbam) =>
					realignedbam = realignBam(sampleName, chr, realparms, chrdedupsortedbam);
					recalibratedbam = recalibrateBam(sampleName, chr, realignedbam,
									 recalparmsindels
									); 
															
				}										       
				else {										  
					/*****
					Recalibration								   
					*****/									  
					// Wait for index to get created before moving on
					samtools_index(vars["SAMTOOLSDIR"], chrdedupsortedbam) => 
					recalibratedbam = recalibrateBam(sampleName, chr,
								         chrdedupsortedbam, recalparmsindels
								        );
				}										       

				int ploidy;									     
				if ( chr=="M" ) { ploidy = 1; }							 
				else { ploidy = 2; }	       

				gvcfSampleChr = HaplotypeCaller(vars["JAVADIR"], vars["GATKDIR"],			
							        strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),    
							        recalibratedbam,					 
							        strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),	
							        string2int(vars["PBSCORES"]), ploidy, chr		
						      	       ) =>
				/*
				The wait signal ('=>') in the previous command must be used as the name of the
				output file is mapped when gvcfSampleChr is created, not when HaplotypeCaller is
				finished running.

				Without this, the gvcfChrCollection array assigns gvcfSampleChr to a location after
				file mapping, not after it is written to
				*/
				gvcfChrCollection[chrIndex] = gvcfSampleChr =>
				rm(recalfiles);
			} // end the loop for all chromosomes

			/*
			Get the names of each of the chromosome split sample files, add a --variants flag that the
			beginning of them, and concatenate them together to get a string that can be fed into
			CombineGVCFs
			
			Example: "--variant sample1.chr1.g.vcf --variant sample1.chr2.g.vcf"
			*/
			string chrSampleArray[];
			foreach f, ind in gvcfChrCollection {
				// Adds a string like "--variant sample2.vcf" to the string array
				chrSampleArray[ind] = strcat("--variant ", filename(f));
                	}
                	// Concatenates everything in the array with a space in between each entry
			string chrSampleVars = join(chrSampleArray, " ");

			/*
			Wait for end of foreach loop, then regroup all of the chromosome files
			and put them in the samplesOut array
			*/
			gvcfSample = CombineGVCFs(vars["JAVADIR"], vars["GATKDIR"], 
						  strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),
						  strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),
						  chrSampleVars
				 		 ) =>
			samplesOut[index] = gvcfSample;
		}
	} // end the loop for all samples
}

/****************************************************************************					      
SECTION ONE: Parse Files
*****************************************************************************/					     

//read-in the runfile with argument --runfile=<runfile_name>							       
argv_accept("runfile"); // return error if user supplies other flagged inputs					      
string configFilename = argv("runfile");

file configFile = input_file(configFilename);									      
string configFileData[] = file_lines(configFile);

string vars[string] = getConfigVariables(configFileData);

file sampleInfoFile = input_file(vars["SAMPLEINFORMATION"]);							       
string sampleLines[] = file_lines(sampleInfoFile);

// Create the output and temporary directories
mkdir(vars["OUTPUTDIR"]);
mkdir(vars["TMPDIR"]);

// Copy the runfile and sampleInfoFile to the docs directory for documentation purposes
file docRunfile < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/docs/", basename_string(configFilename)
			) > = configFile;

file docSampleInfo < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/docs/",				      
			    basename_string(filename(sampleInfoFile))						      
			   ) > = sampleInfoFile;

// Create quality control file											     
//file qcfile < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/docs/QC_test_results.txt") >;		     

/****************************************************************************					      
SECTION TWO: Immediate Error Checking
*****************************************************************************/

/*****														     
ANALYSIS Variable Errors
(The 'callVariants' function relies on this conditional statement being evaluated before its use)
*****/														     

// If the ANALYSIS variable is empty, warn that the default (variant calling without realignment) is being used	    
if (vars["ANALYSIS"] == "") {											      
	printf("WARNING: No analysis choice specified: variant calling with no realignment is the default choice\n");   
}														       
// If anything other than the available options (or nothing) is entered, provide an error and terminate		 
else if (vars["ANALYSIS"] != "ALIGN_ONLY" &&									    
	 vars["ANALYSIS"] != "VC_REALIGN" &&									    
	 vars["ANALYSIS"] != "VC_NO_REALIGN"									    
	) {													     
	string message = strcat("FAILURE: Analysis option <", vars["ANALYSIS"], "> is not a vaild entry.\n",	    
				"\tPlease enter one of the following:\n",					       
				"\t\tFor alignment only: ALIGN_ONLY\n",						 
				"\t\tFor variant calling with realignment: VC_REALIGN\n",			       
				"\t\tFor variant calling without realignment: VC_NO_REALIGN\n"			  
			       );										       
	assert(false, message);											 
}	      

/*														      
ADD SAFETY CHECKS BEFORE THE MAIN LOOP BEGINS MAKING SURE THAT THE NECESSARY INDEX FILES ARE PRESENT FOR THE TYPE OF    
ANALYSIS BEING CONDUCTED												
															
FOR EXAMPLE:													    
	If stage of analysis listed includes realignment, and the reference indel and param files cannot be found,      
	terminate the process now instead of waiting until the realign step to check for the necessary files	    
															
// assert(strlen(recalparmsindels)>1 ,										  
	  strcat("no indels were found for ", chr, " in this folder", vars["REFGENOMEDIR"]/vars["INDELDIR"] ));	 
// assert(strlen(realparms)>1 ,											 
	  strcat("no indels were found for ", chr, "in this folder", vars["REFGENOMEDIR"]/vars["INDELDIR"] ));	  
															
*/	

/****************************************************************************
SECTION THREE: The Main Loop and Joint Genotyping
****************************************************************************/

// Will be true only if we only want the alignments
boolean alignOnlyCheck = vars["ANALYSIS"] == "ALIGN_ONLY"; 

/******************************
Main Loop on Samples is Called
*******************************/
// If alignment only, this array will actually be empty

file samplesOutput[] = mainSampleLoop(alignOnlyCheck);

if (!alignOnlyCheck) {
	wait(samplesOutput) {
		
		/***************
		Joint Genotyping
		****************/
		file jointVCF < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"],
				       "/jointVCFs/jointVCFcalled.vcf"
				      ) >;
		mkdir(strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/jointVCFs"));

		/*
		Construct a string with '--variant' flags before each sample to be in the joint genotyping
		Example: "--variant sample1.g.vcf --variant sample2.g.vcf"
		*/
		string variantSampleArray[];
		foreach item, index in samplesOutput {
			// Adds a string like "--variant sample1.vcf" to the string array
			variantSampleArray[index] = strcat("--variant ", filename(item));
		}
		// Concatenates everything in the array with a space in between each entry
		string variantFlags = join(variantSampleArray, " ");

		jointVCF = GenotypeGVCFs(vars["JAVADIR"], vars["GATKDIR"], strcat(vars["REFGENOMEDIR"],
					 "/", vars["REFGENOME"]), variantFlags
					); 
	}
}


