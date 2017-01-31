// The way to run this script is:
// swift-t -r $PWD/pipelinefunctions VCcallingPipeline.swift --runfile=<Insert Runfile Location>

/* 
Note:
You may need additional logging, for example, how compilation works (add `-L <stc log file name>` to the command
above); or how implementaion is mapped onto servers (you can enable turbine logging by `export TURBINE_LOG=1` before
calling the command above)
*/

/**********************************************************************************************************************
********************************************   SUMMARY OF WORKFLOW   **************************************************
***********************************************************************************************************************

*****************
General functions
*****************

() checkBam(file bamFile, string sampleName, string workflowStep);						      
	- Checks whether the given bam file is empty. If so, throw an error.   

*****************
Helper functions  (Simplifies calls to the command-line tools of the pipeline)
*****************

(file outputSam) alignReads(string sampleName, string read1, string read2, string rgheader);
	- Performs an alignment on the input files using the aligner specifed in the runfile
	- The stderr from the alignment is put into a log file in the same location as the alignment output

(file dedupSortedBam) markDuplicates(string sampleName, file alignedSam, file alignedBam);
	- Sorts and marks the duplicates of the input bam file using the deduplication to specified in the runfile 
	- Writes log file(s) from dedup and sorting steps (Novosort usage only produces one log file, because sorting
	    and marking duplicates is a single step)

(file realignedbam, file targetLog, file intervals) realignBam(string sampleName, string realparms[], file inputBam);
	- Performs realignment utilizing GATK's RealignerTargetCreator and IndelRealigner			       
	- For each step, log files are written									  
															
(file outBam, file recalLog, file printLog, file report) recalibrateBam(string sampleName, file inputBam,
									string recalparmsindels[]
								       );      
	- Performs recalibration utilizing GATK's BaseRecalibrator and PrintReads				       
	- For each step, log files are written

(file gvcfFile) callVariants(string sampleName, file inputBam);							 
	- Call variants using HaplotypeCaller (Without needing chromosome information)

(file gvcfFile) callChrVariants(string sampleName, file inputBam, string chr);
	- Call variants using HaplotypeCaller (With chromosome information) 

*****************
Wrapper functions (Wraps conditional logic that is used multiple times in the workflow)
*****************

(file recalibratedBam) recalibrationWrapper(file prefix) {
	- Create file handles for recalibration step (prefix is at the beginning of file handles)

	if (variant calling with realignment) {
		- Create realignment specific file handles
		- realignBam
		- recalibrateBam
	}
	else {
		- recalibrateBam
	}

}

***************************												
Splitting and Processing by Chromosome (This loop is nested within the main loop)					  
***************************												
															   
(file gvcfSample) callVariantsByChromosome() {

	foreach chromosome in sample {									     
		- Create output file handles								       
		- Separate dedupped and sorted bam by chromosome						  
		- Gather list of reference recalibration files and retrieve info from them

		if (Variant calling with realignment) {									 
			- realignBam											    
			- recalibrateBam
		}
		else {													  
			- recalibrateBam
			- Call variants using HaplotypeCaller
		}													  

	Combine the variants called for each chromosome of the sample into a single file after each		
	  chromosome has been processed
}      

*****************************
Main loop through all samples
*****************************

// Because Swift/T syntax requires that control flow be specified with wait statements, waiting requires an output
     from the function to be waited on, and the workflow requires joint genotyping to be completed only when
     individual variant calling of all the samples has finished, it makes sense to nest the loop through all of the
     samples within a function with an output.

// Joint Genotyping only runs when there are no other possible additions to the samplesOut array, which is the output
//   of the mainSampleLoop

// The following is pseudocode that outlines the inner workings of the main loop through the samples:

(file samplesOut[]) mainSampleLoop(boolean alignOnly) {
	foreach sample in Samples {
		
		- Create sample output directories
		- Create output file handles
		- alignReads (and verify alignment success)
		- markDuplicates (and verify dedup success)
		- QC on deduplicated file

		if (!alignOnly) {   // If the runfile did not specify ANALYSIS=ALIGN_ONLY
			
			- Create array to store the variants for each chromosome of the sample

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

			Combine the variants called for each chromosome of the sample into a single file after each
			  chromosome has been processed
		}
	}
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

/********
Alignment
*********/

(file outputSam) alignReads(string sampleName, string read1, string read2, string rgheader) {
	/*
	 This function returns a .sam file because samblaster requires it
	 To minimize memory usage, delete the .sam file after a .bam file is made from it
	*/

	string AlignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/align/");

	// Log file
	file alignedLog < strcat(AlignDir, sampleName, "_Alignment.log") >;

	// Use the specified alignment tool
	if (vars["ALIGNERTOOL"] == "BWAMEM") {
		// Directly return the .sam file created from bwa_mem
		outputSam, alignedLog = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], 
				    [vars["BWAMEMPARAMS"]], string2int(vars["PBSCORES"]), rgheader
				   );
	} 
	else { // Novoalign is the default aligner
		// Directly return the .sam file created from novoalign
		outputSam, alignedLog = novoalign(vars["NOVOALIGNDIR"], read1, read2, vars["NOVOALIGNINDEX"],
				      [vars["NOVOALIGNPARAMS"]], string2int(vars["PBSCORES"]), rgheader
				     );
	}
}

/**************
Mark Duplicates
***************/

(file dedupSortedBam) markDuplicates(string sampleName, file alignedSam, file alignedBam) {
	/*
	The input file used depends on the dedup program being used:
		alignedSam => Samblaster
		alignedBam => Picard or Novosort
	*/

	string AlignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/align/");

	if (vars["MARKDUPLICATESTOOL"] == "SAMBLASTER") {
		file dedupsam < strcat(vars["TMPDIR"], "/align/", sampleName, ".wDedups.sam") >;
		file dedupbam < strcat(AlignDir, sampleName, ".wDedups.bam") >;
		file samLog < strcat(AlignDir, sampleName, "_SamblasterDedup.log") >; 
		file sortLog < strcat(AlignDir, sampleName, "_Sort.log") >;       	

		// Mark Duplicates
		dedupsam, samLog = samblaster(vars["SAMBLASTERDIR"], alignedSam);
		dedupbam = samtools_view(vars["SAMTOOLSDIR"], dedupsam, string2int(vars["PBSCORES"]), ["-u"]);
		// Delete the dedupsam file once dedupbam has been created
		rm(dedupsam);
		
		// Sort
		dedupSortedBam, sortLog = novosort(vars["NOVOSORTDIR"], dedupbam, vars["TMPDIR"],
					  string2int(vars["PBSCORES"]), ["--compression", "1"]
					 );
	}
	else if (vars["MARKDUPLICATESTOOL"] == "PICARD") {
		// Picard is unique in that it has a metrics file
		file metricsfile < strcat(AlignDir, sampleName, ".picard.metrics") >;
		file alignedsortedbam < strcat(AlignDir, sampleName, ".noDedups.sorted.bam") >;
		file picardLog < strcat(AlignDir, sampleName, "_PicardDedup.log") >;
		file sortLog < strcat(AlignDir, sampleName, "_Sort.log") >;
	
		// Sort
		alignedsortedbam, sortLog = novosort(vars["NOVOSORTDIR"], alignedBam, vars["TMPDIR"],
						     string2int(vars["PBSCORES"]), []
						    );
		// Mark Duplicates
		dedupSortedBam, picardLog, metricsfile = picard(vars["JAVADIR"], vars["PICARDDIR"],
							 	vars["TMPDIR"], alignedsortedbam
							       ); 
	}
	else {	//Novosort is the default duplicate marker
		file novoLog < strcat(AlignDir, sampleName, "_NovosortDedup.log") >;

		// Sort and Mark Duplicates in one step
		dedupSortedBam, novoLog = novosort(vars["NOVOSORTDIR"], alignedBam, vars["TMPDIR"],
					  string2int(vars["PBSCORES"]), ["--markDuplicates"]
					 );
	}
}

/**********
Realignment
***********/

(file realignedbam) realignBam(string prefix, string sampleName, string realparms[], file inputBam) {
	// Prefix is something like: RealignDir, sampleName, ".", chr

	// Log files												    
	file targetLog < strcat(prefix, "_RealignTargetCreator.log") >;				 
	file realignLog < strcat(prefix, "_IndelRealigner.log") >;
	file intervals < strcat(prefix, ".realignTargetCreator.intervals") >;			   
															
	// The inputBam should be indexed before this function is called						
	intervals, targetLog = RealignerTargetCreator(vars["JAVADIR"], vars["GATKDIR"],				 
					   strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),			
					   inputBam, string2int(vars["PBSCORES"]), realparms			    
					  );									    
	realignedbam, realignLog = IndelRealigner(vars["JAVADIR"], vars["GATKDIR"],				     
				      strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),			     
				      inputBam, realparms, intervals						    
				     );										 
	checkBam(realignedbam, sampleName, "realignment");							      
}			   


/************
Recalibration
*************/

(file outBam) recalibrateBam(string prefix, string sampleName, file inputBam, string recalparmsindels[]) {
	// Log files
	file recalLog < strcat(prefix, "_BaseRecalibrator.log") >;
	file printLog < strcat(prefix, "_PrintReads.log") >;
	file report < strcat(prefix, ".recal_report.grp") >;

	// The inputBam should be indexed before this function is called
	report, recalLog = BaseRecalibrator(vars["JAVADIR"], vars["GATKDIR"],
				       strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]), inputBam,
				       string2int(vars["PBSCORES"]), recalparmsindels,
				       strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"])
				      );
	outBam, printLog = PrintReads(vars["JAVADIR"], vars["GATKDIR"],
				      strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]), inputBam,
				      string2int(vars["PBSCORES"]), report
				     );
	checkBam(outBam, sampleName, "recalibration");
}

/*************
VariantCalling (for split chromosome path)
*************/
(file outVCF) callChrVariants(string sampleName, file inputBam, string chr) {

	int ploidy;
	if ( chr == "M" || chr == "chrM" ) {
		ploidy = 1;
	}
	else {
		ploidy = 2;
	}

	string VarcallDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/");
	
	// Log file
	file haploLog < strcat(VarcallDir, sampleName, ".", chr, "_HaplotypeCaller.log") >;

	outVCF, haploLog = HaplotypeCaller(vars["JAVADIR"], vars["GATKDIR"],	     
					   strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),   
					   inputBam,					
					   strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),       
					   string2int(vars["PBSCORES"]), ploidy, chr	       
					  );
}

/*************
VariantCalling (Without splitting chromosomes)
**************/
(file outVCF) callVariants(string sampleName, file inputBam) {

	string VarcallDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/");

	// Log file
	file haploLog < strcat(VarcallDir, sampleName, "_HaplotypeCaller.log") >;

	outVCF, haploLog = HaplotypeCaller(vars["JAVADIR"],
					   vars["GATKDIR"],
					   strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),
					   inputBam,
					   strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),
					   string2int(vars["PBSCORES"])
					  );
}

/********************
BAM File Verification
*********************/

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

/********************
Recalibration Wrapper (If no chromosome splitting, chr input string should be "")
*********************/

(file recalibratedbam) recalibrationWrapper(string prefix, string sampleName, file inputBam, string chr) {
	// Prefix is either RealignDir, sampleName, ".", chr or RealignDir, sampleName"

	/*****
	Gather list of all reference recalibration files,			       
	then retrieve indels and parameters from them				   
	*****/									  
	// Temporary file
	file recalfiles < strcat(vars["TMPDIR"], "/", chr, ".recal_foundfiles.txt") >;
	string recalparmsindels[];
	string realparms[];

	recalfiles = find_files(strcat(vars["REFGENOMEDIR"], "/", vars["INDELDIR"], "/"), strcat("*", chr, ".vcf" ));
	recalparmsindels = split(trim(replace_all(read(sed(recalfiles, "s/^/--knownSites /g")), "\n", " ", 0)), " ");
	realparms = split(trim(replace_all(read(sed(recalfiles, "s/^/-known /g")), "\n", " ", 0)), " ") =>
	rm(recalfiles);	

	if (vars["ANALYSIS"] == "VC_REALIGN") {
		/*****
		Realignment and Recalibration						   
		*****/
		
		// Realign file handles							 
		file realignedbam < strcat(prefix, ".realigned.bam") >;			 

		// Wait for index to get created before moving on			       
		samtools_index(vars["SAMTOOLSDIR"], inputBam) =>		       
		realignedbam = realignBam(prefix, sampleName, realparms, inputBam);
		recalibratedbam = recalibrateBam(prefix, sampleName, realignedbam, recalparmsindels);
	}
	else {
		/*****									  
		Recalibration Only								   
		*****/
		
		// Wait for index to get created before moving on			       
		samtools_index(vars["SAMTOOLSDIR"], inputBam) =>		       
		recalibratedbam = recalibrateBam(prefix, sampleName, inputBam, recalparmsindels);   
	}
}

/****************
Main Sample Loop
*****************/

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
		string VarcallDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/");
		string RealignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/");

		mkdir(AlignDir);
		mkdir(VarcallDir);
		mkdir(RealignDir);

		/*****
		Alignment
		*****/
		file alignedsam < strcat(vars["TMPDIR"], "/align/", sampleName, ".noDedups.sam") >;
		file alignedbam < strcat(AlignDir, sampleName, ".noDedups.bam") >;		

		alignedsam = alignReads(sampleName, read1, read2, rgheader);
		alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, string2int(vars["PBSCORES"]), ["-u"]);

		// Verify alignment was successful
		checkBam(alignedbam, sampleName, "alignment");
															
		/*****
		Deduplication
		*****/
		file dedupsortedbam < strcat(AlignDir, sampleName, ".wDedups.sorted.bam") >;
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
		file flagstats < strcat(AlignDir, sampleName, ".wDedups.sorted.bam", ".flagstats") >;
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
			
			if (vars["SPLIT"] == "YES" || vars["SPLIT"] == "Yes" || vars["SPLIT"] == "yes" ||
			    vars["SPLIT"] == "Y" || vars["SPLIT"] == "y"
			   ) {
			/****************************************************************************
			Split by Chromosome
			****************************************************************************/

				indices = split(vars["CHRNAMES"], ":");							 

				// List of chrgvcf files
				file gvcfChrCollection[];
															
				foreach chr, chrIndex in indices {
					/*****										  
					Separate dedupped and sorted bam by chromosome 
					*****/										  
					file chrdedupsortedbam < strcat(RealignDir, sampleName, ".", chr,	       
									".wDedups.sorted.bam"			   
								       ) >;
					chrdedupsortedbam = samtools_view(vars["SAMTOOLSDIR"], dedupsortedbam,		  
									  string2int(vars["PBSCORES"]), [strcat(chr)]	   
									 ) =>						   

					// Verify chromosome separation worked						  
					checkBam(chrdedupsortedbam, sampleName, "realignment");  
	
					/***************************************************************
					Realignment and/or Recalibration
					****************************************************************/	   
					string prefix = strcat(RealignDir, sampleName, ".", chr);
				
					file recalibratedbam < strcat(prefix, ".recalibrated.bam") >;
					recalibratedbam = recalibrationWrapper(prefix, sampleName, 
									       chrdedupsortedbam, chr
									      );

					/**************
					Variant Calling
					***************/ 

					// Put gvcfvariant handle in the chrCompleted array			     
					file gvcfSampleChr < strcat(VarcallDir, sampleName, ".", chr, ".raw.g.vcf") >;

					gvcfSampleChr = callChrVariants(sampleName, recalibratedbam, chr) =>
					// Add this vcf file to the foreach loop list
					gvcfChrCollection[chrIndex] = gvcfSampleChr;
				} // end the loop for all chromosomes
				
				/*******************************************************
				Combine Variant Calls for each Chromosome
				*******************************************************/

				/*
				Get the names of each of the chromosome split sample files, add a "--variants" string
				before them, and feed the array to CombineGVCF
				*/
				string chrSampleArray[];
				foreach f, ind in gvcfChrCollection {
					/*
					For each sample.vcf added to the chrSampleArray needs "--variant" added in
					front of it. Note: simply concatenating "--variant " with the sample name
					then adding to the array produces an error in the command line call
					*/
					int varFlagPos = ind * 2;
					int namePos = varFlagPos + 1;
				
					chrSampleArray[varFlagPos] = "--variant";
					chrSampleArray[namePos] = filename(f);
				}

				// CombineGVCF output file
				file gvcfSample < strcat(VarcallDir, sampleName, ".raw.g.vcf") >;

				// Log file for HaplotypeCaller							 
				file combineLog < strcat(VarcallDir, sampleName, "_CombineGVCFs.log") >;
	       			/*
				Wait for end of the foreach loop, then regroup all of the chromosome files
				and put them in the samplesOut array
				*/
				gvcfSample, combineLog = CombineGVCFs(vars["JAVADIR"], vars["GATKDIR"], 
							  strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),
							  strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),
							  chrSampleArray
					 		 ) =>
				samplesOut[index] = gvcfSample;
			}
			else {
			/*********************************************
			Call Variants without Splitting by Chromosome 
			**********************************************/

				/*******************************		
				Realignment and/or Recalibration						
				********************************/
				string prefix = strcat(RealignDir, sampleName);

				file recalibratedbam < strcat(prefix, ".recalibrated.bam") >;	
				// The chr entry for the wrapper is an empty string	   
				recalibratedbam = recalibrationWrapper(prefix, sampleName, dedupsortedbam, "");

				/**************
				Variant Calling
				***************/			     
				file gvcfSample < strcat(VarcallDir, sampleName, ".raw.g.vcf") >;  
				gvcfSample = callVariants(sampleName, recalibratedbam) =>	    
				
				// Add to the output array
				samplesOut[index] = gvcfSample;
			}
		} // End of align only check
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
		// Log file for Joint Genotyping
		file jointLog < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/jointVCFs/jointVCF.log") >;
		
		mkdir(strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/jointVCFs"));

		string variantSampleArray[];

		foreach item, index in samplesOutput {
				/*
				For each sample.vcf added to the variantSampleArray needs "--variant" in front of it  

				Note: simply concatenating "--variant " with the sample name then adding to the array   
				produces an error in the command line call
				*/											 
				int varFlagPos = index * 2; 
				int namePos = varFlagPos + 1;
				variantSampleArray[varFlagPos] = "--variant";
				variantSampleArray[namePos] = filename(item);						  
			}

		jointVCF, jointLog = GenotypeGVCFs(vars["JAVADIR"], vars["GATKDIR"], strcat(vars["REFGENOMEDIR"],
					 "/", vars["REFGENOME"]), variantSampleArray
					); 
	}
}


