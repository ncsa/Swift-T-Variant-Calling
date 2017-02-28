/**********************************************************
 Summary
***********************************************************
Input:
	- Runfile
Output:
	- vcf file for each of the input bam files

************************
Pseudocode for main loop
************************
foreach sample {													
	if not split by chromosome {										    
		get recal files											 
		recalWrapper											    
		callVariants											    
	}													       
	else {													  
		foreach chromosome {										    
			get recal files for chromosome								  
			recalWrapper										    
			callChrVariants										 
		}												       
	}													      
}

*/

import files;
import string;
import sys;

import bioapps.realign_varcal;
import generalfunctions.general;

/**********
Realignment
***********/

(file realignedbam) realignBam(string sampleName, string chr, string var[string], string realparms[], file inputBam
			      ) {
	string prePrefix = strcat(var["OUTPUTDIR"],"/", sampleName, "/realign/", sampleName, ".wDedups.sorted.", chr);
	string preLogPrefix = strcat(var["OUTPUTDIR"],"/", sampleName, "/logs/", sampleName, ".wDedups.sorted.", chr);	
	
	// Removes extra '.' when there is no split by chromosome
	string prefix = replace(prePrefix, "..", ".", 0);
	string logPrefix = replace(preLogPrefix, "..", ".", 0);

	int threads = string2int(var["PBSCORES"]) %/ string2int(var["PROCPERNODE"]);

	// Log files												    
	file targetLog < strcat(logPrefix, "_RealignTargetCreator.log") >;				 
	file realignLog < strcat(logPrefix, "_IndelRealigner.log") >;

	file intervals < strcat(prefix, ".realignTargetCreator.intervals") >;			   
															
	// The inputBam should be indexed before this function is called						
	intervals, targetLog = RealignerTargetCreator(var["JAVADIR"], var["GATKDIR"],				 
					   strcat(var["REFGENOMEDIR"], "/", var["REFGENOME"]),			
					   inputBam, threads, realparms			    
					  );									    
	realignedbam, realignLog = IndelRealigner(var["JAVADIR"], var["GATKDIR"],				     
				      strcat(var["REFGENOMEDIR"], "/", var["REFGENOME"]),			     
				      inputBam, realparms, intervals						    
				     );										 
	checkBam(var, realignedbam);
}			   

/************
Recalibration
*************/

(file outBam) recalibrateBam(string sampleName, string chr, string var[string],
			     file inputBam, string recalparmsindels[]
			    ) {

	string prePrefix = strcat(var["OUTPUTDIR"],"/", sampleName, "/realign/", sampleName, ".wDedups.sorted.", chr);
	string preLogPrefix = strcat(var["OUTPUTDIR"],"/", sampleName, "/logs/", sampleName, ".wDedups.sorted.", chr);

	// Removes extra '.' when there is no split by chromosome
	string prefix = replace(prePrefix, "..", ".", 0);
	string logPrefix = replace(preLogPrefix, "..", ".", 0);

	int threads = string2int(var["PBSCORES"]) %/ string2int(var["PROCPERNODE"]);

	// Log files
	file recalLog < strcat(logPrefix, "_BaseRecalibrator.log") >;
	file printLog < strcat(logPrefix, "_PrintReads.log") >;
	file report < strcat(prefix, ".recal_report.grp") >;

	// The inputBam should be indexed before this function is called
	report, recalLog = BaseRecalibrator(var["JAVADIR"], var["GATKDIR"],
				       strcat(var["REFGENOMEDIR"], "/", var["REFGENOME"]), inputBam,
				       threads, recalparmsindels,
				       strcat(var["REFGENOMEDIR"], "/", var["DBSNP"])
				      );
	outBam, printLog = PrintReads(vars["JAVADIR"], vars["GATKDIR"],
				      strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]), inputBam,
				      threads, report
				     );
	checkBam(vars, outBam);
}

/**************************************
 Realignment and recalibration wrapper
***************************************/

(file recalibratedbam) recalibrationWrapper(string sampleName, string chr, string var[string],
					    file inputBam, string realparms[], string recalparmsindels[]
					   ) {

	string prePrefix = strcat(var["OUTPUTDIR"], "/", sampleName, "/realign/", sampleName, ".", chr);
	// If no chr, there will be an extra '.'
	string prefix = replace(prePrefix, "..", ".", 0);

	if (vars["ANALYSIS"] == "VC_REALIGN") {
		/*****
		 Realignment and Recalibration
		*****/

		// Realign file handles
		file realignedbam < strcat(prefix, ".realigned.bam") >;

		// Wait for index to get created before moving on						
		samtools_index(var["SAMTOOLSDIR"], inputBam) =>						 
		realignedbam = realignBam(sampleName, chr, var, realparms, inputBam);				    
		recalibratedbam = recalibrateBam(sampleName, chr, var, realignedbam, recalparmsindels);		 
	}												       
	else {												  
		/*****
		 Recalibration Only									     
		*****/

		// Wait for index to get created before moving on						
		samtools_index(var["SAMTOOLSDIR"], inputBam) =>						 
		recalibratedbam = recalibrateBam(sampleName, chr, var, inputBam, recalparmsindels);		      
	}
}

/******************************************
VariantCalling (for split chromosome path)
*******************************************/
(file outVCF) callChrVariants(string sampleName, file inputBam, string chr) {

	int threads = string2int(vars["PBSCORES"]) %/ string2int(vars["PROCPERNODE"]);

	int ploidy;
	if ( chr == "M" || chr == "chrM" ) { ploidy = 1; }
	else { ploidy = 2; }

	string LogDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/logs/");
	
	// Log file
	file haploLog < strcat(LogDir, sampleName, ".", chr, "_HaplotypeCaller.log") >;

	outVCF, haploLog = HaplotypeCaller(vars["JAVADIR"], vars["GATKDIR"],	     
					   strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),   
					   inputBam,					
					   strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),       
					   threads, ploidy, chr	       
					  );
}

/**********************************************
VariantCalling (Without splitting chromosomes)
***********************************************/
(file outVCF) callVariants(string sampleName, file inputBam) {

	int threads = string2int(vars["PBSCORES"]) %/ string2int(vars["PROCPERNODE"]);

	string LogDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/logs/");

	// Log file
	file haploLog < strcat(LogDir, sampleName, "_HaplotypeCaller.log") >;

	outVCF, haploLog = HaplotypeCaller(vars["JAVADIR"],
					   vars["GATKDIR"],
					   strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),
					   inputBam,
					   strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),
					   threads
					  );
}

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

/************
 Main Loop
*************/

foreach sample in sampleLines {
	string sampleName = split(sample, " ")[0];
	
	string splitVar = vars["SPLIT"];
	
	/**************************************
	 If we are NOT splitting by chromosome
	***************************************/
	if (splitVar != "YES" && splitVar != "Yes" && splitVar != "yes" && splitVar != "Y" && splitVar != "y") {
		/*************************************
		 Gather the recalibration index files		  
		**************************************/
		
		// Temporary file
		file recalfiles < strcat(vars["TMPDIR"], "/", sampleName, ".recal_foundfiles.txt") >;
		recalfiles = find_files(strcat(vars["REFGENOMEDIR"], "/", vars["INDELDIR"], "/"),
					strcat("*", ".vcf")
				       ) =>
		// Get the realign parameters									      
		string recalparmsindels[] = split(
			trim(replace_all(read(sed(recalfiles, "s/^/--knownSites /g")), "\n", " ", 0)), " "		 
						 ) =>
		string realparms[] = split(
			trim(replace_all(read(sed(recalfiles, "s/^/-known /g")), "\n", " ", 0)), " "
					  ) =>
		rm(recalfiles);
		
		/***************************
		 Realign and/or recalibrate
		****************************/
		// Find the input file
		file inputBam = input(strcat(vars["OUTPUTDIR"], "/", sampleName, "/align/", sampleName,
				      ".wDedups.sorted.bam")
				     );

		file recalibratedbam < strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/", sampleName,
					      ".wDedups.sorted.recalibrated.bam"
					     ) >;
		recalibratedbam = recalibrationWrapper(sampleName, "", vars, inputBam, realparms, recalparmsindels);
		
		/**************
		 Call variants
		***************/
		file gvcfVariants < strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/", sampleName,
					   ".wDedups.sorted.recalibrated.g.vcf"
					  ) >;
		gvcfVariants = callVariants(sample, recalibratedbam);
	}  
	/***************************
	 If splitting by chromosome
	****************************/
	else {
		// Get the chromosomes
		string indices[] = split(vars["CHRNAMES"], ":");

		foreach chr in indices {
			/*************************************
			 Gather the recalibration index files
			**************************************/

			// Temporary file
			file recalfiles < strcat(vars["TMPDIR"], "/", sampleName, ".", chr, ".recal_foundfiles.txt") >;
			recalfiles = find_files(strcat(vars["REFGENOMEDIR"], "/", vars["INDELDIR"], "/"),
						strcat("*", chr, ".vcf" )
				       ) =>

			// Get the realign parameters
			string recalparmsindels[] = split(
				trim(replace_all(read(sed(recalfiles, "s/^/--knownSites /g")), "\n", " ", 0)), " "
							 ) =>
			string realparms[] = split(
				trim(replace_all(read(sed(recalfiles, "s/^/-known /g")), "\n", " ", 0)), " "
						  ) =>
			rm(recalfiles);

			/***************************
			 Realign and/or recalibrate
			****************************/
			// Find the input file
			file inputBam = input(strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/", sampleName,
					      ".wDedups.sorted.", chr, ".bam")
					     );

			file recalibratedbam < strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/", sampleName,
						      ".wDedups.sorted.recalibrated.", chr, ".bam" 
						     ) >;
			recalibratedbam = recalibrationWrapper(sampleName, chr, vars, inputBam,
							       realparms, recalparmsindels
							      );

			/**************
			 Call variants
			***************/
			file gvcfVariants < strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/", sampleName,
						   ".wDedups.sorted.recalibrated.", chr, ".g.vcf"
					  	  ) >;
			gvcfVariants = callChrVariants(sample, recalibratedbam, chr);			
		}
	}
}
