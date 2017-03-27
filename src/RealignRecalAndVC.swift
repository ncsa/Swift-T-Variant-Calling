/*

*****************************
 Pseudocode of Main Functions
*****************************

Note:
	This stage has two versions of the main function
		One for variant calling when splitting by chromosome
		The other for when variants are called on a per sample basis

	This duplication of main functions is caused because they must return different
	  data structures
	  	If splitting by chromosome, the output will be a 2D-matrix of VCF files
	  	Otherwise, it will just be an array of VCF files

********************************************************************************

******************************
 Pseudocode for No Split Main
******************************

(file VCFList[]) VCNoSplitMain(file inputBams[]) {
	foreach sample in samples {
		if (VC_STAGE variable == "Y") {
			**************************
			*** EXECUTE THIS STAGE ***
			**************************

			- Find the sample name
			- Gather the recalibration index files
			- Realign and/or recalibrate
			- Call variants
			- Add the output vcf file to the output array

		} else {
			**********************************
			*** THIS STAGE IS NOT EXECUTED ***
			**********************************
			if (sample VCF file is found) {
				*** SUCCESS ***
				- Add the output vcf file to the output array
			} else {
				*** FAILURE ***
				- Write an error message to the fail log

			}
		}
	}
}

********************************************************************************

***************************
 Pseudocode for Split Main
***************************

Note: although the input matrix is in the form inputMatrix[chromosome][sample],
  the output matrix is in the form outputMatrix[sample][chromosome]


(file VCF_list[][]) VCSplitMain(file inputBams[][]) {
	foreach chrSet in inputBams {
		- Get chromosome name

		foreach inputBam in chrSet {
			if (VC_STAGE variable == "Y") {
				**************************
				*** EXECUTE THIS STAGE ***
				**************************

				- Get the sample name
				- Gather the recalibration index files
				- Realign and/or recalibrate
				- Call variants
				- Add the output vcf file to the output array

			} else {
				if (the output vcf file is found) {
					*** SUCCESS ***
					- Add the output vcf file to the output array
				} else {
					*** FAILURE ***
					- Write an error message to the fail log
				}
			}
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
	file report < strcat(prefix, "recal_report.grp") >;

	// The inputBam should be indexed before this function is called
	report, recalLog = BaseRecalibrator(var["JAVADIR"], var["GATKDIR"],
				       strcat(var["REFGENOMEDIR"], "/", var["REFGENOME"]), inputBam,
				       threads, recalparmsindels,
				       strcat(var["REFGENOMEDIR"], "/", var["DBSNP"])
				      );
	outBam, printLog = PrintReads(var["JAVADIR"], var["GATKDIR"],
				      strcat(var["REFGENOMEDIR"], "/", var["REFGENOME"]), inputBam,
				      threads, report
				     );
	checkBam(var, outBam);
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

	if (var["ANALYSIS"] == "VC_REALIGN") {
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
(file outVCF) callChrVariants(string vars[string], string sampleName, file inputBam, string chr) {

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
(file outVCF) callVariants(string vars[string], string sampleName, file inputBam) {

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

/***********************
 Main Functions
************************/
(file VCF_list[]) VCNoSplitMain(string vars[string], file inputBams[], file failLog) {
	foreach sample, index in inputBams {
		
		string baseName = basename(sample); 
		string sampleName = substring(baseName, 0, strlen(baseName) - 23);  // Verified

		if (vars["VC_STAGE"] == "Y") {

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
			recalibratedbam = recalibrationWrapper(sampleName, "", vars, inputBam,
							       realparms, recalparmsindels
							      );    

			/**************
			 Call variants
			***************/
			file gvcfVariants < strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/", sampleName,
					   ".wDedups.sorted.recalibrated.g.vcf"
					  ) >;
			gvcfVariants = callVariants(vars, sampleName, recalibratedbam);
			VCF_list[index] = gvcfVariants;
		}
		else {
			/************************
			 Gather the output files
			*************************/
			string vcfLocation = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/", sampleName,
						    ".wDedups.sorted.recalibrated.g.vcf"
						   );


			// If the vcf file does not exist were it should be, write an error to the failLog
			if ( file_exists(vcfLocation) ) {
				VCF_list[index] = input(vcfLocation);
			}
			else {
				append(failLog, strcat("ERROR: ", vcfLocation, " not found. Did you set VC_STAGE to",
						       " 'N' by accident?\n"
						      )
				      );
			}
		}
	}
}

(file VCF_list[][]) VCSplitMain(string vars[string], file inputBams[][], file failLog) {
	foreach chrSet, chrIndex in inputBams {
		// Input files will have names in the form 'prefix.chrA.bam'
		// This will grab the chr name from the first sample in that chromosome list
		string base = basename(chrSet[0]);
		string trimmed = substring(base, 0, strlen(base) - 4);  // gets rid of 'bam' extension
		string pieces[] = split(trimmed, ".");		// Splits the string by '.'
		string chr = pieces[size(pieces) - 1];		// Grabs the last part, which is the chromosome part

		// Removes '.chr' part of the sample's name
		string sampleName = substring(base, 0, strlen(base) - strlen(chr) - 1);

		foreach inputBam, sampleIndex in chrSet {
			if (vars["VC_STAGE"] == "Y") {
 				/*************************************							  
				 Gather the recalibration index files							   
				**************************************/							 
															
				// Temporary file									       
				file recalfiles < strcat(vars["TMPDIR"], "/", sampleName, ".", chr, 
							 ".recal_foundfiles.txt"
							) >; 
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
				file recalibratedbam < strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/",
							      sampleName, ".wDedups.sorted.recalibrated.", chr, ".bam"		      
							     ) >;							       
				recalibratedbam = recalibrationWrapper(sampleName, chr, vars, inputBam,			 
								       realparms, recalparmsindels			      
								      );
			
				/**************										 
				 Call variants										  
				***************/								       
				file gvcfVariants < strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/",
							   sampleName, ".wDedups.sorted.recalibrated.", chr, ".g.vcf"		       
						 	  ) >;								  
				gvcfVariants = callChrVariants(vars, sampleName, recalibratedbam, chr);
				VCF_list[sampleIndex][chrIndex] = gvcfVariants;
			}
			// If this stage of processing was skipped
			else {
				/************************
			 	 Gather the output files
				*************************/
				string vcfFileLocation = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/",               
								sampleName, ".wDedups.sorted.recalibrated.",
								chr, ".g.vcf"
							       );

				// If the vcf file does not exist were it should be, write an error to the failLog
				if (file_exists(vcfFileLocation)) {
					VCF_list[sampleIndex][chrIndex] = input(vcfFileLocation);
				}
				else {
					append(failLog, strcat("ERROR: ", vcfFileLocation, " not found. Did you set ",
							       "VC_STAGE to 'N' by accident?\n"
							      )
					      );
				}
			}
		}
	}
}
