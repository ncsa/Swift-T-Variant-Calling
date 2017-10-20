/*

*****************************
 Pseudocode of Run Functions
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
 Pseudocode for No Split Run
******************************

(file VCFList[]) VCNoSplitRun(file inputBams[]) {
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
 Pseudocode for Split Run
***************************

Note: although the input matrix is in the form inputMatrix[chromosome][sample],
  the output matrix is in the form outputMatrix[sample][chromosome]


(file VCF_list[][]) VCSplitRun(file inputBams[][]) {
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
import io;

import bioapps.realign_varcal;
import bioappsLoggingFunctions.realign_varcal_logging;
import generalfunctions.general;

/**********
Realignment
***********/

(file realignedbam) realignBam(string sampleName, string chr, string var[string], string realparms[], file inputBam) {

	exec_check(var["GATKJAR"], "GATKJAR"); //Check the proper jar is specified in the runfile
	exec_check(var["JAVAEXE"], "JAVAEXE");

	string prePrefix = strcat(var["OUTPUTDIR"],"/", sampleName, "/realign/", sampleName, ".wDedups.sorted.", chr);
	string preLogPrefix = strcat(var["OUTPUTDIR"],"/", sampleName, "/logs/", sampleName, ".wDedups.sorted.", chr);	
	
	// Removes extra '.' when there is no split by chromosome
	string prefix = replace(prePrefix, "..", ".", 0);
	string logPrefix = replace(preLogPrefix, "..", ".", 0);
	
	int threads = string2int(var["CORES_PER_NODE"]) %/ string2int(var["PROGRAMS_PER_NODE"]);

	// Log files												    
	file targetLog < strcat(logPrefix, "_RealignTargetCreator.log") >;				 
	file realignLog < strcat(logPrefix, "_IndelRealigner.log") >;

	string tmpLogDir = strcat(var["TMPDIR"], "/timinglogs/" );
	file tmptargetLog < strcat(tmpLogDir, sampleName, ".", chr, "_RealignTargetCreator.log")>;
	file tmprealignLog < strcat(tmpLogDir, sampleName, ".", chr, "_IndelRealigner.log")>;


	file intervals < strcat(prefix, ".realignTargetCreator.intervals") >;			   
															
	// The inputBam should be indexed before this function is called						
	intervals, targetLog, tmptargetLog = RealignerTargetCreator_logged (var["JAVAEXE"], var["JAVA_MAX_HEAP_SIZE"], var["GATKJAR"],				 
					   strcat(var["REFGENOMEDIR"], "/", var["REFGENOME"]),
					   inputBam, threads, realparms, sampleName, chr 
					  );									    
	realignedbam, realignLog, tmprealignLog = IndelRealigner_logged (var["JAVAEXE"], var["JAVA_MAX_HEAP_SIZE"], var["GATKJAR"],				     
				      strcat(var["REFGENOMEDIR"], "/", var["REFGENOME"]),			     
				      inputBam, realparms, intervals, sampleName, chr 
				     ) ; 

	checkBam(var, realignedbam);
}			   

/************
Recalibration
*************/

(file outBam) recalibrateBam(string sampleName, string chr, string var[string],
			     file inputBam, string recalparmsindels[]
			    ) {
	
	exec_check(var["GATKJAR"], "GATKJAR"); // Check the proper jar file is specified in the runfile
	exec_check(var["JAVAEXE"], "JAVAEXE");

	string prePrefix = strcat(var["OUTPUTDIR"],"/", sampleName, "/realign/", sampleName, ".wDedups.sorted.", chr);
	string preLogPrefix = strcat(var["OUTPUTDIR"],"/", sampleName, "/logs/", sampleName, ".wDedups.sorted.", chr);

	// Removes extra '.' when there is no split by chromosome
	string prefix = replace(prePrefix, "..", ".", 0);
	string logPrefix = replace(preLogPrefix, "..", ".", 0);

	int threads = string2int(var["CORES_PER_NODE"]) %/ string2int(var["PROGRAMS_PER_NODE"]);

	// Log files
	file recalLog < strcat(logPrefix, "_BaseRecalibrator.log") >;
	file printLog < strcat(logPrefix, "_PrintReads.log") >;
	file report < strcat(prefix, ".recal_report.grp") >;

	string tmpLogDir = strcat(var["TMPDIR"], "/timinglogs/" );
	file tmprecalLog < strcat(tmpLogDir, sampleName, ".", chr, "_BaseRecalibrator.log")>;
	file tmpprintLog < strcat(tmpLogDir, sampleName, ".", chr, "_PrintReads.log")>;
  
	// The inputBam should be indexed before this function is called
	report, recalLog, tmprecalLog = BaseRecalibrator_logged (var["JAVAEXE"], var["JAVA_MAX_HEAP_SIZE"], var["GATKJAR"],
				       strcat(var["REFGENOMEDIR"], "/", var["REFGENOME"]), inputBam,
				       threads, recalparmsindels,
				       strcat(var["REFGENOMEDIR"], "/", var["DBSNP"]), sampleName, chr
				      );
	outBam, printLog, tmpprintLog = PrintReads_logged (var["JAVAEXE"], var["JAVA_MAX_HEAP_SIZE"], var["GATKJAR"],
				      strcat(var["REFGENOMEDIR"], "/", var["REFGENOME"]), inputBam,
				      threads, report, sampleName, chr
				     ); 

	checkBam(var, outBam);
}

/**************************************
 Realignment and recalibration wrapper
***************************************/

(file recalibratedbam) recalibrationWrapper(string sampleName, string chr, string var[string],
					    file inputBam, string realparms[], string recalparmsindels[]
					   ) {

	string prePrefix = strcat(var["OUTPUTDIR"], "/", sampleName, "/realign/", sampleName, ".wDedups.sorted.", chr);
	// If no chr, there will be an extra '.'
	string prefix = replace(prePrefix, "..", ".", 0);

	if (var["REALIGN"] == "YES" ||
	    var["REALIGN"] == "Yes" ||
	    var["REALIGN"] == "yes" ||
	    var["REALIGN"] == "Y" ||
	    var["REALIGN"] == "y" ||
	    var["REALIGN"] == "TRUE" ||
	    var["REALIGN"] == "True" ||
	    var["REALIGN"] == "true" ||
	    var["REALIGN"] == "T" ||
	    var["REALIGN"] == "t"
	   ) {
		/*****
		 Realignment and Recalibration
		*****/

		// Realign file handles
		file realignedbam < strcat(prefix, ".realigned.bam") >;
		trace("recalibrationWrapper call\t", sampleName, "\t", chr);

		// Wait for index to get created before moving on
		samtools_index(var["SAMTOOLSEXE"], inputBam) =>
		realignedbam = realignBam(sampleName, chr, var, realparms, inputBam) =>
		recalibratedbam = recalibrateBam(sampleName, chr, var, realignedbam, recalparmsindels);
	}
	else {
		/*****
		 Recalibration Only
		*****/

		// Wait for index to get created before moving on
		samtools_index(var["SAMTOOLSEXE"], inputBam) =>
		recalibratedbam = recalibrateBam(sampleName, chr, var, inputBam, recalparmsindels);
	}
}

/******************************************
VariantCalling (for split chromosome path)
*******************************************/
(file outVCF) callChrVariants(string vars[string], string sampleName, file inputBam, string chr) {

	int threads = string2int(vars["CORES_PER_NODE"]) %/ string2int(vars["PROGRAMS_PER_NODE"]);

	int ploidy;
	if ( chr == "M" || chr == "chrM" ) { ploidy = 1; }
	else { ploidy = 2; }

	string LogDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/logs/");
	
	// Log file
	file haploLog < strcat(LogDir, sampleName, ".", chr, "_HaplotypeCaller.log") >;

	string tmpLogDir = strcat(vars["TMPDIR"], "/timinglogs/" );
	file tmphaploLog < strcat(tmpLogDir, sampleName, ".", chr, "_HaplotypeCaller.log")>;

	outVCF, haploLog, tmphaploLog = HaplotypeCaller_logged (vars["JAVAEXE"], var["JAVA_MAX_HEAP_SIZE"], vars["GATKJAR"],	     
					   strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),   
					   inputBam,					
					   strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),       
					   threads, ploidy, chr, sampleName	       
					  );
}

/**********************************************
 VariantCalling (Without splitting chromosomes)
***********************************************/
(file outVCF) callVariants(string vars[string], string sampleName, file inputBam) {

	int threads = string2int(vars["CORES_PER_NODE"]) %/ string2int(vars["PROGRAMS_PER_NODE"]);

	string LogDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/logs/");

	// Log file
	file haploLog < strcat(LogDir, sampleName, "_HaplotypeCaller.log") >;

	string tmpLogDir = strcat(vars["TMPDIR"], "/timinglogs/" );
	file tmphaploLog < strcat(tmpLogDir, sampleName, "_HaplotypeCaller.log")>;

	outVCF, haploLog, tmphaploLog = HaplotypeCaller_logged (vars["JAVAEXE"], vars["JAVA_MAX_HEAP_SIZE"], 
					   vars["GATKJAR"],
					   strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),
					   inputBam,
					   strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"]),
					   threads, sampleName
					  );
}

/***********************
 Run Functions
************************/
(file VCF_list[]) VCNoSplitRun(string vars[string], file inputBams[], file failLog) {
	foreach sample, index in inputBams {
		
		string baseName = basename_string(filename(sample)); 
		string sampleName = substring(baseName, 0, strlen(baseName) - 19);  // Verified on Aug 14th, 2017 

		if (vars["VC_STAGE"] == "Y" ||
		    vars["VC_STAGE"] == "y" ||
		    vars["VC_STAGE"] == "YES" ||
		    vars["VC_STAGE"] == "yes" ||
		    vars["VC_STAGE"] == "Yes" ||
		    vars["VC_STAGE"] == "End" ||
		    vars["VC_STAGE"] == "end" ||
		    vars["VC_STAGE"] == "E" ||
		    vars["VC_STAGE"] == "e" ||
		    vars["VC_STAGE"] == "END"
		   ) {

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
						  );

			//rm(recalfiles);

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
							      ) =>

			/**************
			 Call variants
			***************/
			file gvcfVariants < strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/", sampleName,
					   ".wDedups.sorted.recalibrated.g.vcf"
					  ) >;
			gvcfVariants = callVariants(vars, sampleName, recalibratedbam) =>
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
				string message = strcat("ERROR: ", vcfLocation, " not found. Did you set VC_STAGE to",
							" 'N' by accident?\n"
						       );
				append(failLog, message) =>
				exitIfFlagGiven(vars, message);
			}
		}
	}
}

(file VCF_list[][]) VCSplitRun(string vars[string], file inputBams[][], file failLog) {
	foreach chrSet, chrIndex in inputBams {
		foreach inputBam, sampleIndex in chrSet {
			if (vars["VC_STAGE"] == "Y" ||
			    vars["VC_STAGE"] == "y" ||
			    vars["VC_STAGE"] == "YES" ||
			    vars["VC_STAGE"] == "yes" ||
			    vars["VC_STAGE"] == "Yes" ||
			    vars["VC_STAGE"] == "End" ||
			    vars["VC_STAGE"] == "end" ||
			    vars["VC_STAGE"] == "E" ||
			    vars["VC_STAGE"] == "e" ||
			    vars["VC_STAGE"] == "END"
			   ) {
				// Input files will have names in the form 'prefix.chrA.bam'
        // This will grab the chr name from the first sample in that chromosome list
       	string base = basename_string(filename(chrSet[sampleIndex]));
        string trimmed = substring(base, 0, strlen(base) - 4);  // gets rid of 'bam' extension
				string pieces[] = split(trimmed, ".");          // Splits the string by '.'
				string chr = pieces[size(pieces) - 1];// Grabs the last part: the chromosome name

				string nameBase = basename_string(filename(inputBam));
				// Removes the .wDedups.sorted.<chr>.bam portion of the trimmed string,
				//   leaving only the sample name
				string sampleName = substring(nameBase, 0, strlen(nameBase) - strlen(chr) - 1 - 19);

        			trace("VCNoSplitRun\t", sampleName, "\t", chr);                         
        
 				/*************************************
				 Gather the recalibration index files
				**************************************/

				// Temporary file
				file recalfiles < strcat(vars["TMPDIR"], "/", sampleName, ".", chr,
							 ".recal_foundfiles.txt"
							) >;
				recalfiles = find_files(strcat(vars["REFGENOMEDIR"], "/", vars["INDELDIR"], "/"), 
							strcat(chr, ".vcf" )// changed from strcat("*",chr, ".vcf")
					       ) =>
				// Get the realign parameters
				string recalparmsindels[] = split(							      
					trim(replace_all(read(sed(recalfiles, "s/^/--knownSites /g")), "\n", " ", 0)), " "      
								 ) =>							   
				string realparms[] = split(								     
					trim(replace_all(read(sed(recalfiles, "s/^/-known /g")), "\n", " ", 0)), " "	    
							  );								  
				//rm(recalfiles);

				/***************************
				 Realign and/or recalibrate
				****************************/
				file recalibratedbam < strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/",
							      sampleName, ".wDedups.sorted.recalibrated.", chr, ".bam" 
							     ) >;
				recalibratedbam = recalibrationWrapper(sampleName, chr, vars, inputBam,
								       realparms, recalparmsindels
								      ) =>
			
				/**************										 
				 Call variants										  
				***************/								       
				file gvcfVariants < strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/",
							   sampleName, ".wDedups.sorted.recalibrated.", chr, ".g.vcf"
						 	  ) >;
				gvcfVariants = callChrVariants(vars, sampleName, recalibratedbam, chr) =>
				VCF_list[sampleIndex][chrIndex] = gvcfVariants;
			}
			// If this stage of processing was skipped
			else {

				// Input files will have names in the form 'prefix.chrA.bam'                               
        // This will grab the chr name from the first sample in that chromosome list               
        string base = basename_string(filename(chrSet[sampleIndex]));                              
        string trimmed = substring(base, 0, strlen(base) - 4);  // gets rid of 'bam' extension  
        string pieces[] = split(trimmed, ".");          // Splits the string by '.'                
        string chr = pieces[size(pieces) - 1];// Grabs the last part: the chromosome name

				/************************
			 	 Gather the output files
				*************************/

				string nameBase = basename_string(filename(inputBam));                                     
                                // Removes the .wDedups.sorted.<chr>.bam portion of the trimmed string,                    
                                //   leaving only the sample name                                                          
                                string sampleName = substring(nameBase, 0, strlen(nameBase) - strlen(chr) - 1 - 19);

				string vcfFileLocation = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/",
								sampleName, ".wDedups.sorted.recalibrated.",
								chr, ".g.vcf"
							       );

				// If the vcf file does not exist were it should be, write an error to the failLog
				if (file_exists(vcfFileLocation)) {
					VCF_list[sampleIndex][chrIndex] = input(vcfFileLocation);
				}
				else {
					string m = strcat("ERROR: ", vcfFileLocation, " not found. Did you set ",
							  "VC_STAGE to 'N' by accident?\n"
							 );
					append(failLog, m) =>
					exitIfFlagGiven(vars, m);
				}
			}
		}
	}
}
