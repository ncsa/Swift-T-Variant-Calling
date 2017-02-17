/**********************************************************
 Summary
***********************************************************
Input:
	- Runfile
	- Directory where bam file locations are
		- if chromosome splitting, the directory will be piping_files/chr_split_files
		- Without chromosome splitting, this will just be piping_files
Output:
	- vcf file for each of the input bam files
	- File with list of the output vcf files

************************
Pseudocode for main loop
************************

foreach chrFile in Directory {
	if ( chromosome splitting ) {
		gather recalibration files by chr
	}
	else {
		gather all recalibration files	
	}

	foreach bamFile in chrFile {
		if ( variant calling with realignment ) {
			realign
			recalibrate
		}
		else {
			recalibrate
		}

		if ( chromosome splitting ) {
			call variants in chromosome separated way (have chromosome name in file name)
		}
		else {
			call variants on the whole sample at once
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

(file realignedbam) realignBam(string prefix, string variables[string], string realparms[], file inputBam) {
	
	int threads = string2int(vars["PBSCORES"]) %/ string2int(vars["PROCPERNODE"]);

	// Prefix is something like: RealignDir, sampleName, ".", chr
	// Log files												    
	file targetLog < strcat(prefix, "_RealignTargetCreator.log") >;				 
	file realignLog < strcat(prefix, "_IndelRealigner.log") >;
	file intervals < strcat(prefix, ".realignTargetCreator.intervals") >;			   
															
	// The inputBam should be indexed before this function is called						
	intervals, targetLog = RealignerTargetCreator(vars["JAVADIR"], vars["GATKDIR"],				 
					   strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),			
					   inputBam, threads, realparms			    
					  );									    
	realignedbam, realignLog = IndelRealigner(vars["JAVADIR"], vars["GATKDIR"],				     
				      strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]),			     
				      inputBam, realparms, intervals						    
				     );										 
	checkBam(variables, realignedbam);							      
}			   

/************
Recalibration
*************/

(file outBam) recalibrateBam(string prefix, string variables[string], file inputBam, string recalparmsindels[]) {

	int threads = string2int(vars["PBSCORES"]) %/ string2int(vars["PROCPERNODE"]);

	// Log files
	file recalLog < strcat(prefix, "_BaseRecalibrator.log") >;
	file printLog < strcat(prefix, "_PrintReads.log") >;
	file report < strcat(prefix, ".recal_report.grp") >;

	// The inputBam should be indexed before this function is called
	report, recalLog = BaseRecalibrator(vars["JAVADIR"], vars["GATKDIR"],
				       strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]), inputBam,
				       threads, recalparmsindels,
				       strcat(vars["REFGENOMEDIR"], "/", vars["DBSNP"])
				      );
	outBam, printLog = PrintReads(vars["JAVADIR"], vars["GATKDIR"],
				      strcat(vars["REFGENOMEDIR"], "/", vars["REFGENOME"]), inputBam,
				      threads, report
				     );
	checkBam(variables, outBam);
}

/*************
VariantCalling (for split chromosome path)
*************/
(file outVCF) callChrVariants(string sampleName, file inputBam, string chr) {

	int threads = string2int(vars["PBSCORES"]) %/ string2int(vars["PROCPERNODE"]);

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
					   threads, ploidy, chr	       
					  );
}

/*************
VariantCalling (Without splitting chromosomes)
**************/
(file outVCF) callVariants(string sampleName, file inputBam) {

	int threads = string2int(vars["PBSCORES"]) %/ string2int(vars["PROCPERNODE"]);

	string VarcallDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/");

	// Log file
	file haploLog < strcat(VarcallDir, sampleName, "_HaplotypeCaller.log") >;

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

argv_accept("runfile", "input"); // return error if other flags are given
string directory = argv("input");
string configFilename = argv("runfile");

/*
 Get bam file locations										     
*/
string matchPattern = "*bam_list.txt"; // Might need \ before * !!!!!!!!					        
													   
// If there was no chromosome splitting, there will only be one file: "bam_list.txt"				
string chrFileList[] =  file_lines(find_files(directory, matchPattern));

/**************
 Parse runfile
***************/

file configFile = input_file(configFilename);									   
string configFileData[] = file_lines(configFile);

// Gather variables
string vars[string] = getConfigVariables(configFileData);

/*****												 
 Gather list of all reference recalibration files,							   
   then retrieve indels and parameters from them							       
*****/

/*************
 Main Loop
**************/

foreach chrFile in chrFileList {
	// Get bam file locations from the chromosome file
	string bamFileList[] = file_lines(input(chrFile));

	string recalparmsindels[];
	string realparms[];

	/*
	 Define where the vcf location file will be and whether it needs a chr name on it
	   The location of each output vcf file is added to this file as they are generated
	*/

	/*
	 Since we cannot declare a file without knowing its name in advance, and the name of this file depends on
	   whether chromosome splitting occurred (which is determined below), we simply record the name of the
	   file as a string, then initialize the file right before it is used at the end
	*/
	string vcfFileName;

	if (vars["SPLIT"] == "Yes" ||								   
	   vars["SPLIT"] == "Y" ||								     
	   vars["SPLIT"] == "y" ||								     
	   vars["SPLIT"] == "yes" ||
	    vars["SPLIT"] == "YES"								      
	  ) {
		// Get chromosome name
		string splitPath[] = split(chrFile, "/");
		string chr = split(splitPath[size(splitPath) - 1], "_")[0];

		vcfFileName = strcat(directory, "/", chr, "_vcf_list.txt");
		
		// Temporary file
		file recalfiles < strcat(vars["TMPDIR"], "/", chr, ".recal_foundfiles.txt") >;
		recalfiles = find_files(strcat(vars["REFGENOMEDIR"], "/", vars["INDELDIR"], "/"),
					strcat("*", chr, ".vcf" )
				       );
		
		// Get the realign parameters									     
       		recalparmsindels = split(
			trim(replace_all(read(sed(recalfiles, "s/^/--knownSites /g")), "\n", " ", 0)), " "
		);
        	realparms = split(trim(replace_all(read(sed(recalfiles, "s/^/-known /g")), "\n", " ", 0)), " ") =>
        	rm(recalfiles);
	}
	else {
		vcfFileName = strcat(directory, "/vcf_list.txt");

		// Temporary file
		file recalfiles < strcat(vars["TMPDIR"], "/.recal_foundfiles.txt") >;
		recalfiles = find_files(strcat(vars["REFGENOMEDIR"], "/", vars["INDELDIR"], "/"), strcat("*", ".vcf"));
	
		// Get the realign parameters
		recalparmsindels = split(
			trim(replace_all(read(sed(recalfiles, "s/^/--knownSites /g")), "\n", " ", 0)), " "	     
					);
		realparms = split(trim(replace_all(read(sed(recalfiles, "s/^/-known /g")), "\n", " ", 0)), " ") =>
		rm(recalfiles); 
	}

	foreach bamLocation in bamFileList {
		file bamFile = input(bamLocation);

		/*
		 Gather info about file
		*/

		// Gets everything in the filename but the '.bam' extension
		prefix = substring(bamLocation, 0, strlen(bamLocation) - 4);
		
		// Declare recalibration variable
		file recalibratedbam < strcat(prefix, ".recalibrated.bam") > ;		

		if (vars["ANALYSIS"] == "VC_REALIGN") {								 
			/*****
			Realignment and Recalibration								   
	       		*****/

			// Realign file handles									
			file realignedbam < strcat(prefix, ".realigned.bam") >;

			// Wait for index to get created before moving on						 
			samtools_index(vars["SAMTOOLSDIR"], bamFile) =>						  
			realignedbam = realignBam(prefix, vars, realparms, bamFile);				 
			recalibratedbam = recalibrateBam(prefix, vars, realignedbam, recalparmsindels);		 
		}												   
        	else {
			/*****										        
			 Recalibration Only									     
			*****/										        

			// Wait for index to get created before moving on						 
			samtools_index(vars["SAMTOOLSDIR"], bamFile) =>						  
			recalibratedbam = recalibrateBam(prefix, vars, bamFile, recalparmsindels);
		}

		/**************
		Variant Calling
		***************/
		
		// Get the sampleName from the file name						  
		string slashSplit[] = split(prefix, "/");						   
		// Get the third to last item, as this will be the sample name			       
		string sampleName = slashSplit[size(slashSplit) - 3];

		string VarcallDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/"); 

		// Map the output vcf list file we determined the name of earlier
		file vcfLocationFile < vcfFileName > = write("");

		// If split by chromosome
		if (vars["SPLIT"] == "Yes" ||
		    vars["SPLIT"] == "Y" ||
		    vars["SPLIT"] == "y" ||
		    vars["SPLIT"] == "yes" ||
		    vars["SPLIT"] == "YES"
		   ) {
			// Get the chromosome and sampleName from the file name
			// "prefix" will be something like ". . . /realign/Test_sample1.wDedups.sorted.chr1"
			string splitPrefix[] = split(prefix, ".");
			// Get the last item in the split prefix, which is certain to be the chr name
			string chr = splitPrefix[size(splitPrefix) - 1];
			
			// Call variants
			file gvcfSample < strcat(VarcallDir, sampleName, ".", chr, ".raw.g.vcf") >;  
			gvcfSample = callChrVariants(sampleName, recalibratedbam, chr) =>
			// Add output vcf to the list of outputs
			append(vcfLocationFile, strcat(filename(gvcfSample), "\n"));		
		}
		else {
			// Call variants
			file gvcfSample < strcat(VarcallDir, sampleName, ".raw.g.vcf") >;
			gvcfSample = callVariants(sampleName, recalibratedbam) => 
			// Add output vcf to the list of outputs
			append(vcfLocationFile, strcat(filename(gvcfSample), "\n"));
		}
	}
}

