// The way to run this script is:
// swift-t -r $PWD/pipelinefunctions VCcallingPipeline.swift --runfile=HgG0.lowcoverage.chr20.parameters-azza
// Note:
// You may need additional logging, for example, how compilation works (add `-L <stc log file name>` to the command above); or how implementaion is mapped onto servers (you can enable turbine logging by `export TURBINE_LOG=1` before calling the command above)

import sys;
import files;
import string;
import unix;
import assert;
import pipelinefunctions.align_dedup;
import pipelinefunctions.realign_varcall_by_chr;
import pipelinefunctions.merge_vcf;
import pipelinefunctions.joint_vcf;
import pipelinefunctions.miscellaneous;

/////////////////////////// Pipeline begins ////////////////////////////////////

//read-in the runfile with argument --runfile=<runfile_name>

argv_accept("runfile"); // return error if user supplies other flagged inputs
string configFilename = argv("runfile");

file configFile = input_file(configFilename);
string configFileData[] = file_lines(configFile);

string vars[string] = getConfigVariables(configFileData);

file sampleInfoFile = input_file(vars["SAMPLEINFORMATION"]);
string sampleLines[] = file_lines(sampleInfoFile);
int samples_processing_done;

// Copy the runfile and sampleInfoFile to the docs directory for documentation purposes
file docRunfile < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/docs/", basename_string(configFilename)) > = configFile;
file docSampleInfo < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/docs/", basename_string(filename(sampleInfoFile))) > = sampleInfoFile;  

////////////////////// Main loop begins //////////////////////////////////////
foreach sample in sampleLines{

	string sampleInfo[] = split(sample, " ");
	string sampleName = sampleInfo[0];
	string read1 = sampleInfo[1];
	string read2 = sampleInfo[2];
	
	string rgheader = sprintf("@RG\tID:%s\tLB:%s\tPL:%s\tPU:%s\tSM:%s\tCN:%s", sampleName, vars["SAMPLELB"], vars["SAMPLEPL"], sampleName, sampleName, vars["SAMPLECN"]);

	///////////////// Alignment and deduplication (per sample):

	string AlignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/align/");
	string VarcallDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/variant/");
	string RealignDir = strcat(vars["OUTPUTDIR"], "/", sampleName, "/realign/");
 
	mkdir(AlignDir);
	mkdir(VarcallDir);
	mkdir(RealignDir);
	
	file dedupsortedbam < strcat(AlignDir, sampleName, ".wdups.sorted.bam") >;
	file outbam < strcat(RealignDir, sampleName, ".recalibrated.bam") >;
	file rawvariant < strcat(VarcallDir, sampleName, ".GATKCombineGVCF.raw.vcf") >;

	file qcfile <strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/docs/QC_test_results.txt") >;
	file mergedbam < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/", sampleName, ".recalibrated.bam") > ;
	file mergedvariant < strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/", sampleName, ".GATKCombineGVCF.raw.vcf") >;

	// These are not specifically defined!
	file flagstats < strcat(AlignDir, sampleName, ".wdups.sorted.bam", ".flagstats") >;

	// These are temporary files: If piping is implemented, they would not be needed!
	file alignedsam < strcat(vars["TMPDIR"], "/align/", sampleName, ".nodups.sam") >;
	file chr_bamListfile < strcat(vars["TMPDIR"], "/", sampleName, ".chr_bamList.txt") >;
	file chr_vcfListfile < strcat(vars["TMPDIR"], "/", sampleName,".chr_vcfList.txt") >;

	if (vars["MARKDUPLICATESTOOL"] == "SAMBLASTER") {
		file dedupsam < strcat(vars["TMPDIR"], "/align/", sampleName, ".wdups.sam") >;
		file dedupbam < strcat(AlignDir, sampleName, ".wdups.bam") >;	
		///////// Pipe the stages below up to the dedupbam! /////////
		alignedsam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], [vars["BWAMEMPARAMS"]], string2int(vars["PBSCORES"]), rgheader) =>
		dedupsam = samblaster(vars["SAMBLASTERDIR"], alignedsam) =>
		dedupbam = samtools_view(vars["SAMTOOLSDIR"], dedupsam, string2int(vars["PBSCORES"]), ["-u"]) => 
		int numAlignments_dedup = samtools_view2(vars["SAMTOOLSDIR"], filename(dedupbam));
		if (numAlignments_dedup==0) { qcfile = echo(strcat(sampleName, "\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for ", filename(dedupbam), "\n"));	}
		assert (numAlignments_dedup > 0, strcat("bwa mem command did not produce alignments for ", filename(dedupbam), " alignment failed"));
		dedupsortedbam = novosort(vars["NOVOSORTDIR"], dedupbam,vars["TMPDIR"], string2int(vars["PBSCORES"]), ["--compression", "1"]) => 
		int numAlignments_dedupsortedbam = samtools_view2(vars["SAMTOOLSDIR"], filename(dedupsortedbam));
		if (numAlignments_dedupsortedbam==0) { qcfile = echo(strcat(sampleName, "\tALIGNMENT\tFAIL\t novosort command did not produce alignments for ", filename(dedupsortedbam), "\n"));	}
		assert (numAlignments_dedupsortedbam > 0, strcat("novosort command did not produce alignments for ", filename(dedupsortedbam), " Sorting bam failed"));
	} else { 
		file alignedbam <strcat(AlignDir, sampleName, ".nodups.bam")>;
		if  (vars["MARKDUPLICATESTOOL"] == "PICARD") {
			file alignedsortedbam <strcat(AlignDir, sampleName, ".nodups.sorted.bam")>;
			file metricsfile <strcat(AlignDir, sampleName, ".picard.metrics")>;
			alignedsam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], [vars["BWAMEMPARAMS"]], string2int(vars["PBSCORES"]), rgheader) =>
			alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, string2int(vars["PBSCORES"]), ["-u"]) =>
			int numAlignments_aligned = samtools_view2(vars["SAMTOOLSDIR"], filename(alignedbam));
			if (numAlignments_aligned==0) { qcfile = echo(strcat(sampleName, "\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for ", filename(alignedbam), "\n"));	}
			assert (numAlignments_aligned > 0, strcat("bwa mem command did not produce alignments for ", filename(alignedbam), " alignment failed"));
			alignedsortedbam = novosort(vars["NOVOSORTDIR"], alignedbam, vars["TMPDIR"], string2int(vars["PBSCORES"]), []) =>
			dedupsortedbam, metricsfile= picard(vars["JAVADIR"], vars["PICARDDIR"], vars["TMPDIR"], alignedsortedbam ) => 
			int numAlignments_dedupsortedbam = samtools_view2(vars["SAMTOOLSDIR"], filename(dedupsortedbam));
			if (numAlignments_dedupsortedbam==0) { qcfile = echo(strcat(sampleName, "\tALIGNMENT\tFAIL\t picard command did not produce alignments for ", filename(dedupsortedbam), "\n"));	}
			assert (numAlignments_dedupsortedbam > 0, strcat("picard command did not produce alignments for ", filename(dedupsortedbam), " deduplication failed"));

		} else { 			//Default markduplicates tool, novosort
			if  (vars["ALIGNERTOOL"] == "BWAMEM") {
				alignedsam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], [vars["BWAMEMPARAMS"]], string2int(vars["PBSCORES"]), rgheader) =>
				alignedbam = samtools_view(vars["SAMTOOLSDIR"],alignedsam, string2int(vars["PBSCORES"]), ["-u"]);
	                } else {
				alignedsam = novoalign(vars["NOVOALIGNDIR"], read1, read2, vars["NOVOALIGNINDEX"], [vars["NOVOALIGNPARAMS"]], string2int(vars["PBSCORES"]), rgheader) =>
				alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, string2int(vars["PBSCORES"]), ["-u"]);
			}
			int numAlignments_aligned = samtools_view2(vars["SAMTOOLSDIR"], filename(alignedbam));
			if (numAlignments_aligned==0) { qcfile = echo(strcat(sampleName, "\tALIGNMENT\tFAIL\tbwa mem or novoalign command did not produce alignments for ", filename(alignedbam), "\n"));	}
			assert (numAlignments_aligned > 0, strcat("bwa mem command did not produce alignments for ", filename(alignedbam), " alignment failed"));
			wait (alignedbam) {
				dedupsortedbam = novosort(vars["NOVOSORTDIR"], alignedbam,vars["TMPDIR"], string2int(vars["PBSCORES"]), ["--markDuplicates"] ) =>
				int numAlignments_dedupsortedbam = samtools_view2(vars["SAMTOOLSDIR"], filename(dedupsortedbam));
				if (numAlignments_dedupsortedbam==0) { qcfile = echo(strcat(sampleName, "\tALIGNMENT\tFAIL\tnovosort command did not produce alignments for ", filename(dedupsortedbam), "\n"));	}
				assert (numAlignments_dedupsortedbam > 0, strcat("novosort command did not produce alignments for ", filename(dedupsortedbam), "sorting and deduplication failed"));

			} // End wait(alignedbam)
		} //End if-else (dedup tool selection: picard or novosort)
	} //End if-else (dedup tool selection: samblaster or other)

	flagstats = samtools_flagstat(vars["SAMTOOLSDIR"], dedupsortedbam);

	string stat[] = file_lines(flagstats);
	tot_mapped =  split(stat[4], " ")[0];
	tot_reads = split(stat[0], " ")[0];
	tot_dups = split(stat[3], " ")[0];

	perc_dup= string2float(tot_dups)*100/string2float(tot_reads) ;
	perc_mapped= string2float(tot_mapped)*100/string2float(tot_reads) ;

	if ( perc_dup < string2float(vars["DUP_CUTOFF"]) ) {
		if ( perc_mapped > string2float(vars["MAP_CUTOFF"]) ) {
			qcfile = echo(sampleName + "\tQCTEST\tPASS\n\trule1 ok: percent_duplication=" +perc_dup+ "<? duplication_cutoff=" +vars["DUP_CUTOFF"]+ "\n\trule2 ok: percent_mapped=" +perc_mapped+ ">? mapping_cutoff=" +vars["MAP_CUTOFF"]);
		} else {
			qcfile = echo(sampleName + "\tQCTEST\tFAIL\n\trule1 ok: percent_duplication=" +perc_dup+ "<? duplication_cutoff=" +vars["DUP_CUTOFF"]+ "\n\trule2 not ok: percent_mapped=" +perc_mapped+ ">? mapping_cutoff=" +vars["MAP_CUTOFF"]);
		}	
	} else {
//	I need to append to the qcfile if it contains some data from above. How to do so?
//	qcfile = echo(sampleName + "\tQCTEST\tFAIL\n\trule1 not ok: percent_duplication=" +perc_dup+ "<? duplication_cutoff=" +vars["DUP_CUTOFF"]+ "\n\trule2 not evaluated: percent_mapped=" +perc_mapped+ ">? mapping_cutoff=" +vars["MAP_CUTOFF"]);
	}

	////// email findings to redmine and email: 
	////// seems best to be done in tcl, as a typical mail bash is:
	////// echo $MSG | mail -s $SUBJECT $RECIPIENTS
	////// I'm yet to find something similar to pipes here!
	//MSG="ALIGNMENT-DEDUPLICATION for $SampleName finished successfully"
	//echo -e "program=$scriptfile at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"

	indices = split(vars["CHRNAMES"], ":") ;
	int chromosomes_processing_done;

	wait(dedupsortedbam) {
		foreach chr in indices {
			////// Map the output files from this stage! (line 79 in realign_var_call_by_chr.sh onwards!)
			file chrdedupsortedbam <strcat(RealignDir, sampleName, ".", chr, ".wdups.sorted.bam")>;

			file realignedbam <strcat(RealignDir, sampleName, ".", chr, ".realigned.bam")>;
			file recalibratedbam <strcat(RealignDir, sampleName, ".", chr, ".recalibrated.bam")>;
			file intervals < strcat(RealignDir, sampleName, ".", chr, ".realignTargetCreator.intervals") >;
			file recalreport < strcat(RealignDir, sampleName, ".", chr, ".recal_report.grp") >;
			file gvcfvariant <strcat(VarcallDir, sampleName, ".", chr, ".raw.g.vcf")>;
			// These are temporary files only:
			file recalfiles < strcat(vars["TMPDIR"], "/", sampleName, ".", chr, ".recal_foundfiles.txt") >;
			file realfiles < strcat(vars["TMPDIR"], "/", sampleName, ".", chr, ".real_foundfiles.txt") >;

			int ploidy;
			if ( chr=="M" ) {ploidy = 1;} else {ploidy = 2;}
			chrdedupsortedbam = samtools_view(vars["SAMTOOLSDIR"], dedupsortedbam, string2int(vars["PBSCORES"]), [strcat(chr)]) =>
			int numAlignments_chrdedupsortedbam = samtools_view2(vars["SAMTOOLSDIR"], filename(chrdedupsortedbam));
			
      if (numAlignments_chrdedupsortedbam == 0) {
				qcfile = echo(strcat(sampleName,
						     "\tREALIGNMENT\tFAIL\tsamtools command did not produce alignments for ",
						     filename(chrdedupsortedbam), "\n"
						    )
					     );
			}
			
			assert (numAlignments_chrdedupsortedbam > 0,
				strcat("samtools command did not produce alignments for ",
				       filename(chrdedupsortedbam), "splitting by chromosome failed"
				      )
			       );		

			recalfiles = find_files(strcat(vars["REFGENOMEDIR"], "/", vars["INDELDIR"], "/"), 
						strcat("*", chr, ".vcf" )
					       );
			
			string recalparmsindels[] = split(trim(replace_all(read(sed(
				recalfiles, "s/^/--knownSites /g")), "\n", " ", 0)), " "
							 );
	
			realfiles = find_files( strcat(vars["REFGENOMEDIR"], "/", vars["INDELDIR"], "/"),
					       strcat("*", chr, ".vcf")
					      );			
			
			string realparms[] = split(trim(replace_all(read(sed(
				recalfiles, "s/^/-known /g")), "\n", " ", 0)), " "
						  );

	//		assert( strlen(recalparmsindels)>1 , strcat("no indels were found for ", chr, " in this folder",vars["REFGENOMEDIR"]/vars["INDELDIR"] ));
	//		assert( strlen(realparms)>1 , strcat("no indels were found for ", chr, "in this folder",vars["REFGENOMEDIR"]/vars["INDELDIR"] ));
		
			// The void output is necessary so "intervals" will have an "output" to wait for before executing
			void indexed = samtools_index(vars["SAMTOOLSDIR"], chrdedupsortedbam) =>
			intervals = RealignerTargetCreator (vars["JAVADIR"], vars["GATKDIR"], vars["REFGENOMEDIR"]/vars["REFGENOME"], chrdedupsortedbam, string2int(vars["PBSCORES"]), realparms);
			realignedbam = IndelRealigner (vars["JAVADIR"], vars["GATKDIR"], vars["REFGENOMEDIR"]/vars["REFGENOME"], chrdedupsortedbam, realparms, intervals) =>
			int numAlignments_realignedbam = samtools_view2(vars["SAMTOOLSDIR"], filename(realignedbam));
			if (numAlignments_realignedbam==0) { qcfile = echo(strcat(sampleName, "\tREALIGNMENT\tFAIL\tGATK IndelRealigner command did not produce alignments for ", filename(realignedbam), "\n"));	}
			assert (numAlignments_realignedbam > 0, strcat("GATK IndelRealigner command did not produce alignments for ", filename(realignedbam) ));
			recalreport = BaseRecalibrator (vars["JAVADIR"], vars["GATKDIR"], vars["REFGENOMEDIR"]/vars["REFGENOME"], realignedbam, string2int(vars["PBSCORES"]), recalparmsindels, vars["REFGENOMEDIR"]/vars["DBSNP"]) => 
			recalibratedbam = PrintReads (vars["JAVADIR"], vars["GATKDIR"], vars["REFGENOMEDIR"]/vars["REFGENOME"], realignedbam, string2int(vars["PBSCORES"]), recalreport) =>
			int numAlignments_recalibratedbam = samtools_view2(vars["SAMTOOLSDIR"], filename(recalibratedbam));
			if (numAlignments_recalibratedbam==0) { qcfile = echo(strcat(sampleName, "\tRECALIBRATION\tFAIL\tGATK BQSR Recalibration command did not produce alignments for ", filename(recalibratedbam), "\n"));	}
			assert (numAlignments_recalibratedbam > 0, strcat("GATK BQSR Recalibrator command did not produce alignments for ", filename(recalibratedbam) ));
			gvcfvariant = HaplotypeCaller (vars["JAVADIR"], vars["GATKDIR"], vars["REFGENOMEDIR"]/vars["REFGENOME"], recalibratedbam, vars["REFGENOMEDIR"]/vars["DBSNP"], string2int(vars["PBSCORES"]), ploidy, chr) =>
			if ( size(indices) == size(glob(strcat(VarcallDir, sampleName, ".*.raw.g.vcf"))) ) { chromosomes_processing_done = 1; };
		
		} // end the loop for all chromosomes

	} // end the wait for dedupsortedbam
	
	wait (chromosomes_processing_done ) {
		chr_bamListfile = find_files (RealignDir, strcat(sampleName, ".*.recalibrated.bam") );
		chr_vcfListfile = find_files (VarcallDir, strcat(sampleName, ".*.raw.g.vcf") );
	}


	string chr_bamList[] = split(trim(replace_all(read(chr_bamListfile), "\n", " ", 0)), " ") =>
	outbam = novosort (vars["NOVOSORTDIR"], chr_bamList, vars["TMPDIR"], string2int(vars["PBSCORES"]), []) =>
	mergedbam = cp(outbam) =>
	int numAlignments_mergedbam = samtools_view2(vars["SAMTOOLSDIR"], filename(mergedbam));
	if (numAlignments_mergedbam==0) { qcfile = echo(strcat(sampleName, "\tMERGE\tFAIL\tnovosort command did not produce alignments for ", filename(mergedbam), "\n"));	}
	assert (numAlignments_mergedbam > 0, strcat("novosort command did not produce alignments for ", filename(mergedbam) ));
	
	string chr_vcfList[] = split(trim(replace_all(read(sed(chr_vcfListfile, "s/^/--variant /g")), "\n", " ", 0)), " ") =>
	rawvariant = CombineGVCFs (vars["JAVADIR"], vars["GATKDIR"],  vars["REFGENOMEDIR"]/vars["REFGENOME"], vars["REFGENOMEDIR"]/vars["DBSNP"], chr_vcfList ) =>
	mergedvariant = cp (rawvariant) =>
	if ( size(sampleLines) == size(glob(strcat(vars["OUTPUTDIR"]/vars["DELIVERYFOLDER"], "/*.GATKCombineGVCF.raw.vcf"))) ) { samples_processing_done = 1;  };

} /// End the loop throgh all samples



file jointVCF < strcat(vars["OUTPUTDIR"]/vars["DELIVERYFOLDER"], "/jointVCFs/", "jointVCFcalled.vcf") >;
file variantFiles < strcat(vars["TMPDIR"],"/variantFiles.txt") >;
mkdir(strcat(vars["OUTPUTDIR"]/vars["DELIVERYFOLDER"], "/jointVCFs"));
wait(samples_processing_done) {
	variantFiles = find_files (vars["OUTPUTDIR"]/vars["DELIVERYFOLDER"], strcat("*.GATKCombineGVCF.raw.vcf")) =>
	string varlist[] = split(trim(replace_all(read(sed(variantFiles, "s/^/--variant /g")), "\n", " ", 0)), " ") =>
	jointVCF = GenotypeGVCFs (vars["JAVADIR"], vars["GATKDIR"], vars["REFGENOMEDIR"]/vars["REFGENOME"], varlist );
}
