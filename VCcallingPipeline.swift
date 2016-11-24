// The way to run this script is:
// swift-t VCcallingPipeline.swift --runfile=HgG0.lowcoverage.chr20.parameters-azza

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

//read-in the runfile with argument --runfile=<filename>

argv_accept("runfile"); // return error if user supplies other flagged inputs
string configFilename = argv("runfile");

file configFile = input_file(configFilename);
string configFileData[] = file_lines(configFile);

string vars[string] = getConfigVariables(configFileData);

file sampleInfoFile = input_file(vars["SAMPLEINFORMATION"]);
string sampleLines[] = file_lines(sampleInfoFile);
int samples_processing_done ;

////////////////////// Main loop begins //////////////////////////////////////
foreach sample in sampleLines{

	string sampleInfo[] = split(sample, " ");
	string sampleName = sampleInfo[0];
	string read1 = sampleInfo[1];
	string read2 = sampleInfo[2];
	
	string rgheader = sprintf("@RG\tID:%s\tLB:%s\tPL:%s\tPU:%s\tSM:%s\tCN:%s", sampleName, vars["SAMPLELB"], vars["SAMPLEPL"], sampleName, sampleName, vars["SAMPLECN"]);

	///////////////// Alignment and deduplication (per sample):

	string VarcallDir = strcat(vars["OUTPUTDIR"]/sampleName, "/variant/");
	string RealignDir = strcat(vars["OUTPUTDIR"]/sampleName, "/realign/");

	mkdir(VarcallDir);
	
	file dedupsortedbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".wdups.sorted.bam")>;
	file outbam < strcat(RealignDir, sampleName, ".recalibrated.bam")> ;
	file rawvariant < strcat(VarcallDir, sampleName, ".GATKCombineGVCF.raw.vcf") >;

	file qcfile <strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/docs/QC_test_results.txt") >;
	file mergedbam < strcat(vars["OUTPUTDIR"]/vars["DELIVERYFOLDER"]/sampleName, ".recalibrated.bam")> ;
	file mergedvariant < strcat(vars["OUTPUTDIR"]/vars["DELIVERYFOLDER"]/sampleName, ".GATKCombineGVCF.raw.vcf") >;

	// These are not specifically defined!
	file flagstats <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".wdups.sorted.bam", ".flagstats")>;

	// These are temporary files:
	file alignedsam <strcat(vars["TMPDIR"],"/align/", sampleName, ".nodups.sam")>;
	file chr_bamListfile < strcat(vars["TMPDIR"]/sampleName,".chr_bamList.txt") >;
	file chr_vcfListfile < strcat(vars["TMPDIR"]/sampleName,".chr_vcfList.txt") >;

	if (vars["MARKDUPLICATESTOOL"] == "SAMBLASTER") {
		trace("##CASE1: dedup tool is ## SAMBLASTER ##. We use a single command for align-deduplication ##");
		file dedupsam <strcat(vars["TMPDIR"], "/align/", sampleName, ".wdups.sam")>;
		file dedupbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".wdups.bam")>;	
		///////// Pipe the stages below up to the dedupbam! /////////
		alignedsam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], [vars["BWAMEMPARAMS"]], string2int(vars["PBSCORES"]), rgheader) =>
		dedupsam = samblaster(vars["SAMBLASTERDIR"], alignedsam) =>
		dedupbam = samtools_view(vars["SAMTOOLSDIR"], dedupsam, string2int(vars["PBSCORES"]), ["-u"]) =>
		//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit! The code would be something similar to:
		//int numAlignments_sam = samtools_view(vars["SAMTOOLSDIR"], dedupbam);
		//if (numAlignments_sam==0) { qcfile = echo(strcat(sampleName, "\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for ", filename(dedupbam), "\n"));	}
		//assert (numAlignments_sam > 0, strcat("bwa mem command did not produce alignments for ", filename(dedupbam), " alignment failed"));
		dedupsortedbam = novosort(strcat(vars["NOVOCRAFTDIR"],"/","novosort"), dedupbam,vars["TMPDIR"], string2int(vars["PBSCORES"]), ["--compression", "1"]) ;
		//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
	
	} else { 
		file alignedbam <strcat(vars["OUTPUTDIR"],"/align/", sampleName, ".nodups.bam")>;
		if  (vars["MARKDUPLICATESTOOL"] == "PICARD") {
			trace("##CASE 2: dedup tool is ## PICARD ##. One command for align, one for sort, one for dedup ##");
			file alignedsortedbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".nodups.sorted.bam")>;
			file metricsfile <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".picard.metrics")>;
			alignedsam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], [vars["BWAMEMPARAMS"]], string2int(vars["PBSCORES"]), rgheader) =>
			alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, string2int(vars["PBSCORES"]), ["-u"]) =>
			//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
			alignedsortedbam = novosort(strcat(vars["NOVOCRAFTDIR"],"/","novosort"), alignedbam, vars["TMPDIR"], string2int(vars["PBSCORES"]), []) =>
			dedupsortedbam, metricsfile= picard(vars["JAVADIR"], vars["PICARDDIR"], vars["TMPDIR"], alignedsortedbam ) ;
			//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!

		} else {
			trace("##CASE DEFAULT: dedup tool is ## NOVOSORT ##. We use one command for dup-sort ##");
			if  (vars["ALIGNERTOOL"] == "BWAMEM") {
				alignedsam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], [vars["BWAMEMPARAMS"]], string2int(vars["PBSCORES"]), rgheader) =>
				alignedbam = samtools_view(vars["SAMTOOLSDIR"],alignedsam, string2int(vars["PBSCORES"]), ["-u"]);
	                } else {
				alignedsam = novoalign(strcat(vars["NOVOCRAFTDIR"],"/","novoalign"), read1, read2, vars["NOVOALIGNINDEX"], [vars["NOVOALIGNPARAMS"]], string2int(vars["PBSCORES"]), rgheader) =>
				alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, string2int(vars["PBSCORES"]), ["-u"]);
			}
			//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
			wait (alignedbam) {
				dedupsortedbam = novosort(strcat(vars["NOVOCRAFTDIR"],"/","novosort"), alignedbam,vars["TMPDIR"], string2int(vars["PBSCORES"]), ["--markDuplicates"] );	
			}
			//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
		}
	}

	trace("#############     END ALIGNMENT-DEDUPLICATION BLOCK                   ############");
	trace("########   ALIGNMENT QC TEST   FOR SAMPLE", sampleName, "             ###############");
	trace("########   QC rule1: duplication cutoff <= ", vars["DUP_CUTOFF"], "          ###############");
	trace("########   QC rule2: mapped_reads cutoff >= ", vars["MAP_CUTOFF"], "           ###############");

	flagstats = samtools_flagstat(vars["SAMTOOLSDIR"], dedupsortedbam);

	string stat[] = file_lines(flagstats);
	tot_mapped =  split(stat[4], " ")[0];
	tot_reads = split(stat[0], " ")[0];
	tot_dups = split(stat[3], " ")[0];

	perc_dup= string2float(tot_dups)*100/string2float(tot_reads) ;
	perc_mapped= string2float(tot_mapped)*100/string2float(tot_reads) ;

	if ( perc_dup < string2float(vars["DUP_CUTOFF"]) ) {
		trace("#####  " + sampleName + " passed first filter percent_duplicates with value" +perc_dup + ", maximum cutoff is " + vars["DUP_CUTOFF"]);
		if ( perc_mapped > string2float(vars["MAP_CUTOFF"]) ) {
			trace("#####  " + sampleName + " passed second filter map_cutoff with value" +perc_mapped + ", minimum cutoff is " + vars["MAP_CUTOFF"]);		
			qcfile = echo(sampleName + "\tQCTEST\tPASS\n\trule1 ok: percent_duplication=" +perc_dup+ "<? duplication_cutoff=" +vars["DUP_CUTOFF"]+ "\n\trule2 ok: percent_mapped=" +perc_mapped+ ">? mapping_cutoff=" +vars["MAP_CUTOFF"]);
		} else {
			trace("#####  " + sampleName + " DID NOT pass second filter map_cutoff with value" +perc_mapped + ", minimum cutoff is " + vars["MAP_CUTOFF"]);
			qcfile = echo(sampleName + "\tQCTEST\tFAIL\n\trule1 ok: percent_duplication=" +perc_dup+ "<? duplication_cutoff=" +vars["DUP_CUTOFF"]+ "\n\trule2 not ok: percent_mapped=" +perc_mapped+ ">? mapping_cutoff=" +vars["MAP_CUTOFF"]);
		}	
	} else {
	trace("#####  " + sampleName + " DID NOT pass first filter percent_duplicates with value" +perc_dup + ", maximum cutoff is " + vars["DUP_CUTOFF"]);
	qcfile = echo(sampleName + "\tQCTEST\tFAIL\n\trule1 not ok: percent_duplication=" +perc_dup+ "<? duplication_cutoff=" +vars["DUP_CUTOFF"]+ "\n\trule2 not evaluated: percent_mapped=" +perc_mapped+ ">? mapping_cutoff=" +vars["MAP_CUTOFF"]);
	}

	trace("#############   WRAP UP the align_dedup stage   ############");
	////// email findins to redmine and email: 
	////// seems best to be done in tcl, as a typical mail bash is:
	////// echo $MSG | mail -s $SUBJECT $RECIPIENTS
	////// I'm yet to find something similar to pipes here!
	//MSG="ALIGNMENT-DEDUPLICATION for $SampleName finished successfully"
	//echo -e "program=$scriptfile at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"

	indices = split(vars["CHRNAMES"], ":") ;
	int chromosomes_processing_done;

	wait(dedupsortedbam) {
		foreach chr in indices {
			trace("####   Realign-Vcall script for SAMPLE " +sampleName+ " chr=" +chr+ "                             #######");
			////// Map the output files from this stage! (line 79 in realign_var_call_by_chr.sh onwards!)
			file chrdedupsortedbam <strcat(RealignDir, sampleName, ".", chr, ".wdups.sorted.bam")>;
			file realignedbam <strcat(RealignDir, sampleName, ".", chr, ".realigned.bam")>;
			file recalibratedbam <strcat(RealignDir, sampleName, ".", chr, ".recalibrated.bam")>;
			file intervals < strcat(RealignDir, sampleName, ".", chr, ".realignTargetCreator.intervals") >;
			file recalreport < strcat(RealignDir, sampleName, ".", chr, ".recal_report.grp") >;
			file gvcfvariant <strcat(VarcallDir, sampleName, ".", chr, ".raw.g.vcf")>;
			// These are temporary files only:
			file recalfiles < strcat(vars["TMPDIR"]/sampleName, ".", chr, ".recal_foundfiles.txt") >;
			file realfiles < strcat(vars["TMPDIR"]/sampleName, ".", chr, ".real_foundfiles.txt") >;

			trace("#######   REALIGN-RECALIBRATION BLOCK STARTS HERE   " + sampleName + chr + "        ######");
			int ploidy;
			if ( chr=="M" ) {ploidy = 1;} else {ploidy = 2;}
			chrdedupsortedbam = samtools_view(vars["SAMTOOLSDIR"], dedupsortedbam, string2int(vars["PBSCORES"]), [strcat(chr)]);
			//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
			samtools_index(vars["SAMTOOLSDIR"], chrdedupsortedbam);

			recalfiles = find_files( strcat(vars["REFGENOMEDIR"]/vars["INDELDIR"], "/"), strcat("*", chr, ".vcf" ) );
			string recalparmsindels[] = split(trim(replace_all(read(sed(recalfiles, "s/^/--knownSites /g")), "\n", " ", 0)), " ");
	
			realfiles = find_files( strcat(vars["REFGENOMEDIR"]/vars["INDELDIR"], "/"), strcat("*", chr, ".vcf" ) );
			string realparms[] = split(trim(replace_all(read(sed(recalfiles, "s/^/-known /g")), "\n", " ", 0)), " ");

	//		assert( strlen(recalparmsindels)>1 , strcat("no indels were found for ", chr, " in this folder",vars["REFGENOMEDIR"]/vars["INDELDIR"] ));
	//		assert( strlen(realparms)>1 , strcat("no indels were found for ", chr, "in this folder",vars["REFGENOMEDIR"]/vars["INDELDIR"] ));
		
			intervals = RealignerTargetCreator (vars["JAVADIR"], strcat(vars["GATKDIR"]/"GenomeAnalysisTK.jar"), vars["REFGENOMEDIR"]/vars["REFGENOME"], chrdedupsortedbam, string2int(vars["PBSCORES"]), realparms);
			realignedbam = IndelRealigner (vars["JAVADIR"], strcat(vars["GATKDIR"]/"GenomeAnalysisTK.jar"), vars["REFGENOMEDIR"]/vars["REFGENOME"], chrdedupsortedbam, realparms, intervals);
			//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
			recalreport = BaseRecalibrator (vars["JAVADIR"], strcat(vars["GATKDIR"]/"GenomeAnalysisTK.jar"), vars["REFGENOMEDIR"]/vars["REFGENOME"], realignedbam, string2int(vars["PBSCORES"]), recalparmsindels, vars["REFGENOMEDIR"]/vars["DBSNP"]) ;
			recalibratedbam = PrintReads (vars["JAVADIR"], strcat(vars["GATKDIR"]/"GenomeAnalysisTK.jar"), vars["REFGENOMEDIR"]/vars["REFGENOME"], realignedbam, string2int(vars["PBSCORES"]), recalreport);
			//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!

			trace("#############    GATK VARIANT CALLING   FOR SAMPLE " + sampleName + " : " +  chr + "   ###########");
			gvcfvariant = HaplotypeCaller (vars["JAVADIR"], strcat(vars["GATKDIR"]/"GenomeAnalysisTK.jar"), vars["REFGENOMEDIR"]/vars["REFGENOME"], recalibratedbam, vars["REFGENOMEDIR"]/vars["DBSNP"], string2int(vars["PBSCORES"]), ploidy, chr) =>
			if ( size(indices) == size(glob(strcat(VarcallDir, sampleName, ".*.raw.g.vcf"))) ) { chromosomes_processing_done = 1; };
		
		} // end the loop for all chromosomes

	} // end the wait for dedupsortedbam
	
	wait (chromosomes_processing_done ) {
		chr_bamListfile = find_files (RealignDir, strcat(sampleName, ".*.recalibrated.bam") );
		chr_vcfListfile = find_files (VarcallDir, strcat(sampleName, ".*.raw.g.vcf") );
	}

	trace("#######   MERGE BAMS BLOCK STARTS HERE  FOR             " + sampleName + "      ######");

	string chr_bamList[] = split(trim(replace_all(read(chr_bamListfile), "\n", " ", 0)), " ") =>
	outbam = novosort (strcat(vars["NOVOCRAFTDIR"],"/","novosort"), chr_bamList, vars["TMPDIR"], string2int(vars["PBSCORES"]), []) =>
	mergedbam = cp(outbam);
	//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!

	trace("#######   MERGE VCFs BLOCK STARTS HERE  FOR             " + sampleName + "       ######");
	
	string chr_vcfList[] = split(trim(replace_all(read(sed(chr_vcfListfile, "s/^/--variant /g")), "\n", " ", 0)), " ") =>
	rawvariant = CombineGVCFs (vars["JAVADIR"], strcat(vars["GATKDIR"]/"GenomeAnalysisTK.jar"),  vars["REFGENOMEDIR"]/vars["REFGENOME"], vars["REFGENOMEDIR"]/vars["DBSNP"], chr_vcfList ) =>
	mergedvariant = cp (rawvariant) =>
	if ( size(sampleLines) == size(glob(strcat(vars["OUTPUTDIR"]/vars["DELIVERYFOLDER"], "/*.GATKCombineGVCF.raw.vcf"))) ) { samples_processing_done = 1;  };

} /// End the loop throgh all samples


trace("####    Now launching joint_genotyping script for all SAMPLEs: each 200 together        ##########");

file jointVCF < strcat(vars["OUTPUTDIR"]/vars["DELIVERYFOLDER"], "/jointVCFs/", "jointVCFcalled.vcf") >;
file variantFiles < strcat(vars["TMPDIR"],"/variantFiles.txt") >;
mkdir(strcat(vars["OUTPUTDIR"]/vars["DELIVERYFOLDER"], "/jointVCFs"));
wait(samples_processing_done) {
	variantFiles = find_files (vars["OUTPUTDIR"]/vars["DELIVERYFOLDER"], strcat("*.GATKCombineGVCF.raw.vcf")) =>
	string varlist[] = split(trim(replace_all(read(sed(variantFiles, "s/^/--variant /g")), "\n", " ", 0)), " ") =>
	jointVCF = GenotypeGVCFs (vars["JAVADIR"], strcat(vars["GATKDIR"]/"GenomeAnalysisTK.jar"), vars["REFGENOMEDIR"]/vars["REFGENOME"], varlist );
}
