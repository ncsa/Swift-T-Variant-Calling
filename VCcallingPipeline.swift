// The way to run this script is:
// swift-t VCcallingPipeline.swift --runfile=HgG0.lowcoverage.chr20.parameters-azza

import sys;
import files;
import string;
import unix;
import assert;
import pipelinefunctions.align_dedup;
import pipelinefunctions.realign_varcall_by_chr;
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

////////////////////// Main loop begins //////////////////////////////////////
foreach sample in sampleLines{

	string sampleInfo[] = split(sample, " ");
	string sampleName = sampleInfo[0];
	string read1 = sampleInfo[1];
	string read2 = sampleInfo[2];
	
	string rgheader = sprintf("@RG\tID:%s\tLB:%s\tPL:%s\tPU:%s\tSM:%s\tCN:%s", sampleName, vars["SAMPLELB"], vars["SAMPLEPL"], sampleName, sampleName, vars["SAMPLECN"]);

	///////////////// Alignment and deduplication (per sample):

	//Needed files (same order as in align_dedup.sh script):
	// outputdir=$rootdir/$SampleName= OUTPUTDIR/SAMPLENAME !!
	file qcfile <strcat(vars["OUTPUTDIR"], "/", vars["DELIVERYFOLDER"], "/docs/QC_test_results.txt") >;

	file alignedsam <strcat(vars["TMPDIR"],"/align/", sampleName, ".nodups.sam")>;
	
	file dedupsortedbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".wdups.sorted.bam")>;

	// These are not specifically defined!
	file flagstats <strcat(filename(dedupsortedbam), ".flagstats")>;

	if (vars["MARKDUPLICATESTOOL"] == "SAMBLASTER") {
		trace("##CASE1: dedup tool is ## SAMBLASTER ##. We use a single command for align-deduplication ##");
		///////// Pipe everything up to the dedupbam! /////////
		file dedupsam <strcat(vars["TMPDIR"], "/align/", sampleName, ".wdups.sam")>;
		file dedupbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".wdups.bam")>;
		alignedsam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], vars["BWAMEMPARAMS"], string2int(vars["PBSCORES"]), rgheader);
		dedupsam = samblaster(vars["SAMBLASTERDIR"], alignedsam);
		dedupbam = samtools_view(vars["SAMTOOLSDIR"], dedupsam, string2int(vars["PBSCORES"]),"-u");
		//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
		dedupsortedbam = novosort(strcat(vars["NOVOCRAFTDIR"],"/","novosort"), dedupbam,vars["TMPDIR"], string2int(vars["PBSCORES"]), "--compression 1");
		//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
	
	} else { 
		file alignedbam <strcat(vars["OUTPUTDIR"],"/align/", sampleName, ".nodups.bam")>;
		if  (vars["MARKDUPLICATESTOOL"] == "PICARD") {
			trace("##CASE 2: dedup tool is ## PICARD ##. One command for align, one for sort, one for dedup ##");
			file alignedsortedbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".nodups.sorted.bam")>;
			file metricsfile <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".picard.metrics")>;
			alignedsam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], vars["BWAMEMPARAMS"], string2int(vars["PBSCORES"]), rgheader);
			alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, string2int(vars["PBSCORES"]), "-u");
			//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
			alignedsortedbam = novosort(strcat(vars["NOVOCRAFTDIR"],"/","novosort"), alignedbam, vars["TMPDIR"], string2int(vars["PBSCORES"]), "\"\"");                  
			dedupsortedbam, metricsfile= picard(vars["JAVADIR"], vars["PICARDDIR"], vars["TMPDIR"], alignedsortedbam ) ; 
			//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!

		} else {
			trace("##CASE DEFAULT: dedup tool is ## NOVOSORT ##. We use one command for dup-sort ##");
			if  (vars["ALIGNERTOOL"] == "BWAMEM") {
				alignedsam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], vars["BWAMEMPARAMS"], string2int(vars["PBSCORES"]), rgheader);
				alignedbam = samtools_view(vars["SAMTOOLSDIR"],alignedsam, string2int(vars["PBSCORES"]), "-u");
	                } else {
				alignedsam = novoalign(strcat(vars["NOVOCRAFTDIR"],"/","novoalign"), read1, read2, vars["NOVOALIGNINDEX"], vars["NOVOALIGNPARAMS"], string2int(vars["PBSCORES"]), rgheader);
				alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, string2int(vars["PBSCORES"]),"-u");
			}
			//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
			dedupsortedbam = novosort(strcat(vars["NOVOCRAFTDIR"],"/","novosort"), alignedbam,vars["TMPDIR"], string2int(vars["PBSCORES"]), strcat("--markDuplicates") );	
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

	foreach chr in indices {
		trace("####   Realign-Vcall script for SAMPLE " +sampleName+ " chr=" +chr+ "                             #######");
		////// Map the output files from this stage! (line 79 in realign_var_call_by_chr.sh onwards!)
		// qcfile is the same as in the alignment stage
		// inputbam is really the dedupsortedbam file
		
		//RealignDir=outputdir/realign=rootdir/$SampleName/realign
		//			      =vars["OUTPUTDIR"]/$SampleName/realig
		file chrdedupsortedbam <strcat(vars["OUTPUTDIR"], "/",sampleName, "/realign/",sampleName, ".", chr, ".wdups.sorted.bam")>;
		file realignedbam <strcat(sampleName, ".", chr, ".realigned.bam")>;
		file recalibratedbam <strcat(sampleName, ".", chr, ".recalibrated.bam")>;
		file rawvariant <strcat(sampleName, ".", chr, ".raw.g.vcf")>;
		
		trace("#######   REALIGN-RECALIBRATION BLOCK STARTS HERE   " + sampleName + chr + "        ######");
		if ( chr=="M" ) {
			ploidy = 1;
		} else {
			ploidy = 2;
		}
		chrdedupsortedbam = samtools_view(vars["SAMTOOLSDIR"], dedupsortedbam, string2int(vars["PBSCORES"]), strcat(chr));
		//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
		file chrdedupsortedbamindex <strcat(sample,".",chr,".wdups.sorted.bam")>;
/*		chrdedupsortedbamindex = samtools_index(vars["SAMTOOLSDIR"], chrdedupsortedbam);

		file recalfiles < strcat(vars["TMPDIR"],"/","recalparms2foundfiles.txt")> ;
	//	recalfiles = find_files( strcat(vars["REFGENOMEDIR"], "/", vars["INDELDIR"], "/"), strcat("*", chr, ".vcf" ) );
		string recalparms2 = replace_all(read(input_file(sed(recalfiles, "s/^/ --knownSites /g"))), "\n", " ", 0);
		trace("\n\n" + recalparms2 + "\n\n");
*/		
	} 

}

