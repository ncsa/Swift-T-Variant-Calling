// The way to run this script is:
// swift-t VCcallingPipeline.swift --runfile=HgG0.lowcoverage.chr20.parameters-azza

import sys;
import files;
import string;
import unix;
import assert;
import align_dedup_test;
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

	mkdir(VarcallDir);
	
	file dedupsortedbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".wdups.sorted.bam")>;


	// These are not specifically defined!

	// These are temporary files:
	file alignedsam <strcat(vars["TMPDIR"],"/align/", sampleName, ".nodups.sam")>;

		trace("##CASE1: dedup tool is ## SAMBLASTER ##. We use a single command for align-deduplication ##");
		file dedupsam <strcat(vars["TMPDIR"], "/align/", sampleName, ".wdups.sam")>;
		file dedupbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".wdups.bam")>;	
		///////// Pipe the stages below up to the dedupbam! /////////
		alignedsam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], vars["BWAMEMPARAMS"], string2int(vars["PBSCORES"]), rgheader) =>
		dedupsam = samblaster(vars["SAMBLASTERDIR"], alignedsam) =>
		dedupbam = samtools_view(vars["SAMTOOLSDIR"], dedupsam, string2int(vars["PBSCORES"]),"-u") ;
		//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit! The code would be something similar to:
		//int numAlignments_sam = samtools_view(vars["SAMTOOLSDIR"], dedupbam);
		//if (numAlignments_sam==0) { qcfile = echo(strcat(sampleName, "\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for ", filename(dedupbam), "\n"));	}
		//assert (numAlignments_sam > 0, strcat("bwa mem command did not produce alignments for ", filename(dedupbam), " alignment failed"));
		dedupsortedbam = novosort(strcat(vars["NOVOCRAFTDIR"],"/","novosort"), dedupbam,vars["TMPDIR"], string2int(vars["PBSCORES"])  ) ;
		//////// At this stage, check numAlignments, and report if alignment has failed (in the qcfile), and exit!
	

}
