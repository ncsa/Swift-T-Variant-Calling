// The way to run this script is:
// swift-t -r $PWD/pipelinefunctions VCcallingPipeline.swift --runfile=HgG0.lowcoverage.chr20.parameters-azza

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

	// These are temporary files: If piping is implemented, they would not be needed!
	file alignedsam <strcat(vars["TMPDIR"],"/align/", sampleName, ".nodups.sam")>;
	file chr_bamListfile < strcat(vars["TMPDIR"]/sampleName,".chr_bamList.txt") >;
	file chr_vcfListfile < strcat(vars["TMPDIR"]/sampleName,".chr_vcfList.txt") >;

	file alignedbam <strcat(vars["OUTPUTDIR"],"/align/", sampleName, ".nodups.bam")>;
			file alignedsortedbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".nodups.sorted.bam")>;
			file metricsfile <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".picard.metrics")>;
			alignedsam = bwa_mem(vars["BWADIR"], read1, read2, vars["BWAINDEX"], [vars["BWAMEMPARAMS"]], string2int(vars["PBSCORES"]), rgheader) =>
			alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, string2int(vars["PBSCORES"]), ["-u"]) =>
//			int numAlignments_aligned = samtools_view2(vars["SAMTOOLSDIR"], filename(alignedbam));
//			if (numAlignments_aligned==0) { qcfile = echo(strcat(sampleName, "\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for ", filename(alignedbam), "\n"));	}
//			assert (numAlignments_aligned > 0, strcat("bwa mem command did not produce alignments for ", filename(alignedbam), " alignment failed"));
			alignedsortedbam = novosort(strcat(vars["NOVOCRAFTDIR"],"/","novosort"), alignedbam, vars["TMPDIR"], string2int(vars["PBSCORES"]), []) =>
			dedupsortedbam, metricsfile= picard(vars["JAVADIR"], vars["PICARDDIR"], vars["TMPDIR"], alignedsortedbam ) => 
			int numAlignments_dedupsortedbam = samtools_view2(vars["SAMTOOLSDIR"], filename(dedupsortedbam));
//			if (numAlignments_dedupsortedbam==0) { qcfile = echo(strcat(sampleName, "\tALIGNMENT\tFAIL\t picard command did not produce alignments for ", filename(dedupsortedbam), "\n"));	}
//			assert (numAlignments_dedupsortedbam > 0, strcat("picard command did not produce alignments for ", filename(dedupsortedbam), " deduplication failed"));

} /// End the loop throgh all samples


