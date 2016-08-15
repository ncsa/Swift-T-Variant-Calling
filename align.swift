// The way to run this script is:
// swift-t align.swift --runfile=HgG0.lowcoverage.chr20.parameters-azza


import sys;
import files;
import string;

/////////////////////////////////// Pipeline functions definition

/////// Alignment functions:
//@dispatch=WORKER
app (file output) bwa (string bwadir, string read1, string read2, string INDEX, string bwamemparams, int PBSCORES,  string rgheader){
	bwadir "mem" "-M" bwamemparams "-t" PBSCORES "-R" rgheader  INDEX read1 read2 @stdout=output;
}

//@dispatch=WORKER
app (file output) novoalign (string novoaligndir, string read1, string read2, string INDEX, string novoalignparams, int PBSCORES) {
	novoaligndir novoalignparams "-c" PBSCORES "-d" INDEX "-f" read1 read2 @stdout=output;
}

//@dispatch=WORKER
app (file output) samblaster(string samblasterdir, string inputFilename){
	samblasterdir "-M" inputFilename @stdout=output;
}

//@dispatch=WORKER
app (file output) samtools_view(string samtoolsdir, file inputFilename, int thr, string u){
	samtoolsdir "view" "-@" thr "-bS" u inputFilename @stdout=output;
}

//@dispatch=WORKER
app (file output) samtools_flagstat(string samtoolsdir, string inputFilename){
	samtoolsdir "flagstat" inputFilename @stdout=output;
}

//@dispatch=WORKER
app () samtools_index(string samtoolsdir, string inputFilename) {
	samtoolsdir "index" inputFilename
}

//@dispatch=WORKER
app (file output) novosort (string novosortdir, file inputFilename, string tmpdir, int thr, string sortoptions){
	novosortdir "--index" "--tmpdir" tmpdir "--threads" thr inputFilename "-o" output sortoptions;
}

//@dispatch=WORKER
app (file outputfile, file metricsfile) picard (string javadir, string picarddir, string tmpdir, file inputFilename ){
        javadir "-jar" picarddir "MarkDuplicates"  "INPUT=" inputFilename "OUTPUT=" outputfile "METRICS_FILE=" metricsfile;

// This is the command line that is working:
// $javadir -Xmx8g -jar $picarddir MarkDuplicates INPUT=/home/azza/swift-project/Results/align/HG00108.lowcoverage.chr20.smallregion.nodups.bam
}

app (file output) touch (string input){
        "touch" input;
}

app (file output) find (string dir, string pattern){
	"/usr/bin/find" dir "-name" pattern @stdout=output;
}

///////////////////////////////// Reading the runfile parameters:
(string data[string]) getConfigVariables(string lines[])
{
	foreach line in lines
	{
		string keyValuePair[] = split(line, "=");
		string name = keyValuePair[0];
		string value = keyValuePair[1];
		data[name] = value;
	}
}

////////////////////////////////// Running the pipeline:
//get Configuration filename with argument --params [filename]

argv_accept("runfile"); // return error if user supplies other flagged inputs
string configFilename = argv("runfile");

file configFile = input_file(configFilename);
string configFileData[] = file_lines(configFile);

string vars[string] = getConfigVariables(configFileData);

file sampleInfoFile = input_file(vars["SAMPLEINFORMATION"]);
string sampleLines[] = file_lines(sampleInfoFile);

foreach sample in sampleLines{
	string sampleInfo[] = split(sample, " ");
	string sampleName = sampleInfo[0];
	string read1 = sampleInfo[1];
	string read2 = sampleInfo[2];

	///////////////// Alignment and deduplication (per sample):

	string rgheader = sprintf("@RG\tID:%s\tLB:%s\tPL:%s\tPU:%s\tSM:%s\tCN:%s\t", sampleName, vars["SAMPLELB"], vars["SAMPLEPL"], sampleName, sampleName, vars["SAMPLECN"]);

	file alignedsam <strcat(vars["TMPDIR"],"/align/", sampleName, ".nodups.sam")>;
	file alignedbam <strcat(vars["OUTPUTDIR"],"/align/", sampleName, ".nodups.bam")>;
	file dedupsam <strcat(vars["TMPDIR"], "/align/", sampleName, ".wdups.sam")>;
	file dedupbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".wdups.bam")>;
	file dedupsortedbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".wdups.sorted.bam")>;
	file alignedsortedbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".nodups.sorted.bam")>;
	file metricsfile <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".picard.metrics")>;
//	file flagstats <strcat(filename(dedupsortedbam), ".flagstats")>;


/*	if (vars["MARKDUPLICATESTOOL"] == "SAMBLASTER") {
			alignedsam = bwa(vars["BWADIR"], read1, read2, vars["BWAINDEX"], vars["BWAMEMPARAMS"], string2int(vars["PBSCORES"]), rgheader);
			dedupsam = samblaster(vars["SAMBLASTERDIR"], filename(alignedsam));
			dedupbam = samtools_view(vars["SAMTOOLSDIR"], dedupsam, string2int(vars["PBSCORES"]),"-u");
			dedupsortedbam = novosort(strcat(vars["NOVOCRAFTDIR"],"/","novosort"), filename(dedupbam),vars["TMPDIR"], string2int(vars["PBSCORES"]), "--compression 1");
	} else { */
	       if (vars["MARKDUPLICATESTOOL"] == "NOVOSORT") {
	       		if  (vars["ALIGNERTOOL"] == "BWAMEM"){
				alignedsam = bwa(vars["BWADIR"], read1, read2, vars["BWAINDEX"], vars["BWAMEMPARAMS"], string2int(vars["PBSCORES"]), rgheader);
				alignedbam = samtools_view(vars["SAMTOOLSDIR"],alignedsam, string2int(vars["PBSCORES"]), "-u");
	                } else {
				alignedsam = novoalign(strcat(vars["NOVOCRAFTDIR"],"/","novoalign"), read1, read2, strcat(vars["REFGENOMEDIR"],"/",vars["NOVOALIGNINDEX"]), vars["NOVOALIGNPARAMS"], string2int(vars["PBSCORES"]));
				alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, string2int(vars["PBSCORES"])," ");
			}
			dedupsortedbam = novosort(strcat(vars["NOVOCRAFTDIR"],"/","novosort"), alignedbam,vars["TMPDIR"], string2int(vars["PBSCORES"]), strcat("--markDuplicates -r ", "\"" ,rgheader, "\"") );
		} else { 
			if  (vars["MARKDUPLICATESTOOL"] == "PICARD") { 

			alignedsam = bwa(vars["BWADIR"], read1, read2, vars["BWAINDEX"], vars["BWAMEMPARAMS"], string2int(vars["PBSCORES"]), rgheader);
			alignedbam = samtools_view(vars["SAMTOOLSDIR"], alignedsam, string2int(vars["PBSCORES"]), "-u");
			alignedsortedbam = novosort(strcat(vars["NOVOCRAFTDIR"],"/","novosort"), alignedbam, vars["TMPDIR"], string2int(vars["PBSCORES"]), "\"\"");                  
			dupsortedbam, metricsfile= picard(vars["JAVADIR"], vars["PICARDDIR"], vars["TMPDIR"], alignedsortedbam ) ; 
		}  
	}
  //    }

	///////////////// Alignment QC (per sample):
	//start from line 568 in align_dedup.sh
//	flagstats = samtools_flagstat(vars["SAMTOOLSDIR"], filename(dedupsortedbam));
	// parse the file to get the 'mapped', 'in total' and 'duplicates'
//	tot_mapped = 
//	tot_reads = 
//	tot_dups =

	

/*	foreach chr in vars["CHRNAMES"] {	
		file chrdedupsortedbam <strcat(sample,".",chr,".wdups.sorted.bam")>;
		// file chrdedupsortedbamindex <strcat(sample,".",chr,".wdups.sorted.bam")>;
		chrdedupsortedbam = samtools_view(vars["SAMTOOLSDIR"], filename(dedupsortedbam), string2int(vars["PBSCORES"]), strcat("-u -h"," ", chr));
				    samtools_index(vars["SAMTOOLSDIR"], filename(chrdedupsortedbam));
		
		recalparms2=$( find ${PWD} -name "${chr}.*.vcf" | sed "s/^/ --knownSites /g" | tr "\n" " " )

		string recalparms2= find(vars["REFGENOMEDIR"]/vars["INDELDIR"], strcat(chr,".*.vcf"))
		strcat("--knownSites",)





	} */

}

