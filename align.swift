import sys;
import files;
import string;
import io;
/////////////////////////////////// Pipeline functions definition

/////// Alignment functions:
app (file output) bwa (string read1, string read2, string INDEX, string bwamemparams, int PBSCORES,  string rgheader){
	"/usr/bin/bwa" "mem" "-M" bwamemparams "-t" PBSCORES "-R" rgheader  INDEX read1 read2 @stdout=output;
}

app (file output) novoalign (string read1, string read2, string INDEX, string novoalignparams, int PBSCORES,  string rgheader) {
	"/usr/local/src/novocraft/novoalign" novoalignparams "-c" PBSCORES "-d" INDEX "-f" read1 read2 @stdout=output;
}

app (file output) samblaster(string inputFilename){
	"samblaster" "-M" inputFilename @stdout=output;
}

app (file output) samtools(string inputFilename, int thr, string u){
	"samtools" "view" "-@" thr "-bS" u inputFilename @stdout=output;
}

app (file output) novosort (string inputFilename, string tmpdir, int thr, string sortoptions){
	"novosort" "--index" "--tmpdir" tmpdir "--threads" thr sortoptions inputFilename "-o" output;
}


app () picard (string java, string picard, string tmpdir, string inputFilename, string outputfile, string metricsfile){
//      java "-Xmx8g" "-Djava.io.tmpdir" tmpdir "-jar" picard "MarkDuplicates"  "INPUT" inputFilename "OUTPUT" output "TMP_DIR" tmpdir "ASSUME_
        java "-jar" picard "MarkDuplicates"  inputFilename outputfile metricsfile ;

// This is the command line that is working:
// $javadir -Xmx8g -jar $picarddir MarkDuplicates INPUT=/home/azza/swift-project/Results/align/HG00108.lowcoverage.chr20.smallregion.nodups.bam
}

app (file output) touch (string input){
        "touch" input;
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

	string rgheader = sprintf("@RG\tID:%s\tLB:%s\tPL:%s\tPU:%s\tSM:%s\tCN:%s\t", sampleName, vars["SAMPLELB"], vars["SAMPLEPL"], sampleName, sampleName, vars["SAMPLECN"]);

	file alignedsam <strcat(vars["TMPDIR"],"/align/", sampleName, ".nodups.sam")>;
	file alignedbam <strcat(vars["OUTPUTDIR"],"/align/", sampleName, ".nodups.bam")>;
	file dedupsam <strcat(vars["TMPDIR"], "/align/", sampleName, ".wdups.sam")>;
	file dedupbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".wdups.bam")>;
	file dedupsortedbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".wdups.sorted.bam")>;
	file alignedsortedbam <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".nodups.sorted.bam")>;


	//file sortedbam<SingleFileMapper; file=strcat(parameters.OUTPUTDIR,"/align/", sampleName, ".nodups.sorted.bam")>;


		switch(vars["MARKDUPLICATESTOOL"]){
		case "SAMBLASTER":
			alignedsam = bwa(read1, read2, vars["BWAINDEX"], vars["BWAMEMPARAMS"], string2int(vars["PBSCORES"]), rgheader);
			dedupsam = samblaster(filename(alignedsam));
			dedupbam = samtools(filename(dedupsam), string2int(vars["PBSCORES"]),"-u");
			dedupsortedbam = novosort(filename(dedupbam),vars["TMPDIR"], string2int(vars["PBSCORES"]), "--compression 1");

		case "NOVOSORT":
			switch (vars["ALIGNERTOOL"])
			{
			case "BWAMEM":
				alignedsam = bwa(read1, read2, vars["BWAINDEX"], vars["BWAMEMPARAMS"], string2int(vars["PBSCORES"]), rgheader);
				alignedbam = samtools(filename(alignedsam), string2int(vars["PBSCORES"]), "-u");
	                case "NOVOALIGN":
				alignedsam = novoalign(read1, read2, strcat(vars["REFGENOMEDIR"],"/",vars["NOVOALIGNINDEX"]), vars["NOVOALIGNPARAMS"], string2int(vars["PBSCORES"]));
				alignedbam = samtools(filename(alignedsam), string2int(vars["PBSCORES"])," ");
			}
			dedupsortedbam = novosort(filename(alignedbam),vars["TMPDIR"], string2int(vars["PBSCORES"]), strcat("--markDuplicates -r",rgheader) );

		case "PICARD":
			file dedupsortedbam =  touch(filename(dedupsortedbam));
                        file metrics <strcat(vars["OUTPUTDIR"], "/align/", sampleName, ".picard.metrics")>;
                        metricsfile = touch(filename(metrics));	

			alignedsam = bwa(read1, read2, vars["BWAINDEX"], vars["BWAMEMPARAMS"], string2int(vars["PBSCORES"]), rgheader);
			alignedbam = samtools(filename(alignedsam), string2int(vars["PBSCORES"]), "-u");
			alignedsortedbam = novosort(filename(alignedbam),vars["TMPDIR"], string2int(vars["PBSCORES"])," ");                  
			picard(vars["JAVADIR"], vars["PICARDIR"], vars["TMPDIR"], strcat("INPUT=",filename(alignedsortedbam)), strcat("OUTPUT=",filename(dedupsortedbam)), strcat("METRICS_FILE=",filename(metricsfile))   ) ;
 
	}

start from line 568 in align_dedup.sh

}





