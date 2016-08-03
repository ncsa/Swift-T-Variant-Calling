import sys;
import files;
import string;
import io;
/////////////////////////////////// Pipeline functions definition

/////// Alignment functions:
app (file output) bwa (string read1, string read2, string INDEX, string bwamemparams, int PBSCORES,  string rgheader){
	"/usr/bin/bwa" "mem" bwamemparams "-t" PBSCORES "-R" rgheader  INDEX read1 read2 @stdout=output;
}

////////////////////////////////
app (file output) samtools(string inputFilename, int thr){
	"samtools" "view" "-@" thr "-bSu" inputFilename @stdout=output;
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

//	file alignedsam <strcat(vars["OUTPUTDIR"],"/align/", sampleName, ".nodups.sam")>;

//	alignedsam = bwa(read1, read2, vars["BWAINDEX"], vars["BWAMEMPARAMS"], string2int(vars["PBSCORES"]), rgheader);

//	alignedsam = input(strcat(vars["OUTPUTDIR"],"/align/", sampleName, ".nodups.sam"));
//	file alignedbam <strcat(vars["OUTPUTDIR"],"/align/", sampleName, ".nodups.bam")>;
//	alignedbam = samtools(filename(alignedsam), string2int(vars["PBSCORES"]));

	//file sortedbam<SingleFileMapper; file=strcat(parameters.OUTPUTDIR,"/align/", sampleName, ".nodups.sorted.bam")>;
	//sortedbam = novosort(parameters.TMPDIR, parameters.PBSCORES, alignedbam)";
}





