import sys;
import files;
import string;

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

app (file output) bwa (string params[string], int PBSCORES, string read1, string read2, string rgheader){
	 "bwa" "mem" params["BWAMEMPARAMS"] "-t" PBSCORES "-R" rgheader  params["BWAINDEX"] read1 read2 @stdout=output;
}

app (file output) samtools(string inputFilename, int thr){
	"samtools" "view" "-@" thr "-bSu" inputFilename @stdout=output;
}

//get Configuration filename with argument --params [filename]

string configFilename = argv("params", "runfile");
file configFile = input_file(configFilename);
string configFileData[] = file_lines(configFile);
string vars[string] = getConfigVariables(configFileData);

file sampleInfoFile = input_file(vars["SAMPEINFORMATION"]);
string sampleLines[] = file_lines(sampleInfoFile);

foreach sample in sampleLines{

	string sampleInfo[] = split(sample, " ");
	string sampleName = sampleInfo[0];
	string read1 = sampleInfo[1];
	string read2 = sampleInfo[2];

	string rgheader = sprintf("@RG\tID:%s\tLB:%s\tPL:%s\tPU:%s\tSM:%s\tCN:%s\t", sampleName, vars["SAMPLELB"], vars["SAMPLEPL"], sampleName, sampleName, vars["SAMPLECN"]);

	file alignedsam <strcat(vars["OUTPUTDIR"],"/align/", sampleName, ".nodups.sam")>;
	alignedsam = bwa(vars, string2int(vars["PBSCORES"]), read1, read2, rgheader);

	file alignedbam <strcat(vars["OUTPUTDIR"],"/align/", sampleName, ".nodups.bam")>;
	alignedbam = samtools(filename(alignedsam), string2int(vars["PBSCORES"]));

	//file sortedbam<SingleFileMapper; file=strcat(parameters.OUTPUTDIR,"/align/", sampleName, ".nodups.sorted.bam")>;
	//sortedbam = novosort(parameters.TMPDIR, parameters.PBSCORES, alignedbam)";
}
