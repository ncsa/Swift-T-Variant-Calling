type file;

type configData {
	string SAMPLEINFORMATION;
	string OUTPUTDIR;
	string TMPDIR;
	string SCRIPTDIR;
	string EMAIL;
	int REPORTTICKET;
	string RUNMETHOD;
	string ANALYSIS;
	string INPUTFORMAT;
	int PAIRED;
	int READLENGTH;
	string MULTISAMPLE;
	string PROVENANCE;
	string SAMPLEID;
	string SAMPLELB;
	string SAMPLEPL;
	string SAMPLEPU;
	string SAMPLESM;
	string SAMPLECN;
	string ALIGNERTOOL;
	string SORTMERGETOOL;
	string MARKDUPLICATESTOOL;
	string BWAMEMPARAMS;
	string CHRNAMES;
	string BWAINDEX;
	int MAP_CUTOFF;
	int DUP_CUTOFF;
	string MARKDUP;
	string REFGENOMEDIR;
	string REFGENOME;
	string DBSNP;
	string INDELDIR;
	string OMNI;
	string SAMBLASTERDIR;
	string PICARDIR;
	string GATKDIR;
	string SAMDIR;
	string BWAMEMDIR;
	string TABIXDIR;
	string JAVADIR;
	string NOVOCRAFTDIR;
	string VCFTOOLSDIR;
	string FASTQCDIR;
	string DELIVERYFOLDER;
	string PBSPROJECTID;
	int PBSNODES;
	int PBSCORES;
	string PBSQUEUE;
	string PBSWALLTIME;
}
type sampleInfo {
	string SampleName;
	string read1;
	string read2;
}

app (file alignedsam) bwa (configData info, sampleInfo sample, string rgheader){
	 bwa "mem" info.BWAMEMPARAMS "-t" info.PBSCORES "-R" rgheader  info.BWAINDEX sample.read1 sample.read2 stdout=filename(alignedsam);
}
app (file alignedbam) samtools(file input, int thr){
	samtools "view" "-@" thr "-bSu" filename(input) stdout=filename(alignedbam);
}
#app (file dedupsortedbam) novosort(string tempDir, int threads, file input){
#	 novosort "--index" "--tmpdir" tempDir "--threads" threads input stdout=filename(dedupsortedbam)
#}

string parametersFilename = arg("params", "runfile");
file configFile<SingleFileMapper; file=parametersFilename>;
configData parameters = readStructured(filename(configFile));
file sampleInfoFile<SingleFileMapper; file = parameters.SAMPLEINFORMATION>;
sampleInfo[] samples = readData(sampleInfoFile);

foreach sample in samples{
	string rgheader = sprintf("@RG\tID:%s\tLB:%s\tPL:%s\tPU:%s\tSM:%s\tCN:%s\t", sample.SampleName, parameters.SAMPLELB, parameters.SAMPLEPL, sample.SampleName, sample.SampleName, parameters.SAMPLECN);
	file alignedsam<SingleFileMapper; file=strcat(parameters.OUTPUTDIR,"/align/", sample.SampleName, ".nodups.sam")>;
	alignedsam = bwa(parameters, sample, rgheader);
	file alignedbam<SingleFileMapper; file=strcat(parameters.OUTPUTDIR,"/align/", sample.SampleName, ".nodups.bam")>;
	alignedbam = samtools(alignedsam, parameters.PBSCORES);
	#file sortedbam<SingleFileMapper; file=strcat(parameters.OUTPUTDIR,"/align/", sample.SampleName, ".nodups.sorted.bam")>;
	#sortedbam = novosort(parameters.TMPDIR, parameters.PBSCORES, alignedbam)";
}
