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

file configFile<"runfile">;
configData parameters = readStructured(filename(configFile));
file sampleInfoFile<SingleFileMapper; file = parameters.SAMPLEINFORMATION>;
sampleInfo[] samples = readData(sampleInfoFile);



