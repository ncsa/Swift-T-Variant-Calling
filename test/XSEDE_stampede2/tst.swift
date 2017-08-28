import files;
import string;

string chr = "2222";

app f() { "false" ; }
f();


file bam = input("/home/azza/swift-project/Results/run2/HG00108.lowcoverage.chr20.smallregion/align/HG00108.lowcoverage.chr20.smallregion.wDedups.sorted.bam");

string prefix = substring( filename(bam), 0, strlen(filename(bam)) - 4 );
string fileName = basename_string(prefix);
string sampleName = substring( fileName, 0, strlen(fileName) - 15);

string outputName = strcat("OUTPUTDIR", "/", sampleName, "/realign/", fileName, ".",
                                                                        chr, ".bam"
                                                                        );
trace("*********************\t sampleName:" + sampleName + "\t****************************");
trace("*********************\t fileName:" + fileName + "\t****************************");
trace("*********************\t outputName:" + outputName + "\t****************************");

