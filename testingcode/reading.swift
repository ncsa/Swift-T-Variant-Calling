import files;
import unix;
import string;
variantFiles = find_files ("/home/azza/swift-project/Results/delivery", strcat("*.GATKCombineGVCF.raw.vcf")) ;

string varlist[] = split(trim(replace_all(read(variantFiles), "\n", " ", 0)), " ");

foreach v in varlist {
	trace ("\n**********\t\t"+ v+"\t\t***********");
}

//jointVCF = GenotypeGVCFs (strcat("/usr/local/src/gatk/GenomeAnalysisTK.jar"), "/home/azza/swift-project/Workshop_low_pass/ref/human_g1k_v37_chr20.fa", varlist );



@dispatch=WORKER
app (file outputfile) GenotypeGVCFs (string gatk, string reference, string variants[] ) {
        "/usr/bin/java" "-Xmx8g" "-jar" gatk "-T" "GenotypeGVCFs" "-R" reference  variants "-o" outputfile;
}

app (file output) find_files (string dir, string pattern){
        "find" dir "-name" pattern @stdout=output;
}

