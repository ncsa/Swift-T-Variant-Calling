@dispatch=WORKER
app (file outputfile) GenotypeGVCFs (string javadir, string gatkdir, string reference, string variants[] ) {
        javadir "-Xmx8g" "-jar" gatkdir "-T" "GenotypeGVCFs" "-R" reference  variants "-o" outputfile;
}


