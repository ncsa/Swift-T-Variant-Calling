@dispatch=WORKER
app (file outputfile, file outLog) GenotypeGVCFs (string javadir, string gatkdir, string reference, string variants[] ) {
        javadir "-Xmx8g" "-jar" gatkdir "-T" "GenotypeGVCFs" "-R" reference  variants "-o" outputfile @stderr=outLog;
}


