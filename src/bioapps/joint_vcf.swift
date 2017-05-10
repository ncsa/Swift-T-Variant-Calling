@dispatch=WORKER
app (file outputfile, file outLog) GenotypeGVCFs (string javadir, string gatkdir, string reference, string variants[], string threads) {
        javadir "-Xmx8g" "-jar" gatkdir "-T" "GenotypeGVCFs" "-R" reference  variants "-o" outputfile "-nt" threads  @stderr=outLog;
}


