@dispatch=WORKER
app (file outputfile, file outLog) GenotypeGVCFs (string javaexe, string gatkdir, string reference, string variants[], string threads) {
        javaexe "-Xmx8g" "-jar" gatkdir "-T" "GenotypeGVCFs" "-R" reference  variants "-o" outputfile "-nt" threads  @stderr=outLog;
}


