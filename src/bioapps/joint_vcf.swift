@dispatch=WORKER
app (file outputfile, file outLog) GenotypeGVCFs (string javaexe, string gatkjar, string reference, string variants[], string threads) {
        javaexe "-Xmx8g" "-jar" gatkjar "-T" "GenotypeGVCFs" "-R" reference  variants "-o" outputfile "-nt" threads  @stderr=outLog;
}


