@dispatch=WORKER
app (file outputfile, file outLog) GenotypeGVCFs (string javaexe, string java_heap, string gatkjar, string reference, string variants[], string threads) {
        javaexe java_heap "-jar" gatkjar "-T" "GenotypeGVCFs" "-R" reference  variants "-o" outputfile "-nt" threads  @stderr=outLog;
}


