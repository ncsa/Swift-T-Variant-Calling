@dispatch=WORKER
app (file outputfile, file outLog) CombineGVCFs (string javaexe, string java_heap, string gatkjar, string reference, string dbsnp, string variants[] ) {
        javaexe java_heap "-jar" gatkjar "-T" "CombineGVCFs" "-R" reference "--dbsnp" dbsnp variants "-o" outputfile @stderr=outLog;
}


