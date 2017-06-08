@dispatch=WORKER
app (file outputfile, file outLog) CombineGVCFs (string javaexe, string gatkdir, string reference, string dbsnp, string variants[] ) {
        javaexe "-Xmx8g" "-jar" gatkdir "-T" "CombineGVCFs" "-R" reference "--dbsnp" dbsnp variants "-o" outputfile @stderr=outLog;
}


