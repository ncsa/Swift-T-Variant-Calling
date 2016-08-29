app (file output) find (string dir, string pattern){
	"/usr/bin/find" dir "-name" pattern @stdout=output;
}

//@dispatch=WORKER
app () samtools_index(string samtoolsdir, string inputFilename) {
	samtoolsdir "index" inputFilename
}

