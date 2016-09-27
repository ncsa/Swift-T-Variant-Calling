app (file output) find (string dir, string pattern){
	"/usr/bin/find" dir "-name" pattern @stdout=output;
}

//@dispatch=WORKER
app (file output) samtools_index(string samtoolsdir, file inputFilename) {
	samtoolsdir "index" inputFilename @stdout=output;
}

