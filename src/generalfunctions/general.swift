// These are functions necessary for the pipeline, but hard to classify into a job.

import unix;
import files;
import string;
import assert;

import bioapps.align_dedup;

(boolean signal) append(file f, string s) {
	string fName = filename(f);
	string original = read(f);
	string appended = strcat(original, s);
	file newFile < fName >;
	newFile = write(appended) =>
	signal = true;
}

exitIfFlagGiven(string vars[string], string message) {
	// If user doesn't want to kill job if a sample fails
	if (vars["EXIT_ON_ERROR"] == "F" ||
	    vars["EXIT_ON_ERROR"] == "f" ||
	    vars["EXIT_ON_ERROR"] == "False" ||
	    vars["EXIT_ON_ERROR"] == "false" ||
	    vars["EXIT_ON_ERROR"] == "FALSE" ||
	    vars["EXIT_ON_ERROR"] == "N" ||
	    vars["EXIT_ON_ERROR"] == "n" ||
	    vars["EXIT_ON_ERROR"] == "No" ||
	    vars["EXIT_ON_ERROR"] == "NO" ||
	    vars["EXIT_ON_ERROR"] == "no"
	   ) {
		// Do nothing. The user elected to not kill upon failure
	} else {
		assert(false, message);
	}
}

// Reading the runfile parameters:
(string data[string]) getConfigVariables(string lines[])
{
	foreach line in lines
	{
		string keyValuePair[] = split(line, "=");
		string name = keyValuePair[0];
		string value = keyValuePair[1];
		data[name] = value;
	}
}

app (file output) find_files (string dir, string pattern){
	"find" dir "-name" pattern @stdout=output;
}

app (void v) rm(file f) {
	"rm" f;
}

app (void v) rm(file f[]) { 
	"rm" f;													 
}

() logging (string tmpdir, file timeLog){
        file tmplogs[] = glob(strcat(tmpdir, "/timinglogs/*"));
        append(timeLog, read(cat(tmplogs))) =>
        rm(tmplogs);
}

// Convert an array of files to an array of strings with those file names
(string filenames[]) filesToFileNames(file files[]) {
	foreach f, index in files {
		fName = filename(f);
		filenames[index] = fName;
	}
}

(boolean success) checkBam(string vars[string], file bamFile) {						       
	// Swift/T would try to use samtools_view2 even before bamFile is read. Added explicit wait command		
	wait (bamFile) {
		int alignNum = samtools_view2(vars["SAMTOOLSEXE"], filename(bamFile));				
		success = alignNum > 0;										    
	}
} 
