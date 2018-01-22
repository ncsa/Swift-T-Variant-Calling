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

/*
 * Takes info from temporary log files and creates the final timing log file; also deletes the temporary ones
 */
() logging (string tmpdir, file timeLog){
        file tmplogs[] = glob(strcat(tmpdir, "/timinglogs/*"));
        append(timeLog, read(cat(tmplogs))) =>
        foreach i in tmplogs {
		rm(i);
	}
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

// Modifying the pure swift/t function to accept strings as well
@pure
(string t) file_type(string f)
"turbine" "0.0.2"
[ "set <<t>> [ file type [ lindex <<f>> 0 ] ]" ];


(boolean exec_ok) exec_check (string exec, string parameter){
        file_exists(exec) =>
	string fileType = file_type(exec) =>
	assert(fileType == "file" || fileType == "link", 
	       strcat("The executable: \n\t", exec,
		      "\n Referred to by the parameter: \n\t", parameter, 
		      " is not properly specified in your runfile!"
		     )
	      );
        exec_ok = true;
}

