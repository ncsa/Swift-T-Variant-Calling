// These are functions necessary for the pipeline, but hard to classify into a job.

import unix;
import files;
import string;

append(file f, string s) {
  string fName = filename(f);
  string original = read(f);
  string appended = strcat(original, s);
  file newFile < fName >;
  newFile = write(appended);
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

// Convert an array of files to an array of strings with those file names
(string filenames[]) filesToFileNames(file files[]) {
	foreach f, index in files {
		fName = filename(f);
		filenames[index] = fName;
	}
}

(boolean success) checkBam(string variables[string], file bamFile) {						       
	// Swift/T would try to use samtools_view2 even before bamFile is read. Added explicit wait command		
	wait (bamFile) {
		int alignNum = samtools_view2(variables["SAMTOOLSDIR"], filename(bamFile));				
		success = alignNum > 0;										    
	}
} 
