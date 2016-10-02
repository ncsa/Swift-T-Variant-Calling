import files;
import string;
import unix;

app (file output) find (string dir, string pattern){
        "find" dir "-name" pattern @stdout=output;
}

file fining <"f.txt">  = find(strcat("/home/azza/swift-project/Workshop_low_pass/ref", "/", "IndelsByChr", "/"), "*20.vcf");

file new <"n.txt">  = sed(fining, "s/^/ --knownSites /g");

string result = replace_all(read(input_file(filename(new))), "\n", " ", 0);
trace("\n\n" + result + "\n\n");



/*

foundfiles = glob(strcat("/home/azza/swift-project/Workshop_low_pass/ref", "/", "IndelsByChr", "/","*20.vcf"));

//parmsfiles = sed(foundfiles[], "s/^/ --knownSites/g");

string ff[];
for (int i=0; i<1; i=i+1) {
	f= filename(foundfiles[i]);
	ff[i] = strcat(" --knownSites ", f);
}

foreach file in foundfiles {
	ff = filename(file);
//	ff  strcat(" --knownSites ", f);
}


recalparms2 = join(ff," ");
trace(recalparms2);

//recalparms2 = replace_all(foundfiles, "^", " --knownSites", 0);
//trace("\n\n", recalparms2);


*/
