import "stdlib.v2";

type file;

string config = arg("fileName", "config");
file configFile <"H3A_NextGen_assessment.Chr1_50X.set3.runfile">;

app (file trimmed) fileTrimmer (file untrimmed){
	trim @untrimmed stdout=fileName(trimmed);
}

file trimmedConfig = fileTrimmer(configFile);

string fileData = read(f = trimmedConfig, format = "CSV");
