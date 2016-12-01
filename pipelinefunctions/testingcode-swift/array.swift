import files;
import unix;
 
string arr [] = ["adele", "book", "clock"];
string b[];
int c;
for (int i=0; i<size(arr); i=i+1) {
	trace ("#########\t\t" + arr[i] + "#########");
	b[i] = arr[i];
}

foreach l in arr {
	trace ("*********\t\t" + l + "***********");
	c = 1;
	c = 4;
}
file tmp <"tmp.swift">;
tmp = touch ();

if (size(b) == size(glob("*.swift"))) {
	trace ("---------\t" +"Size match"+ "\t-----------");
	trace ("---------\t" + size(b) + "\t-----------");
	trace ("---------\t" + size(arr) + "\t-----------");
}


