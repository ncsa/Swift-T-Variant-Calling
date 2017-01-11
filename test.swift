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

file example <"/ui/ncsa/jacobrh/example.txt">;

example = write("line1\n");
append(example, "Line2\n");
