
//app (file out) sort (file in){
//	"sort" in @stdout=out;
//}

app (file out) ls (){
	"ls" "-htr"  @stdout=out;
}


file t <"test.txt">;
t= ls();// =>
	//sort();
//ls -tr1 | sort
