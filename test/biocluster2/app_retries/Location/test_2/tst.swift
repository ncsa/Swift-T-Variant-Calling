import location;
import string;

// Defining an app that always fails:
app f() { "false" ; }

trace("*********************\t :" + "Testing app retries on (Biocluster 2)" +"\t****************************");


/* type struct location {int rank; LocationStrictness strictness; LocationAccuracy accuracy} */


string hosts[] = hostmapList();

trace("Number of hosts: \t" + size(hosts));

trace("Host_number\t"+ "Host_name\t" + "Rank\t" + "Strictness\t" + "Accuracy" )=>
foreach i in [0: size(hosts)-1]{
	loc = hostmapOne(hosts[i]) ;
	trace(i + "\t"+ hosts[i] + "\t" + loc.rank + "\t" + loc.strictness + "\t" + loc.accuracy);
}

trace("Host_number\t"+ "WorkerRank\t" + "Rank\t" + "Strictness\t" + "Accuracy" ) =>
foreach i in [0: size(hosts)-1]{
	rank = hostmapOneWorkerRank(hosts[i]); 
	loc = rank2location(rank);
	trace(i + "\t"+ hosts[i] + "\t" + rank + "\t" + loc.rank + "\t" + loc.strictness + "\t" + loc.accuracy);
}



