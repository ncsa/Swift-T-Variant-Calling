import assert;
import files;
import io;
import string;
import location;

app (file o) hostname() {
  "hostname" @stdout=o;
}

(string o) extract_hostname(file f) {
  o = trim(read(f));
}

main {
  foreach i in [1:500] {
	    string host1 = extract_hostname(hostname());
	    int rank = hostmapOneRank(host1);   //change to hostmapOneWorkerRank to gaurantee smooth testing
	    location loc = location(rank, HARD, NODE);
	    string host2 = extract_hostname(@location=loc hostname());
	    printf("\t%i: rank: %i Hostname: %s, Dispatched host: %s, identical: %i", i, rank, host1, host2, host1==host2);

    /* Check that apps run on the correct host! */
    //assertEqual(host2, host1, "hosts");
  }
}

