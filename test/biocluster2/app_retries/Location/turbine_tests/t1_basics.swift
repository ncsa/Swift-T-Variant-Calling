import random;
import location;
import assert;
import io;
import files;

@dispatch=WORKER
(int rank) f(int i) "turbine" "0.0.1" [
  "set <<rank>> [ adlb::rank ]; after 5"
];

app (file o) hostname() {
	  "hostname" @stdout=o;
}

int N = 50;

int actual_ranks[];
file hosts[];

foreach i in [1:N] {
    int target_rank = randomWorkerRank();
    location target = locationFromRank(target_rank);

    actual_rank = @location=target f(0);
    host = @location=target hostname();

    trace(i + ":Target rank:\t" + target_rank + "\tActual rank:\t" + actual_rank + "\tActual host:\t" + read(host));

    hosts[i] = host;
    actual_ranks[i] = actual_rank;
  }

