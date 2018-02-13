
import files;
import io;
import string;
import location;

/**
   Run as TURBINE_SRAND=<seed> turbine ... 568.tcl
*/

app (file o) hostname()
{
  "hostname" @stdout=o;
}

@dispatch=WORKER
f(int i) "turbine" "0.0.1" [
  "puts \"HELLO <<i>>\""
  ];

main
{
  file tmp<"tmp.txt"> = @location=locationFromRank(12) hostname();
  string name = trim(read(tmp));
  printf("\t\tname: %s", name);
  location rank = hostmap_one(name);
  printf("\t\trank: %i", rank.rank);
  @location=locationFromRank(12) f(12);
  location worker_rank_2 = hostmap_one_worker(name);
  printf("\t\tworker_rank_2: %i", worker_rank_2.rank);
}
