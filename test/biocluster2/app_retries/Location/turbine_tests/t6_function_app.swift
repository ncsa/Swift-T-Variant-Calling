// Test app location dispatch

import assert;
import files;
import io;
import string;
import location;

@dispatch=WORKER
app (file o) hostname() {
  "hostname" @stdout=o;
}

@par
(string o) extract_hostname(file f) {
  o = trim(read(f));
}

main {
  foreach i in [1:10] {
    string host1 = extract_hostname(hostname());
    // Run on same host
    string host2 = @par=1 extract_hostname(@location=hostmap_one_worker(host1) hostname());
    string host3 = extract_hostname(@location=locationFromRank(i) hostname());
    printf("\t%d\t\tRandom host:\t%s\tHost by hostmap:\t%s\tHost by rank:\t%s", i, host1, host2, host3);
   assertEqual(host1, host2, sprintf("Check hostnames same trial %i", i));
  }
}
