// Test app location dispatch

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
    string host2 = extract_hostname(@soft_location=hostmap_one_worker(host1) hostname());
    printf("\t\tHostname %i: %s\t Soft host: %s\t%i", i, host1, host2, host1!=host2);
  }
}
