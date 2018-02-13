// THIS-TEST-SHOULD-NOT-COMPILE

import io;
import location;
import string;

//@dispatch=WORKER
@par
() f_init(int n) {
	printf("\tHello: %i", n);
}

@dispatch=WORKER
(int rank) f(int i) "turbine" "0.0.1" [
  "set <<rank>> [ adlb::rank ]; after 5"
  ];

// Integer not location
@location=locationFromRank(0)  
 f(10);
@par=4 f_init(40);
