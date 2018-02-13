import io;
import location;

type wrapped_rank int;
type my_type float;

main {
  location A[];
  foreach i in [1:10] {
    A[i] = random_worker();
  }

  foreach j in [1:10] {
	  printf("\t\t: Default code:\t %d: %d", j, wrapped_rank(A[j].rank));
	  printf("\t\t: My code, no type:\t %d: %d", j, A[j].rank);
  }

  trace(wrapped_rank(5));
  trace(my_type(25));
}
