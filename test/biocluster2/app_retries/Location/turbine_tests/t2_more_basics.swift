import location;
//THIS-TEST-SHOULD-NOT-COMPILE
// Check that we can't give location to local task


// Should be local by default
@dispatch=WORKER 
f(int i) "turbine" "0.0.1" [
  "puts \"HELLO <<i>>\""
];

main {
  @location=location_from_rank(10)
  f(10);

  @location=location_from_rank(11)
  f(11);

  @location=randomWorker()
  f(13);
  @location=randomWorker()
  f(12);


  f(1000);

  f(2000);
  f(43974590);


}
