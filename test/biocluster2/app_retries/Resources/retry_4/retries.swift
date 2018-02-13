// Defining an app that always fails:
app f() { "false" ; }



trace("*********************\t :" + "Testing app retries on (Biocluster 2)" +"\t****************************");


// Calling the always falling app
f();


