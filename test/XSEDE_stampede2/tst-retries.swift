// Defining an app that always fails:
app f() { "false" ; }



trace("*********************\t :" + "Testing app retries on XSESE (stampede2)" +"\t****************************");


// Calling the always falling app
f();


