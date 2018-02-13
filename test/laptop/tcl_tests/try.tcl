set f [open tst a]
try {
	    puts $f "some message"
	        
} finally {
     	close $f
}
