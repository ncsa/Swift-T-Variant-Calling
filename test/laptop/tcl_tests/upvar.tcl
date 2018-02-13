proc top { topArg } {
	set localArg [expr $topArg+1]
	puts "Before calling bottom localArg is: $localArg"
	bottom localArg
	puts "After calling bottom, localArg is: $localArg"
}
proc bottom { bottomArg } {
	upvar $bottomArg arg
	puts "bottom is passed $bottomArg with a value of $arg"
	incr arg
}
top 2
