proc errorProc { first second } {
	global errorInfo

	# $fail will be non−zero if $first is non−numeric.
	set fail [catch { expr 5 * $first } result]
	puts "First argument check: $fail"

	# if $fail is set, generate an error
	if { $fail } {
		error "Bad first argument"
	}

	# This will fail if $second is non−numeric or 0
	set fail [catch { expr $first/$second } dummy]
	puts "Second argument check: $fail"

	if { $fail } {
		error "Bad second argument" \
		"second argument fails math test\cback \n
		$errorInfo"
	}

	error "errorProc always fails" "evaluating error" \
	[list USER {123} { Non−Standard User−Defined Error}]
}
# Example Script

set fail [catch { errorProc X 0} returnString]

puts "Overall check: $fail"

if { $fail } {
	puts "Failed in errorProc"
	puts "Return string: $returnString"
	puts "Error Info: $errorInfo\n"
}

puts "call errorProc with a 0 second argument"
if {[catch { errorProc 1 0} returnString]} {
	puts "Failed in errorProc"
	puts "Return string: $returnString"
	puts "Error Info: $errorInfo\n"
}

puts "call errorProc with valid arguments"
set fail [catch { errorProc 1 1} returnString]
if { $fail } {
	if {[string first USER $errorCode] == 0} {
		puts "errorProc failed as expected"
		puts "returnString is: $returnString"
		puts "errorInfo: $errorInfo"
} else {
	puts "errorProc failed for an unknown reason"
}
}
