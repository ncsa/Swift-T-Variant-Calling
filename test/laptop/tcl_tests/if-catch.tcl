

set someFile "tst"
if { ! [catch {open $someFile r} msg options] } {
	puts "SUCCESS"
}

puts "MESSAGE is: $msg"
puts "OPTIONS is: $options"


if {[catch {exec grep -q pippo /etc/passwd} result ]} {
	    # non-zero exit status, get it:
  set status [lindex $errorCode 2]
} else {
  # exit status was 0
  # result contains the result of your command
  set status 0
}

puts $status
puts $result
puts $errorCode
