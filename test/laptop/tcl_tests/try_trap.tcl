try {
    set f [open /usr/src w]
} trap { POSIX EISDIR } { msg } {
    puts "failed to open tst: it's a directory"
} trap { POSIX ENOENT } { msg } {
    puts "failed to open tst: it doesn't exist"
}

puts "\nExecution result from the try script \n \t $msg"
