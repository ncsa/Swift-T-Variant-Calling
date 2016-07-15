#!/bin/bash

# run this file as: fileparser.sh <some_file>
# the output would be a file stripped from empty lines and comments

cp $1 cleanFile
#remove empty lines from input:
sed -i '/^\s*$/d' cleanFile

#remove files starting with a comment:
sed -i '/#/d' cleanFile

cat cleanFile
