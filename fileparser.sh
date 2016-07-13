#!/bin/bash

# run this file as: fileparser.sh <some_file>
# the output would be a file stripped from empty lines and comments


#remove empty lines from input:
sed -i '/^\s*$/d' $1

#remove files starting with a comment:
sed -i '/#/d' $1
