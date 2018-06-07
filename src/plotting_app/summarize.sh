#!/bin/bash

#$1 is either a directory (of small timinglogs), or a timinglog.txt

raw_input=$1 

if [[ -d $raw_input ]]; then
    cd $raw_input
    find . -name '*.txt' -exec cat {} \; > partial_run_timing.log
    raw_logging_file=partial_run_timing.log

elif [[ -f $raw_input ]]; then
    raw_logging_file=$raw_input
else
    echo -e "The given arguement, **$raw_input** is NOT valid. \n\nPlease specify either:\n 1. A *Timing.log* file, or \n 2. The path to the *TMPDIR* of partial logs \n"
    echo "Exiting now!"
    exit 1
fi

echo 
Rscript summarize.R $raw_logging_file
