#!/bin/bash

#$1 is either a directory (of small timinglogs), or a timinglog.txt

raw_input=$1 
home=$(pwd)
echo $home
if [[ -d $raw_input ]]; then
    cd $raw_input
    find . -name '*.log' -exec cat {} \; > $(pwd)/partial_run_timing.txt
    raw_logging_file=$(pwd)/partial_run_timing.txt
    cd $home

elif [[ -f $raw_input ]]; then
    raw_logging_file=$raw_input

else
    echo -e "The given arguement, **$raw_input** is NOT valid. \n\nPlease specify either:\n 1. A *Timing.log* file, or \n 2. The path to the *TMPDIR* of partial logs \n"
    echo "Exiting now!"
    exit 1
fi

Rscript summarize.R $raw_logging_file
