#!/bin/bash
set -eu

# UNUSED alternative implementation of logging()

TMPDIR=$1
TIMELOG=$2
TOOLDIRNAME=$3

TMPLOGS=$( ls $TMPDIR/timinglogs/$TOOLDIRNAME )
if (( ${#TMPLOGS} ))
then
  cat $TMPLOGS >> $TIMELOG
  rm -v $TMPLOGS
fi
echo "Updated $TIMELOG"

# file tmplogs[] = glob(strcat(tmpdir, "/timinglogs/", toolDirName, "/*"));
# if (size(tmplogs) > 0) {
# 	append(timeLog, read(cat(tmplogs))); //=>

# 	/*
# 	These rm calls caused bugs

# 	*/
# 	/*foreach i in tmplogs {
# 		rm(i);
# 	}*/
# }
