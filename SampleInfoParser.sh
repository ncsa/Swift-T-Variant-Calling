#!/bin/bash

line=0

if [ -s cleaned.sampleinfo ]; then
	truncate -s 0 cleaned.sampleinfo
fi

while read sampledetails; 
do
	if [ ! `expr ${#sampledetails}` -lt 1 ]; then 

		SampleName=${sampledetails%% *}
		read1=${sampledetails% *}
		read2=${sampledetails##* }

		echo "[${line}]SampleName=${SampleName}" >> cleaned.sampleinfo
		echo "[${line}]read1=${read1}" >> cleaned.sampleinfo
		echo "[${line}]read2=${read2}" >> cleaned.sampleinfo
		let line+=1
	fi
done < "$1"
