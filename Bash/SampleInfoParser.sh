#!/bin/bash

line=0

touch "cleaned.sampleinfo"
while read sampledetails;
do
	if [ ! `expr ${#sampledetails}` -lt 1 ]; then

		SampleName=${sampledetails%% *}
		read1=${sampledetails% *}
		read2=${sampledetails##* }

		echo "[${line}].SampleName=${SampleName}"
		echo "[${line}].read1=${read1}"
		echo "[${line}].read2=${read2}"
		let line+=1
	fi
done < "$1"
