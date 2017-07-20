#!/usr/bin/env tclsh

package provide align 0.2

namespace eval alignment {
	proc samtools_view {samtoolsexe inputFile} {
	# set numAlignments [ open "|$samtoolsexe view -c $inputFile" w+ ]
	set numAlignments [ exec $samtoolsexe view -c $inputFile ]

	}
}
