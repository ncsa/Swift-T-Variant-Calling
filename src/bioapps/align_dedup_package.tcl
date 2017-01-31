#!/usr/bin/env tclsh

package provide align 0.2

namespace eval alignment {
	proc samtools_view {samtoolsdir inputFile} {
	# set numAlignments [ open "|$samtoolsdir view -c $inputFile" w+ ]
	set numAlignments [ exec $samtoolsdir view -c $inputFile ]

	}
}
