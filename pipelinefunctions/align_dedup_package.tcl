#!/usr/bin/env tclsh

package provide align 0.2

namespace eval alignment {
	proc samtools_view {samtoolsdir inputFile} {
#		set numAlignments [ open "|$samtoolsdir view -c $inputFile" w+ ]
	set numAlignments [ exec $samtoolsdir view -c $inputFile ]

	}
}

#set samtoolsdir {/usr/local/bin/samtools}
#set input /home/azza/swift-project/Swift-Variant-Calling/pipelinefunctions/tmp.sam

#set num [ alignment::samtools_view $samtoolsdir $input ]
#puts "num = $num"
#app (file output) novosort (string novosortdir, file inputFile, string tmpdir, int thr, string sortoptions[]){
	#// processing a single file (sorting and indexing input)
#	novosortdir "--index" "--tmpdir" tmpdir "--threads" thr inputFile "-o" output sortoptions; 
#	// novosort has dual function to also mark duplicates
#}

#@dispatch=WORKER
#app (file output) novosort (string novosortdir, string inputFile[], string tmpdir, int thr, string sortoptions[]){
#	// processing multi-input files together (merging files)
#	novosortdir "--index" "--tmpdir" tmpdir "--threads" thr inputFile "-o" output; 
#	// novosort has dual function to also mark duplicates
#}



