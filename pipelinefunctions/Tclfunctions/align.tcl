#!/usr/bin/env tclsh

package provide align 0.0

proc bwa {bwaDir Index r1 r2 rgheader args} {
	exec $bwaDir mem $Index $r1 $r2 -R $rgheader $args >tmp.sam 
}

