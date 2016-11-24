package provide mytcl 0.0

namespace eval mytclsort {
	proc novosort { novosortdir, inputFile, tmpdir, thr, output, args } {
		 $novosortdir "--index" "--tmpdir" $tmpdir "--threads" $thr $inputFile "-o" $output $args;
	}
}
	
