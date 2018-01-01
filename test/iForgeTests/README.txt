This was written by Matt Kendzior on Nov 27th, 2017 to describe the iForge testing process:

I started testing the Swift/T workflow using two samples, just to validate that the functionalities and variant calling tools worked on iForge. The TestCases.txt file in the test directory lists the functionalists and tools for validation. The runfiles and qsub files in this test directory are named after successful tests. I was able to debug by first checking the Swift/T stdout and stderr logs produced by the PBSTorque to validate that the features worked. The error messages were sufficient to realize what was wrong with the Swift/T workflow source code and to fix it. If the error was caused by a tool failing on a sample, this could then be debugged looking at the error logs for the tool in the logs directory of the samples output directory.

_________ 2 sample paired end reads  (BWA-mem as aligner, Picard for marking duplicates)_________ 
This test run was with two samples sequenced with paired end reads, so that the workflow could run to completion, from the alignment/deduplication stage through the joint genotyping stage, but still be run on one node. Two sample processes per node was the recommended amount on iForge. BWA-mem was used as the aligner and Picard was used for marking duplicates The workflow failed on several steps. This test led us find and debug an error with Picard, the VC_NO_REALIGN function, errors running to completion with the Y set variables for each stage. The runfile after dubugging these errors is: 5,8c,19_SUCCESS.runfile.

_________2 sample paired end reads (Novoalign as aligner, samblaster for marking duplicates)_________
The test was done with the same two paired end samples on one node as the first run, but just to test Novoalign as the aligner, and samblaster for marking duplicates. I also validated that the workflow could successfully perform the alignment stage only. Both of these tools originally failed, but the logs provided in the sample logs output directory were sufficient for debugging. I also validated that the workflow could successfully perform the alignment stage only. The runfile after dubugging these errors is: 2,4,8a_SUCCESS.runfile. 

_________2 sample single end reads_________
This test was done to validate the workflow could handle samples sequenced with single ended reads.

_________Validate that workflow will not run if no samples pass the quality control check_________
This test was run on 16 samples, with two samples per node, to test whether the workflow will stop if no samples pass quality control checks. For this test we set the MAP_CUTOFF to 101, which lead each sample to fail the quality control and the workflow to successfully exit. The runfile from this test is 12_SUCCESS.runfile

_________Validate that workflow will run if one or more samples fail, but not all._________
This test used 12 samples, with 2 samples per node, to validate that if the EXIT_ON_ERROR flag was set to N, the workflow will continue running through for the samples that did not fail.
The runfile from this test is 11_SUCCESS.runfile

_________Validate all other features and tools on a 6 sample run_________
I then tested out all of the other various flags and tools with different parameters using 6 samples, 2 per node. This test run validated many different features, and these are listed in the name of the runfile:
1,3b,6,7,8b,9,13,15,16,17a,18_SUCCESS.runfile.

