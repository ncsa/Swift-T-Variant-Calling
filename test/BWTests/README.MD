To get the details on the testing on Blue Waters, look at the summary.txt files in each of the folders in this directory

# Summary of Results as of April 23rd, 2018

Between October and December of 2017, Jacob was able to get successful runs of the pipeline on up to 100 samples using 110 nodes. We attempted to analyze 600 samples on 601 nodes, but it never succeeded.

Here is a summary of the October-December tests:

## Tmp dir flag problem

When trying to execute the pipeline on more than one node, Jacob quickly found that the default tmp space (/tmp) is not accessable by all nodes, and Swift/T would die as a result.

To solve this, one needs to set the swift tmp space to a custom location. The way Jacob did this was to have Swift/T automatically generate a qsub file, kill the resulting job, and manually add the following lines to the generated qsub file (This process is described in detail on the README.md):

```
#PBS -V

# Note: Make sure this directory is created before running the workflow, and make sure it is not just '/tmp'

# All three are used because it wasn't clear which one of them was the one that worked
export SWIFT_TMP=/path/to/tmp_dir
export TMPDIR=/path/to/tmp_dir
export TMP=/path/to/tmp_dir
```

After talking to Justin, it turns out that there is an `-e` flag in which one can pass in these values (as `-e <key>=<value>`), although this was not tested at this time.

## 2 Nodes
With two nodes, the run was only successful when running with a single program per node, i.e., samples cannot be doubled up on a node.

## 110 Nodes
(While this is only looking at 100 samples, Jacob added more nodes to make room for more ADLB Servers, although after talking to Justin, a single ADLB Server should be sufficient for this kind of job.)

For a long time, these tests were failing because Jacob accidentally had duplicates in the sampleinfo file.

After this was fixed, Jacob tested runs with different programs marking duplicates (Picard, Samblaster, and Novosort). 

The run with 100 nodes that finished successfully was the run where there was 1 sample per node and novosort was used to mark duplicates.

## 601 Nodes

(All of the runs below used novosort to mark duplicates, as that was the tool that worked successfully in the 110 node tests)

None of these runs ended up working.

### Tried the following tests

#### 1200 sample run (2 samples per node)
Failed in HaplotypeCaller

When inspected, GATK complains that Haplotype Caller does not have enough heap memory (this is controlled when calling Java with the -Xmx flag. It is set to 16g in the runfile, which should be plenty of memory for the data sets we are testing against.

(This is an example of the same type of error message, not from this specific run, i.e. the sample name is probably different from that in the log file reported above. But generally the error messages have been the same)

(The relevant part of the GATK error message)

\##### ERROR MESSAGE: An error occurred because you did not provide enough memory to run this program. You can use the -Xmx argument (before the -jar argument) to adjust the maximum heap size provided to Java. Note that this is a JVM argument, not a GATK argument.

#### 600 sample run (1 samples per node)

Same error as above

#### 600 sample run (1 samples per node; increased heap size from 16G to 32G (should be way more than enough memory))

Same error as above

#### Tried 300 samples instead

Still failed
Appears to fail in a different spot
It may be failing when trying to delete a temporary sam file. Not very clear from the Swift/T log

## Conclusion

This is where Jacob's testing ended.

In early January, Justin started experimenting with the workflow on BlueWaters. There are a few examples of how he changed Jacob's pipeline starting shell script to make it easier for users to use it as a template.
