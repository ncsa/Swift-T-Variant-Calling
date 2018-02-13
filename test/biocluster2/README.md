These tests took place after testing on XSEDE, so for general info; you may refer to [XSEDE tests report](https://github.com/ncsa/Swift-T-Variant-Calling/tree/master/test/XSEDE_stampede2). 

Two sets of sets were done on Biocluster (which uses slurm as a resource manager):

1. H3ABioNet accreditation practice data](http://h3data.cbio.uct.ac.za/assessments/NextGenVariantCalling/practice/) of a single chromosome at 30X

2. Custom synthetic data at 30,50 and 70X generated for comparison of runtime among samples based on read depth. Both the generation scripts and variant calling scripts are provided; along with the resulting `Timing.log` file. A visualization of this log is below.

<img src=./BioclusterTiming.log.png width="700">

