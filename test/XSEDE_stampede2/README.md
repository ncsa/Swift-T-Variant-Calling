This file aims to list a few details on the tests performed on XSEDE resources and/or explain what the files represent. Note that Stampede2 runs Slurm as a resource manager.

## Job submission:
This can be done by composing a slurm script and then submitting the job explicitly to the scheduler. An example is:

```
$ vi submit.tst.swift.slurm
# A standard sbatch file (starting with SBATCH declarations: e.g #SBATCH -J tstSwift), then followed by the command to launch Swift/T

$ sbatch submit.tst.swift.slurm
# Submitting this command returns the job id, and also, depending on the SBATCH directives in the submission script,
# in this case it  creates 2 files in the current directory: one for the errors (myjob.e###) and another for the output 
# (myjob.o###); where ### is the jobid
```

Another alternative, is to call swift/t directly:

```
$ swift-t -l -n 5 tst.swift
# In this case, the job is submitted to the cluster right away.
```

Note that in either case, the resources allocated to the Swift/T job are **NOT** those specified by the SBATCH directives. Swift/T job resources need to be explicitly specified via standard Swift/T `export` variables and/or command line flags

## Environment settings:

These files are inside the `data_preparation` folder, and essentially generate some of the needed resources for running this Swift/T pipeline on Stampede environment. They are intended to be submitted to slurm directly, not necessarily through swift/t jobs (just to ease porting to other slurm based environments)
In particular, one file (`prepare_reference_genome.slurm.sbatch.sh`) does the indexing and needed preparations for the reference genome (as per the GATK recommendations), and the other (`split_vcf_by_chromosome.slurm.sbatch.sh`) splits up the 1000 Genomes and Mills files from the GATK bundle per chromosome, to faciliate splitting the stages of Realignment/Recalibration/HaplotypeCalling per chromosome while running this pipeline.

## GIAB experiments:
Tests on XSEDE resources used the following scenarios (note that some of these tests date back to April 2017. Since then, some parameter names in the runfile have changed throughout development; and also that some of these tests detected errors with how apps are being called. For example, the Read Group format in bwa 0.7.16 has changed from previous releases [For example](https://github.com/ncsa/Swift-T-Variant-Calling/issues/18).  ). 

1. ["ANALYSIS"]=<ALIGN>, ["SPLIT"]=<YES>, ["PAIRED"]=<1>, ["EXIT_ON_ERROR"]=<YES>, ["ALIGN_DEDUP_STAGE"]=<Y>, ["CHR_SPLIT_STAGE"]=<Y>, ["VC_STAGE"]=<Y>, ["COMBINE_VARIANT_STAGE"]=<Y>, ["JOINT_GENOTYPING_STAGE"]=<Y>, ["ALIGNERTOOL"]=<NOVOALIGN>, ["MARKDUPLICATESTOOL"]=<PICARD>, ["CHRNAMES"]=<1:2:3> ==> 12 MPI ranks, 3 nodes

2. ["ANALYSIS"]=<ALIGN>, ["SPLIT"]=<YES>, ["PAIRED"]=<1>, ["EXIT_ON_ERROR"]=<YES>, ["ALIGN_DEDUP_STAGE"]=<Y>, ["CHR_SPLIT_STAGE"]=<Y>, ["VC_STAGE"]=<Y>, ["COMBINE_VARIANT_STAGE"]=<Y>, ["JOINT_GENOTYPING_STAGE"]=<Y>, ["ALIGNERTOOL"]=<BWAMEM>, ["MARKDUPLICATESTOOL"]=<PICARD>, ["CHRNAMES"]=<1:2:3> ==> 12 MPI ranks, 4 nodes

3. ["ANALYSIS"]=<ALIGN>, ["SPLIT"]=<YES>, ["PAIRED"]=<1>, ["EXIT_ON_ERROR"]=<YES>, ["ALIGN_DEDUP_STAGE"]=<Y>, ["CHR_SPLIT_STAGE"]=<Y>, ["VC_STAGE"]=<Y>, ["COMBINE_VARIANT_STAGE"]=<Y>, ["JOINT_GENOTYPING_STAGE"]=<Y>,  ["ALIGNERTOOL"]=<BWAMEM>, ["MARKDUPLICATESTOOL"]=<PICARD>, ["CHRNAMES"]=<1:2:3> ==> 4 MPI ranks, 2 nodes

4. ["REALIGN"]=<NO>, ["SPLIT"]=<YES>, ["EXIT_ON_ERROR"]=<YES>, ["ALIGN_STAGE"]=<Y>, ["DEDUP_SORT_STAGE"]=<Y>, ["CHR_SPLIT_STAGE"]=<Y>, ["VC_STAGE"]=<Y>, ["COMBINE_VARIANT_STAGE"]=<Y>, ["JOINT_GENOTYPING_STAGE"]=<Y>, ["PAIRED"]=<1>, ["ALIGNERTOOL"]=<BWAMEM>, ["MARKDUPLICATESTOOL"]=<PICARD>, ["CHRNAMES"]=<1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:23> ==> 3 MPI ranks, 1 node

5. ["REALIGN"]=<NO>, ["SPLIT"]=<YES>, ["EXIT_ON_ERROR"]=<YES>, ["ALIGN_STAGE"]=<Y>, ["DEDUP_SORT_STAGE"]=<Y>, ["CHR_SPLIT_STAGE"]=<Y>, ["VC_STAGE"]=<Y>, ["COMBINE_VARIANT_STAGE"]=<Y>, ["JOINT_GENOTYPING_STAGE"]=<Y>, ["PAIRED"]=<1>, ["ALIGNERTOOL"]=<NOVOALIGN>, ["MARKDUPLICATESTOOL"]=<PICARD>, ["CHRNAMES"]=<1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:23> ==> 3 MPI ranks, 1 node


6. ["REALIGN"]=<NO>, ["SPLIT"]=<YES>, ["EXIT_ON_ERROR"]=<YES>, ["ALIGN_STAGE"]=<Y>, ["DEDUP_SORT_STAGE"]=<Y>, ["CHR_SPLIT_STAGE"]=<Y>, ["VC_STAGE"]=<Y>, ["COMBINE_VARIANT_STAGE"]=<Y>, ["JOINT_GENOTYPING_STAGE"]=<Y>, ["PAIRED"]=<1>, ["ALIGNERTOOL"]=<BWAMEM>, ["MARKDUPLICATESTOOL"]=<PICARD>, ["CHRNAMES"]=<1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22> ==> 3 MPI ranks, 1 node

Note: the `runfile`, `sampleinfo` and `settings.sh` files provided in the `GIAB` directory are an example configuration options to a pipeline run as in the above scenarios (used in test scenario 4 above). The corresponding files: `turbine-directory.txt` is created automatically to point to where the swift-t log files are; whereas the file `log.swift.GIAB.run` contains a more elaborate log containing job submission to Slurm details in addition to swift-t pipeline compilation details. In both senarios 4 and 5 above, note the user error of specifying chromosome names up to 23. Since splitting by chromosome was enabled, and there is no chromosome 23 in the human reference genome, the corresponding bam file was not there leading to terminating the workflow (the complete log is not provided here). This error was corrected in the last scenario 6 as can be seen above.

## Miscelleneous:

The `settings.sh` file contains the variables needed to allocate resources for a Swift/T pipeline run. I was experimenting with these variables; namely: PROCS, NODES and PPN to see how to translate Swift/T options to Slurm options (See: https://github.com/swift-lang/swift-t/issues/134)

