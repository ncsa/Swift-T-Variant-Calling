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


## Miscelleneous:

The `settings.sh` file contains the variables needed to allocate resources for a Swift/T pipeline run. I was experimenting with these variables; namely: PROCS, NODES and PPN to see how to translate Swift/T options to Slurm options (See: https://github.com/swift-lang/swift-t/issues/134)

