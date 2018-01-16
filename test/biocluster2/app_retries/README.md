This README serves to report back on efforts to enable robust  Swift/T app retries. The testing environment is Biocluster2, which is slurm based. Below are the directions for tests thus far, and lessons learnt:

# Slurm resources allocation:

The main parameters to allocate resources for a Swift/T run and their respective usage is as follows:
`PPN`: Processes-per-node        ==> `#SBATCH --ntasks-per-node=getenv(PPN)`: [default = 1](https://github.com/swift-lang/swift-t/blob/master/turbine/code/scripts/submit/run-init.zsh#L85)
`PROCS`: Number of MPI processes ==> *Not directly translated to slurm parameter*: [default = 0?](https://github.com/swift-lang/swift-t/blob/master/turbine/code/scripts/submit/run-init.zsh#L111)
`NODES`:                         ==> `#SBATCH --nodes=getenv(NODES)`

These are defined in the script [run-init.zsh  ](https://github.com/swift-lang/swift-t/blob/master/turbine/code/scripts/submit/run-init.zsh), before propogating to the respective slurm submission scripts [turbine-slurm-run.zsh](https://github.com/swift-lang/swift-t/blob/master/turbine/code/scripts/submit/slurm/turbine-slurm-run.zsh) and [turbine-slurm.sh.m4](https://github.com/swift-lang/swift-t/blob/master/turbine/code/scripts/submit/slurm/turbine-slurm.sh.m4i#L46)

The key peace of code that govrens these variables is below:

```
export NODES=$(( PROCS/PPN ))
(( PROCS % PPN )) && (( NODES++ )) || true
export TURBINE_WORKERS=$(( PROCS - ADLB_SERVERS ))
declare NODES PROCS PPN 
```
Below, are a summary of how these have been set, and what resources where actually reserved. Actual files from eact test are within the [Resources directory here](https://github.com/ncsa/Swift-T-Variant-Calling/tree/master/test/biocluster2/app_retries/Resources), including the email_report recieved from the scheduler upon the start of a given job.


|Folder_name|PPN| PROCS| NODES|Resulting_Resources_Allocated (total_mpi_ranks (=PPN\*NODES) and nodes (=NODEs))|
|-------|----|------|------|------|
|retry_2|3|12|5| 15 MPI ranks, on 5 nodes (logical as per this table heading definition)|
|retry_1|3|12|-| 6\* MPI ranks, on 2\* nodes |
|retry_3|3|-|5|  15 MPI ranks, on 5 nodes (logical as per this table heading definition)|
|retry_4|-|12|5| 5 MPI ranks, on 5 nodes(logical because default PPN=1)|
|retry_5|-|-|5| 5 MPI ranks, on 5 nodes (logical because default PPN=1)|
|retry_6|-|5|-| 2\* MPI ranks, on 2\* nodes\* |
|retry_7|5|-|-| 6\* MPI ranks, on 2\* nodes+|


The star (\*) denotes strange value.

**These tests make me question what the effect of the `PROCS` parameter (defined as the Number of MPI processes) is. **

# Location tests: 

Probably, the main lesson learnt here is that a location object is a data structure with 3 components: `rank`, `strictness`, and `accuracy`; where:
`rank`: is the mpi-rank
`strictness`: is the location specification constratint, which can be either `HARD` or `SOFT`
`Accuracy`: is the choice to run a given task on either a specific `RANK` or `NODE`

# Other observations:

- Random number generation doesn't really seem random. I noticed that when calling `rand_worker` funciton, I would always get a worker with rank 0. This seemed strage, but then I simplified the `tests_location/1631*`  scricpt, and 


# App retries cases: 
Retries applies to leaf function failures, not to programatic errors. For example, exceeding array limits would cause the program to fail, but no retry apply in this case (this is obviious, but it helps to write out)
