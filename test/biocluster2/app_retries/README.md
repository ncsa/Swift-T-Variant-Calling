This README serves to report back on efforts to enable robust  Swift/T app retries. The testing environment is Biocluster2, which is slurm based. Below are the directions for tests thus far, and lessons learnt:

# Slurm resources allocation

SLURM jobs and resources are scheduled according to the file /code/scripts/submit/slurm/turbine-slurm.sh.m4; where resources allocations (at least NODES and PPN) are obtained from these environemnt variables as per the  /turbine/code/scripts/submit/run-init.zsh  

The logic implement from line 180 (below), is not in alignment with how many NODES are reserved- not sure why. See below
for more detail:

```
# Round NODE up for extra processes
	export NODES=$(( PROCS/PPN ))
	(( PROCS % PPN )) && (( NODES++ )) || true
	export TURBINE_WORKERS=$(( PROCS - ADLB_SERVERS ))
	declare NODES PROCS PPN 
	```


## Resources folder:

Here, I play with values as follows:

|Retry_i|PPN | PROCS| NODES|Remarks|
|-------|----|------|------|------|
|retry_1|3|12| -| NODES are mysteriously calculated(2)|
|retry_2|3|12|5| For some reason, the output.txt file suggests 15 MPI ranks|
|retry_3|3|-|5| I get in total 15 MPI ranks (the question is therefore what is the effect of PROCS parameter)|
|retry_4|-|12|5| There are 5 MPI ranks and 5 nodes (meaning that PPN=1 was automatically set).|
|retry_5|-|-|5| I get the same resources as in retry_4; so this is a good indication that PROCS is more or less a moot variable?|


# Important note: 
Retries applies to leaf function failures, not to programatic errors. For example, exceeding array limits would cause the program to fail, but no retry apply in this case (this is obviious, but it helps to write out)

# Location directive 

Probably, the main lesson learnt here is that a location object is a data structure with 3 components: `rank`, `strictness`, and `accuracy`; where:
`rank`: is the mpi-rank
`strictness`: is the location specification constratint, which can be either `HARD` or `SOFT`
`Accuracy`: is the choice to run a given task on either a specific `RANK` or `NODE`


