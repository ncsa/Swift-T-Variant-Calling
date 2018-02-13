export PROCS=12 #total number of MPI processes
export PPN=3 #Number of processes per node 
export NODES=5
export QUEUE=normal
export WALLTIME=00:30:00

#export TURBINE_LOG=1
export ADLB_DEBUG_RANKS=1

export TEST_ADLB_WORKERS=2
export MAIL_ENABLED=0
export MAIL_ADDRESS=azzaea@gmail.com
export TURBINE_APP_RETRIES=10
export SWIFT_TMP=/home/groups/hpcbio_shared/azza/swift_T_project/results/synthetic_data/tmp_swift
export TURBINE_SRAND=$( date +%s )
