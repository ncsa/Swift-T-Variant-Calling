export PROCS=6 #total number of MPI processes
export PPN=2 #Number of processes per node 
export NODES=3
export QUEUE=normal
export TURBINE_LOG=1
export MAIL_ENABLED=1
export MAIL_ADDRESS=azzaea@gmail.com
export WALLTIME=100:00:00
export TURBINE_SBATCH_ARGS="--mem=100g"
#export ADLB_PERF_COUNTERS=1
# the variables below didn't really chip in much!
#export  ADLB_REPORT_LEAKS=1
#export  ADLB_TRACE=true
#export ADLB_DEBUG=true 
export TURBINE_LOG=1
export ADLB_DEBUG_RANKS=1
export SWIFT_TMP=/home/groups/hpcbio_shared/azza/swift_T_project/results/synthetic_data/tmp_swift


