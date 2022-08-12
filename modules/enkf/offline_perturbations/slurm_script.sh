#!/bin/bash
#SBATCH --account=nn2993k  #ACCOUNT_NUMBER
#SBATCH --job-name=perturbations
#SBATCH --time=0-1:0:0    #WALL_TIME_DAYS-WALL_TIME_HOURS:WALL_TIME_MINUTES:0
##SBATCH --nodes=1         # for debug, run members in sequential. 
##SBATCH --qos=devel
#SBATCH --nodes=40         # run members in parallel
#SBATCH --ntasks-per-node=128
#SBATCH --output=slurm-%j.log         # Stdout
#SBATCH --error=slurm-%j.log          # Stderr

work_path=$1
ensemble_size=$2
for (( i=1; i<=${ensemble_size}; i++ )); do
    mem_path=${work_path}/result/mem$i
    mkdir -p $mem_path
    cd $mem_path
    cp $work_path/bin/p_pseudo2D .
    cp $work_path/pseudo2D.nml .
    srun --nodes=1 --mpi=pmi2 -n1 ./p_pseudo2D & #2>&1 | tee task.log & # add '&' for parallel mode or in sequential mode for debug # 
    echo " $i over ${ensemble_size} members"
done
wait