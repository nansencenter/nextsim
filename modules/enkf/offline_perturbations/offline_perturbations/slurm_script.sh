#!/bin/bash 
#SBATCH --account=nn2993k   #ACCOUNT_NUMBER
#SBATCH --job-name=perturbations
#SBATCH --time=0-0:15:0     #WALL_TIME_DAYS-WALL_TIME_HOURS:WALL_TIME_MINUTES:0
#SBATCH --nodes=40           #NUM_NODES
#SBATCH --ntasks-per-node=128 

work_path=$1
ensemble_size=$2
for (( i=1; i<=${ensemble_size}; i++ )); do
    mem_path=${work_path}/result/mem$i
    mkdir -p $mem_path
    cd $mem_path
    cp $work_path/bin/p_pseudo2D . 
    cp $work_path/pseudo2D.nml .
    ./p_pseudo2D
    echo "Finish $i over ${ensemble_size} members"
done