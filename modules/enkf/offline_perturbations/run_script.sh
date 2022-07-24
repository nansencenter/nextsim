#!/bin/bash
#!/bin/bash 
#SBATCH --account=nn2993k   #ACCOUNT_NUMBER
#SBATCH --job-name=perturbations
#SBATCH --time=0-1:0:0     #WALL_TIME_DAYS-WALL_TIME_HOURS:WALL_TIME_MINUTES:0
#SBATCH --nodes=40           #NUM_NODES
#SBATCH --ntasks-per-node=128 
#SBATCH --output=slurm.perturbation.%j.log         # Stdout
#SBATCH --error=slurm.perturbation.%j.log          # Stderr
work_path=/cluster/work/users/chengsukun/offline_perturbations
# create random perturbations in ensemble
ensemble_size=40
# # - In ./src/main_pseudo2D.F90,     
# #     - Set the length of a sequential perturbations for one member to variable i_step.
# #     - Set domain size to xdim = 1024, ydim = 1024,xy_full = xdim*ydim.  It uses FFT, which generates faster when using power of 2.
# sed -i "s/^integer, parameter:: length_of_perturbation=.*$/integer, parameter:: length_of_perturbation=$((240*4)), xdim=1024,ydim=1024/g" \
#     ${work_path}/src/main_pseudo2D.F90

# # compile the code 
# [ ! -d bin ] && mkdir bin
# [ ! -d objs ] && mkdir objs
# rm -rf $work_path/result/*
# cd $work_path/src
# make clean; make 

#
rm -rf ${work_path}/result/*
# # option 1: run on login nodes
# for (( i=1; i<=${ensemble_size}; i++ )); do
#     mem_path=${work_path}/result/mem$i
#     mkdir -p $mem_path
#     cd $mem_path
#     cp $work_path/bin/p_pseudo2D . 
#     cp $work_path/pseudo2D.nml .
#     ./p_pseudo2D 
#     echo "Finish $i over ${ensemble_size} members"
# done

# option 2: run on slurm
ulimit -s 2000000 # set sufficient stack  https://stackoverflow.com/questions/66034666/ulimit-stack-size-through-slurm-script
# number_of_nodes=4
# sbatch --nodes=$number_of_nodes --qos=devel $work_path/slurm_script.sh $work_path $ensemble_size $number_of_nodes 

for (( i=1; i<=${ensemble_size}; i++ )); do
    mem_path=${work_path}/result/mem$i
    mkdir -p $mem_path
    cd $mem_path
    cp $work_path/bin/p_pseudo2D . 
    cp $work_path/pseudo2D.nml .
    srun --nodes=1 --mpi=pmi2 -n1 ./p_pseudo2D  & # 2>&1 | tee task.log
    echo " $i over ${ensemble_size} members"
done
wait