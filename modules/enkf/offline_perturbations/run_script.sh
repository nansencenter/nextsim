#!/bin/bash
work_path=$(cd `dirname $BASH_SOURCE`;pwd)
# - In ./src/main_pseudo2D.F90,     
#     - Set the length of a sequential perturbations for one member to variable i_step.
#     - Set domain size to xdim = 1024, ydim = 1024,xy_full = xdim*ydim.  It uses FFT, which generates faster when using power of 2.
sed -i "s/^integer, parameter:: length_of_perturbation=.*$/integer, parameter:: length_of_perturbation=$((240*4)), xdim=1024,ydim=1024/g" \
    ${work_path}/src/main_pseudo2D.F90

# compile the code 
cd $work_path/src
make clean; make 
cd $work_path

# create random perturbations in ensemble
ensemble_size=40
rm -rf $work_path/result/*
ulimit -s 2000000 # set sufficient stack  https://stackoverflow.com/questions/66034666/ulimit-stack-size-through-slurm-script
sbatch  $work_path/slurm_script.sh $work_path $ensemble_size 
