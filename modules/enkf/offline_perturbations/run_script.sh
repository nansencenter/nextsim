#!/bin/bash

work_path=$(cd `dirname $0`;pwd)
echo "work_path=" $work_path
[ ! -d bin ] && mkdir bin
[ ! -d objs ] && mkdir objs
rm -rf $work_path/result/*
cd $work_path/src
make clean; make 
#
ensemble_size=40
for (( i=1; i<=${ensemble_size}; i++ )); do
    mem_path=${work_path}/result/mem$i
    mkdir -p $mem_path
    cd $mem_path
    cp $work_path/bin/p_pseudo2D . 
    cp $work_path/pseudo2D.nml .
    ./p_pseudo2D
    echo "Finish $i over ${ensemble_size} members"
done

# ulimit -s 2000000 # error of insufficient stack  https://stackoverflow.com/questions/66034666/ulimit-stack-size-through-slurm-script
# sbatch  $work_path/slurm_script.sh $work_path $ensemble_size 
