#!/bin/bash
err_report() {
    echo "Error on line $1"
}
trap 'err_report $LINENO' ERR

# ------------------------------------------
work_path=$(cd `dirname $0`;pwd)
echo "work_path=" $work_path
rm -rf $work_path/result/*
cd $work_path/src
make clean; make 

ensemble_size=1
for (( i=1; i<=${ensemble_size}; i++ )); do
    mem_path=${work_path}/result/mem$i
    mkdir -p $mem_path
    cd $mem_path
    cp $work_path/bin/p_pseudo2D . 
    cp $work_path/pseudo2D.nml .
    ./p_pseudo2D
    echo "Finish $i over ${ensemble_size} members"
done
