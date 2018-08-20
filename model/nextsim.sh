#!/bin/bash
CONFIG=${1}     # name of config file
CPUS=${2-'4'}   # number of CPUS
export NEXTSIMDIR=/src
export SIMDATADIR=/src
export SIMFORECASTDIR=/simforecast
mpirun --allow-run-as-root -np $CPUS /usr/local/bin/nextsim.exec --config-files=$CONFIG
