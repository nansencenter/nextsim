#!/bin/bash
CONFIG=${1}     # name of config file
CPUS=${2-'4'}   # number of CPUS
export NEXTSIMDIR=/nextsim
export SIMDATADIR=/nextsim
export SIMFORECASTDIR=/simforecast
mpirun --allow-run-as-root -np $CPUS /usr/local/bin/nextsim.exec --config-files=$CONFIG
