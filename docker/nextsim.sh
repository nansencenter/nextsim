#!/bin/bash
export NEXTSIMDIR=/nextsimdir
export SIMDATADIR=/simdatadir
export SIMFORECASTDIR=/simforecastdir
mpirun --allow-run-as-root -np 4 /nextsim/model/bin/nextsim.exec -mat_mumps_icntl_14 60 --config-files=/nextsimdir/docker/test.cfg
