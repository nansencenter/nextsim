#!/bin/bash
export NEXTSIMDIR=/nextsimdir
export SIMDATADIR=/nextsimdir
export SIMFORECASTDIR=/simforecastdir
mpirun --allow-run-as-root -np 5 /nextsim/model/bin/nextsim.exec --config-files=/nextsimdir/docker/test.cfg
