#!/bin/bash
docker run --rm -it \
    --cap-add=SYS_PTRACE \
    --security-opt seccomp=unconfined \
    --network=host \
    -e DISPLAY=$DISPLAY \
    -v /home/antonk/nextsim/nextsim:/nextsimdir \
    -v /Data/sim/data/:/Data/sim/data/ \
    -v /Data/nextsimf/data/:/Data/nextsimf/data/ \
    -v /Data/nextsimf/src/nextsim-cxx/nextsim-develop-new/data/:/nextsimdir/data \
    -v /home/antonk/nextsim/nextsimf_forecasts/:/simforecastdir/ \
    $1
