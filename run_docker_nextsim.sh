#!/bin/bash
docker run --rm -it \
    --security-opt seccomp=unconfined \
    -v /home/antonk/nextsim/nextsim:/nextsim \
    -v /Data/sim/data/:/Data/sim/data/ \
    -v /Data/nextsimf/data/:/Data/nextsimf/data/ \
    -v /Data/nextsimf/src/nextsim-cxx/nextsim-develop-new/data/:/simdata/data \
    -v /home/antonk/nextsim/nextsimf_forecasts/:/simforecast/ \
    nextsim /nextsim/test.cfg 7
