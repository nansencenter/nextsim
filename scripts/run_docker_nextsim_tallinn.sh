#!/bin/bash
docker run --rm -it \
    --security-opt seccomp=unconfined \
    -v /home/antonk/nextsim/nextsim/data:/data \
    -v /home/antonk/nextsim/nextsim/mesh:/mesh \
    -v /Data/sim/data:/Data/sim/data \
    -v /Data/nextsimf/data:/Data/nextsimf/data \
    -v /home/antonk/nextsim/forecasts:/forecasts \
    nextsim /forecasts/test.cfg 7
