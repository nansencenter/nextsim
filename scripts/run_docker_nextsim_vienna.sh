#!/bin/bash

docker run --rm -it \
    --security-opt seccomp=unconfined \
    -v /Data/sim/data/:/Data/sim/data/:ro \
    -v /Data/nextsimf/data/:/Data/nextsimf/data/:ro \
    -v /mnt/10.11.12.231/sim/:/mnt/10.11.12.231/sim/:ro \
    -v /mnt/10.11.12.232/nextsimf/:/mnt/10.11.12.232/nextsimf/:ro \
    -v /home/timill/Github-Repos/nextsim-cxx/nextsim/mesh:/nextsim_mesh_dir/:ro \
    -v /home/timill/Github-Repos/nextsim-cxx/nextsim/data:/nextsim_data_dir/:ro \
    -v /home/timill/docker_run/:/docker_io/ \
    nextsim /docker_io/nextsim_model_options_parallel.cfg 2
    #nextsim
