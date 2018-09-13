#!/bin/bash

VOLUMES=()

# mount some directories as themselves so that various links to data work (read-only)
# NB some data is linked from the nextsim source directory
# eg .mpp files, .geo files
VOLUMES+=("-v /Data/sim/data/:/Data/sim/data/:ro")
VOLUMES+=("-v /Data/nextsimf/data/:/Data/nextsimf/data/:ro")
VOLUMES+=("-v /mnt/10.11.12.231/sim/:/mnt/10.11.12.231/sim/:ro")
VOLUMES+=("-v /mnt/10.11.12.232/nextsimf/:/mnt/10.11.12.232/nextsimf/:ro")
VOLUMES+=("-v /home/timill/Github-Repos/nextsim-cxx/nextsim:/home/timill/Github-Repos/nextsim-cxx/nextsim:ro")

# $NEXTSIM_MESH_DIR (read-only)
VOLUMES+=("-v /home/timill/Github-Repos/nextsim-cxx/nextsim/mesh/mesh_links:/mesh:ro")

# $NEXTSIM_DATA_DIR (read-only)
VOLUMES+=("-v /home/timill/Github-Repos/nextsim-cxx/nextsim/data/data_links:/data:ro")

# where docker can find inputs and outputs (read-write)
# NB could separate input and output if you preferred
VOLUMES+=("-v /home/timill/docker_run/:/docker_io/")


ARGS=()
# arguments for run_in_docker.sh
# - comment this line to just enter the container without running the model
ARGS+=(/docker_io/nextsim_model_options_parallel.cfg 2)

docker run --rm -it \
    --security-opt seccomp=unconfined \
    ${VOLUMES[@]} nextsim ${ARGS[@]}
