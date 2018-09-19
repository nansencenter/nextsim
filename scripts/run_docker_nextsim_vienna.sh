#!/bin/bash
if [ $# -eq 0 ]
then
   echo Usage:
   echo "$0 [run_model]" 
   echo "If run_model==0, just go into bash"
   echo "- otherwise run the model"
   exit 1
fi
RUN_MODEL=$1

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
if [ "$RUN_MODEL" != 0 ]
then
   # arguments for run_in_docker.sh
   # - input any argument to script to run the model
   ARGS+=(/docker_io/nextsim_model_options_parallel.cfg 2)
fi

docker run --rm -it \
    --security-opt seccomp=unconfined \
    ${VOLUMES[@]} nextsim ${ARGS[@]}
