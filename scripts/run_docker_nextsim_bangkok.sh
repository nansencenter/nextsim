#!/bin/bash

. neXtSIM.env
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
VOLUMES+=("-v ${FORCLOC}:${FORCREM}:ro")

# $NEXTSIM_MESH_DIR (read-only)
VOLUMES+=("-v ${MSHPATH}:/mesh:ro")

# $NEXTSIM_DATA_DIR (read-only)
VOLUMES+=("-v ${REMLINK}:/data:ro")

# where docker can find inputs and outputs (read-write)
# NB could separate input and output if you preferred
VOLUMES+=("-v ${MEMPATH}:/docker_io")

echo "Ensemble Member: " ${MEMPATH}

ARGS=()
if [ "$RUN_MODEL" != 0 ]
then
   # arguments for run_in_docker.sh
   # - input any argument to script to run the model
   ARGS+=(/docker_io/nextsim.cfg ${NPROC})
fi

echo ${VOLUMES[@]} ${ARGS[@]}

docker run --rm \
    --name=${MEMNAME} \
    --security-opt seccomp=unconfined \
    ${VOLUMES[@]} nextsim ${ARGS[@]}
