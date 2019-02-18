#! /bin/bash
export NEXTSIM_MESH_DIR=/mesh
export NEXTSIM_DATA_DIR=/data
export PATH=$PATH:$NEXTSIMDIR/model/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NEXTSIMDIR/lib
source /root/.nextsimrc

cd /docker_io
# run bash if parameters are not given
if [ $# -eq 0 ]
then
    /bin/bash
    exit 0
fi

# show help if only one parameter is given
if [ $# -eq 1 ]
then
    echo "Usage: $0 [CONFIG_FILE NUM_CPUS]"
    echo ""
    echo "Run bash or parallel neXtSIM using MPI"
    echo "    CONFIG_FILE, name of the config file for neXtSIM"
    echo "       NUM_CPUS, number of CPUs. Must be >= 2."
    echo ""
    echo "If CONFIG_FILE NUM_CPUS are not given /bin/bash is executed"
    exit 0
fi

# set config file name
config=$1
ncpu=$2
if [ $ncpu -lt 2 ]
then
   echo "Error: num_cpus cannot be less than 2"
   exit 2
fi

# Run the nextsim model usin MPI
mpirun --allow-run-as-root -np $ncpu nextsim.exec -mat_mumps_icntl_14 60 --config-files=$config
