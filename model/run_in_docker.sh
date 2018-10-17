#! /bin/bash
export NEXTSIM_MESH_DIR=/mesh
export NEXTSIM_DATA_DIR=/data
export PATH=$PATH:$NEXTSIMDIR/model/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NEXTSIMDIR/lib
source /root/.nextsimrc

# run bash if parameters are not given
if [ $# -eq 0 ]
then
    /bin/bash
    exit 0
fi

# show help if only one parameter is given
if [ $# -eq 1 ] || [ $# -ge 4 ]
then
    echo "Usage:"
    echo "$0"
    echo "to just run bash, or"
    echo "$0 [CONFIG_FILE] [NUM_CPUS]"
    echo "to run the  model with log output on the screen, or"
    echo "$0 [CONFIG_FILE] [NUM_CPUS] [LOG_FILE]"
    echo "to run the  model with log output saved to LOG_FILE"
    echo "NB LOGFILE should be the path inside docker and directory containing it"
    echo "needs to be mounted with -v for it to be accessible outside docker"
    echo "NB we need NUM_CPUS>= 2."
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

# Run the nextsim model using MPI
RUN_NEXTSIM_COMMAND="mpirun --allow-run-as-root -np $ncpu nextsim.exec -mat_mumps_icntl_14 60 --config-files=$config"
if [ $# -ge 3 ]
then
   # output log to a text file
   logdir=`dirname $3`
   if [ ! -d $logdir ]
   then
      echo "Invalid path to output log file: $3"
      exit 1
   fi
   $RUN_NEXTSIM_COMMAND &> $3
else
   # output log to the screen
   $RUN_NEXTSIM_COMMAND
fi
