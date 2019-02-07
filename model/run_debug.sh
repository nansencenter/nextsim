#! /bin/bash

# script to run code to find memory errors
# run in a screen (be careful of possibly changed enviroment variables)
# or use nohup as valgrind is very slow
# NB code needs to compiled with env variable NEXTSIM_BUILD_TYPE=debug

export MPIWRAP_DEBUG=warn

if [ "$1" = "" -o "$2" = "" ]
then
    echo "Usage: $0 [CONFIG_FILE NUM_CPUS]"
    echo "Or: $0 [CONFIG_FILE NUM_CPUS MUMPS_MEM]"
    echo "where MUMPS_MEM is the memory reserved for the solver (%)"
    exit 1
fi

CONFIG=$1
NCPU=$2
# Memory reserved for solver
MUMPS_MEM=${3-200}

# record changes from last git commit:
# file gets moved from current dir to "output_directory" inside nextsim code
# NB want file paths relative to $NEXTSIMDIR
P=`pwd`
cd $NEXTSIMDIR
git diff > $P/git_changes.txt
cd $P

# valgrind options
vopts[0]="--log-file=vlgrnd.log"
vopts[1]="--leak-check=full"  # see details of leaked memory
vopts[2]="--track-origins=yes" # see where uninitialised values come from

# nextsim options
nsopts[0]="--config-file=$CONFIG"
nsopts[1]="--debugging.maxiteration=1" # just run nextsim for 1 time step

prog=bin/nextsim.exec
if [ `pwd` != $NEXTSIMDIR/model ]
then
   # make a local copy of executable
   # (so can recompile code and run somewhere else)
   mkdir -p bin
   cp $NEXTSIMDIR/model/$prog $prog
fi

# settings required for mac
kernel=$(uname -s)
if [ $kernel == "Darwin" ]
then
    # mac
    export DYLD_LIBRARY_PATH=$NEXTSIMDIR/lib:$BOOST_LIBDIR
    export LD_PRELOAD=/opt/local/lib/valgrind/libmpiwrap-amd64-darwin.so
fi

# settings required for johansen
if [ $(hostname) == "johansen.ad.nersc.no" ]
then
	mpi_opts="--mca pml ob1 --mca btl self,tcp"
	export LD_PRELOAD=/usr/lib64/valgrind/libmpiwrap-amd64-linux.so
fi

# Run the nextsim model
mpirun $mpi_opts -np $NCPU valgrind ${vopts[@]} $prog -mat_mumps_icntl_23 $MUMPS_MEM ${nsopts[@]} 2>&1 | tee simdebug.log
