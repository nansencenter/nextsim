#! /bin/bash

# script to run code to find memory errors
# run in a screen (be careful of possibly changed enviroment variables)
# or use nohup as valgrind is very slow
# NB code needs to compiled with env variable NEXTSIM_BUILD_TYPE=debug

export MPIWRAP_DEBUG=warn

if [ $# -lt 3 ]
then
        echo "Usage: $0 config1.cfg config1.cfg num_cpus"
        echo "Usage: $0 config1.cfg config1.cfg num_cpus envfile"
        exit 1
fi

configs=($1 $2)
ncpu=$3
if [ $# -ge 4 ]
then
   echo "source $4"
   source $4
fi

# record changes from last git commit:
# file gets moved from current dir to "output_directory" inside nextsim code
# NB want file paths relative to $NEXTSIMDIR
P=`pwd`
cd $NEXTSIMDIR
git diff > $P/git_changes.txt
cd $P

# valgrind options
vopts+=("--log-file=vlgrnd.log")
vopts+=("--leak-check=full")  # see details of leaked memory
vopts+=("--track-origins=yes") # see where uninitialised values come from

# nextsim options
nsopts=()
nsopts+=("--config-files=${configs[@]}")
nsopts+=("--debugging.maxiteration=3") # just run nextsim for 1 time step
nsopts+=("--nextwim.couplingfreq=1") # just run wim for as short a time as possible

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
mpirun $mpi_opts -np $ncpu valgrind ${vopts[@]} $prog -mat_mumps_icntl_14 60 ${nsopts[@]} 2>&1 | tee wimdebug.log
