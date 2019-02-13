#! /bin/bash

# script to run code to find memory errors
# run in a screen (be careful of possibly changed enviroment variables)
# or use nohup as valgrind is very slow
# NB code needs to compiled with env variable NEXTSIM_BUILD_TYPE=debug

export MPIWRAP_DEBUG=warn

if [ "$1" = "" -o "$2" = "" ]
then
        echo "Usage: $0 config_file.cfg num_cpus"
        echo "or: $0 config_file.cfg num_cpus maxits"
        echo "where maxits=None (run for simul.duration) or a positive integer (run for maxits timesteps)"
        exit 1
fi

config=$1
ncpu=$2
maxits=${3-"None"}

# record changes from last git commit:
# file gets moved from current dir to "output_directory" inside nextsim code
# NB want file paths relative to $NEXTSIMDIR
P=`pwd`
cd $NEXTSIMDIR
git diff > $P/git_changes.txt
cd $P

# valgrind options
vopts=()
vopts+=("--log-file=vlgrnd.log")
vopts+=("--leak-check=full")  # see details of leaked memory
vopts+=("--track-origins=yes") # see where uninitialised values come from

# nextsim options
nsopts=()
nsopts+=("--config-file=$config")
if [ $maxits != "None" ]
then
   if [ $maxits -gt 0 ] 2>/dev/null 
   then
      # no error from -gt if integer
      nsopts+=("--debugging.maxiteration=$maxits") # just run nextsim for 1 time step
   else
      echo "maxits argument needs to be 'None' or an integer >0"
      exit 1
   fi
fi

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
mpirun $mpi_opts -np $ncpu valgrind ${vopts[@]} $prog -mat_mumps_icntl_14 60 ${nsopts[@]} 2>&1 | tee simdebug.log
