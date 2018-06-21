#! /bin/bash

# script to run code to find memory errors
# run in a screen (be careful of possibly changed enviroment variables)
# or use nohup as valgrind is very slow
# NB code needs to compiled with env variable NEXTSIM_BUILD_TYPE=debug

if [ "$1" ]
then
        config=$1
else
        config=nextsim.cfg
fi

# valgrind options
vopts[0]="--log-file=vlgrnd.log"
vopts[1]="--leak-check=full"  # see details of leaked memory
vopts[2]="--track-origins=yes" # see where uninitialised values come from

# nextsim options
nsopts[0]="--config-file=$config"
nsopts[1]="--debugging.max_iteration=1" # just run nextsim for 1 time step

prog=bin/nextsim.exec
if [ `pwd` != $NEXTSIMDIR/model ]
then
   # make a local copy of executable
   # (so can recompile code and run somewhere else)
   rm -rf bin #make sure any old executable is deleted
   mkdir -p bin
   cp $NEXTSIMDIR/model/$prog $prog
fi

# extra settings required for mac
kernel=$(uname -s)
if [ $kernel == "Darwin" ]
then
    # mac
    export DYLD_LIBRARY_PATH=$NEXTSIMDIR/lib:$BOOST_LIBDIR
fi

# Run the nextsim model
valgrind ${vopts[@]} $prog ${nsopts[@]} > simdebug.log


# Run the CPU profiler (google perftools)
nbt=`echo ${NEXTSIM_BUILD_TYPE,,}`
# this is now lower case
if [ "$nbt" == "debug" ]
then
	echo "============================================"
	echo "CPU profiling starts..."
	CMDSTRPROF="pprof --pdf --functions --focus=run --cum --drop_negative --nodecount=50 bin/nextsim.exec profile.log > profile.pdf"
	echo $CMDSTRPROF
	eval $CMDSTRPROF
	echo "CPU profiling done"
	echo "============================================"
	echo "Run the following command to analyze the CPU profile:"
	echo "open profile.pdf"
fi
