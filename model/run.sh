#! /bin/bash

if [ $# -ne 1 ]
then
   echo "not enough arguments"
   echo "run.sh [config file name]"
   exit
fi
config=$1

# record changes from last git commit:
# file gets moved from current dir to "output_directory" inside nextsim code
# NB want file paths relative to $NEXTSIMDIR
P=`pwd`
cd $NEXTSIMDIR
git diff > $P/git_changes.txt
cd $P

prog=bin/nextsim.exec
if [ `pwd` != $NEXTSIMDIR/model ]
then
   # make a local copy of executable
   # (so can recompile code and run somewhere else)
   rm -rf bin #make sure any old executable is deleted
   mkdir -p bin
   echo "cp $NEXTSIMDIR/model/$prog $prog"
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
$prog --config-files=$config


# Run the CPU profiler (google perftools)
# nbt=`echo ${NEXTSIM_BUILD_TYPE,,}`
nbt=`echo $NEXTSIM_BUILD_TYPE  | tr '[:upper:]' '[:lower:]'`
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
