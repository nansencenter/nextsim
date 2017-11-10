#! /bin/bash

if [ "$1" ]
then
        config=$1
else
        config=nextsim.cfg
fi

# record changes from last git commit:
# file gets moved to "output_directory"
P=`pwd`
cd $NEXTSIMDIR
git diff > $P/git_changes.txt
cd $P

prog=bin/nextsim.exec
if [ `pwd` != $NEXTSIMDIR ]
then
   # make a local copy of executable
   # (so can recompile code and run somewhere else)
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
$prog --config-files=$config


# do git diff to record changes to code
# default directory
outdir=$NEXTSIMDIR/matlab
lin=`grep "output_directory" $config`
outdir=${lin#*=}
git diff > $outdir/git_changes.txt


# Run the CPU profiler (google perftools)
if [ \( "$NEXTSIM_BUILD_TYPE" == "DEBUG" \) -o \( "$NEXTSIM_BUILD_TYPE" == "Debug" \) -o \( "$NEXTSIM_BUILD_TYPE" == "debug" \) ]; then
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
