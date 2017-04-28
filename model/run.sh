#! /bin/bash

if [ "$1" ]
then
        config=$1
else
        config=nextsim.cfg
fi

kernel=$(uname -s)
if [ $kernel == "Darwin" ]
then
    export DYLD_LIBRARY_PATH=$NEXTSIMDIR/lib:$BOOST_DIR/lib
    # Run the nextsim model
    bin/nextsim.exec --config-files=$config
else
    # Run the nextsim model
    bin/nextsim.exec --config-files=$config
fi

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
