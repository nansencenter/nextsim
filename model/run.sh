#! /bin/sh

# Run the nextsim model
bin/nextsim.exec --config-files=nextsim.cfg

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
