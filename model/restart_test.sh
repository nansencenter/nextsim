#! /bin/bash

if [ "$1" == "--help" ]
then
	echo "Usage: $0 [config.cfg [executable.exec]]"
	exit 0
fi

# Test wether we restart correctly

# The first input can be the config file - but use nextsim.cfg if the first input is empty
cfgfile=$1
: ${cfgfile:="nextsim.cfg"}

# The second input can be the executable file - but use bin/nextsim.exec if the second input is empty
execfile=$2
: ${execfile:="bin/nextsim.exec"}

# Run from start
$execfile --setup.use_restart=false --simul.duration=0.5 --setup.write_restart=true --setup.restart_time_step=0.25 --simul.output_per_day=4 --config-file=$cfgfile #|| exit 5

# Save the old final files
mkdir -p ../matlab/restart_test
mv ../matlab/*_1000.??? ../matlab/restart_test/

# Run with a restart
$execfile --setup.use_restart=true -setup.step_nb=1 --simul.duration=0.5 --setup.write_restart=false --simul.ouput_per_day=4 --config-file=$cfgfile #|| exit 5

# Test for diff
results=0
diff -qs ../matlab/restart_test/field_1000.bin ../matlab/field_1000.bin
(( results=results+$? ))
diff -qs ../matlab/restart_test/mesh_1000.bin ../matlab/mesh_1000.bin
(( results=results+$? ))

if [ $results -gt 0 ]
then
	echo "Restart DOES NOT produce bit-wise identical results"
else
	echo "Restart produces bit-wise identical results"
fi


