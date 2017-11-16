#! /bin/bash

export DYLD_LIBRARY_PATH=$NEXTSIMDIR/lib:$BOOST_LIBDIR

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

# Create a temporary directory
tmpdir_start=$(mktemp -d)
tmpdir_restart=$(mktemp -d)

# Run from start
$execfile --setup.exporter_precision=double --simul.output_directory=$tmpdir_start --setup.use_restart=false --simul.duration=0.5 --setup.write_restart=true --setup.restart_time_step=0.25 --simul.output_per_day=4 --config-files=$cfgfile || exit 5

# Copy the restart files to ../restart so they can be read by the model
cp $tmpdir_start/restart/* ../restart/

# Run with a restart
$execfile --setup.exporter_precision=double --simul.output_directory=$tmpdir_restart --setup.use_restart=true --setup.step_nb=1 --simul.duration=0.5 --setup.write_restart=false --simul.ouput_per_day=4 --config-files=$cfgfile || exit 6

# Test for diff
results=0
diff -qs $tmpdir_start/field_1000.bin $tmpdir_restart/field_1000.bin
(( results=results+$? ))
diff -qs $tmpdir_start/mesh_1000.bin $tmpdir_restart/mesh_1000.bin
(( results=results+$? ))

if [ $results -gt 0 ]
then
    echo "Restart DOES NOT produce bit-wise identical results"
    echo "Leaving temporary files in place:"
    echo "First run is in: $tmpdir_start"
    echo "Restarted run is in: $tmpdir_restart"
    exit 1
else
    echo "Restart produces bit-wise identical results"
    echo "Removing temporary files"
    rm -rf $tmpdir_start $tmpdir_restart
    exit 0
fi


