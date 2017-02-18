#! /bin/bash

if [ "$1" == "--help" ]
then
	echo "Usage: $0 [config.cfg [executable.exec]]"
	exit 0
fi

# Test wether we get the same results with a different number of threads

# The first input can be the config file - but use nextsim.cfg if the first input is empty
cfgfile=$1
: ${cfgfile:="nextsim.cfg"}

# The second input can be the executable file - but use bin/nextsim.exec if the second input is empty
execfile=$2
: ${execfile:="bin/nextsim.exec"}

# Create temporary directories
tmpdir1=$(mktemp -d)
tmpdir8=$(mktemp -d)

# Run with one thread
OMP_NUM_THREADS=1 $execfile --simul.output_directory=$tmpdir1 --simul.duration=0.25 --simul.output_per_day=4 --config-files=$cfgfile #|| exit 5

# Run with 8 threads
OMP_NUM_THREADS=8 $execfile --simul.output_directory=$tmpdir8 --simul.duration=0.25 --simul.output_per_day=4 --config-files=$cfgfile #|| exit 5

# Test for diff
results=0
diff -qs $tmpdir1/field_1000.bin $tmpdir8/field_1000.bin
(( results=results+$? ))
diff -qs $tmpdir1/mesh_1000.bin $tmpdir8/mesh_1000.bin
(( results=results+$? ))

if [ $results -gt 0 ]
then
    echo "Different threads DO NOT produce bit-wise identical results"
    echo "Leaving temporary files in place:"
    echo "Run using 1 thread is in: $tmpdir1"
    echo "Run using 8 threads is in: $tmpdir8"
    exit 1
else
    echo "Different threads produce bit-wise identical results"
fi


