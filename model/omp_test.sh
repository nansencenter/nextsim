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

# Run with one thread
OMP_NUM_THREADS=1 $execfile --simul.duration=0.25 --simul.output_per_day=4 --config-file=$cfgfile #|| exit 5

# Save the old final files
mkdir -p ../matlab/omp_test
mv ../matlab/*_1000.??? ../matlab/omp_test/

# Run with 8 threads
OMP_NUM_THREADS=8 $execfile --simul.duration=0.25 --simul.ouput_per_day=4 --config-file=$cfgfile #|| exit 5

# Test for diff
results=0
diff -qs ../matlab/omp_test/field_1000.bin ../matlab/field_1000.bin
(( results=results+$? ))
diff -qs ../matlab/omp_test/mesh_1000.bin ../matlab/mesh_1000.bin
(( results=results+$? ))

if [ $results -gt 0 ]
then
	echo "Different threads DO NOT produce bit-wise identical results"
else
	echo "Different threads produce bit-wise identical results"
fi


