#! /bin/bash

if [ "$1" = "" -o "$2" = "" ]
then
   echo "Usage: $0 config_file.cfg num_cpus"
   echo "Usage: $0 config_file.cfg num_cpus envfile"
   echo envfile is a file to be sourced to set some environment variables
   exit 1
fi

config=$1
ncpu=$2

if [ $ncpu -lt 2 ]
then
   echo "Error: num_cpus cannot be less than 2"
   exit 2
fi

if [ $# -ge 3 ]
then
   echo "source $3"
   source $3
fi

# record changes from last git commit:
# file gets moved from current dir to "output_directory" inside nextsim code
# NB want file paths relative to $NEXTSIMDIR
P=`pwd`
cd $NEXTSIMDIR
git diff > $P/git_changes.txt
cd $P

prog=bin/nextsim.exec
if [ `pwd` != "$NEXTSIMDIR/model" ]
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

# bugfix for johansen
if [ $(hostname) == "johansen.ad.nersc.no" ]
then
   opts="--mca pml ob1 --mca btl self,tcp"
fi

# Run the nextsim model
logfile=$(basename $config .cfg).log
mpirun $opts -np $ncpu $prog -mat_mumps_icntl_23 400 --config-files=$config 2>&1 | tee $logfile

if [ "$MAILTO" != "" ]
then
   # send email at finish of run
   # export MAILTO=... to activate
   mail -s "neXtSIM run finished ($config)" -a $logfile $MAILTO
fi

# Run the CPU profiler (google perftools)
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
