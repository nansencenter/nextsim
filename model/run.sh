#! /bin/bash

if [ "$1" = "" -o "$2" = "" ]
then
   echo "Usage: $0 CONFIG_FILE NUM_CPUS"
   echo "Usage: $0 CONFIG_FILE NUM_CPUS MUMPS_MEM"
   echo "where MUMPS_MEM is the memory reserved for the solver (%)"
   exit 1
fi

CONFIG=$1
NCPU=$2
if [ $NCPU -lt 2 ]
then
   echo "Error: input NUM_CPUS cannot be less than 2"
   exit 2
fi

# Memory reserved for solver
MUMPS_MEM=${3-400}

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
logfile=$(basename $CONFIG .cfg).log
mpirun $opts -np $NCPU $prog -mat_mumps_icntl_23 $MUMPS_MEM --config-files=$CONFIG 2>&1 | tee $logfile

if [ "$MAILTO" != "" ]
then
   # send email at finish of run
   # export MAILTO=... to activate
   mail -s "neXtSIM run finished ($CONFIG)" -a $logfile $MAILTO
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
