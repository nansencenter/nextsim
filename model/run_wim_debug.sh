# script to run code to find memory errors
# run in a screen (be careful of possibly changed enviroment variables)
# or use nohup as valgrind is very slow
# NB code needs to compiled with env variable NEXTSIM_BUILD_TYPE=debug

if [ $# -lt 2 ]
then
   echo "not enough arguments"
   echo "run_wim.sh cfg1 cfg2"
   echo or
   echo "run_wim.sh cfg1 cfg2 [file to source environment variables from]"
   echo or
   echo "run_wim.sh cfg1 cfg2 envfile max_iterations"
   echo or
   echo "run_wim.sh cfg1 cfg2 envfile max_iterations wim_coupling_freq"
   exit 1
else
   configs=($1 $2)
fi

if [ $# -ge 3 ]
then
   echo source $3
   source $3
   echo " "
fi

maxits=2
if [ $# -ge 4 ]
then
   maxits=$4
fi

wimfq=1
if [ $# -ge 5 ]
then
   wimfq=$5
fi

# valgrind options
vopts=()
vopts+=("--log-file=vlgrnd.log")
vopts+=("--leak-check=full")     # see details of leaked memory
vopts+=("--track-origins=yes")   # see where uninitialised values come from

# nextsim options
nsopts=()
nsopts+=("--config-files=${configs[@]}")
nsopts+=("--simul.maxiteration=$maxits")  # just run nextsim for 1 or a few time steps
nsopts+=("--nextwim.couplingfreq=$wimfq")      # run wim for only one time step of nextsim

prog=bin/nextsim.exec
if [ `pwd` != $NEXTSIMDIR/model ]
then
   # make a local copy of executable
   # (so can recompile code and run somewhere else)
   rm -rf bin #make sure any old executable is deleted
   mkdir -p bin
   cp $NEXTSIMDIR/model/$prog $prog
fi

valgrind ${vopts[@]} $prog ${nsopts[@]} > wimdebug.log
