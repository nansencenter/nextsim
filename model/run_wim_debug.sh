# script to run code to find memory errors
# run in a screen (be careful of possibly changed enviroment variables)
# or use nohup as valgrind is very slow
# NB code needs to compiled with env variable NEXTSIM_BUILD_TYPE=debug

# valgrind options
vopts[0]="--log-file=vlgrnd.log"
vopts[1]="--leak-check=full"  # see details of leaked memory
vopts[2]="--track-origins=yes" # see where uninitialised values come from

# nextsim options
nsopts[0]="--config-files=coupling_wim.cfg wim.cfg"
nsopts[1]="--simul.maxiteration=1" # just run nextsim for 1 time step
nsopts[2]="--nextwim.couplingfreq=1" # run wim for only one time step of nextsim

prog=bin/nextsim.exec
if [ `pwd` != $NEXTSIMDIR/model ]
then
   # make a local copy of executable
   # (so can recompile code and run somewhere else)
   mkdir -p bin
   cp $NEXTSIMDIR/model/$prog $prog
fi

valgrind ${vopts[@]} $prog ${nsopts[@]} > wimdebug.log
