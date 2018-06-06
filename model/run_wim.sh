#! /bin/sh

if [ $# -lt 2 ]
then
   echo "not enough arguments"
   echo "run_wim.sh cfg1 cfg2"
   echo or
   echo "run_wim.sh cfg1 cfg2 [file to source environment variables from]"
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
   rm -rf bin #make sure any old executable is deleted
   mkdir -p bin
   echo "cp $NEXTSIMDIR/model/$prog $prog"
   cp $NEXTSIMDIR/model/$prog $prog
fi

# extra settings needed for mac
kernel=$(uname -s)
if [ "$kernel" == "Darwin" ]
then
    # mac
    export DYLD_LIBRARY_PATH=$NEXTSIMDIR/lib:$BOOST_LIBDIR
fi

# set the log name - don't want to accidentally overwrite it
logfile=${configs[0]}
logfile=${logfile%.*}
logfile=${logfile}-${configs[1]}
logfile=${logfile%.*}
n=1
lf0=$logfile
while [ -f $logfile.log ]
do
   logfile=${lf0}-${n}
   (( n++ ))
done
logfile=${logfile}.log

# Run the nextsim model coupled with wim
echo "$prog --config-files=${configs[@]} > $logfile"
$prog --config-files=${configs[@]} > $logfile 2>&1

exit 0
# print info on plotting on grid
scr=$WIM2D_PATH/fortran/tools/plot_prog.sh

echo ""
echo "Make plots/movie with:"
echo $scr
$scr
