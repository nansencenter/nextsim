#! /bin/sh

# record changes from last git commit:
# file gets moved to "output_directory"
P=`pwd`
cd $NEXTSIMDIR
git diff > $P/git_changes.txt
cd $P

prog=bin/nextsim.exec
if [ `pwd` != $NEXTSIMDIR/model ]
then
   # make a local copy of executable
   # (so can recompile code and run somewhere else)
   mkdir -p bin
   cp $NEXTSIMDIR/model/$prog $prog
fi

# extra settings needed for mac
kernel=$(uname -s)
if [ $kernel == "Darwin" ]
then
    # mac
    export DYLD_LIBRARY_PATH=$NEXTSIMDIR/lib:$BOOST_DIR/lib
fi

# Run the nextsim model coupled with wim
$prog --config-files=coupling_wim.cfg wim.cfg

# print info on plotting on grid
scr=$WIM2D_PATH/fortran/tools/plot_prog.sh

echo ""
echo "Make plots/movie with:"
echo $scr
$scr
