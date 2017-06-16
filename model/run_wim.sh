#! /bin/sh

# record changes from last git commit:
# file gets moved to "output_directory"
P=`pwd`
cd $NEXTSIMDIR
git diff > $P/git_changes.txt
cd $P

# Run the nextsim model coupled with wim
kernel=$(uname -s)
if [ $kernel == "Darwin" ]
then
    # mac
    export DYLD_LIBRARY_PATH=$NEXTSIMDIR/lib:$BOOST_DIR/lib
    $NEXTSIMDIR/model/bin/nextsim.exec --config-files=coupling_wim.cfg wim.cfg
else
    # linux
    $NEXTSIMDIR/model/bin/nextsim.exec --config-files=coupling_wim.cfg wim.cfg
fi

# print info on plotting on grid
scr=$WIM2D_PATH/fortran/tools/plot_prog.sh

echo ""
echo "Make plots/movie with:"
echo $scr
$scr
