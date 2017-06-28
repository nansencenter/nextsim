# export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$NEXTSIMDIR/lib # not needed on mac
bin/WIM2d_single_call_cpp.exec --config-file=wim.cfg

echo " "
echo "To plot results, run"
echo "\$WIM2D_PATH/fortran/tools/plot_prog.sh"
$WIM2D_PATH/fortran/tools/plot_prog.sh
