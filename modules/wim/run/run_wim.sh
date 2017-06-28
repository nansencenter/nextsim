kernel=`uname -s`
if [ "$kernel" == "Linux" ]
then
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NEXTSIMDIR/lib:$BOOST_DIR/lib
fi
bin/WIM2d_single_call_cpp.exec --config-file=wim.cfg

echo " "
echo "To plot results, run"
echo "\$WIM2D_PATH/fortran/tools/plot_prog.sh"
$WIM2D_PATH/fortran/tools/plot_prog.sh
