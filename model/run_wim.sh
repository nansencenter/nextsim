#! /bin/sh

# Run the nextsim model coupled with wim
bin/nextsim.exec --config-files=coupling_wim.cfg wim.cfg

scr=$WIM2D_PATH/fortran/tools/plot_prog.sh

echo ""
echo "Make plots/movie with:"
echo $scr
$scr
