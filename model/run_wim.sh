#! /bin/sh
# bin/nextwim.exec --config-files=nextsim_wim.cfg wim.cfg
bin/nextsim.exec --config-files=coupling_wim.cfg wim.cfg

scr=$WIM2D_PATH/fortran/tools/plot_prog.sh

echo ""
echo "Make plots/movie with:"
echo $scr
$scr
