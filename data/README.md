# README #

### Introduction ###
The nextsim model expects a directory called $NEXTSIM_DATA_DIR to have:
* links (or copies) to the netcdf files needed for the simulations
* links (or copies) to other files (eg drifter text files) needed for the simulations

### How do I get set up? ###
Choose a directory in which to put $NEXTSIM_DATA_DIR:
cd [path to the directory]
nextsim/data/scripts/populate_with_data_links.sh [root data directory]
nextsim/data/scripts/populate_with_forecast_data_links.sh [root forcing data directory]
export NEXTSIM_DATA_DIR=[current directory]/data_links

where [root data directory] is the directory with copies of the meshes (eg /Data/sim/data on johansen)
and [root forcing data directory] is the dir with the ice-ocean forecast (TOPAZ)
forcing files (eg /Data/nextsimf/data on johansen).
