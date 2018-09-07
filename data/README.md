# README #

### Introduction ###
To access forcing data, the nextsim model expects a directory called
$NEXTSIM_DATA_DIR to have 2 subdirectories:
* netcdf_data_links:
  - links (or copies) to the netcdf files needed for the simulations
    (no subfolders)
* other_data_links:
  - links (or copies) to other files (eg drifter text files) needed for the simulations
    (no subfolders)

### How do I get set up? ###
To use a directory (including this one) as $NEXTSIM_DATA_DIR
cd [path to the directory]
nextsim/data/scripts/populate_with_data_links.sh [dir1]
nextsim/data/scripts/populate_with_forecast_data_links.sh [dir2]
export NEXTSIM_DATA_DIR=[path to the directory]

where dir1 is the dir with most of the forcing files (eg /Data/sim/data on johansen),
and dir2 is the dir with the ocean forecast forcing files (eg /Data/nextsimf/data on johansen).
