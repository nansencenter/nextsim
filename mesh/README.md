# README #

### Introduction ###
To access the initial meshes, .mpp files for stereographic projections, and .geo files,
the nextsim model expects a directory called $NEXTSIM_MESH_DIR to have 3 subdirectories:
* meshfile_links:
  - links (or copies) to the initial meshes
* mpp_files
  - .mpp files for stereographic projections
* geo_files
  - .geo files

### How do I get set up? ###
To use a directory (including this one) as $NEXTSIM_MESH_DIR
cd [path to the directory]
nextsim/mesh/scripts/populate_with_mesh_links.sh [dir1]
export NEXTSIM_MESH_DIR=[path to the directory]

where dir1 is the dir with copies of the meshes (eg /Data/sim/data/mesh on johansen)
