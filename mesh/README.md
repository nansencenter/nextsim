# README #

### Introduction ###
The nextsim model expects a directory called $NEXTSIM_MESH_DIR to have:
- links (or copies) to the initial meshes
- links to .mpp files for stereographic projections
- links to .geo files (TODO remove this from model)

### How do I get set up? ###
Choose a directory in which to put $NEXTSIM_MESH_DIR:
cd [path to the directory]
nextsim/mesh/scripts/populate_with_mesh_links.sh [root mesh directory]
export NEXTSIM_MESH_DIR=[current directory]/mesh_links

where [root mesh directory] is the directory with copies of the meshes (eg /Data/sim/data/mesh on johansen)
