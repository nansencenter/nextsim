# README #

### Introduction ###
The nextsim model expects a directory called $NEXTSIM_MESH_DIR to have:
- links (or copies) to the initial meshes
- links to .mpp files for stereographic projections

### How do I get set up? ###
1. Simple option:
   * copy the mesh you need to nextsim/mesh
   * `export NEXTSIM_MESH_DIR=[full path to nextsim/mesh]`
2. Option to make it easier to change region
   * Choose a directory in which to put $NEXTSIM_MESH_DIR:
   * `cd [path to the directory]`
   * Get the mpp files if you are outside nextsim/mesh:
     `ln -s [path to nextsim/mesh]/*.mpp .`
   * link all the meshes into that directory:
     `nextsim/mesh/scripts/populate_with_mesh_links.sh [root mesh directory]`
     where [root mesh directory] is the directory with copies of the meshes (eg /Data/sim/data/mesh on johansen)
   * `export NEXTSIM_MESH_DIR=[current directory]`

