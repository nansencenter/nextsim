#! /bin/bash
#
# Make a directory called "mesh_links" in the current directory
# and link all the mesh-related files needed by the model into this directory
# These files are the meshes themselves (*.msh files), the .mpp files defining different
# stereographic projections, and the .geo files (TODO not really needed
# - remove the ref's in the code)
#
# The .msh files are linked from ROOT_MESH_DIR (1st command line input to script)
# which has all the needed data (eg /Data/sim/data/mesh on johansen)"
# 
# Usage: `basename $thisfile` [ROOT_MESH_DIR]"
# Where ROOT_MESH_DIR is the directory where the data will be linked from"

CWD=`pwd`

if [ $# -eq 0 ]
then
   echo "Usage: `basename $0` [ROOT_MESH_DIR]"
   echo "Where ROOT_MESH_DIR is the directory where the data will be linked from"
   echo "(eg /Data/sim/data/mesh on johansen)"
   exit 1
fi

ROOT_MESH_DIR="$1"
if [ ! -d "$ROOT_MESH_DIR" ]
then
   echo "Input directory << $ROOT_MESH_DIR >> does not exist"
   exit 1
fi

# link the meshes
echo "We link the meshes from $ROOT_MESH_DIR"; 
ln -s $ROOT_MESH_DIR/unref/*.msh .
ln -s $ROOT_MESH_DIR/split/*.msh .

# link WIM grids
ln -s $ROOT_MESH_DIR/wim_grids/*.a .
ln -s $ROOT_MESH_DIR/wim_grids/*.b .
