#! /bin/bash
#
# Link the files from ROOT_MESH_DIR (1st command line input to script)
# which has all the needed data (eg /Data/sim/data/mesh on johansen)"
# 
# Usage: `basename $thisfile` [ROOT_MESH_DIR]"
# Where ROOT_MESH_DIR is the directory where the data will be linked from"

thisfile=$0
thisdir=`dirname $thisfile`
CWD=`pwd`

if [ $# -eq 0 ]
then
   echo "Usage: `basename $thisfile` [ROOT_MESH_DIR]"
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

MESH_LINKS_DIR=$CWD/meshfile_links
mkdir -p $MESH_LINKS_DIR
cd $MESH_LINKS_DIR

# link the meshes
echo "We link the meshes from $ROOT_MESH_DIR"; 
ln -s $ROOT_MESH_DIR/unref/*.msh .
ln -s $ROOT_MESH_DIR/split/*.msh .

# link WIM grids
mkdir -p wim_grids
cd wim_grids
ln -s $ROOT_MESH_DIR/wim_grids/*.* .

# link to .mpp files if they are not present
NEXTSIM_MESH_SRC=`readlink -f $thisdir/..`
cd $CWD
if [ ! -d mpp_files  ]
then
   mkdir -p mpp_files
   cd mpp_files
   ln -s $NEXTSIM_MESH_SRC/mpp_files/* .
fi

# link to .geo files if they are not present
cd $CWD
if [ ! -d geo_files  ]
then
   mkdir -p geo_files
   cd geo_files
   ln -s $NEXTSIM_MESH_SRC/geo_files/* .
fi
