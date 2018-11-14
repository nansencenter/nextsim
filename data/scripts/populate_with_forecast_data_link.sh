#! /bin/bash
#
# Link all the ice-ocean forecast files (TOPAZ) needed by the model into the current directory
#
# The links are made from ROOT_DATA_DIR (1st command line input to script)
# which has all the needed data (eg /Data/nextsimf/data on johansen)"
# 
# Usage: `basename $thisfile` [ROOT_DATA_DIR]"
# Where ROOT_DATA_DIR is the directory where the data will be linked from"

thisfile=$0

if [ $# -eq 0 ]
then
   echo "Usage: `basename $thisfile` [ROOT_DATA_DIR]"
   echo "Where ROOT_DATA_DIR is the directory where the data will be linked from"
   echo "(eg /Data/nextsimf/data on johansen)"
   exit 1
fi

ROOT_DATA_DIR="$1"
if [ ! -d "$ROOT_DATA_DIR" ]
then
   echo "Input directory << $ROOT_DATA_DIR >> does not exist"
   exit 1
fi

# link the data
ln -s $ROOT_DATA_DIR/TOPAZ4RC_daily .
