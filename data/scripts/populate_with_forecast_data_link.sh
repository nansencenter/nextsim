#! /bin/bash
#
# Link the ocean forecast files from ROOT_DATA_DIR (1st command line input to script)
# which has all the needed data (eg /Data/nextsimf/data on johansen)"
# 
# Usage: `basename $thisfile` [ROOT_DATA_DIR]"
# Where ROOT_DATA_DIR is the directory where the data will be linked from"

thisfile=$0
thisdir=`dirname $thisfile`
CWD=`pwd`

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
NETCDF_LINKS_DIR=$CWD/netcdf_data_links
mkdir -p $NETCDF_LINKS_DIR
cd $NETCDF_LINKS_DIR
ln -s $ROOT_DATA_DIR/TOPAZ4RC_daily/*.nc .

# lot of links so update the list of files
ls -lh > FileList.txt
