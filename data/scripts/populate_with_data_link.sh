#! /bin/bash
#
# Make a directory called "data_links" in the current directory
# and link all the input-data-related files needed by the model into this directory
# These files are for forcings, initial conditions and topography (*.nc files),
# and text files for the drifters

# Link the files from ROOT_DATA_DIR (1st command line input to script)
# which has all the needed data (eg /Data/sim/data on johansen)"
# 
# Usage: `basename $thisfile` [ROOT_DATA_DIR]"
# Where ROOT_DATA_DIR is the directory where the data will be linked from"

thisfile=`readlink -f $0`
thisdir=`dirname $thisfile`
CWD=`pwd`

if [ $# -eq 0 ]
then
   echo "Usage: `basename $thisfile` [ROOT_DATA_DIR]"
   echo "Where ROOT_DATA_DIR is the directory where the data will be linked from"
   echo "(eg /Data/sim/data on johansen)"
   exit 1
fi

ROOT_DATA_DIR="$1"
if [ ! -d "$ROOT_DATA_DIR" ]
then
   echo "Input directory << $ROOT_DATA_DIR >> does not exist"
   exit 1
fi

DATA_LINKS_DIR=$CWD/data_links
mkdir -p $DATA_LINKS_DIR
cd $DATA_LINKS_DIR
echo "We link the data from $ROOT_DATA_DIR"; 

# link the netcdf data
ln -s $ROOT_DATA_DIR/BATHYMETRY/*.nc .
ln -s $ROOT_DATA_DIR/TOPAZ4/198910_201512/*.nc .
ln -s $ROOT_DATA_DIR/ERAI/*.nc .
ln -s $ROOT_DATA_DIR/ASR_FINAL/*.nc .
ln -s $ROOT_DATA_DIR/CFSR/*.nc .
ln -s $ROOT_DATA_DIR/CS2_SMOS_v13/*.nc .
ln -s $ROOT_DATA_DIR/AMSR2_ice_conc/*.nc .
ln -s $ROOT_DATA_DIR/NIC_ice_charts/*.nc .
ln -s $ROOT_DATA_DIR/NIC_weekly_ice_charts/*.nc .
ln -s $ROOT_DATA_DIR/SIT_data/icesat_filled_10prods/*.nc .
ln -s $ROOT_DATA_DIR/currents_from_altimeter/*.nc .
ln -s $ROOT_DATA_DIR/RGPS_ice_drift/trajectories/*.txt .
ln -s $ROOT_DATA_DIR/WW3arctic_RT/*.nc .
ln -s $ROOT_DATA_DIR/ECMWF_forecast_arctic/*.nc .
ln -s $ROOT_DATA_DIR/ERAI_waves_1deg/*.nc .
ln -s $ROOT_DATA_DIR/WW3arctic_RT/*.nc .

# OSISAF products
ln -s $ROOT_DATA_DIR/OSISAF_ice_conc/polstere/*_nh_polstere/*.nc .
ln -s $ROOT_DATA_DIR/OSISAF_ice_type/*/*/*.nc .
ln -s $ROOT_DATA_DIR/OSISAF_ice_drift/*/*/*.nc .

# - example files needed for OSISAF drifters grids
ln -s $ROOT_DATA_DIR/OSISAF_ice_drift/2017/07/ice_drift_nh_polstere-625_multi-oi_201707101200-201707121200.nc ice_drift_nh_polstere-625_multi-oi.nc
ln -s $ROOT_DATA_DIR/OSISAF_ice_drift/ice_drift*refined*.nc .

# For the Moorings
ln -s $ROOT_DATA_DIR/AMSR2_ice_conc/Arc_20170710_res3.125_pyres.nc Arc_res3.125_pyres.nc
ln -s $ROOT_DATA_DIR/AMSRE_ice_conc/*/*.nc .

# CS2-SMOS ice thickness
# script links weekly files to 7 daily files
$thisdir/CS2_SMOS_pp.sh $ROOT_DATA_DIR/CS2_SMOS_v13/*.nc

# Other links (just so they don't get buried by all the .nc files)
# - eg drifters
ln -s $ROOT_DATA_DIR/IABP/IABP_buoys*.txt .

# Link some other files (possibly for testing) (TODO needed?)
ln -s $thisdir/../misc/* .

# lot of links so update the list of files
ls -lh > FileList.txt
