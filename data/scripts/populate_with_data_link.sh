#! /bin/bash
#
# Link the files from ROOT_DATA_DIR (1st command line input to script)
# which has all the needed data (eg /Data/sim/data on johansen)"
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
   echo "(eg /Data/sim/data on johansen)"
   exit 1
fi

ROOT_DATA_DIR="$1"
if [ ! -d "$ROOT_DATA_DIR" ]
then
   echo "Input directory << $ROOT_DATA_DIR >> does not exist"
   exit 1
fi

NETCDF_LINKS_DIR=$CWD/netcdf_data_links
mkdir -p $NETCDF_LINKS_DIR
OTHER_LINKS_DIR=$CWD/other_data_links
mkdir -p $OTHER_LINKS_DIR
echo "We link the data from $ROOT_DATA_DIR"; 

# link the netcdf data
cd $NETCDF_LINKS_DIR
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

# For the OSISAF drifters
ln -s $ROOT_DATA_DIR/OSISAF_ice_drift/2017/07/ice_drift_nh_polstere-625_multi-oi_201707101200-201707121200.nc ice_drift_nh_polstere-625_multi-oi.nc

# For the Moorings
ln -s $ROOT_DATA_DIR/AMSR2_ice_conc/Arc_20170710_res3.125_pyres.nc Arc_res3.125_pyres.nc

for year in {2002..2011}
do
     ln -s $ROOT_DATA_DIR/AMSRE_ice_conc/$year/*.nc .
done
for year in {2005..2017}
do
     ln -s $ROOT_DATA_DIR/OSISAF_ice_type/$year/*/*.nc .
done
for year in {2009..2017}
do
     ln -s $ROOT_DATA_DIR/OSISAF_ice_drift/$year/*/*.nc .
done
for year in {2005..2017}
do
   TMP_DATADIR=$ROOT_DATA_DIR/OSISAF_ice_conc/polstere/${year}_nh_polstere
   
   echo $TMP_DATADIR
   if [ -d "$TMP_DATADIR" ]; then
         # Control will enter here if $DIRECTORY exists.
     ln -s $TMP_DATADIR/*.nc .
   fi
done

# CS2-SMOS ice thickness
# script links weekly files to 7 daily files
$thisdir/CS2_SMOS_pp.sh $ROOT_DATA_DIR/CS2_SMOS_v13/*.nc

# lot of links so update the list of files
ls -lh > FileList.txt

# Other links (just so they don't get buried by all the .nc files)
# - eg drifters
cd $OTHER_LINKS_DIR
ln -s $ROOT_DATA_DIR/IABP/IABP_buoys*.txt .
