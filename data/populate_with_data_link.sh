#! /bin/bash

if [ -z "$LOCALDATADIR" ]; then 
	echo "LOCALDATADIR is unset, we try to link withe the data from /Data/sim/data or /net/sverdrup-1/vol/sim/data"; 

	TMP_NEXTSIM_DATADIR=/net/sverdrup-1/vol/sim/data
	if [ -d "$TMP_NEXTSIM_DATADIR" ]; then
  		# Control will enter here if $DIRECTORY exists.
		export NEXTSIM_DATADIR=$TMP_NEXTSIM_DATADIR
	fi

	TMP_NEXTSIM_DATADIR=/Data/sim/data
	if [ -d "$TMP_NEXTSIM_DATADIR" ]; then
  		# Control will enter here if $DIRECTORY exists.
		export NEXTSIM_DATADIR=$TMP_NEXTSIM_DATADIR
	fi

else 
	echo "data are linked from '$LOCALDATADIR'"; 
	
	TMP_NEXTSIM_DATADIR=$LOCALDATADIR
	if [ -d "$TMP_NEXTSIM_DATADIR" ]; then
  		# Control will enter here if $DIRECTORY exists.
		export NEXTSIM_DATADIR=$TMP_NEXTSIM_DATADIR
	fi
fi


if [ -z "$NEXTSIM_DATADIR" ]; then 
	echo "Not able to define NEXTSIM_DATADIR, some data may not be linked"; 
else
	echo "We link the data from $NEXTSIM_DATADIR"; 
# link the data
ln -s $NEXTSIM_DATADIR/BATHYMETRY/*.nc .
ln -s $NEXTSIM_DATADIR/TOPAZ4/198910_201512/*.nc .
ln -s $NEXTSIM_DATADIR/ERAI/*.nc .
ln -s $NEXTSIM_DATADIR/ASR_FINAL/*.nc .
ln -s $NEXTSIM_DATADIR/CFSR/*.nc .
ln -s $NEXTSIM_DATADIR/CS2_SMOS_v13/*.nc .
ln -s $NEXTSIM_DATADIR/AMSR2_ice_conc/*.nc .
ln -s $NEXTSIM_DATADIR/NIC_ice_charts/*.nc .
ln -s $NEXTSIM_DATADIR/NIC_weekly_ice_charts/*.nc .
ln -s $NEXTSIM_DATADIR/SIT_data/icesat_filled_10prods/*.nc .
ln -s $NEXTSIM_DATADIR/currents_from_altimeter/*.nc .
ln -s $NEXTSIM_DATADIR/RGPS_ice_drift/trajectories/*.txt .
ln -s $NEXTSIM_DATADIR/WW3arctic_RT/*.nc .
ln -s $NEXTSIM_DATADIR/ECMWF_forecast_arctic/*.nc .

# For the OSISAF drifters
ln -s $NEXTSIM_DATADIR/OSISAF_ice_drift/2017/07/ice_drift_nh_polstere-625_multi-oi_201707101200-201707121200.nc ice_drift_nh_polstere-625_multi-oi.nc

# For the Moorings
ln -s $NEXTSIM_DATADIR/AMSR2_ice_conc/Arc_20170710_res3.125_pyres.nc Arc_res3.125_pyres.nc

for year in {2002..2011}
	do
		ln -s $NEXTSIM_DATADIR/AMSRE_ice_conc/$year/*.nc .
	done
for year in {2005..2017}
	do
		ln -s $NEXTSIM_DATADIR/OSISAF_ice_type/$year/*/*.nc .
	done
for year in {2009..2017}
	do
		ln -s $NEXTSIM_DATADIR/OSISAF_ice_drift/$year/*/*.nc .
	done
for year in {2005..2017}
	do
		TMP_DATADIR=$NEXTSIM_DATADIR/OSISAF_ice_conc/polstere/${year}_nh_polstere
        	
		echo $TMP_DATADIR
		if [ -d "$TMP_DATADIR" ]; then
                # Control will enter here if $DIRECTORY exists.
			ln -s $TMP_DATADIR/*.nc .
        	fi

	done
fi
