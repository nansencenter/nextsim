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


# link the data
ln -s $NEXTSIM_DATADIR/BATHYMETRY/*.nc .
ln -s $NEXTSIM_DATADIR/TOPAZ4/198910_201512/*.nc .
ln -s $NEXTSIM_DATADIR/ERAI/*.nc .
ln -s $NEXTSIM_DATADIR/ASR_FINAL/*.nc .
ln -s $NEXTSIM_DATADIR/CFSR/*.nc .
ln -s $NEXTSIM_DATADIR/CS2_SMOS_thickness/*.nc .
ln -s $NEXTSIM_DATADIR/AMSR2_ice_conc/*.nc .
ln -s $NEXTSIM_DATADIR/SIT_data/icesat_filled_10prods/*.nc .
ln -s $NEXTSIM_DATADIR/currents_from_altimeter/*.nc .
ln -s $NEXTSIM_DATADIR/RGPS_ice_drift/trajectories/*.txt .

for year in {2002..2011}
	do
		ln -s $NEXTSIM_DATADIR/AMSRE_ice_conc/$year/*.nc .
	done
for year in {2005..2017}
	do
		ln -s $NEXTSIM_DATADIR/OSISAF_ice_type/$year/*/*.nc .
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
