#! /bin/bash

if [ -z "$LOCALFORECASTDATADIR" ]; then 
	echo "LOCALFORECASTDATADIR is unset, we try to link withe the data from /Data/nextsim/data"; 

	TMP_NEXTSIM_FORECAST_DATADIR=/Data/nextsimf/data
	if [ -d "$TMP_NEXTSIM_FORECAST_DATADIR" ]; then
		echo "Found a directory $TMP_NEXTSIM_FORECAST_DATADIR"
  		# Control will enter here if $DIRECTORY exists.
		export NEXTSIM_FORECAST_DATADIR=$TMP_NEXTSIM_FORECAST_DATADIR
	fi

else 
	echo "data are linked from '$LOCALFORECASTDATADIR'"; 
	
	TMP_NEXTSIM_DATADIR=$LOCALFORECASTDATADIR
	if [ -d "$TMP_NEXTSIM_FORECAST_DATADIR" ]; then
  		# Control will enter here if $DIRECTORY exists.
		export NEXTSIM_FORECAST_DATADIR=$TMP_NEXTSIM_FORECAST_DATADIR
	fi
fi


if [ -z "$NEXTSIM_FORECAST_DATADIR" ]; then 
	echo "Not able to define LOCALFORECASTDATADIR, some data may not be linked"; 
else
	# link the data
	ln -s $NEXTSIM_FORECAST_DATADIR/TOPAZ4RC_daily/*.nc .
fi
