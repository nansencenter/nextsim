#! /bin/bash

if [ -z "$LOCALDATADIR" ]; then 
	echo "LOCALDATADIR is unset, we try to link with the data from /Data/sim/data or /net/sverdrup-1/vol/sim/data"; 

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
	echo "Not able to define NEXTSIM_DATADIR, some mesh may not be linked"; 
else
	# link the meshes
	echo "We link the meshes from $NEXTSIM_DATADIR"; 
	ln -s $NEXTSIM_DATADIR/mesh/*.msh .
	ln -s $NEXTSIM_DATADIR/mesh/Arctic_meshes_for_neXtSIM/*.msh .
fi
