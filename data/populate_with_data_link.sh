#! /bin/sh

export NEXTSIM_DATADIR=/net/sverdrup-1/vol/sim/data

# link the data
ln -s $NEXTSIM_DATADIR/BATHYMETRY/*.nc .
ln -s $NEXTSIM_DATADIR/TOPAZ4/198910_201512/*.nc .
ln -s $NEXTSIM_DATADIR/ERAI/*.nc .
ln -s $NEXTSIM_DATADIR/ASR_FINAL/*.nc .
ln -s $NEXTSIM_DATADIR/AMSR2_ice_conc/*.nc .
for year in {2002..20011..1}
	do
		ln -s $NEXTSIM_DATADIR/AMSRE_ice_conc/year/*.nc .
	done
