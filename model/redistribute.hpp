/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   thin_ice_redistribute.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Tue Nov 17 11:46:10 2015
 */

#ifndef __ThinIce_HPP
#define __ThinIce_HPP 1

#include <constants.hpp>

namespace Nextsim
{
	using namespace physical;

/* Compute the redistribution of thin ice. */
/* Returns the change in volume and concentration of the thick ice as well as
 * the change in volume of thin ice. This is called by thermo_ow_mex as well as
 * thin_ice_redistribute_mex */

void thin_ice_redistribute(double v_thin, double vs_thin, double newice, double c, double tanalpha, double rtanalpha, double hi_thin_max, \
		double *v_thin_new, double *v_newice, double *del_c_new, double *del_vs)
{
	/* Physical constants */

	/* Inputs:
	 * v_thin:	Volume of thin ice
	 * vs_thin:	Volume of snow on thin ice
	 * newice:	Thickness of newly formed ice
	 * c:		Concentration of thick ice
	 * tanalpha:	simul_in.h_thin_max/simul_in.c_thin_max
	 * rtanalpha:	1/tanalpha
	 * hi_thin_max:	Maximum thickness of thin ice */

	/* Outputs:
	 * v_thin_new:	Volume of thin ice
	 * v_newice:	Volume of newly formed (thick) ice
	 * del_c_new:	Change in concentration of thick ice
	 * del_vs:	Volume of snow transfered from the thin ice to the thick */


	double c_thin, h0, del_c, hs_thin;

	/* -----------------------------------------------------------------------
	 * Compute the redistribution of ice volume from the thin ice class to
	 * the tick ice.
	 * ----------------------------------------------------------------------- */

	/* c_thin needs to change as the volume changes */
	c_thin = fmin(v_thin/hmin, sqrt(2.*v_thin*rtanalpha));
	if ( c_thin > 0. )
	{
		hs_thin = vs_thin/c_thin;
	} else {
		hs_thin = 0.;
	}
	v_thin = v_thin + newice*fmax(0., 1.-c-c_thin);
	c_thin = fmin(v_thin/hmin, sqrt(2.*v_thin*rtanalpha));

	/* First we check if the cell is full of thick ice, in which case we
	 * just add add the thin ice to the thick ice and exit */
	if ( c == 1. )
	{
		*del_c_new  = 0.;
		*v_newice   = v_thin;
		*v_thin_new = 0.;
        *del_vs = vs_thin;
		return;
	}

	/* Two cases: Thin ice fills the cell or not */
	if ( c_thin <= 1.-c )
	{
		h0 = sqrt(2.*v_thin*tanalpha);
	} else {
		h0 = v_thin/(1.-c) + 0.5*(1.-c)*tanalpha;
		c_thin = 1.-c;
	}

	/* Concentration and volum changes for the thick ice */
	/* Two cases again: The lead closes or not */
	del_c      = fmax( 0., h0-hi_thin_max )*rtanalpha;
	if ( del_c <= 1.-c )
	{
		*del_c_new = del_c;
		newice     = 0.5*del_c*(h0+hi_thin_max);
	} else {
		*del_c_new = 1.-c;
		newice     = v_thin;
	}
	*v_newice  = newice;

	/* New volume for the thin ice */
	*v_thin_new = v_thin - *v_newice;

	/* Change in snow volume */
	*del_vs = fmin( vs_thin,del_c*hs_thin);
}

} // Nextsim
#endif
