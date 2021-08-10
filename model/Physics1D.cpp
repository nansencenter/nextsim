/*
 * Physics1D.cpp
 *
 *  Created on: 6 Aug 2021
 *      Author: Tim Spain, <timothy.spain@nersc.no>
 */

#include "Physics1D.hpp"

namespace Nextsim {
Physics1D::Physics1D():
fe(0), vm(0), timer(0)
{
} // Physics1D::Physics1D()

void Physics1D::setFE(FiniteElement& fe) {
	this->fe = fe;

	// Copy across settings
	settings.iceCategoryType = fe.M_ice_cat_type;
	settings.oceanType = fe.M_ocean_type;
	settings.thermoType = fe.M_thermo_type;
	settings.weldingType = fe.M_welding_type;
	settings.oceanHeatFluxScheme = fe.M_Qio_type;
	settings.freezingPointType = fe.M_freezingpoint_type;

	// Reference the FiniteElement timer instance
	timer = fe.M_timer;

} // void Physics1D::setFE(FiniteElement&)

Physics1D::Settings::Settings(FiniteElement& fe, po::variables_map& vm) :
		iceCategoryType(fe.M_ice_cat_type),
		oceanType(fe.M_ocean_type),
		thermoType(fe.M_thermo_type),
		weldingType(fe.M_welding_type),
		oceanHeatFluxScheme(fe.M_Qio_type),
		freezingPointType(fe.M_freezingpoint_type),
		drag_ocean_t(vm["thermo.drag_ocean_t"].as<double>()),
		drag_ocean_q(vm["thermo.drag_ocean_q"].as<double>()),
		ocean_albedo(vm["thermo.albedoW"].as<double>())
{ }


void Physics1D::thermo(int dt, int i) {
    //! 2) Calculate atmospheric fluxes

    //! Calculate the ocean-atmosphere fluxes
	timer.tick("fluxes");
	timer.tick("ow_fluxes");
} // void Physics1D::thermo(int dt, int i)

void Physics1D::OWBulkFluxes(double& Qow, // scalar versions of the arguments
		double& Qlw,         // of the FiniteElement version
		double& Qsw,
		double& Qlh,
		double& Qsh,
		double& evap,
		double& tau,         // variables read directly from
		double sst,          // FiniteElement arrays
		double t_air,
		double mslp,
		double Qsw_in
		) {
#ifdef AEROBULK
    if ( M_ocean_bulk_formula != aerobulk::algorithm::OTHER ) {
    	// aerobulk
    } else {
#else
    {
#endif
    	// not aerobulk
    }
    }
}

} // namespace Nextsim
