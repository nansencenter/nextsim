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

void Physics1D::setVariablesFromMap(po::variables_map& vm) {
	this->vm = vm;
} // void Physics1D::setVariablesFromMap(po::variables_map&)

void Physics1D::thermo(int dt, int i) {
    //! 2) Calculate atmospheric fluxes

    //! Calculate the ocean-atmosphere fluxes
	timer.tick("fluxes");
	timer.tick("ow_fluxes");
} // void Physics1D::thermo(int dt, int i)

} // namespace Nextsim
