/*
 * Physics1D.cpp
 *
 *  Created on: 6 Aug 2021
 *      Author: Tim Spain, <timothy.spain@nersc.no>
 */

#include "Physics1D.hpp"

namespace Nextsim {
Physics1D::Physics1D():
fe(0)
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
} // void Physics1D::setFE(FiniteElement&)


} // namespace Nextsim
