/*
 * Physics1D.hpp
 *
 *  Created on: 6 Aug 2021
 *      Author: Tim Spain, <timothy.spain@nersc.no>
 */

#ifndef MODEL_PHYSICS1D_HPP_
#define MODEL_PHYSICS1D_HPP_

#include "finiteelement.hpp"
#include "enums.hpp"

namespace Nextsim {
class Physics1D {
public:
	Physics1D();

	void setFE(FiniteElement & fe);


private:
	FiniteElement& fe;

	class Settings {
	public:
		setup::IceCategoryType iceCategoryType;
		setup::OceanType oceanType;
		setup::ThermoType thermoType;
		setup::WeldingType weldingType;
		setup::OceanHeatfluxScheme oceanHeatFluxScheme;
		setup::FreezingPointType freezingPointType;
	} settings;
}; // class Physics1D
} // namespace Nextsim
#endif /* MODEL_PHYSICS1D_HPP_ */
