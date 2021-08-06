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
#include "timer.hpp"

namespace po = boost::program_options;

namespace Nextsim {
class Physics1D {
public:
	Physics1D();

	void setFE(FiniteElement & fe);
	void setVariablesFromMap(po::variables_map&);

	// Calculate the physics for a single element
	void thermo(int dt, int i);
private:
	FiniteElement& fe;
	po::variables_map& vm;
	Timer& timer;

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
