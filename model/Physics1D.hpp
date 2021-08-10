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
#include "constants.hpp"

namespace po = boost::program_options;

namespace Nextsim {

class PhysicsSettings;

class Physics1D {
public:
	Physics1D();

	void setFE(FiniteElement& fe);

	// Calculate the physics for a single element
	void thermo(int dt, int i);

	// Argument monster version of the 1D version of FiniteElement::OWBulkFluxes
	// TODO: supersede this with something better
	void OWBulkFluxes(double& Qow, // scalar versions of the arguments
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
			);

	// Argument monster version of incomingLongwave
	// TODO: supersede this with something better
	double incomingLongwave(double t_air_centigrade, double t_cc);

	// Argument monster version of the first half of specificHumidity (air)
	// TODO: supersede this with something better
	double specificHumidityAir(); // TODO: fill in the argument list

	// Argument monster version of the windSpeedElement
	// TODO: supersede this with something better
	double windSpeed(); // TODO: fill in the argument list

private:
	FiniteElement& fe;
	po::variables_map& vm;
	Timer& timer;

private:
	class Settings {
	public:
		Settings(FiniteElement& fe, po::variables_map& vm);

		double out();

		const setup::IceCategoryType iceCategoryType;
		const setup::OceanType oceanType;
		const setup::ThermoType thermoType;
		const setup::WeldingType weldingType;
		const setup::OceanHeatfluxScheme oceanHeatFluxScheme;
		const setup::FreezingPointType freezingPointType;

		const double drag_ocean_t;
		const double drag_ocean_q;
		const double ocean_albedo;

	} settings;
}; // class Physics1D
} // namespace Nextsim
#endif /* MODEL_PHYSICS1D_HPP_ */
