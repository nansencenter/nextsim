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
	// Requires the wind speed calculation to be done externally and be passed in.
	// TODO: supersede this with something better
	void OWBulkFluxes(double& Qow, // scalar versions of the arguments
			double& Qlw,         // of the FiniteElement version
			double& Qsw,
			double& Qlh,
			double& Qsh,
			double& evap,
			double& tau,         // variables read directly from
			double sst,          // FiniteElement arrays
			double sss,
			double t_air,
			double t_cc,
			double mslp,
			double Qsw_in,
			double windspeed
			);

	// Argument monster version of a Aerobulk wrapper
	// TODO: supersede this with something better
	void aerobulkWrapper(double& Qlh,
			double& Qsh,
			double& evap,
			double& tau,
			double sst,
			double sss,
			double t_air,
			double t_cc,
			double mslp,
			double Qsw);

	void nonRadiativeFluxes(double& Qlh,
			double& Qsh,
			double& evap,
			double& tau,
			double sst,
			double sss,
			double t_air,
			double mslp,
			double windSpeed);

	// Argument monster version of incomingLongwave
	// TODO: supersede this with something better
	double incomingLongwave(double t_air_centigrade, double t_cc);

	// Air density
	double airDensity(double mslp, double t_air, double sphuma);

	// Argument monster version of the first half of specificHumidity (air)
	// TODO: supersede this with something better
	double specificHumidityAir(double t_air, double slp);

	double specificHumidityWater(double sst, double sss, double slp);
	double specificHumidityIce(double t_ice, double slp);
	double dSH_dT(double t_ice, double slp);

	double sensibleHeatFlux(double density, double specificHumidity, double windSpeed, double temperatureDifference);
	double latentHeatVaporization(double sst);
	double evaporation(double density, double windSpeed, double specificHumidityDifference);

	double oceanDrag(double density, double windSpeed);
private:
	FiniteElement& fe;
	po::variables_map& vm;
	Timer& timer;

	double generalHumidity(double T, double S, double p, bool isIce);

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
#ifdef AEROBULK
		const aerobulk::algorithm oceanBulkFormula;
#else
		const int oceanBulkFormula;
#endif
		const double drag_ocean_t;
		const double drag_ocean_q;
		const double ocean_albedo;

	} settings;
}; // class Physics1D
} // namespace Nextsim
#endif /* MODEL_PHYSICS1D_HPP_ */
