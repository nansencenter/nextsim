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
#ifdef AEROBULK
		oceanBulkFormula(fe.M_ocean_bulk_formula),
#else
		oceanBulkFormula(0),
#endif
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
		double sss,
		double t_air,
		double t_cc,
		double mslp,
		double Qsw_in,
		double qsrml,
		double windSpeed
		) {
#ifdef AEROBULK
    if ( settings.oceanBulkFormula != aerobulk::algorithm::OTHER ) {
    	// aerobulk
    	aerobulkWrapper(Qlh, Qsh, evap, tau, sst, sss, t_air, t_cc, mslp, Qsw);
    } else {
    	nonRadiativeFluxes(Qlh, Qsh, evap, tau, sst, sss, t_air, mslp, windSpeed);
    }
#else
	nonRadiativeFluxes(Qlh, Qsh, evap, tau, sst, sss, t_air, mslp, windSpeed);
#endif
	// radiative fluxes
	radiativeFluxes(Qlw, Qsw, sst, t_air, t_cc, Qsw_in);

	// Total outward flux
	Qow = Qlw + Qsh + Qlh + Qsw;

#ifdef OASIS
    if ( M_ocean_type == setup::OceanType::COUPLED )
    	// Correct for the amount of shortwave flux absorbed in the top layer
    	// of the ocean, derived from the ocean model. Equivalent to replacing
    	// Qsw above with Qsw * qsrml.
    	Qow -= (1 - qsrml) * Qsw
#endif
}  // void Physics1D::OWBulkFluxes(...)


void Physics1D::nonRadiativeFluxes(
		double& Qlh, double& Qsh, double& evap, double& tau,
		double sst, double sss, double t_air, double mslp, double windSpeed) {

	// not aerobulk: Heat, momentum and vapour fluxes
	double sphuma = specificHumidityAir(t_air, mslp);
	double sphumw = specificHumidityWater(sst, sss, mslp);
	double rhoair = airDensity(mslp, t_air, sphuma);
	Qsh = sensibleHeatFlux(rhoair, sphuma, windSpeed, sst - t_air);
	evap = evaporation(rhoair, windSpeed, sphumw - sphuma);
	Qlh = evap * latentHeatVaporization(sst);
	tau = oceanDrag(rhoair, windSpeed);
} // void Physics1D::nonRadiativeFluxes(...)

void Physics1D::aerobulkWrapper(double& Qlh,
			double& Qsh,
			double& evap,
			double& tau,
			double sst,
			double sss,
			double t_air,
			double t_cc,
			double mslp,
			double Qsw) {
	double sst_kelvin = sst + physical::tfrwK;
	double t2m_kelvin = t_air + physical::tfrwK;
	double Qsw_in = Qsw;
	double Qlw_in = incomingLongwave(t_air, t_cc);
	double sphuma = specificHumidityAir(sst, mslp);
	double windSpeed = windSpeed();
	// Aerobulk does not support a scalar interface :-(
	// vectors to receive the calculated values
	std::vector<double> v_Qlh(1);
	std::vector<double> v_Qsh(1);
	std::vector<double> v_tau(1);
	std::vector<double> dummy(1);
	std::vector<double> T_s(1);
	std::vector<double> lvap(1);
	// Heights of measurements in metres
	const double HEIGHT_OF_TEMPERATURE = 2.;
	const double HEIGHT_OF_WINDS = 10.;
	// Transform scalars into 1-element std::vectors
#ifdef AEROBULK
	aerobulk::model(settings.oceanBulkFormula,
			HEIGHT_OF_TEMPERATURE, HEIGHT_OF_WIND,
			std::vector<double>({sst_kelvin}),
			std::vector<double>({t2m_kelvin}),
			std::vector<double>({sphuma}),
			std::vector<double>({windSpeed}),
			std::vector<double>({0.}),
			std::vector<double>({mslp}),
			Qlh, Qsh, tau, dummy,
			std::vector<double>({Qsw_in}),
			std::vector<double>({Qlw_in}),
			T_s);
	lvap = aerobulk::lvap(sst);
#endif
	// Transfer data back to the current element with post-processing:
	// Change sign on the fluxes, divide tau with wind speed, and calculate
	// evaporation
	Qlh = -v_Qlh[0];
	Qsh = -v_Qsh[0];
	tau = v_tau[0] / (windSpeed * windSpeed);
	evap = Qlh / lvap[0];

} // void Physics1D::aerobulkWrapper(...)

double Physics1D::incomingLongwave(double t_air_celsius, double t_cc) {
    // S. B. Idso & R. D. Jackson, Thermal radiation from the atmosphere, J. Geophys. Res. 74, 5397-5403, (1969)
	// σT^4{1 - c exp [−d (273 - T)^2]}
	double c = 0.261;
	double d = 7.77e-4;
	// Stefan-Boltzmann
	double lwr = stefanBoltzmanInCelsius(t_air_celsius);
	// Idso and Jackson emittance correction
	lwr *= (1 - c * exp(-d * pow(t_air_celsius, 2)));
	// Cloud cover correction from TODO: reference
	double f = 0.275;
	lwr *= (1 + f * t_cc);
	return lwr;
} // double Physics1D::incomingLongwave(double t_air_centigrade, double t_cc)

double Physics1D::stefanBoltzmanInCelsius(double temperature) {
	// Absolute temperature
	double t = temperature + physical::tfrwK;
	return physical::sigma_sb * t * t * t * t;
} // double Physics1D::stefanBoltzman(double temperature)

double Physics1D::airDensity(double mslp, double t_c, double sphuma) {
	double t_k = t_c + physical::tfrwK;
	double rhoDry = mslp / (physical::Ra_dry * t_k);
	double humidityCorrection = 1 - sphuma * (1 - physical::Ra_vap/physical::Ra_dry);

	return rhoDry * humidityCorrection;
} // double Physics1D::airDensity(...)


// Specific Humidity coefficients for ice
const static double ICE_A = 2.2e-4;
const static double ICE_B = 3.83e-6;
const static double ICE_C = 6.4e-10;

const static double ICE_a = 6.1115e2;
const static double ICE_b = 23.036;
const static double ICE_c = 279.82;
const static double ICE_d = 333.7;

const static double WATER_A = 7.2e-4;
const static double WATER_B = 3.20e-6;
const static double WATER_C = 5.9e-10;

const static double WATER_a = 6.1121e2;
const static double WATER_b = 18.729;
const static double WATER_c = 257.87;
const static double WATER_d = 227.3;

const static double ALPHA = 0.62197;
const static double BETA = 1. - ALPHA;

// Helper functions
double f(double A, double B, double C, double p, double T) {

	const double mbFromPa = 1e-2;
	return 1 + A + p * mbFromPa * ( B + C * T );
}

double f(double p, double T, bool isIce) {
	return (isIce) ? f(ICE_A, ICE_B, ICE_C, p, T) : f(WATER_A, WATER_B, WATER_C, p, T);
}

double est(double a, double b, double c, double d, double T, double S) {
	const double salinity_coefficient = 5.37e-4;
	return a * exp( (b-T/d) * T / (T + c)) * (1 - salinity_coefficient * S);
}

double est(double T, double S, bool isIce) {
	return (isIce) ? est(ICE_a, ICE_b, ICE_c, ICE_d, T, 0.) : est(WATER_a, WATER_b, WATER_c, WATER_d, T, S);
}

double Physics1D::generalHumidity(double T, double S, double p, bool isIce) {

	double f = f(p, T, isIce);
	double est = est(T, S, isIce);
	double specificHumidity = (ALPHA * f * est) / (p - BETA * f * est);

	return specificHumidity;
} // double generalHumidity(double T, double S, double p, bool isIce)

double Physics1D::specificHumidityAir(double t_air, double slp) {
	return generalHumidity(t_air, 0., slp, false);
}

double Physics1D::specificHumidityWater(double sst, double sss, double slp) {
	return generalHumidity(sst, sss, slp, false);
}

double Physics1D::specificHumidityIce(double t_ice, double slp) {
	return generalHumidity(t_ice, 0., slp, true);
}

double Physics1D::dSH_dT(double t_ice, double slp) {

	double f		  = f(slp, t_ice, true);
	double est		  = est(t_ice, 0., true);
    double dfdT       = 2. * ICE_C * ICE_B * t_ice;
    double destdT     = ( ICE_b * ICE_c * ICE_d - t_ice * (2. * ICE_c + t_ice) )/( ICE_d*std::pow(ICE_c+t_ice,2) )*est(t_ice, slp, true);
    double dsphumdT   = ALPHA * slp * ( f*destdT + est*dfdT ) /
    						std::pow(slp - BETA*est*f,2);
    return dsphumdT;
} // double Physics1D::dSH_dT(double t_ice, double slp)

double Physics1D::sensibleHeatFlux(double density, double specificHumidity, double windSpeed, double temperatureDifference) {
	return settings.drag_ocean_t * density * (physical::cpa + specificHumidity * physical::cpv) * windSpeed * temperatureDifference;
} // double Physics1D::sensibleHeatFlux(double rho, double sh, double wind, double t_diff)

double Physics1D::latentHeatVaporization(double sst) {
	return physical::Lv0 +
			sst * (-2.36418e3 +
					sst * (1.58927 +
							sst * 6.14342e-2));
} // double Physics1D::latentHeatVaporization(double sst)

double Physics1D::evaporation(double density, double wind, double sh_diff) {
	return settings.drag_ocean_q * density * wind *sh_diff;
} // double Physics1D::evaporation(double rho, double wind, double sh_diff)

double Physics1D::oceanDrag(double density, double windSpeed) {
	// Drag coefficient from Gill(1982) / Smith (1980)
	double drag_ocean_m = 1e-2 * std::max(1., std::min(2., 0.61 + 0.063 * windSpeed));
	return density * drag_ocean_m;
} // double Physics1D::oceanDrag(double density, double windSpeed)

void Physics1D::radiativeFluxes(double& Qlw, double& Qsw,
		double sst, double t_air, double t_cc, double Qsw_in) {
	// Shortwave flux. Positive is outward.
	Qsw = -Qsw_in * (1. - settings.ocean_albedo);
	// Emitted longwave flux. Positive is outward.
	double Qlw_out = stefanBoltzmanInCelsius(t_air);
	// Net outward longwave flux.
	Qlw = Qlw_out - incomingLongwave(t_air, t_cc);

} // void Physics1D::radiativeFluxes(...)

} // namespace Nextsim
