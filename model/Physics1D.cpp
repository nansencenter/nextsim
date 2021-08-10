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
    	double sst_kelvin = sst + physical::tfrwK;
    	double t2m_kelvin = t_air + physical:trfwK;
    	double Qlw_in = incomingLongwave();
    	double sphuma = specificHumidityAir();
    	double windSpeed = windSpeed();
    	// Aerobulk does not support a scalar interface :-(
    } else {
#else
    {
#endif
    	// not aerobulk
    }
} // void Physics1D::OWBulkFluxes(...)

double Physics1D::incomingLongwave(double t_air_centigrade, double t_cc) {
	// Convert temperatures to kelvin
	double t_air_kelvin = t_air_centigrade + physical::tfrwK;
    // S. B. Idso & R. D. Jackson, Thermal radiation from the atmosphere, J. Geophys. Res. 74, 5397-5403, (1969)
	// σT^4{1 - c exp [−d (273 - T)^2]}
	double c = 0.261;
	double d = 7.77e-4;
	// Stefan-Boltzmann
	double lwr = physical::sigma_sb * pow(t_air_kelvin, 4);
	// Idso and Jackson emittance correction
	lwr *= (1 - c * exp(-d * pow(t_air_centigrade, 2)));
	// Cloud cover correction from TODO: reference
	double f = 0.275;
	lwr *= (1 + f * t_cc);
	return lwr;
} // double Physics1D::incomingLongwave(double t_air_centigrade, double t_cc)

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
}

} // namespace Nextsim
