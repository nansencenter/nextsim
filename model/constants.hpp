/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   constants.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Thu Sep 10 10:45:25 2015
 */

#ifndef __Constants_HPP
#define __Constants_HPP 1

/* Physical constants for the ice as well as not-so-constant values that we
 * keep constant anyway */

namespace physical
{
/* Heat capacity of ice (for Winton model, excluding internal melt) [J/K/kg] */
const double C = 2100;

/* Minimum ice concentration allowed [0 1] */
const double cmin = 0.001;

/* Speciffic heat of air [J/K/kg] */
const double cpa = 1004.64;

/* Speciffic heat of sea water [J/K/kg] */
const double cpw = 4186.84;

/* Emisivity of ice */
const double eps = 0.996;

/* Gravitational acceleration [m/s^2] */
const double g = 9.8;

/* Minimum ice thickness allowed [m] */
const double hmin = 0.01;

/* Heat conductivity of ice [W/K/m] */
const double ki = 2.0334;

/* Heat conductivity of snow [W/K/m] */
const double ks = 0.3096;

/* Latent heat of fusion [J/Kg] */
const double Lf = 334e3;

/* Latent heat of evaporation at 0degC [J/Kg] */
const double Lv0 = 2500.79e3;

/* Proportionalty constant between salinity and freezing temperature of sea
 * water [C] */
const double mu = 0.055;

/* Gas constant for dry air [J/kg/K] */
const double Ra = 286.9;

/* Density of ice (same as in NEMO-LIM) */
const double rhoi = 917.;

/* Density of fresh water */
const double rhofw = 1000.;

/* Density of ocean water */
const double rhow = 1025.;

/* Density of snow (same as in NEMO-LIM) */
const double rhos = 330.;

/* Salinity of ice */
const double si = 5.;

/* Stephan-Boltzmann constant [W/m2/K4] */
const double sigma_sb = 5.67E-8;

/* Freezing temperature of water in Kelvin [K] */
const double tfrwK = 273.15;

/* van Karman constant */
const double vonKarman = 0.4;
}

#endif
