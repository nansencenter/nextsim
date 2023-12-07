/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   constants.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Thu Sep 10 10:45:25 2015
 */

#ifndef __Constants_HPP
#define __Constants_HPP 1

//! Physical constants for the ice as well as not-so-constant values that we keep constant anyway

namespace physical
{
    //! Heat capacity of ice (for Winton model, excluding internal melt) [J/K/kg]
    const double C = 2100;

    //! Minimum ice concentration allowed [0 1]
    const double cmin = 1e-12;

    //! Specific heat of air [J/K/kg]
    const double cpa = 1000.5;

    //! Specific heat of water vapour [J/K/kg]
    const double cpv = 1860.;

    //! Specific heat of sea water [J/K/kg]
    const double cpw = 4186.84;

    //! Emissivity of ice
    const double eps = 0.996;

    //! Gravitational acceleration [m/s^2]
    const double g = 9.8;

    //! Minimum ice thickness allowed [m]
    const double hmin = 0.01;

    //! Heat conductivity of ice [W/K/m]
    const double ki = 2.0334;

    //! Latent heat of fusion [J/Kg]
    const double Lf = 333.55e3;

    //! Latent heat of evaporation at 0degC [J/Kg]
    const double Lv0 = 2.5e6;

    //! Proportionality constant between salinity and freezing temperature of sea water [C kg/g]
    const double mu = 0.055;

    //! Gas constant for dry air [J/kg/K]
    const double Ra_dry = 287.058;

    //! Gas constant for water vapour [J/kg/K]
    const double Ra_vap = 461.5;

    //! Density of ice (same as in NEMO-LIM) [kg/m3]
    const double rhoi = 917.;

    //! Density of fresh water [kg/m3]
    const double rhofw = 1000.;

    //! Density of ocean water [kg/m3]
    const double rhow = 1025.;

    //! Density of snow (same as in NEMO-LIM) [kg/m3]
    const double rhos = 330.;

    //! Salinity of ice [g/kg]
    const double si = 5.;

    //! Stephan-Boltzmann constant [W/m2/K4]
    const double sigma_sb = 5.67E-8;

    //! Freezing temperature of water in Kelvin [K]
    const double tfrwK = 273.15;

    //! Von Karman constant
    const double vonKarman = 0.4;

    //! Constant of gravity [m/s2]
    const double gravity = 9.80616;

    //! Rotational velocity [rad/s]
    const double omega = 7.292e-5;

    //! Density of air [kg/m3]
    const double rhoa = 1.22;

    //! Temperature of freezing of sea water [degrees Celsius]
    const double ocean_freezing_temp = -1.8;
}

#endif
