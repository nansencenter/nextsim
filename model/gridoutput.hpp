/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gridoutput.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Thu Aug  4 09:47:27 CEST 2016
 */

#ifndef __GridOutput_H
#define __GridOutput_H 1

#include <environment.hpp>
#include <assert.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/format.hpp>
#include <gmshmesh.hpp>
#include <gmshmeshseq.hpp>
#include <Bamgx.h>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromMeshToGridx.h>
#include <ConservativeRemapping.hpp>
#include <BamgTriangulatex.h>
#include <netcdf>
#include <debug.hpp>

/**
 * @class GridOutput
 * @brief Manage output of gridded data
 *
 * @see
 *
 */

namespace Nextsim
{
class GridOutput
{
public:

    //////////////////////////////////////////////////////////////////////
    // Misc. enumerations and type definitions
    //////////////////////////////////////////////////////////////////////

    enum fileLength
    {
        // How much data to put into each output file
        inf     = 0, // Never start a new file
        daily   = 1,
        weekly  = 2,
        monthly = 3,
        yearly  = 4
    };

    enum variableKind
    {
        nodal      =  0,
        elemental  =  1
    };

    enum interpMethod
    {
        meshToMesh = 0,
        conservative = 1
    };

    typedef struct Grid
    {
        Grid()
            : loaded(false), defined(false)
        {}

        Grid(std::string file, std::string lat, std::string lon, bool transp=false)
            : gridFile(file), latName(lat), lonName(lon), transpose(transp),
                thetaName(""), interp_method(interpMethod::meshToMesh),
                loaded(false), defined(true)
        {}

        Grid(std::string file, std::string lat, std::string lon, std::string theta, bool transp=false)
            : gridFile(file), latName(lat), lonName(lon), thetaName(theta), transpose(transp),
                interp_method(interpMethod::meshToMesh),
                loaded(false), defined(true)
        {}

        Grid(std::string file, std::string lat, std::string lon, std::string theta, interpMethod method, bool transp=false)
            : gridFile(file), latName(lat), lonName(lon), thetaName(theta), transpose(transp), interp_method(method),
                loaded(false), defined(true)
        {}

        interpMethod interp_method;

        std::string gridFile;
        std::string latName;
        std::string lonName;
        std::string thetaName;

        bool loaded;
        bool defined;
        bool transpose;

        std::string dimNameX;
        std::string dimNameY;
        std::vector<double> gridLAT;
        std::vector<double> gridLON;
        std::vector<double> gridX;
        std::vector<double> gridY;

        std::vector<double> gridTheta;
        std::vector<double> gridCornerLat;
        std::vector<double> gridCornerLon;
        std::vector<double> gridCornerX;
        std::vector<double> gridCornerY;

    } Grid;

    //////////////////////////////////////////////////////////////////////
    // Output variable definitions
    //////////////////////////////////////////////////////////////////////
    /*
     * We need to define the variable ID and the associated name,
     * longName, stdName and units. Try to follow the CF table of
     * standard names whenever possible
     * (http://cfconventions.org/standard-names.html).
     */

    enum variableID
    {
        // Land-sea mask
        lsm =  0,

        // Prognostic variables
        conc         =  1,
        thick        =  2,
        snow         =  3,
        damage       =  4,
        ridge_ratio  =  5,
        VT_x         =  6,
        VT_y         =  7,
        tsurf        =  8,
        sst          =  9,
        sss          = 10,
        tsurf_ice    = 11,
        t1           = 12,
        t2           = 13,
        h_young      = 14,
        hs_young     = 15,
        conc_young   = 16,
        fyi_fraction = 17,
        age_d        = 18,
        age          = 19,
        conc_upd     = 20,
        sigma_11     = 21,
        sigma_22     = 22,
        sigma_12     = 23,

        // Diagnostic variables
        Qa           = 100,
        Qsw          = 101,
        Qlw          = 102,
        Qsh          = 103,
        Qlh          = 104,
        Qo           = 105,
        delS         = 106,
        rain         = 107,
        evap         = 108,
        d_crit       = 109,
        vice_melt    = 110,
        del_hi       = 111,
        del_hi_young = 112,
        newice       = 113,
        snow2ice     = 114,
        mlt_top      = 115,
        mlt_bot      = 116,
        del_vi_young = 117,
        sigma_n      = 118,
        sigma_s      = 119,
        divergence   = 120,
        albedo       = 121,

        // Forcing variables
        tair     = 200,
        sphuma   = 201,
        mixrat   = 202,
        d2m      = 203,
        mslp     = 204,
        Qsw_in   = 205,
        Qlw_in   = 206,
        tcc      = 207,
        snowfall = 208,
        precip   = 209,
        snowfr   = 210,
        wind_x   = 211,
        wind_y   = 212,
        wspeed   = 213,
        tau_ax   = 214,
        tau_ay   = 215,

        // WIM variables
        dmax        = 300,
        dmean       = 301,

        // Coupling variables not already covered elsewhere
        taux       = 901,
        tauy       = 902,
        taumod     = 903,
        fwflux     = 904,
        QNoSw      = 905,
        QSwOcean   = 906,
        saltflux   = 907,
        fwflux_ice = 908,
        tauwix     = 909,
        tauwiy     = 910,

        // Non-output variables - all negative
        proc_mask = -1,
        ice_mask  = -2
    };


    typedef struct Variable
    {
        Variable(){}

        Variable(variableID id)
            : Variable(id, false)
        {}

        Variable(variableID id, bool masked)
            : varID(id), mask(masked), data_mesh(), data_grid()
        {
            switch (varID)
            {
                // Land-sea mask
                case (variableID::lsm):
                    name     = "lsm";
                    longName = "Land Sea Mask";
                    stdName  = "land_sea_mask";
                    Units    = "1";
                    cell_methods = "area: point";
                    break;

                // Prognostic variables
                case (variableID::conc):
                    name     = "sic";
                    longName = "Sea Ice Concentration";
                    stdName  = "sea_ice_area_fraction";
                    Units    = "1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::thick):
                    name     = "sit";
                    longName = "Sea Ice Thickness";
                    stdName  = "sea_ice_thickness";
                    Units    = "m";
                    cell_methods = "area: mean";
                    break;
                case (variableID::damage):
                    name     = "damage";
                    longName = "Sea Ice Damage";
                    stdName  = "sea_ice_damage";
                    Units    = "1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::ridge_ratio):
                    name     = "ridge_ratio";
                    longName = "Sea Ice Ridge Ratio";
                    stdName  = "sea_ice_ridge_ratio";
                    Units    = "1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::snow):
                    name     = "snt";
                    longName = "Surface Snow Thickness";
                    stdName  = "surface_snow_thickness";
                    Units    = "m";
                    cell_methods = "area: mean";
                    break;
                case (variableID::VT_x):
                    name     = "siu";
                    longName = "Sea Ice X Velocity";
                    stdName  = "sea_ice_x_velocity";
                    Units    = "m s-1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::VT_y):
                    name     = "siv";
                    longName = "Sea Ice Y Velocity";
                    stdName  = "sea_ice_y_velocity";
                    Units    = "m s-1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::tsurf):
                    name     = "ts";
                    longName = "Surface Temperature";
                    stdName  = "surface_temperature";
                    Units    = "degree_Celsius";
                    cell_methods = "area: mean";
                    // CF cannonical units are K, but we can use C also
                    break;
                case (variableID::sst):
                    name     = "sst";
                    longName = "Sea Surface Temperature";
                    stdName  = "sea_surface_temperature";
                    Units    = "degree_Celsius";
                    cell_methods = "area: mean";
                    // CF cannonical units are K, but we can use C also
                    break;
                case (variableID::sss):
                    name     = "sss";
                    longName = "Sea Surface Salinity";
                    stdName  = "sea_surface_salinity";
                    Units    = "1e-3";
                    cell_methods = "area: mean";
                    break;
                case (variableID::tsurf_ice):
                    name     = "tsi";
                    longName = "Sea Ice Surface Temperature";
                    stdName  = "sea_ice_surface_temperature";
                    Units    = "degree_Celsius";
                    cell_methods = "area: mean";
                    // CF cannonical units are K, but we can use C also
                    break;
                case (variableID::t1):
                    name     = "t1";
                    longName = "Ice Temperature 1";
                    stdName  = "ice_temperature_1";
                    Units    = "degree_Celsius";
                    cell_methods = "area: mean";
                    // CF cannonical units are K, but we can use C also
                    break;
                case (variableID::t2):
                    name     = "t2";
                    longName = "Ice Temperature 2";
                    stdName  = "ice_temperature_2";
                    Units    = "degree_Celsius";
                    cell_methods = "area: mean";
                    // CF cannonical units are K, but we can use C also
                    break;
                case (variableID::conc_young):
                    name     = "sic_young";
                    longName = "Young Ice Concentration";
                    stdName  = "young_ice_area_fraction";
                    Units    = "1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::h_young):
                    name     = "sit_young";
                    longName = "Young Ice Thickness";
                    stdName  = "young_ice_thickness";
                    Units    = "m";
                    cell_methods = "area: mean";
                    break;
                case (variableID::hs_young):
                    name     = "snt_young";
                    longName = "Surface Snow Thickness on young ice";
                    stdName  = "surface_snow_thickness_on_young_ice";
                    Units    = "m";
                    cell_methods = "area: mean";
                    break;
                case (variableID::fyi_fraction):
                    name     = "fyi_fraction";
                    longName = "First Year Ice Fraction";
                    stdName  = "fyi_fraction";
                    Units    = "1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::age_d):
                    name     = "sia_det";
                    longName = "Detectable sea ice age";
                    stdName  = "det_sea_ice_age";
                    Units    = "s";
                    cell_methods = "area: mean";
                    break;
                case (variableID::age):
                    name     = "sia";
                    longName = "Sea ice age";
                    stdName  = "aea_ice_age";
                    Units    = "s";
                    cell_methods = "area: mean";
                    break;
                case (variableID::conc_upd):
                    name     = "conc_upd";
                    longName = "conc_upd";
                    stdName  = "conc_upd";
                    Units    = "1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::sigma_11):
                    name     = "sigma_11";
                    longName = "Stress tensor 11";
                    stdName  = "stress_tensor_11";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;
                case (variableID::sigma_22):
                    name     = "sigma_22";
                    longName = "Stress tensor 22";
                    stdName  = "stress_tensor_22";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;
                case (variableID::sigma_12):
                    name     = "sigma_12";
                    longName = "Stress tensor 12";
                    stdName  = "stress_tensor_12";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;

                // Diagnostic variables
                case (variableID::Qa):
                    name     = "hfs";
                    longName = "Surface Upward Heat Flux In Air";
                    stdName  = "surface_upward_heat_flux_in_air";
                    Units    = "W m-2";
                    cell_methods = "area: mean";
                    break;
                case (variableID::Qsw):
                    name     = "rss";
                    longName = "Surface Net Upward Shortwave Flux";
                    stdName  = "surface_net_upward_shortwave_flux";
                    Units    = "W m-2";
                    cell_methods = "area: mean";
                    break;
                case (variableID::Qlw):
                    name     = "rls";
                    longName = "Surface Net Upward Longwave Flux";
                    stdName  = "surface_net_upward_longwave_flux";
                    Units    = "W m-2";
                    cell_methods = "area: mean";
                    break;
                case (variableID::Qsh):
                    name     = "hfss";
                    longName = "Surface Upward Sensible Heat Flux";
                    stdName  = "surface_upward_sensible_heat_flux";
                    Units    = "W m-2";
                    cell_methods = "area: mean";
                    break;
                case (variableID::Qlh):
                    name     = "hfsl";
                    longName = "Surface Upward Latent Heat Flux";
                    stdName  = "surface_upward_latent_heat_flux";
                    Units    = "W m-2";
                    cell_methods = "area: mean";
                    break;
                case (variableID::Qo):
                    name     = "hfos";
                    longName = "Surface Upward Heatflux In Ocean";
                    stdName  = "surface_upward_heatflux_in_ocean";
                    Units    = "W m-2";
                    cell_methods = "area: mean";
                    break;
                case (variableID::rain):
                    name     = "rain";
                    longName = "Surface Ocean Liquid Precipitation";
                    stdName  = "surface_ocean_liquid_precipitation";
                    Units    = "kg m-2 s-1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::delS):
                    name     = "sfo";
                    longName = "Downward Saltflux In Ocean";
                    stdName  = "downward_slatflux_in_ocean";
                    Units    = "g m-2 day-1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::evap):
                    name     = "evap";
                    longName = "Surface Ocean Evaporation";
                    stdName  = "surface_ocean_evaporation";
                    Units    = "kg m-2 s-1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::d_crit):
                    name     = "d_crit";
                    longName = "Distance_To_Yield_Criterion";
                    stdName  = "distance_to_yield_criterion";
                    Units    = "1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::del_vi_young):
                    name     = "del_vi_young";
                    longName = "Young Ice Volume Melted or Formed per Day per Surface Area";
                    stdName  = "young_ice_volume_melted_or_formed_per_day_per_surface_area";
                    Units    = " m/day";
                    cell_methods = "area: mean";
                    break;
                case (variableID::vice_melt):
                    name     = "vice_melt";
                    longName = "Ice Volume Melted or Formed per Day per Surface Area";
                    stdName  = "ice_volume_melted_or_formed_per_day_per_surface_area";
                    Units    = " m/day";
                    cell_methods = "area: mean";
                    break;
                case (variableID::del_hi):
                    name     = "del_hi";
                    longName = "Growth-melt rate of (thick) ice";
                    stdName  = "growth_melt_rate_of_thick_ice";
                    Units    = " m/day";
                    cell_methods = "area: mean";
                    break;
                case (variableID::del_hi_young):
                    name     = "del_hi_young";
                    longName = "Growth-melt rate of young ice";
                    stdName  = "growth_melt_rate_of_young_ice";
                    Units    = " m/day";
                    cell_methods = "area: mean";
                    break;
                case (variableID::newice):
                    name     = "newice";
                    longName = "Ice formed in open water by supercooling";
                    stdName  = "ice_formed_in_open_water_by_supercooling";
                    Units    = " m/day";
                    cell_methods = "area: mean";
                    break;
                case (variableID::mlt_bot):
                    name     = "mlt_bot";
                    longName = "Ice melted at bottom";
                    stdName  = "ice_melted_at_bottom";
                    Units    = " m/day";
                    cell_methods = "area: mean";
                    break;
                case (variableID::mlt_top):
                    name     = "mlt_top";
                    longName = "Ice melted at top";
                    stdName  = "ice_melted_at_top";
                    Units    = " m/day";
                    cell_methods = "area: mean";
                    break;
                case (variableID::snow2ice):
                    name     = "snow2ice";
                    longName = "Ice formed from snow by flooding";
                    stdName  = "ice_formed_from_snow_by_flooding";
                    Units    = " m/day";
                    cell_methods = "area: mean";
                    break;
                case (variableID::albedo):
                    name     = "albedo";
                    longName = "Surface albedo";
                    stdName  = "surface_albedo";
                    Units    = "";
                    cell_methods = "area: mean";
                    break;
                case (variableID::sigma_n):
                    name     = "sigma_n";
                    longName = "Normal internal stress";
                    stdName  = "normal_internal_stress";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;
                case (variableID::sigma_s):
                    name     = "sigma_s";
                    longName = "Shear internal stress";
                    stdName  = "shear_internal_stress";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;
                case (variableID::tau_ax):
                    name     = "tau_ax";
                    longName = "Eastward Stress at Ice Surface";
                    stdName  = "eastward_stress_at_ice_surface";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;
                case (variableID::tau_ay):
                    name     = "tau_ay";
                    longName = "Northward Stress at Ice Surface";
                    stdName  = "northward_stress_at_ice_surface";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;
                case (variableID::divergence):
                    name     = "divergence";
                    longName = "Ice Velocity Divergence";
                    stdName  = "ice_velocity_divergence";
                    Units    = "1/s";
                    cell_methods = "area: mean";
                    break;

                // Coupling variables
                case (variableID::taux):
                    name     = "taux";
                    longName = "Eastward Stress at Ocean Surface";
                    stdName  = "eastward_stress_at_ocean_surface";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;
                case (variableID::tauy):
                    name     = "tauy";
                    longName = "Northward Stress at Ocean Surface";
                    stdName  = "northward_stress_at_ocean_surface";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;
                case (variableID::taumod):
                    name     = "taumod";
                    longName = "Downward Stress Magnitude at Ocean Surface";
                    stdName  = "downward_stress_magnitude_at_ocean_surface";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;
                case (variableID::fwflux):
                    name     = "fwflux";
                    longName = "Surface Net Downward Ocean Freshwater Flux";
                    stdName  = "surface_net_downward_ocean_freshwater_flux";
                    Units    = "kg m-2 s-1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::fwflux_ice):
                    name     = "fwflux_ice";
                    longName = "Surface Ice Originating Downward Ocean Freshwater Flux";
                    stdName  = "surface_ice_originating_downward_ocean_freshwater_flux";
                    Units    = "kg m-2 s-1";
                    cell_methods = "area: mean";
                    break;
                case (variableID::QNoSw):
                    name     = "rsnos";
                    longName = "Surface Net Downward Nonsolar Heatflux";
                    stdName  = "surface_net_downward_nonsolar_heatflux";
                    Units    = "W m-2";
                    cell_methods = "area: mean";
                    break;
                case (variableID::QSwOcean):
                    name     = "rsso";
                    longName = "Ocean Net Downward Short Wave Flux";
                    stdName  = "ocean_net_downward_short_wave_flux";
                    Units    = "W m-2";
                    cell_methods = "area: mean";
                    break;
                case (variableID::saltflux):
                    name     = "sfi";
                    longName = "Downward Sea Ice Basal Salt Flux";
                    stdName  = "downward_sea_ice_basal_salt_flux";
                    Units    = "W m-2";
                    cell_methods = "area: mean";
                    break;
                case (variableID::tauwix):
                    name     = "tauwix";
                    longName = "Eastward Stress waves on ice";
                    stdName  = "eastward_stress_waves_on_ice";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;
                case (variableID::tauwiy):
                    name     = "tauwiy";
                    longName = "Northward Stress waves on ice";
                    stdName  = "northward_stress_waves_on_ice";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;

                //WIM variables
                case (variableID::dmax):
                    name     = "dmax";
                    longName = "Maximum floe size";
                    stdName  = "maximum_floe_size";
                    Units    = "m";
                    cell_methods = "area: mean where sea_ice";
                    break;
                
                case (variableID::dmean):
                    name     = "dmean";
                    longName = "Mean floe size";
                    stdName  = "mean_floe_size";
                    Units    = "m";
                    cell_methods = "area: mean where sea_ice";
                    break;

                //forcing variables
                case (variableID::tair):
                    name     = "t2m";
                    longName = "2 metre air temperature";
                    stdName  = "2_metre_air_temperature";
                    Units    = "C";
                    cell_methods = "area: mean";
                    break;

                case (variableID::sphuma):
                    name     = "hus";
                    longName = "specific humidity";
                    stdName  = "specific_humidity";
                    Units    = "kg/kg";
                    cell_methods = "area: mean";
                    break;

                case (variableID::mixrat):
                    name     = "mixrat";
                    longName = "humidity mixing ratio";
                    stdName  = "humidity_mixing_ratio";
                    Units    = "1";
                    cell_methods = "area: mean";
                    break;

                case (variableID::d2m):
                    name     = "d2m";
                    longName = "dew point temperature";
                    stdName  = "dew_point_temperature";
                    Units = "C";
                    cell_methods = "area: mean";
                    break;

                case (variableID::mslp):
                    name     = "psl";
                    longName = "pressure at sea level";
                    stdName  = "pressure_at_sea_level";
                    Units    = "Pa";
                    cell_methods = "area: mean";
                    break;

                case (variableID::Qsw_in):
                    name     = "ssrd";
                    longName = "downward shortwave radiation flux";
                    stdName  = "surface_downwelling_shortwave_flux_in_air";
                    Units    = "W/m^2";
                    cell_methods = "area: mean";
                    break;

                case (variableID::Qlw_in):
                    name     = "strd";
                    longName = "downward thermal radiation flux";
                    stdName  = "surface_downwelling_longwave_flux_in_air";
                    Units    = "W/m^2";
                    cell_methods = "area: mean";
                    break;

                case (variableID::tcc):
                    name     = "tcc";
                    longName = "total cloud cover";
                    stdName  = "cloud_area_fraction";
                    Units    = "1";
                    cell_methods = "area: mean";
                    break;

                case (variableID::snowfall):
                    name     = "sf";
                    longName = "snowfall rate";
                    stdName  = "snowfall_rate";
                    Units    = "kg/m^2/s";
                    cell_methods = "area: mean";
                    break;

                case (variableID::snowfr):
                    name     = "snowfr";
                    longName = "fraction of precipitation that is snow";
                    stdName  = "snow_fraction_of_precipitation";
                    Units    = "1";
                    cell_methods = "area: mean";
                    break;

                case (variableID::precip):
                    name     = "tp";
                    longName = "total precipitation rate";
                    stdName  = "total_precipitation_rate";
                    Units    = "kg/m^2/s";
                    cell_methods = "area: mean";
                    break;

                case (variableID::wind_x):
                    name     = "wndx";
                    longName = "Wind X velocity";
                    stdName  = "wind_x_velocity";
                    Units    = "m/s";
                    cell_methods = "area: mean";
                    break;

                case (variableID::wind_y):
                    name     = "wndy";
                    longName = "Wind Y velocity";
                    stdName  = "wind_y_velocity";
                    Units    = "m/s";
                    cell_methods = "area: mean";
                    break;

                case (variableID::wspeed):
                    name     = "wspeed";
                    longName = "Wind speed";
                    stdName  = "wind_speed";
                    Units    = "m/s";
                    cell_methods = "area: mean";
                    break;

                // Non-output variables
                case (variableID::proc_mask):
                    name     = "proc_mask";
                    longName = "MPI Processor Mask";
                    stdName  = "mpi_proc_mask";
                    Units    = "1";
                    cell_methods = "";
                    break;

                case (variableID::ice_mask):
                    name     = "ice_mask";
                    longName = "Ice Mask";
                    stdName  = "ice_mask";
                    Units    = "1";
                    cell_methods = "";
                    break;

                default:
                    throw std::logic_error("GridOutput::Grid: variableID not defined: "+std::to_string(varID)+"\n");
            }
        }

        variableID varID;
        std::string name;
        std::string longName;
        std::string stdName;
        std::string Units;
        std::string cell_methods;
        bool mask;

        std::vector<double> data_mesh;
        std::vector<double> data_grid;

#ifdef OASIS
        int cpl_id;
#endif

    } Variable;

    typedef struct Vectorial_Variable
    {
        Vectorial_Variable() {}

        Vectorial_Variable(std::pair<int,int> id)
            : components_Id(id)
        {}

        std::pair<int,int> components_Id;
    } Vectorial_Variable;

    typedef typename GmshMesh::bimap_type bimap_type;

    ///////////////////////////////////////////////////////////////////////
    // Constructors (and destructor)
    ///////////////////////////////////////////////////////////////////////
    GridOutput();

    GridOutput(BamgMesh* bamgmesh, int nb_local_el, int ncols, int nrows, double mooring_spacing, double xmin, double ymin, std::vector<Variable> variables, variableKind kind,
            double averaging_period, bool false_easting);

    GridOutput(BamgMesh* bamgmesh, int nb_local_el, Grid grid, std::vector<Variable> variables, variableKind kind,
            double averaging_period, bool false_easting,
        BamgMesh* bamgmesh_root = NULL,
        bimap_type const & transfer_map = boost::bimaps::bimap<int,int>(),
        Communicator const & comm = Environment::comm());

    GridOutput(BamgMesh* bamgmesh, int nb_local_el, int ncols, int nrows, double mooring_spacing, double xmin, double ymin, std::vector<Variable> variables, variableKind kind, std::vector<Vectorial_Variable> vectorial_variables,
            double averaging_period, bool false_easting);

    GridOutput(BamgMesh* bamgmesh, int nb_local_el, Grid grid, std::vector<Variable> variables, variableKind kind, std::vector<Vectorial_Variable> vectorial_variables,
            double averaging_period, bool false_easting,
        BamgMesh* bamgmesh_root = NULL,
        bimap_type const & transfer_map = boost::bimaps::bimap<int,int>(),
        Communicator const & comm = Environment::comm());

    GridOutput(BamgMesh* bamgmesh, int nb_local_el, int ncols, int nrows, double mooring_spacing, double xmin, double ymin, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables,
            double averaging_period, bool false_easting);

    GridOutput(BamgMesh* bamgmesh, int nb_local_el, Grid grid, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables,
            double averaging_period, bool false_easting,
        BamgMesh* bamgmesh_root = NULL,
        bimap_type const & transfer_map = boost::bimaps::bimap<int,int>(),
        Communicator const & comm = Environment::comm());

    GridOutput(BamgMesh* bamgmesh, int nb_local_el, int ncols, int nrows, double mooring_spacing, double xmin, double ymin, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables,
            double averaging_period, bool false_easting);

    GridOutput(BamgMesh* bamgmesh, int nb_local_el, Grid grid, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables,
            double averaging_period, bool false_easting,
        BamgMesh* bamgmesh_root = NULL,
        bimap_type const & transfer_map = boost::bimaps::bimap<int,int>(),
        Communicator const & comm = Environment::comm());

    ~GridOutput();

    void setLSM(BamgMesh* bamgmesh);

    void updateGridMean(BamgMesh* bamgmesh, int nb_local_el,
            std::vector<double> const& UM);
    void resetGridMean();
    void resetMeshMean(BamgMesh* bamgmesh, bool regrid, int nb_local_el,
            const std::vector<int>& gridP, const std::vector<std::vector<int>>& triangles, const std::vector<std::vector<double>>& weights);
    void resetMeshMean(BamgMesh* bamgmesh,
            bool regrid = false,
            int nb_local_el = 0,
            bimap_type const & transfer_map = boost::bimaps::bimap<int,int>(),
            BamgMesh* bamgmesh_root = NULL);
    std::string initNetCDF(const std::string file_prefix, const fileLength file_length,
            const double current_time, const bool append=false);
    void createProjectionVariable(netCDF::NcFile &dataFile);
    void appendNetCDF(std::string filename, double timestamp);

    // Return a mask
    std::vector<int> const &getMask() const { return M_lsm; }

    int M_ncols;
    int M_nrows;
    double M_mooring_spacing = 0;
    bool M_is_regular_grid;
    double M_averaging_period;
    int M_grid_size;
    bool M_false_easting;
    Grid M_grid;
    std::vector<Variable> M_nodal_variables;
    std::vector<Variable> M_elemental_variables;
    std::vector<Vectorial_Variable> M_vectorial_variables;

    double M_miss_val = -1e+14; // Must be smaller than any expected result

    std::vector<int> const &getGridP() const { return M_gridP; }
    std::vector<std::vector<int>> const &getTriangles() const { return M_triangles; }
    std::vector<std::vector<double>> const &getWeights() const { return M_weights; }

    void info();

private:

    LogLevel M_log_level;
    bool M_log_all;

    double M_xmin;
    double M_ymax;

    std::vector<int> M_lsm;
    bool M_use_lsm;

    int M_ice_mask_indx;

    GridOutput(std::vector<Variable> variables, variableKind kind, double averaging_period, bool false_easting);

    GridOutput(std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, double averaging_period, bool false_easting);

    GridOutput(std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables,
            double averaging_period, bool false_easting);

    void initRegularGrid(BamgMesh* bamgmesh, int nb_local_el, int ncols, int nrows, double mooring_spacing, double xmin, double ymin);

    void initArbitraryGrid(BamgMesh* bamgmesh, int nb_local_el, Grid& grid, Communicator const & comm,
            BamgMesh* bamgmesh_root = NULL,
            bimap_type const & transfer_map = boost::bimaps::bimap<int,int>());

    void initMask();

    void applyLSM();

    void updateGridMeanWorker(BamgMesh* bamgmesh, std::vector<double> const& UM,
            variableKind kind, interpMethod method, std::vector<Variable>& variables,
            double miss_val);

    void rotateVectors(Vectorial_Variable const& vectorial_variable, int nb_var, double* &interp_out, double miss_val);

    std::vector<int> M_gridP;
    std::vector<std::vector<int>> M_triangles;
    std::vector<std::vector<double>> M_weights;

    void setProcMask(BamgMesh* bamgmesh, int nb_local_el,
            std::vector<double> const& UM);
    std::vector<double> M_proc_mask;

    Communicator M_comm;

    void broadcastWeights(std::vector<int>& gridP,
            std::vector<std::vector<int>>& triangles,
            std::vector<std::vector<double>>& weights );

    void stereographicProjection(std::vector<double> const & x, std::vector<double> const & y, std::vector<double> const & z,
        std::vector<double>& ps_x, std::vector<double>& ps_y);
};
} // Nextsim
#endif // __GridOutput_H
