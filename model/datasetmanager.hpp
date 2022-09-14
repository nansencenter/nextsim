/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   meshhandler.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 17 Aug 2022 07:59:06 CEST
 */

#ifndef __DataSetManager_HPP
#define __DataSetManager_HPP 1

#include <date.hpp>
#include <constants.hpp>
#include <enums.hpp>
#include <optionhandler.hpp>
#include <dataset.hpp>
#include <externaldata.hpp>
#include <timer.hpp>
#include <meshhandler.hpp>

namespace Nextsim
{

    class DataSetManager : public MeshHandler
    {
        // Types
        protected:
            typedef DataSet Dataset;
            typedef ExternalData external_data;
            typedef typename std::vector<external_data*> external_data_vec ;
            typedef typename std::vector<Dataset*> Dataset_vec ;

            typedef boost::ptr_vector<external_data> externaldata_ptr_vector;

        // Methods
        public:
            DataSetManager(Timer::timer* timer)
                : vm(Environment::vm())
            {
                this->initOptAndParam();
                M_timer = timer;
            }

            void initExternalData();

        protected:
            void setReloadDatasets();

            void checkReloadDatasets(external_data_vec const& ext_data_vec,
                            double const CRtime, std::vector<double> &RX, std::vector<double> &RY);
            void checkReloadMainDatasets(double const CRtime);

        private:
            void initOptAndParam();
            void initDatasets();
            void forcingAtmosphere();
            void forcingOcean();
            void initBathymetry();
#ifdef OASIS
            void forcingWaves();
#endif

        protected:
            setup::AtmosphereType M_atmosphere_type;
            setup::OceanType M_ocean_type;

            // Datasets
            Dataset M_atmosphere_nodes_dataset;
            Dataset M_atmosphere_elements_dataset;
            Dataset M_atmosphere_bis_elements_dataset;
            Dataset M_ocean_nodes_dataset;
            Dataset M_ocean_elements_dataset;
            Dataset M_bathymetry_elements_dataset;
#ifdef OASIS
            Dataset M_wave_nodes_dataset;
            Dataset M_wave_elements_dataset;
#endif

            Dataset M_ice_topaz_elements_dataset;
            Dataset M_ice_icesat_elements_dataset;
            Dataset M_ice_piomas_elements_dataset;
            Dataset M_ice_amsre_elements_dataset;
            Dataset M_ice_osisaf_elements_dataset;
            Dataset M_ice_osisaf_type_elements_dataset;
            Dataset M_ice_amsr2_elements_dataset;
            Dataset M_ice_nic_elements_dataset;
            Dataset M_ice_nic_weekly_elements_dataset;
            Dataset M_ice_cs2_smos_elements_dataset;
            Dataset M_ice_smos_elements_dataset;

            // Datasets for nesting from outer domain with coarse resolution
            Dataset M_nesting_ocean_elements_dataset;
            Dataset M_nesting_nodes_dataset;
            Dataset M_nesting_ice_elements_dataset;
            Dataset M_nesting_distance_nodes_dataset;
            Dataset M_nesting_distance_elements_dataset;
            Dataset M_nesting_dynamics_elements_dataset;

            Dataset_vec M_datasets_regrid;

            // Thermodynamic and dynamic forcing
            // Atmosphere
            external_data M_wind;         // Surface wind [m/s]
            external_data M_tair;         // 2 m temperature [C]
            external_data M_mixrat;       // Mixing ratio
            external_data M_mslp;         // Atmospheric pressure [Pa]
            external_data M_Qsw_in;       // Incoming short-wave radiation [W/m2]
            external_data M_Qlw_in;       // Incoming long-wave radiation [W/m2]
            external_data M_tcc;          // Incoming long-wave radiation [W/m2]
            external_data M_precip;       // Total precipitation [m]
            external_data M_snowfall;     // Snowfall rate [kg/m^2/s]
            external_data M_snowfr;       // Fraction of precipitation that is snow
            external_data M_dair;         // 2 m dew point [C]
            external_data M_sphuma;       // Speciffic humidity of the atmosphere [kg/kg]

            // Ocean
            external_data M_ocean;        // "Geostrophic" ocean currents [m/s]
            external_data M_ssh;          // Sea surface elevation [m]

            external_data M_ocean_temp;   // Ocean temperature in top layer [C]
            external_data M_ocean_salt;   // Ocean salinity in top layer [C]
            external_data M_mld;          // Mixed-layer depth [m]

            external_data M_qsrml;        // Fraction of short wave radiation absorbed by the mixed layer

            // Nesting
            external_data M_nesting_dist_elements; // Distance to the nearest open boundaries
            external_data M_nesting_dist_nodes; // Distance to the nearest open boundaries
            external_data M_nesting_conc; // sea_ice_area_fraction from the outer domain
            external_data M_nesting_thick; // sea_ice_thickness from the outer domain
            external_data M_nesting_snow_thick; // surface_snow_thickness from the outer domain
            external_data M_nesting_h_young ; // young_ice_thickness from the outer domain
            external_data M_nesting_conc_young ; // young_ice_area_fraction from the outer domain
            external_data M_nesting_hs_young ; // surface_snow_thickness_on_young_ice from the outer domain
            external_data M_nesting_damage; // damage from the outer domain
            external_data M_nesting_ridge_ratio; // ridge_ratio from the outer domain
            external_data M_nesting_VT1; // X-velocity from the outer domain
            external_data M_nesting_VT2; // Y-velocity from the outer domain
            external_data M_nesting_sigma1; // 1st component stress tensor from the outer domain
            external_data M_nesting_sigma2; // 2nd component stress tensor from the outer domain
            external_data M_nesting_sigma3; // 3rd component stress tensor from the outer domain

            // Bathymetry
            external_data M_element_depth;

            external_data_vec M_external_data_elements;
            std::vector<std::string> M_external_data_elements_names;//list of names for debugging and exporting
            external_data_vec M_external_data_nodes;
            std::vector<std::string> M_external_data_nodes_names;//list of names for debugging and exporting

        private:
            double time_init;
            double M_spinup_duration;
            int M_ensemble_member;

            setup::BathymetryType M_bathymetry_type;

            // =============================================================================
            // variables needed for nesting
            bool M_use_nesting;
            bool M_use_ocean_nesting;
            std::string M_nest_outer_mesh;
            std::string M_nest_inner_mesh;
            std::string M_nest_method;
            std::string M_nudge_function;
            double M_nudge_timescale;
            double M_nudge_lengthscale;
            bool M_nest_dynamic_vars;
            // =============================================================================

            Timer::timer* M_timer;

            po::variables_map vm;
    };
}
#endif
