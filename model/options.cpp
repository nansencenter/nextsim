/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   options.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue Jul 14 13:23:45 2015
 */

#include <boost/program_options.hpp>
#include <constants.hpp>

#if defined (WAVES)
#include <options_wim.hpp>
#endif

namespace po = boost::program_options;

double const days_in_sec = 24.0*3600.0;
namespace Nextsim
{

    po::options_description
    descrOptions()

    {
        po::options_description desc("Options");

        desc.add_options()
            ("help,h", "Print help messages")
            ("config-file", po::value<std::string>(), "specify a configuration file")
            ("config-files", po::value<std::vector<std::string> >()->multitoken(), "specify a list of configuration files" )

            /*
             *-----------------------------------------------------------------------------------
             * SIMULATION PARAMETERS
             * -----------------------------------------------------------------------------------
             */

            // - basics
            ("simul.time_init", po::value<std::string>()->default_value( "2008-Mar-05" ), "")
            ("simul.duration", po::value<double>()->default_value( 1. ), "")
            ("simul.timestep", po::value<double>()->default_value( 200. ), "")
            ("simul.spinup_duration", po::value<double>()->default_value( 1. ), "")

            /*
             *-----------------------------------------------------------------------------------
             * DEBUGGING OPTIONS
             * -----------------------------------------------------------------------------------
             */
            ("debugging.bamg_verbose", po::value<int>()->default_value( 7 ), "")
            ("debugging.log-level", po::value<std::string>()->default_value( "info" ), "")
            ("debugging.max_iteration", po::value<int>()->default_value( -1 ), "")

            /*
             *-----------------------------------------------------------------------------------
             * NUMERICS
             * -----------------------------------------------------------------------------------
             */

            // remeshing
            ("simul.regrid", po::value<std::string>()->default_value( "bamg" ), "No-regridding or bamg")
            ("simul.regrid_output_flag", po::value<bool>()->default_value( false ), "")
            ("simul.regrid_angle", po::value<double>()->default_value( 10. ), "")
            ("simul.interp_with_cavities", po::value<bool>()->default_value( true ), "")

            // advection scheme
            // - ALE_smoothing_step_nb<0 is the eulerian case where M_UM is not changed and then =0.
            // - ALE_smoothing_step_nb=0 is the purely Lagrangian case where M_UM is updated with M_VT
            // - ALE_smoothing_step_nb>0 is the ALE case where M_UM is updated with a smoothed version of M_VT
            ("simul.ALE_smoothing_step_nb", po::value<int>()->default_value( 0 ), "")

            // solver
            ("solver.ksp-type", po::value<std::string>()->default_value( "preonly" ), "")
            ("solver.pc-type", po::value<std::string>()->default_value( "cholesky" ), "")
            ("solver.mat-package-type", po::value<std::string>()->default_value( "cholmod" ), "")
            ("solver.ksp-view", po::value<bool>()->default_value( false ), "")
            ("solver.ksp-convergence-info", po::value<bool>()->default_value( true ), "")

            /*
             *-----------------------------------------------------------------------------------
             * SETUP
             * -----------------------------------------------------------------------------------
             */

            // setup
            ("setup.atmosphere-type", po::value<std::string>()->default_value( "asr" ), "")
            ("setup.ocean-type", po::value<std::string>()->default_value( "constant" ), "")
            ("setup.ice-type", po::value<std::string>()->default_value( "constant" ), "")
            ("setup.bathymetry-type", po::value<std::string>()->default_value( "etopo" ), "")
            ("setup.basal_stress-type", po::value<std::string>()->default_value( "lemieux" ), "")
            ("setup.use_assimilation", po::value<bool>()->default_value( false ), "")
            ("setup.dynamics-type", po::value<std::string>()->default_value( "default" ), "")
            ("setup.thermo-type", po::value<std::string>()->default_value( "winton" ), "")

            // mesh
            ("mesh.path", po::value<std::string>()->default_value( "nextsimdir" ), "nextsimdir or simdatadir")
            ("mesh.filename", po::value<std::string>()->default_value( "medium_Arctic_10km.msh" ), "")
            ("mesh.mppfile", po::value<std::string>()->default_value( "NpsNextsim.mpp" ), "")
            //not used: ("mesh.hsize", po::value<double>()->default_value( 0.01 ), "") // to be checked

            // moorings
            ("moorings.use_moorings", po::value<bool>()->default_value( false ), "")
            ("moorings.snapshot", po::value<bool>()->default_value( false ), "")
            ("moorings.file_length", po::value<std::string>()->default_value( "inf" ), "")
            ("moorings.spacing", po::value<double>()->default_value( 10 ), "km")
            ("moorings.output_timestep", po::value<double>()->default_value( 1 ), "days")
            ("moorings.variables", po::value<std::vector<std::string>>()->multitoken()->default_value(
                        std::vector<std::string>
                            {"conc", "thick", "snow", "conc_thin", "h_thin", "hs_thin", "velocity_xy"},
                             "conc    thick    snow    conc_thin    h_thin    hs_thin    velocity_xy"
                    )->composing(), "list of variable names (put on separate lines in config file)")
            ("moorings.grid_file", po::value<std::string>()->default_value( "" ), "") // It must be a netcdf file having x y as dimensions and latitude longitude as variables

            // drifters
            ("drifters.concentration_limit", po::value<double>()->default_value( 0.15 ), "")

            ("drifters.use_iabp_drifters", po::value<bool>()->default_value( false), "")

            ("drifters.use_osisaf_drifters", po::value<bool>()->default_value( false ), "")
            ("drifters.use_refined_osisaf_grid", po::value<bool>()->default_value( false ),
                "true: if using OSISAF drifters, use grid refined by a factor of 9, so averaged model results can be compared to the data; false: use same grid as OSISAF drift dataset")

            ("drifters.use_equally_spaced_drifters", po::value<bool>()->default_value( false), "")
            ("drifters.equally_spaced_drifters_output_time_step", po::value<double>()->default_value( 1. ), "days (must be a multiple of 0.5)")
            ("drifters.spacing", po::value<double>()->default_value( 10 ), "spacing of equally spaced drifters in km")

            ("drifters.use_rgps_drifters", po::value<bool>()->default_value( false), "")
            ("drifters.rgps_drifters_output_time_step", po::value<double>()->default_value( 0.5 ), "days (must be a multiple of 0.5)")
            ("drifters.RGPS_time_init", po::value<std::string>()->default_value( "2007-12-01" ), "")

            ("drifters.use_sidfex_drifters", po::value<bool>()->default_value( false), "are we using SIDFEX drifters")
            ("drifters.sidfex_drifters_output_time_step", po::value<double>()->default_value( 0.5 ), "days (must be a multiple of 0.5)")
            ("drifters.sidfex_filename", po::value<std::string>()->default_value( "" ), "text file with initial buoy positions")


            // restart
            // - inputs
            ("restart.start_from_restart", po::value<bool>()->default_value( false ), "")
            ("restart.restart_string", po::value<std::string>()->default_value( "" ), "")
            ("restart.step_nb", po::value<int>()->default_value( 0 ), "")
            ("restart.input_path", po::value<std::string>()->default_value( "" ),
                    "where to find restarts (default is $NEXTSIMDIR/restart)")
            ("restart.restart_at_rest", po::value<bool>()->default_value( false ), "")
            ("restart.reset_time_counter", po::value<bool>()->default_value( false ),
                "false: simulation starts at simul.time_init eg for forecast; true: simulation starts at simul.time_init+pcpt*simul.timestep eg to restart interrupted simulation")

            // - outputs
            ("restart.write_restart", po::value<bool>()->default_value( false ), "")
            ("restart.output_time_step", po::value<double>()->default_value( 15 ), "days")
            ("restart.debugging", po::value<bool>()->default_value( false ),
                "save restart every time step for debugging (only with build type DEBUG)")

            // general outputs
            ("output.output_per_day", po::value<int>()->default_value( 4 ), "")
            ("output.logfile", po::value<std::string>()->default_value( "" ), "")
            ("output.save_forcing_fields", po::value<bool>()->default_value( false ), "")
            ("output.save_diagnostics", po::value<bool>()->default_value( false ), "")

            // exporter
            ("output.datetime_in_filename", po::value<bool>()->default_value( false ), "filename outputs are eg [mesh,field]_20180101T000000Z.[bin,dat]")
            ("output.exporter_path", po::value<std::string>()->default_value( "" ), "")
            ("output.exporter_precision", po::value<std::string>()->default_value("float"),
                    "float (default) or double (almost only for testing)")

            /*
             *-----------------------------------------------------------------------------------
             IDEAL CONFIGURATIONS
             * -----------------------------------------------------------------------------------
             */

            // - if setup.bathymetry-type = constant
            ("ideal_simul.constant_bathymetry", po::value<double>()->default_value( 200. ), "")

            // - if setup.ice-type = constant
            ("ideal_simul.init_thickness", po::value<double>()->default_value( 1.0 ), "")
            ("ideal_simul.init_concentration", po::value<double>()->default_value( 1.0 ), "")
            ("ideal_simul.init_thin_conc", po::value<double>()->default_value( 0. ), "")
            ("ideal_simul.init_snow_thickness", po::value<double>()->default_value( 0. ), "")

            // - if atmosphere-type=constant
            // -- thermodynamics
            ("ideal_simul.constant_tair", po::value<double>()->default_value( -25. ), "")
            ("ideal_simul.constant_dair", po::value<double>()->default_value( -1. ), "")
            ("ideal_simul.constant_mixrat", po::value<double>()->default_value( 0.001 ), "")
            ("ideal_simul.constant_mslp", po::value<double>()->default_value( 1013e2 ), "")
            ("ideal_simul.constant_Qsw_in", po::value<double>()->default_value( 50. ), "")
            ("ideal_simul.constant_Qlw_in", po::value<double>()->default_value( 250. ), "")
            ("ideal_simul.constant_precip", po::value<double>()->default_value( 1e-5 ), "")
            ("ideal_simul.constant_snowfr", po::value<double>()->default_value( 0.9 ), "")
            ("ideal_simul.constant_Qdw", po::value<double>()->default_value( 0. ), "")
            ("ideal_simul.constant_Fdw", po::value<double>()->default_value( 0. ), "")
            ("ideal_simul.constant_mld", po::value<double>()->default_value( 9. ), "")

            // - if atmosphere-type=constant
            // -- dynamics
            ("ideal_simul.constant_wind_u", po::value<double>()->default_value( 0. ), "")
            ("ideal_simul.constant_wind_v", po::value<double>()->default_value( -10. ), "")

            // - if ocean-type=constant
            ("ideal_simul.constant_ocean_u", po::value<double>()->default_value( 0. ), "")
            ("ideal_simul.constant_ocean_v", po::value<double>()->default_value( 0. ), "")
            ("ideal_simul.constant_ssh", po::value<double>()->default_value( 0. ), "")
            // -----------------------------------------------------------------------------------


            /*
             *-----------------------------------------------------------------------------------
             * DYNAMICS
             * -----------------------------------------------------------------------------------
             */

            // - internal stresses
            ("dynamics.alea_factor", po::value<double>()->default_value( 0. ), "")
            ("dynamics.young", po::value<double>()->default_value( 5.49e+9 ), "Pa") // 5.49e+9 is a more reasonable than 9GPa, and same as used in WIM paper.
            ("dynamics.cfix", po::value<double>()->default_value( 40e+3 ), "Pa")
            ("dynamics.nu0", po::value<double>()->default_value( 0.3 ), "")
            ("dynamics.tan_phi", po::value<double>()->default_value( 0.7 ), "")
            ("dynamics.tract_coef", po::value<double>()->default_value( 5./6 ), "")
            ("dynamics.compr_strength", po::value<double>()->default_value( 750e+3 ), "Pa")
            ("dynamics.ridging_exponent", po::value<double>()->default_value( -20. ), "")

            // - c,h limits for where to use MEB rheology and where to use the Laplacian free drift thing
            ("dynamics.min_h", po::value<double>()->default_value( 0.05 ), "")
            ("dynamics.min_c", po::value<double>()->default_value( 0.01 ), "")

            // - Ratio of ridged ice cohesion and compressive strength compared to level ice (1. does nothing)
            ("dynamics.ridge_to_normal_cohesion_ratio", po::value<double>()->default_value( 1. ), "")
            // - Scaling of cohesion w.r.t. ice thickness (normalisation factor = 1 and exponent = 0 does nothing)
            ("dynamics.cohesion_thickness_normalisation", po::value<double>()->default_value( 1. ), "")
            ("dynamics.cohesion_thickness_exponent", po::value<double>()->default_value( 1. ), "")
            // - scaling with respect to horizontal resolution
            ("dynamics.scale_coef", po::value<double>()->default_value( 0.1 ), "")

            ("dynamics.use_temperature_dependent_healing", po::value<bool>()->default_value( false ), "")
            ("dynamics.time_relaxation_damage", po::value<double>()->default_value( 25. ), "days")
            ("dynamics.deltaT_relaxation_damage", po::value<double>()->default_value( 20. ), "Kelvin")
            ("dynamics.undamaged_time_relaxation_sigma", po::value<double>()->default_value( 1e7 ), "seconds")
                // from V. Dansereau et al.: A Maxwell elasto-brittle rheology for sea ice modelling
            ("dynamics.exponent_relaxation_sigma", po::value<double>()->default_value( 4. ), "")
                // from V. Dansereau et al.: A Maxwell elasto-brittle rheology for sea ice modelling

            // - water and air drag parameterizations
            ("dynamics.ERAi_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
            ("dynamics.ECMWF_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
            ("dynamics.ASR_quad_drag_coef_air", po::value<double>()->default_value( 0.0049 ), "")
            ("dynamics.CFSR_quad_drag_coef_air", po::value<double>()->default_value( 0.0023 ), "")
                // Updated value, based on comparison with OSISAF drift in the free drift case
            ("dynamics.lin_drag_coef_air", po::value<double>()->default_value( 0. ), "")
            ("dynamics.quad_drag_coef_water", po::value<double>()->default_value( 0.0055 ), "")
            ("dynamics.lin_drag_coef_water", po::value<double>()->default_value( 0. ), "")

            // - Coriolis force
            ("dynamics.use_coriolis", po::value<bool>()->default_value( true ), "")
            ("dynamics.oceanic_turning_angle", po::value<double>()->default_value( 25. ), "")

            // - basal stress parameterization
            ("dynamics.Lemieux_basal_k1", po::value<double>()->default_value( 10. ), "")
            ("dynamics.Lemieux_basal_k2", po::value<double>()->default_value( 15. ), "")
            ("dynamics.Lemieux_basal_Cb", po::value<double>()->default_value( 20. ), "")
            ("dynamics.Lemieux_basal_u_0", po::value<double>()->default_value( 5e-5 ), "")
            ("dynamics.Lemieux_basal_u_crit", po::value<double>()->default_value( 5e-4 ), "")

            /*
             *-----------------------------------------------------------------------------------
             * THERMODYNAMICS
             * -----------------------------------------------------------------------------------
             */
            ("thermo.Qio-type", po::value<std::string>()->default_value( "basic" ), "")
            ("thermo.use_thermo_forcing", po::value<bool>()->default_value( true ), "")
            ("thermo.albedoW", po::value<double>()->default_value( 0.07 ), "")
            ("thermo.alb_scheme", po::value<int>()->default_value( 3 ), "")
            ("thermo.flooding", po::value<bool>()->default_value( true ), "")
            ("thermo.alb_ice", po::value<double>()->default_value( 0.63 ), "")
            ("thermo.alb_sn", po::value<double>()->default_value( 0.88 ), "")
            ("thermo.I_0", po::value<double>()->default_value( 0.17 ), "")
            ("thermo.Qdw", po::value<double>()->default_value( 0.5 ), "")
            ("thermo.Fdw", po::value<double>()->default_value( 0. ), "")
            ("thermo.newice_type", po::value<int>()->default_value( 4 ), "4: THIN_ICE; else CLASSIC")
            ("thermo.melt_type", po::value<int>()->default_value( 1 ), "")
            ("thermo.hnull", po::value<double>()->default_value( 0.25 ), "")
            ("thermo.PhiF", po::value<double>()->default_value( 4. ), "")
            ("thermo.PhiM", po::value<double>()->default_value( 0.5 ), "")
            ("thermo.h_thin_max", po::value<double>()->default_value( 0.5 ), "")
            ("thermo.h_thin_min", po::value<double>()->default_value( 0.05 ), "")

            ("thermo.drag_ice_t", po::value<double>()->default_value( 1.3e-3 ), "")
            ("thermo.drag_ocean_u", po::value<double>()->default_value( 1.1e-3 ), "")
            ("thermo.drag_ocean_t", po::value<double>()->default_value( 0.83e-3 ), "")
            ("thermo.drag_ocean_q", po::value<double>()->default_value( 1.5e-3 ), "")

            // - diffusivity
            ("thermo.diffusivity_sss", po::value<double>()->default_value( 100. ), "") //[m^2/s]
            ("thermo.diffusivity_sst", po::value<double>()->default_value( 100. ), "") //[m^2/s]

            // - relaxation of slab ocean to ocean forcing
            ("thermo.ocean_nudge_timeT", po::value<double>()->default_value( 30*days_in_sec), "")
            ("thermo.ocean_nudge_timeS", po::value<double>()->default_value( 30*days_in_sec), "")

            /*
             *-----------------------------------------------------------------------------------
             * NESTING
             * -----------------------------------------------------------------------------------
             */

            ("nesting.use_nesting", po::value<bool>()->default_value( false ), "")
            ("nesting.use_ocean_nesting", po::value<bool>()->default_value( false ), "")
            ("nesting.outer_mesh",po::value<std::string>()->default_value( "" ),"")
            ("nesting.inner_mesh",po::value<std::string>()->default_value( "" ),"")
            ("nesting.method", po::value<std::string>()->default_value( "nudging" ), "")
            ("nesting.nudge_function", po::value<std::string>()->default_value( "exponential" ), "")
            ("nesting.nudge_timescale", po::value<double>()->default_value((1./2.)*days_in_sec), "")
            ("nesting.nudge_lengthscale", po::value<double>()->default_value(10.), "")
            ("nesting.nest_dynamic_vars", po::value<bool>()->default_value( false ), "")

            /*
             *-----------------------------------------------------------------------------------
             * FORECASTING
             * -----------------------------------------------------------------------------------
             */
            ("forecast.air_temperature_correction", po::value<double>()->default_value( 0. ),
                "for use in BADA (Bias-Aware-Data-Assimilation)")
            ("forecast.true_forecast", po::value<bool>()->default_value( false ),
                "if false, use <<analysis/best est>> - forecast that started on the same day as the model time; else use the one that began on simul.time_init, or  forecast.time_init_atm_fc/forecast.time_init_ocean_fc (if provided)")
            ("forecast.time_init_atm_fc", po::value<std::string>()->default_value( "" ),
                "if(forecast.true_forecast), get atmospheric forecast starting from this date as opposed to simul.time_init (eg if usual one is absent)")
            ("forecast.time_init_ocean_fc", po::value<std::string>()->default_value( "" ),
                "if(forecast.true_forecast), get ocean forecast starting from this date as opposed to simul.time_init (eg if usual one is absent)")

            /*
             *-----------------------------------------------------------------------------------
             * MODULES
             * -----------------------------------------------------------------------------------
             */
#if defined(OASIS)
            ("coupler.timestep", po::value<double>()->default_value( 3600. ), "") // The coupling time step
            ("coupler.with_ocean", po::value<bool>()->default_value( false ), "") // Do we couple with an ocean model
            ("coupler.atm_from_ocean", po::value<bool>()->default_value( false ), "") // Do we get atmospheric state from the ocean model
            ("coupler.with_waves", po::value<bool>()->default_value( false ), "") // Do we couple with a wave model
            ("coupler.with_atm", po::value<bool>()->default_value( false ), "")  // Do we couple with an atmospheric model
#endif

#if defined(WAVES)
        ;
        return desc.add( Wim::descrWimOptions() );
#else
        ;
        return desc;
#endif
    }

} // Nextsim
