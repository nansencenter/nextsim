/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   options.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue Jul 14 13:23:45 2015
 */

#include <boost/program_options.hpp>
#include <constants.hpp>

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
             * NUMERICS
             * -----------------------------------------------------------------------------------
             */

            // solver
            ("solver.ksp-type", po::value<std::string>()->default_value( "preonly" ), "")
            ("solver.pc-type", po::value<std::string>()->default_value( "cholesky" ), "")
            ("solver.mat-package-type", po::value<std::string>()->default_value( "cholmod" ), "")
            ("solver.ksp-view", po::value<bool>()->default_value( false ), "")
            ("solver.ksp-convergence-info", po::value<bool>()->default_value( true ), "")
            ("solver.ksp-reuse-prec", po::value<bool>()->default_value( false ), "")
            ("solver.ksp-monitor", po::value<bool>()->default_value( false ), "")

            // debugging options
            ("simul.verbose", po::value<int>()->default_value( 7 ), "")
            ("simul.log-level", po::value<std::string>()->default_value( "info" ), "")
            ("simul.maxiteration", po::value<int>()->default_value( 1e+8 ), "")

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
            ("simul.basal_stress-type", po::value<std::string>()->default_value( "lemieux" ), "")
            ("simul.expansion_factor", po::value<double>()->default_value( 0.05 ),
                    "used in DataSet::loadGrid - makes loaded grid slightly larger than actually needed (default is by 5%)")

            // mesh
            ("mesh.path", po::value<std::string>()->default_value( "nextsimdir" ), "nextsimdir or simdatadir")
            ("mesh.filename", po::value<std::string>()->default_value( "medium_Arctic_10km.msh" ), "")
            ("mesh.fileformat", po::value<std::string>()->default_value( "binary" ), "")
            ("mesh.mppfile", po::value<std::string>()->default_value( "NpsNextsim.mpp" ), "")
            ("mesh.hsize", po::value<double>()->default_value( 0.01 ), "")
            ("setup.domain-type", po::value<std::string>()->default_value( "unref" ), "")
            ("mesh.partitioner", po::value<std::string>()->default_value( "metis" ), "mesh partitioner: chaco or metis")
            ("mesh.partition-space", po::value<std::string>()->default_value( "disk" ), "")

            // simul
            ("simul.time_init", po::value<std::string>()->default_value( "2008-Mar-05" ), "")
            ("simul.duration", po::value<double>()->default_value( 1. ), "")
            ("simul.timestep", po::value<double>()->default_value( 200. ), "")
            ("simul.spinup_duration", po::value<double>()->default_value( 1. ), "")

            // frequency for info output
            ("simul.ptime_per_day", po::value<int>()->default_value( 12 ), "")

            // moorings
            ("simul.use_moorings", po::value<bool>()->default_value( false ), "")
            ("simul.mooring_snapshot", po::value<bool>()->default_value( false ), "")
            ("simul.mooring_file_length", po::value<std::string>()->default_value( "inf" ), "")
            ("simul.mooring_spacing", po::value<double>()->default_value( 10 ), "km")
            ("simul.mooring_output_timestep", po::value<double>()->default_value( 1 ), "days")
            ("simul.mooring_names", po::value<std::vector<std::string>>()->multitoken()->default_value(std::vector<std::string>(),
                    "conc thick snow conc_thin h_thin hs_thin velocity_xy")->composing(), "")
            ("simul.mooring_grid_file", po::value<std::string>()->default_value( "" ), "") // It must be a netcdf file having x y as dimensions and latitude longitude as variables
            ("simul.mooring_parallel_output", po::value<bool>()->default_value( false ), "")

            // drifters
            ("simul.use_iabp_drifters", po::value<bool>()->default_value( false), "")
            ("simul.equallyspaced_drifters_output_time_step", po::value<double>()->default_value( 0. ), "days") // must be a multiple of 0.5
            ("simul.rgps_drifters_output_time_step", po::value<double>()->default_value( 0. ), "days") // must be a multiple of 0.5
            ("simul.use_osisaf_drifters", po::value<bool>()->default_value( false ), "")
            ("simul.drifter_climit", po::value<double>()->default_value( 0.15 ), "")
            ("simul.drifter_spacing", po::value<double>()->default_value( 10 ), "")
            ("simul.RGPS_time_init", po::value<std::string>()->default_value( "2008-Jan-01" ), "")

            // bathymetry
            ("simul.constant_bathymetry", po::value<double>()->default_value( 200. ), "")

            // if ice-type=constant
            ("simul.init_thickness", po::value<double>()->default_value( 1.0 ), "")
            ("simul.init_concentration", po::value<double>()->default_value( 1.0 ), "")
            ("simul.init_thin_conc", po::value<double>()->default_value( 0. ), "")
            ("simul.init_snow_thickness", po::value<double>()->default_value( 0. ), "")

            //inputs (restart)
            ("setup.use_restart", po::value<bool>()->default_value( false ), "")
            ("setup.step_nb", po::value<int>()->default_value( 0 ), "")
            ("setup.restart_path", po::value<std::string>()->default_value( "" ),
                    "where to find restarts (default is $NEXTSIMDIR/restart)")

            // outputs (restart)
            ("setup.write_restart", po::value<bool>()->default_value( false ), "")
            ("setup.restart_at_rest", po::value<bool>()->default_value( false ), "")
            ("setup.restart_time_step", po::value<double>()->default_value( 20 ), "days")

            // other outputs
            ("simul.output_per_day", po::value<int>()->default_value( 24 ), "")
            ("simul.logfile", po::value<std::string>()->default_value( "" ), "")
            ("simul.save_forcing_field", po::value<bool>()->default_value( false ), "")
            ("simul.save_diagnostics", po::value<bool>()->default_value( false ), "")

            // exporter
            ("exporter.path", po::value<std::string>()->default_value( "" ), "")
            ("exporter.precision", po::value<std::string>()->default_value("float"),
                    "float (default) or double (almost only for testing)")

            /*
             *-----------------------------------------------------------------------------------
             * DYNAMICS
             * -----------------------------------------------------------------------------------
             */
            ("setup.dynamics-type", po::value<std::string>()->default_value( "default" ), "")

            // internal stresses
            ("simul.use_coriolis", po::value<bool>()->default_value( true ), "")
            ("simul.alea_factor", po::value<double>()->default_value( 0. ), "")
            ("simul.young", po::value<double>()->default_value( 5.49e+9 ), "Pa") // 5.49e+9 is a more reasonable than 9GPa, and same as used in WIM paper.
            ("simul.cfix", po::value<double>()->default_value( 40e+3 ), "Pa")
            ("simul.nu0", po::value<double>()->default_value( 0.3 ), "")
            ("simul.tan_phi", po::value<double>()->default_value( 0.7 ), "")
            ("simul.tract_coef", po::value<double>()->default_value( 5./6 ), "")
            ("simul.compr_strength", po::value<double>()->default_value( 750e+3 ), "Pa")
            ("simul.ridging_exponent", po::value<double>()->default_value( -20. ), "")

            // - Ratio of ridged ice cohesion and compressive strength compared to level ice (1. does nothing)
            ("simul.ridge_to_normal_cohesion_ratio", po::value<double>()->default_value( 1. ), "")
            // - Scaling of cohesion w.r.t. ice thickness (normalisation factor = 1 and exponent = 0 does nothing)
            ("simul.cohesion_thickness_normalisation", po::value<double>()->default_value( 1. ), "")
            ("simul.cohesion_thickness_exponent", po::value<double>()->default_value( 1. ), "")
            // - scaling with respect to horizontal resolution
            ("simul.scale_coef", po::value<double>()->default_value( 0.1 ), "")

            ("simul.use_temperature_dependent_healing", po::value<bool>()->default_value( false ), "")
            ("simul.time_relaxation_damage", po::value<double>()->default_value( 25. ), "days")
            ("simul.deltaT_relaxation_damage", po::value<double>()->default_value( 20. ), "Kelvin")
            ("simul.undamaged_time_relaxation_sigma", po::value<double>()->default_value( 1e7 ), "seconds")
                // from V. Dansereau et al.: A Maxwell elasto-brittle rheology for sea ice modelling
            ("simul.exponent_relaxation_sigma", po::value<double>()->default_value( 4. ), "")
                // from V. Dansereau et al.: A Maxwell elasto-brittle rheology for sea ice modelling

            // water and air drag parameterizations
            ("simul.ERAi_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
            ("simul.ECMWF_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
            ("simul.ASR_quad_drag_coef_air", po::value<double>()->default_value( 0.0049 ), "")
            ("simul.CFSR_quad_drag_coef_air", po::value<double>()->default_value( 0.0023 ), "")
                // Updated value, based on comparison with OSISAF drift in the free drift case
            ("simul.lin_drag_coef_air", po::value<double>()->default_value( 0. ), "")
            ("simul.quad_drag_coef_water", po::value<double>()->default_value( 0.0055 ), "")
            ("simul.lin_drag_coef_water", po::value<double>()->default_value( 0. ), "")

            // basal stress parameterization
            ("simul.Lemieux_basal_k1", po::value<double>()->default_value( 10. ), "")
            ("simul.Lemieux_basal_k2", po::value<double>()->default_value( 15. ), "")
            ("simul.Lemieux_basal_Cb", po::value<double>()->default_value( 20. ), "")
            ("simul.Lemieux_basal_u_0", po::value<double>()->default_value( 5e-5 ), "")
            ("simul.Lemieux_basal_u_crit", po::value<double>()->default_value( 5e-4 ), "")

            /*
             *-----------------------------------------------------------------------------------
             * THERMODYNAMICS
             * -----------------------------------------------------------------------------------
             */
            ("setup.thermo-type", po::value<std::string>()->default_value( "winton" ), "")
            ("simul.Qio-type", po::value<std::string>()->default_value( "basic" ), "")
            ("simul.use_thermo_forcing", po::value<bool>()->default_value( true ), "")
            ("simul.albedoW", po::value<double>()->default_value( 0.07 ), "")
            ("simul.alb_scheme", po::value<int>()->default_value( 3 ), "")
            ("simul.flooding", po::value<bool>()->default_value( true ), "")
            ("simul.alb_ice", po::value<double>()->default_value( 0.63 ), "")
            ("simul.alb_sn", po::value<double>()->default_value( 0.88 ), "")
            ("simul.I_0", po::value<double>()->default_value( 0.17 ), "")
            ("simul.Qdw", po::value<double>()->default_value( 0.5 ), "")
            ("simul.Fdw", po::value<double>()->default_value( 0. ), "")
            ("simul.newice_type", po::value<int>()->default_value( 4 ), "4: THIN_ICE; else CLASSIC")
            ("simul.melt_type", po::value<int>()->default_value( 1 ), "")
            ("simul.hnull", po::value<double>()->default_value( 0.25 ), "")
            ("simul.PhiF", po::value<double>()->default_value( 4. ), "")
            ("simul.PhiM", po::value<double>()->default_value( 0.5 ), "")
            ("simul.h_thin_max", po::value<double>()->default_value( 0.5 ), "")
            ("simul.h_thin_min", po::value<double>()->default_value( 0.05 ), "")

            ("simul.drag_ice_t", po::value<double>()->default_value( 1.3e-3 ), "")
            ("simul.drag_ocean_u", po::value<double>()->default_value( 1.1e-3 ), "")
            ("simul.drag_ocean_t", po::value<double>()->default_value( 0.83e-3 ), "")
            ("simul.drag_ocean_q", po::value<double>()->default_value( 1.5e-3 ), "")

            // diffusivity
            ("simul.diffusivity_sss", po::value<double>()->default_value( 100. ), "") //[m^2/s]
            ("simul.diffusivity_sst", po::value<double>()->default_value( 100. ), "") //[m^2/s]

            // if atmosphere-type=constant
            ("simul.constant_tair", po::value<double>()->default_value( -25. ), "")
            ("simul.constant_dair", po::value<double>()->default_value( -1. ), "")
            ("simul.constant_mixrat", po::value<double>()->default_value( 0.001 ), "")
            ("simul.constant_mslp", po::value<double>()->default_value( 1013e2 ), "")
            ("simul.constant_Qsw_in", po::value<double>()->default_value( 50. ), "")
            ("simul.constant_Qlw_in", po::value<double>()->default_value( 250. ), "")
            ("simul.constant_precip", po::value<double>()->default_value( 1e-5 ), "")
            ("simul.constant_snowfr", po::value<double>()->default_value( 0.9 ), "")
            ("simul.constant_Qdw", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_Fdw", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_mld", po::value<double>()->default_value( 9. ), "")

            /*
             *-----------------------------------------------------------------------------------
             * EXTERNAL FORCINGS
             * -----------------------------------------------------------------------------------
             */
            ("simul.ocean_nudge_timeT", po::value<double>()->default_value( 30*days_in_sec), "")
            ("simul.ocean_nudge_timeS", po::value<double>()->default_value( 30*days_in_sec), "")
            ("simul.constant_wind_u", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_wind_v", po::value<double>()->default_value( -10. ), "")
            ("simul.constant_ocean_u", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_ocean_v", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_ssh", po::value<double>()->default_value( 0. ), "")
            ("simul.oceanic_turning_angle", po::value<double>()->default_value( 25. ), "")

            /*
             *-----------------------------------------------------------------------------------
             * FORECASTING
             * -----------------------------------------------------------------------------------
             */
            ("simul.forecast", po::value<bool>()->default_value( false ), "")
            ("simul.use_ship", po::value<bool>()->default_value( false ), "")
            ("simul.min_h", po::value<double>()->default_value( 0.05 ), "")
            ("simul.min_c", po::value<double>()->default_value( 0.01 ), "")

            /*
             *-----------------------------------------------------------------------------------
             * OTHERS (TO BE SORTED)
             * -----------------------------------------------------------------------------------
             */
            //only used in commented-out-test: ("simul.drift_limit_concentration", po::value<double>()->default_value( 0.05 ), "")
            //not used: ("simul.ERAcorr2T", po::value<std::vector<std::string>>()->multitoken()->zero_tokens()->composing(), "")
            //not used: ("simul.close_all_boundaries", po::value<bool>()->default_value( false ), "")
            //not used: ("simul.use_stability_drag", po::value<double>()->default_value( 0. ), "")

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
            ("simul.use_wim", po::value<bool>()->default_value( false ), "")
            ("simul.wim_grid", po::value<bool>()->default_value( false ), "")
            ("setup.wave-type", po::value<std::string>()->default_value( "constant" ), "constant, ww3a")
            ("setup.wave-time-interp-option", po::value<std::string>()->default_value( "step" ), "step, linear")

            ("simul.constant_significant_wave_height", po::value<double>()->default_value( 2. ), "")
            ("simul.constant_wave_mean_direction", po::value<double>()->default_value( 90. ), "")
            ("simul.constant_wave_peak_frequency", po::value<double>()->default_value( 0.05 ), "")

            ("simul.export_after_wim_call", po::value<bool>()->default_value( false ), "")
            ;

        desc.add( Wim::descrWimOptions() );
#endif
        ;

        return desc;
    }

} // Nextsim
