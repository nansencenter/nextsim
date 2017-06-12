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
#include <wimoptions.hpp>
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
             * NUMERICS
             * -----------------------------------------------------------------------------------
             */

            // solver
            ("solver.ksp-type", po::value<std::string>()->default_value( "preonly" ), "")
            ("solver.pc-type", po::value<std::string>()->default_value( "cholesky" ), "")
            ("solver.mat-package-type", po::value<std::string>()->default_value( "cholmod" ), "")
            ("solver.ksp-view", po::value<bool>()->default_value( false ), "")
            ("solver.ksp-convergence-info", po::value<bool>()->default_value( true ), "")

            // debugging options
            ("simul.verbose", po::value<int>()->default_value( 7 ), "")
            ("simul.log-level", po::value<std::string>()->default_value( "info" ), "")
            ("simul.maxiteration", po::value<int>()->default_value( -1 ), "")


            // remesher
            ("simul.regrid", po::value<std::string>()->default_value( "bamg" ), "No-regridding or bamg")
            ("simul.regrid_output_flag", po::value<bool>()->default_value( false ), "")
            ("simul.regrid_angle", po::value<double>()->default_value( 10. ), "")
            ("simul.interp_with_cavities", po::value<bool>()->default_value( true ), "")
            ("simul.hsize", po::value<double>()->default_value( 0.01 ), "") // to be checked


            // advection scheme
            ("simul.ALE_smoothing_step_nb", po::value<int>()->default_value( 0 ), "")
            // ALE_smoothing_step_nb<0 is the eulerian case where M_UM is not changed and then =0.
            // ALE_smoothing_step_nb=0 is the purely Lagrangian case where M_UM is updated with M_VT
            // ALE_smoothing_step_nb>0 is the ALE case where M_UM is updated with a smoothed version of M_VT

            // TODO: partitioner (for the parallel version)

            /*
             *-----------------------------------------------------------------------------------
             * SETUP
             * -----------------------------------------------------------------------------------
             */

            // setup
            ("setup.atmosphere-type", po::value<std::string>()->default_value( "asr" ), "")
            ("setup.ocean-type", po::value<std::string>()->default_value( "constant" ), "")

            ("setup.ice-type", po::value<std::string>()->default_value( "constant" ), "")
            // if ice-type=constant
            // not used: ("simul.init_damage", po::value<std::string>()->default_value( "No-damage" ), "")
            ("simul.init_thickness", po::value<double>()->default_value( 1.0 ), "")
            ("simul.init_concentration", po::value<double>()->default_value( 1.0 ), "")
            ("simul.init_thin_conc", po::value<double>()->default_value( 0. ), "")
            ("simul.init_snow_thickness", po::value<double>()->default_value( 0. ), "")
            // not used: ("simul.init_snow_thin_thickness", po::value<double>()->default_value( 0. ), "")

            // outputs
            ("setup.use_restart", po::value<bool>()->default_value( false ), "")
            ("setup.write_restart", po::value<bool>()->default_value( false ), "")
            ("setup.restart_at_rest", po::value<bool>()->default_value( false ), "")
            ("setup.restart_time_step", po::value<double>()->default_value( 15 ), "days")
            ("setup.step_nb", po::value<int>()->default_value( 0 ), "")
            ("simul.output_per_day", po::value<int>()->default_value( 4 ), "")
            ("simul.output_directory", po::value<std::string>()->default_value( "" ), "")
            ("simul.logfile", po::value<std::string>()->default_value( "" ), "")
            // not used: ("simul.diagnostics_frequency", po::value<int>()->default_value( 1 ), "")
            ("simul.save_forcing_field", po::value<bool>()->default_value( false ), "")
            // not used: ("simul.slim_bin.use", po::value<bool>()->default_value( false ), "")
            // not used: ("simul.slim_bin.names", po::value<std::vector<std::string>>()->multitoken()->zero_tokens()->composing(), "")
            // not used: ("simul.slim_bin.output_timestep", po::value<double>()->default_value( std::numeric_limits<double>::infinity() ), "")
            ("setup.exporter_precision", po::value<std::string>()->default_value("float"), "float (default) or double (almost only for testing)")

            // simulation
            ("simul.mesh_filename", po::value<std::string>()->default_value( "bigarctic10km.msh" ), "")
            ("simul.mesh_path", po::value<std::string>()->default_value( "nextsimdir" ), "nextsimdir or simdatadir")
            ("simul.proj_filename", po::value<std::string>()->default_value( "NpsNextsim.mpp" ), "")
            ("simul.time_init", po::value<std::string>()->default_value( "2008-Mar-05" ), "")
            ("simul.duration", po::value<double>()->default_value( 1. ), "")
            ("setup.bathymetry-type", po::value<std::string>()->default_value( "etopo" ), "")
            ("simul.basal_stress-type", po::value<std::string>()->default_value( "lemieux" ), "")
            ("simul.constant_bathymetry", po::value<double>()->default_value( 200. ), "")
            ("simul.timestep", po::value<double>()->default_value( 200. ), "")
            ("simul.spinup_duration", po::value<double>()->default_value( 1. ), "")

            // TODO: moorings and drifters
            ("simul.use_moorings", po::value<bool>()->default_value( false ), "")
            ("simul.mooring_snapshot", po::value<bool>()->default_value( false ), "")
            ("simul.mooring_file_length", po::value<std::string>()->default_value( "inf" ), "")
            // not used: ("simul.mooring_names", po::value<std::vector<std::string>>()->multitoken()->zero_tokens()->composing(), "")
            ("simul.mooring_spacing", po::value<double>()->default_value( 10 ), "km")
            // not used: ("simul.mooring_grid_file", po::value<std::string>()->default_value( "" ), "")
            ("simul.mooring_output_timestep", po::value<double>()->default_value( 1 ), "days")

            ("simul.use_iabp_drifters", po::value<bool>()->default_value( false), "")
            ("simul.equallyspaced_drifters_output_time_step", po::value<double>()->default_value( 0. ), "days") // must be a multiple of 0.5
            ("simul.rgps_drifters_output_time_step", po::value<double>()->default_value( 0. ), "days") // must be a multiple of 0.5
            ("simul.use_osisaf_drifters", po::value<bool>()->default_value( false ), "") 
            
            ("simul.drifter_climit", po::value<double>()->default_value( 0.15 ), "")
            ("simul.drifter_spacing", po::value<double>()->default_value( 10 ), "")

            ("simul.RGPS_time_init", po::value<std::string>()->default_value( "2007-12-01" ), "")
               
            /*
             *-----------------------------------------------------------------------------------
             * DYNAMICS
             * -----------------------------------------------------------------------------------
             */

            // internal stresses
            // not used: ("simul.fixed_cohesion", po::value<bool>()->default_value( true ), "")
            ("simul.use_coriolis", po::value<bool>()->default_value( true ), "")
            ("simul.alea_factor", po::value<double>()->default_value( 0. ), "")
            ("simul.young", po::value<double>()->default_value( 5.49e+9 ), "Pa") // 5.49e+9 is a more reasonable than 9GPa, and same as used in WIM paper.
            ("simul.cfix", po::value<double>()->default_value( 40e+3 ), "Pa")
            ("simul.nu0", po::value<double>()->default_value( 0.3 ), "")
            ("simul.tan_phi", po::value<double>()->default_value( 0.7 ), "")
            ("simul.tract_coef", po::value<double>()->default_value( 5./6 ), "")
            ("simul.compr_strength", po::value<double>()->default_value( 750e+3 ), "Pa")
            ("simul.ridging_exponent", po::value<double>()->default_value( -20. ), "")
            ("simul.ridge_to_normal_cohesion_ratio", po::value<double>()->default_value( 1. ), "") // Ratio of ridged ice cohesion and compressive strength compared to level ice (1. does nothing)
            // Scaling of cohesion w.r.t. ice thickness (normalisation factor = 1 and exponent = 0 does nothing)
            ("simul.cohesion_thickness_normalisation", po::value<double>()->default_value( 1. ), "")
            ("simul.cohesion_thickness_exponent", po::value<double>()->default_value( 1. ), "")

            ("simul.time_relaxation_damage", po::value<double>()->default_value( 25. ), "days")
            ("simul.deltaT_relaxation_damage", po::value<double>()->default_value( 40. ), "Kelvin")
            
            ("simul.undamaged_time_relaxation_sigma", po::value<double>()->default_value( 1e7 ), "seconds") // from V. Dansereau et al.: A Maxwell elasto-brittle rheology for sea ice modelling
            ("simul.exponent_relaxation_sigma", po::value<double>()->default_value( 4. ), "")              // from V. Dansereau et al.: A Maxwell elasto-brittle rheology for sea ice modelling

            // scaling with respect to horizontal resolution
            ("simul.scale_coef", po::value<double>()->default_value( 0.2 ), "")
            // not used: ("simul.rheology_scaling", po::value<bool>()->default_value( false ), "")
            // not used: ("simul.rheology_scaling_ref_dx", po::value<double>()->default_value( 10. ), "")
            // not used: ("simul.rheology_scaling_factor", po::value<double>()->default_value( -0.0638 ), "")


            // water and air drag parameterizations
            ("simul.ERAi_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
            ("simul.ECMWF_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
            ("simul.ASR_quad_drag_coef_air", po::value<double>()->default_value( 0.0049 ), "")
            ("simul.CFSR_quad_drag_coef_air", po::value<double>()->default_value( 0.0023 ), "") // Updated value, based on comparison with OSISAF drift in the free drift case
            ("simul.lin_drag_coef_air", po::value<double>()->default_value( 0. ), "")
            ("simul.quad_drag_coef_water", po::value<double>()->default_value( 0.0055 ), "")
            ("simul.lin_drag_coef_water", po::value<double>()->default_value( 0. ), "")

            // basal stress parameterization
            //("simul.Lemieux_basal_gamma", po::value<double>()->default_value( 10. ), "")
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

            ("setup.thermo-type", po::value<std::string>()->default_value( "zero-layer" ), "")
            ("simul.Qio-type", po::value<std::string>()->default_value( "basic" ), "")
            ("simul.use_thermo_forcing", po::value<bool>()->default_value( true ), "")
            ("simul.albedoW", po::value<double>()->default_value( 0.07 ), "")
            // not used: ("simul.read_hw", po::value<double>()->default_value( 1. ), "")
            ("simul.alb_scheme", po::value<int>()->default_value( 1 ), "")
            ("simul.flooding", po::value<bool>()->default_value( true ), "")
            ("simul.alb_ice", po::value<double>()->default_value( 0.64 ), "")
            ("simul.alb_sn", po::value<double>()->default_value( 0.85 ), "")
            ("simul.I_0", po::value<double>()->default_value( 0.17 ), "")
            ("simul.Qdw", po::value<double>()->default_value( 0.5 ), "")
            ("simul.Fdw", po::value<double>()->default_value( 0. ), "")
            ("simul.newice_type", po::value<int>()->default_value( 1 ), "")
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

            ("simul.diffusivity_sss", po::value<double>()->default_value( 100. ), "") //[m^2/s]
            ("simul.diffusivity_sst", po::value<double>()->default_value( 100. ), "") //[m^2/s]

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

            ("simul.drift_limit_concentration", po::value<double>()->default_value( 0.05 ), "")
            ("simul.ERAcorr2T", po::value<std::vector<std::string>>()->multitoken()->zero_tokens()->composing(), "")
            ("simul.use_stability_drag", po::value<double>()->default_value( 0. ), "")
            ("simul.close_all_boundaries", po::value<bool>()->default_value( false ), "")

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

            ("simul.constant_significant_wave_height", po::value<double>()->default_value( 2. ), "")
            ("simul.constant_wave_mean_direction", po::value<double>()->default_value( 90. ), "")
            ("simul.constant_wave_peak_frequency", po::value<double>()->default_value( 0.05 ), "")
            ;

        return desc.add( Wim::descrWimOptions() );
#else
        ;
        return desc;
#endif
    }

} // Nextsim
