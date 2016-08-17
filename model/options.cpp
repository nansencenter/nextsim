/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   options.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue Jul 14 13:23:45 2015
 */


#include <boost/program_options.hpp>
#include <constants.hpp>
#include <wimoptions.hpp>

namespace po = boost::program_options;

// namespace Wim
// {
//     po::options_description descrWimOptions();
// }

double const days_in_sec = 24.0*3600.0;
namespace Nextsim
{
    //po::options_description _test = descrWimOptions();

    po::options_description
    descrOptions()
    {
        po::options_description desc("Options");

        desc.add_options()
            ("help,h", "Print help messages")
            ("config-file", po::value<std::string>(), "specify a configuration file")
            ("config-files", po::value<std::vector<std::string> >()->multitoken(), "specify a list of configuration files" )

            ("solver.ksp-type", po::value<std::string>()->default_value( "preonly" ), "")
            ("solver.pc-type", po::value<std::string>()->default_value( "cholesky" ), "")
            ("solver.mat-package-type", po::value<std::string>()->default_value( "cholmod" ), "")
            ("solver.ksp-view", po::value<bool>()->default_value( false ), "")
            ("solver.ksp-convergence-info", po::value<bool>()->default_value( true ), "")

            ("setup.atmosphere-type", po::value<std::string>()->default_value( "asr" ), "")
            ("setup.ocean-type", po::value<std::string>()->default_value( "constant" ), "")
            ("setup.ice-type", po::value<std::string>()->default_value( "constant" ), "")
            ("setup.bathymetry-type", po::value<std::string>()->default_value( "etopo" ), "")
            ("setup.drifter-type", po::value<std::string>()->default_value( "none" ), "")
            ("setup.use_restart", po::value<bool>()->default_value( false ), "")
            ("setup.write_restart", po::value<bool>()->default_value( false ), "")
            ("setup.restart_time_step", po::value<double>()->default_value( 20 ), "days")
            ("setup.step_nb", po::value<int>()->default_value( 0 ), "")

            ("simul.mesh_filename", po::value<std::string>()->default_value( "bigarctic10km.msh" ), "")
            ("simul.hsize", po::value<double>()->default_value( 0.01 ), "")
            ("simul.verbose", po::value<int>()->default_value( 7 ), "")
            ("simul.log-level", po::value<std::string>()->default_value( "info" ), "")
            ("simul.time_init", po::value<std::string>()->default_value( "2008-Mar-05" ), "")
            ("simul.duration", po::value<double>()->default_value( 1. ), "")
            ("simul.spinup_duration", po::value<double>()->default_value( 1. ), "")
            ("simul.output_per_day", po::value<int>()->default_value( 24 ), "")
            ("simul.diagnostics_frequency", po::value<int>()->default_value( 1 ), "")
            ("simul.save_forcing_field", po::value<bool>()->default_value( false ), "")
            ("simul.forecast", po::value<bool>()->default_value( false ), "")
            ("simul.use_ship", po::value<bool>()->default_value( false ), "")
            ("simul.use_moorings", po::value<bool>()->default_value( false ), "")
            ("simul.mooring_names", po::value<std::vector<std::string>>()->multitoken()->zero_tokens()->composing(), "")
            ("simul.mooring_spacing", po::value<double>()->default_value( 10 ), "km")
            ("simul.mooring_output_timestep", po::value<double>()->default_value( 1 ), "days")
            ("simul.slim_bin.use", po::value<bool>()->default_value( false ), "")
            ("simul.slim_bin.names", po::value<std::vector<std::string>>()->multitoken()->zero_tokens()->composing(), "")
            ("simul.slim_bin.output_timestep", po::value<double>()->default_value( std::numeric_limits<double>::infinity() ), "")
            ("simul.drifters_output_timestep", po::value<double>()->default_value( std::numeric_limits<double>::infinity() ), "")
            ("simul.drifter_climit", po::value<double>()->default_value( 0.15 ), "")
            ("simul.reduce_domain_to_sea_ice", po::value<bool>()->default_value( true ), "")
            ("simul.transport_scheme", po::value<std::string>()->default_value( "lagrangian" ), "")
            ("simul.regrid", po::value<std::string>()->default_value( "bamg" ), "No-regridding or bamg")
            ("simul.regrid_output_flag", po::value<bool>()->default_value( false ), "")
            ("simul.regrid_angle", po::value<double>()->default_value( 10. ), "")
            ("simul.interp_with_cavities", po::value<bool>()->default_value( true ), "")
            ("simul.interp_forcing_every_timestep", po::value<bool>()->default_value( false ), "")
            ("simul.timestep", po::value<double>()->default_value( 200. ), "")
            ("simul.scale_coef", po::value<double>()->default_value( 0.1 ), "")
            ("simul.alea_factor", po::value<double>()->default_value( 0. ), "")
            ("simul.ridging_exponent", po::value<double>()->default_value( -40. ), "")
            ("simul.fixed_cohesion", po::value<bool>()->default_value( true ), "")
            ("simul.divergence_min", po::value<double>()->default_value( 0.01 ), "")
            ("simul.compression_factor", po::value<double>()->default_value( 12000. ), "")
            ("simul.exponent_compression_factor", po::value<double>()->default_value( 2. ), "")
            ("simul.init_damage", po::value<std::string>()->default_value( "No-damage" ), "")
            ("simul.time_relaxation_damage", po::value<double>()->default_value( 28. ), "")
            ("simul.deltaT_relaxation_damage", po::value<double>()->default_value( 40. ), "")
            ("simul.young", po::value<double>()->default_value( 9000e+6 ), "")
            ("simul.cfix", po::value<double>()->default_value( 40e+3 ), "")
            ("simul.nu0", po::value<double>()->default_value( 0.3 ), "")
            ("simul.tan_phi", po::value<double>()->default_value( 0.7 ), "")
            ("simul.tract_coef", po::value<double>()->default_value( 5./6 ), "")
            ("simul.compr_strength", po::value<double>()->default_value( 750e+3 ), "")
            ("simul.rheology_scaling", po::value<bool>()->default_value( false ), "")
            ("simul.rheology_scaling_ref_dx", po::value<double>()->default_value( 10. ), "")
            ("simul.rheology_scaling_factor", po::value<double>()->default_value( -0.0638 ), "")
            ("simul.Lemieux_basal_k1", po::value<double>()->default_value( 10. ), "")
            ("simul.Lemieux_basal_k2", po::value<double>()->default_value( 15. ), "")
            ("simul.Lemieux_drag_coef_air", po::value<double>()->default_value( 0.0012 ), "")
            ("simul.Lemieux_basal_Cb", po::value<double>()->default_value( 20. ), "")
            ("simul.Lemieux_basal_u_0", po::value<double>()->default_value( 5e-5 ), "")
            ("simul.Lemieux_basal_u_crit", po::value<double>()->default_value( 5e-4 ), "")
            ("simul.albedoW", po::value<double>()->default_value( 0.07 ), "")
            ("simul.read_hw", po::value<double>()->default_value( 1. ), "")
            ("simul.hw", po::value<double>()->default_value( 9. ), "")
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
            ("simul.c_thin_max", po::value<double>()->default_value( 1. ), "")
            ("simul.ERAcorr2T", po::value<std::vector<std::string>>()->multitoken()->zero_tokens()->composing(), "")
            ("simul.gravity", po::value<double>()->default_value( 9.80616 ), "")
            ("simul.omega", po::value<double>()->default_value( 7.292e-5 ), "")
            ("simul.rho_air", po::value<double>()->default_value( 1.3 ), "")
            ("simul.use_stability_drag", po::value<double>()->default_value( 0. ), "")
            ("simul.drag_ice_t", po::value<double>()->default_value( 1.3e-3 ), "")
            ("simul.drag_ocean_u", po::value<double>()->default_value( 1.1e-3 ), "")
            ("simul.drag_ocean_t", po::value<double>()->default_value( 0.83e-3 ), "")
            ("simul.drag_ocean_q", po::value<double>()->default_value( 1.5e-3 ), "")
            ("simul.use_ocean_forcing", po::value<bool>()->default_value( false ), "")
            ("simul.ASR_quad_drag_coef_air", po::value<double>()->default_value( 0.0076 ), "")
            ("simul.ERAi_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
            ("simul.ECMWF_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
            ("simul.ASRFINAL_quad_drag_coef_air", po::value<double>()->default_value( 0.0037 ), "")
            ("simul.lin_drag_coef_air", po::value<double>()->default_value( 0. ), "")
            ("simul.rho_water", po::value<double>()->default_value( 1025. ), "")
            ("simul.quad_drag_coef_water", po::value<double>()->default_value( 0.0055 ), "")
            ("simul.lin_drag_coef_water", po::value<double>()->default_value( 0. ), "")
            ("simul.oceanic_turning_angle", po::value<double>()->default_value( 25. ), "")
            ("simul.drift_limit_concentration", po::value<double>()->default_value( 0.05 ), "")
            ("simul.min_h", po::value<double>()->default_value( 0. ), "")
            ("simul.min_c", po::value<double>()->default_value( 0. ), "")
            ("simul.ridge_h", po::value<double>()->default_value( 3. ), "")
            ("simul.scale_coef", po::value<double>()->default_value( 0.1 ), "")
            ("simul.close_all_boundaries", po::value<bool>()->default_value( false ), "")
            ("simul.init_thickness", po::value<double>()->default_value( 1.0 ), "")
            ("simul.init_concentration", po::value<double>()->default_value( 1.0 ), "")
            ("simul.init_thin_max_thickness", po::value<double>()->default_value( 0.1 ), "")
            ("simul.init_snow_thickness", po::value<double>()->default_value( 0. ), "")
            ("simul.init_snow_thin_thickness", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_wind_u", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_wind_v", po::value<double>()->default_value( -10. ), "")
            ("simul.constant_tair", po::value<double>()->default_value( -25. ), "")
            ("simul.constant_dair", po::value<double>()->default_value( -1. ), "")
            ("simul.constant_mixrat", po::value<double>()->default_value( 0.001 ), "")
            ("simul.constant_mslp", po::value<double>()->default_value( 1013e2 ), "")
            ("simul.constant_Qsw_in", po::value<double>()->default_value( 50. ), "")
            ("simul.constant_Qlw_in", po::value<double>()->default_value( 250. ), "")
            ("simul.constant_precip", po::value<double>()->default_value( 1e-5 ), "")
            ("simul.constant_snowfr", po::value<double>()->default_value( 0.9 ), "")
            ("simul.ocean_nudge", po::value<bool>()->default_value( false ), "")
            ("simul.ocean_nudge_timeT", po::value<double>()->default_value( 30*days_in_sec), "")
            ("simul.ocean_nudge_timeS", po::value<double>()->default_value( 30*days_in_sec), "")
            ("simul.constant_Qdw", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_Fdw", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_ocean_u", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_ocean_v", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_ssh", po::value<double>()->default_value( 0. ), "")
            ("simul.constant_ocean_temp", po::value<double>()->default_value( -1.8 ), "")
            ("simul.constant_ocean_salt", po::value<double>()->default_value( -1.8/physical::mu ), "")
            ("simul.constant_mld", po::value<double>()->default_value( 9. ), "")
            ("simul.constant_bathymetry", po::value<double>()->default_value( 200. ), "")
            ("simul.use_thermo_forcing", po::value<bool>()->default_value( true ), "")
            ("simul.use_wim", po::value<bool>()->default_value( false ), "")
            ("simul.wim_grid", po::value<bool>()->default_value( false ), "")
            ("simul.maxiteration", po::value<int>()->default_value( 1e+8 ), "")
            ;
        return desc.add( Wim::descrWimOptions() );
    }

} // Nextsim
