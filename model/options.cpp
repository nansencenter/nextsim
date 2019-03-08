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
    //! Describes the options to be set in a simulation (parameters, setup, dynamics, thermodynamics, numerics, ...)
    po::options_description
    descrOptions()

    {
        po::options_description desc("Options");

        desc.add_options()
            ("help,h", "Print help messages")
            ("config-file", po::value<std::string>(), "specify a configuration file")
            ("config-files", po::value<std::vector<std::string> >()->multitoken(), "specify a list of configuration files" )

        
            //-----------------------------------------------------------------------------------
            //! - Simulation parameters
            //-----------------------------------------------------------------------------------
        

            // Basics
            ("simul.time_init", po::value<std::string>()->default_value( "" ),
                "Start date/time of simulation. Formats: yyyy-mm-dd, yyyy-mm-dd HH:MM:SS; can also use 3 letter month name for 'mm' eg Mar for March")
            ("simul.duration", po::value<double>()->default_value( -1. ),
                "Length of simulation in days.")
            ("simul.timestep", po::value<int>()->default_value( 200 ), "Model timestep in seconds.")
            ("simul.thermo_timestep", po::value<int>()->default_value( 3600 ), "Thermodynamic timestep in seconds.")
            ("simul.spinup_duration", po::value<double>()->default_value( 1. ), "Spinup duration in days over which the forcing is linearly increased from 0 to its correct value.")

        
             //-----------------------------------------------------------------------------------
             //! - Debugging options
             // -----------------------------------------------------------------------------------
        
            ("debugging.bamg_verbose", po::value<int>()->default_value( 6 ),
                 "Bamg verbose mode: 0 is not verbose, 6 is very verbose")
            ("debugging.log-level", po::value<std::string>()->default_value( "info" ),
                "Nextsim printouts. Options: debug, info, warning, error")
            ("debugging.ptime_per_day", po::value<int>()->default_value( 12 ), "frequency of info printouts.")
            ("debugging.maxiteration", po::value<int>()->default_value( -1 ),
                "Stop simulation after this number of model time steps (overrides simul.duration)")
            ("debugging.check_fields", po::value<bool>()->default_value( false ),
                "call checkFields")
            ("debugging.test_proc_number", po::value<int>()->default_value( -1 ),
                "print out fields during checkFields() if on this processor number (M_rank) (do nothing if <0)")
            ("debugging.test_element_number", po::value<int>()->default_value( -1 ),
                "print out fields during checkFields() at this element number (local to M_rank = debugging.test_proc_number) (do nothing if <0)")
        
             //-----------------------------------------------------------------------------------
             //! - Numerics
             // -----------------------------------------------------------------------------------
        

            // remeshing
            ("numerics.regrid", po::value<std::string>()->default_value( "bamg" ),
                "Options for regridding: No-regridding or bamg")
            ("numerics.regrid_angle", po::value<double>()->default_value( 10. ),
                "Minimum value that any angle in an element can have.")

            // Hotfix for issue #53 - we only have pure Lagrangian now.
            // advection scheme
            // - diffusive Eulerian case where M_UM is kept as 0
            // - purely Lagrangian case where M_UM is updated with M_VT
            // - ALE case where M_UM is updated with a smoothed version of M_VT
            // ("numerics.advection_scheme", po::value<std::string>()->default_value( "Lagrangian" ), "Options: Lagrangian, ALE, Eulerian")
            // ("numerics.ALE_smoothing_step_nb", po::value<int>()->default_value( 2 ),
            //     "Number of time steps to average over when smoothing in ALE scheme. 0: pure Lagrangian; <0: pure Eulerian")

            // solver
            ("solver.ksp-type", po::value<std::string>()->default_value( "preonly" ), "")
            ("solver.pc-type", po::value<std::string>()->default_value( "cholesky" ), "")
            ("solver.mat-package-type", po::value<std::string>()->default_value( "cholmod" ), "")
            ("solver.ksp-view", po::value<bool>()->default_value( false ), "")
            ("solver.ksp-convergence-info", po::value<bool>()->default_value( true ), "")
            ("solver.ksp-reuse-prec", po::value<bool>()->default_value( false ), "")
            ("solver.ksp-monitor", po::value<bool>()->default_value( false ), "")

        
             //-----------------------------------------------------------------------------------
             //! - Setup
             //-----------------------------------------------------------------------------------


            // - Setup
            ("setup.atmosphere-type", po::value<std::string>()->default_value( "asr" ), "which atmospheric forcing?")
            ("setup.ocean-type", po::value<std::string>()->default_value( "constant" ), "which oceanic forcing?")
            ("setup.ice-type", po::value<std::string>()->default_value( "constant" ), "ice initialisation or assimilation option")
            ("setup.bathymetry-type", po::value<std::string>()->default_value( "etopo" ), "bathymetry option")
            ("setup.basal_stress-type", po::value<std::string>()->default_value( "lemieux" ), "type of basal stress model")
            ("setup.use_assimilation", po::value<bool>()->default_value( false ), "use assimilation or not")
            ("setup.dynamics-type", po::value<std::string>()->default_value( "default" ), "type of dynamics")
            ("setup.thermo-type", po::value<std::string>()->default_value( "winton" ), "which thermodynamics model")

            // mesh
            // -- input
            ("mesh.filename", po::value<std::string>()->default_value( "medium_Arctic_10km.msh" ), "name of .msh file to use")
            ("mesh.mppfile", po::value<std::string>()->default_value( "NpsNextsim.mpp" ), "name of .mpp projection file to use")
            ("mesh.partitioner", po::value<std::string>()->default_value( "metis" ),
                "mesh partitioner: chaco or metis")
            ("mesh.partitioner-fileformat", po::value<std::string>()->default_value( "binary" ),
                "Format for saving partitioned mesh. Options: ascii, binary")
            ("mesh.partitioner-space", po::value<std::string>()->default_value( "memory" ),
                "where the partitioned mesh is kept (disk/memory)")
            //not used: ("mesh.hsize", po::value<double>()->default_value( 0.01 ), "") // to be checked
            ("mesh.type", po::value<std::string>()->default_value( "from_unref" ),
                "from_unref (implies constant vertice length) or from_split (implies variable vertice length)")


            // -- moorings
            ("moorings.use_moorings", po::value<bool>()->default_value( false ), "do we use moorings (netcdf output to grid)?")
#ifdef OASIS
            ("moorings.grid_type", po::value<std::string>()->default_value( "coupled" ),
                "[coupled|regular|from_file] for using the coupling grid, a regular spaced grid, or a grid read in from the file moorings.grid_file (default: coupled)")
#else
            ("moorings.grid_type", po::value<std::string>()->default_value( "regular" ),
                "[regular|from_file] for regular spaced grid or grid read in from the file moorings.grid_file (default: regular)")
#endif
            ("moorings.snapshot", po::value<bool>()->default_value( false ), "do we output snapshots in time or do we use time-averaging?")
            ("moorings.file_length", po::value<std::string>()->default_value( "inf" ), "daily, weekly, monthly, or yearly mooring files; or inf (single file)")
            ("moorings.spacing", po::value<double>()->default_value( 10 ), "spacing between grid points (km), regular grid in the model's stereographic projection")
            ("moorings.output_timestep", po::value<double>()->default_value( 1 ), "time interval between mooring records (days)")
            ("moorings.output_time_step_units", po::value<std::string>()->default_value("days"),
                "units of moorings.output_time_step: days or time_steps")
            ("moorings.variables", po::value<std::vector<std::string>>()->multitoken()->default_value(
                        std::vector<std::string>
                            {"conc", "thick", "snow", "conc_thin", "h_thin", "hs_thin", "velocity"},
                             "conc    thick    snow    conc_thin    h_thin    hs_thin    velocity"
                    )->composing(), "list of variable names (put on separate lines in config file)")
            ("moorings.grid_file", po::value<std::string>()->default_value( "" ),
                "Grid file with locations for moorings output. It must be a netcdf file with two dimensional lat and lon")
            ("moorings.grid_latitude", po::value<std::string>()->default_value( "latitude" ), "The name of the latitude variable in the mooring_grid_file")
            ("moorings.grid_longitude", po::value<std::string>()->default_value( "longitude" ), "The name of the longitude variable in the mooring_grid_file")
            ("moorings.grid_transpose", po::value<bool>()->default_value( false ), "If true we assume the first dimension is y and the second x.")
            ("moorings.false_easting", po::value<bool>()->default_value( true ),
                "true: we output vectors relative to the output grid; false: we give their north-south components")
            ("moorings.parallel_output", po::value<bool>()->default_value( false ), "")


            // -- drifters
            ("drifters.concentration_limit", po::value<double>()->default_value( 0.15 ), "lower concentration limit for ignoring a drifter")
            ("drifters.use_iabp_drifters", po::value<bool>()->default_value( false), "Use IABP drifters or not")
            ("drifters.iabp_drifters_output_time_step", po::value<double>()->default_value( .5 ),
                 "interval between IABP drifter outputs (days): 0.5/2n, n=1,2,... down to timestep")

            ("drifters.use_osisaf_drifters", po::value<bool>()->default_value( false ), "Use OSISAF drifters?")
            ("drifters.osisaf_drifters_output_time_step", po::value<double>()->default_value( 2. ),
                 "interval between OSISAF drifter outputs (days): 2.0/n, n=1,2,... down to timestep")
            ("drifters.use_refined_osisaf_grid", po::value<bool>()->default_value( false ),
                "true: if using OSISAF drifters, use grid refined by a factor of 9, so averaged model results can be compared to the data; false: use same grid as OSISAF drift dataset")

            ("drifters.use_equally_spaced_drifters", po::value<bool>()->default_value( false), "use equally spaced drifters?")
            ("drifters.equally_spaced_drifters_output_time_step", po::value<double>()->default_value( .5 ),
                 "interval between equally-spaced drifter outputs (days): integer multiple of timestep")
            ("drifters.spacing", po::value<double>()->default_value( 10 ), "spacing of equally spaced drifters in km")

            ("drifters.use_rgps_drifters", po::value<bool>()->default_value( false), "")
            ("drifters.rgps_drifters_output_time_step", po::value<double>()->default_value( 0.5 ),
                 "interval between RGPS drifter outputs (days): integer multiple of timestep")
            ("drifters.RGPS_time_init", po::value<std::string>()->default_value( "2007-12-01" ),
                "time to init RGPS drifters: date format yyyy-mm-dd or yyyy-mmm-dd (eg 2008-Mar-05); can also add time with HH:MM:SS (eg 2008-Mar-05 00:00:00)")

             ("drifters.use_sidfex_drifters", po::value<bool>()->default_value( false), "are we using SIDFEX drifters")
             ("drifters.sidfex_drifters_output_time_step", po::value<double>()->default_value( 0.5 ),
                  "interval between SIDFEX drifter outputs (days): integer multiple of timestep")
             ("drifters.sidfex_filename", po::value<std::string>()->default_value( "" ), "text file with initial buoy positions")

            // - Restart
            // -- inputs
            ("restart.start_from_restart", po::value<bool>()->default_value( false ),
                "are we starting from a restart file?")
            ("restart.input_path", po::value<std::string>()->default_value( "" ),
                    "where to find restart files")
            ("restart.input_filename", po::value<std::string>()->default_value( "" ),
                "if we are starting from a restart file, the field files' names will be [restart.input_path]/[restart.filename].[bin,dat]")
            ("restart.reset_time_counter", po::value<bool>()->default_value( false ),
                "true: simulation starts at simul.time_init eg for forecast; false: simulation starts at simul.time_init+pcpt*simul.timestep eg to restart interrupted simulation")

            // -- outputs
            ("restart.write_final_restart", po::value<bool>()->default_value( false ),
                "To write a out a restart file at the end of the run")
            ("restart.write_interval_restart", po::value<bool>()->default_value( false ),
                "To write a out a restart file at an interval given by restart.restart_output_interval")
            ("restart.write_initial_restart", po::value<bool>()->default_value( false ),
                "To write a out a restart file at the start of the run")
            ("restart.output_interval", po::value<double>()->default_value( 15 ), "How often to write restarts (with restart.write_interval_restart), in days or time_steps, according to restart.output_interval_units")

            // -- restart debugging
            ("restart.output_interval_units", po::value<std::string>()->default_value("days"),
                "Units of restart.output_interval: days or time_steps")
            ("restart.restart_at_rest", po::value<bool>()->default_value( false ),
                "reset ice velocity to zero if starting from restart")
            ("restart.output_before_regrid", po::value<bool>()->default_value( false ),
                "if true, export results before regrid")
            ("restart.output_after_regrid", po::value<bool>()->default_value( false ),
                "if true, export results after regrid")

            // -- general outputs
            ("output.output_per_day", po::value<int>()->default_value( 4 ), "")
            ("output.logfile", po::value<std::string>()->default_value( "" ), "")
            ("output.save_forcing_fields", po::value<bool>()->default_value( false ), "")
            ("output.save_diagnostics", po::value<bool>()->default_value( false ), "")
#if 0
            //TODO issue193 uncomment these lines to set export variables using config file (finish another time)
            ("output.variables", po::value<std::vector<std::string>>()->multitoken()
                ->default_value( std::vector<std::string>{}, "")->composing(),
                "list of variable names (put on separate lines in config file)")
#endif

            // --exporter
            ("output.datetime_in_filename", po::value<bool>()->default_value( false ),
                "filename outputs are eg [mesh,field]_20180101T000000Z.[bin,dat]")
            ("output.exporter_path", po::value<std::string>()->default_value( "nextsim_outputs" ),
                "Path where results should be exported")
            ("output.exporter_precision", po::value<std::string>()->default_value("float"),
                    "float (default) or double (almost only for testing)")

        
             //-----------------------------------------------------------------------------------
             //! - Ideal configurations (constant initial bathymetry, ice, atmosphere and ocean)
             //-----------------------------------------------------------------------------------
        

            // - if setup.bathymetry-type = constant
            ("ideal_simul.constant_bathymetry", po::value<double>()->default_value( 200. ), "")

            // - if setup.ice-type = constant
            ("ideal_simul.init_thickness", po::value<double>()->default_value( 1.0 ), "")
            ("ideal_simul.init_concentration", po::value<double>()->default_value( 1.0 ), "")
            ("ideal_simul.init_thin_conc", po::value<double>()->default_value( 0. ), "")
            ("ideal_simul.init_snow_thickness", po::value<double>()->default_value( 0. ), "")
            ("ideal_simul.init_SST_limit", po::value<double>()->default_value( 2. ), "")

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
            ("ideal_simul.constant_wind_v", po::value<double>()->default_value( 0. ), "")

            // - if ocean-type=constant
            ("ideal_simul.constant_ocean_u", po::value<double>()->default_value( 0. ), "")
            ("ideal_simul.constant_ocean_v", po::value<double>()->default_value( 0. ), "")
            ("ideal_simul.constant_ssh", po::value<double>()->default_value( 0. ), "")
            // -----------------------------------------------------------------------------------


        
             //-----------------------------------------------------------------------------------
             //! - Dynamics
             //-----------------------------------------------------------------------------------
        

            // - Internal stresses
            ("dynamics.alea_factor", po::value<double>()->default_value( 0.0 ), "")     // Fraction of C_fix that will be added to C_fix as some alea on the cohesion
            ("dynamics.young", po::value<double>()->default_value( 5.9605e+08 ), "Pa")  // 5.3645e+09 gives an elastic wave speed of 1500 m/s and td0 = 6.666 s for resolution of 10 km
                                                                                        // 2.3842e+09 gives an elastic wave speed of 1000 m/s and td0 = 10 s for of 10 km
                                                                                        // 5.9605e+08 gives an elastic wave speed of 500 m/s and td0 = 20 s for resolution of 10 km
            ("dynamics.C_lab", po::value<double>()->default_value( 6.8465e+6 ), "Pa")   // Cohesion value at the lab scale (10^6 Pa is the order of magnitude determined by Schulson).
            ("dynamics.nu0", po::value<double>()->default_value( 0.3 ), "")
            ("dynamics.tan_phi", po::value<double>()->default_value( 0.7 ), "")
            ("dynamics.ridging_exponent", po::value<double>()->default_value( -20. ), "")

            // - C,h limits for where to use MEB rheology and where to use the Laplacian free drift thing
            ("dynamics.min_h", po::value<double>()->default_value( 0.05 ), "")
            ("dynamics.min_c", po::value<double>()->default_value( 0.01 ), "")

            // - Ratio of ridged ice cohesion and compressive strength compared to level ice (1. does nothing)
            ("dynamics.ridge_to_normal_cohesion_ratio", po::value<double>()->default_value( 1. ), "")

            ("dynamics.use_temperature_dependent_healing", po::value<bool>()->default_value( false ), "")
            ("dynamics.time_relaxation_damage", po::value<double>()->default_value( 25. ), "days")
            ("dynamics.deltaT_relaxation_damage", po::value<double>()->default_value( 20. ), "Kelvin")
            ("dynamics.undamaged_time_relaxation_sigma", po::value<double>()->default_value( 1e7 ), "seconds")
            // from V. Dansereau et al.: A Maxwell elasto-brittle rheology for sea ice modelling
            ("dynamics.exponent_relaxation_sigma", po::value<double>()->default_value( 5. ), "")
                // from V. Dansereau et al.: A Maxwell elasto-brittle rheology for sea ice modelling

            // - Water and air drag parameterizations
            ("dynamics.ERAi_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
            ("dynamics.ECMWF_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
            ("dynamics.ASR_quad_drag_coef_air", po::value<double>()->default_value( 0.0049 ), "")
            ("dynamics.CFSR_quad_drag_coef_air", po::value<double>()->default_value( 0.0023 ), "")
            // Updated value, based on comparison with OSISAF drift in the free drift case
            ("dynamics.lin_drag_coef_air", po::value<double>()->default_value( 0. ), "")
            ("dynamics.quad_drag_coef_water", po::value<double>()->default_value( 0.0055 ), "")
            ("dynamics.lin_drag_coef_water", po::value<double>()->default_value( 0. ), "")

            // - Coriolis parameters
            ("dynamics.use_coriolis", po::value<bool>()->default_value( true ), "")
            ("dynamics.oceanic_turning_angle", po::value<double>()->default_value( 25. ), "")

            // - Basal stress parameters
            ("dynamics.Lemieux_basal_k1", po::value<double>()->default_value( 10. ), "")
            ("dynamics.Lemieux_basal_k2", po::value<double>()->default_value( 15. ), "")
            ("dynamics.Lemieux_basal_Cb", po::value<double>()->default_value( 20. ), "")
            ("dynamics.Lemieux_basal_u_0", po::value<double>()->default_value( 5e-5 ), "")
            ("dynamics.Lemieux_basal_u_crit", po::value<double>()->default_value( 5e-4 ), "")
        
            // - Damage equation discretization
            //   disc_scheme is either : explicit, implicit, recursive
            //   td_type is either : fixed or damage_dependent
            ("damage.disc_scheme", po::value<std::string>()->default_value( "explicit" ), "which discretization scheme for the damage equation?")
            ("damage.td_type", po::value<std::string>()->default_value( "fixed" ), "is the char. time for damage fixed or damage dependent?")
        
        
             //-----------------------------------------------------------------------------------
             //! - Thermodynamics
             //-----------------------------------------------------------------------------------
        
            ("thermo.use_thermo_forcing", po::value<bool>()->default_value( true ), "")
            ("thermo.Qio-type", po::value<std::string>()->default_value( "basic" ), "")
            ("thermo.freezingpoint-type", po::value<std::string>()->default_value( "linear" ), "How to calculate the freezing point of sea water, either linear or non-linear")
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

            // -- diffusivity
            ("thermo.diffusivity_sss", po::value<double>()->default_value( 100. ), "") //[m^2/s]
            ("thermo.diffusivity_sst", po::value<double>()->default_value( 100. ), "") //[m^2/s]

            // -- relaxation of slab ocean to ocean forcing
            ("thermo.ocean_nudge_timeT", po::value<double>()->default_value( 30*days_in_sec),
                "relaxation time of slab ocean temperature to ocean forcing")
            ("thermo.ocean_nudge_timeS", po::value<double>()->default_value( 30*days_in_sec),
                "relaxation time of slab ocean salinity to ocean forcing")

            // -- relating to thermodynamic forcing
            ("thermo.use_parameterised_long_wave_radiation", po::value<bool>()->default_value(false),
                "True: use total cloud cover parameterisation of long wave incoming radiation - only works if dataset has QLW_IN. False: use forcing from atmospheric datasets - only works if dataset has TCC")

#ifdef AEROBULK
            ("thermo.ocean_bulk_formula", po::value<std::string>()->default_value( "nextsim" ), "Bulk formula to calculate ocean-atmosphere fluxes [ nextsim (default) | coare | coare3.5 | ncar | ecmwf ]")
#endif
        
             //-----------------------------------------------------------------------------------
             //! - Nesting
             // -----------------------------------------------------------------------------------

            ("nesting.use_nesting", po::value<bool>()->default_value( false ),
                "Use nesting at/near boundaries")
            ("nesting.use_ocean_nesting", po::value<bool>()->default_value( false ),
                "Use slab ocean from outer neXtSIM model as forcing")
            ("nesting.outer_mesh",po::value<std::string>()->default_value( "" ),
                "Filenames are nesting_[outer_mesh]_[yyyymmdd].nc")
            ("nesting.inner_mesh",po::value<std::string>()->default_value( "" ),
                "Grid filename is nsting_grid_[inner_mesh].nc")
            ("nesting.method", po::value<std::string>()->default_value( "nudging" ),
                "Options: nudging")
            ("nesting.nudge_timescale", po::value<double>()->default_value((1./2.)*days_in_sec),
                "relaxation timescale for nudging at boundary")
            ("nesting.nudge_function", po::value<std::string>()->default_value( "exponential" ),
                "Functional form for nudging frequency as a function of distance to boundary. Options: exponential, linear. Depends on nudge_length_scale")
            ("nesting.nudge_lengthscale", po::value<double>()->default_value(10.),
                "Length scale (units = mesh resolution) of the nudging function")
            ("nesting.nest_dynamic_vars", po::value<bool>()->default_value( false ), "")

        
             //-----------------------------------------------------------------------------------
             //! - Forecasting
             //-----------------------------------------------------------------------------------
        
            ("forecast.air_temperature_correction", po::value<double>()->default_value( 0. ),
                "for use in BADA (Bias-Aware-Data-Assimilation)")
            ("forecast.true_forecast", po::value<bool>()->default_value( false ),
                "if false, use <<analysis/best est>> - forecast that started on the same day as the model time; else use the one that began on simul.time_init, or  forecast.time_init_atm_fc/forecast.time_init_ocean_fc (if provided)")
            ("forecast.time_init_atm_fc", po::value<std::string>()->default_value( "" ),
                "if(forecast.true_forecast), get atmospheric forecast starting from this date as opposed to simul.time_init (eg if usual one is absent)")
            ("forecast.time_init_ocean_fc", po::value<std::string>()->default_value( "" ),
                "if(forecast.true_forecast), get ocean forecast starting from this date as opposed to simul.time_init (eg if usual one is absent)")

        
             //-----------------------------------------------------------------------------------
             //! - Modules
             //-----------------------------------------------------------------------------------
        
#if defined(OASIS)
            ("coupler.component_name", po::value<std::string>()->default_value( "nxtsim" ), "Component name (6 characters) same as in the namcouple")
            ("coupler.timestep", po::value<int>()->default_value( 3600 ), "Coupling time step")
            ("coupler.exchange_grid_file", po::value<std::string>()->default_value( "coupler/NEMO.nc" ), "File containing neccesary grid information for coupling.")
            // ("coupler.with_ocean", po::value<bool>()->default_value( false ), "Do we couple with an ocean model?")
            // ("coupler.with_waves", po::value<bool>()->default_value( false ), "Do we couple with a wave model?")
            // ("coupler.with_atm", po::value<bool>()->default_value( false ), "Do we couple with an atmospheric model?")
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
