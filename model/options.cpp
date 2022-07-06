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
            ("simul.spinup_duration", po::value<double>()->default_value( 1. ), "Spinup duration in days over which the forcing is linearly increased from 0 to its correct value.")

             //-----------------------------------------------------------------------------------
             //! - Debugging options
             // -----------------------------------------------------------------------------------

            ("debugging.bamg_verbose", po::value<int>()->default_value( 0 ),
                 "Bamg verbose mode: 0 is not verbose, 6 is very verbose")
            ("debugging.gmsh_verbose", po::value<int>()->default_value( 0 ),
                 "Gmsh verbose mode: 0: silent except for fatal errors, 1: +errors, 2: +warnings, 3: +direct, 4: +information, 5: +status, 99: +debug")
            ("debugging.log-level", po::value<std::string>()->default_value( "info" ),
                "Nextsim printouts. Options: debug, info, warning, error")
            ("debugging.log-all", po::value<bool>()->default_value( false ),
                "Whether printouts from debugging.log-level are to be done on all processors [true|false(default)].")
            ("debugging.ptime_percent", po::value<int>()->default_value( 5 ), "frequency of info printouts in % of runtime.")
            ("debugging.maxiteration", po::value<int>()->default_value( -1 ),
                "Stop simulation after this number of model time steps (overrides simul.duration)")
            ("debugging.check_fields", po::value<bool>()->default_value( false ),
                "call checkFields")
            ("debugging.test_proc_number", po::value<int>()->default_value( -1 ),
                "print out fields during checkFields() if on this processor number (M_rank) (do nothing if <0)")
            ("debugging.test_element_number", po::value<int>()->default_value( -1 ),
                "print out fields during checkFields() at this element number (local to M_rank = debugging.test_proc_number) (do nothing if <0)")
            ("debugging.check_velocity_fields", po::value<bool>()->default_value( false ),
                "If check_velocity_fields is true: find outlier nodes with extreme velocities printed to DEBUG")
            ("debugging.check_fields_fast", po::value<bool>()->default_value( true ),
                "Do a quick sanity check on select fields. Export binary files and stop model if the check fails.")

             //-----------------------------------------------------------------------------------
             //! - Numerics
             // -----------------------------------------------------------------------------------


            // remeshing
            ("numerics.regrid", po::value<std::string>()->default_value( "bamg" ),
                "Options for regridding: No-regridding or bamg")
            ("numerics.regrid_angle", po::value<double>()->default_value( 10. ),
                "Minimum value that any angle in an element can have.")
            ("numerics.nit_ow", po::value<int>()->default_value( 50. ),
                "Number of iterations taken to smooth velocity into open water")

            // Hotfix for issue #53 - we only have pure Lagrangian now.
            // advection scheme
            // - diffusive Eulerian case where M_UM is kept as 0
            // - purely Lagrangian case where M_UM is updated with M_VT
            // - ALE case where M_UM is updated with a smoothed version of M_VT
            // ("numerics.advection_scheme", po::value<std::string>()->default_value( "Lagrangian" ), "Options: Lagrangian, ALE, Eulerian")
            // ("numerics.ALE_smoothing_step_nb", po::value<int>()->default_value( 2 ),
            //     "Number of time steps to average over when smoothing in ALE scheme. 0: pure Lagrangian; <0: pure Eulerian")

            //-----------------------------------------------------------------------------------
            //! - Setup
            //-----------------------------------------------------------------------------------


            // - Setup
            ("setup.atmosphere-type", po::value<std::string>()->default_value( "asr" ), "which atmospheric forcing?")
            ("setup.ocean-type", po::value<std::string>()->default_value( "constant" ), "which oceanic forcing?")
            ("setup.ice-type", po::value<std::string>()->default_value( "constant" ), "ice initialisation or assimilation option")
            ("setup.bathymetry-type", po::value<std::string>()->default_value( "etopo" ), "bathymetry option")
            ("setup.bathymetry-file", po::value<std::string>()->default_value( "ETOPO_Arctic_2arcmin.nc" ), "Bathymetry file for basal stress calculations (ETOPO_Arctic_2arcmin.nc)")
            ("setup.basal_stress-type", po::value<std::string>()->default_value( "lemieux" ), "type of basal stress model")
            ("setup.use_assimilation", po::value<bool>()->default_value( false ), "use assimilation or not")
            ("setup.dynamics-type", po::value<std::string>()->default_value( "bbm" ), "type of dynamics [ bbm | no_motion | evp | mevp | free_drift ] ")
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
            ("moorings.use_conservative_remapping", po::value<bool>()->default_value( false ),
                 "Use the conservative remapping scheme to interpolate onto the moorings grid. Requires mooring.grid_type=from_file, mesh.type=from_split and .msh and .nc files from the mkCplMesh.m script in nextsim-env (default: false)")
#endif
            ("moorings.snapshot", po::value<bool>()->default_value( false ), "do we output snapshots in time or do we use time-averaging?")
            ("moorings.file_length", po::value<std::string>()->default_value( "inf" ), "daily, weekly, monthly, or yearly mooring files; or inf (single file)")
            ("moorings.spacing", po::value<double>()->default_value( 10 ), "spacing between grid points (km), regular grid in the model's stereographic projection")
            ("moorings.output_timestep", po::value<double>()->default_value( 1 ), "time interval between mooring records (days)")
            ("moorings.output_time_step_units", po::value<std::string>()->default_value("days"),
                "units of moorings.output_time_step: days or time_steps")
            ("moorings.variables", po::value<std::vector<std::string>>()->multitoken()->default_value(
                        std::vector<std::string>
                            {"conc", "thick", "snow", "conc_young", "h_young", "hs_young", "velocity"},
                             "conc    thick    snow    conc_young    h_young    hs_young    velocity"
                    )->composing(), "list of variable names (put on separate lines in config file)")
            ("moorings.grid_file", po::value<std::string>()->default_value( "" ),
                "Grid file with locations for moorings output. It must be a netcdf file with two dimensional lat and lon")
            ("moorings.grid_latitude", po::value<std::string>()->default_value( "latitude" ), "The name of the latitude variable in the mooring_grid_file")
            ("moorings.grid_longitude", po::value<std::string>()->default_value( "longitude" ), "The name of the longitude variable in the mooring_grid_file")
            ("moorings.grid_transpose", po::value<bool>()->default_value( false ), "If true we assume the first dimension is x and the second y (non-standard ordering).")
            ("moorings.false_easting", po::value<bool>()->default_value( true ),
                "true: we output vectors relative to the output grid; false: we give their north-south components. NB only implemented for grid_type=regular")
            ("moorings.parallel_output", po::value<bool>()->default_value( false ), "")


            // -- drifters
            ("drifters.concentration_limit", po::value<double>()->default_value( 0.15 ), "lower concentration limit for ignoring a drifter")
            ("drifters.use_iabp_drifters", po::value<bool>()->default_value( false), "Use IABP drifters or not")
            ("drifters.iabp_drifters_output_time_step", po::value<double>()->default_value( .5 ),
                 "interval between IABP drifter outputs (days): 0.5/2n, n=1,2,... down to timestep")
            ("drifters.iabp_ignore_restart", po::value<bool>()->default_value( false),
                "do not load IABP drifters from restart even if present")

            ("drifters.use_osisaf_drifters", po::value<bool>()->default_value( false ), "Use OSISAF drifters?")
            ("drifters.osisaf_drifters_output_time_step", po::value<double>()->default_value( 2. ),
                 "interval between OSISAF drifter outputs (days): 2.0/n, n=1,2,... down to timestep")
            ("drifters.use_refined_osisaf_grid", po::value<bool>()->default_value( false ),
                "true: if using OSISAF drifters, use grid refined by a factor of 9, so averaged model results can be compared to the data; false: use same grid as OSISAF drift dataset")

            ("drifters.use_equally_spaced_drifters", po::value<bool>()->default_value( false), "use equally spaced drifters?")
            ("drifters.equally_spaced_drifters_output_time_step", po::value<double>()->default_value( .5 ),
                 "interval between equally-spaced drifter outputs (days): integer multiple of timestep")
            ("drifters.spacing", po::value<double>()->default_value( 10 ), "spacing of equally spaced drifters in km")
            ("drifters.equally_spaced_ignore_restart", po::value<bool>()->default_value( false),
                "do not load equally spaced drifters from restart even if present")

            ("drifters.use_rgps_drifters", po::value<bool>()->default_value( false), "")
            ("drifters.rgps_drifters_output_time_step", po::value<double>()->default_value( 0.5 ),
                 "interval between RGPS drifter outputs (days): integer multiple of timestep")
            ("drifters.RGPS_time_init", po::value<std::string>()->default_value( "2007-12-01" ),
                "time to init RGPS drifters: date format yyyy-mm-dd or yyyy-mmm-dd (eg 2008-Mar-05); can also add time with HH:MM:SS (eg 2008-Mar-05 00:00:00)")

             ("drifters.use_sidfex_drifters", po::value<bool>()->default_value( false),
                "are we using SIDFEX drifters")
             ("drifters.sidfex_drifters_output_time_step", po::value<double>()->default_value( 0.5 ),
                  "interval between SIDFEX drifter outputs (days): integer multiple of timestep")
             ("drifters.sidfex_filename", po::value<std::string>()->default_value( "" ), "text file with initial buoy positions")
             ("drifters.sidfex_time_init", po::value<std::string>()->default_value( "" ),
                  "time to init SIDFEx drifters: date format yyyy-mm-dd or yyyy-mmm-dd (eg 2008-Mar-05); can also add time with HH:MM:SS (eg 2008-Mar-05 00:00:00)")
             ("drifters.sidfex_ignore_restart", po::value<bool>()->default_value( false),
                "do not load SIDFEx buoys from restart even if present")

            // - Restart
            // -- inputs
            ("restart.start_from_restart", po::value<bool>()->default_value( false ),
                "are we starting from a restart file?")
            ("restart.check_restart", po::value<bool>()->default_value( false ),
                "check restart file at init time?")
            ("restart.input_path", po::value<std::string>()->default_value( "" ),
                    "where to find restart files")
            ("restart.basename", po::value<std::string>()->default_value( "" ),
                "The base of a restart file name. If we are starting from restart files, the files' names will be (restart.input_path)/{field|mesh}_(restart.basename).{bin,dat}")
            ("restart.type", po::value<std::string>()->default_value( "extend" ),
                "Restart type: [extend|continue]. Extend (default): M_time_init is taken as the time inside the restart file and simul.duration is added to that. Continue: M_time_init is read from the configuration file (simul.time_init) and duration is added to that.")

            // -- outputs
            ("restart.write_final_restart", po::value<bool>()->default_value( false ),
                "To write a out a restart file at the end of the run")
            ("restart.write_interval_restart", po::value<bool>()->default_value( false ),
                "To write a out a restart file at an interval given by restart.restart_output_interval")
            ("restart.write_initial_restart", po::value<bool>()->default_value( false ),
                "To write a out a restart file at the start of the run")
            ("restart.output_interval", po::value<double>()->default_value( 15 ), "How often to write restarts (with restart.write_interval_restart), in days or time_steps, according to restart.output_interval_units")
            ("restart.datetime_in_filename", po::value<bool>()->default_value( true ),
                "filename outputs are eg [mesh,field]_20180101T000000Z.[bin,dat]")

            // -- restart debugging
            ("restart.output_interval_units", po::value<std::string>()->default_value("days"),
                "Units of restart.output_interval: days or time_steps")
            ("restart.restart_at_rest", po::value<bool>()->default_value( false ),
                "reset ice velocity to zero if starting from restart")
            ("restart.write_restart_before_regrid", po::value<bool>()->default_value( false ),
                "if true, write restart before regrid")
            ("restart.write_restart_after_regrid", po::value<bool>()->default_value( false ),
                "if true, write restart after regrid")

            // -- general outputs
            ("output.output_per_day", po::value<int>()->default_value( 4 ),
               "Positive integer specifies number of outputs per day, Zero cancels output, Negative integer forces ouput at each timestep")
            ("output.save_forcing_fields", po::value<bool>()->default_value( false ), "")
            ("output.save_diagnostics", po::value<bool>()->default_value( false ), "")
            ("output.export_before_regrid", po::value<bool>()->default_value( false ),
                "if true, export results before regrid")
            ("output.export_after_regrid", po::value<bool>()->default_value( false ),
                "if true, export results after regrid")
#if 0
            //TODO issue193 uncomment these lines to set export variables using config file (finish another time)
            ("output.variables", po::value<std::vector<std::string>>()->multitoken()
                ->default_value( std::vector<std::string>{}, "")->composing(),
                "list of variable names (put on separate lines in config file)")
#endif

            // --exporter
            ("output.datetime_in_filename", po::value<bool>()->default_value( true ),
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
            ("ideal_simul.init_young_conc", po::value<double>()->default_value( 0. ), "")
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
            ("dynamics.C_lab", po::value<double>()->default_value( 2.0e6 ), "Pa")       // Cohesion value at the lab scale (10^6 Pa is the order of magnitude determined by Schulson).
            ("dynamics.nu0", po::value<double>()->default_value( 1./3. ), "")
            ("dynamics.tan_phi", po::value<double>()->default_value( 0.7 ), "")
            ("dynamics.compr_strength", po::value<double>()->default_value( 1e10 ), "Pa")
            ("dynamics.compaction_param", po::value<double>()->default_value( -20. ), "")

            // - C,h limits for where to use MEB rheology and where to use the Laplacian free drift thing
            ("dynamics.min_h", po::value<double>()->default_value( 0.05 ), "")
            ("dynamics.min_c", po::value<double>()->default_value( 0.01 ), "")

            ("dynamics.use_temperature_dependent_healing", po::value<bool>()->default_value( false ), "")
            ("dynamics.time_relaxation_damage", po::value<double>()->default_value( 25. ), "days")
            ("dynamics.deltaT_relaxation_damage", po::value<double>()->default_value( 20. ), "Kelvin")
            ("dynamics.undamaged_time_relaxation_sigma", po::value<double>()->default_value( 1e7 ), "seconds")
            // from V. Dansereau et al.: A Maxwell elasto-brittle rheology for sea ice modelling
            ("dynamics.exponent_relaxation_sigma", po::value<double>()->default_value( 5. ), "")
                // from V. Dansereau et al.: A Maxwell elasto-brittle rheology for sea ice modelling

            // - Water and air drag parameterizations
            ("dynamics.ERAi_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
            ("dynamics.ERA5_quad_drag_coef_air", po::value<double>()->default_value( 0.0020 ), "")
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

            // - Pressure term parameters
            ("dynamics.exponent_compression_factor", po::value<double>()->default_value( 1.5 ), "Power of ice thickness in the pressure term")
            ("dynamics.compression_factor", po::value<double>()->default_value( 10e3 ), "Max pressure for damaged converging ice")

            // - EVP!
            ("dynamics.substeps", po::value<int>()->default_value( 120 ),
                "Number of explicit sub-steps (default 120)")
            ("dynamics.evp.e", po::value<double>()->default_value( 2. ),
                "Ellipse ratio (default 2)")
            ("dynamics.evp.Pstar", po::value<double>()->default_value( 27.5e3 ),
                "P* (default 27.5e3)")
            ("dynamics.evp.C", po::value<double>()->default_value( 20 ),
                "Compaction parameter (C, default 20)")
            ("dynamics.evp.dmin", po::value<double>()->default_value( 1e-9 ),
                "Minimum delta (default 1e-9)")

            // - mEVP!
            ("dynamics.mevp.alpha", po::value<double>()->default_value( 500 ), "Alpha of the mEVP method (default 500)")
            ("dynamics.mevp.beta",  po::value<double>()->default_value( 500 ), "Beta of the mEVP method (default 500)")

             //-----------------------------------------------------------------------------------
             //! - Thermodynamics
             //-----------------------------------------------------------------------------------

            ("thermo.use_thermo_forcing", po::value<bool>()->default_value( true ), "")
            ("thermo.Qio-type", po::value<std::string>()->default_value( "basic" ), "")
            ("thermo.freezingpoint-type", po::value<std::string>()->default_value( "linear" ), "How to calculate the freezing point of sea water, either linear or unesco formula")
            ("thermo.albedoW", po::value<double>()->default_value( 0.07 ), "")
            ("thermo.alb_scheme", po::value<int>()->default_value( 3 ), "")
            ("thermo.flooding", po::value<bool>()->default_value( true ), "")
            ("thermo.alb_ice", po::value<double>()->default_value( 0.63 ), "")
            ("thermo.alb_sn", po::value<double>()->default_value( 0.88 ), "")
            ("thermo.I_0", po::value<double>()->default_value( 0.17 ), "")
            ("thermo.Qdw", po::value<double>()->default_value( 0.5 ), "")
            ("thermo.Fdw", po::value<double>()->default_value( 0. ), "")
            ("thermo.newice_type", po::value<int>()->default_value( 4 ), "4: YOUNG_ICE; else CLASSIC")
            ("thermo.melt_type", po::value<int>()->default_value( 1 ), "")
            ("thermo.hnull", po::value<double>()->default_value( 0.25 ), "")
            ("thermo.PhiF", po::value<double>()->default_value( 4. ), "")
            ("thermo.PhiM", po::value<double>()->default_value( 0.5 ), "")
            ("thermo.h_young_max", po::value<double>()->default_value( 0.5 ), "")
            ("thermo.h_young_min", po::value<double>()->default_value( 0.05 ), "")
            ("thermo.snow_cond", po::value<double>()->default_value( 0.3096 ),
                "snow conductivity (W/(K m)")

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

            // -- assimilation compensating flux
            ("thermo.use_assim_flux", po::value<bool>()->default_value(false),
             "Add a heat flux that compensates for assimilation of concentration")
            ("thermo.assim_flux_exponent", po::value<double>()->default_value(1.0),
             "Exponent of factor for heat flux that compensates for assimilation of concentration")

#ifdef AEROBULK
            ("thermo.ocean_bulk_formula", po::value<std::string>()->default_value( "coare" ), "Bulk formula to calculate ocean-atmosphere fluxes [ nextsim | coare (default) | coare3.5 | ncar | ecmwf ]")
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
            ("forecast.ec2_time_res_hours", po::value<double>()->default_value( 6. ),
                "specify the time resolution in hours here if want to change from 6")


             //-----------------------------------------------------------------------------------
             //! - Modules
             //-----------------------------------------------------------------------------------

#if defined(OASIS)
            ("coupler.component_name", po::value<std::string>()->default_value( "nxtsim" ), "Component name (6 characters) same as in the namcouple")
            ("coupler.timestep", po::value<int>()->default_value( 3600 ), "Coupling time step")
            ("coupler.exchange_grid_file", po::value<std::string>()->default_value( "coupler/NEMO.nc" ),
                 "File containing neccesary grid information for coupling. Name is relative to $NEXTSIM_DATA_DIR")
            // ("coupler.with_ocean", po::value<bool>()->default_value( false ), "Do we couple with an ocean model?")
            ("coupler.with_waves", po::value<bool>()->default_value( false ), "Do we couple with a wave model?")
            // ("coupler.with_atm", po::value<bool>()->default_value( false ), "Do we couple with an atmospheric model?")
            ("coupler.BGC_active", po::value<bool>()->default_value( false ), "Activate sending of MSLP and wind speed for the BGC model.")
            ("coupler.rcv_first_layer_depth", po::value<bool>()->default_value( false ), "Activate/deactivate (default) receiving of the depth of the first ocean layer (I_MLD in namcouple).")


            //-----------------------------------------------------------------------------------
            //!wave_coupling
            //-----------------------------------------------------------------------------------
            ("wave_coupling.receive_wave_stress", po::value<bool>()->default_value( true ),
             "Do we receive the wave stress from the wave model?")
            // FSD related
            ("wave_coupling.num_fsd_bins", po::value<int>()->default_value( 0 ), "Select a number of bins for FSD")
            ("wave_coupling.fsd_type", po::value<std::string>()->default_value("constant_size"), "Type of FSD bin width : constant_size or constant_area")
            ("wave_coupling.fsd_bin_cst_width", po::value<double>()->default_value( 10 ), "Select a number of constant width for FSD bins")
            ("wave_coupling.fsd_min_floe_size", po::value<double>()->default_value( 10 ), "Select a number of constant width for FSD bins")
            // Rheology related (careful to agree with WW3 !!)
            ("wave_coupling.floes_flex_strength", po::value<double>()->default_value( 0.27e6 ), "Value used in Williams et al. 2013")
            ("wave_coupling.floes_flex_young", po::value<double>()->default_value( 5.49e9 ), "Value used in Williams et al. 2013")
            ("wave_coupling.welding_type", po::value<std::string>()->default_value( "none" ), "none -> no welding, roach -> welding following Roach et al. 2018")
            ("wave_coupling.welding_kappa", po::value<double>()->default_value( 0.01 ), "Coagulation rate : values in range 0.001->0.01 (Roach et al., 2018)")
            ("wave_coupling.fsd_welding_use_scaled_area", po::value<bool>()->default_value( false ), "Roach et al., 2018-> True")
            ("wave_coupling.dmax_c_threshold", po::value<double>()->default_value( 0.1 ), "Sea ice concentration thereshold to determine dmax value.             By default, dmax is the size associated with the 9th decile of sea ice")
            ("wave_coupling.fsd_unbroken_floe_size", po::value<double>()->default_value( 1000. ), "Floe size assumed for unbroken sea ice")
            ("wave_coupling.fsd_damage_type", po::value<int>()->default_value(0), "Type of relationship between damage and FSD: 0 none, 1 related to unbroken floes concentration, 2 : related to break-up prob.")
            ("wave_coupling.fsd_damage_max", po::value<double>()->default_value(0.99), "Max. value for damage due to break-up")
            // Break-up parameters 
            ("wave_coupling.breakup_thick_min", po::value<double>()->default_value( 0.), "Minimum sea ice thickness assumed in breakup determination")
            ("wave_coupling.breakup_prob_type", po::value<int>()->default_value(0), "0: Default, use tau_w as a time scale for breakup ")
            ("wave_coupling.breakup_cell_average_thickness", po::value<bool>()->default_value(false), "Are we using M_thick (volume, if false) or sea ice thickness averaged on sea ice cover only (if true) in break-up determination ?")
            ("wave_coupling.breakup_timescale_tuning", po::value<double>()->default_value( 1.), "Time scale tuning for breakup_prob_type=0 (tau_w)")
            ("wave_coupling.breakup_type", po::value<std::string>()->default_value( "uniform_size" ), "Type of break-up redistribution : none, uniform in size,(uniform_size) or followin Zhang et al. 2015 (zhang)")
            ("wave_coupling.breakup_coef1", po::value<double>()->default_value( 0.5), "Coef1 is used to tune the size D of fsd regime transition : D=lambda_wave*coef1")
            ("wave_coupling.breakup_coef2", po::value<double>()->default_value( 1.), "Coef2 is used to tune the size range on which the fsd transition is seen : range=lambda_wave*coef2")
            ("wave_coupling.breakup_coef3", po::value<double>()->default_value( 1.), "Coef3 is used to tune the size range on which the fsd transition is seen : range=dflex*coef3")
            ("wave_coupling.breakup_prob_cutoff", po::value<double>()->default_value( 0.0015), "Minimum value of break-up probability for which break-up occurs (if prob<prob_cut_off -> prob=0.)")

            ("wave_coupling.distinguish_mech_fsd", po::value<bool>()->default_value( true ), "Do we distinguish a mech. FSD from the real FSD ?")
            ("wave_coupling.debug_fsd", po::value<bool>()->default_value( false ), "Do we check ice area conservation in FSD each time it is modified ?")
#endif

            // for ensemble forcing
            ("statevector.ensemble_member", po::value<int>()->default_value(0),
                "id of ensemble member (NB starts from 1)")

            //-----------------------------------------------------------------------------------
            //!Age settings
            //-----------------------------------------------------------------------------------
            ("age.reset_date", po::value<std::string>()->default_value( "0915" ), "Select the date which resets all ice to multiyear ice")
            ("age.reset_by_date", po::value<bool>()->default_value( false ), "Choose whether to reset myi on a date (true) or by an amount of melt days (false)")
            ("age.include_young_ice", po::value<bool>()->default_value( true ), "If ice-type is young ice, choose whether to include it when resetting multiyear on reset_date")
            ("age.reset_freeze_days", po::value<double>()->default_value( 3. ), "If reset by freeze days, this is number of consecutive freezing days before ice becomes myi")
            ("age.equal_ridging", po::value<bool>()->default_value( false ), "When ridging, if fyi is present, ridge fyi preferentially (false) or both myi and fyi evenly (true)")
            ("age.equal_melting", po::value<bool>()->default_value( true ), "When ridging, if fyi is present, melt fyi preferentially (false) or both myi and fyi evenly (true)")

#if defined(WAVES)
        ;
        return desc.add( Wim::descrWimOptions() );
#else
        ;
        return desc;
#endif
    }

} // Nextsim
