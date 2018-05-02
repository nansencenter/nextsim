/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

#include <options_wim.hpp>

namespace Wim
{
    po::options_description
    descrWimOptions()
    {
        po::options_description desc("Options");

        desc.add_options()
            ("help", "Print help messages")
            ("config-file", po::value<std::string>(),
                  "specify .cfg file")

            // [wimgrid]
            // - to read in a grid file
            ("wimgrid.gridfilename", po::value<std::string>()->default_value( "" ),
                  "wim grid binary filename - default is to create the grid from parameters in cfg file (uncoupled code), or to use the nextsim mesh to determine the parameters (coupled code)")

            // - to create a new regular grid
            ("wimgrid.nx", po::value<int>()->default_value( 150 ),
                  "Record length in x direction")
            ("wimgrid.ny", po::value<int>()->default_value( 10 ),
                  "Record length in y direction")
            ("wimgrid.dx", po::value<double>()->default_value( 4e+3 ),
                  "Resolution in x direction [m]")
            ("wimgrid.dy", po::value<double>()->default_value( 4e+3 ),
                  "Resolution in y direction [m]")
            ("wimgrid.xmin", po::value<double>()->default_value( -298.e+3 ),
                  "xmin [m]")
            ("wimgrid.ymin", po::value<double>()->default_value( -60e+3 ),
                  "ymin [m]")

            // - extra grid options
            ("wimgrid.landon3edges", po::value<bool>()->default_value( false ),
                  "Add land on upper,lower and RH edges")
            ("wimgrid.useregulargridtools", po::value<bool>()->default_value( false ),
                  "use regular grid tools if possible (else triangulate grid and use mesh tools)")

            // [wimsetup]
            ("wimsetup.initialtime", po::value<std::string>()->default_value( "2015-01-01 00:00:00" ),
                  "Initial time")
            ("wimsetup.duration", po::value<double>()->default_value( 43200.0 ),
                  "length of simulation [s]")

            // - wim spectral discretisation
            ("wimsetup.tmin", po::value<double>()->default_value( 2.5 ),
                  "Minimum wave period in a spectrum [s]")
            ("wimsetup.tmax", po::value<double>()->default_value( 25. ),
                  "Maximum wave period in a spectrum [s]")
            ("wimsetup.nwavefreq", po::value<int>()->default_value( 1 ),
                  "Number of wave frequencies")
            ("wimsetup.nwavedirn", po::value<int>()->default_value( 16 ),
                  "Number of wave directions")

            // - forcing (external wave model)
            ("wimsetup.wave-type", po::value<std::string>()->default_value( "set_in_wim" ),
                "set_in_wim, ww3a, eraiw_1deg")
            ("wimsetup.wave-time-interp-option", po::value<std::string>()->default_value( "step" ), "step, linear")

            // [wim]
            // - wave atten options
            ("wim.atten", po::value<bool>()->default_value( true ),
                  "Do attenuation")
            ("wim.scatmod", po::value<std::string>()->default_value( "dissipated" ),
                  "Scattered energy is dissipated (=dissipated), distributed isotropically (=isotropic)")
            ("wim.young", po::value<double>()->default_value( 5.49e+9 ),
                  "Young's modulus [Pa]")
            ("wim.dragrp", po::value<double>()->default_value( 13. ),
                  "Robinson-Palmer drag coefficient [Pa.s/m]")

            // - numerics
            ("wim.advopt", po::value<std::string>()->default_value( "y-periodic" ),
                  "Not periodic (=notperiodic), periodic in y only (=y-periodic), periodic in both x,y (=xy-periodic)")
            ("wim.advdim", po::value<int>()->default_value( 2 ),
                  "Dimension of advection scheme (1 or 2)")
            ("wim.steady", po::value<bool>()->default_value( true ),
                  "Steady-state (=true), or not steady-state (=false)")
            ("wim.cfl", po::value<double>()->default_value( 0.7 ),
                  "CFL number")

            // - breaking and FSD
            ("wim.breaking", po::value<bool>()->default_value( true ),
                  "Do breaking (=true), or turn off breaking (=false)")
            ("wim.fsdopt", po::value<std::string>()->default_value( "PowerLawSmooth" ),
                  "FSD parameterisation: 'PowerLawSmooth' or 'RG'")

            // -- numerical parameters
            ("wim.dfloemin", po::value<double>()->default_value( 20. ),
                  "Minimum floe size [m]")
            ("wim.cicemin", po::value<double>()->default_value( 0.05 ),
                  "Minimum ice conc considered by WIM")
            ("wim.dfloepackthresh", po::value<double>()->default_value( 400. ),
                  "Don't let Dmax grow above this value [m]")

            // - other bool param's
            ("wim.refhsice", po::value<bool>()->default_value( false ),
                  "Inside ice, Hs corresponds to water (=false) or ice (=true) displacement")
            ("wim.useicevel", po::value<bool>()->default_value( false ),
                  "Inside ice, use <<correct>> group velocity (=true), or water group velocity (=false)")

            //initial conditions in idealised simulation
            ("wim.hsinc", po::value<double>()->default_value( 3. ),
                  "Incident significant wave height [m]")
            ("wim.tpinc", po::value<double>()->default_value( 12. ),
                  "Incident peak period [s]")
            ("wim.mwdinc", po::value<double>()->default_value( -90. ),
                  "Incident mean wave-from direction [deg]")
            ("wim.unifc", po::value<double>()->default_value( 0.7 ),
                  "Initial const conc")
            ("wim.unifh", po::value<double>()->default_value( 1. ),
                  "Initial const thickness [m]")
            ("wim.dfloepackinit", po::value<double>()->default_value( 300. ),
                  "Initial value in pack (unbroken) ice [m]")

            // [wimdiag]
            // outputs/diagnostics of WIM
            ("wimdiag.checkinit", po::value<bool>()->default_value( true ),
                  "Do/don't dump initial states of each call to WIM.run() to binary files  (true/false)")
            ("wimdiag.checkprog", po::value<bool>()->default_value( true ),
                  "Do/don't dump intermediate states to binary files (true/false)")
            ("wimdiag.checkfinal", po::value<bool>()->default_value( true ),
                  "Do/don't dump final states after each call to WIM.run() to binary files (true/false)")
            ("wimdiag.checkincwaves", po::value<bool>()->default_value( true ),
                  "Do/don't dump input wave fields after each call to WIM.run() to binary files (true/false)")
            ("wimdiag.savelog", po::value<bool>()->default_value( true ),
                  "Do/don't save diagnostic file after each call to WIM.run() to text file (true/false)")
            ("wimdiag.dumpfreq", po::value<int>()->default_value( 10 ),
                  "frequency of dumping (# WIM timesteps)")
            ("wimdiag.outparentdir", po::value<std::string>()->default_value( "out_cpp" ),
                  "Parent directory for the output files")
            ("wimdiag.itest", po::value<int>()->default_value( -1 ),
                  "row index for detailed diagnostics (<0 for no diagnostics)")
            ("wimdiag.jtest", po::value<int>()->default_value( -1 ),
                  "column index for detailed diagnostics (<0 for no diagnostics)")


            // [nextwim]
            // coupling of wim to nextsim

            ("nextwim.use_wim", po::value<bool>()->default_value( false ),
                "use the WIM")

            ("nextwim.couplingfreq", po::value<int>()->default_value( 20 ),
                  "Coupling frequency between neXtSIM and WIM (# neXtSIM time-steps)")
            ("nextwim.coupling-option", po::value<std::string>()->default_value( "break_on_mesh" ),
                  "Coupling option: naive->interp nfloes onto mesh after exiting WIM; break_on_mesh->import mesh and do breaking on mesh in parallel; run_on_mesh->run WIM on neXtSIM mesh")
            ("nextwim.wim_damage_mesh", po::value<bool>()->default_value( true ),
                  "If ice is broken by waves, increase damage to nextwim.wim_damage_value (if nextwim_coupling_option=breaking_on_mesh)")
            ("nextwim.wim_damage_value", po::value<double>()->default_value( 0.999 ),
                  "If ice is broken by waves, increase damage to this value (if wim_damage_mesh=true)")

            ("nextwim.applywavestress", po::value<bool>()->default_value( true ),
                  "Use wave stress from WIM in neXtSIM momentum equations")

            // - exporting
            ("nextwim.exportresults", po::value<bool>()->default_value( true ),
                  "Export results in coupled mode")
            ("nextwim.export_diags_mesh", po::value<bool>()->default_value( false ),
                  "export wave diagnostics on mesh or not")
            ("nextwim.export_after_wim_call", po::value<bool>()->default_value( false ), "")
            ("nextwim.test_and_exit", po::value<bool>()->default_value( false ),
                  "Stop the run after 1 call to the wim")

            ;
        return desc;
    }

} // WIMOPT
