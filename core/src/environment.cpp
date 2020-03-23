/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   environment.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Jul  6 17:31:33 2015
 */

#include <environment.hpp>

namespace Nextsim
{

Environment::Environment( int& argc, char** &argv )
    :
    mpienv(argc, argv)
{
    mpicomm = Communicator::commSelf();

    this->setEnvironmentVariables();

    int ierr = 0;
    ierr = PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
    CHKERRABORT( mpicomm, ierr );
}

Environment::Environment( int& argc, char** &argv, po::options_description desc )
    :
    mpienv(argc, argv)
{
    //! - 1) Set NEXTSIM_[MESH,DATA]_DIR for inputs
    this->setEnvironmentVariables();

    //! - 2) Read the config files to get the run-time options
    try
    {
        //po::store(po::parse_command_line(argc, argv, desc),vmenv);

        po::store(po::command_line_parser(argc, argv)
                  .options(desc)
                  .style(po::command_line_style::unix_style)// | po::command_line_style::allow_long_disguise)
                  .allow_unregistered()
                  .run(),
                  vmenv);


        if ( vmenv.count( "config-file" ) || vmenv.count( "config-files" ) )
        {
            if ( vmenv.count( "config-file" ) )
                configFiles = {vmenv["config-file"].as<std::string>()};
            else if ( vmenv.count( "config-files" ) )
                configFiles = vmenv["config-files"].as<std::vector<std::string>>();
        }

        if (configFiles.size() == 0)
        {
            std::stringstream msg;
            msg << "Please provide (a) config file(s) with command line options\n";
            msg << "--config-file=CFGFILE or --config-files=CFGFILE1 CFGFILE2\n";
            throw std::runtime_error(msg.str());
        }

        // loop over config files
        for ( std::string cfgfile : configFiles )
        {
            if ( fs::exists( cfgfile ) )
            {
                if (Communicator::commSelf().rank()==0)
                    std::cout << "Reading " << cfgfile << "...\n";
                std::ifstream ifs( cfgfile.c_str() );
                // 3rd argument of parse_config_file: true for ignoring unknown options (false else)
                po::store( parse_config_file( ifs, desc, false ), vmenv );
            }
            else
                std::cout << "Cannot find " << "config-file `" << cfgfile <<"`\n";
        }

        po::notify(vmenv);
    }

    // catches program_options exceptions

    catch (po::error& e)
    {
        throw std::runtime_error(std::string(e.what()));
    }

    catch ( std::exception& e )
    {
        std::cout << "Application option parsing: unknown option:" << e.what()
                  << " (the .cfg file or some options may not have been read properly)\n";

        throw std::runtime_error(std::string(e.what()));
    }

    catch ( ... )
    {
        std::cout << "Application option parsing: unknown exception triggered"
                  << "(the .cfg file or some options may not have been read properly)\n";

        throw std::runtime_error("...");
    }

    //! -3) Initialise communicator, PETSc, and OASIS (if compiled in)
#ifdef OASIS
    // For OASIS we need to get the local communicator first

    MPI_Comm localComm, cplComm; // local and couple MPI communicators
    //const std::string comp_name = std::string("nxtsim");  // Component name (6 characters) same as in the namcouple
    const std::string comp_name = vmenv["coupler.component_name"].as<std::string>();

    // Initialise OASIS
    int comp_id;
    int ierror = OASIS3::init_comp(&comp_id, comp_name);
    if (ierror != 0) {
        std::cout << "oasis_init_comp abort by nextsim with error code " << ierror << std::endl;
        OASIS3::abort(comp_id, comp_name, "Problem calling OASIS3::init_comp");
    }

    // Get the local communicator
    ierror = OASIS3::get_localcomm(&localComm);
    if (ierror != 0) {
        std::cout << "oasis_get_localcomm abort by nextsim with error code " << ierror << std::endl;
        OASIS3::abort(comp_id, comp_name, "Problem calling OASIS3::get_localcomm");
    }

    // create a Nextsim::Communicator from the OASIS communicator
    mpicomm = Communicator(localComm);

    // Tell PETSc to use the right communicator
    PETSC_COMM_WORLD = localComm;

    // Create the coupler communicator - only root communicates with the coupler
    ierror = OASIS3::create_couplcomm(mpicomm.rank()==0, &localComm, &cplComm);
    if (ierror != 0) {
        std::cout << "oasis_create_couplcomm abort by nextsim with error code " << ierror << std::endl;
        OASIS3::abort(comp_id, comp_name, "Problem calling OASIS3::create_couplcomm");
    }
#else
    mpicomm = Communicator::commSelf();
#endif

    int ierr = 0;
    ierr = PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
    CHKERRABORT( mpicomm, ierr );

    //! -4) set other useful variables it would be convenient to have access to
    //! across multiple classes
    //! * nextsim .mppfile
    nextsim_mppfile = (boost::format( "%1%/%2%" )
            % this->nextsimMeshDir().string()
            % this->vm()["mesh.mppfile"].as<std::string>()
            ).str();

    const boost::unordered_map<const std::string, LogLevel> str2log = boost::assign::map_list_of
        ("debug", DEBUG)
        ("verbose", VERBOSE)
        ("info", INFO)
        ("warning", WARNING)
        ("error", ERROR);
    std::string const s = vmenv["debugging.log-level"].as<std::string>();
    if(str2log.count(s)==0)
        throw std::runtime_error("bad value for option debugging.log-level " + s);
    log_level = str2log.find(s)->second;
    log_all = vmenv["debugging.log-all"].as<bool>();
}


Environment::~Environment()
{
    int ierr = 0;
    ierr = PetscFinalize();
}


void
Environment::setEnvironmentVariables()
{
    //! get and check $NEXTSIM_MESH_DIR and $NEXTSIM_DATA_DIR
    //! and set them to nextsim_mesh_dir_env and nextsim_data_dir_env

    char * senv = ::getenv( "NEXTSIM_MESH_DIR" );
    if ( (senv == NULL) || (senv[0] == '\0') )
    {
        std::cout << "you must define 'NEXTSIM_MESH_DIR' environment variable"
            << " for directory with meshes and mppfiles"<<"\n";
        throw std::logic_error("invalid environment variable");
    }
    nextsim_mesh_dir_env = fs::path(std::string(senv));

	senv = ::getenv( "NEXTSIM_DATA_DIR" );
    if ( (senv == NULL) || (senv[0] == '\0') )
    {
        std::cout << "you must define 'NEXTSIM_DATA_DIR' environment variable"
            << " for directory with netcdf input files and other input files"
            <<" (eg drifter text files)\n";
        throw std::logic_error("invalid environment variable");
    }
    nextsim_data_dir_env = fs::path(std::string(senv));
}//setEnvironmentVariables


Communicator Environment::mpicomm;
po::variables_map Environment::vmenv;
fs::path Environment::nextsim_data_dir_env;
fs::path Environment::nextsim_mesh_dir_env;
std::string Environment::nextsim_mppfile;
std::vector<std::string> Environment::configFiles;
LogLevel Environment::log_level;
bool Environment::log_all;

MemoryUsage
Environment::logMemoryUsage(std::string const& message)
{
    MemoryUsage mem;

    PetscMemoryGetCurrentUsage( &mem.memory_usage );
    std::cout << message
              << " PETSC get current memory usage (resident memory): ["
              << mem.memory_usage/1e3
              << " KB]  ["
              << mem.memory_usage/1e6
              << " MB]  ["
              << mem.memory_usage/1e9
              << " GB]\n" ;
    //PetscMemoryGetMaximumUsage( &mem );
    //LOG(INFO) << logMessage << " PETSC get maximum memory usag (resident memory): " << mem/1e6 << "  MB " << mem/1e9 << " GB" ;

    PetscMallocGetCurrentUsage( &mem.petsc_malloc_usage );
    std::cout << message
              << " PETSC get current PETSC Malloc usage: ["
              << mem.petsc_malloc_usage/1e3
              << " KB]  ["
              << mem.petsc_malloc_usage/1e6
              << " MB]  ["
              << mem.petsc_malloc_usage/1e9
              << " GB]\n" ;

    PetscMallocGetMaximumUsage( &mem.petsc_malloc_maximum_usage );
    std::cout << message
              << " PETSC get maximum PETSC Malloc usage(largest memory ever used so far): ["
              << mem.petsc_malloc_maximum_usage/1e3
              << " KB]  ["
              << mem.petsc_malloc_maximum_usage/1e6
              << " MB]  ["
              << mem.petsc_malloc_maximum_usage/1e9
              << " GB]\n" ;

    return mem;
}//logMemoryUsage

} // Nextsim
