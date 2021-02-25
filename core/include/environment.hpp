/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   environment.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Thu Jul  2 17:50:41 2015
 */


#ifndef Environment_HPP
#define Environment_HPP 1

#include <communicator.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/unordered_map.hpp>
#include <boost/assign/list_of.hpp>
#include "debug.hpp"
#if defined OASIS
#include<oasis_cpp_interface.h>
#endif

namespace Nextsim
{
namespace fs = boost::filesystem;
namespace po = boost::program_options;

class Environment
{
public:

    Environment( int& argc, char** &argv );

    Environment( int& argc, char** &argv, po::options_description desc );

    ~Environment();

    static Communicator const& comm() { return mpicomm; }

    static po::variables_map const& vm() { return vmenv; }

    //! return $NEXTSIM_MESH_DIR
    static fs::path const& nextsimMeshDir() { return nextsim_mesh_dir_env; }

    //! return $NEXTSIM_DATA_DIR
	static fs::path const& nextsimDataDir() { return nextsim_data_dir_env; }

    //! return neXtSIM .mpp file
	static std::string const& nextsimMppfile() { return nextsim_mppfile; }

    //! return neXtSIM config files
	static std::vector<std::string> const& nextsimConfigFiles() { return configFiles; }
    
    //! return log level
	static LogLevel const& logLevel() { return log_level; }

    //! return bool to log on all processors
	static bool const& logAll() { return log_all; }

    //! get $NEXTSIM_DATA_DIR and $NEXTSIM_MESH_DIR variables,
    //! and some sub-directories
    void setEnvironmentVariables();

private:

    boost::mpi::environment mpienv;
    static Communicator mpicomm;
    static po::variables_map vmenv;

    //! $NEXTSIM_DATA_DIR
    static fs::path nextsim_data_dir_env;

    //! $NEXTSIM_MESH_DIR
    static fs::path nextsim_mesh_dir_env;

    //! .mpp file for nextsim
    static std::string nextsim_mppfile;


    //! neXtSIM config files
    static std::vector<std::string> configFiles;

    //! log level
    static LogLevel log_level;

    //! log on all cpus
    static bool log_all;
};

} // Nextsim
#endif // ENVIRONMENT_HPP
