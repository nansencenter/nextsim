/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   environment.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Thu Jul  2 17:50:41 2015
 */


#ifndef Environment_HPP
#define Environment_HPP 1

#include <communicator.hpp>
#include <petscsys.h>
#include <petsc.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

namespace Nextsim
{
namespace fs = boost::filesystem;
namespace po = boost::program_options;

struct MemoryUsage
{
    MemoryUsage()
        :
        memory_usage(0),
        petsc_malloc_usage(0),
        petsc_malloc_maximum_usage(0)
    {}

    MemoryUsage(MemoryUsage const& m )
        :
        memory_usage(m.memory_usage),
        petsc_malloc_usage(m.petsc_malloc_usage),
        petsc_malloc_maximum_usage(m.petsc_malloc_maximum_usage)
    {}

    MemoryUsage& operator=(MemoryUsage const& m )
    {
        if ( this != &m )
        {
            memory_usage = m.memory_usage;
            petsc_malloc_usage = m.petsc_malloc_usage;
            petsc_malloc_maximum_usage = m.petsc_malloc_maximum_usage;
        }
        return *this;
    }

    PetscLogDouble memory_usage;
    PetscLogDouble petsc_malloc_usage;
    PetscLogDouble petsc_malloc_maximum_usage;

};

class Environment
{
public:

    Environment( int& argc, char** &argv );

    Environment( int& argc, char** &argv, po::options_description desc, Communicator const& comm = Communicator::commSelf());

    ~Environment();

    static Communicator const& comm() { return mpicomm; }

    static po::variables_map const& vm() { return vmenv; }

    //! return $NEXTSIM_MESH_DIR
    static fs::path const& nextsimMeshDir() { return nextsim_mesh_dir_env; }

    //! return $NEXTSIM_DATA_DIR
	static fs::path const& nextsimDataDir() { return nextsim_data_dir_env; }

    //! get $NEXTSIM_DATA_DIR and $NEXTSIM_MESH_DIR variables,
    //! and some sub-directories
    void setEnvironmentVariables();

    static MemoryUsage logMemoryUsage(std::string const& message);

private:

    boost::mpi::environment mpienv;
    static Communicator mpicomm;
    static po::variables_map vmenv;

    //! $NEXTSIM_DATA_DIR
    static fs::path nextsim_data_dir_env;

    //! $NEXTSIM_MESH_DIR
    static fs::path nextsim_mesh_dir_env;
};

} // Nextsim
#endif // ENVIRONMENT_HPP
