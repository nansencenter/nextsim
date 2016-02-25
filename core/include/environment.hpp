/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

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

    Environment( int& argc, char** &argv, po::options_description desc );

    ~Environment();

    static Communicator const& comm() { return mpicomm; }

    static po::variables_map const& vm() { return vmenv; }

    static fs::path const& nextsimDir() { return nextsimdirenv; }
	
    static fs::path const& simdataDir() { return simdatadirenv; }
	
	static fs::path const& simforecastDir() { return simforecastdirenv; }

    static MemoryUsage logMemoryUsage(std::string const& message);

private:

    boost::mpi::environment mpienv;
    static Communicator mpicomm;
    static po::variables_map vmenv;
    static fs::path nextsimdirenv;
    static fs::path simdatadirenv;
    static fs::path simforecastdirenv;

};

} // Nextsim
#endif // ENVIRONMENT_HPP
