/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   environment.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Thu Jul  2 17:50:41 2015
 */


#ifndef FEELPP_ENVIRONMENT_HPP
#define FEELPP_ENVIRONMENT_HPP 1

#include <communicator.hpp>
#include <petscsys.h>
#include <petsc.h>

namespace Nextsim
{

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

    ~Environment();

    static Communicator const& comm() { return mpicomm; }

    static MemoryUsage logMemoryUsage(std::string const& message);

private:

    boost::mpi::environment mpienv;
    static Communicator mpicomm;

};

} // Nextsim
#endif // FEELPP_ENVIRONMENT_HPP
