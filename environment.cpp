/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

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
    //mpicomm()
{
    mpicomm = Communicator::commSelf();
    int ierr = 0;
    ierr = PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
    CHKERRABORT( mpicomm, ierr );
}

Environment::~Environment()
{
    int ierr = 0;
    ierr = PetscFinalize();
}

Communicator Environment::mpicomm;


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
}


} // Nextsim
