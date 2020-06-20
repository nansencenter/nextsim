/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   petsc.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Jul  6 13:47:19 2015
 */

#ifndef __NextsimPETSc_H
#define __NextsimPETSc_H 1

#include <petsc.h>
#include <petscerror.h>

#if !defined(PETSC_VERSION_LT) || !defined(PETSC_VERSION_GT) || !defined(PETSC_VERSION_GE)
#define PETSC_VERSION_LESS_THAN(major,minor,subminor)                   \
	((PETSC_VERSION_MAJOR < (major) ||                                  \
	  (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR < (minor) || \
	                                      (PETSC_VERSION_MINOR == (minor) && \
	                                       PETSC_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)

#define PETSC_VERSION_GREATER_THAN(major,minor,subminor)                \
	((PETSC_VERSION_MAJOR > (major) ||                                  \
	  (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR > (minor) || \
	                                      (PETSC_VERSION_MINOR == (minor) && \
	                                       PETSC_VERSION_SUBMINOR > (subminor))))) ? 1 : 0)

#define PETSC_VERSION_GREATER_OR_EQUAL_THAN(major,minor,subminor)       \
	((PETSC_VERSION_MAJOR > (major) ||                                  \
	  (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR > (minor) || \
	                                      (PETSC_VERSION_MINOR == (minor) && \
	                                       PETSC_VERSION_SUBMINOR >= (subminor))))) ? 1 : 0)
#else
#define PETSC_VERSION_LESS_THAN(major,minor,subminor)                   \
	PETSC_VERSION_LT(major,minor,subminor)

#define PETSC_VERSION_GREATER_THAN(major,minor,subminor)                \
	PETSC_VERSION_GT(major,minor,subminor)

#define PETSC_VERSION_GREATER_OR_EQUAL_THAN(major,minor,subminor)       \
	PETSC_VERSION_GE(major,minor,subminor)
#endif


#endif // __NextsimPETSc_H
