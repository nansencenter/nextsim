/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   nextsimgmsh.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Wed May  1 13:36:05 2019
 */

#ifndef __GMSHVERSION_HPP
#define __GMSHVERSION_HPP 1

//#include <Gmsh.h>
#include <GmshVersion.h>

#define GMSH_VERSION_LESS_THAN(major,minor,subminor)                   \
    ((GMSH_MAJOR_VERSION < (major) ||                                  \
      (GMSH_MAJOR_VERSION == (major) && (GMSH_MINOR_VERSION < (minor) || \
                                          (GMSH_MINOR_VERSION == (minor) && \
                                           GMSH_PATCH_VERSION < (subminor))))) ? 1 : 0)

#define GMSH_VERSION_GREATER_THAN(major,minor,subminor)                \
    ((GMSH_MAJOR_VERSION > (major) ||                                  \
      (GMSH_MAJOR_VERSION == (major) && (GMSH_MINOR_VERSION > (minor) || \
                                          (GMSH_MINOR_VERSION == (minor) && \
                                           GMSH_PATCH_VERSION > (subminor))))) ? 1 : 0)

#define GMSH_VERSION_GREATER_OR_EQUAL_THAN(major,minor,subminor)       \
    ((GMSH_MAJOR_VERSION > (major) ||                                  \
      ((GMSH_MAJOR_VERSION == (major)) && ((GMSH_MINOR_VERSION > (minor)) || \
                                         ((GMSH_MINOR_VERSION == (minor)) && \
                                          ( GMSH_PATCH_VERSION >= (subminor))) ))) ? 1 : 0)
#endif // __GMSHVERSION_HPP
