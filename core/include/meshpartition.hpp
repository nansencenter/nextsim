/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   meshpartition.hpp
 * @author Abdoulaye Samake <abdama@beijing.ad.nersc.no>
 * @date   Mon Aug 15 16:44:48 2016
 */

#ifndef __MeshPartition_HPP
#define __MeshPartition_HPP 1

//#include <meshPartitionOptions.h>
//#include <meshPartition.h>
#include <Context.h>
//#include <GModel.h>

namespace Nextsim
{
namespace mesh
{
    enum class Partitioner
    {
        CHACO = 1,
        METIS = 2
    };

    enum class PartitionSpace
    {
        MEMORY = 0,
        DISK   = 1
    };
} // mesh
} // Nextsim
#endif
