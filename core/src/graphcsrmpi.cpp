/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   graphcsrmpi.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Feb 22 11:48:16 2016
 */

#include <graphcsrmpi.hpp>

namespace Nextsim
{
GraphCSRMPI::GraphCSRMPI()
	:
	M_dnnz(),
    M_onnz(),
    M_global()
{}

GraphCSRMPI::GraphCSRMPI(std::vector<int> const& d_nnz,
                         std::vector<int> const& o_nnz,
                         std::vector<int> const& global)
	:
	M_dnnz(d_nnz),
    M_onnz(o_nnz),
    M_global(global)
{}


} // Nextsim
