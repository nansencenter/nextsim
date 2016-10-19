/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

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
    M_global_without_ghost(),
    M_global_with_ghost()
{}

GraphCSRMPI::GraphCSRMPI(std::vector<int> const& d_nnz,
                         std::vector<int> const& o_nnz,
                         std::vector<int> const& global_without_ghost,
                         std::vector<int> const& global_with_ghost)
	:
	M_dnnz(d_nnz),
    M_onnz(o_nnz),
    M_global_without_ghost(global_without_ghost),
    M_global_with_ghost(global_with_ghost)
{}

GraphCSRMPI::GraphCSRMPI(GraphCSRMPI const& g)
    :
    M_dnnz(g.M_dnnz),
    M_onnz(g.M_onnz),
    M_global_without_ghost(g.M_global_without_ghost),
    M_global_with_ghost(g.M_global_with_ghost)
{}

} // Nextsim
