/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   graphcsrmpi.hpp
 * @author Abdoulaye Samake <abdama@beijing.ad.nersc.no>
 * @date   Mon Feb 22 11:49:50 2016
 */

#ifndef __GraphCSRMPI_HPP
#define __GraphCSRMPI_HPP 1

#include <vector>

namespace Nextsim
{
class GraphCSRMPI
{
public:

    GraphCSRMPI();

	GraphCSRMPI(std::vector<int> const& d_nnz,
                std::vector<int> const& o_nnz,
                std::vector<int> const& global_without_ghost,
                std::vector<int> const& global_with_ghost);

	std::vector<int> const& nNzOnProc() const {return M_dnnz;}
    std::vector<int> const& nNzOffProc() const {return M_onnz;}
    std::vector<int> const& globalIndicesWithoutGhost() const {return M_global_without_ghost;}
    std::vector<int> const& globalIndicesWithGhost() const {return M_global_with_ghost;}

private:

	std::vector<int> M_dnnz;
    std::vector<int> M_onnz;
    std::vector<int> M_global_without_ghost;
    std::vector<int> M_global_with_ghost;
};

} // Nextsim
#endif
