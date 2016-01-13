/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   graphcsr.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Fri Jan  8 15:21:09 2016
 */

#include <graphcsr.hpp>

namespace Nextsim
{
GraphCSR::GraphCSR()
	:
	M_nnz(),
	M_ia(),
	M_ja(),
	M_a()
{}

GraphCSR::GraphCSR(std::vector<int> const& nnz,
                   std::vector<int> const& ia,
                   std::vector<int> const& ja,
                   std::vector<double> const& a)
	:
	M_nnz(nnz),
	M_ia(ia),
	M_ja(ja),
	M_a(a)
{}


} // Nextsim
