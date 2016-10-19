/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   graphcsr.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Fri Jan  8 15:02:02 2016
 */

#ifndef __GraphCSR_HPP
#define __GraphCSR_HPP 1

#include <vector>

namespace Nextsim
{
class GraphCSR
{
public:

    GraphCSR();

	GraphCSR(std::vector<int> const& nnz,
	         std::vector<int> const& ia,
	         std::vector<int> const& ja,
	         std::vector<double> const& a);

	std::vector<int> const& nNz() const {return M_nnz;}
	std::vector<int> const& ia() const {return M_ia;}
	std::vector<int> const& ja() const {return M_ja;}
	std::vector<double> const& a() const {return M_a;}

private:

	std::vector<int> M_nnz;
	std::vector<int> M_ia;
	std::vector<int> M_ja;
	std::vector<double> M_a;

};

} // Nextsim
#endif
