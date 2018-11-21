/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

#ifndef __aerobulk_HPP
#define __aerobulk_HPP 1

#include <vector>
#include <string>
#include <cassert>
#include <cstdarg>
#include <stdexcept>

namespace aerobulk
{

    enum class algorithm
    {
        neXtSIM = 0,
        COARE   = 1,
        COARE35 = 2,
        NCAR    = 3,
        ECMWF   = 4
    };

    std::string algorithm_to_string(algorithm algo);

    // To check the size of the inputs
    int check_sizes(int count, ...);

    // Interface to calculate fluxes and drag
    /* rad_sw and rad_lw default to an empty vector, in which case we don't use
     * the cool-skin warm-layer scheme */
    void model( algorithm algo, const double z_t, const double z_u,
            const std::vector<double>& sst,  const std::vector<double>& t_zt, const std::vector<double>& q_zt,
            const std::vector<double>& U_zu, const std::vector<double>& slp,
            std::vector<double>& QL, std::vector<double>& QH, std::vector<double>& Cd_rho_U,
            const std::vector<double>& rad_sw = std::vector<double>(),
            const std::vector<double>& rad_lw = std::vector<double>());

}

#endif
