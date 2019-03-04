/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

#include "aerobulk.hpp"

extern "C"
{
    void aerobulk_nextsim_skin( const char *, const double *, const double *, const double *, const double *,
            const double *, const double *, const double *,
            double *, double *, double *, double *,
            int *, int *,
            const double *, const double *);

    void aerobulk_nextsim_no_skin( const char *, const double *, const double *, const double *, const double *,
            const double *, const double *, const double *,
            double *, double *, double *, double *,
            int *, int *);
}

std::string aerobulk::algorithm_to_string(aerobulk::algorithm algo)
{
    std::string return_value = std::string("unknown");
    switch (algo)
    {
        case aerobulk::algorithm::neXtSIM:
            throw std::logic_error("You cannot call aerobulk routines with algorithm == neXtSIM\n");
        case aerobulk::algorithm::COARE:
            return_value = std::string("coare");
        case aerobulk::algorithm::COARE35:
            return_value = std::string("coare35");
        case aerobulk::algorithm::NCAR:
            return_value = std::string("ncar");
        case aerobulk::algorithm::ECMWF:
            return_value = std::string("ecmwf");
    }
    return return_value;
}

// Chekc input sizes
int aerobulk::check_sizes(int count, ...)
{
    va_list ap;

    va_start(ap, count); // Requires the last fixed parameter (to get the address)
    int size = va_arg(ap, int); // Get the first size - the one we compare everyone else's with

    for (int i=1; i<count; i++) // Start from 1 because we already called var_arg once
        assert( size == va_arg(ap, int) ); // Increments ap to the next argument

    va_end(ap);

    return size;
}

// Interface to calculate fluxes and drag
void aerobulk::model( algorithm algo, double z_t, double z_u,
        std::vector<double>& sst, std::vector<double>& t_zt, std::vector<double>& q_zt,
        std::vector<double>& U_zu, std::vector<double>& slp,
        std::vector<double>& QL, std::vector<double>& QH, std::vector<double>& Cd_rho_U,
        std::vector<double>& evap,
        const std::vector<double>& rad_sw, const std::vector<double>& rad_lw)
{
    // Algorithm type
    std::string calgo = algorithm_to_string(algo);
    int l = calgo.size();

    // Check the length of the inputs and record it
	int m = aerobulk::check_sizes(5, sst.size(), t_zt.size(), q_zt.size(), U_zu.size(), slp.size());

    // Check the length of the optional inputs and record it
	int n = aerobulk::check_sizes(2, rad_sw.size(), rad_lw.size());

    // Set the size of the outputs
    QL.resize(m);
    QH.resize(m);
    Cd_rho_U.resize(m);
    evap.resize(m);

    // The actual function call - we need to send the addresses/pointer because it's a C interface to a Fortran routine
    if ( m == n )
        aerobulk_nextsim_skin( calgo.c_str(), &z_t, &z_u, &sst[0], &t_zt[0], &q_zt[0], &U_zu[0], &slp[0],
                &QL[0], &QH[0], &Cd_rho_U[0], &evap[0],
                &l, &m,
                &rad_sw[0], &rad_lw[0]);
    else
        aerobulk_nextsim_no_skin( calgo.c_str(), &z_t, &z_u, &sst[0], &t_zt[0], &q_zt[0], &U_zu[0], &slp[0],
                &QL[0], &QH[0], &Cd_rho_U[0], &evap[0],
                &l, &m);
}
