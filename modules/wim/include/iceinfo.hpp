/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   wimtypes.hpp
 * @author Timothy Williams <timothy.williams@nersc.no>
 * @date   Mon Aug  10 2017
 */

#ifndef __ICEINFO_H
#define __ICEINFO_H 1

#include <cmath>
#include <iostream>
#include <boost/program_options.hpp>
#include <omp.h>

#ifndef PI
#define PI M_PI
#endif

namespace Wim
{

namespace po = boost::program_options;

template<typename T=float> class IceInfo
{

    typedef T value_type;
    typedef typename std::vector<value_type>        value_type_vec;
    typedef typename std::vector<value_type_vec*>   value_type_vec_ptrs;
    po::variables_map vm;
    int M_num_elements = 0;
    int M_max_threads  = 0;

    // ===============================================================
    // parameters
    value_type M_gravity          = 9.81;
    value_type M_rhowtr           = 1025.;
    value_type M_rhoice           = 922.5;
    value_type M_poisson          = 0.3;
    value_type M_dmin             = 20.;
    value_type M_xi               = 2.;
    value_type M_fragility        = 0.9;
    value_type M_dfloe_miz_thresh = 200.;
    value_type M_vbf              = .1;
    // ===============================================================


    // ===============================================================
    // can reset with vm
    value_type M_drag_rp            = 0.;
    value_type M_young              = 5.49e9;
    value_type M_cice_min           = .05;
    value_type M_dfloe_pack_init    = 300.;
    value_type M_dfloe_pack_thresh  = 350.;
    std::string M_fsdopt            = "PowerLawSmooth";
    // ===============================================================


    // ===============================================================
    //dependent parameters - set with init()
    value_type M_epsc;
    value_type M_sigma_c;
    // ===============================================================


    //private functions
    void init();

public:

    //destructor
    ~IceInfo();

    //constructors
    IceInfo();
    IceInfo(po::variables_map const& vm);

    //explicit copy operator
    //IceInfo(IceInfo<value_type> const& tmp);

    // ===============================================================
    // public functions
    void update(po::variables_map const& vm);
    void clearFields();
    void setFields();
    void setFields(value_type_vec const& conc,
            value_type_vec const& vol,
            value_type_vec const& nfloes);

    void doBreaking(value_type_vec const& mom0,value_type_vec const& mom2,
            value_type_vec const& var_strain);

    value_type_vec dfloeToNfloes(value_type_vec const& dfloe_in,
            value_type_vec const& conc_in);
    value_type dfloeToNfloes(value_type const& dfloe_in,
            value_type const& conc_in);
    value_type nfloesToDfloe(value_type const& nfloes_in,
            value_type const& conc_in);
    value_type_vec nfloesToDfloe(value_type_vec const& nfloes_in,
            value_type_vec const& conc_in);

    void floeScaling(
            value_type const& dmax,int const& moment, value_type& dave);
    void floeScalingSmooth(
            value_type const& dmax,int const& moment, value_type& dave);

    void getDave(value_type_vec &dave,int const& moment);
    std::vector<double> getAttenParams();

    value_type gravity() const { return M_gravity; }
    value_type rhowtr() const { return M_rhowtr; }
    value_type rhoice() const { return M_rhoice; }
    value_type poisson() const { return M_poisson; }
    value_type dmin() const { return M_dmin; }
    value_type xi() const { return M_gravity; }
    value_type fragility() const { return M_fragility; }
    value_type dfloeMizThresh() const { return M_dfloe_miz_thresh; }
    value_type vbf() const { return M_vbf; }
    value_type dragRp() const { return M_drag_rp; }
    value_type young() const { return M_young; }
    value_type ciceMin() const { return M_cice_min; }
    value_type epsc() const { return M_epsc; }
    value_type sigmac() const { return M_sigma_c; }
    value_type dfloePackInit() const { return M_dfloe_pack_init; }
    value_type dfloePackThresh() const { return M_dfloe_pack_thresh; }
    std::string fsdopt() const { return M_fsdopt; }
    // ===============================================================


    // ===============================================================
    //public variables
    value_type_vec M_conc;
    value_type_vec M_vol;
    value_type_vec M_broken;
    value_type_vec M_thick;
    value_type_vec M_dfloe;
    value_type_vec M_nfloes;
    value_type_vec M_mask;
    // ===============================================================
};

}//namespace Wim
#endif
