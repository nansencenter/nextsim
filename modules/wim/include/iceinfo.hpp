/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   wimtypes.hpp
 * @author Timothy Williams <timothy.williams@nersc.no>
 * @date   Mon Aug  10 2017
 */

/**
 * Name conventions:
 * M_... = global variable
 * T_... = custom typename
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

template<typename T=float> class IceParams
{

    typedef T T_val;

    //private functions
    void init();

public:

    // ===============================================================
    // parameters
    T_val M_gravity          = 9.81;
    T_val M_rhowtr           = 1025.;
    T_val M_rhoice           = 922.5;
    T_val M_poisson          = 0.3;
    T_val M_dmin             = 20.;
    T_val M_xi               = 2.;
    T_val M_fragility        = 0.9;
    T_val M_dfloe_miz_thresh = 200.;
    T_val M_vbf              = .1;
    // ===============================================================


    // ===============================================================
    // can reset with vm
    T_val M_drag_rp            = 0.;
    T_val M_young              = 5.49e9;
    T_val M_cice_min           = .05;
    T_val M_dfloe_pack_init    = 300.;
    T_val M_dfloe_pack_thresh  = 350.;
    std::string M_fsdopt       = "PowerLawSmooth";
    // ===============================================================


    // ===============================================================
    //dependent parameters - set with init()
    T_val M_epsc;
    T_val M_sigma_c;
    // ===============================================================


    // ===============================================================
    //destructor
    ~IceParams() {}

    //constructors
    IceParams() {}
    IceParams(po::variables_map const& vm) { this->update(vm); }
    // ===============================================================


    // ===============================================================
    // public functions
    void update(po::variables_map const& vm);
    std::vector<double> getAttenParams();

    // get parameters
    T_val gravity()         const { return M_gravity; }
    T_val rhowtr()          const { return M_rhowtr; }
    T_val rhoice()          const { return M_rhoice; }
    T_val poisson()         const { return M_poisson; }
    T_val dmin()            const { return M_dmin; }
    T_val xi()              const { return M_gravity; }
    T_val fragility()       const { return M_fragility; }
    T_val dfloeMizThresh()  const { return M_dfloe_miz_thresh; }
    T_val vbf()             const { return M_vbf; }
    T_val dragRp()          const { return M_drag_rp; }
    T_val young()           const { return M_young; }
    T_val ciceMin()         const { return M_cice_min; }
    T_val epsc()            const { return M_epsc; }
    T_val sigmac()          const { return M_sigma_c; }
    T_val dfloePackInit()   const { return M_dfloe_pack_init; }
    T_val dfloePackThresh() const { return M_dfloe_pack_thresh; }
    std::string fsdopt()    const { return M_fsdopt; }

};//class IceParams


template<typename T=float> class IceInfo
{

    typedef T T_val;
    typedef typename std::vector<T_val> T_val_vec;
    typedef IceParams<T_val>*           T_icep_ptr;
    T_icep_ptr M_ice_params;//pointer to IceParams object - will be wim.M_ice_params

    int M_num_elements = 0;
    int M_max_threads  = 0;

    //private functions
    void init();

public:

    // ===============================================================
    //destructor
    ~IceInfo() {}

    //constructors
    IceInfo() {}
    IceInfo(T_icep_ptr ice_params)
        :
        M_ice_params(ice_params)
    {
        this->init();
    }
    // ===============================================================

    // ===============================================================
    // public functions
    void clearFields();
    void setFields();
    void setFields(T_val_vec const& conc,
            T_val_vec const& vol,
            T_val_vec const& nfloes);

    void doBreaking(T_val_vec const& mom0,T_val_vec const& mom2,
            T_val_vec const& var_strain);

    T_val_vec dfloeToNfloes(T_val_vec const& dfloe_in,
            T_val_vec const& conc_in);
    T_val dfloeToNfloes(T_val const& dfloe_in,
            T_val const& conc_in);
    T_val nfloesToDfloe(T_val const& nfloes_in,
            T_val const& conc_in);
    T_val_vec nfloesToDfloe(T_val_vec const& nfloes_in,
            T_val_vec const& conc_in);

    void floeScaling(
            T_val const& dmax,int const& moment, T_val& dave);
    void floeScalingSmooth(
            T_val const& dmax,int const& moment, T_val& dave);

    void getDave(T_val_vec &dave,int const& moment);
    // ===============================================================


    // ===============================================================
    //public variables
    T_val_vec M_conc;
    T_val_vec M_vol;
    T_val_vec M_broken;
    T_val_vec M_thick;
    T_val_vec M_dfloe;
    T_val_vec M_nfloes;
    T_val_vec M_mask;
    // ===============================================================

};//class IceInfo

}//namespace Wim
#endif
