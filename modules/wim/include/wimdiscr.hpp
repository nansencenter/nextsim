/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   wimdiscr.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug  3 11:53:19 2015
 */


#ifndef __WIMDISCR_H
#define __WIMDISCR_H 1

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/any.hpp>
#include <boost/format.hpp>
#include <boost/mpi/timer.hpp>
#include <iomanip>
#include <omp.h>

#ifdef __cplusplus
extern "C"
{
#endif
#include<RTparam_outer.h>
#ifdef __cplusplus
}
#endif


#define PI M_PI

namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace Wim
{
template<typename T=float> class WimDiscr
{
	typedef T value_type;
    typedef size_t size_type;
	typedef boost::multi_array<value_type, 2> array2_type;
    typedef boost::multi_array<value_type, 3> array3_type;
    typedef boost::multi_array<value_type, 4> array4_type;
    typedef typename array2_type::index index;

public:

    WimDiscr(po::variables_map const& vmin)
        :
        vm(vmin),
        nx(vm["wim.nx"].template as<int>()),
        ny(vm["wim.ny"].template as<int>())
    {}

    void gridProssessing();
	void readGridFromFile(std::string const& filein);
    void readDataFromFile(std::string const& filein);
    void exportResults(size_type const& timestp, value_type const& t_out) const;
    void init();

    void assign(std::vector<value_type> const& ice_c = std::vector<value_type>(),
                std::vector<value_type> const& ice_h = std::vector<value_type>(),
                std::vector<value_type> const& n_floes = std::vector<value_type>());

    void timeStep();

    void run(std::vector<value_type> const& ice_c = std::vector<value_type>(),
             std::vector<value_type> const& ice_h = std::vector<value_type>(),
             std::vector<value_type> const& n_floes = std::vector<value_type>());

    void floeScaling(value_type const& dmax, value_type& dave);
    void advAttenSimple(array3_type& Sdir, array2_type& Sfreq,array2_type& taux_omega,array2_type& tauy_omega, array2_type const& ag2d_eff);
    void advAttenIsotropic(array3_type& Sdir, array2_type& Sfreq,array2_type& taux_omega,array2_type& tauy_omega, array2_type const& ag2d_eff);
    void waveAdvWeno(array2_type& h, array2_type const& u, array2_type const& v);
    void weno3pdV2(array2_type const& gin, array2_type const& u, array2_type const& v, array2_type const& scuy,
                   array2_type const& scvx, array2_type const& scp2i, array2_type const& scp2, array2_type& saoout);

    void padVar(array2_type const& u, array2_type& upad);
    void calcMWD();

    array2_type getX() const { return X_array; }
    array2_type getY() const { return Y_array; }
    array2_type getSCUY() const { return SCUY_array; }
    array2_type getSCVX() const { return SCVX_array; }
    array2_type getSCP2() const { return SCP2_array; }
    array2_type getSCP2I() const { return SCP2I_array; }
    array2_type getLANDMASK() const { return LANDMASK_array; }

    std::vector<value_type> getTaux() const { return tau_x; }
    std::vector<value_type> getTauy() const { return tau_y; }
    std::vector<value_type> getNFloes() const { return nfloes; }

private:

    po::variables_map const& vm;
    int nx, ny, nxext, nyext, nbdy, nbdx;
    array2_type X_array, Y_array, SCUY_array, SCVX_array, SCP2_array, SCP2I_array, LANDMASK_array;

    value_type cfl, dom, guess, Hs_inc, Tp_inc, mwd_inc, Tmin, Tmax, gravity, om;
    value_type xmax, ym, x0, y0, dx, dy, x_edge, unifc, unifh, dfloe_pack_init, dfloe_pack_thresh, amin, amax, dt;
    value_type rhowtr, rhoice, poisson, dmin, xi, fragility, young, visc_rp, kice, kwtr, int_adm, modT, argR, argT, rhoi, rho, rhow;
    value_type fmin, fmax, df, epsc, sigma_c, vbf, vb, flex_rig_coeff;

    int nwavedirn, nwavefreq, advdim, ncs;
    bool ref_Hs_ice, atten, icevel, steady, breaking;
    std::string scatmod, advopt;
    std::vector<value_type> wavedir, wt_simp, wt_om, freq_vec, vec_period, wlng, ag, ap;

    array2_type steady_mask, wave_mask, ice_mask, wtr_mask, icec, iceh, atten_dim, damp_dim, ag2d_eff_temp, mwd, Hs, Tp;
    array3_type ag_eff, ap_eff, wlng_ice, atten_nond, damping, disp_ratio, sdf3d_dir_temp;
    array4_type sdf_dir, sdf_inc;

    array2_type S_freq, taux_om, tauy_om;
    array2_type hp;
    array2_type Fdmax, Ftaux, Ftauy, Fhs, Ftp;

    //array2_type dfloe, nfloes, tau_x, tau_y;
    std::vector<value_type> dfloe, nfloes, tau_x, tau_y;

    boost::mpi::timer chrono;

};

} // namespace Wim

#endif
