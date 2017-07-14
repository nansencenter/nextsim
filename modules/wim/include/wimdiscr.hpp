/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

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
#include <InterpFromGridToMeshx.h>
#include <date.hpp>
#include <iomanip>
#include <omp.h>

#ifdef __cplusplus
extern "C"
{
#endif
#include<RTparam_outer.h>
#include <mapx.h>
#ifdef __cplusplus
}
#endif

#define PI M_PI

namespace Wim
{

namespace po = boost::program_options;
namespace fs = boost::filesystem;

template<typename T=float> class WimDiscr
{
	typedef T value_type;
    typedef typename std::vector<value_type> value_type_vec;
    typedef size_t size_type;
	typedef boost::multi_array<value_type, 2> array2_type;
    typedef boost::multi_array<value_type, 3> array3_type;
    typedef boost::multi_array<value_type, 4> array4_type;
    typedef typename array2_type::index index;

public:

    typedef struct BreakInfo
    {
        // information needed for breaking
        value_type conc;   // concentration
        value_type thick;  // thickness
        value_type mom0;
        value_type mom2;
        value_type var_strain;
        mutable value_type dfloe;
        mutable bool broken;
    } BreakInfo;

    typedef struct WimGrid
    {
        // information describing wim grid
        int nx;
        int ny;
        value_type dx;
        value_type dy;
        std::vector<value_type> X;
        std::vector<value_type> Y;
        std::vector<value_type> x;
        std::vector<value_type> y;
    } WimGrid;


public:

    WimDiscr()
        :
        vm(),
        nx(),
        ny()
    {}

    WimDiscr(po::variables_map const& vmIn)
        :
        vm(vmIn),
        nx(vm["wim.nx"].template as<int>()),
        ny(vm["wim.ny"].template as<int>())
    {}

    void gridProcessing();
    void saveGrid();
    void readGridFromFile();
    void readFromBinary(std::fstream &in, array2_type& in_array, int off = 0, std::ios_base::seekdir direction = std::ios::beg, int addx = 0, int addy = 0);
    void readDataFromFile(std::string const& filein);
    void exportResults(std::string const& output_type, value_type const& t_out) const;
    void testInterp(std::string const& output_type,
                    value_type const& t_out,
                    std::vector<std::vector<value_type>> const& vectors,
                    std::vector<std::string> const& names
                    ) const;
    void saveLog(value_type const& t_out) const;

    void init(int const nextsim_cpt=0);
    void assign();

    void update(std::vector<value_type> const& icec_in = std::vector<value_type>(),
                std::vector<value_type> const& iceh_in = std::vector<value_type>(),
                std::vector<value_type> const& nfloes_in = std::vector<value_type>(),
                std::vector<value_type> const& swh_in = std::vector<value_type>(),
                std::vector<value_type> const& mwp_in = std::vector<value_type>(),
                std::vector<value_type> const& mwd_in = std::vector<value_type>());

    void updateWaveMedium();

    void timeStep();

    void doBreaking(BreakInfo const& breakinfo);

    void setMesh(std::vector<value_type> const& m_rx,
            std::vector<value_type> const& m_ry,
                 std::vector<value_type> const& m_conc,
                 std::vector<value_type> const& m_thick,
                 std::vector<value_type> const& m_nfloes,
                 std::string const& units="km");

    value_type nfloesToDfloe(
                 value_type const& m_nfloes,
                 value_type const& m_conc);
    std::vector<value_type> nfloesToDfloe(
                 std::vector<value_type> const& m_nfloes,
                 std::vector<value_type> const& m_conc);
    
    value_type dfloeToNfloes(
                 value_type const& m_dfloe,
                 value_type const& m_conc);
    std::vector<value_type> dfloeToNfloes(
                 std::vector<value_type> const& m_dfloe,
                 std::vector<value_type> const& m_conc);

    std::vector<bool> getBrokenMesh() const {return mesh_broken;}
    std::vector<value_type> getNfloesMesh();

    void clearMesh();

    WimGrid wimGrid(std::string const& units="m");

    void run(std::vector<value_type> const& icec_in = std::vector<value_type>(),
             std::vector<value_type> const& iceh_in = std::vector<value_type>(),
             std::vector<value_type> const& nfloes_in = std::vector<value_type>(),
             std::vector<value_type> const& swh_in = std::vector<value_type>(),
             std::vector<value_type> const& mwp_in = std::vector<value_type>(),
             std::vector<value_type> const& mwd_in = std::vector<value_type>());

    //===========================================================================
    //FSD: Dmax -> <D^moment> conversion
    void floeScaling(
          value_type const& dmax, int const& moment, value_type& dave);
    void floeScalingSmooth(
          value_type const& dmax, int const& moment, value_type& dave);
    //===========================================================================

    //===========================================================================
    //advection/attenuation
    void advAttenSimple(
          array3_type& Sdir, array2_type& Sfreq,
          array2_type& taux_omega,array2_type& tauy_omega,
          array2_type& sdx_omega,array2_type& sdy_omega,
          array2_type const& ag2d_eff);
    void advAttenIsotropic(array3_type& Sdir, array2_type& Sfreq,
          array2_type& taux_omega,array2_type& tauy_omega,
          array2_type& sdx_omega,array2_type& sdy_omega,
          array2_type const& ag2d_eff);
    void waveAdvWeno(
          array2_type& h, array2_type const& u, array2_type const& v);
    void weno3pdV2(
          array2_type const& gin, array2_type const& u, array2_type const& v,
          array2_type const& scuy, array2_type const& scvx,
          array2_type const& scp2i, array2_type const& scp2,
          array2_type& saoout);
    void padVar(array2_type const& u, array2_type& upad,std::string const& advopt_);
    //===========================================================================


    void calcMWD();
    void idealWaveFields(value_type const xfac);
    void idealIceFields (value_type const xfac);
    void inputWaveFields(value_type_vec const& swh_in,
            value_type_vec const& mwp_in,
            value_type_vec const& mwd_in);
    void setIncWaveSpec();
    void inputIceFields(value_type_vec const& icec_in,
            value_type_vec const& iceh_in,
            value_type_vec const& nfloes_in);
    //void getWimCenters(value_type& x,value_type& y,value_type const& rotangle);

    value_type thetaDirFrac(value_type const& th1_, value_type const& dtheta_, value_type const& mwd_);
    value_type thetaInRange(value_type const& th_, value_type const& th1, bool const& close_on_right=false);

    array2_type getX() const { return X_array; }
    array2_type getY() const { return Y_array; }
    array2_type getSCUY() const { return SCUY_array; }
    array2_type getSCVX() const { return SCVX_array; }
    array2_type getSCP2() const { return SCP2_array; }
    array2_type getSCP2I() const { return SCP2I_array; }
    array2_type getLANDMASK() const { return LANDMASK_array; }

    std::string getWimGridFilename() const { return wim_gridfile; }
    //std::vector<int> getWimShape();

    std::vector<value_type> getTaux() const { return tau_x; }
    std::vector<value_type> getTauy() const { return tau_y; }
    std::vector<value_type> getStokesDriftx() const { return stokes_drift_x; }
    std::vector<value_type> getStokesDrifty() const { return stokes_drift_y; }
    std::vector<value_type> getNfloes() const { return nfloes; }


private:

    po::variables_map vm;
    int nx, ny, nxext, nyext, nbdy, nbdx, nghost;
    int wim_itest, wim_jtest;
    array2_type X_array, Y_array, SCUY_array, SCVX_array,
                SCP2_array, SCP2I_array, LANDMASK_array;
    std::vector<value_type> x_col,y_row;

    value_type cfl, dom, guess, Hs_inc, Tp_inc, mwd_inc, Tmin, Tmax, gravity, om;
    value_type xmax, ym, x0, y0, dx, dy, x_edge, unifc, unifh,
               dfloe_pack_init, dfloe_pack_thresh, amin, amax;
    value_type rhowtr, rhoice, poisson, dmin, xi, fragility,
               young, visc_rp, kice, kwtr, int_adm, modT, argR, argT, rhoi, rho, rhow;
    value_type fmin, fmax, df, epsc, sigma_c, vbf, vb, flex_rig_coeff;
    value_type dt,duration;

    int nwavedirn, nwavefreq, advdim, ncs ,nt;
    bool ref_Hs_ice, atten, useicevel, steady, breaking, dumpDiag;
    bool docoupling;
    std::string scatmod, advopt, fsdopt;
    std::string wim_gridfile;
    std::vector<value_type> wavedir, wt_simp, wt_om, freq_vec, vec_period, wlng, ag, ap;
    std::vector<value_type> Hs,Tp,mwd,wave_mask;

    array2_type steady_mask, ice_mask, wtr_mask,
                icec, iceh, swh_in_array,mwp_in_array,mwd_in_array,
                dave, atten_dim, damp_dim, ag2d_eff_temp;
    array3_type ag_eff, ap_eff, wlng_ice, atten_nond, damping, disp_ratio, sdf3d_dir_temp;
    array4_type sdf_dir, sdf_inc;

    array2_type S_freq, taux_om, tauy_om,
                stokes_drift_x_om, stokes_drift_y_om;
    array2_type hp;
    array2_type Fdmax, Ftaux, Ftauy, Fhs, Ftp;

    std::vector<value_type> dfloe, nfloes, tau_x, tau_y,stokes_drift_x,stokes_drift_y;//row-major order (C)
    std::vector<value_type> mesh_x, mesh_y, mesh_conc, mesh_thick, mesh_dfloe;
    std::vector<bool> mesh_broken;
    bool break_on_mesh;
    int mesh_num_elements;

    boost::mpi::timer chrono;
    std::string init_time_str;
    value_type restart_time_shift;
    int cpt;

};

} // namespace Wim

#endif
