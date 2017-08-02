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
#include <boost/unordered_map.hpp>
#include <InterpFromGridToMeshx.h>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromMeshToGridx.h>
#include <BamgTriangulatex.h>
#include <Bamgx.h>
#include <iomanip>
#include <omp.h>
#include <gmshmesh.hpp>


#define PI M_PI

namespace Wim
{

namespace po = boost::program_options;
namespace fs = boost::filesystem;

template<typename T=float> class WimDiscr
{
    // ==========================================================================================
	typedef T value_type;
    typedef typename std::vector<value_type> value_type_vec;
    typedef typename std::vector<std::vector<value_type>*> value_type_vec_ptrs;
    typedef size_t size_type;
	typedef boost::multi_array<value_type, 2> array2_type;
    typedef boost::multi_array<value_type, 3> array3_type;
    typedef boost::multi_array<value_type, 4> array4_type;
    typedef typename array2_type::index index;
    
    //gmsh types
    typedef typename Nextsim::GmshMesh::point_type point_type;
    typedef typename Nextsim::GmshMesh::element_type element_type;
    typedef Nextsim::GmshMesh mesh_type;

    typedef struct IceInfo
    {
        // information describing ice fields
        // - ie that needed for interpMeshToGrid, interpMeshToMesh
        value_type_vec conc;        // total ice conc
        value_type_vec vol;         // total ice volume (conc*thickness)
        value_type_vec thick;       // ice thickness
        value_type_vec nfloes;      // N_floes = conc/Dfloe^2
        value_type_vec dfloe;       // max floe size
        value_type_vec mask;        // 1 if ice present, else 0
        value_type_vec broken;      // 1 if broken during this call to the wim
    } IceInfo;

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


    // ==========================================================================================

public:

    // ====================================================================================
    // public types
    typedef boost::unordered_map<std::string,value_type_vec>  unord_map_vecs_type;

    typedef struct MeshInfo
    {
        // information describing nextsim mesh
        // - only basic info
        // - ie that needed for interpMeshToGrid, interpMeshToMesh
        bool initialised;                       // initialised yet?
        int num_nodes;                          // number of nodes (not needed for structured grids)
        int num_elements;                       // number of elements 
        std::vector<int> index;                 // indices of nodes corresponding to the elements
        std::vector<int> element_connectivity;  // indices of nodes corresponding to the elements
        value_type_vec nodes_x;                 // x-coords of nodes (not needed for structured grids)
        value_type_vec nodes_y;                 // y-coords of nodes (not needed for structured grids)
        value_type_vec elements_x;              // x-coords of elements
        value_type_vec elements_y;              // y-coords of elements
        value_type_vec surface;                 // y-coords of elements
    } MeshInfo;

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
    // ====================================================================================


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

    //make a grid
    void gridProcessing();
    void gridProcessing(mesh_type const &mesh);
    void gridFromParameters();
    void gridPostProcessing();

    void saveGrid();
    void readGridFromFile();
    void readFromBinary(std::fstream &in, value_type_vec& in_array, int off = 0, std::ios_base::seekdir direction = std::ios::beg,
            int addx = 0, int addy = 0);
    void readDataFromFile(std::string const& filein);
    void exportResults(std::string const& output_type, value_type const& t_out) const;
    void testInterp(std::string const& output_type,
                    value_type const& t_out,
                    std::vector<std::vector<value_type>> const& vectors,
                    std::vector<std::string> const& names
                    ) const;
    void saveLog(value_type const& t_out) const;
    void saveOptionsLog();

    void init(int const nextsim_cpt=0);
    void init(mesh_type const &mesh,int const nextsim_cpt=0);

    void initConstant(int const nextsim_cpt);
    void assign();

    void update();

    void updateWaveMedium();

    void timeStep(int &lcpt);
    void doBreaking(BreakInfo const& breakinfo);

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

    //std::vector<bool> getBrokenMesh() const {return mesh_broken;}
    void getFsdMesh(value_type_vec &nfloes_out,value_type_vec &dfloe_out,value_type_vec &broken);
    void getFsdMesh(value_type_vec &nfloes_out,value_type_vec &dfloe_out, value_type_vec &broken,
            value_type_vec const & conc_tot, mesh_type const &mesh_in,value_type_vec const &um_in);

    value_type getModelTime(int lcpt=0) const {return M_update_time+lcpt*M_timestep;}

    // ========================================================================
    // breaking on mesh
    void setMesh( mesh_type const &mesh);
    void setMesh( mesh_type const &mesh,value_type_vec const &um);
    void setMesh( mesh_type const &mesh,value_type_vec const &um,BamgMesh *bamgmesh);
    void resetMesh( mesh_type const &mesh);
    void resetMesh( mesh_type const &mesh,value_type_vec const &um);
    value_type_vec relativeMeshDisplacement(mesh_type const &mesh_in,value_type_vec const &um_in) const;
    void setIceFields( std::vector<value_type> const& m_conc,  // conc
                       std::vector<value_type> const& m_vol, // ice vol or effective thickness (conc*thickness)
                       std::vector<value_type> const& m_nfloes,// Nfloes=conc/Dmax^2
                       bool const pre_regrid);
    void transformIce(IceInfo &ice_info);

    void clearMeshFields();
    void gridToPoints(
        value_type_vec_ptrs &output_data,
        value_type_vec_ptrs const &input_data,
        value_type_vec &Rx, value_type_vec &Ry);
    void meshToGrid(
        value_type_vec_ptrs &output_data,       //output data
        value_type_vec_ptrs const &input_data); //input data
    void meshToPoints(
        value_type_vec_ptrs &output_data,       //output data
        value_type_vec_ptrs const &input_data,  //input data
        value_type_vec &Rx,                     //location of output data (x-coord)
        value_type_vec &Ry);                    //location of output data (y-coord)

    unord_map_vecs_type returnFieldsElements(std::vector<std::string> const&fields,
            value_type_vec &xel, value_type_vec &yel);
    unord_map_vecs_type returnFieldsElements(std::vector<std::string> const&fields,
            mesh_type const &mesh_in,value_type_vec const &um_in);
    unord_map_vecs_type returnFieldsNodes(std::vector<std::string> const&fields,
            value_type_vec &xnod, value_type_vec &ynod);
    unord_map_vecs_type returnFieldsNodes(std::vector<std::string> const&fields,
            mesh_type const &mesh_in,value_type_vec const &um_in);

    void returnWaveStress(value_type_vec &M_tau, value_type_vec &xnod, value_type_vec &ynod);
    void returnWaveStress(value_type_vec &M_tau, mesh_type const &mesh_in,value_type_vec const &um_in);

    // ========================================================================

    WimGrid wimGrid(std::string const& units="m");

    void run();

    //===========================================================================
    //FSD: Dmax -> <D^moment> conversion
    void floeScaling(
          value_type const& dmax, int const& moment, value_type& dave);
    void floeScalingSmooth(
          value_type const& dmax, int const& moment, value_type& dave);
    //===========================================================================

    //===========================================================================
    //advection/attenuation
    void advectDirections( array2_type& Sdir, value_type_vec const& ag2d_eff);
    void advectDirectionsMesh( array2_type& Sdir, value_type_vec& ag2d_eff);
    void attenSimple(
          array2_type& Sdir, value_type_vec& Sfreq,
          value_type_vec& taux_omega,value_type_vec& tauy_omega,
          value_type_vec& sdx_omega,value_type_vec& sdy_omega,
          value_type_vec const& ag2d_eff);
    void attenIsotropic(array2_type& Sdir, value_type_vec& Sfreq,
          value_type_vec& taux_omega,value_type_vec& tauy_omega,
          value_type_vec& sdx_omega,value_type_vec& sdy_omega,
          value_type_vec const& ag2d_eff);
    void waveAdvWeno(
          value_type_vec& h, value_type_vec const& u, value_type_vec const& v);
    void weno3pdV2(
          value_type_vec const& gin, value_type_vec const& u, value_type_vec const& v,
          value_type_vec const& scuy, value_type_vec const& scvx,
          value_type_vec const& scp2i, value_type_vec const& scp2,
          value_type_vec& saoout);
    void padVar(value_type_vec const& u, value_type_vec& upad,std::string const& advopt_);
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

    value_type_vec getX() const { return X_array; }
    value_type_vec getY() const { return Y_array; }
    value_type_vec getSCUY() const { return SCUY_array; }
    value_type_vec getSCVX() const { return SCVX_array; }
    value_type_vec getSCP2() const { return SCP2_array; }
    value_type_vec getSCP2I() const { return SCP2I_array; }
    value_type_vec getLANDMASK() const { return LANDMASK_array; }

    std::string getWimGridFilename() const { return wim_gridfile; }
    //std::vector<int> getWimShape();

    std::vector<value_type> getTaux() const { return tau_x; }
    std::vector<value_type> getTauy() const { return tau_y; }
    std::vector<value_type> getStokesDriftx() const { return stokes_drift_x; }
    std::vector<value_type> getStokesDrifty() const { return stokes_drift_y; }
    std::vector<value_type> getNfloes() const { return wim_ice.nfloes; }


private:

    std::vector<int> iv_tmp = {};
    value_type_vec   vtv_tmp = {};
    MeshInfo mesh_info_tmp = {
            initialised             : false,
            num_nodes               : 0,
            num_elements            : 0,
            index                   : iv_tmp,
            element_connectivity    : iv_tmp,
            nodes_x                 : vtv_tmp,
            nodes_y                 : vtv_tmp,
            elements_x              : vtv_tmp,
            elements_y              : vtv_tmp,
            surface                 : vtv_tmp
    };

    po::variables_map vm;
    int nx, ny, nxext, nyext, nbdy, nbdx, nghost;
    int num_p_wim,num_q_wim,num_u_wim,num_v_wim;
    int M_itest;
    value_type_vec X_array, Y_array, SCUY_array, SCVX_array,
                SCP2_array, SCP2I_array, LANDMASK_array;
    std::vector<value_type> x_col,y_row;

    value_type cfl, dom, guess, Tmin, Tmax, gravity, om;
    value_type xmax, ym, x0, y0, dx, dy, x_edge, unifc, unifh,
               dfloe_pack_init, dfloe_pack_thresh, amin, amax;
    value_type rhowtr, rhoice, poisson, dmin, xi, fragility, cice_min, dfloe_miz_thresh,
               young, drag_rp, kice, kwtr, int_adm, modT, argR, argT, rhoi, rho, rhow;
    value_type fmin, fmax, df, epsc, sigma_c, vbf, vb, flex_rig_coeff;
    value_type M_timestep,duration;

    int nwavedirn, nwavefreq, advdim, ncs ,nt;
    bool ref_Hs_ice, atten, useicevel, M_steady, breaking;
    bool M_dump_diag;
    bool docoupling;
    std::string scatmod, advopt, fsdopt;
    std::string wim_gridfile;
    value_type_vec wavedir, wt_simp, wt_om, freq_vec, vec_period, wlng, ag, ap;
    value_type_vec Hs,Tp,mwd,wave_mask,M_steady_mask;

    //value_type_vec ice_mask, icec, iceh;
    value_type_vec swh_in_array,mwp_in_array,mwd_in_array,
                dave, atten_dim, damp_dim, ag2d_eff_temp;
    array2_type ag_eff, ap_eff, wlng_ice, atten_nond, damping, disp_ratio, sdf3d_dir_temp;
    array3_type sdf_dir, sdf_inc;

    value_type_vec S_freq, taux_om, tauy_om,
                stokes_drift_x_om, stokes_drift_y_om;
    value_type_vec hp;
    value_type_vec Fdmax, Ftaux, Ftauy, Fhs, Ftp;

    //std::vector<value_type> dfloe, nfloes;
    std::vector<value_type> tau_x, tau_y,stokes_drift_x,stokes_drift_y;//row-major order (C)
    //std::vector<value_type> mesh_x, mesh_y, mesh_conc, mesh_thick, mesh_dfloe;
    //std::vector<bool> mesh_broken;
    bool break_on_mesh;

    boost::mpi::timer chrono;
    std::string init_time_str;
    value_type M_restart_time = 0.;//time of restarting, relative to init_time
    value_type M_update_time  = 0;//time of start of call to wim.run(), relative to init_time
    value_type M_time_mesh_set = 0;//time of start of call to wim.run(), relative to init_time
    int M_cpt = 0;//global counter
    int M_num_elements;
    int max_threads;
    bool M_regular = false;
    bool M_initialised_ice = false;
    bool M_initialised_waves = false;

    MeshInfo M_wim_triangulation,nextsim_mesh,nextsim_mesh_old;
    IceInfo wim_ice, nextsim_ice;
    std::vector<int> wet_indices;

    // =========================================================================
    //to run WIM on nextsim mesh
    bool M_wim_on_mesh = false;
    value_type_vec M_UM;//displacement of mesh nodes between calls to wim.run()
                        //- for correction to group velocity at advection time
    // =========================================================================

};

} // namespace Wim

#endif
