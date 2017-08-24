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
#include <iceinfo.hpp>

#ifdef PI
#undef PI
#define PI M_PI
#endif

namespace Wim
{

namespace po = boost::program_options;
namespace fs = boost::filesystem;

template<typename T=float> class WimDiscr
{
    // ==========================================================================================
    typedef T value_type;
    typedef typename std::vector<value_type>        value_type_vec;
    typedef typename std::vector<value_type_vec*>   value_type_vec_ptrs;
    typedef typename std::vector<value_type_vec>    value_type_vec2d;//vector of vectors
    typedef typename std::vector<value_type_vec2d>  value_type_vec3d;//vector of vectors of vectors
    typedef size_t size_type;
    typedef boost::multi_array<value_type, 2> array2_type;
    typedef boost::multi_array<value_type, 3> array3_type;
    typedef boost::multi_array<value_type, 4> array4_type;
    typedef typename array2_type::index index;
    
    //gmsh types
    typedef typename Nextsim::GmshMesh::point_type point_type;
    typedef typename Nextsim::GmshMesh::element_type element_type;
    typedef Nextsim::GmshMesh mesh_type;

    typedef Wim::IceInfo<double> ice_type;
#if 0
    typedef struct IceInfo
    {
        // information describing ice fields
        // - ie that needed for interpMeshToGrid, interpMeshToMesh
        mutable value_type_vec conc;        // total ice conc
        mutable value_type_vec vol;         // total ice volume (conc*thickness)
        mutable value_type_vec thick;       // ice thickness
        mutable value_type_vec nfloes;      // N_floes = conc/Dfloe^2
        mutable value_type_vec dfloe;       // max floe size
        mutable value_type_vec mask;        // 1 if ice present, else 0
        mutable value_type_vec broken;      // 1 if broken during this call to the wim
    } IceInfo;
#endif

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
    typedef boost::unordered_map<std::string,value_type_vec*> unord_map_vec_ptrs_type;
    typedef boost::unordered_map<std::string,std::string>     unord_map_type;

    typedef struct MeshInfo
    {
        // information describing nextsim mesh
        // - only basic info
        // - ie that needed for interpMeshToGrid, interpMeshToMesh
        bool initialised;                           // initialised yet?
        int num_nodes;                              // number of nodes (not needed for structured grids)
        int num_elements;                           // number of elements 
        int max_node_el_conn;                       // max number of elements connected to a node
        std::vector<int> index;                     // indices of nodes corresponding to the elements
        std::vector<int> element_connectivity;      // indices of neighbouring elements
        std::vector<int> node_element_connectivity; // indices of elements connected to each node
        std::vector<bool> mask_dirichlet;           // is the node on a coast (edge is coastal if this is true for both nodes)?
        value_type_vec nodes_x;                     // x-coords of nodes (not needed for structured grids)
        value_type_vec nodes_y;                     // y-coords of nodes (not needed for structured grids)
        value_type_vec elements_x;                  // x-coords of elements
        value_type_vec elements_y;                  // y-coords of elements
        value_type_vec surface;                     // surface area of elements
        std::vector<int> id;                        // id's of nodes
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

    void exportResults(std::string const& output_type);
    void exportResultsGrid(unord_map_vec_ptrs_type& extract_fields,
            std::vector<std::string> const& strings);
    void exportResultsMesh(unord_map_vec_ptrs_type & extract_fields,
            std::vector<std::string> const &filenames,
            bool export_mesh=true, bool export_fields=true);
    void exportMesh(std::string const &filename);
    void testMesh();

    void testInterp(std::string const& output_type,
                    value_type const& t_out,
                    std::vector<std::vector<value_type>> const& vectors,
                    std::vector<std::string> const& names
                    ) const;
    void saveLog(value_type const& t_out) const;
    void saveOptionsLog();

    void init1(int const& nextsim_cpt);
    void init(int const& nextsim_cpt=0);
    void init(mesh_type const &mesh,int const& nextsim_cpt=0);
    void init(mesh_type const &mesh,BamgMesh* bamgmesh,int const& flag_fix,int const& nextsim_cpt=0);

    void init2();
    void initConstant(int const& nextsim_cpt);
    void assign();
    void assignSpatial();

    void update();

    void updateWaveMedium();

    void timeStep();
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

    void getFsdMesh(value_type_vec &nfloes_out,value_type_vec &dfloe_out,value_type_vec &broken);
    void getFsdMesh(value_type_vec &nfloes_out,value_type_vec &dfloe_out, value_type_vec &broken,
            value_type_vec const & conc_tot, mesh_type const &mesh_in,value_type_vec const &um_in);
    void getFsdMesh(value_type_vec &nfloes_out,value_type_vec &dfloe_out, value_type_vec &broken,
            value_type_vec const & conc_tot, mesh_type const &mesh_in);

    //time in seconds from init_time
    value_type getModelTime(int lcpt=0) const {return M_update_time+lcpt*M_timestep;}

    //time in days from ref time (1900-1-1)
    value_type getNextsimTime() const;

    // ========================================================================
    // breaking on mesh
    void setMesh( mesh_type const &mesh);
    void setMesh( mesh_type const &mesh,value_type_vec const &um);
    void setMesh( mesh_type const &mesh,value_type_vec const &um,BamgMesh* bamgmesh,int const& flag_fix,bool const& regridding=false);
    void setMesh( mesh_type const &mesh,BamgMesh* bamgmesh,int const& flag_fix,bool const& regridding=false);

    void resetMesh( mesh_type const &mesh);
    void resetMesh( mesh_type const &mesh,value_type_vec const &um);

    void setIceFields( std::vector<value_type> const& m_conc,  // conc
                       std::vector<value_type> const& m_vol, // ice vol or effective thickness (conc*thickness)
                       std::vector<value_type> const& m_nfloes,// Nfloes=conc/Dmax^2
                       bool const pre_regrid);
    //void transformIce(IceInfo &ice_info);

    value_type_vec getRelativeMeshDisplacement(mesh_type const &mesh_in) const;
    value_type_vec getRelativeMeshDisplacement(mesh_type const &mesh_in,value_type_vec const &um_in) const;
    void updateWaveSpec( mesh_type const &mesh);
    void updateWaveSpec( mesh_type const &mesh,value_type_vec const &um);

    void clearMeshFields() {nextsim_ice.clearFields();}
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
    void elementsToNodes(
        value_type_vec_ptrs &output_data,       //output data
        value_type_vec_ptrs const &input_data); //input data

    unord_map_vecs_type returnFieldsElements(std::vector<std::string> const&fields,
            value_type_vec &xel, value_type_vec &yel, value_type_vec const&surface_fac);
    unord_map_vecs_type returnFieldsElements(std::vector<std::string> const&fields,
            mesh_type const &mesh_in,value_type_vec const &um_in);
    unord_map_vecs_type returnFieldsElements(std::vector<std::string> const&fields,
            mesh_type const &mesh_in);

    unord_map_vecs_type returnFieldsNodes(std::vector<std::string> const&fields,
            value_type_vec &xnod, value_type_vec &ynod);
    unord_map_vecs_type returnFieldsNodes(std::vector<std::string> const&fields,
            mesh_type const &mesh_in,value_type_vec const &um_in);
    unord_map_vecs_type returnFieldsNodes(std::vector<std::string> const&fields,
            mesh_type const &mesh_in);

    value_type_vec getSurfaceFactor(mesh_type const &mesh_in);

    void returnWaveStress(value_type_vec &M_tau, value_type_vec &xnod, value_type_vec &ynod);
    void returnWaveStress(value_type_vec &M_tau, mesh_type const &mesh_in,value_type_vec const &um_in);
    void returnWaveStress(value_type_vec &M_tau, mesh_type const &mesh_in);

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
    void advectDirections( value_type_vec2d& Sdir, value_type_vec const& ag2d_eff);
    void advectDirectionsMesh( value_type_vec2d& Sdir, value_type_vec& ag2d_eff,
            value_type_vec const& boundary_vals);
    void attenSimple(
            value_type_vec2d& Sdir, value_type_vec& Sfreq,
            value_type_vec& taux_omega,value_type_vec& tauy_omega,
            value_type_vec& sdx_omega,value_type_vec& sdy_omega,
            value_type_vec const& ag2d_eff);
    void attenIsotropic(
            value_type_vec2d& Sdir, value_type_vec& Sfreq,
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
    void padVar(value_type_vec const& u, value_type_vec& upad,
            std::string const& advopt_,bool const&steady=false);
    //===========================================================================

    void intWaveSpec();//TODO test this
    void intDirns(value_type_vec2d const& Sdir, value_type_vec& Sfreq,
            value_type_vec& sdx_omega, value_type_vec& sdy_omega);

    void idealWaveFields(value_type const xfac);
    void idealIceFields (value_type const xfac);
    void inputWaveFields(value_type_vec const& swh_in,
            value_type_vec const& mwp_in,
            value_type_vec const& mwd_in);
    void setIncWaveSpec(value_type_vec const& wave_mask);
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

    //for use at regridding time (M_wim_on_mesh)
    value_type_vec getMeshDisplacement() const { return M_UM; }
    value_type_vec3d getWaveSpec() const { return M_sdf_dir; }
    void setRelativeMeshDisplacement(value_type_vec const&um_in) { M_UM = um_in; return; }
    void setWaveSpec(value_type_vec3d const&sdf_in)
    {
        // reset wave spectrum after regrid;
        // also integrate if now so that initial wave fields are correct
        // and so they are ready for export;
        M_sdf_dir = sdf_in;
        this->intWaveSpec();
        return;
    }

    std::string getWimGridFilename() const { return wim_gridfile; }
    //std::vector<int> getWimShape();


private:

    po::variables_map vm;
    int nx, ny, nxext, nyext, nbdy, nbdx, nghost;
    int num_p_wim,num_q_wim,num_u_wim,num_v_wim;
    int M_itest;
    value_type_vec X_array, Y_array, SCUY_array, SCVX_array,
                SCP2_array, SCP2I_array, LANDMASK_array;
    std::vector<value_type> x_col,y_row;

    value_type M_cfl, M_length_cfl, M_max_cg;
    value_type M_current_time;
    value_type Tmin, Tmax;// gravity;
    value_type xmax, ym, x0, y0, dx, dy, x_edge;
    value_type M_dfloe_pack_init;
    //value_type dfloe_pack_thresh;
    //value_type rhowtr, rhoice, poisson, dmin, xi, fragility, cice_min, dfloe_miz_thresh,
               //young, drag_rp;
    //value_type epsc, sigma_c, vbf;
    //value_type kice, kwtr, int_adm, modT, argR, argT, rhoi, rho, rhow;//TODO these are temporary?
    value_type M_timestep,M_duration;
    int M_num_timesteps;

    int nwavedirn, nwavefreq, M_advdim;
    bool M_ref_Hs_ice, M_atten, M_useicevel, M_steady, M_breaking;
    bool M_dump_diag;
    bool docoupling;
    std::string M_scatmod, M_advopt;// fsdopt;
    std::string wim_gridfile;

    //dimension of wavedir
    value_type_vec M_wavedir, M_quadrature_wt_dir;
    
    //dimension of freq
    value_type_vec M_quadrature_wt_freq, M_freq_vec;
    value_type_vec M_wlng_wtr, M_ag_wtr, M_ap_wtr;

    //dimension of space
    value_type_vec M_Hs,M_Tp,M_mwd;
    value_type_vec M_swh_in,M_mwp_in,M_mwd_in;
    value_type_vec M_dave;

    //depend on freq and position
    value_type_vec2d M_ag_eff, M_agnod_eff, M_ap_eff, M_wlng_ice, M_atten_nond, M_damping, M_disp_ratio;

    //depend on freq, dirn and position
    value_type_vec2d M_open_boundary_vals;
    value_type_vec3d M_sdf_dir;
    value_type_vec3d M_sdf_dir_inc;//
    value_type_vec M_steady_mask;//TODO for grid advection, make steady a boundary condition

    //these are only temporary vectors, but they are global in order to
    //save creating and destroying them extremely often
    value_type_vec Mtmp_sdf_freq;                                   //for mom0 integral
    value_type_vec Mtmp_taux_om, Mtmp_tauy_om;                      //for taux,tau_y integrals
    value_type_vec Mtmp_stokes_drift_x_om, Mtmp_stokes_drift_y_om;  //for mwd_x,mwd_y, stokes_drift_x,stokes_drift_y integrals
    value_type_vec Mtmp_mom0,Mtmp_mom2,Mtmp_var_strain;
    value_type_vec Mtmp_mom0w,Mtmp_mom2w;
    value_type_vec Mtmp_atten_dim, Mtmp_damp_dim;

    value_type_vec M_mwd_x, M_mwd_y, M_tau_x, M_tau_y,M_stokes_drift_x,M_stokes_drift_y;//row-major order (C)

    int M_max_threads;
    boost::mpi::timer chrono;

    std::string M_init_time_str;
    value_type M_restart_time = 0.;// time of restarting, relative to init_time
    value_type M_update_time  = 0; // time of start of call to wim.run(), relative to init_time
    value_type M_time_mesh_set = 0;// time of start of call to wim.run(), relative to init_time
    int M_cpt = 0;//global counter
    int M_num_elements;
    bool M_regular = false;
    bool M_initialised_ice = false;
    bool M_initialised_waves = false;

    int M_nb_export_nextwim = 0;
    int M_nb_export_inc     = 0;
    int M_nb_export_init    = 0;
    int M_nb_export_prog    = 0;
    int M_nb_export_final   = 0;
    int M_nb_mesh_test      = 0;

    MeshInfo M_wim_triangulation,nextsim_mesh,nextsim_mesh_old;

    
    ice_type wim_ice, nextsim_ice;
    std::vector<int> wet_indices;

    bool M_break_on_mesh = false;// do breaking on nextsim mesh as well as on grid
    bool M_wim_on_mesh   = false;// to run WIM on nextsim mesh
    
    value_type_vec M_UM;// displacement of mesh nodes between calls to wim.run()
                        // - only used if running WIM on nextsim mesh
                        // - for correction to group velocity at advection time
    // =========================================================================

    MeshInfo mesh_info_tmp = {
            initialised               : false,
            num_nodes                 : 0,
            num_elements              : 0,
            max_node_el_conn          : 0,
            index                     : {},
            element_connectivity      : {},
            node_element_connectivity : {},
            mask_dirichlet            : {},
            nodes_x                   : {},
            nodes_y                   : {},
            elements_x                : {},
            elements_y                : {},
            surface                   : {},
            id                        : {}
    };

};

} // namespace Wim

#endif
