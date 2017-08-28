/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   wimdiscr.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Timothy Williams <timothy.williams@nersc.no>
 * @date   Mon Aug  3 11:53:19 2015
 */

/**
 * Name conventions:
 * M_... = global variable
 * T_... = custom typename
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
    typedef T                                 T_val;
    typedef typename std::vector<T_val>       T_val_vec;
    typedef typename std::vector<T_val_vec*>  T_val_vec_ptrs;
    typedef typename std::vector<T_val_vec>   T_val_vec2d;//vector of vectors
    typedef typename std::vector<T_val_vec2d> T_val_vec3d;//vector of vectors of vectors
    
    //gmsh types
    typedef Nextsim::GmshMesh T_gmsh;

    //other types
    typedef Wim::IceParams<T_val>                         T_icep;
    typedef Wim::IceInfo<T_val>                           T_ice;
    typedef std::vector<T_ice>                            T_vec_ice;
    typedef boost::unordered_map<std::string,std::string> T_map;
    typedef boost::unordered_map<std::string,T_val_vec*>  T_map_vec_ptrs;

    enum IceType//distinguish between wim ice and sim ice
    {
        wim = 0,
        sim = 1
    };
    // ==========================================================================================

public:

    // ====================================================================================
    // public types
    typedef boost::unordered_map<std::string,T_val_vec> T_map_vec;

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
        T_val_vec nodes_x;                     // x-coords of nodes (not needed for structured grids)
        T_val_vec nodes_y;                     // y-coords of nodes (not needed for structured grids)
        T_val_vec elements_x;                  // x-coords of elements
        T_val_vec elements_y;                  // y-coords of elements
        T_val_vec surface;                     // surface area of elements
        std::vector<int> id;                        // id's of nodes
    } MeshInfo;

    typedef struct WimGrid
    {
        // information describing wim grid
        int nx;
        int ny;
        T_val dx;
        T_val dy;
        std::vector<T_val> X;
        std::vector<T_val> Y;
        std::vector<T_val> x;
        std::vector<T_val> y;
    } WimGrid;
    // ====================================================================================


public:

    // ====================================================================================
    //constructors
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
    // ====================================================================================


    // ====================================================================================
    // functions

    //make a grid
    void gridProcessing();
    void gridProcessing(T_gmsh const &mesh);
    void gridFromParameters();
    void gridPostProcessing();

    void saveGrid();
    void readGridFromFile();
    void readFromBinary(std::fstream &in, T_val_vec& in_array, int off = 0, std::ios_base::seekdir direction = std::ios::beg,
            int addx = 0, int addy = 0);
    void readDataFromFile(std::string const& filein);

    void exportResults(std::string const& output_type);
    void exportResultsGrid(T_map_vec_ptrs& extract_fields,
            std::vector<std::string> const& strings);
    void exportResultsMesh(T_map_vec_ptrs & extract_fields,
            std::vector<std::string> const &filenames,
            bool export_mesh=true, bool export_fields=true);
    void exportMesh(std::string const &filename);
    void testMesh();

    void saveLog(T_val const& t_out) const;
    void saveOptionsLog();

    void init1(int const& nextsim_cpt);
    void init(int const& nextsim_cpt=0);
    void init(T_gmsh const &mesh,int const& nextsim_cpt=0);//init WIM from gmsh mesh
    void init(T_gmsh const &mesh,BamgMesh* bamgmesh,int const& flag_fix,int const& nextsim_cpt=0);

    void init2();
    void initConstant(int const& nextsim_cpt);
    void assign();
    void assignSpatial();

    void update();

    void updateWaveMedium();

    void timeStep();

    T_val nfloesToDfloe(
                 T_val const& m_nfloes,
                 T_val const& m_conc);
    T_val_vec nfloesToDfloe(
                 std::vector<T_val> const& m_nfloes,
                 std::vector<T_val> const& m_conc);
    
    T_val dfloeToNfloes(
                 T_val const& m_dfloe,
                 T_val const& m_conc);
    T_val_vec dfloeToNfloes(
                 std::vector<T_val> const& m_dfloe,
                 std::vector<T_val> const& m_conc);

    void getFsdMesh(T_val_vec &nfloes_out,T_val_vec &dfloe_out,T_val_vec &broken);
    void getFsdMesh(T_val_vec &nfloes_out,T_val_vec &dfloe_out, T_val_vec &broken,
            T_val_vec const & conc_tot, T_gmsh const &mesh_in,T_val_vec const &um_in);
    void getFsdMesh(T_val_vec &nfloes_out,T_val_vec &dfloe_out, T_val_vec &broken,
            T_val_vec const & conc_tot, T_gmsh const &mesh_in);

    //time in seconds from init_time
    T_val getModelTime(int lcpt=0) const {return M_update_time+lcpt*M_timestep;}

    //time in days from ref time (1900-1-1)
    T_val getNextsimTime() const;

    // ========================================================================
    // breaking on mesh
    void setMesh( T_gmsh const &mesh);
    void setMesh( T_gmsh const &mesh,T_val_vec const &um);
    void setMesh( T_gmsh const &mesh,T_val_vec const &um,BamgMesh* bamgmesh,int const& flag_fix,bool const& regridding=false);
    void setMesh( T_gmsh const &mesh,BamgMesh* bamgmesh,int const& flag_fix,bool const& regridding=false);

    void resetMesh( T_gmsh const &mesh);
    void resetMesh( T_gmsh const &mesh,T_val_vec const &um);

    void setIceFields( std::vector<T_val> const& m_conc,  // conc
                       std::vector<T_val> const& m_vol, // ice vol or effective thickness (conc*thickness)
                       std::vector<T_val> const& m_nfloes,// Nfloes=conc/Dmax^2
                       bool const pre_regrid);

    T_val_vec getRelativeMeshDisplacement(T_gmsh const &mesh_in) const;
    T_val_vec getRelativeMeshDisplacement(T_gmsh const &mesh_in,T_val_vec const &um_in) const;
    void updateWaveSpec( T_gmsh const &mesh);
    void updateWaveSpec( T_gmsh const &mesh,T_val_vec const &um);

    void clearMeshFields() {M_ice[IceType::sim].clearFields();}
    void gridToPoints(
        T_val_vec_ptrs &output_data,
        T_val_vec_ptrs const &input_data,
        T_val_vec &Rx, T_val_vec &Ry);
    void meshToGrid(
        T_val_vec_ptrs &output_data,       //output data
        T_val_vec_ptrs const &input_data); //input data
    void meshToPoints(
        T_val_vec_ptrs &output_data,       //output data
        T_val_vec_ptrs const &input_data,  //input data
        T_val_vec &Rx,                     //location of output data (x-coord)
        T_val_vec &Ry);                    //location of output data (y-coord)
    void elementsToNodes(
        T_val_vec_ptrs &output_data,       //output data
        T_val_vec_ptrs const &input_data); //input data

    T_map_vec returnFieldsElements(std::vector<std::string> const&fields,
            T_val_vec &xel, T_val_vec &yel, T_val_vec const&surface_fac);
    T_map_vec returnFieldsElements(std::vector<std::string> const&fields,
            T_gmsh const &mesh_in,T_val_vec const &um_in);
    T_map_vec returnFieldsElements(std::vector<std::string> const&fields,
            T_gmsh const &mesh_in);

    T_map_vec returnFieldsNodes(std::vector<std::string> const&fields,
            T_val_vec &xnod, T_val_vec &ynod);
    T_map_vec returnFieldsNodes(std::vector<std::string> const&fields,
            T_gmsh const &mesh_in,T_val_vec const &um_in);
    T_map_vec returnFieldsNodes(std::vector<std::string> const&fields,
            T_gmsh const &mesh_in);

    T_val_vec getSurfaceFactor(T_gmsh const &mesh_in);

    void returnWaveStress(T_val_vec &M_tau, T_val_vec &xnod, T_val_vec &ynod);
    void returnWaveStress(T_val_vec &M_tau, T_gmsh const &mesh_in,T_val_vec const &um_in);
    void returnWaveStress(T_val_vec &M_tau, T_gmsh const &mesh_in);

    // ========================================================================

    WimGrid wimGrid(std::string const& units="m");

    void run();

    //===========================================================================
    //FSD: Dmax -> <D^moment> conversion
    void floeScaling(
          T_val const& dmax, int const& moment, T_val& dave);
    void floeScalingSmooth(
          T_val const& dmax, int const& moment, T_val& dave);
    //===========================================================================

    //===========================================================================
    //advection/attenuation
    void advectDirections( T_val_vec2d& Sdir, T_val_vec const& ag2d_eff);
    void advectDirectionsMesh( T_val_vec2d& Sdir, T_val_vec& ag2d_eff,
            T_val_vec const& boundary_vals);
    void attenSimple(
            T_val_vec2d& Sdir, T_val_vec& Sfreq,
            T_val_vec& taux_omega,T_val_vec& tauy_omega,
            T_val_vec& sdx_omega,T_val_vec& sdy_omega,
            T_val_vec const& ag2d_eff);
    void attenIsotropic(
            T_val_vec2d& Sdir, T_val_vec& Sfreq,
            T_val_vec& taux_omega,T_val_vec& tauy_omega,
            T_val_vec& sdx_omega,T_val_vec& sdy_omega,
            T_val_vec const& ag2d_eff);
    void waveAdvWeno(
            T_val_vec& h, T_val_vec const& u, T_val_vec const& v);
    void weno3pdV2(
            T_val_vec const& gin, T_val_vec const& u, T_val_vec const& v,
            T_val_vec const& scuy, T_val_vec const& scvx,
            T_val_vec const& scp2i, T_val_vec const& scp2,
            T_val_vec& saoout);
    void padVar(T_val_vec const& u, T_val_vec& upad,
            std::string const& advopt_,bool const&steady=false);
    //===========================================================================

    void intWaveSpec();
    void intDirns(T_val_vec2d const& Sdir, T_val_vec& Sfreq,
            T_val_vec& sdx_omega, T_val_vec& sdy_omega);

    void idealWaveFields(T_val const xfac);
    void idealIceFields (T_val const xfac);
    void inputWaveFields(T_val_vec const& swh_in,
            T_val_vec const& mwp_in,
            T_val_vec const& mwd_in);
    void setIncWaveSpec(T_val_vec const& wave_mask);

    T_val thetaDirFrac(T_val const& th1_, T_val const& dtheta_, T_val const& mwd_);
    T_val thetaInRange(T_val const& th_, T_val const& th1, bool const& close_on_right=false);

    T_val_vec getX() const { return X_array; }
    T_val_vec getY() const { return Y_array; }
    T_val_vec getSCUY() const { return SCUY_array; }
    T_val_vec getSCVX() const { return SCVX_array; }
    T_val_vec getSCP2() const { return SCP2_array; }
    T_val_vec getSCP2I() const { return SCP2I_array; }
    T_val_vec getLANDMASK() const { return LANDMASK_array; }

    //for use at regridding time (M_wim_on_mesh)
    T_val_vec getMeshDisplacement() const { return M_UM; }
    T_val_vec3d getWaveSpec() const { return M_sdf_dir; }
    void setRelativeMeshDisplacement(T_val_vec const&um_in) { M_UM = um_in; return; }
    void setWaveSpec(T_val_vec3d const&sdf_in)
    {
        // reset wave spectrum after regrid;
        // also integrate if now so that initial wave fields are correct
        // and so they are ready for export;
        M_sdf_dir = sdf_in;
        this->intWaveSpec();
        return;
    }

    std::string getWimGridFilename() const { return wim_gridfile; }


private:

    po::variables_map vm;
    int nx, ny, nxext, nyext, nbdy, nbdx, nghost;
    int num_p_wim,num_q_wim,num_u_wim,num_v_wim;
    int M_itest;
    T_val_vec X_array, Y_array, SCUY_array, SCVX_array,
                SCP2_array, SCP2I_array, LANDMASK_array;
    T_val_vec x_col,y_row;

    T_val M_cfl, M_length_cfl, M_max_cg;
    T_val M_current_time;
    T_val Tmin, Tmax;
    T_val xmax, ym, x0, y0, dx, dy, x_edge;
    T_val M_dfloe_pack_init;
    T_val M_timestep = 0.;
    T_val M_duration = 0.;
    int M_num_timesteps;

    int nwavedirn, nwavefreq, M_advdim;
    bool M_ref_Hs_ice, M_atten, M_useicevel, M_steady, M_breaking;
    bool M_dump_diag;
    bool docoupling;
    std::string M_scatmod, M_advopt;// fsdopt;
    std::string wim_gridfile;

    //dimension of wavedir
    T_val_vec M_wavedir, M_quadrature_wt_dir;
    
    //dimension of freq
    T_val_vec M_quadrature_wt_freq, M_freq_vec;
    T_val_vec M_wlng_wtr, M_ag_wtr, M_ap_wtr;

    //dimension of space
    T_val_vec M_Hs,M_Tp,M_mwd;
    T_val_vec M_swh_in,M_mwp_in,M_mwd_in;
    T_val_vec M_dave;

    //depend on freq and position
    T_val_vec2d M_ag_eff, M_agnod_eff, M_ap_eff, M_wlng_ice, M_atten_nond, M_damping, M_disp_ratio;

    //depend on freq, dirn and position
    T_val_vec2d M_open_boundary_vals;
    T_val_vec3d M_sdf_dir;
    T_val_vec3d M_sdf_dir_inc;//
    T_val_vec M_steady_mask;//TODO for grid advection, make steady a boundary condition

    //these are only temporary vectors, but they are global in order to
    //save creating and destroying them extremely often
    T_val_vec Mtmp_sdf_freq;                                   //for mom0 integral
    T_val_vec Mtmp_taux_om, Mtmp_tauy_om;                      //for taux,tau_y integrals
    T_val_vec Mtmp_stokes_drift_x_om, Mtmp_stokes_drift_y_om;  //for mwd_x,mwd_y, stokes_drift_x,stokes_drift_y integrals
    T_val_vec Mtmp_mom0,Mtmp_mom2,Mtmp_var_strain;
    T_val_vec Mtmp_mom0w,Mtmp_mom2w;
    T_val_vec Mtmp_atten_dim, Mtmp_damp_dim;

    T_val_vec M_mwd_x, M_mwd_y, M_tau_x, M_tau_y,M_stokes_drift_x,M_stokes_drift_y;//row-major order (C)

    int M_max_threads;
    boost::mpi::timer chrono;

    std::string M_init_time_str;
    T_val M_restart_time     = 0.;// time of restarting, relative to init_time
    T_val M_update_time      = 0; // time of start of call to wim.run(), relative to init_time
    T_val M_time_mesh_set    = 0;// time of start of call to wim.run(), relative to init_time
    int M_cpt                = 0;//global counter
    int M_num_elements       = 0;
    bool M_regular           = false;
    bool M_initialised_ice   = false;
    bool M_initialised_waves = false;

    int M_nb_export_nextwim = 0;
    int M_nb_export_inc     = 0;
    int M_nb_export_init    = 0;
    int M_nb_export_prog    = 0;
    int M_nb_export_final   = 0;
    int M_nb_mesh_test      = 0;

    MeshInfo M_wim_triangulation,nextsim_mesh,nextsim_mesh_old;
    std::vector<int> wet_indices;

    T_icep M_ice_params;
    //T_map_ice M_ice;
    T_vec_ice M_ice;

    bool M_break_on_mesh = false;// do breaking on nextsim mesh as well as on grid
    bool M_wim_on_mesh   = false;// to run WIM on nextsim mesh
    
    T_val_vec M_UM;// displacement of mesh nodes between calls to wim.run()
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
