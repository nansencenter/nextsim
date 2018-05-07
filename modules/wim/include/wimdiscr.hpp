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
//#include <InterpFromGridToMeshx.h>
//#include <InterpFromMeshToMesh2dx.h>
//#include <InterpFromMeshToGridx.h>
#include <Bamgx.h>
#include <iomanip>
#include <omp.h>
#include <gmshmesh.hpp>
#include <iceinfo.hpp>
#include <meshtools.hpp>
#include <gridinfo.hpp>

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
    typedef Wim::MeshInfo<T_val>                          T_mesh;
    typedef Wim::GridInfo<T_val>                          T_grid;
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
    // ====================================================================================


public:

    // ====================================================================================
    //constructors
    ~WimDiscr() {}
    WimDiscr() {}
    WimDiscr(po::variables_map const& vmIn);
    // ====================================================================================


    // ====================================================================================
    // init
    void initStandAlone();
    void initCoupled(int const& nextsim_cpt,
            T_gmsh const &movedmesh,
            BamgMesh* bamgmesh,
            int const& flag_fix);
    void initConstant(int const& nextsim_cpt=0);
    void initRemaining();
    void assign();
    void assignSpatial();

    // update
    void update();
    void updateWaveMedium();

    // main routines
    void run();
    void timeStep();
    // ====================================================================================


    // ====================================================================================
    // export
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
    // ====================================================================================


    // ========================================================================
    // set mesh

    // wim on its own grid
    // - need mesh for interpolation etc
    void setMeshSimple( T_gmsh const &mesh);

    // wim on nextsim mesh
    void setMeshFull( T_gmsh const &mesh,
            BamgMesh* bamgmesh,
            int const& flag_fix,
            bool const& regridding=false);

    // wrapper for setMeshSimple and setMeshFull
    void setMesh( T_gmsh const &mesh,
            BamgMesh* bamgmesh,
            int const& flag_fix,
            bool const& regridding=false);

    T_val_vec getRelativeMeshDisplacement(T_gmsh const &movedmesh) const;
    T_val_vec getMeshDisplacement() const { return M_UM; }
    void setRelativeMeshDisplacement(T_val_vec const&um_in) { M_UM = um_in; return; }

    void updateWaveSpec( T_gmsh const &mesh);
    T_val_vec getSurfaceFactor(T_gmsh const &movedmesh);

    T_val_vec3d getWaveSpec() const { return M_sdf_dir; }
    void setWaveSpec(T_val_vec3d const&sdf_in)
    {
        // reset wave spectrum after regrid;
        // also integrate if now so that initial wave fields are correct
        // and so they are ready for export;
        M_sdf_dir = sdf_in;
        this->intWaveSpec();
        return;
    }
    // ==============================================================================


    // ==============================================================================
    // set ice
    void idealIceFields (T_val const xfac=0.7);
    void setIceFields( std::vector<T_val> const& conc,  // conc
                       std::vector<T_val> const& vol, // ice vol or effective thickness (conc*thickness)
                       std::vector<T_val> const& nfloes);// Nfloes=conc/Dmax^2
                       //bool const pre_regrid);
    void clearMeshFields() { M_ice[IceType::sim].clearFields(); }

    // set waves
    void idealWaveFields(T_val const xfac=0.8);
    void setWaveFields(T_val_vec const& swh_in,
            T_val_vec const& mwp_in,
            T_val_vec const& mwd_in);
    void setIncWaveSpec(T_val_vec const& wave_mask);
    // ==============================================================================


    // ==============================================================================
    // output fields, wave stress,...
    T_map_vec returnFieldsElements(std::vector<std::string> const&fields,
            T_val_vec &xel, T_val_vec &yel, T_val_vec const&surface_fac);
    T_map_vec returnFieldsElements(std::vector<std::string> const&fields,
            T_gmsh const &movedmesh);

    T_map_vec returnFieldsNodes(std::vector<std::string> const&fields,
            T_val_vec &xnod, T_val_vec &ynod);
    T_map_vec returnFieldsNodes(std::vector<std::string> const&fields,
            T_gmsh const &movedmesh);

    T_val_vec returnWaveStress(T_val_vec &xnod, T_val_vec &ynod);
        //base interface - pass in the nodes for interpolation to these points
    T_val_vec returnWaveStress(T_gmsh const &movedmesh);
        //pass in the moved mesh, then get the nodes
    T_val_vec returnWaveStress();
        //if running on mesh, have stresses on elements: interp them to nodes

    T_val_vec combineVectorComponents(T_val_vec const &vec_x, T_val_vec const &vec_y);
        //combine x,y components of vector into one vector:
        // {vec_x[0], ..., vec_x[N-1], vec_y[0], ..., vec_y[N-1]}
    // ========================================================================


    // ========================================================================
    // fsd stuff
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
            T_val_vec const & conc_tot, T_gmsh const &movedmesh);
    // ========================================================================


    // ==========================================================================
    // advection/attenuation
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

    // integrate the wave spectrum
    void intWaveSpec();
    void intDirns(T_val_vec2d const& Sdir, T_val_vec& Sfreq,
            T_val_vec& sdx_omega, T_val_vec& sdy_omega);
    // ==========================================================================


    // ==========================================================================
    // other functions

    // interpolation
    void interpIceMeshToGrid();

    // time in seconds from init_time
    T_val getModelTime(int lcpt=0) const { return M_update_time+lcpt*M_timestep; }

    // time in days from ref time (1900-1-1)
    T_val getNextsimTime() const;

    T_val thetaDirFrac(T_val const& th1_, T_val const& dtheta_, T_val const& mwd_) const;
    T_val thetaInRange(T_val const& th_, T_val const& th1, bool const& close_on_right=false) const;

    T_val_vec getX() const
    { 
        if(!M_wim_on_mesh)
            return M_grid.getX();
        else
            return M_mesh.getX();
    }
    T_val_vec getY() const
    { 
        if(!M_wim_on_mesh)
            return M_grid.getY();
        else
            return M_mesh.getY();
    }
    int getNumElements()
    { 
        if(M_wim_on_mesh)
            return M_mesh.M_num_elements;
        else
            return M_grid.M_num_elements;
    }
    void getRangeXY(T_val &xmin,T_val &xmax, T_val &ymin, T_val &ymax) const;
    void printRange(std::string const &name, T_val_vec const &vec, int const & prec=0) const;
    void getRange(T_val_vec const &vec, T_val &xmin, T_val &xmax) const;

    std::string getWimGridFilename() const { return M_grid.M_gridfile; }
    // ==========================================================================


private:

    po::variables_map vm;
    int M_itest;

    T_val M_cfl, M_length_cfl, M_max_cg;
    T_val Tmin, Tmax;
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
    T_val M_current_time     = 0.;// current model time (s), relative to init_time (updated every timestep)
    T_val M_restart_time     = 0.;// time (s) of restarting, relative to init_time
    T_val M_update_time      = 0; // time (s) of start of call to wim.run(), relative to init_time
    T_val M_time_mesh_set    = 0; // time (s) of start of call to wim.run(), relative to init_time
    int M_cpt                = 0; // global counter
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

    // mesh, grid objects
    T_mesh M_mesh,M_mesh_old;
    T_grid M_grid;

    T_icep M_ice_params;// ice parameters
    T_vec_ice M_ice;    // vector of ice objects

    bool M_break_on_mesh = false;// do breaking on nextsim mesh as well as on grid
    bool M_wim_on_mesh   = false;// to run WIM on nextsim mesh
    bool M_do_coupling   = false;// run WIM coupled to nextsim
    
    T_val_vec M_land_mask;
    T_val_vec M_UM;// displacement of mesh nodes between calls to wim.run()
                        // - only used if running WIM on nextsim mesh
                        // - for correction to group velocity at advection time
    // =========================================================================
};

} // namespace Wim

#endif
