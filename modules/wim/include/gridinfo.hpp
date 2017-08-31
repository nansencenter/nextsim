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

#ifndef __GRIDINFO_H
#define __GRIDINFO_H 1

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
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
//#include <Bamgx.h>
//#include <iomanip>
#include <omp.h>
//#include <gmshmesh.hpp>
#include <meshtools.hpp>

#ifdef PI
#undef PI
#define PI M_PI
#endif

namespace Wim
{

namespace po = boost::program_options;
namespace fs = boost::filesystem;

template<typename T=float> class GridInfo
{
    // ==========================================================================================
    typedef T                                 T_val;
    typedef typename std::vector<T_val>       T_val_vec;
    typedef typename std::vector<T_val_vec*>  T_val_vec_ptrs;
    typedef typename std::vector<T_val_vec>   T_val_vec2d;//vector of vectors
    typedef typename std::vector<T_val_vec2d> T_val_vec3d;//vector of vectors of vectors
    
    //gmsh types
    typedef Nextsim::GmshMesh T_gmsh;
    typedef Wim::MeshInfo<T_val> T_mesh;

    //other types
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
    //constructors
    ~GridInfo() {}
    GridInfo() {}
    GridInfo(po::variables_map const& vmIn);
    GridInfo(po::variables_map const& vmIn,T_mesh & mesh_in);
    // ====================================================================================


    // ====================================================================================
    // functions

    //make a grid
    void gridFromParameters();
    void gridPostProcessing();
    void setupAdvection();

    void saveGrid();
    void readGridFromFile();
    void readFromBinary(std::fstream &in, T_val_vec& in_array, int off = 0, std::ios_base::seekdir direction = std::ios::beg,
            int addx = 0, int addy = 0);

    void interpToPoints(
        T_val_vec_ptrs &output_data,
        T_val_vec_ptrs const &input_data,
        T_val_vec &Rx, T_val_vec &Ry);


    //===========================================================================
    //advection/attenuation
    void waveAdvWeno(
            T_val_vec& h, T_val_vec const& u, T_val_vec const& v, T_val const& timestep);
    void weno3pdV2(
            T_val_vec const& gin, T_val_vec const& u, T_val_vec const& v,
            T_val_vec const& scuy, T_val_vec const& scvx,
            T_val_vec const& scp2i, T_val_vec const& scp2,
            T_val_vec& saoout, T_val const& timestep);
    void padVar(T_val_vec const& u, T_val_vec& upad,
            std::string const& advopt_,bool const&steady=false);
    //===========================================================================

    T_val_vec getX() const { return M_px; }
    T_val_vec getY() const { return M_py; }
    bool isRegular() const { return (M_use_regular&&M_regular); }
    std::string getWimGridFilename() const { return M_gridfile; }


    // public variables
    // - scalars
    int M_num_px, M_num_py;
    T_val M_xmax, M_xmin, M_ymin, M_ymax, M_dx, M_dy;
    int M_num_elements;
    T_val M_resolution;
    bool M_initialised = false;
    std::string M_gridfile = "";

    // - arrays
    T_val_vec M_px, M_py, M_scuy, M_scvx,
                M_scp2, M_scp2i, M_land_mask;

private:

    po::variables_map vm;
    int M_num_px_ext, M_num_py_ext, M_nbdy_y, M_nbdy_x;
    int M_num_p, M_num_q, M_num_u, M_num_v;

    T_val_vec M_px_vec,M_py_vec;


    int M_max_threads;
    boost::mpi::timer chrono;

    std::string M_advopt;
    int M_advdim;
    bool M_use_regular = false;//use interp routines for regular grid if possible
    bool M_regular     = false;//is the grid regular

    T_mesh M_wim_triangulation;
    std::vector<int> M_wet_indices;

    //tmp variables in advection
    T_val_vec Mtmp_sao,Mtmp_hp,Mtmp_hp_tmp;
    T_val_vec Mtmp_u_pad, Mtmp_v_pad, Mtmp_scp2_pad, Mtmp_scp2i_pad, Mtmp_scuy_pad, Mtmp_scvx_pad, Mtmp_h_pad;
    T_val_vec Mtmp_ful, Mtmp_fuh, Mtmp_fvl, Mtmp_fvh, Mtmp_gt;
};

} // namespace Wim

#endif
