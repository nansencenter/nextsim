/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   tools.hpp
 * @author Timothy Williams <timothy.williams@nersc.no>
 * @date   Mon Jul 3 2017
 * created to move some of the finiteelement functions
 * to shorten file and also make them more accessible 
 * (eg so they can be called from the WIM)
 */

#ifndef __MESHTOOLS_H
#define __MESHTOOLS_H 1

#include <Bamgx.h>
#include <MemOps.h>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromMeshToGridx.h>
#include <gmshmesh.hpp>
#include <gmshmeshseq.hpp>
#include <omp.h>
//
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

//extern "C"
//{
//#include <mapx.h>
//}
//
//#include <wimdiscr.hpp>

namespace Wim
{

template<typename T=float,typename FEMeshType=Nextsim::GmshMesh> class MeshInfo
{

    typedef T T_val;
    typedef std::vector<T_val> T_val_vec;
    typedef typename std::vector<T_val_vec*>  T_val_vec_ptrs;
    //typedef Nextsim::GmshMesh T_gmsh;
    typedef FEMeshType T_gmsh;

    //private variables    
    int  M_max_node_el_conn = 0;
    T_val_vec M_element_connectivity;
    T_val_vec M_node_element_connectivity;
    T_val_vec M_mask_dirichlet;

public:
    int M_num_nodes    = 0;
    int M_num_elements = 0;
    T_val M_resolution;
    std::vector<int> M_index;
    std::vector<int> M_id;
    T_val_vec M_nodes_x;
    T_val_vec M_nodes_y;
    T_val_vec M_elements_x;
    T_val_vec M_elements_y;
    T_val_vec M_surface;

    enum class E_mesh_type {
        uninitialised = 0,//not initialised
        simple        = 1,//simple
        full          = 2//full info - eg can use advect
    };
    E_mesh_type M_mesh_type = E_mesh_type::uninitialised;

    // ===============================================================
    //destructor
    ~MeshInfo() {}

    // constructors
    // - default
    MeshInfo() {}

    // - nodes_only
    MeshInfo(std::vector<int> const &index,T_val_vec const &nodes_x,T_val_vec const &nodes_y);
    MeshInfo(T_val_vec const &nodes_x,T_val_vec const &nodes_y);

    // - simple
    MeshInfo(std::vector<int> const &index,T_val_vec const &nodes_x,T_val_vec const &nodes_y,
       T_val_vec const &elements_x,T_val_vec const &elements_y);
    MeshInfo( T_gmsh const &mesh);

    // - full
    MeshInfo( T_gmsh const &mesh,BamgMesh* bamgmesh,int const& flag_fix);
    // ===============================================================


    // ===============================================================
    //functions
    bool initialised() const
    { return ( !(M_mesh_type == E_mesh_type::uninitialised) ); }

    void initSimple();
    void setElements();
    void setFields( T_gmsh const &mesh);
    T_val lengthCfl() const { return M_resolution/12.; }
    void advect(T_val** interp_elt_out_ptr, // pointer to pointer to output data
            T_val* interp_elt_in,           // pointer to input data
            T_val* VC_in,                   // pointer to convective velocities (len = 2*num_nodes)
            int* interp_method,             // pointer to interp methods for each variable
            int nb_var,                     // number of variables
            T_val time_step,                // time step (s)
            T_val* inc_values);             // values incoming from open boundaries
                                            // - currently just specify a constant (ie spatially invariant) value for each variable
    void interpToPoints(
            T_val_vec_ptrs &output_data,            //output data
            T_val_vec_ptrs const &input_data,       //input data
            T_val_vec &Rx,                          //location of output data (x-coord)
            T_val_vec &Ry,                          //location of output data (y-coord)
            std::vector<int> const& wet_indices={}, //optional: say which cells of input_data are wet
            bool  const &use_default=false,         //optional: use a default value for points outside the mesh (else use the nearest element)
            T_val const &default_val=0.);           //optional: default value for points outside the mesh
    void elementsToNodes(
            T_val_vec_ptrs &output_data,       //output data
            T_val_vec_ptrs const &input_data); //input data

    void interpToGrid(
            T_val_vec_ptrs &output_data,            //output data
            T_val_vec_ptrs const &input_data,       //input data
            T_val const &xmin, T_val const &ymax,   //grid parameters
            int const&nx, int const&ny,
            T_val const &dx, T_val const &dy,
            T_val const &default_value=0.);

    T_val_vec getX() const { return M_elements_x; }
    T_val_vec getY() const { return M_elements_y; }

    T_val jacobian(T_val const x0, T_val const y0,
            T_val const x1, T_val const y1,
            T_val const x2, T_val const y2) const
    {
        //signed area of a triangle with vertices (going anti-clockwise)
        //(x0,y0), (x1,y1), (x2,y2)
        return (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
    }


};//class MeshInfo

} // namespace Wim
#endif
