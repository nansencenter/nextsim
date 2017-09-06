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

template<typename T=float> class MeshInfo
{

    typedef T T_val;
    typedef std::vector<T_val> T_val_vec;
    typedef typename std::vector<T_val_vec*>  T_val_vec_ptrs;
    typedef Nextsim::GmshMesh T_gmsh;

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
    void initSimple();
    void setElements();
    void setFields( T_gmsh const &mesh);
    T_val lengthCfl() const { return 0.1*M_resolution; }
    void advect(T_val** interp_elt_out_ptr, // pointer to pointer to output data
            T_val* interp_elt_in,           // pointer to input data
            T_val* VC_in,                   // pointer to convective velocities (len = 2*num_nodes)
            int* interp_method,             // pointer to interp methods for each variable
            int nb_var,                     // number of variables
            T_val time_step,                // time step (s)
            T_val* inc_values);             // values incoming from open boundaries
                                            // - currently just specify a constant (ie spatially invariant) value for each variable
    void interpToPoints(
            T_val_vec_ptrs &output_data,          //output data
            T_val_vec_ptrs const &input_data,     //input data
            T_val_vec &Rx,                        //location of output data (x-coord)
            T_val_vec &Ry,                        //location of output data (y-coord)
            std::vector<int> const& wet_indices={}); //optional: say which cells of input_data are wet

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
    // ===============================================================

};//class MeshInfo


namespace MeshTools
{

// ==========================================================================================
//typenames
typedef typename Nextsim::GmshMesh::element_type T_gmsh_el;
typedef Nextsim::GmshMesh T_gmsh;
//typedef typename Wim::WimDiscr<double>::MeshInfo mesh_info_type_dbl;
// ==========================================================================================


// ==========================================================================================
// mesh functions
double jacobian(double const x0, double const y0, double const x1, double const y1,double const x2, double const y2);
double jacobian(T_gmsh_el const& element, T_gmsh const& mesh);
double jacobian(T_gmsh_el const& element, T_gmsh const& mesh,
                std::vector<double> const& um, double factor = 1.);
std::vector<double> sides(T_gmsh_el const& element, T_gmsh const& mesh);
std::vector<double> minMaxSide(T_gmsh const& mesh);
//void movedMesh(std::vector<double> const& um, double factor = 0);
double measure(double const x0, double const y0, double const x1, double const y1,double const x2, double const y2);
double measure(T_gmsh_el const& element, T_gmsh const& mesh);
double measure(T_gmsh_el const& element, T_gmsh const& mesh,
               std::vector<double> const& um, double factor = 1.);

double minAngles(T_gmsh_el const& element, T_gmsh const& mesh);
double minAngle(T_gmsh const& mesh);

double minAngle(T_gmsh const& mesh, std::vector<double> const& um, double factor);

bool flip(T_gmsh const& mesh, std::vector<double> const& um, double factor);

double resolution(T_gmsh const& mesh);

std::vector<double> hminVertices(T_gmsh const& mesh, BamgMesh const* bamg_mesh);
std::vector<double> hmaxVertices(T_gmsh const& mesh, BamgMesh const* bamg_mesh);

std::vector<double> AllMinAngle(T_gmsh const& mesh, std::vector<double> const& um, double factor);

//void nodesToElements(double const* depth, std::vector<double>& v);
// ==========================================================================================

#if 0
void advect(double** interp_elt_out_ptr, // pointer to pointer to output data
double* interp_elt_in,               // pointer to input data
mesh_info_type_dbl* mesh_info,       // pointer to structure with mesh info: positions of nodes and elements,
                                     //  index (maps elements to nodes), element connectivity
double* VC_in,                       // pointer to convective velocities (len = 2*num_nodes)
int* interp_method,                  // pointer to interp methods for each variable
int nb_var,                          // number of variables
double time_step,                    // time step (s)
double* inc_values);                 // pointer to values coming in from open boundaries (len = nb_var)
#endif
//void diffuse

// ==========================================================================================
// other functions
std::string gitRevision();
std::string system(std::string const& command);
std::string getEnv(std::string const& envname);
// ==========================================================================================

} // namespace MeshTools

} // namespace Wim
#endif
