/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   tools.hpp
 * @author Timothy Williams <timothy.williams@nersc.no>
 * @date   Mon Jul 3 2017
 * created to move some of the finiteelement functions
 * to shorten file and also make them more accessible 
 * (eg so they can be called from the WIM)
 */


//#include <solverpetsc.hpp>
//#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/version.hpp>
#include <boost/format.hpp>
//#include <BamgConvertMeshx.h>
#include <Bamgx.h>
//#include <InterpFromMeshToMesh2dCavities.h>
//#include <InterpFromMeshToMesh2dx.h>
//#include <InterpFromGridToMeshx.h>
//#include <InterpFromMeshToGridx.h>
#include <gmshmesh.hpp>
//#include <graphcsr.hpp>
//#include <debug.hpp>
#include <omp.h>
//#include <externaldata.hpp>
//#include <gridoutput.hpp>
//#include <dataset.hpp>
//#include <drifters.hpp>
//
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

extern "C"
{
#include <mapx.h>
}

#include <wimdiscr.hpp>

namespace NextsimTools
{

    // ==========================================================================================
    //typenames
    typedef typename Nextsim::GmshMesh::point_type point_type;
    typedef typename Nextsim::GmshMesh::element_type element_type;
    typedef Nextsim::GmshMesh mesh_type;
    typedef typename Wim::WimDiscr<double>::MeshInfo mesh_info_type_dbl;
    // ==========================================================================================


    // ==========================================================================================
    // mesh functions
    double jacobian(element_type const& element, mesh_type const& mesh);
    double jacobian(element_type const& element, mesh_type const& mesh,
                    std::vector<double> const& um, double factor = 1.);
    std::vector<double> sides(element_type const& element, mesh_type const& mesh);
    std::vector<double> minMaxSide(mesh_type const& mesh);
    //void movedMesh(std::vector<double> const& um, double factor = 0);
    double measure(element_type const& element, mesh_type const& mesh);
    double measure(element_type const& element, mesh_type const& mesh,
                   std::vector<double> const& um, double factor = 1.);

    double minAngles(element_type const& element, mesh_type const& mesh);
    double minAngle(mesh_type const& mesh);

    double minAngle(mesh_type const& mesh, std::vector<double> const& um, double factor);

    bool flip(mesh_type const& mesh, std::vector<double> const& um, double factor);

    double resolution(mesh_type const& mesh);

    std::vector<double> hminVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh);
    std::vector<double> hmaxVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh);

    std::vector<double> AllMinAngle(mesh_type const& mesh, std::vector<double> const& um, double factor);

    //void nodesToElements(double const* depth, std::vector<double>& v);
    // ==========================================================================================

#if 1
    void advect(double** interp_elt_out_ptr, // pointer to pointer to output data
        double* interp_elt_in,      // pointer to input data
        mesh_info_type_dbl* mesh_info,        // pointer to structure with mesh info: positions of nodes and elements,
                                    //  index (maps elements to nodes), element connectivity
        double* VC_in,              // pointer to convective velocities (len = 2*num_nodes)
        int* interp_method,         // pointer to interp methods for each variable
        int nb_var,                 // number of variables
        double time_step);          // time step (s)
#endif
    //void diffuse

    // ==========================================================================================
    // other functions
    std::string gitRevision();
    std::string system(std::string const& command);
    std::string getEnv(std::string const& envname);
    // ==========================================================================================

} // namespace NextsimTools
