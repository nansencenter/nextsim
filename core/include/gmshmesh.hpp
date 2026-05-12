/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gmshmesh.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Wed Jul 29 12:51:51 2015
 */

#ifndef __GmshMesh_HPP
#define __GmshMesh_HPP 1

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <assert.hpp>
#include <boost/bimap.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/mpi/collectives.hpp>

#include <environment.hpp>
#include <entities.hpp>
//#include <MElement.h>
//#include <GModel.h>
#include "debug.hpp"

extern "C"
{
#include <mapx.h>
}

// from Gmsh
void SwapBytes(char *array, int size, int n);

namespace Nextsim
{

class GmshMesh
{
public:

    typedef Nextsim::entities::GMSHPoint point_type;
    typedef Nextsim::entities::GMSHElement element_type;

    typedef boost::bimap<int,int> bimap_type;
    typedef bimap_type::value_type position;


    GmshMesh( Communicator const& comm = Environment::comm() );

    GmshMesh( GmshMesh const& mesh );

    // GmshMesh(std::vector<point_type> const& nodes,
    //          std::vector<element_type> const& edges,
    //          std::vector<element_type> const& triangles,
    //          Communicator const& comm = Environment::comm());

    void readFromFile(std::string const& filename, std::string const& format="ascii");
    void readFromFileBinary(std::ifstream& ifs);
    void readFromFileASCII(std::ifstream& ifs);
    void writeToFile(std::string const& filename);
    void move(std::vector<double> const& um, double factor);
    void update(std::vector<point_type> const& nodes, std::vector<element_type> const& triangles, int numElements);
    void allGather(std::vector<int> const& field_in, std::vector<std::vector<int> >& field_out, int& acc_size);
    void nodalGrid();

    Communicator const& comm() const { return M_comm; }
    std::string const& version() const {return M_version;}
    std::string const& ordering() const {return M_ordering;}
    std::string const& mppfile() const {return M_mppfile;}
    std::map<int, point_type > const& nodes() const {return M_nodes;}
    std::vector<element_type> const& triangles() const {return M_triangles;}
    std::vector<element_type> const& edges() const {return M_edges;}

    std::map<std::string, std::vector<int> > markerNames() const {return M_marker_names;}

    int numGlobalNodes() const {return M_global_num_nodes;}
    int numGlobalNodesFromSarialMesh() const {return M_global_num_nodes_from_serial;}
    int numNodes() const {return M_num_nodes;}
    int numTriangles() const {return M_num_triangles;}
    int numEdges() const {return M_num_edges;}

    int numLocalNodesWithoutGhost() const {return M_nldof_without_ghost;}
    int numLocalNodesWithGhost() const {return M_nldof_with_ghost;}
    int numLocalGhost() const {return M_nlghost;}

    int numGlobalElements() const {return M_global_num_elements;}
    int numGlobalElementsFromSarialMesh() const {return M_global_num_elements_from_serial;}
    int numTrianglesWithoutGhost() const {return M_num_triangles_without_ghost;}

    void setCommunicator(Communicator const& comm) {M_comm=comm;}
    void setOrdering(std::string const& order) {M_ordering=order;}

    void setNodes(std::map<int, point_type > const& nodes) {M_nodes=nodes;}
    void setEdges(std::vector<element_type> const& edges) {M_edges=edges;}
    void setTriangles(std::vector<element_type> const& triangles) {M_triangles=triangles;}

    void setNumNodes(int const& nnodes) {M_num_nodes=nnodes;}
    void setNumEdges(int const& nlns) {M_num_edges=nlns;}
    void setNumTriangles(int const& ntrs) {M_num_triangles=ntrs;}

    void stereographicProjection();

    std::vector<int> const& localDofWithoutGhost() const {return M_local_dof_without_ghost;}
    std::vector<int> const& localDofWithGhost() const {return M_local_dof_with_ghost;}
    std::vector<int> const& localGhost() const {return M_local_ghost;}

    std::vector<int> const& localDofWithGhostInit() const {return M_local_dof_with_ghost_init;}
    std::vector<int> const& trianglesIdWithGhost() const {return M_triangles_id_with_ghost;}

    bimap_type const& transferMapInit() const {return M_transfer_map;}
    bimap_type const& transferMap() const {return M_transfer_map_reordered;}

    bimap_type const& transferMapElt() const {return M_transfer_map_elt;}
    //bimap_type const& transferMapReordered() const {return M_transfer_map_reordered;}

    //bimap_type const& mapNodes() const {return M_reorder_map_nodes;}
    //std::map<int,int> const& mapNodes() const {return M_reorder_map_nodes;}
    //bimap_type const& mapElements() const {return M_reorder_map_elements;}
    //std::map<int,int> const& mapElements() const {return M_reorder_map_elements;}

    std::vector<int> const& mapNodes() const {return M_map_nodes;}
    std::vector<int> const& mapElements() const {return M_map_elements;}

    std::vector<int> indexTr() const;

    std::vector<double> coordX() const;
    std::vector<double> coordY() const;

    std::vector<double> coordX(double const& rotangle) const;
    std::vector<double> coordY(double const& rotangle) const;

    std::vector<double> bCoordX() const;
    std::vector<double> bCoordY() const;

    std::vector<double> bCoordX(double const& rotangle) const;
    std::vector<double> bCoordY(double const& rotangle) const;

    std::vector<std::vector<double>>
        vertices(std::vector<int> const& indices) const;
    std::vector<std::vector<double>>
        vertices(std::vector<int> const& indices,
                std::vector<double> const& um, double factor) const;

    std::vector<double> meanLat() const;
    std::vector<double> meanLon() const;

    std::vector<double> lat() const;
    std::vector<double> lon() const;

    std::vector<int> indexTrPartition() const;
    std::vector<double> coordXPartition() const;
    std::vector<double> coordYPartition() const;

    void setId(std::vector<int> const& newid);
    std::vector<int> id() const;

private:

    Communicator M_comm;
    std::string M_version;
    std::string M_ordering;
    std::string M_mppfile;
    LogLevel M_log_level;
    bool M_log_all;

    std::vector<point_type> M_nodes_vec;
    std::map<int, point_type > M_nodes;
    std::vector<element_type> M_triangles;
    std::vector<element_type> M_edges;

    std::vector<int> M_local_dof_with_ghost;
    std::vector<int> M_local_dof_with_ghost_init;
    std::vector<int> M_local_dof_without_ghost;
    std::vector<int> M_local_ghost;
    std::vector<int> M_triangles_id_with_ghost;

    int M_global_num_nodes;
    int M_global_num_nodes_from_serial;
    int M_num_nodes;
    int M_num_triangles;
    int M_num_edges;

    int M_nldof_with_ghost;
    int M_nldof_without_ghost;
    int M_nlghost;

    int M_global_num_elements;
    int M_global_num_elements_from_serial;
    int M_num_triangles_without_ghost;


    bimap_type M_transfer_map;
    bimap_type M_transfer_map_reordered;
    bimap_type M_transfer_map_elt;

    // container for storing the mesh marker names
    std::map<std::string, std::vector<int> > M_marker_names;

    //bimap_type M_reorder_map_nodes;
    //std::map<int,int> M_reorder_map_nodes; (not stored in memory)
    //bimap_type M_reorder_map_elements;
    //std::map<int,int> M_reorder_map_elements; (not stored in memory)

    std::vector<int> M_map_elements;
    std::vector<int> M_map_nodes;

    std::map<std::string,std::pair<boost::mpi::timer,double> > timer;
};

} // Nextsim
#endif
