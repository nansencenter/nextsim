/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

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

extern "C"
{
#include <mapx.h>
}


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

    void readFromFile(std::string const& filename);
    void writeTofile(std::string const& filename);
    void move(std::vector<double> const& um, double factor);
    void nodalGrid();

    Communicator const& comm() const { return M_comm; }
    std::string const& version() const {return M_version;}
    std::string const& ordering() const {return M_ordering;}
    std::map<int, point_type > const& nodes() const {return M_nodes;}
    std::vector<element_type> const& triangles() const {return M_triangles;}
    std::vector<element_type> const& edges() const {return M_edges;}

    int numGlobalNodes() const {return M_global_num_nodes;}
    int numNodes() const {return M_num_nodes;}
    int numTriangles() const {return M_num_triangles;}
    int numEdges() const {return M_num_edges;}

    int numLocalNodesWithoutGhost() const {return M_nldof_without_ghost;}
    int numLocalNodesWithGhost() const {return M_nldof_with_ghost;}
    int numLocalGhost() const {return M_nlghost;}

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

    //bimap_type const& transferMap() const {return M_transfer_map;}
    bimap_type const& transferMap() const {return M_transfer_map_reordered;}
    bimap_type const& transferMapReordered() const {return M_transfer_map_reordered;}

    bimap_type const& mapNodes() const {return M_reorder_map_nodes;}
    bimap_type const& mapElements() const {return M_reorder_map_elements;}

    std::vector<int> indexTr() const;

    std::vector<double> coordX() const;
    std::vector<double> coordY() const;

    std::vector<double> coordX(double const& rotangle) const;
    std::vector<double> coordY(double const& rotangle) const;

    std::vector<double> bCoordX() const;
    std::vector<double> bCoordY() const;

    std::vector<double> bCoordX(double const& rotangle) const;
    std::vector<double> bCoordY(double const& rotangle) const;

    std::vector<double> meanLat() const;
    std::vector<double> meanLon() const;

    std::vector<int> indexTrPartition() const;
    std::vector<double> coordXPartition() const;
    std::vector<double> coordYPartition() const;


private:

    Communicator M_comm;
    std::string M_version;
    std::string M_ordering;

    std::vector<point_type> M_nodes_vec;
    std::map<int, point_type > M_nodes;
    std::vector<element_type> M_triangles;
    std::vector<element_type> M_edges;

    int M_global_num_nodes;
    int M_num_nodes;
    int M_num_triangles;
    int M_num_edges;

    int M_nldof_with_ghost;
    int M_nldof_without_ghost;
    int M_nlghost;

    int M_num_triangles_without_ghost;

    std::vector<int> M_local_dof_with_ghost;
    std::vector<int> M_local_dof_without_ghost;
    std::vector<int> M_local_ghost;

    bimap_type M_transfer_map;
    bimap_type M_transfer_map_reordered;

    bimap_type M_reorder_map_nodes;
    bimap_type M_reorder_map_elements;
};

} // Nextsim
#endif
