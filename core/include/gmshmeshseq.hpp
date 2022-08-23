/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gmshmeshseq.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Tue Jun 21 14:26:34 2016
 */

#ifndef __GmshMeshSeq_HPP
#define __GmshMeshSeq_HPP 1

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <assert.hpp>
#include <boost/bimap.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <environment.hpp>
#include <entities.hpp>
#include <meshpartition.hpp>
#include "debug.hpp"
#include <meshmath.hpp>

extern "C"
{
#include <mapx.h>
}

int PartitionMesh( GModel *const model, meshPartitionOptions &options );

namespace Nextsim
{
class GmshMeshSeq
{
public:

    typedef Nextsim::entities::GMSHPoint point_type;
    typedef Nextsim::entities::GMSHElement element_type;

    typedef boost::bimap<int,int> bimap_type;
    typedef bimap_type::value_type position;


    GmshMeshSeq();

    GmshMeshSeq(std::vector<point_type> const& nodes,
                std::vector<element_type> const& edges,
                std::vector<element_type> const& triangles);

    GmshMeshSeq(std::vector<point_type> const& nodes,
                std::vector<element_type> const& triangles);

    GmshMeshSeq(GmshMeshSeq const& mesh);

    //~GmshMeshSeq();

    void readFromFile(std::string const& filename);
    void writeToFile(std::string const& filename);

    void partition(std::string const& filename,
                   mesh::Partitioner const& partitioner=mesh::Partitioner::METIS,
                   mesh::PartitionSpace const& space=mesh::PartitionSpace::MEMORY,
                   std::string const& format="ascii");

    void move(std::vector<double> const& um, double factor);
    void reorder(bimap_type const& rmap_nodes, bimap_type const& rmap_elements);

    std::string const& version() const {return M_version;}
    std::string const& ordering() const {return M_ordering;}
    std::string const& mppfile() const {return M_mppfile;}

    std::vector<point_type> const& nodes() const {return M_nodes;}
    std::vector<element_type> const& triangles() const {return M_triangles;}
    std::vector<element_type> const& edges() const {return M_edges;}

    std::map<std::string, std::vector<int> > const& markerNames() const {return M_marker_names;}

    int numNodes() const {return M_num_nodes;}
    int numTriangles() const {return M_num_triangles;}
    int numEdges() const {return M_num_edges;}

    void setOrdering(std::string const& order) {M_ordering=order;}

    void setNodes(std::vector<point_type> const& nodes) {M_nodes=nodes;}
    void setEdges(std::vector<element_type> const& edges) {M_edges=edges;}
    void setTriangles(std::vector<element_type> const& triangles) {M_triangles=triangles;}

    void setNumNodes(int const& nnodes) {M_num_nodes=nnodes;}
    void setNumEdges(int const& nlns) {M_num_edges=nlns;}
    void setNumTriangles(int const& ntrs) {M_num_triangles=ntrs;}

    void update(std::vector<point_type> const& nodes,
                std::vector<element_type> const& triangles,
                std::string const& ordering = "gmsh");


    void stereographicProjection();

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

    std::vector<std::vector<double>>
        vertices(std::vector<int> const& indices) const;
    std::vector<std::vector<double>>
        vertices(std::vector<int> const& indices,
                std::vector<double> const& um, double factor) const;

    std::vector<double> lat() const;
    std::vector<double> lon() const;

    void setId(std::vector<int> const& newid);
    std::vector<int> id() const;

    void initGModel();
    void writeToGModel();
    void clear();

    double jacobian(element_type const& element) const
    { return MeshMath::jacobian(GmshMeshSeq::vertices(element.indices)); }

    double jacobian(element_type const& element,
                    std::vector<double> const& um, double factor = 1.) const
    { return MeshMath::jacobian(GmshMeshSeq::vertices(element.indices, um, factor)); }

    //------------------------------------------------------------------------------------------------------
    //! Calculates the area of triangular mesh elements.
    //! Called by the advect() and other functions.
    double measure(element_type const& element) const
    { return (1./2)*std::abs(GmshMeshSeq::jacobian(element)); }

    double measure(element_type const& element,
            std::vector<double> const& um, double factor = 1.) const
    { return (1./2)*std::abs(GmshMeshSeq::jacobian(element,um,factor)); }

private:

    void partitionMemory(std::string const& filename,
                         mesh::Partitioner const& partitioner=mesh::Partitioner::METIS,
                         std::string const& format="ascii");

    void partitionDisk(std::string const& filename,
                       mesh::Partitioner const& partitioner=mesh::Partitioner::METIS,
                       std::string const& format="ascii");


private:

    std::string M_version;
    std::string M_ordering;
    std::string M_mppfile;
    std::vector<point_type> M_nodes;
    std::vector<element_type> M_triangles;
    std::vector<element_type> M_edges;
    int M_num_nodes;
    int M_num_triangles;
    int M_num_edges;

    Communicator M_comm;
    LogLevel M_log_level;
    bool M_log_all;

    // container for storing the mesh marker names
    std::map<std::string, std::vector<int> > M_marker_names;

    //meshPartitionOptions M_partition_options;
    GModel*  M_gmodel;
    std::map<std::string,std::pair<boost::mpi::timer,double> > timer;
};

} // Nextsim
#endif
