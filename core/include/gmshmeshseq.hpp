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

#include <environment.hpp>
#include <entities.hpp>
#include <meshpartition.hpp>

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
    void writeTofile(std::string const& filename);

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

    std::vector<double> bcoordX() const;
    std::vector<double> bcoordY() const;

    std::vector<double> bcoordX(double const& rotangle) const;
    std::vector<double> bcoordY(double const& rotangle) const;

    std::vector<double> meanLat() const;
    std::vector<double> meanLon() const;

    void setId(std::vector<int> const& newid);
    std::vector<int> id() const;

    void initGModel();
    void writeToGModel(std::string const& filename);
    void clear();

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

    meshPartitionOptions M_partition_options;
    GModel*  M_gmodel;
    std::map<std::string,std::pair<boost::mpi::timer,double> > timer;
};

} // Nextsim
#endif
