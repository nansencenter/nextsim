/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

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

extern "C"
{
#include <mapx.h>
}


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

    void readFromFile(std::string const& filename);
    void writeTofile(std::string const& filename);
    void partition(std::string const& filename, std::string const& partitioner = "chaco");
    void move(std::vector<double> const& um, double factor);
    void reorder(bimap_type const& rmap_nodes, bimap_type const& rmap_elements);

    std::string const& version() const {return M_version;}
    std::string const& ordering() const {return M_ordering;}

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

private:

    std::string M_version;
    std::string M_ordering;
    std::vector<point_type> M_nodes;
    std::vector<element_type> M_triangles;
    std::vector<element_type> M_edges;
    int M_num_nodes;
    int M_num_triangles;
    int M_num_edges;
};

} // Nextsim
#endif
