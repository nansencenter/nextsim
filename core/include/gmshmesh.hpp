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
#include <algorithm>
#include <assert.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <environment.hpp>
#include <MElement.h>

extern "C"
{
#include <mapx.h>
}


namespace Nextsim
{

namespace entities
{
class GMSHPoint
{
public:

    GMSHPoint()
        :
        id(0),
        coords()
    {}

    std::vector<double> coords;
    int id;

};

class GMSHElement
{
public:

    GMSHElement()
        :
        number( 0 ),
        type( MSH_PNT ),
        physical( 0 ),
        elementary( 0 ),
        numVertices(0),
        indices()
    {}

    GMSHElement( int n,
                 int t,
                 int p,
                 int e,
                 int _numVertices,
                 std::vector<int> const& _indices )
        :
        number( n ),
        type( t ),
        physical( p ),
        elementary( e ),
        numVertices( _numVertices ),
        indices( _indices )
    {}

    int number;
    int type;
    int physical;
    int elementary;

    // vertices
    int numVertices;
    std::vector<int> indices;
};
}

class GmshMesh
{
public:

    typedef Nextsim::entities::GMSHPoint point_type;
    typedef Nextsim::entities::GMSHElement element_type;

    GmshMesh();

    GmshMesh(std::vector<point_type> const& nodes,
             std::vector<element_type> const& edges,
             std::vector<element_type> const& triangles);

    GmshMesh(std::vector<point_type> const& nodes,
             std::vector<element_type> const& triangles);

    GmshMesh(GmshMesh const& mesh);

    void readFromFile(std::string const& filename);
    void writeTofile(std::string const& filename);
    void writeGeometry(std::string const& geofile, int nx, int ny, double xmin, double ymin, double dx, double dy);
    void move(std::vector<double> const& um, double factor);
    void setId(std::vector<int> const& new_id);

    std::string const& version() const {return M_version;}
    std::string const& ordering() const {return M_ordering;}
    std::string const& projfile() const {return M_projection_file;}
    std::vector<point_type> const& nodes() const {return M_nodes;}
    //std::vector<element_type> const& elements() const {return M_elements;}
    std::vector<element_type> const& triangles() const {return M_triangles;}
    std::vector<element_type> const& edges() const {return M_edges;}

    int numNodes() const {return M_num_nodes;}
    //int numElements() const {return M_num_elements;}
    int numTriangles() const {return M_num_triangles;}
    int numEdges() const {return M_num_edges;}

    void setOrdering(std::string const& order) {M_ordering=order;}
    void setNodes(std::vector<point_type> const& nodes) {M_nodes=nodes;}
    //void setElements(std::vector<element_type> const& elements) {M_elements=elements;}
    void setEdges(std::vector<element_type> const& edges) {M_edges=edges;}
    void setTriangles(std::vector<element_type> const& triangles) {M_triangles=triangles;}
    void setNumNodes(int const& nnodes) {M_num_nodes=nnodes;}
    //void setNumElements(int const& nelts) {M_num_elements=nelts;}
    void setNumEdges(int const& nlns) {M_num_edges=nlns;}
    void setNumTriangles(int const& ntrs) {M_num_triangles=ntrs;}
    void stereographicProjection();

    void update(std::vector<point_type> const& nodes,
                std::vector<element_type> const& triangles);

    std::vector<int> indexTr() const;
    std::vector<double> coordX() const;
    std::vector<double> coordY() const;
    std::vector<int> id() const;

    std::vector<double> coordX(double const& rotangle) const;
    std::vector<double> coordY(double const& rotangle) const;

    std::vector<double> bcoordX() const;
    std::vector<double> bcoordY() const;

    std::vector<double> bcoordX(double const& rotangle) const;
    std::vector<double> bcoordY(double const& rotangle) const;

    std::vector<double> meanLat() const;
    std::vector<double> meanLon() const;

    std::vector<double> lat() const;
    std::vector<double> lon() const;

private:

    std::string M_version;
    std::string M_ordering;
    std::string M_projection_file;
    std::vector<point_type> M_nodes;
    //std::vector<element_type> M_elements;
    std::vector<element_type> M_triangles;
    std::vector<element_type> M_edges;
    int M_num_nodes;
    //int M_num_elements;
    int M_num_triangles;
    int M_num_edges;
};

} // Nextsim
#endif
