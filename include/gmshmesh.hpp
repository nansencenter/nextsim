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
#include <algorithm>
#include <assert.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <GModel.h>
#include <MElement.h>

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

    //~GmshMesh();

    void readFromFile(std::string const& filename);

    std::string const& version() const {return M_version;}
    std::map<int, point_type > const& nodes () const {return M_nodes;}
    std::map<int, element_type > const& elements () const {return M_elements;}

    int numNodes() const {return M_num_nodes;}
    int numElements() const {return M_num_elements;}

private:

    std::map<int, point_type > M_nodes;
    std::map<int, element_type > M_elements;
    std::string M_version;
    int M_num_nodes;
    int M_num_elements;
};

} // Nextsim
#endif
