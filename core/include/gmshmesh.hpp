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
//#include <GModel.h>
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
        numPartitions( 1 ),
        partition( 0 ),
        ghosts(),
        is_on_processor( false ),
        is_ghost( false ),
        ghost_partition_id( -1 ),
        numVertices(0),
        indices()
    {}

    GMSHElement( int n,
                 int t,
                 int p,
                 int e,
                 int _numPartitions,
                 int _partition,
                 std::vector<int> const& _ghosts,
                 int _numVertices,
                 std::vector<int> const& _indices,
                 int worldcommrank,
                 int worldcommsize)
        :
        number( n ),
        type( t ),
        physical( p ),
        elementary( e ),
        numPartitions( _numPartitions ),
        partition( (_partition % worldcommsize) ),
        ghosts( _ghosts ),
        is_on_processor( false ),
        is_ghost( false ),
        ghost_partition_id( partition ),
        numVertices( _numVertices ),
        indices( _indices )
    {
        setPartition(worldcommrank,worldcommsize);
    }


    bool isOnProcessor() const { return is_on_processor; }
    bool isGhost() const { return is_ghost; }
    int ghostPartitionId() const { return ghost_partition_id; }

    void setPartition(int worldcommrank, int worldcommsize)
    {
        // maybe proc id not start to 0
        for ( auto _itghost=ghosts.begin(),_enghost=ghosts.end() ; _itghost!=_enghost ; ++_itghost )
            *_itghost = ( (*_itghost) % worldcommsize);

        if ( worldcommsize == 1 )
        {
            is_on_processor = true;
            is_ghost = false;
        }
        else if ( worldcommrank == partition )
        {
            is_on_processor = true;
            is_ghost = false;
        }
        else
        {
            // is the element a ghost cell
            // look into ghosts if 'partition' is present
            auto it = std::find( ghosts.begin(), ghosts.end(), worldcommrank );
            if ( it != ghosts.end() )
            {
                is_on_processor = true;
                is_ghost = true;
                ghost_partition_id = partition;
            }
        }
    }

    int number;
    int type;
    int physical;
    int elementary;

    // partitioning info
    int numPartitions;
    int partition;
    std::vector<int> ghosts;
    bool is_on_processor;
    bool is_ghost;
    int ghost_partition_id;

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

    typedef boost::bimap<int,int> bimap_type;
    typedef bimap_type::value_type position;


    GmshMesh( Communicator const& comm = Environment::comm() );

    GmshMesh(std::vector<point_type> const& nodes,
             std::vector<element_type> const& edges,
             std::vector<element_type> const& triangles,
             Communicator const& comm = Environment::comm());

    void readFromFile(std::string const& filename);
    void writeTofile(std::string const& filename);
    void move(std::vector<double> const& um, double factor);
    void nodalGrid();

    Communicator const& comm() const { return M_comm; }
    std::string const& version() const {return M_version;}
    std::string const& ordering() const {return M_ordering;}
    std::vector<point_type> const& nodes() const {return M_nodes;}
    //std::vector<element_type> const& elements() const {return M_elements;}
    std::vector<element_type> const& triangles() const {return M_triangles;}
    std::vector<element_type> const& edges() const {return M_edges;}

    int numNodes() const {return M_num_nodes;}
    //int numElements() const {return M_num_elements;}
    int numTriangles() const {return M_num_triangles;}
    int numEdges() const {return M_num_edges;}

    void setCommunicator(Communicator const& comm) {M_comm=comm;}
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

    std::vector<int> const& localDofWithoutGhost() const {return M_local_dof_without_ghost;}
    std::vector<int> const& localDofWithGhost() const {return M_local_dof_with_ghost;}
    std::vector<int> const& localGhost() const {return M_local_ghost;}
    //bimap_type const& transferMap() const {return M_transfer_map;}
    bimap_type const& transferMap() const {return M_transfer_map_reordered;}
    bimap_type const& transferMapReordered() const {return M_transfer_map_reordered;}
    bool isValid(element_type const& elt) const;

    std::vector<int> indexTr() const;
    std::vector<double> coordX() const;
    std::vector<double> coordY() const;

    std::vector<double> coordX(double const& rotangle) const;
    std::vector<double> coordY(double const& rotangle) const;

    std::vector<double> bCoordX() const;
    std::vector<double> bCoordY() const;
    std::vector<double> meanLat() const;
    std::vector<double> meanLon() const;

private:

    Communicator M_comm;
    std::string M_version;
    std::string M_ordering;
    std::vector<point_type> M_nodes;
    //std::vector<element_type> M_elements;
    std::vector<element_type> M_triangles;
    std::vector<element_type> M_edges;
    int M_num_nodes;
    //int M_num_elements;
    int M_num_triangles;
    int M_num_edges;

    std::vector<int> M_local_dof_with_ghost;
    std::vector<int> M_local_dof_without_ghost;
    std::vector<int> M_local_ghost;
    bimap_type M_transfer_map;
    bimap_type M_transfer_map_reordered;

    // std::vector<int> local_node_without_ghost;
    // std::vector<int> local_node_with_ghost;
    // std::vector<int> local_ghost;
};

} // Nextsim
#endif
