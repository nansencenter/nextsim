/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   entities.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Tue Jun 21 14:52:48 2016
 */

#ifndef __Entities_HPP
#define __Entities_HPP 1

#include <vector>
#include <algorithm>
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

    GMSHElement( int n,
                 int t,
                 int p,
                 int e,
                 int _numPartitions,
                 int _partition,
                 std::vector<int> const& _ghosts,
                 std::vector<bool> const& _ghostNodes,
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
        ghostNodes( _ghostNodes ),
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
    std::vector<bool> ghostNodes;
    bool is_on_processor;
    bool is_ghost;
    int ghost_partition_id;

    // vertices
    int numVertices;
    std::vector<int> indices;
};
}

} // Nextsim
#endif
