/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gmshmesh.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Wed Jul 29 15:26:47 2015
 */

#include <gmshmesh.hpp>
#include <boost/mpi.hpp>
#include <boost/format.hpp>
#include <boost/serialization/string.hpp>
#include <vector>
#include <string>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <numeric>

namespace Nextsim
{
GmshMesh::GmshMesh(Communicator const& comm)
    :
    M_comm(comm),
    M_version("2.2"),
    M_ordering("gmsh"),
    M_nodes(),
    M_triangles(),
    M_edges(),
    M_nodes_vec(),
    M_local_dof_with_ghost(),
    M_local_dof_without_ghost(),
    M_local_ghost(),
    M_triangles_id_with_ghost(),
    M_num_nodes(0),
    M_num_triangles(0),
    M_num_edges(0),
    M_nldof_with_ghost(),
    M_nldof_without_ghost(),
    M_nlghost(),
    M_num_triangles_without_ghost(),
    M_transfer_map(),
    M_transfer_map_reordered(),
    M_transfer_map_elt(),
    //M_reorder_map_nodes(),
    //M_reorder_map_elements(),
    M_map_nodes(),
    M_map_elements(),
    timer(),
    M_mppfile(Environment::nextsimMppfile()),
    M_log_level(Environment::logLevel()),
    M_log_all(Environment::logAll())
{}

GmshMesh::GmshMesh(GmshMesh const& mesh)
    :
    M_comm(mesh.M_comm),
    M_version(mesh.M_version),
    M_ordering(mesh.M_ordering),
    M_mppfile(mesh.M_mppfile),
    M_log_level(mesh.M_log_level),
    M_log_all(mesh.M_log_all),
    M_nodes(mesh.M_nodes),
    M_triangles(mesh.M_triangles),
    M_edges(mesh.M_edges),
    M_nodes_vec(mesh.M_nodes_vec),
    M_local_dof_with_ghost(mesh.M_local_dof_with_ghost),
    M_local_dof_without_ghost(mesh.M_local_dof_without_ghost),
    M_local_ghost(mesh.M_local_ghost),
    M_triangles_id_with_ghost(mesh.M_triangles_id_with_ghost),
    M_num_nodes(mesh.M_num_nodes),
    M_num_triangles(mesh.M_num_triangles),
    M_num_edges(mesh.M_num_edges),
    M_nldof_with_ghost(mesh.M_nldof_with_ghost),
    M_nldof_without_ghost(mesh.M_nldof_without_ghost),
    M_nlghost(mesh.M_nlghost),
    M_num_triangles_without_ghost(mesh.M_num_triangles_without_ghost),
    M_transfer_map(mesh.M_transfer_map),
    M_transfer_map_reordered(mesh.M_transfer_map_reordered),
    //M_reorder_map_nodes(mesh.M_reorder_map_nodes),
    //M_reorder_map_elements(mesh.M_reorder_map_elements)
    M_map_nodes(mesh.M_map_nodes),
    M_map_elements(mesh.M_map_elements)
{}


void
GmshMesh::readFromFile(std::string const& gmshmshfile, std::string const& format)
{
    LOG(DEBUG)<<"Reading Msh file "<< gmshmshfile <<"\n";

    std::ifstream ifs ( gmshmshfile.c_str() );

    if ( !ifs.is_open() )
    {
        std::ostringstream ostr;
        std::cout << "Invalid file name " << gmshmshfile << " (file not found)\n";
        ostr << "Invalid file name " << gmshmshfile << " (file not found)\n";
        throw std::invalid_argument( ostr.str() );
    }

    if (format == "binary")
        this->readFromFileBinary(ifs);
    else if (format == "ascii")
        this->readFromFileASCII(ifs);
    else
    {
        std::cout << "invalid mesh file format"<<"\n";
        throw std::logic_error("invalid mesh file format");
    }

    // create nodal partitions
    timer["in.nodal"].first.restart();
    if (M_comm.size() > 1)
        this->nodalGrid();

    LOG(DEBUG)<<"-------------------INSIDE: NODALGRID done in "<< timer["in.nodal"].first.elapsed() <<"s\n";
}//readFromFile

void
GmshMesh::readFromFileASCII(std::ifstream& ifs)
{
    char buf[256];
    ifs >> buf;

    std::string theversion;
    double version = 2.2;

    if (std::string( buf ) == "$MeshFormat")
    {
        int format, size;
        ifs >> theversion >> format >> size;

        LOG(DEBUG) << "GMSH mesh file version : " << theversion
                  << " format: " << (format?"binary":"ascii")
                  << " size of double: " << size << "\n";

        ASSERT(boost::lexical_cast<double>( theversion ) >= 2, "Nextsim supports only Gmsh version >= 2");

        version = boost::lexical_cast<double>( theversion );

        ifs >> buf;

        ASSERT(std::string( buf ) == "$EndMeshFormat","invalid file format entry");

        ifs >> buf;

        LOG(DEBUG) << "[gmshmesh::reading] " << buf << " (expect $PhysicalNames)\n";

        if ( std::string( buf ) == "$PhysicalNames" )
        {
            int nnames;
            ifs >> nnames;

            for ( int n = 0; n < nnames; ++n )
            {
                int id, topodim;
                std::string name;

                ifs >> topodim >> id >> name;

                boost::trim( name );
                boost::trim_if( name,boost::is_any_of( "\"" ) );

                LOG(DEBUG) << "[gmshmesh::reading] topodim: "  << topodim << " id: " << id << " name: " << name << "\n";

                std::vector<int> marker_data = {id, topodim};
                M_marker_names.insert(std::make_pair(name,marker_data));
            }

            ifs >> buf;
            ASSERT(std::string( buf ) == "$EndPhysicalNames","invalid file format entry");

            ifs >> buf;
        }

    }

    // Read NODES

    if ( !( std::string( buf ) == "$NOD" ||
            std::string( buf ) == "$Nodes" ||
            std::string( buf ) == "$ParametricNodes") )
    {
        LOG(WARNING)<< "invalid nodes string '" << buf << "' in gmsh importer. It should be either $Nodes.\n";
    }

    bool has_parametric_nodes = ( std::string( buf ) == "$ParametricNodes" );
    unsigned int __n;
    ifs >> __n;

    M_num_nodes = __n;

    M_global_num_nodes_from_serial = M_num_nodes;

    LOG(DEBUG) << "Reading "<< __n << " nodes\n";

    M_nodes_vec.resize(__n);
    std::vector<double> coords(3,0);

    for ( unsigned int __i = 0; __i < __n; ++__i )
    {
        int id = 0;

        ifs >> id
             >> coords[0]
             >> coords[1]
             >> coords[2];

        M_nodes_vec[id-1].id = id;
        M_nodes_vec[id-1].coords = coords;
    }

    ifs >> buf;

    ASSERT(std::string( buf ) == "$EndNodes","invalid end nodes string");

    // Read ELEMENTS

    ifs >> buf;

    ASSERT(std::string( buf ) == "$Elements","invalid elements string");

    int numElements;
    ifs >> numElements;

    M_global_num_elements_from_serial = numElements;

    LOG(DEBUG) << "Reading " << numElements << " elements...\n";
    std::map<int,int> __gt;

    int cpt_edge = 0;
    int cpt_triangle = 0;
    int num_edge = 0;
    int num_edge_diff = 0;
    bool first_triangle = true;

    for(int i = 0; i < numElements; i++)
    {
        int number, type, physical = 0, elementary = 0;
        int const numVertices = 3;// only use triangular elements here
        std::vector<int> ghosts;
        std::vector<bool> ghostNodes;
        int numTags;
        int partition = (this->comm().size()>1)?this->comm().rank():0;

        ifs >> number  // elm-number
             >> type // elm-type
             >> numTags; // number-of-tags

        // Skip non-triangular elements (type != 2)
        if (type != 2)
        {
            ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            ++num_edge;
            continue;
        }

        int numPartitions = 1;

        for(int j = 0; j < numTags; j++)
        {
            int tag;
            ifs >> tag;
            if(j == 0) physical = tag;
            else if(j == 1) elementary = tag;
            else if((j == 2) && (numTags > 3)) numPartitions = tag;
            else if(j == 3) partition = tag-1;
            else if((j >= 4) && (j < 4 + numPartitions - 1)) ghosts.push_back((-tag)-1);
        }

        std::vector<int> indices(numVertices);
        for(int j = 0; j < numVertices; j++)
        {
            ifs >> indices[j];
        }

        if (M_ordering=="bamg")
        {
            std::next_permutation(indices.begin()+1,indices.end());
        }

        if (first_triangle)
        {
            if (num_edge == 0)
            {
                num_edge_diff = 0;
            }
            else
            {
                num_edge_diff = (number == 1) ? 0 : num_edge;
            }

            first_triangle = false;
        }

        number = number - num_edge_diff;

        Nextsim::entities::GMSHElement gmshElt( number,
                                                type,
                                                physical,
                                                elementary,
                                                numPartitions,
                                                partition,
                                                ghosts,
                                                ghostNodes,
                                                numVertices,
                                                indices,
                                                this->comm().rank(),
                                                this->comm().size());

        if (gmshElt.isOnProcessor() == false)
            continue;

        if (type == 2)
        {
            M_triangles.push_back(gmshElt);
            ++cpt_triangle;
        }
        else if (type == 1)
        {
            M_edges.push_back(gmshElt);
            ++cpt_edge;
        }

        if ( __gt.find( type ) != __gt.end() )
            ++__gt[ type ];
        else
            __gt[type]=1;

    } // element description loop

    M_global_num_elements_from_serial = M_global_num_elements_from_serial - num_edge;

    for ( auto const& it : __gt )
    {
        const char* name = Nextsim::entities::getElementTypeName(it.first);

        if (std::string(name) == "Triangle 3")
            M_num_triangles = it.second;
        else if (std::string(name) == "Line 2")
            M_num_edges = it.second;
    }

    // make sure that we have read everything
    ifs >> buf;

    ASSERT(std::string( buf ) == "$EndElements","invalid end elements string");
}// readFromFileASCII


void
GmshMesh::readFromFileBinary(std::ifstream& ifs)
{
    char buf[256];
    ifs >> buf;

    std::string theversion;
    double version = 2.2;
    bool swap = false;

    if (std::string( buf ) == "$MeshFormat")
    {
        int format, size;
        ifs >> theversion >> format >> size;

        LOG(DEBUG) << "GMSH mesh file version : " << theversion
                  << " format: " << (format?"binary":"ascii")
                  << " size of double: " << size << "\n";

        ASSERT(boost::lexical_cast<double>( theversion ) >= 2, "Nextsim supports only Gmsh version >= 2");

        version = boost::lexical_cast<double>( theversion );

        // ----------------------------------------------------------------------
        char c=ifs.get();
        ASSERT( c == '\n', "Invalid character");

        int one;
        ifs.read(reinterpret_cast<char*>(&one), sizeof(int));

        if(one != 1)
        {
            swap = true;
            LOG(DEBUG) << "one before swap : " << one << "\n";
            if(swap) GmshMesh::SwapBytes(&one, sizeof(int), 1); // Note: removed (char*) cast
            LOG(DEBUG) << "one after swap : " << one << "\n";
            LOG(DEBUG) <<"Swapping bytes from binary file (to be done)\n";
        }
        // ----------------------------------------------------------------------

        ifs >> buf;

        ASSERT(std::string( buf ) == "$EndMeshFormat","invalid file format entry");

        ifs >> buf;

        LOG(DEBUG) << "[gmshmesh::reading] " << buf << " (expect $PhysicalNames)\n";

    }

    // Read NODES

    if ( !( std::string( buf ) == "$NOD" ||
            std::string( buf ) == "$Nodes" ||
            std::string( buf ) == "$ParametricNodes") )
    {
        LOG(WARNING)<< "invalid nodes string '" << buf << "' in gmsh importer. It should be either $Nodes.\n";
    }

    bool has_parametric_nodes = ( std::string( buf ) == "$ParametricNodes" );
    unsigned int __n;
    ifs >> __n;

    // eat  '\n' in binary mode otherwise the next binary read will get screwd
    ifs.get();

    M_num_nodes = __n;

    M_global_num_nodes_from_serial = M_num_nodes;

    LOG(DEBUG) << "Reading "<< __n << " nodes\n";

    M_nodes_vec.resize(__n);
    std::vector<double> coords(3,0);

    for ( unsigned int __i = 0; __i < __n; ++__i )
    {
        int id = 0;

        ifs.read(reinterpret_cast<char*>(&id), sizeof(int));
        if(swap) GmshMesh::SwapBytes(&id, sizeof(int), 1);
        ifs.read(reinterpret_cast<char*>(&coords[0]), 3*sizeof(double));
        if(swap) GmshMesh::SwapBytes(&coords[0], sizeof(double), 3);

        M_nodes_vec[id-1].id = id;
        M_nodes_vec[id-1].coords = coords;
    }

    // eat  '\n' in binary mode otherwise the next binary read will get screwd
    ifs.get();

    ifs >> buf;

    // make sure that we have read all the points

    ASSERT(std::string( buf ) == "$EndNodes","invalid end nodes string");

    // Read ELEMENTS

    ifs >> buf;

    ASSERT(std::string( buf ) == "$Elements","invalid elements string");

    int numElements;
    ifs >> numElements;

    // eat  '\n' in binary mode otherwise the next binary read will get screwd
    ifs.get();

    M_global_num_elements_from_serial = numElements;

    LOG(DEBUG) << "Reading " << numElements << " elements...\n";

    std::map<int,int> __gt;

    int cpt_edge = 0;
    int cpt_triangle = 0;
    int num_edge = 0;
    int num_edge_diff = 0;
    bool first_triangle = true;

    int cptii = 0;

    int numElementsPartial = 0;
    while(numElementsPartial < numElements)
    {
        int header[3];

        ifs.read(reinterpret_cast<char*>(&header), 3*sizeof(int));
        if(swap) GmshMesh::SwapBytes(header, sizeof(int), 3);  // Note: removed (char*) cast

        int type = header[0];
        int numElems = header[1];
        int numTags = header[2];

        int numVertices = Nextsim::entities::getNumVerticesForElementType(type);
        if (numVertices == 0)
        {
            std::cout << "Invalid GMSH element type " << type << "\n";
            throw std::logic_error("Invalid GMSH element type");
        }
        const char* name = Nextsim::entities::getElementTypeName(type);

        if ( numVertices <= 0 )
        {
            std::cout << "Unsupported element type " << name << "\n";
            throw std::logic_error("Unsupported element type");
        }

        unsigned int n = 1 + numTags + numVertices;

        if (type != 2)
        {
            int off = sizeof(int)*(1 + numTags + numVertices);
            ifs.seekg(off, std::ios::cur); // skip from the direction (beginning/current/end) position of the file

            numElementsPartial += numElems;
            ++num_edge;
            continue;
        }

        std::vector<int> data(n);
        std::vector<int> indices(numVertices);
        std::vector<int> ghosts;


        for(int i = 0; i < numElems; i++)
        {
            ghosts.clear();
            std::vector<bool> ghostNodes;

            ifs.read(reinterpret_cast<char*>(data.data()), n*sizeof(int));
            if(swap) GmshMesh::SwapBytes(data.data(), sizeof(int), n);

            int number = data[0];
            int physical = (numTags > 0) ? data[1] : 0;
            int elementary = (numTags > 1) ? data[2] : 0;
            int numPartitions = (version >= 2.2 && numTags > 3) ? data[3] : 1;
            int partition = (version < 2.2 && numTags > 2) ? data[3]-1 : (version >= 2.2 && numTags > 3) ? data[4]-1 : 0;

            if(numPartitions > 1)
            {
                for(int j = 0; j < numPartitions - 1; j++)
                {
                    ghosts.push_back( (-data[5 + j]) -1 );
                }
            }

            std::copy( &data[numTags + 1], &data[numTags + 1]+numVertices, indices.begin() );

            if (M_ordering=="bamg")
            {
                std::next_permutation(indices.begin()+1,indices.end());
            }

            if (first_triangle)
            {
                if (num_edge == 0)
                    num_edge_diff = 0;
                else
                    num_edge_diff = (number == 1) ? 0 : num_edge;
                first_triangle = false;
            }
            number = number - num_edge_diff;

            Nextsim::entities::GMSHElement gmshElt( number,
                                                    type,
                                                    physical,
                                                    elementary,
                                                    numPartitions,
                                                    partition,
                                                    ghosts,
                                                    ghostNodes,
                                                    numVertices,
                                                    indices,
                                                    this->comm().rank(),
                                                    this->comm().size());

            if (gmshElt.isOnProcessor() == false)
                continue;

            if (type == 2)
            {
                M_triangles.push_back(gmshElt);
                ++cpt_triangle;
            }
            else if (type == 1)
            {
                M_edges.push_back(gmshElt);
                ++cpt_edge;
            }

            if ( __gt.find( type ) != __gt.end() )
            {
                ++__gt[ type ];
            }
            else
            {
                __gt[type]=1;
            }
        }

        numElementsPartial += numElems;

    } // while

    M_global_num_elements_from_serial = M_global_num_elements_from_serial - num_edge;

    for ( auto const& it : __gt )
    {
        const char* name = Nextsim::entities::getElementTypeName(it.first);
        if (std::string(name) == "Triangle 3")
            M_num_triangles = it.second;
        else if (std::string(name) == "Line 2")
            M_num_edges = it.second;
    }

    ifs >> buf;

    // make sure that we have read everything
    ifs >> buf;

    ASSERT(std::string( buf ) == "$EndElements","invalid end elements string");
}// readFromFileBinary


void
GmshMesh::move(std::vector<double> const& um, double factor)
{
    if ((um.size() != 0) && (factor != 0))
    {
        ASSERT(2*M_nodes.size()==um.size(),"invalid size of displacement vector");

        int cpt = 0;
        for (auto it=M_nodes.begin(), en=M_nodes.end(); it!=en; ++it)
        {
            it->second.coords[0] += factor*um[cpt];
            it->second.coords[1] += factor*um[cpt+M_num_nodes];

            ++cpt;
        }
    }
}

void
GmshMesh::update(std::vector<point_type> const& nodes,
                 std::vector<element_type> const& triangles,
                 int numTrianglesGlobal,
                 int numNodesGlobal)
{
#ifndef NDEBUG
    // Sanity check: every triangle index must be a valid 1-based subscript
    // into `nodes`. This is the invariant nodalGrid() silently depends on
    // via M_nodes[k+1] = M_nodes_vec[local_dof_with_ghost[k]-1].
    for (auto const& tri : triangles)
        for (int idx : tri.indices)
            ASSERT(idx >= 1 && idx <= (int)nodes.size(),
                    "GmshMesh::update: triangle index out of range of node array "
                    "(did you pass a local node subset where a global one was expected?)");
#endif

    M_nodes_vec = nodes;
    M_triangles = triangles;
    M_num_nodes = nodes.size();
    M_num_triangles = triangles.size();
    M_global_num_elements_from_serial = numTrianglesGlobal;

    // Explicit global node count, when known (parallel/MMG path).
    // Falls back to nodes.size() for the serial/root-mesh case, where
    // that value is already correct (no duplication across ranks to
    // resolve, since there's only one rank/mesh).
    M_global_num_nodes_from_serial = (numNodesGlobal >= 0) ? numNodesGlobal : M_num_nodes;
}


void
GmshMesh::stereographicProjection()
{
    // polar stereographic projection
    mapx_class *map;
    std::vector<char> str(M_mppfile.begin(), M_mppfile.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    int cpt = 0;
    for (auto it=M_nodes.begin(), en=M_nodes.end(); it!=en; ++it)
    {
        //compute latitude and longitude from cartesian coordinates
        double _x = it->second.coords[0];
        double _y = it->second.coords[1];
        double _z = it->second.coords[2];

        // compute radius
        double radius = std::sqrt(std::pow(_x,2.)+std::pow(_y,2.)+std::pow(_z,2.));

        double latitude = std::asin(_z/radius)*(180./PI);
        double longitude = std::atan2(_y,_x);

        longitude = longitude-2*PI*std::floor(longitude/(2*PI));
        longitude = longitude*(180./PI);

        double x_, y_;
        int status = forward_mapx(map,latitude,longitude,&x_,&y_);

        it->second.coords[0] = x_;
        it->second.coords[1] = y_;
        it->second.coords[2] = 0.0;

#if 0
        if (cpt < 10)
        {
            std::cout<<"latitude= "<< latitude <<"\n";
            std::cout<<"longitude= "<< longitude <<"\n";

            std::cout<<"        xcart= "<< it->second.coords[0] <<"\n";
            std::cout<<"        ycart= "<< it->second.coords[1] <<"\n";
            std::cout<<"        zcart= "<< it->second.coords[2] <<"\n";
        }
#endif

        ++cpt;
    }

    close_mapx(map);
}


void
GmshMesh::nodalGrid()
{
    // --- Step 1: identify nodes touched by ghost triangles (interface nodes) ---
    std::vector<int> ghosts_nodes_f;

    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
        if (it->is_ghost)
            for (int const& index : it->indices)
                ghosts_nodes_f.push_back(index);

    // --- Step 2: build this rank's tentative "owned" node set from its non-ghost triangles ---
    // A node is claimed here if either:
    //  (a) it never appears in a ghost triangle (unambiguously local), or
    //  (b) it does appear in a ghost triangle elsewhere, but THIS triangle
    //      (which is non-ghost) itself borders a ghost partition (it->ghosts
    //      non-empty) -- i.e. this is an interface node visited via a
    //      boundary-adjacent triangle.
    // NB: interface nodes can be claimed this way by MULTIPLE ranks
    // simultaneously (once per rank whose local triangles touch them).
    // That over-claiming is intentional here and is resolved later via
    // allGather + cross-rank de-duplication (see below).
    std::sort(ghosts_nodes_f.begin(), ghosts_nodes_f.end());
    ghosts_nodes_f.erase(std::unique( ghosts_nodes_f.begin(), ghosts_nodes_f.end() ), ghosts_nodes_f.end());
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        if (!it->is_ghost)
        {
            for (int const& index : it->indices)
            {
                bool const touches_ghost = std::binary_search(
                        ghosts_nodes_f.begin(), ghosts_nodes_f.end(), index);
                bool const claims_node = !touches_ghost || (it->ghosts.size() > 0);
                if (claims_node)
                    M_local_dof_without_ghost.push_back(index);
            }
        }
    }

    std::sort(M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end());
    M_local_dof_without_ghost.erase(std::unique( M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end() ), M_local_dof_without_ghost.end());

    // --- Step 3: seed M_local_ghost with nodes seen on triangles from
    // higher-or-equal-ranked partitions but not already claimed as owned ---
    std::vector<int> all_local_nodes;

    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        if (M_comm.rank() <= it->partition) // to check
        {
            bool is_found = false;

            for (int i=0; i<3; ++i)
            {
                if (std::binary_search(M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end(),it->indices[i]))
                {
                    is_found = true;
                    break;
                }
            }

            if (!is_found)
                continue;

            for (int i=0; i<3; ++i)
            {
                all_local_nodes.push_back(it->indices[i]);
            }
        }
    }

    std::sort(all_local_nodes.begin(), all_local_nodes.end());
    all_local_nodes.erase(std::unique( all_local_nodes.begin(), all_local_nodes.end() ), all_local_nodes.end());

    std::set_difference(all_local_nodes.begin(), all_local_nodes.end(),
                        M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end(),
                        std::back_inserter(M_local_ghost));

    // --- Step 4: cross-rank reconciliation of ownership ---
    // `num_nodes` is the RAW sum of each rank's tentative "owned" count
    // (from Step 2), i.e. BEFORE removing cross-rank duplicate claims on
    // interface nodes. It is therefore >= the true global node count
    // whenever any interface node was claimed by more than one rank.
    std::vector<std::vector<int> > renumbering;
    int num_nodes;

    this->allGather(M_local_dof_without_ghost, renumbering, num_nodes);

    LOG(DEBUG)<<"num_nodes = "<< num_nodes <<"\n";

    // TODO: once GmshMesh::update() explicitly tracks the true global node
    // count (rather than overloading M_num_nodes for both the MMG global-
    // array case and the readFromFile per-rank case), replace M_num_nodes
    // below with that explicit value (e.g. M_global_num_nodes_from_serial).

    // Resolve duplicate ownership: for each pair of ranks (ii, jj<ii),
    // any node both claim is stripped from the higher-ranked renumbering[ii]
    // and reassigned as a ghost on that rank instead. Lower rank always
    // wins ownership. Applied in increasing `ii` order, this correctly
    // collapses n-way shared nodes (not just pairs) to a single owner.
    int const expected_global_num_nodes = M_global_num_nodes_from_serial;
    if (expected_global_num_nodes < num_nodes)
    {
        LOG(DEBUG)<<"---------------------------------------Post-processing needed for nodal mesh partition: "<< expected_global_num_nodes <<" < "<< num_nodes <<"\n";

        for (int ii=0; ii<renumbering.size(); ++ii)
        {
            for (int jj=0; jj<ii; ++jj)
            {
                std::vector<int> duplicated_dofs;

                std::set_intersection(renumbering[ii].begin(),renumbering[ii].end(),
                                      renumbering[jj].begin(),renumbering[jj].end(),
                                      std::back_inserter(duplicated_dofs));

                for (int kk=0; kk<duplicated_dofs.size(); ++kk)
                {
                    auto it = std::lower_bound(renumbering[ii].begin(), renumbering[ii].end(), duplicated_dofs[kk]);
                    if (it != renumbering[ii].end() && *it == duplicated_dofs[kk]) renumbering[ii].erase(it);

                    if (M_comm.rank() == ii) M_local_ghost.push_back(duplicated_dofs[kk]);
                }
            }
        }

        M_local_dof_without_ghost = renumbering[M_comm.rank()];
    } // check the nodal partition end

    std::sort(M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end());
    std::sort(M_local_ghost.begin(), M_local_ghost.end());

    std::copy_n(M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.size(), std::back_inserter(M_local_dof_with_ghost));
    std::copy_n(M_local_ghost.begin(), M_local_ghost.size(), std::back_inserter(M_local_dof_with_ghost));

    // Post-dedup, authoritative global node count: the sum of each rank's
    // FINAL (post-reconciliation) owned-node count. This is what M_ndof /
    // M_mesh.numGlobalNodes() downstream is expected to equal.
    int size = 0;
    for (int ii=0; ii<M_comm.size(); ++ii) size += renumbering[ii].size();
    std::vector<int> reorder(size+M_num_nodes+1);
    M_map_nodes.resize(size);

    // Build a globally-contiguous renumbering of dofs (needed for PETSc):
    // rank 0's owned nodes get global ids [1..n0], rank 1's get
    // [n0+1..n0+n1], etc. Each node gets two contiguous slots (u,v
    // velocity components), offset by `size`.
    int cpts = 0;
    int cpts_dom = 0;

    for (int ii=0; ii<M_comm.size(); ++ii)
    {
        int sr = renumbering[ii].size();
        for (int jj=0; jj<sr; ++jj)
        {
            // add first component (u) of velocity
            reorder[renumbering[ii][jj]] = cpts+1+cpts_dom;
            // add second component (v) of velocity
            reorder[renumbering[ii][jj]+M_num_nodes] = cpts+1+sr+cpts_dom;
            M_map_nodes[renumbering[ii][jj]-1] = cpts;
            ++cpts;
        }
        cpts_dom += renumbering[ii].size();
    }

    M_local_dof_with_ghost_init = M_local_dof_with_ghost;
    auto local_dof_with_ghost = M_local_dof_with_ghost;

    M_nldof_with_ghost = M_local_dof_with_ghost.size();
    M_nldof_without_ghost = M_local_dof_without_ghost.size();
    M_nlghost = M_local_ghost.size();

    M_local_dof_with_ghost.resize(2*M_nldof_with_ghost);
    M_local_dof_without_ghost.resize(2*M_nldof_without_ghost);

    for (int k=0; k<local_dof_with_ghost.size(); ++k)
    {
        // mapping from old global numbering to local numbering
        M_transfer_map.insert(position(local_dof_with_ghost[k],k+1));

        int rdof = reorder[local_dof_with_ghost[k]];
        int rdofv = reorder[local_dof_with_ghost[k]+M_num_nodes];

        // mapping from new global numbering to local numbering
        M_transfer_map_reordered.insert(position(rdof,k+1));

        M_local_dof_with_ghost[k] = rdof;
        M_local_dof_with_ghost[k+M_nldof_with_ghost] = rdofv;

        M_nodes[k+1] = M_nodes_vec[local_dof_with_ghost[k]-1];

        if (k < M_nldof_without_ghost)
        {
            M_local_dof_without_ghost[k] = rdof;
            M_local_dof_without_ghost[k+M_nldof_without_ghost] = rdofv;
        }
        else
        {
            M_local_ghost[k-M_nldof_without_ghost] = rdof;
        }
    }

    M_nodes_vec.clear();
    M_nodes_vec.shrink_to_fit();

    std::sort(M_local_ghost.begin(), M_local_ghost.end());
    M_global_num_nodes = size; // <- must be `size` (post-dedup), NOT M_num_nodes
                               //    or the raw `num_nodes` -- see history/bug
                               //    where using the wrong one caused
                               //    inconsistent M_ndof and downstream
                               //    segfaults in gatherNodalField/move().
    M_num_nodes = M_nodes.size();

    std::vector<int> triangles_num_without_ghost;
    std::vector<element_type> _triangles = M_triangles;
    M_triangles.resize(0);

    for (auto it=_triangles.begin(), end=_triangles.end(); it!=end; ++it)
    {
        if ((M_comm.rank() <= it->partition))
        {
            bool _test = false;

            for (int i=0; i<3; ++i)
            {
                if ((M_transfer_map.left.find(it->indices[i]) == M_transfer_map.left.end()))
                {
                    _test = true;
                    break;
                }
            }

            if (_test)
                continue;

            it->ghostNodes.assign(3,false);

            // new add
            for (int i=0; i<3; ++i)
            {
                int rdof = reorder[it->indices[i]];

                it->indices[i] = M_transfer_map.left.find(it->indices[i])->second;

                if (std::binary_search(M_local_ghost.begin(),M_local_ghost.end(),rdof))
                {
                    it->ghostNodes[i] = true;
                }
            }
            // end

            M_triangles.push_back(*it);

            if ((M_comm.rank() == it->partition))
            {
                triangles_num_without_ghost.push_back(it->number);
            }
        }
    }

    M_num_triangles = M_triangles.size();

    // check the nodal partitions
    int elt_size = triangles_num_without_ghost.size();
    int num_elements = boost::mpi::all_reduce(M_comm, elt_size, std::plus<int>());

    if (M_global_num_elements_from_serial != num_elements)
    {
        LOG(DEBUG)<<"---------------------------------------Post-processing needed for element mesh partition: "<< M_global_num_elements_from_serial <<" != "<< num_elements <<"\n";
    }


    // move renumbering of triangles here (previously at the end of this function)
    // --------------------------------BEGINNING-------------------------
    int num_trls;

    std::vector<int> diff_trs;

    if (M_global_num_elements_from_serial != num_elements)
    {
        renumbering.resize(0);
        this->allGather(triangles_num_without_ghost, renumbering, num_trls);

        cpts = 0;
        cpts_dom = 0;

        std::vector<int> all_trs;

        for (int ii=0; ii<M_comm.size(); ++ii)
        {
            for (int jj=0; jj<renumbering[ii].size(); ++jj)
            {
                all_trs.push_back(renumbering[ii][jj]);
                ++cpts;
            }

            cpts_dom += renumbering[ii].size();
        }

        std::sort(all_trs.begin(), all_trs.end());

        std::vector<int> global_trs(all_trs.size());
        std::iota(global_trs.begin(), global_trs.end(), 1);

        std::set_difference(global_trs.begin(), global_trs.end(),
                            all_trs.begin(), all_trs.end(),
                            std::back_inserter(diff_trs));

        LOG(DEBUG)<<"---------------------------------------MISSING ELEMENTS= \n";

        for (int i=0; i<diff_trs.size(); ++i)
        {
            LOG(DEBUG)<<"                                                         ---IDS["<< i <<"]= "<< diff_trs[i] <<"\n";
        }
    }
    // --------------------------------END-------------------------------


    // --------------------------------BEGINNING-------------------------

    // elements in partition first and ghost at the end
    _triangles = M_triangles;
    M_triangles.resize(0);
    M_num_triangles_without_ghost = 0;

    for (auto it=_triangles.begin(), end=_triangles.end(); it!=end; ++it)
    {
        // treatment of missing elements
        for (int i=0; i<diff_trs.size(); ++i)
        {
            if (it->number == diff_trs[i])
            {
                auto ghosts_ = it->ghosts;
                int min_rank = *std::min_element(ghosts_.begin(), ghosts_.end());

                if ((M_comm.rank() == min_rank))
                {
                    triangles_num_without_ghost.push_back(it->number);
                    it->partition = M_comm.rank();
                }
            }
        }

        if (M_comm.rank() == it->partition)
        {
            M_triangles.push_back(*it);
            M_triangles_id_with_ghost.push_back(it->number);
            ++M_num_triangles_without_ghost;
        }
    }

    for (auto it=_triangles.begin(), end=_triangles.end(); it!=end; ++it)
    {
        if (M_comm.rank() != it->partition)
        {
            M_triangles.push_back(*it);
            M_triangles_id_with_ghost.push_back(it->number);
        }
    }

    for (int k=0; k<M_triangles_id_with_ghost.size(); ++k)
    {
        // mapping from old global numbering to local numbering
        M_transfer_map_elt.insert(position(M_triangles_id_with_ghost[k],k+1));
    }

    // --------------------------------END-------------------------------

    // reorder edge nodes
    for (auto it=M_edges.begin(), end=M_edges.end(); it!=end; ++it)
    {
        it->ghostNodes.assign(2,false);

        for (int i=0; i<2; ++i)
        {
            int rdof = reorder[it->indices[i]];

            it->indices[i] = M_transfer_map.left.find(it->indices[i])->second;

            if (std::binary_search(M_local_ghost.begin(),M_local_ghost.end(),rdof))
            {
                it->ghostNodes[i] = true;
            }
        }
    }

    // --------------------------------BEGINNING-------------------------
    std::sort(triangles_num_without_ghost.begin(), triangles_num_without_ghost.end());
    renumbering.resize(0);
    //int num_trls;
    allGather(triangles_num_without_ghost, renumbering, num_trls);

    cpts = 0;
    size = 0;
    for (int ii=0; ii<M_comm.size(); ++ii) size += renumbering[ii].size();
    M_global_num_elements = size;
    M_map_elements.resize(M_global_num_elements);

    for (int ii=0; ii<M_comm.size(); ++ii)
    {
        for (int jj=0; jj<renumbering[ii].size(); ++jj)
        {
            M_map_elements[renumbering[ii][jj]-1] = cpts++;
        }
    }
    // --------------------------------END-------------------------------

}//nodalGrid

void
GmshMesh::allGather(std::vector<int> const& field_in, std::vector<std::vector<int> >& field_out, int& acc_size)
{
    int fd_size = field_in.size();

    std::vector<int> container_size;
    boost::mpi::all_gather(M_comm, fd_size, container_size);
    int num_elts = std::accumulate(container_size.begin(),container_size.end(),0);
    acc_size = num_elts;

    std::vector<int> field_gather(num_elts);

    std::vector<int> displs(M_comm.size(), 0);
    for (int k = 1; k < M_comm.size(); ++k) {
        displs[k] = displs[k - 1] + container_size[k - 1];
    }

    int ier = MPI_Allgatherv(&field_in[0], fd_size, MPI_INT, &field_gather[0], &container_size[0], &displs[0], MPI_INT, MPI_Comm(M_comm));

    field_out.resize(M_comm.size());

    int global_indexing = 0;
    for (int ii=0; ii<M_comm.size(); ++ii)
    {
        int current_size = container_size[ii];
        field_out[ii].resize(current_size);

        for (int jj=0; jj<current_size; ++jj)
        {
            field_out[ii][jj] = field_gather[global_indexing+jj];
        }
        global_indexing += current_size;
    }
}


std::vector<int>
GmshMesh::indexTr() const
{
    std::vector<int> index;
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        for (int i=0; i<3; ++i)
        {
            index.push_back(it->indices[i]);
        }
    }

    return index;
}

std::vector<double>
GmshMesh::coordX() const
{
    std::vector<double> x(M_num_nodes);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        x[cpt] = it->second.coords[0];
        ++cpt;
    }

    return x;
}

std::vector<double>
GmshMesh::coordY() const
{
    std::vector<double> y(M_num_nodes);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        y[cpt] = it->second.coords[1];
        ++cpt;
    }

    return y;
}

std::vector<double>
GmshMesh::coordX(double const& rotangle) const
{
    std::vector<double> x(M_num_nodes);
    int cpt = 0;
    double cos_rotangle = std::cos(rotangle);
    double sin_rotangle = std::sin(rotangle);
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        x[cpt] = cos_rotangle*(it->second.coords[0]) + sin_rotangle*(it->second.coords[1]);
        ++cpt;
    }

    return x;
}

std::vector<double>
GmshMesh::coordY(double const& rotangle) const
{
    std::vector<double> y(M_num_nodes);
    int cpt = 0;
    double cos_rotangle=std::cos(rotangle);
    double sin_rotangle=std::sin(rotangle);
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        y[cpt] = -sin_rotangle*(it->second.coords[0]) + cos_rotangle*(it->second.coords[1]);
        ++cpt;
    }

    return y;
}

std::vector<double>
GmshMesh::bCoordX() const
{
    std::vector<double> node(M_num_nodes);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        node[cpt] = it->second.coords[0];
        ++cpt;
    }

    std::vector<double> bcoord_x(M_num_triangles);
    cpt = 0;
    double x = 0.;
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        x = 0.;

        for (int i=0; i<3; ++i)
        {
            x += node[it->indices[i]-1];
        }

        bcoord_x[cpt] = x/3.;

        ++cpt;
    }

    return bcoord_x;
}

std::vector<double>
GmshMesh::bCoordY() const
{
    std::vector<double> node(M_num_nodes);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        node[cpt] = it->second.coords[1];
        ++cpt;
    }

    std::vector<double> bcoord_y(M_num_triangles);
    cpt = 0;
    double y = 0.;
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        y = 0.;

        for (int i=0; i<3; ++i)
        {
            y += node[it->indices[i]-1];
        }

        bcoord_y[cpt] = y/3.;

        ++cpt;
    }

    return bcoord_y;
}

std::vector<double>
GmshMesh::bCoordX(double const& rotangle) const
{
    std::vector<double> bcoord_x(M_num_triangles);
    double cos_rotangle=std::cos(rotangle);
    double sin_rotangle=std::sin(rotangle);
    int cpt = 0;
    double x = 0.;
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        x = 0.;

        for (int i=0; i<3; ++i)
        {
            x += cos_rotangle*(M_nodes.find(it->indices[i])->second.coords[0]) + sin_rotangle*(M_nodes.find(it->indices[i])->second.coords[1]);
        }

        bcoord_x[cpt] = x/3.;

        ++cpt;
    }

    return bcoord_x;
}

std::vector<double>
GmshMesh::bCoordY(double const& rotangle) const
{
    std::vector<double> bcoord_y(M_num_triangles);
    double cos_rotangle=std::cos(rotangle);
    double sin_rotangle=std::sin(rotangle);
    int cpt = 0;
    double y = 0.;
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        y = 0.;

        for (int i=0; i<3; ++i)
        {
            y += -sin_rotangle*(M_nodes.find(it->indices[i])->second.coords[0]) + cos_rotangle*(M_nodes.find(it->indices[i])->second.coords[1]);
        }

        bcoord_y[cpt] = y/3.;

        ++cpt;
    }

    return bcoord_y;
}

std::vector<double>
GmshMesh::meanLon() const
{
    mapx_class *map;
    std::vector<char> str(M_mppfile.begin(), M_mppfile.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    std::vector<double> mean_lon(M_num_triangles);
    double lat = 0.;
    double lon = 0.;

    std::vector<double> X = this->bCoordX();
    std::vector<double> Y = this->bCoordY();

    for (int elt=0; elt<M_num_triangles; ++elt)
    {
        int status = inverse_mapx(map,X[elt],Y[elt],&lat,&lon);
        mean_lon[elt] = lon;
    }

    close_mapx(map);

    return mean_lon;
}

std::vector<double>
GmshMesh::meanLat() const
{
    mapx_class *map;
    std::vector<char> str(M_mppfile.begin(), M_mppfile.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    std::vector<double> mean_lat(M_num_triangles);
    double lat = 0.;
    double lon = 0.;

    std::vector<double> X = this->bCoordX();
    std::vector<double> Y = this->bCoordY();

    for (int elt=0; elt<M_num_triangles; ++elt)
    {
        int status = inverse_mapx(map,X[elt],Y[elt],&lat,&lon);
        mean_lat[elt] = lat;
    }

    close_mapx(map);

    return mean_lat;
}

std::vector<double>
GmshMesh::lon() const
{
    mapx_class *map;
    std::vector<char> str(M_mppfile.begin(), M_mppfile.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    std::vector<double> node_lon(M_num_triangles);
    double lat = 0.;
    double lon = 0.;

    std::vector<double> X = this->coordX();
    std::vector<double> Y = this->coordY();

    for (int nod=0; nod<M_num_nodes; ++nod)
    {
        int status = inverse_mapx(map,X[nod],Y[nod],&lat,&lon);
        node_lon[nod] = lon;
    }

    close_mapx(map);

    return node_lon;
}

std::vector<double>
GmshMesh::lat() const
{
    mapx_class *map;
    std::vector<char> str(M_mppfile.begin(), M_mppfile.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    std::vector<double> node_lat(M_num_triangles);
    double lat = 0.;
    double lon = 0.;

    std::vector<double> X = this->coordX();
    std::vector<double> Y = this->coordY();

    for (int nod=0; nod<M_num_nodes; ++nod)
    {
        int status = inverse_mapx(map,X[nod],Y[nod],&lat,&lon);
        node_lat[nod] = lat;
    }

    close_mapx(map);

    return node_lat;
}

std::vector<int>
GmshMesh::indexTrPartition() const
{
    std::vector<int> index;
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        if (!it->is_ghost)
        {
            for (int i=0; i<3; ++i)
            {
                index.push_back(it->indices[i]);
            }
        }
    }

    return index;
}

std::vector<double>
GmshMesh::coordXPartition() const
{
    std::vector<double> x(M_nldof_without_ghost);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        if (cpt < M_nldof_without_ghost)
        {
            x[cpt] = it->second.coords[0];
        }
        else
        {
            break;
        }

        ++cpt;
    }

    return x;
}

std::vector<double>
GmshMesh::coordYPartition() const
{
    std::vector<double> y(M_nldof_without_ghost);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        if (cpt < M_nldof_without_ghost)
        {
            y[cpt] = it->second.coords[1];
        }
        else
        {
            break;
        }

        ++cpt;
    }

    return y;
}

void
GmshMesh::setId(std::vector<int> const& newid)
{
    if ( newid.size() != 0 )
    {
        ASSERT(M_nodes.size()==newid.size(),"invalid size of new_id vector");

        for (int i=0; i<M_nodes.size(); ++i)
            M_nodes[i].id = newid[i];
    }
}

std::vector<int>
GmshMesh::id() const
{
    std::vector<int> mesh_id(M_num_nodes);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        mesh_id[cpt] = it->second.id;
        ++cpt;
    }

    return mesh_id;
}


// ------------------------------------------------
//! return the vertices for a given list of indices
//! called by FiniteElement::shapeCoeff() and FiniteElement::jacobian()
std::vector<std::vector<double>>
GmshMesh::vertices(std::vector<int> const& indices) const
{
    int const nv = indices.size();
    std::vector<std::vector<double>> vertices(nv);
    for(int i=0; i<nv; i++)
        vertices[i] = M_nodes[indices[i]].coords;
    return vertices;
}//vertices


std::vector<std::vector<double>>
GmshMesh::vertices(std::vector<int> const& indices,
        std::vector<double> const& um, double factor) const
{
    int const nv = indices.size();
    auto vertices = this->vertices(indices);
    for(int i=0; i<nv; i++)
        for(int k=0; k<2; k++)
            vertices[i][k] += factor*um[indices[i]-1+k*M_num_nodes];
    return vertices;
}//vertices


void GmshMesh::SwapBytes(void* array, size_t size, size_t n) {
    unsigned char* p = static_cast<unsigned char*>(array);
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < size/2; ++i) {
            std::swap(p[i], p[size-1-i]);
        }
        p += size;
    }
}// SwapBytes


} // Nextsim
