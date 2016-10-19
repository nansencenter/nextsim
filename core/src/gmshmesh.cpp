/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gmshmesh.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Wed Jul 29 15:26:47 2015
 */

#include <gmshmesh.hpp>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <vector>
#include <string>
#include <iterator>
#include <algorithm>
#include <iostream>

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
    M_reorder_map_nodes(),
    M_reorder_map_elements()
{}

GmshMesh::GmshMesh(GmshMesh const& mesh)
    :
    M_comm(mesh.M_comm),
    M_version(mesh.M_version),
    M_ordering(mesh.M_ordering),
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
    M_reorder_map_nodes(mesh.M_reorder_map_nodes),
    M_reorder_map_elements(mesh.M_reorder_map_elements)
{}

#if 0
GmshMesh::GmshMesh(std::vector<point_type> const& nodes,
                   std::vector<element_type> const& edges,
                   std::vector<element_type> const& triangles,
                   Communicator const& comm)
    :
    M_comm(comm),
    M_version("2.2"),
    M_ordering("gmsh"),
    M_nodes_vec(nodes),
    M_triangles(triangles),
    M_edges(edges),
    M_num_nodes(nodes.size()),
    M_num_triangles(triangles.size()),
    M_num_edges(edges.size())
{}
#endif

void
GmshMesh::readFromFile(std::string const& filename)
{
    std::string gmshmshfile = Environment::nextsimDir().string() + "/mesh/" + filename;

    if (M_comm.rank() == 0)
        std::cout<<"Reading Msh file "<< gmshmshfile <<"\n";

    std::ifstream __is ( gmshmshfile.c_str() );

    if ( !__is.is_open() )
    {
        std::ostringstream ostr;
        std::cout << "Invalid file name " << gmshmshfile << " (file not found)\n";
        ostr << "Invalid file name " << gmshmshfile << " (file not found)\n";
        throw std::invalid_argument( ostr.str() );
    }

    char __buf[256];
    __is >> __buf;

    std::string theversion;

    double version = 2.2;

    if (std::string( __buf ) == "$MeshFormat")
    {
        int format, size;
        __is >> theversion >> format >> size;
        std::cout << "GMSH mesh file version : " << theversion << " format: " << (format?"binary":"ascii") << \
            " size of double: " << size << "\n";

        ASSERT(boost::lexical_cast<double>( theversion ) >= 2, "Nextsim supports only Gmsh version >= 2");

        version = boost::lexical_cast<double>( theversion );

        __is >> __buf;

        ASSERT(std::string( __buf ) == "$EndMeshFormat","invalid file format entry");

        __is >> __buf;

        std::cout << "[importergmsh] " << __buf << " (expect $PhysicalNames)\n";

    }

    // Read NODES

    //std::cout << "buf: "<< __buf << "\n";

    if ( !( std::string( __buf ) == "$NOD" ||
            std::string( __buf ) == "$Nodes" ||
            std::string( __buf ) == "$ParametricNodes") )
    {
        std::cout<< "invalid nodes string '" << __buf << "' in gmsh importer. It should be either $Nodes.\n";
    }

    bool has_parametric_nodes = ( std::string( __buf ) == "$ParametricNodes" );
    unsigned int __n;
    __is >> __n;

    M_num_nodes = __n;

    M_global_num_nodes_from_serial = M_num_nodes;

    //std::map<int, Nextsim::entities::GMSHPoint > gmshpts;
    if (M_comm.rank() == 0)
        std::cout << "Reading "<< __n << " nodes\n";

    M_nodes_vec.resize(__n);
    std::vector<double> coords(3,0);

    for ( unsigned int __i = 0; __i < __n; ++__i )
    {
        int id = 0;

        __is >> id
             >> coords[0]
             >> coords[1]
             >> coords[2];

        M_nodes_vec[id-1].id = id;
        M_nodes_vec[id-1].coords = coords;
    }

    __is >> __buf;
    //std::cout << "buf: "<< __buf << "\n";

    // make sure that we have read all the points

    ASSERT(std::string( __buf ) == "$EndNodes","invalid end nodes string");

    // Read ELEMENTS

    __is >> __buf;

    ASSERT(std::string( __buf ) == "$Elements","invalid elements string");

    int numElements;
    __is >> numElements;

    //M_num_elements = numElements;
    M_global_num_elements_from_serial = numElements;
    //M_global_num_elements_from_serial = 0;

    if (M_comm.rank() == 0)
        std::cout << "Reading " << numElements << " elements...\n";
    //std::list<Nextsim::entities::GMSHElement> __et; // tags in each element
    std::map<int,int> __gt;

    int cpt_edge = 0;
    int cpt_triangle = 0;
    int num_edge = 0;

    for(int i = 0; i < numElements; i++)
    {
        int number, type, physical = 0, elementary = 0, numVertices;
        std::vector<int> ghosts;
        std::vector<bool> ghostNodes;
        int numTags;
        int partition = (this->comm().size()>1)?this->comm().rank():0;

        __is >> number  // elm-number
             >> type // elm-type
             >> numTags; // number-of-tags

        if (type == 1)
        {
            // if current element is an edge, stop reading and go to next line
            __is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            num_edge = number;
            continue;
        }

        int numPartitions = 1;

        for(int j = 0; j < numTags; j++)
        {
            int tag;
            __is >> tag;
            if(j == 0) physical = tag;
            else if(j == 1) elementary = tag;
            else if((j == 2) && (numTags > 3)) numPartitions = tag;
            else if(j == 3) partition = tag-1;
            else if((j >= 4) && (j < 4 + numPartitions - 1)) ghosts.push_back((-tag)-1);
        }

        numVertices = MElement::getInfoMSH(type);
        ASSERT(numVertices!=0,"unknown number of vertices for element type");

        std::vector<int> indices(numVertices);
        for(int j = 0; j < numVertices; j++)
        {
            __is >> indices[j];
            // check
            //indices[j] = indices[j]-1;
        }

        if (M_ordering=="bamg")
        {
            std::next_permutation(indices.begin()+1,indices.end());
        }

        //int cpt_elt = (type == 2) ? cpt_triangle : cpt_edge;
        //if (type == 2)
        //std::cout<<"cpt_elt= "<< cpt_elt <<"\n";

        //std::cout<<"On proc "<< this->comm().rank() <<" : Global size= "<< this->comm().size() <<"\n";

        number = number - num_edge;

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


        //M_triangles.insert(std::make_pair(number,gmshElt));


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
        const char* name;
        MElement::getInfoMSH( it.first, &name );
        std::cout<<"["<< M_comm.rank() <<"] " << " Read " << it.second << " " << name << " elements\n";

        if (std::string(name) == "Triangle 3")
            M_num_triangles = it.second;
        else if (std::string(name) == "Line 2")
            M_num_edges = it.second;
    }

    // make sure that we have read everything
    __is >> __buf;

    ASSERT(std::string( __buf ) == "$EndElements","invalid end elements string");

    // we are done reading the MSH file

    //std::cout<<"nodalGrid starts\n";

    // create local dofs
    if (M_comm.size() > 1)
        this->nodalGrid();

    //std::cout<<"nodalGrid done\n";
}

void
GmshMesh::writeTofile(std::string const& filename)
{
    std::string gmshmshfile = Environment::nextsimDir().string() + "/mesh/" + filename;
    std::fstream gmshfile(gmshmshfile, std::ios::out | std::ios::trunc);

    if (gmshfile.is_open())
    {
        gmshfile << "$MeshFormat\n";
        gmshfile << "2.2 0 8\n";
        gmshfile << "$EndMeshFormat\n";

        gmshfile << "$Nodes\n";
        gmshfile << M_num_nodes << "\n";

        int node = 0;
        for (auto it=M_nodes.begin(), en=M_nodes.end(); it!=en; ++it)
        {
            gmshfile << node + 1
                     << "  " << it->second.coords[0]
                     << "  " << it->second.coords[1]
                     << "  0.0\n";

            ++node;
        }
        gmshfile << "$EndNodes\n";


        int element_type = 2;
        int tag_num = 2;
        int tag1 = 1;
        int tag2 = 0;

        gmshfile << "$Elements\n";
        gmshfile << M_num_triangles << "\n";

        //for ( int element = 0; element < nels; element++ )
        int element = 0;
        for (auto it=M_triangles.begin(), en=M_triangles.end(); it!=en; ++it)
        {
            //tag2 = element +1;

            gmshfile << element + 1
                     << "  " << element_type
                     << "  " << tag_num
                     << "  " << tag1
                     << "  " << tag2;

            for (int i = 0; i < 3; i++ )
            {
                gmshfile << "  " << it->indices[i];
            }
            gmshfile << "\n";

            ++element;
        }
        gmshfile << "$EndElements\n";


    }
    else
    {
        std::cout << "Cannot open " << gmshmshfile  << "\n";
        std::cerr << "error: open file " << gmshmshfile << " for output failed!" <<"\n";
        std::abort();
    }
}

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
GmshMesh::stereographicProjection()
{
    // polar stereographic projection
    mapx_class *map;
    std::string filename = Environment::nextsimDir().string() + "/data/NpsNextsim.mpp";
    std::vector<char> str(filename.begin(), filename.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    // std::cout<<"MFILE= "<< std::string(map->mpp_filename) <<"\n";
    // std::cout<<"PROJE= "<< std::string(map->projection_name) <<"\n";

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
    std::vector<int> ghosts_nodes_f;
    std::vector<int> ghosts_nodes_s;
    std::vector<int> ghosts_nodes_t;

    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        if (it->is_ghost)
        {
            for (int const& index : it->indices)
            {
                if ((std::find(ghosts_nodes_f.begin(),ghosts_nodes_f.end(),index) == ghosts_nodes_f.end()))
                {
                    ghosts_nodes_f.push_back(index);
                }
            }
        }

#if 0
        auto it_indices = it->indices;
        if ((std::find(it_indices.begin(),it_indices.end(),5730) != it_indices.end()))
        {
            std::cout<<"---------------------------------------------------\n";
            std::cout<<"it->rank                = "<< M_comm.rank() <<"\n";
            std::cout<<"it->number              = "<< it->number <<"\n";
            std::cout<<"it->type                = "<< it->type <<"\n";
            std::cout<<"it->physical            = "<< it->physical <<"\n";
            std::cout<<"it->elementary          = "<< it->elementary <<"\n";
            std::cout<<"it->numPartitions       = "<< it->numPartitions <<"\n";
            std::cout<<"it->partition           = "<< it->partition <<"\n";
            std::cout<<"it->is_ghost            = "<< it->is_ghost <<"\n";
            std::cout<<"ghost_partition_id      = "<< it->ghost_partition_id <<"\n";

            for (int i=0; i<it->indices.size(); ++i)
                std::cout<<"                          it->indices["<< i <<"]= "<< it->indices[i] <<"\n";

            std::cout<<"it->ghosts              = "<<"\n";
            for (int i=0; i<it->ghosts.size(); ++i)
                std::cout<<"                          it->ghosts["<< i <<"]= "<< it->ghosts[i] <<"\n";
        }
#endif
    }


    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        if (!it->is_ghost)
        {
            for (int const& index : it->indices)
            {
                if ((std::find(ghosts_nodes_f.begin(),ghosts_nodes_f.end(),index) == ghosts_nodes_f.end()))
                {
                    M_local_dof_without_ghost.push_back(index);
                }

                if ((it->ghosts.size() > 0) && (std::find(ghosts_nodes_f.begin(),ghosts_nodes_f.end(),index) != ghosts_nodes_f.end()))
                {
#if 0
                    int neigh = *std::min_element(it->ghosts.begin(),it->ghosts.end());
                    int min_id = std::min(neigh,it->partition);
                    //int max_id = std::max(neigh,it->partition);

                    if (this->comm().rank() == min_id)
                    {
                        //M_local_dof_without_ghost.push_back(index);
                    }
#endif
                    M_local_dof_without_ghost.push_back(index);
                }
            }
        }
    }

    std::sort(ghosts_nodes_f.begin(),ghosts_nodes_f.end());

    std::sort(M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end());
    M_local_dof_without_ghost.erase(std::unique( M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end() ), M_local_dof_without_ghost.end());


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
    std::sort(all_local_nodes.begin(), all_local_nodes.end());

    std::sort(M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end());

    std::set_difference(all_local_nodes.begin(), all_local_nodes.end(),
                        M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end(),
                        std::back_inserter(M_local_ghost));


    // check the nodal partitions
    int nd_size = M_local_dof_without_ghost.size();
    int num_nodes = boost::mpi::all_reduce(M_comm, nd_size, std::plus<int>());

    if (M_num_nodes != num_nodes)
    {
        std::cout<<"---------------------------------------Post-processing needed for nodal mesh partition: "<< M_num_nodes <<" != "<< num_nodes <<"\n";
    }

    // gather operations for global node renumbering
    std::vector<std::vector<int>> renumbering;
    boost::mpi::all_gather(M_comm,
                           M_local_dof_without_ghost,
                           renumbering);

    // check the nodal partition starts
    if (M_num_nodes != num_nodes)
    {
        if (M_num_nodes < num_nodes)
        {
            for (int ii=0; ii<renumbering.size(); ++ii)
            {
                for (int jj=0; jj<renumbering.size(); ++jj)
                {
                    if (ii != jj)
                    {
                        std::vector<int> duplicated_dofs;

                        std::set_intersection(renumbering[ii].begin(),renumbering[ii].end(),
                                              renumbering[jj].begin(),renumbering[jj].end(),
                                              std::back_inserter(duplicated_dofs));

                        if (duplicated_dofs.size() == 0)
                            continue;

                        for (int kk=0; kk<duplicated_dofs.size(); ++kk)
                        {
                            if (jj < ii)
                            {
                                renumbering[ii].erase(std::remove(renumbering[ii].begin(), renumbering[ii].end(), duplicated_dofs[kk]),
                                                      renumbering[ii].end());

                                if (M_comm.rank() == ii)
                                {
                                    M_local_ghost.push_back(duplicated_dofs[kk]);
                                }
                            }
                        }
                    }
                }
            }

            M_local_dof_without_ghost = renumbering[M_comm.rank()];
        }
#if 0
        else
        {
            std::sort(M_local_ghost.begin(), M_local_ghost.end());
            std::vector<std::vector<int>> global_ghosts;
            // gather operations for global ghost numbering
            boost::mpi::all_gather(M_comm,
                                   M_local_ghost,
                                   global_ghosts);


            std::vector<int> global_from_mesh(M_num_nodes);
            std::iota(global_from_mesh.begin(), global_from_mesh.end(), 1);

            std::vector<int> global_from_code;//(num_nodes);

            for (int ii=0; ii<renumbering.size(); ++ii)
            {
                for (int jj=0; jj<renumbering[ii].size(); ++jj)
                {
                    global_from_code.push_back(renumbering[ii][jj]);
                }
            }

            std::sort(global_from_code.begin(), global_from_code.end());

            std::vector<int> diff_trs;
            std::set_difference(global_from_mesh.begin(), global_from_mesh.end(),
                                global_from_code.begin(), global_from_code.end(),
                                std::back_inserter(diff_trs));

            for (int kk=0; kk<diff_trs.size(); ++kk)
            {
                //std::cout<<
                std::cout<<"---------------------------------------MISSING NODES= \n";
                std::cout<<"                                                         ---[" << M_comm.rank() << "]: IDS["<< kk <<"]= "<< diff_trs[kk] <<"\n";
            }

            for (int kk=0; kk<diff_trs.size(); ++kk)
            {
                std::vector<int> dom_ids;

                for (int ii=0; ii<global_ghosts.size(); ++ii)
                {
                    //global_from_code.push_back(renumbering[ii][jj]);
                    //if ((std::find(global_ghosts[ii].begin(),global_ghosts[ii].end(),diff_trs[kk]) != global_ghosts[ii].end()))
                    if (std::binary_search(global_ghosts[ii].begin(), global_ghosts[ii].end(),diff_trs[kk]))
                    {
                        dom_ids.push_back(ii);
                    }
                }

                int valid_id = *std::min_element(dom_ids.begin(), dom_ids.end());

                if (M_comm.rank() == valid_id)
                {
                    std::cout<<"----------------------------------WORKING["<< M_comm.rank() <<"]: "<< diff_trs[kk] <<"\n";

                    M_local_ghost.erase(std::remove(M_local_ghost.begin(), M_local_ghost.end(), diff_trs[kk]),
                                        M_local_ghost.end());

                    M_local_dof_without_ghost.push_back(diff_trs[kk]);
                }
            }
        }
#endif
    } // check the nodal partition end

    std::sort(M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end());
    std::sort(M_local_ghost.begin(), M_local_ghost.end());

    std::copy_n(M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.size(), std::back_inserter(M_local_dof_with_ghost));
    std::copy_n(M_local_ghost.begin(), M_local_ghost.size(), std::back_inserter(M_local_dof_with_ghost));


    //bimap_type reorder;
    std::map<int,int> reorder;

    int cpts = 0;
    int cpts_dom = 0;

    for (int ii=0; ii<M_comm.size(); ++ii)
    {
        int sr = renumbering[ii].size();

        for (int jj=0; jj<renumbering[ii].size(); ++jj)
        {
            //reorder.insert(position(renumbering[ii][jj],cpts+1+cpts_dom));
            reorder[renumbering[ii][jj]] = cpts+1+cpts_dom;

            //reorder_root[renumbering[ii][jj]] = cpts+1;
            M_reorder_map_nodes.insert(position(renumbering[ii][jj],cpts+1));

            // add second component for velocity
            reorder[renumbering[ii][jj]+M_num_nodes] = cpts+1+sr+cpts_dom;

#if 0
            if (M_comm.rank() == 0)
            {
                //std::cout<<"TEST MAPPING["<< cpts+1 <<"]= "<< renumbering[ii][jj] <<"\n";
                //std::cout<<"MAPPING["<< cpts+cpts_dom+1 <<"]= "<< renumbering[ii][jj] <<"\n";
                //std::cout<<"MAPPING["<< cpts+sr+cpts_dom+1 <<"]= "<< renumbering[ii][jj]+M_num_nodes <<"\n";
                //M_nodes_root[cpts] = M_nodes_vec[renumbering[ii][jj]-1];
            }
#endif

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

        int rdof = reorder.find(local_dof_with_ghost[k])->second;
        int rdofv = reorder.find(local_dof_with_ghost[k]+M_num_nodes)->second;

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
            //M_local_ghost[k-M_nldof_without_ghost+M_nlghost] = rdofv;
        }
    }

    std::sort(M_local_ghost.begin(), M_local_ghost.end());
    M_global_num_nodes = M_num_nodes;
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
                int rdof = reorder.find(it->indices[i])->second;

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
        std::cout<<"---------------------------------------Post-processing needed for element mesh partition: "<< M_global_num_elements_from_serial <<" != "<< num_elements <<"\n";
    }


    // move renumbering of triangles here (previously at the end of this function)
    // --------------------------------BEGINNING-------------------------
    renumbering.resize(0);

    // gather operations for global element renumbering
    boost::mpi::all_gather(M_comm,
                           triangles_num_without_ghost,
                           renumbering);

    std::vector<int> diff_trs;

    if (M_global_num_elements_from_serial != num_elements)
    {
        cpts = 0;
        cpts_dom = 0;

        std::vector<int> all_trs;

        for (int ii=0; ii<M_comm.size(); ++ii)
        {
            for (int jj=0; jj<renumbering[ii].size(); ++jj)
            {
                //M_reorder_map_elements.insert(position(renumbering[ii][jj],cpts+1));
                all_trs.push_back(renumbering[ii][jj]);
                ++cpts;
            }

            cpts_dom += renumbering[ii].size();
        }

        std::sort(all_trs.begin(), all_trs.end());

        std::vector<int> global_trs(all_trs.size());//(M_reorder_map_elements.size());
        std::iota(global_trs.begin(), global_trs.end(), 1);

        std::set_difference(global_trs.begin(), global_trs.end(),
                            all_trs.begin(), all_trs.end(),
                            std::back_inserter(diff_trs));

        if (M_comm.rank() == 0)
        {
            std::cout<<"---------------------------------------MISSING ELEMENTS= \n";

            for (int i=0; i<diff_trs.size(); ++i)
            {
                std::cout<<"                                                         ---IDS["<< i <<"]= "<< diff_trs[i] <<"\n";
            }
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
            int rdof = reorder.find(it->indices[i])->second;

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
    // gather operations for global element renumbering
    boost::mpi::all_gather(M_comm,
                           triangles_num_without_ghost,
                           renumbering);

    cpts = 0;
    cpts_dom = 0;

    for (int ii=0; ii<M_comm.size(); ++ii)
    {
        for (int jj=0; jj<renumbering[ii].size(); ++jj)
        {
            M_reorder_map_elements.insert(position(renumbering[ii][jj],cpts+1));
            ++cpts;
        }
    }
    // --------------------------------END-------------------------------

    M_global_num_elements = M_reorder_map_elements.size();

    // std::cout<<"["<< this->comm().rank() << "] M_global_num_elements= "<< M_global_num_elements <<"\n";
    // std::cout<<"["<< this->comm().rank() << "] M_num_elements       = "<< triangles_num_without_ghost.size() <<"\n";
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
    std::vector<double> bcoord_x(M_num_triangles);
    int cpt = 0;
    double x = 0.;
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        x = 0.;

        for (int i=0; i<3; ++i)
        {
            x += M_nodes.find(it->indices[i])->second.coords[0];
        }

        bcoord_x[cpt] = x/3.;

        ++cpt;
    }

    return bcoord_x;
}

std::vector<double>
GmshMesh::bCoordY() const
{
    std::vector<double> bcoord_y(M_num_triangles);
    int cpt = 0;
    double y = 0.;
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        y = 0.;

        for (int i=0; i<3; ++i)
        {
            y += M_nodes.find(it->indices[i])->second.coords[1];
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
    std::string filename = Environment::nextsimDir().string() + "/data/NpsNextsim.mpp";
    std::vector<char> str(filename.begin(), filename.end());
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
    std::string filename = Environment::nextsimDir().string() + "/data/NpsNextsim.mpp";
    std::vector<char> str(filename.begin(), filename.end());
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
    std::string filename = Environment::nextsimDir().string() + "/data/NpsNextsim.mpp";
    std::vector<char> str(filename.begin(), filename.end());
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
    std::string filename = Environment::nextsimDir().string() + "/data/NpsNextsim.mpp";
    std::vector<char> str(filename.begin(), filename.end());
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

} // Nextsim
