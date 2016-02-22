/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   gmshmesh.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Wed Jul 29 15:26:47 2015
 */

#include <gmshmesh.hpp>

namespace Nextsim
{
GmshMesh::GmshMesh(Communicator const& comm)
    :
    M_comm(comm),
    M_version("2.2"),
    M_ordering("gmsh"),
    M_nodes(),
    //M_elements(),
    M_triangles(),
    M_edges(),
    M_num_nodes(0),
    //M_num_elements(0),
    M_num_triangles(0),
    M_num_edges(0)
{}

GmshMesh::GmshMesh(std::vector<point_type> const& nodes,
                   std::vector<element_type> const& edges,
                   std::vector<element_type> const& triangles,
                   Communicator const& comm)
    :
    M_comm(comm),
    M_version("2.2"),
    M_ordering("gmsh"),
    M_nodes(nodes),
    M_triangles(triangles),
    M_edges(edges),
    M_num_nodes(nodes.size()),
    M_num_triangles(triangles.size()),
    M_num_edges(edges.size())
{}

void
GmshMesh::readFromFile(std::string const& filename)
{
    std::string gmshmshfile = Environment::nextsimDir().string() + "/mesh/" + filename;
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

    //std::map<int, Nextsim::entities::GMSHPoint > gmshpts;
    std::cout << "Reading "<< __n << " nodes\n";

    M_nodes.resize(__n);
    std::vector<double> coords(3,0);

    for ( unsigned int __i = 0; __i < __n; ++__i )
    {
        int id = 0;

        __is >> id
             >> coords[0]
             >> coords[1]
             >> coords[2];

        M_nodes[id-1].id = id;
        M_nodes[id-1].coords = coords;
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

    std::cout << "Reading " << numElements << " elements...\n";
    //std::list<Nextsim::entities::GMSHElement> __et; // tags in each element
    std::map<int,int> __gt;

    int cpt_edge = 0;
    int cpt_triangle = 0;

    for(int i = 0; i < numElements; i++)
    {
        int number, type, physical = 0, elementary = 0, numVertices;
        std::vector<int> ghosts;
        int numTags;
        int partition = (this->comm().size()>1)?this->comm().rank():0;

        __is >> number  // elm-number
             >> type // elm-type
             >> numTags; // number-of-tags

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

        int cpt_elt = (type == 2) ? cpt_triangle : cpt_edge;

        //std::cout<<"On proc "<< this->comm().rank() <<" : Global size= "<< this->comm().size() <<"\n";

        Nextsim::entities::GMSHElement gmshElt( number /*cpt_elt*/,
                                                type,
                                                physical,
                                                elementary,
                                                numPartitions,
                                                partition,
                                                ghosts,
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

    for ( auto const& it : __gt )
    {
        const char* name;
        MElement::getInfoMSH( it.first, &name );
        std::cout << "Read " << it.second << " " << name << " elements\n";

        if (std::string(name) == "Triangle 3")
            M_num_triangles = it.second;
        else if (std::string(name) == "Line 2")
            M_num_edges = it.second;
    }

    // make sure that we have read everything
    __is >> __buf;

    ASSERT(std::string( __buf ) == "$EndElements","invalid end elements string");

    // we are done reading the MSH file

    // create local dofs
    this->nodalGrid();
}

void
GmshMesh::writeTofile(std::string const& filename)
{
    //std::string gmshfilename = (boost::format( "../data/arctic10km.msh" ) ).str();
    std::string gmshmshfile = Environment::nextsimDir().string() + "/mesh/" + filename;
    std::fstream gmshfile(gmshmshfile, std::ios::out | std::ios::trunc);

    if (gmshfile.is_open())
    {
        gmshfile << "$MeshFormat\n";
        gmshfile << "2.2 0 8\n";
        gmshfile << "$EndMeshFormat\n";

        gmshfile << "$Nodes\n";
        gmshfile << M_num_nodes << "\n";

        //for ( int node = 0; node < M_num_nodes; node++ )
        int node = 0;
        for (auto it=M_nodes.begin(), en=M_nodes.end(); it!=en; ++it)
        {
            // gmshfile << node + 1
            //          << "  " << it->second.coords[0]
            //          << "  " << it->second.coords[1]
            //          << "  0.0\n";

            gmshfile << node + 1
                     << "  " << it->coords[0]
                     << "  " << it->coords[1]
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
            //<< "  " << element + 1;

            for (int i = 0; i < 3; i++ )
            {
                //gmshfile << "  " << it->second.indices[i];
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

        for (int i=0; i<M_nodes.size(); ++i)
        {
            //std::cout<<"ADDED= "<< factor*um[2*(element.indices[i]-1)] <<"\n";
            M_nodes[i].coords[0] += factor*um[i];
            M_nodes[i].coords[1] += factor*um[i+M_num_nodes];
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

    std::cout<<"MFILE= "<< std::string(map->mpp_filename) <<"\n";
    std::cout<<"PROJE= "<< std::string(map->projection_name) <<"\n";

    int cpt = 0;

    for (auto it=M_nodes.begin(), en=M_nodes.end(); it!=en; ++it)
    {
        //compute latitude and longitude from cartesian coordinates
        double _x = it->coords[0];
        double _y = it->coords[1];
        double _z = it->coords[2];

        // compute radius
        double radius = std::sqrt(std::pow(_x,2.)+std::pow(_y,2.)+std::pow(_z,2.));

        double latitude = std::asin(_z/radius)*(180./PI);
        double longitude = std::atan2(_y,_x);

        longitude = longitude-2*PI*std::floor(longitude/(2*PI));
        longitude = longitude*(180./PI);

        double x_, y_;
        int status = forward_mapx(map,latitude,longitude,&x_,&y_);

        it->coords[0] = x_;
        it->coords[1] = y_;
        it->coords[2] = 0.0;

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
    }

    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        if (!it->is_ghost)
        {
            for (int const& index : it->indices)
            {
                if ((std::find(M_local_dof_without_ghost.begin(),M_local_dof_without_ghost.end(),index) == M_local_dof_without_ghost.end())
                    && (std::find(ghosts_nodes_f.begin(),ghosts_nodes_f.end(),index) == ghosts_nodes_f.end()))
                {
                    M_local_dof_without_ghost.push_back(index);
                }


                if ((it->ghosts.size() > 0) && (std::find(ghosts_nodes_f.begin(),ghosts_nodes_f.end(),index) != ghosts_nodes_f.end()))
                {
                    int neigh = *std::min_element(it->ghosts.begin(),it->ghosts.end());
                    int min_id = std::min(neigh,it->partition);
                    int max_id = std::max(neigh,it->partition);

                    if (this->comm().rank() == min_id)
                    {
                        M_local_dof_without_ghost.push_back(index);
                    }
                    else if (this->comm().rank() == max_id)
                    {
                        ghosts_nodes_s.push_back(index);
                    }
                }
            }
        }
    }

    std::sort(ghosts_nodes_f.begin(),ghosts_nodes_f.end());

    std::sort(ghosts_nodes_s.begin(), ghosts_nodes_s.end());
    ghosts_nodes_s.erase(std::unique( ghosts_nodes_s.begin(), ghosts_nodes_s.end() ), ghosts_nodes_s.end());

    std::sort(M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end());
    M_local_dof_without_ghost.erase(std::unique( M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end() ), M_local_dof_without_ghost.end());


    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        if (((it->is_ghost) && (this->comm().rank() > it->partition)) && (it->ghosts.size() > 0) )
        {
            for (int const& index : it->indices)
            {
                if ((std::find(ghosts_nodes_t.begin(),ghosts_nodes_t.end(),index) == ghosts_nodes_t.end())
                    && (std::find(ghosts_nodes_s.begin(),ghosts_nodes_s.end(),index) == ghosts_nodes_s.end()))
                {
                    ghosts_nodes_t.push_back(index);
                }
            }
        }
    }

    std::sort(ghosts_nodes_t.begin(),ghosts_nodes_t.end());


    std::vector<int> diff_nodes;
    std::set_difference(ghosts_nodes_f.begin(), ghosts_nodes_f.end(),
                        M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.end(),
                        std::back_inserter(diff_nodes));

    std::set_difference(diff_nodes.begin(), diff_nodes.end(),
                        ghosts_nodes_t.begin(), ghosts_nodes_t.end(),
                        std::back_inserter(M_local_ghost));


    std::copy_n(M_local_dof_without_ghost.begin(), M_local_dof_without_ghost.size(), std::back_inserter(M_local_dof_with_ghost));
    std::copy_n(M_local_ghost.begin(), M_local_ghost.size(), std::back_inserter(M_local_dof_with_ghost));

    for (int k=0; k<M_local_dof_with_ghost.size(); ++k)
    {
        M_transfer_map.insert(position(M_local_dof_with_ghost[k],k+1));
    }

#if 1

    std::vector<element_type> _triangles = M_triangles;
    M_triangles.resize(0);

    for (auto it=_triangles.begin(), end=_triangles.end(); it!=end; ++it)
    {
        if ((M_comm.rank() <= it->partition))
        {
            bool _test = false;

            for (int i=0; i<3; ++i)
            {
                if (M_transfer_map.left.find(it->indices[i]) == M_transfer_map.left.end())
                {
                    _test = true;
                    break;
                }
            }

            if (_test)
                continue;

            M_triangles.push_back(*it);
        }
    }

    M_num_triangles = M_triangles.size();
#endif

}

std::vector<int>
GmshMesh::indexTr() const
{
    std::vector<int> index;//(3*M_num_triangles);
    //int cpt = 0;
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        //bool elt_valid = this->isValid(*it);
        //std::cout<<"--------INDEX "<< it->number <<" isvalid= "<< elt_valid <<"\n";
#if 0
        if ((M_comm.rank() <= it->partition)) //&& (elt_valid))
        {
            bool _test = false;

            for (int i=0; i<3; ++i)
            {
                if (M_transfer_map.left.find(it->indices[i]) == M_transfer_map.left.end())
                {
                    _test = true;
                    break;
                }
            }

            if (_test)
                continue;

            std::cout<<"--------INDEX "<< it->number <<"\n";
            for (int i=0; i<3; ++i)
            {
                index.push_back(M_transfer_map.left.find(it->indices[i])->second);
            }
        }
#endif

        for (int i=0; i<3; ++i)
        {
            index.push_back(M_transfer_map.left.find(it->indices[i])->second);
        }

        // index[3*cpt] = it->indices[0];//it->first;
        // index[3*cpt+1] = it->indices[1];
        // index[3*cpt+2] = it->indices[2];
        // ++cpt;
    }

    return index;
}

bool
GmshMesh::isValid(element_type const& elt) const
{
    if (M_comm.rank() == elt.partition)
        return true;

    bool _test = true;

    for (int i=0; i<3; ++i)
    {
        //_test = ((std::binary_search(M_local_ghost.begin(),M_local_ghost.end(),elt.indices[i])) && (_test));
        //_test = ((std::binary_search(M_local_ghost.begin(),M_local_ghost.end(),elt.indices[i])) && _test);
        _test = ((std::binary_search(M_local_dof_with_ghost.begin(),M_local_dof_with_ghost.end(),elt.indices[i])) && _test);
        //_test = ((M_transfer_map.left.find(elt.indices[i]) != M_transfer_map.left.end()) || _test);
        //index.push_back(M_transfer_map.left.find(it->indices[i])->second);

        //std::cout<<"CoordIDX= "<< elt.indices[i] << " test= "<< _test <<"\n";
    }

    return _test;
}


std::vector<double>
GmshMesh::coordX() const
{
    //std::vector<double> x(M_num_nodes);
    std::vector<double> x(M_local_dof_with_ghost.size());
    //int cpt = 0;
    //for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    for (int cpt=0; cpt<M_local_dof_with_ghost.size(); ++cpt)
    {
        x[cpt] = M_nodes[M_local_dof_with_ghost[cpt]-1].coords[0];
        //++cpt;
    }

    return x;
}

std::vector<double>
GmshMesh::coordY() const
{
    //std::vector<double> y(M_num_nodes);
    std::vector<double> y(M_local_dof_with_ghost.size());
    //int cpt = 0;
    //for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    for (int cpt=0; cpt<M_local_dof_with_ghost.size(); ++cpt)
    {
        y[cpt] = M_nodes[M_local_dof_with_ghost[cpt]-1].coords[1];
        //++cpt;
    }

    return y;
}

std::vector<double>
GmshMesh::coordX(double const& rotangle) const
{
    std::vector<double> x(M_num_nodes);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        x[cpt] = std::cos(rotangle)*(it->coords[0]) + std::sin(rotangle)*(it->coords[1]);
        ++cpt;
    }

    return x;
}

std::vector<double>
GmshMesh::coordY(double const& rotangle) const
{
    std::vector<double> y(M_num_nodes);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        y[cpt] = -std::sin(rotangle)*(it->coords[0]) + std::cos(rotangle)*(it->coords[1]);
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
            x += M_nodes[it->indices[i]-1].coords[0];
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
            //y += M_nodes[it->indices[i]-1].coords[1];
            y += M_nodes[it->indices[i]-1].coords[1];
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


} // Nextsim
