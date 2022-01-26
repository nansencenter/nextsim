/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gmshmeshseq.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Tue Jun 21 14:25:38 2016
 */

#include <gmshmeshseq.hpp>
#include <GModel.h>
#include <GmshMessage.h>
#include <boost/format.hpp>

namespace Nextsim
{
GmshMeshSeq::GmshMeshSeq()
    :
    M_version("2.2"),
    M_ordering("gmsh"),
    M_nodes(),
    M_triangles(),
    M_edges(),
    M_num_nodes(0),
    M_num_triangles(0),
    M_num_edges(0),
    M_marker_names(),
    timer(),
    M_mppfile(Environment::nextsimMppfile()),
    M_log_level(Environment::logLevel()),
    M_log_all(Environment::logAll())

{}

GmshMeshSeq::GmshMeshSeq(std::vector<point_type> const& nodes,
                         std::vector<element_type> const& edges,
                         std::vector<element_type> const& triangles)
    :
    M_version("2.2"),
    M_ordering("gmsh"),
    M_nodes(nodes),
    M_triangles(triangles),
    M_edges(edges),
    M_num_nodes(nodes.size()),
    M_num_triangles(triangles.size()),
    M_num_edges(edges.size()),
    M_mppfile(Environment::nextsimMppfile()),
    M_log_level(Environment::logLevel()),
    M_log_all(Environment::logAll())

{}

GmshMeshSeq::GmshMeshSeq(std::vector<point_type> const& nodes,
                         std::vector<element_type> const& triangles)
    :
    M_version("2.2"),
    M_ordering("gmsh"),
    M_nodes(nodes),
    M_triangles(triangles),
    M_num_nodes(nodes.size()),
    M_num_triangles(triangles.size()),
    M_mppfile(Environment::nextsimMppfile()),
    M_log_level(Environment::logLevel()),
    M_log_all(Environment::logAll())

{}

GmshMeshSeq::GmshMeshSeq(GmshMeshSeq const& mesh)
        :
    M_version(mesh.M_version),
    M_ordering(mesh.M_ordering),
    M_mppfile(mesh.M_mppfile),
    M_nodes(mesh.M_nodes),
    M_triangles(mesh.M_triangles),
    M_num_nodes(mesh.M_num_nodes),
    M_num_triangles(mesh.M_num_triangles),
    M_log_level(mesh.M_log_level),
    M_log_all(mesh.M_log_all)

{}

// GmshMeshSeq::~GmshMeshSeq()
// {
//     //M_gmodel->deleteMesh();
//     //M_gmodel->destroy();
//     delete M_gmodel;
// }

void
GmshMeshSeq::readFromFile(std::string const& gmshmshfile)
{
    LOG(VERBOSE)<<"Reading Msh file "<< gmshmshfile <<"\n";

    std::ifstream gmshfile ( gmshmshfile.c_str() );

    if ( !gmshfile.is_open() )
    {
        std::ostringstream ostr;
        LOG(ERROR) << "Invalid file name " << gmshmshfile << " (file not found)\n";
        ostr << "Invalid file name " << gmshmshfile << " (file not found)\n";
        throw std::invalid_argument( ostr.str() );
    }

    char __buf[256];
    gmshfile >> __buf;

    std::string theversion;

    double version = 2.2;

    if (std::string( __buf ) == "$MeshFormat")
    {
        int format, size;
        gmshfile >> theversion >> format >> size;
        LOG(DEBUG) << "GMSH mesh file version : " << theversion << " format: " << (format?"binary":"ascii") << \
            " size of double: " << size << "\n";

        ASSERT(boost::lexical_cast<double>( theversion ) >= 2, "Nextsim supports only Gmsh version >= 2");

        version = boost::lexical_cast<double>( theversion );

        gmshfile >> __buf;

        ASSERT(std::string( __buf ) == "$EndMeshFormat","invalid file format entry");

        gmshfile >> __buf;

        LOG(DEBUG) << "[importergmsh] " << __buf << " (expect $PhysicalNames)\n";

        if ( std::string( __buf ) == "$PhysicalNames" )
        {
            int nnames;
            gmshfile >> nnames;

            for ( int n = 0; n < nnames; ++n )
            {
                int id, topodim;
                std::string name;

                gmshfile >> topodim >> id >> name;

                boost::trim( name );
                boost::trim_if( name,boost::is_any_of( "\"" ) );

                LOG(DEBUG) << "[gmshmesh::reading] topodim: "  << topodim << " id: " << id << " name: " << name << "\n";

                std::vector<int> marker_data = {id, topodim};
                M_marker_names.insert(std::make_pair(name,marker_data));
            }

            gmshfile >> __buf;
            ASSERT(std::string( __buf ) == "$EndPhysicalNames","invalid file format entry");

            gmshfile >> __buf;
        }
    }

    // Read NODES

    //std::cout << "buf: "<< __buf << "\n";

    if ( !( std::string( __buf ) == "$NOD" ||
            std::string( __buf ) == "$Nodes" ||
            std::string( __buf ) == "$ParametricNodes") )
    {
        LOG(ERROR)<< "invalid nodes string '" << __buf << "' in gmsh importer. It should be either" 
        << "$Nodes or $ParametricNodes.\n";
    }

    bool has_parametric_nodes = ( std::string( __buf ) == "$ParametricNodes" );
    unsigned int __n;
    gmshfile >> __n;

    M_num_nodes = __n;

    //std::map<int, Nextsim::entities::GMSHPoint > gmshpts;
    LOG(DEBUG) << "Reading "<< __n << " nodes\n";

    M_nodes.resize(__n);
    std::vector<double> coords(3,0);

    for ( unsigned int __i = 0; __i < __n; ++__i )
    {
        int id = 0;

        gmshfile >> id
             >> coords[0]
             >> coords[1]
             >> coords[2];

        M_nodes[id-1].id = id;
        M_nodes[id-1].coords = coords;
    }

    gmshfile >> __buf;
    //std::cout << "buf: "<< __buf << "\n";

    // make sure that we have read all the points

    ASSERT(std::string( __buf ) == "$EndNodes","invalid end nodes string");

    // Read ELEMENTS

    gmshfile >> __buf;

    ASSERT(std::string( __buf ) == "$Elements","invalid elements string");

    int numElements;
    gmshfile >> numElements;

    //M_num_elements = numElements;

    LOG(DEBUG) << "Reading " << numElements << " elements...\n";
    //std::list<Nextsim::entities::GMSHElement> __et; // tags in each element
    std::map<int,int> __gt;

    int cpt_edge = 0;
    int cpt_triangle = 0;

    for(int i = 0; i < numElements; i++)
    {
        int number, type, physical = 0, elementary = 0, numVertices;
        int numTags;

        gmshfile >> number  // elm-number
             >> type // elm-type
             >> numTags; // number-of-tags

        for(int j = 0; j < numTags; j++)
        {
            int tag;
            gmshfile >> tag;
            if(j == 0) physical = tag;
            else if(j == 1) elementary = tag;
        }

        numVertices = MElement::getInfoMSH(type);

        ASSERT(numVertices!=0,"unknown number of vertices for element type");

        std::vector<int> indices(numVertices);
        for(int j = 0; j < numVertices; j++)
        {
            gmshfile >> indices[j];
            // check
            //indices[j] = indices[j]-1;
        }

        if (M_ordering=="bamg")
        {
            std::next_permutation(indices.begin()+1,indices.end());
        }

        // Nextsim::entities::GMSHElement gmshElt( number,
        //                                         type,
        //                                         physical,
        //                                         elementary,
        //                                         numVertices,
        //                                         indices );

        //__et.push_back( gmshElt );
        //M_elements.insert(std::make_pair(number,gmshElt));

        if (type == 2)
        {
            Nextsim::entities::GMSHElement gmshElt( cpt_triangle,
                                                    type,
                                                    physical,
                                                    elementary,
                                                    numVertices,
                                                    indices );

            //M_triangles.insert(std::make_pair(number,gmshElt));
            M_triangles.push_back(gmshElt);

            ++cpt_triangle;
        }
        else if (type == 1)
        {
            Nextsim::entities::GMSHElement gmshElt( cpt_edge,
                                                    type,
                                                    physical,
                                                    elementary,
                                                    numVertices,
                                                    indices );

            //M_edges.insert(std::make_pair(number,gmshElt));
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
        LOG(DEBUG) << "Read " << it.second << " " << name << " elements\n";

        if (std::string(name) == "Triangle 3")
            M_num_triangles = it.second;
        else if (std::string(name) == "Line 2")
            M_num_edges = it.second;
    }

    // make sure that we have read everything
    gmshfile >> __buf;

    ASSERT(std::string( __buf ) == "$EndElements","invalid end elements string");

    // we are done reading the MSH file
}//readFromFile

void
GmshMeshSeq::writeToFile(std::string const& gmshmshfile)
{
    std::fstream gmshfile(gmshmshfile, std::ios::out | std::ios::trunc);

    if (!gmshfile.is_open())
    {
        LOG(ERROR) << "Cannot open " << gmshmshfile  << "\n";
        std::cerr << "error: open file " << gmshmshfile << " for output failed!" <<"\n";
        std::abort();
    }

    gmshfile << "$MeshFormat\n";
    gmshfile << "2.2 0 8\n";
    gmshfile << "$EndMeshFormat\n";

    gmshfile << "$Nodes\n";
    gmshfile << M_num_nodes << "\n";

    int node = 0;
    for (auto it=M_nodes.begin(), en=M_nodes.end(); it!=en; ++it)
    {
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

    int element = 0;
    for (auto it=M_triangles.begin(), en=M_triangles.end(); it!=en; ++it)
    {
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

}//writeToFile

void
GmshMeshSeq::update(std::vector<point_type> const& nodes,
                    std::vector<element_type> const& triangles,
                    std::string const& ordering)
{
    M_ordering = ordering;
    M_nodes = nodes;
    M_triangles = triangles;
    M_num_nodes = nodes.size();
    M_num_triangles = triangles.size();
}

void
GmshMeshSeq::initGModel()
{
    M_gmodel = new GModel();

    Msg::SetVerbosity(Environment::vm()["debugging.gmsh_verbose"].as<int>());
    CTX::instance()->terminal = 1;
    CTX::instance()->mesh.saveTopology = 0;
    CTX::instance()->mesh.fileFormat = FORMAT_MSH;
    CTX::instance()->mesh.mshFileVersion = 2.2;
    //M_partition_options.num_partitions = Environment::comm().size();
    CTX::instance()->partitionOptions.num_partitions = Environment::comm().size();
}

void
GmshMeshSeq::clear()
{
    delete M_gmodel;
}

void
GmshMeshSeq::writeToGModel()
{
    std::map<int, MVertex*> vertexMap;
    std::vector<std::vector<int> > vertexIndices(M_num_triangles);
    std::vector<int> elementNum(M_num_triangles);
    std::vector<int> elementType(M_num_triangles);
    std::vector<int> physical(M_num_triangles);
    std::vector<int> elementary(M_num_triangles);
    std::vector<int> partition(M_num_triangles);

    int cpt_node = 0;
    for (auto it=M_nodes.begin(), en=M_nodes.end(); it!=en; ++it)
    {
        vertexMap.insert(std::make_pair(cpt_node+1, new MVertex(it->coords[0], it->coords[1], 0.)));
        ++cpt_node;
    }

    int cpt_element = 0;

#if 0
    // uncomment if needed to add the edge elements
    for (auto it=M_edges.begin(), en=M_edges.end(); it!=en; ++it)
    {
        vertexIndices[cpt_element] = it->indices;
        elementNum[cpt_element] = cpt_element+1;
        elementType[cpt_element] = 1;
        physical[cpt_element] = 0;;
        elementary[cpt_element] = 0;;
        partition[cpt_element] = 0;;
        ++cpt_element;
    }
#endif

    for (auto it=M_triangles.begin(), en=M_triangles.end(); it!=en; ++it)
    {
        vertexIndices[cpt_element] = it->indices;
        elementNum[cpt_element] = cpt_element+1;
        elementType[cpt_element] = 2;
        physical[cpt_element] = 0;;
        elementary[cpt_element] = 0;;
        partition[cpt_element] = 0;;
        ++cpt_element;
    }

    M_gmodel = GModel::createGModel(
                                    vertexMap,
                                    elementNum,
                                    vertexIndices,
                                    elementType,
                                    physical,
                                    elementary,
                                    partition
                                    );
}

void
GmshMeshSeq::partition(std::string const& mshfile,
                       mesh::Partitioner const& partitioner,
                       mesh::PartitionSpace const& space,
                       std::string const& format)
{

    if ((partitioner != mesh::Partitioner::CHACO) && (partitioner != mesh::Partitioner::METIS))
        throw std::logic_error("invalid partitioner");

    if (space == mesh::PartitionSpace::MEMORY)
        this->partitionMemory(mshfile, partitioner, format);
    else if (space == mesh::PartitionSpace::DISK)
        this->partitionDisk(mshfile, partitioner, format);
    else
        throw std::logic_error("invalid partition space");
}

void
GmshMeshSeq::partitionMemory(std::string const& mshfile,
                             mesh::Partitioner const& partitioner,
                             std::string const& format)
{
    CTX::instance()->partitionOptions.partitioner =  (int)partitioner;
    CTX::instance()->partitionOptions.algorithm = 2;
    //CTX::instance()->partitionOptions.edge_matching = 3;// comment it after
    //CTX::instance()->partitionOptions.refine_algorithm = 2; // do not use because of non-contiguous mesh partition

    CTX::instance()->partitionOptions.createPartitionBoundaries = false;

    PartitionMesh( M_gmodel, CTX::instance()->partitionOptions);
    M_gmodel->writeMSH(mshfile, 2.2, (format=="binary")?true:false);

    M_gmodel->deleteMesh();
    M_gmodel->destroy();
}

void
GmshMeshSeq::partitionDisk(std::string const& mshfile,
                           mesh::Partitioner const& partitioner,
                           std::string const& format)
{
    if (!fs::exists(mshfile))
    {
        std::string const msg = "Cannot find mesh file" + mshfile;
        throw std::runtime_error(msg);
    }

    CTX::instance()->partitionOptions.partitioner =  (int)partitioner;
    CTX::instance()->mesh.binary = 1;

    std::ostringstream gmshstr;
    gmshstr << BOOST_PP_STRINGIZE( gmsh )
            << " -" << 2
            << " -part " << Environment::comm().size()
            << " -format " << "msh2";

    if (format == "binary")
    {
        gmshstr << " -bin";
    }

    gmshstr << " -string " << "\"Mesh.MshFileVersion="<< 2.2 <<";\""
            << " -string " << "\"Mesh.Partitioner="<< (int)partitioner <<";\""
            << " -string " << "\"Mesh.MetisAlgorithm="<< 2 <<";\"" // 1 = recursive (default), 2 = K-way
            << " -string " << "\"Mesh.MetisRefinementAlgorithm="<< 2 <<";\""
            << " -string " << "\"General.Verbosity="<< Environment::vm()["debugging.gmsh_verbose"].as<int>() << ";\""
            << " " << mshfile;

    LOG(DEBUG) << "[Gmsh::generate] execute '" <<  gmshstr.str() << "'\n";
    auto err = ::system( gmshstr.str().c_str() );
}

void
GmshMeshSeq::move(std::vector<double> const& um, double factor)
{
    if ((um.size() != 0) && (factor != 0))
    {
        ASSERT(2*M_nodes.size()==um.size(),"invalid size of displacement vector");

        for (int i=0; i<M_nodes.size(); ++i)
        {
            M_nodes[i].coords[0] += factor*um[i];
            M_nodes[i].coords[1] += factor*um[i+M_num_nodes];
        }
    }
}

void
GmshMeshSeq::reorder(bimap_type const& rmap_nodes, bimap_type const& rmap_elements)
{
    ASSERT(rmap_nodes.size()==M_num_nodes,"invalid size of reorder mapping (nodes)");
    ASSERT(rmap_elements.size()==M_num_triangles,"invalid size of reorder mapping (nodes)");

    auto _nodes = M_nodes;

    for (int i=0; i<M_nodes.size(); ++i)
    {
        int ri = rmap_nodes.right.find(i+1)->second-1;
        M_nodes[i].coords[0] = _nodes[ri].coords[0];
        M_nodes[i].coords[1] = _nodes[ri].coords[1];
    }

    auto _triangles = M_triangles;

    for (int i=0; i<M_num_triangles; ++i)
    {
        int ri = rmap_elements.right.find(i+1)->second-1;
        M_triangles[i] = _triangles[ri];
    }

    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        for (int i=0; i<3; ++i)
        {
            it->indices[i] = rmap_nodes.left.find(it->indices[i])->second;
        }
    }
}

void
GmshMeshSeq::stereographicProjection()
{
    // polar stereographic projection
    mapx_class *map;
    std::vector<char> str(M_mppfile.begin(), M_mppfile.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    LOG(DEBUG)<<"MFILE= "<< std::string(map->mpp_filename) <<"\n";
    LOG(DEBUG)<<"PROJE= "<< std::string(map->projection_name) <<"\n";

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

        ++cpt;
    }

    close_mapx(map);
}

std::vector<int>
GmshMeshSeq::indexTr() const
{
    std::vector<int> index(3*M_num_triangles);
    int cpt = 0;
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        index[3*cpt] = it->indices[0];//it->first;
        index[3*cpt+1] = it->indices[1];
        index[3*cpt+2] = it->indices[2];
        ++cpt;
    }

    return index;
}

std::vector<double>
GmshMeshSeq::coordX() const
{
    std::vector<double> x(M_num_nodes);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        x[cpt] = it->coords[0];
        ++cpt;
    }

    return x;
}

std::vector<double>
GmshMeshSeq::coordY() const
{
    std::vector<double> y(M_num_nodes);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        y[cpt] = it->coords[1];
        ++cpt;
    }

    return y;
}

std::vector<double>
GmshMeshSeq::coordX(double const& rotangle) const
{
    std::vector<double> x(M_num_nodes);
    int cpt = 0;
    double cos_rotangle=std::cos(rotangle);
    double sin_rotangle=std::sin(rotangle);
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        x[cpt] = cos_rotangle*(it->coords[0]) + sin_rotangle*(it->coords[1]);
        ++cpt;
    }

    return x;
}

std::vector<double>
GmshMeshSeq::coordY(double const& rotangle) const
{
    std::vector<double> y(M_num_nodes);
    int cpt = 0;
    double cos_rotangle=std::cos(rotangle);
    double sin_rotangle=std::sin(rotangle);
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        y[cpt] = -sin_rotangle*(it->coords[0]) + cos_rotangle*(it->coords[1]);
        ++cpt;
    }

    return y;
}

std::vector<double>
GmshMeshSeq::bCoordX() const
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
GmshMeshSeq::bCoordY() const
{
    std::vector<double> bcoord_y(M_num_triangles);
    int cpt = 0;
    double y = 0.;
    for (auto it=M_triangles.begin(), end=M_triangles.end(); it!=end; ++it)
    {
        y = 0.;

        for (int i=0; i<3; ++i)
        {
            y += M_nodes[it->indices[i]-1].coords[1];
        }

        bcoord_y[cpt] = y/3.;

        ++cpt;
    }

    return bcoord_y;
}

std::vector<double>
GmshMeshSeq::bCoordX(double const& rotangle) const
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
            x += cos_rotangle*(M_nodes[it->indices[i]-1].coords[0])+ sin_rotangle*(M_nodes[it->indices[i]-1].coords[1]);
        }

        bcoord_x[cpt] = x/3.;

        ++cpt;
    }

    return bcoord_x;
}

std::vector<double>
GmshMeshSeq::bCoordY(double const& rotangle) const
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
            y += -sin_rotangle*(M_nodes[it->indices[i]-1].coords[0])+ cos_rotangle*(M_nodes[it->indices[i]-1].coords[1]);
        }

        bcoord_y[cpt] = y/3.;

        ++cpt;
    }

    return bcoord_y;
}

std::vector<double>
GmshMeshSeq::meanLon() const
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
GmshMeshSeq::meanLat() const
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
GmshMeshSeq::lon() const
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
GmshMeshSeq::lat() const
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

void
GmshMeshSeq::setId(std::vector<int> const& newid)
{
    if ( newid.size() != 0 )
    {
        ASSERT(M_nodes.size()==newid.size(),"invalid size of new_id vector");

        for (int i=0; i<M_nodes.size(); ++i)
            M_nodes[i].id = newid[i];
    }
}

std::vector<int>
GmshMeshSeq::id() const
{
    std::vector<int> node_id(M_num_nodes);
    int cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        node_id[cpt] = it->id;
        ++cpt;
    }

    return node_id;
}


// ------------------------------------------------
//! return the vertices for a given list of indices
//! called by FiniteElement::shapeCoeff() and FiniteElement::jacobian()
std::vector<std::vector<double>>
GmshMeshSeq::vertices(std::vector<int> const& indices) const
{
    int const nv = indices.size();
    std::vector<std::vector<double>> vertices(nv);
    for(int i=0; i<nv; i++)
        vertices[i] = M_nodes[indices[i]-1].coords;
    return vertices;
}


std::vector<std::vector<double>>
GmshMeshSeq::vertices(std::vector<int> const& indices,
        std::vector<double> const& um, double factor) const
{
    int const nv = indices.size();
    auto vertices = this->vertices(indices);
    for(int i=0; i<nv; i++)
        for(int k=0; k<2; k++)
            vertices[i][k] += factor*um[indices[i]-1+k*M_num_nodes];
    return vertices;
}//vertices

} // Nextsim
