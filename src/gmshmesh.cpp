/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   gmshmesh.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Wed Jul 29 15:26:47 2015
 */

#include <gmshmesh.hpp>

namespace Nextsim
{
	GmshMesh::GmshMesh()
		:
		M_version("2.2"),
        M_ordering("gmsh"),
        M_nodes(),
        M_elements(),
        M_triangles(),
        M_lines(),
        M_num_nodes(0),
        M_num_elements(0),
        M_num_triangles(0),
        M_num_lines(0)
	{}

    void GmshMesh::readFromFile(std::string const& filename)
    {
        std::cout<<"Reading Msh file "<< filename <<"\n";

        std::ifstream __is ( filename.c_str() );

        if ( !__is.is_open() )
        {
            std::ostringstream ostr;
            std::cout << "Invalid file name " << filename << " (file not found)\n";
            ostr << "Invalid file name " << filename << " (file not found)\n";
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

        std::cout << "buf: "<< __buf << "\n";

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

        std::vector<double> coords(3,0);

        for ( unsigned int __i = 0; __i < __n; ++__i )
        {
            int id = 0;

            __is >> id
                 >> coords[0]
                 >> coords[1]
                 >> coords[2];

            M_nodes[id].id = id;
            M_nodes[id].coords = coords;
        }

        __is >> __buf;
        std::cout << "buf: "<< __buf << "\n";

        // make sure that we have read all the points

        ASSERT(std::string( __buf ) == "$EndNodes","invalid end nodes string");

        // Read ELEMENTS

        __is >> __buf;

        ASSERT(std::string( __buf ) == "$Elements","invalid elements string");

        int numElements;
        __is >> numElements;

        M_num_elements = numElements;

        std::cout << "Reading " << numElements << " elements...\n";
        //std::list<Nextsim::entities::GMSHElement> __et; // tags in each element
        std::map<int,int> __gt;

        for(int i = 0; i < numElements; i++)
        {
            int number, type, physical = 0, elementary = 0, numVertices;
            int numTags;

            __is >> number  // elm-number
                 >> type // elm-type
                 >> numTags; // number-of-tags

            for(int j = 0; j < numTags; j++)
            {
                int tag;
                __is >> tag;
                if(j == 0) physical = tag;
                else if(j == 1) elementary = tag;
            }

            numVertices = MElement::getInfoMSH(type);

            ASSERT(numVertices!=0,"unknown number of vertices for element type");

            std::vector<int> indices(numVertices);
            for(int j = 0; j < numVertices; j++)
            {
                __is >> indices[j];
            }

            if (M_ordering=="bamg")
            {
                std::next_permutation(indices.begin()+1,indices.end());
            }

            Nextsim::entities::GMSHElement gmshElt( number,
                                                    type,
                                                    physical,
                                                    elementary,
                                                    numVertices,
                                                    indices );

            //__et.push_back( gmshElt );
            M_elements.insert(std::make_pair(number,gmshElt));

            if (type == 2)
                M_triangles.insert(std::make_pair(number,gmshElt));
            else if (type == 1)
                M_lines.insert(std::make_pair(number,gmshElt));


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
                M_num_lines = it.second;
        }

        // make sure that we have read everything
        __is >> __buf;

        ASSERT(std::string( __buf ) == "$EndElements","invalid end elements string");

        // we are done reading the MSH file
    }

    void GmshMesh::writeTofile(std::string const& filename)
    {
        //std::string gmshfilename = (boost::format( "../data/arctic10km.msh" ) ).str();
        std::fstream gmshfile(filename, std::ios::out | std::ios::trunc);

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
                //<< "  " << element + 1;

                for (int i = 0; i < 3; i++ )
                {
                    gmshfile << "  " << it->second.indices[i] ;
                }
                gmshfile << "\n";

                ++element;
            }
            gmshfile << "$EndElements\n";


        }
        else
        {
            std::cout << "Cannot open " << filename  << "\n";
            std::cerr << "error: open file " << filename << " for output failed!" <<"\n";
            std::abort();
        }
    }

    void GmshMesh::project(std::string const& filename)
    {
        // polar stereographic projection
        mapx_class *map;
        std::vector<char> str(filename.begin(), filename.end());
        str.push_back('\0');

        map = init_mapx(&str[0]);

        std::cout<<"MFILE= "<< std::string(map->mpp_filename) <<"\n";
        std::cout<<"PROJE= "<< std::string(map->projection_name) <<"\n";

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
    }

} // Nextsim
