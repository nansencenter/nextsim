/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   gmshmeshseq.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Tue Jun 21 14:25:38 2016
 */

#include <gmshmeshseq.hpp>

namespace Nextsim
{
GmshMeshSeq::GmshMeshSeq()
    :
    M_nodes(),
    M_triangles(),
    M_edges(),
    M_num_nodes(0),
    M_num_triangles(0),
    M_num_edges(0)
{}

GmshMeshSeq::GmshMeshSeq(std::vector<point_type> const& nodes,
                         std::vector<element_type> const& edges,
                         std::vector<element_type> const& triangles)
    :
    M_nodes(nodes),
    M_triangles(triangles),
    M_edges(edges),
    M_num_nodes(nodes.size()),
    M_num_triangles(triangles.size()),
    M_num_edges(edges.size())
{}

void
GmshMeshSeq::writeTofile(std::string const& filename)
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
    }
    else
    {
        std::cout << "Cannot open " << gmshmshfile  << "\n";
        std::cerr << "error: open file " << gmshmshfile << " for output failed!" <<"\n";
        std::abort();
    }
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
GmshMeshSeq::stereographicProjection()
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
GmshMeshSeq::bcoordX() const
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
GmshMeshSeq::bcoordY() const
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
GmshMeshSeq::bcoordX(double const& rotangle) const
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
GmshMeshSeq::bcoordY(double const& rotangle) const
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
    std::string filename = Environment::nextsimDir().string() + "/data/NpsNextsim.mpp";
    std::vector<char> str(filename.begin(), filename.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    std::vector<double> mean_lon(M_num_triangles);
    double lat = 0.;
    double lon = 0.;

    std::vector<double> X = this->bcoordX();
    std::vector<double> Y = this->bcoordY();

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
    std::string filename = Environment::nextsimDir().string() + "/data/NpsNextsim.mpp";
    std::vector<char> str(filename.begin(), filename.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    std::vector<double> mean_lat(M_num_triangles);
    double lat = 0.;
    double lon = 0.;

    std::vector<double> X = this->bcoordX();
    std::vector<double> Y = this->bcoordY();

    for (int elt=0; elt<M_num_triangles; ++elt)
    {
        int status = inverse_mapx(map,X[elt],Y[elt],&lat,&lon);
        mean_lat[elt] = lat;
    }

    close_mapx(map);

    return mean_lat;
}

} // Nextsim
