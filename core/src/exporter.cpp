/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   exporter.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Mon Dec 28 15:59:43 2015
 */

#include <exporter.hpp>

namespace Nextsim
{
Exporter::Exporter(std::string const& precision)
	:
	M_mrecord(),
    M_frecord(),
    M_type_record(),
    M_name_record(),
    M_precision(precision)
{}

template<typename Type>
void
Exporter::writeContainer(std::fstream& out, std::vector<Type> const& container)
{
    if (out.is_open())
	{
        int fsize = container.size();
        out.write((char*)&fsize, sizeof(fsize)); // write first the record length

        int typesize = sizeof(Type);
        if ((M_precision != "double") && (typesize == sizeof(double)))
        {
            for (int i=0; i<container.size(); ++i)
            {
                // convert to float before writing into binary file
                float fvalue = (float)container[i];
                out.write((char*)&fvalue, sizeof(float));
            }
        }
        else
        {
            // write as original format
            out.write((char*)&container[0], container.size() * sizeof(Type));
        }
	}
	else
	{
		std::cout << "Cannot open " << out  << "\n";
		std::cerr << "error: open file " << out << " for output failed!" <<"\n";
		std::abort();
	}
}

void
Exporter::writeMesh(std::fstream& out, GmshMesh const& Mesh)
{
    std::string description;

    writeContainer(out, Mesh.indexTr());
    description = (boost::format( "%1% %2%" )
                   % "Elements"
                   % "int" ).str();
    M_mrecord.push_back(description);


    writeContainer(out, Mesh.coordX());
    description = (boost::format( "%1% %2%" )
                   % "Nodes_x"
                   % M_precision ).str();
    M_mrecord.push_back(description);

    writeContainer(out, Mesh.coordY());
    description = (boost::format( "%1% %2%" )
                   % "Nodes_y"
                   % M_precision ).str();

    M_mrecord.push_back(description);

    writeContainer(out, Mesh.id());
    description = (boost::format( "%1% %2%" )
                   % "id"
                   % "int" ).str();
    M_mrecord.push_back(description);
}

template<typename Type>
void
Exporter::writeField(std::fstream& out, std::vector<Type> const& field, std::string const& name)
{
    std::string description;
    writeContainer(out, field);

    description = (boost::format( "%1% %2%" )
                   % name
                   % M_precision ).str();

	M_frecord.push_back(description);
}

void
Exporter::writeRecord(std::fstream& out, std::string const& rtype)
{
	if (out.is_open())
	{
        if (rtype == "field")
        {
            for (auto const& str : M_frecord)
            {
                //out << std::setw(15) << std::left << R.name;
                out << str <<"\n";
            }
        }
        else
        {
            for (auto const& str : M_mrecord)
            {
                //out << std::setw(15) << std::left << R.name;
                out << str <<"\n";
            }
        }
	}
	else
	{
		std::cout << "Cannot open " << out  << "\n";
		std::cerr << "error: open file " << out << " for output failed!" <<"\n";
		std::abort();
	}
}

void
Exporter::readRecord(std::ifstream& in)
{
    std::string name;
    std::string type;

    while ( in >> name >> type )
    {
        M_name_record.push_back(name);
        M_type_record.push_back(type);
    }
}

void
Exporter::loadFile(std::fstream& in, boost::unordered_map<std::string, std::vector<int>>& field_map_int, boost::unordered_map<std::string, std::vector<double>>& field_map_dbl)
{
    int reclen;
    for ( int i=0; i<M_type_record.size(); ++i )
    {
        in.read((char*) &reclen, sizeof(reclen));

        if ( M_type_record[i] == "double" )
        {
            std::vector<double> dvec(reclen);
            in.read((char*) &dvec[0], reclen*sizeof(double));
            field_map_dbl.emplace(M_name_record[i], dvec);
        }
        else if ( M_type_record[i] == "int" )
        {
            std::vector<int> ivec(reclen);
            in.read((char*) &ivec[0], reclen*sizeof(int));
            field_map_int.emplace(M_name_record[i], ivec);
        }
        else
        {
            throw std::logic_error("unknown type in file");
        }
    }
}

template void Exporter::writeField<double>(std::fstream&, std::vector<double> const&, std::string const&);
template void Exporter::writeField<int>(std::fstream&, std::vector<int> const&, std::string const&);

} // Nextsim
