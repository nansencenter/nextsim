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
Exporter::writeContainer(std::fstream& out, std::vector<Type> const& container, std::string const precision)
{
    if (out.is_open())
	{
        int fsize = container.size();
        out.write((char*)&fsize, sizeof(fsize)); // write first the record length

        int typesize = sizeof(Type);
        if ((precision != "double") && (typesize == sizeof(double)))
        {
            // if input is double, but only ask for float precision,
            // convert to float before writing into binary file
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
		//std::cout << "Cannot open " << out  << "\n";
		//std::cerr << "error: open file " << out << " for output failed!" <<"\n";
		std::abort();
	}
}

void
Exporter::writeMesh(std::fstream& out, GmshMesh const& Mesh)
{
    auto xnod = Mesh.coordX();
    auto ynod = Mesh.coordY();
    auto elements = Mesh.indexTr();
    auto idnod = Mesh.id();
    this->writeMesh(out, xnod, ynod, idnod, elements);
}

template<typename Type>
void
Exporter::writeMesh(std::fstream& out, std::vector<Type> const &xnod, std::vector<Type> const &ynod,
        std::vector<int> const &idnod, std::vector<int> const &elements)
{
    std::string description;

    std::vector<int> field_int = elements;
    writeContainer(out, field_int, "int");
    description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                   % "Elements"
                   % "int"
                   % field_int.size()
                   % *std::min_element(field_int.begin(), field_int.end())
                   % *std::max_element(field_int.begin(), field_int.end()) ).str();
    M_mrecord.push_back(description);

    field_int = idnod;
    writeContainer(out, field_int, "int");
    description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                   % "id"
                   % "int"
                   % field_int.size()
                   % *std::min_element(field_int.begin(), field_int.end())
                   % *std::max_element(field_int.begin(), field_int.end()) ).str();
    M_mrecord.push_back(description);

    std::string prec = "double";
    if (sizeof(Type)==sizeof(float))
        prec = "float";

    auto field = xnod;
    writeContainer(out, field, prec);
    description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                   % "Nodes_x"
                   % prec
                   % field.size()
                   % *std::min_element(field.begin(), field.end())
                   % *std::max_element(field.begin(), field.end()) ).str();
    M_mrecord.push_back(description);

    field = ynod;
    writeContainer(out, field, prec);
    description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                   % "Nodes_y"
                   % prec
                   % field.size()
                   % *std::min_element(field.begin(), field.end())
                   % *std::max_element(field.begin(), field.end()) ).str();

    M_mrecord.push_back(description);
}


template<typename Type>
void
Exporter::writeField(std::fstream& out, std::vector<Type> const& field, std::string const& name)
{
    //std::cout<<"Writing "<<name<<": len = "<<field.size()<<"\n";
    std::string description;

    std::string precision;
    // We must choose either integer or floating point
    if ( sizeof(Type) == sizeof(int) )
        precision = "int";
    else
        precision = M_precision;

    // Time should always be in double
    if ( name == "Time" )
        precision = "double";

    writeContainer(out, field, precision);
    description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                   % name
                   % precision
                   % field.size()
                   % ( (field.size() > 0) ? *std::min_element(field.begin(), field.end()) : 0 )
                   % ( (field.size() > 0) ? *std::max_element(field.begin(), field.end()) : 0 ) ).str();

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
//		std::cout << "Cannot open " << out  << "\n";
//		std::cerr << "error: open file " << out << " for output failed!" <<"\n";
		std::abort();
	}
}

void
Exporter::readRecord(std::ifstream& in)
{
    std::string name, type, size, min_element, max_element;

    while ( in >> name >> type >> size >> min_element >> max_element)
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
            std::string error_str = "Unknow typ in file: " + M_type_record[i];
            throw std::logic_error(error_str);
        }
    }
}

template void Exporter::writeField<double>(std::fstream&, std::vector<double> const&, std::string const&);
template void Exporter::writeField<int>(std::fstream&, std::vector<int> const&, std::string const&);

template void Exporter::writeMesh<double>(std::fstream&, std::vector<double> const&, std::vector<double> const&,
        std::vector<int> const&,std::vector<int> const&);
template void Exporter::writeMesh<float>(std::fstream&, std::vector<float> const&, std::vector<float> const&,
        std::vector<int> const&,std::vector<int> const&);

} // Nextsim
