/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   exporter.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Mon Dec 28 15:59:43 2015
 */

#include <exporter.hpp>

namespace Nextsim
{
Exporter::Exporter()
	:
	M_mrecord(),
    M_frecord()
{}

template<typename Type>
void
Exporter::writeContainer(std::fstream& out, std::vector<Type> const& container)
{
    int fsize = container.size();

    if (out.is_open())
	{
        //out.write((char*)&name, sizeof(name));
        out.write((char*)&fsize, sizeof(fsize));
		out.write((char*)&container[0], container.size() * sizeof(Type));

        // for (int i=0; i<10; ++i)
        // {
        //     std::cout<<"Concentration["<< i <<"]= "<< container[i] <<"\n";
        // }
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
               % sizeof(int) ).str();
    M_mrecord.push_back(description);


    writeContainer(out, Mesh.coordX());
    description = (boost::format( "%1% %2%" )
               % "Nodes_x"
               % sizeof(double) ).str();
    M_mrecord.push_back(description);

    writeContainer(out, Mesh.coordY());
    description = (boost::format( "%1% %2%" )
               % "Nodes_y"
               % sizeof(double) ).str();
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
               % sizeof(Type) ).str();

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

template void Exporter::writeField<double>(std::fstream&, std::vector<double> const&, std::string const&);
template void Exporter::writeField<int>(std::fstream&, std::vector<int> const&, std::string const&);

} // Nextsim
