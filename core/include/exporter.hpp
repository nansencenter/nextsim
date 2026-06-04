/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   exporter.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Mon Dec 28 14:00:31 2015
 */

#ifndef __Exporter_HPP
#define __Exporter_HPP 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <gmshmesh.hpp>
#include <gmshmeshseq.hpp>
#include <boost/format.hpp>
#include <boost/unordered_map.hpp>

namespace Nextsim
{
struct Record
{
	std::string name;
	int size;
};

class Exporter
{
public:

	Exporter(std::string const& precision = "float");

    template<typename Type>
    void writeContainer(std::fstream& out, std::vector<Type> const& container, std::string const precision);

    void writeMesh(std::fstream& out, GmshMesh const& Mesh);
    void writeMesh(std::fstream& out, GmshMeshSeq const& Mesh);
    template<typename Type>
	void writeMesh(std::fstream& out, std::vector<Type> const& xnod, std::vector<Type> const& ynod,
            std::vector<int> const& idnod, std::vector<int> const& elements);

    template<typename Type>
	void writeField(std::fstream& out, std::vector<Type> const& field, std::string const& name);
	void writeRecord(std::fstream& out, std::string const& rtype = "field");

    void loadFile(std::fstream &in, boost::unordered_map<std::string, std::vector<int>> &field_map_int, boost::unordered_map<std::string, std::vector<double>> &field_map_dbl);
    void readRecord(std::ifstream &in);

private:

    std::vector<std::string> M_mrecord;
    std::vector<std::string> M_frecord;

    //std::vector<int> M_type_record;
    std::vector<std::string> M_type_record;
    std::vector<std::string> M_name_record;

    std::string M_precision;
};
} // Nextsim
#endif
