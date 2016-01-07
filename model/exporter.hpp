/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   exporter.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Dec 28 14:00:31 2015
 */

#ifndef __Exporter_HPP
#define __Exporter_HPP 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <gmshmesh.hpp>
#include <boost/format.hpp>

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

	Exporter();

    template<typename Type>
    void writeContainer(std::fstream& out, std::vector<Type> const& container);
    void writeMesh(std::fstream& out, GmshMesh const& Mesh);
    template<typename Type>
	void writeField(std::fstream& out, std::vector<Type> const& field, std::string const& name);
	void writeRecord(std::fstream& out, std::string const& rtype = "field");

private:

    std::vector<std::string> M_mrecord;
    std::vector<std::string> M_frecord;

};
} // Nextsim
#endif
