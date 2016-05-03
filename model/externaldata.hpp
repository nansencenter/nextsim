/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   externaldata.hpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#ifndef __ExternalData_H
#define __ExternalData_H 1

#include <environment.hpp>
     //#include <boost/filesystem/operations.hpp>
     //#include <boost/filesystem/fstream.hpp>
     //#include <boost/format.hpp>
#include <assert.hpp>
#include <finiteelement.hpp>

/**
 * @class ExternalData
 * @brief Manage the external data
 *
 * @see
 *
 */

namespace Nextsim
{

class ExternalData
{

public:

	typedef std::size_t size_type;
    typedef double value_type;

    typedef FiniteElement::Dataset Dataset;
    typedef FiniteElement::Grid Grid;

    ExternalData(Dataset *dataset, GmshMesh const& Mesh, int VariableId );

	~ExternalData();
#if 0
    void init( const size_type n, bool fast = false );

    void resize(const size_type n, bool fast = false);

	void close();
#endif
    void settime( const double current_time );
    
	value_type operator() (const size_type i);

	void clear();
    
    void loadDataset(Dataset *dataset, GmshMesh const& Mesh);//(double const& u, double const& v);
    
    void loadGrid(Grid *grid);

#if 0
	size_type size() const;

	bool isInitialized() const { return M_is_initialized; }

	bool closed() const { return M_is_closed; }
#endif

public:
    double fdt;
    std::vector<double> fcoeff;
    Dataset *M_dataset;
    GmshMesh M_mesh;

private:
    std::string M_datasetname;
    int M_VariableId;
	bool M_is_initialized;
    double M_current_time;
    

};

} // Nextsim
#endif // __ExternalData_H
