/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   externaldata.hpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#ifndef __ExternalData_H
#define __ExternalData_H 1

#include <environment.hpp>
#include <assert.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <gmshmesh.hpp>
#include <Bamgx.h>
#include <InterpFromGridToMeshx.h>
#include <InterpFromMeshToMesh2dx.h>
#include <BamgTriangulatex.h>
#include <netcdf>
#include <dataset.hpp>


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
    
    typedef DataSet::Dataset Dataset;
    typedef DataSet::Grid Grid;
    typedef DataSet::Dimension Dimension;
    typedef DataSet::Variable Variable;
    typedef DataSet::Vectorial_Variable Vectorial_Variable;

    ExternalData();

    ExternalData(Dataset *dataset, GmshMesh const& mesh, int VariableId, bool is_vector );

    ExternalData(Dataset *dataset, GmshMesh const& mesh, int VariableId, bool is_vector,
        double SpinUpStartingTime, double SpinUpDuration  );

    ExternalData(double ConstantValue );

    ExternalData(double ConstantValue, double ConstantValuebis );

    ExternalData(double ConstantValue,
        double SpinUpStartingTime, double SpinUpDuration  );

    ExternalData(double ConstantValue, double ConstantValuebis,
        double SpinUpStartingTime, double SpinUpDuration  );

	~ExternalData();
#if 0
    void init( const size_type n, bool fast = false );

    void resize(const size_type n, bool fast = false);

	void close();
#endif
    void check_and_reload(GmshMesh const& mesh, const double current_time );

	value_type operator[] (const size_type i);

	void clear();

    void loadDataset(Dataset *dataset, GmshMesh const& mesh);//(double const& u, double const& v);

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

private:
    std::string M_datasetname;
    int M_VariableId;
    double M_current_time;
    bool M_is_constant;
    double M_constant_value;
    double M_constant_valuebis;
    bool M_is_vector;
    int M_VariableIdbis;
    double M_SpinUpStartingTime;
    double M_SpinUpDuration;
    double M_factor;
};

} // Nextsim
#endif // __ExternalData_H
