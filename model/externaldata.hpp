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
//#include <boost/ptr_container/ptr_vector.hpp>
#include <gmshmesh.hpp>
#include <Bamgx.h>
#include <InterpFromGridToMeshx.h>
#include <InterpFromMeshToMesh2dx.h>
#include <BamgTriangulatex.h>
#include <netcdf>



/**
 * @class ExternalData
 * @brief Manage the external data
 *
 * @see
 *
 */

namespace Nextsim
{
    enum InterpolationType
    {
        FromGridToMesh = 0,
        FromMeshToMesh2dx = 1,
        FromMeshToMesh2dCavities = 2,
    };

class ExternalData
{

public:

	typedef std::size_t size_type;
    typedef double value_type;

    typedef struct Dimension
    {
        std::string name;
        int start;
        int end;
        bool cyclic; // if cyclic, then the first value will also be used for interpolation after the last value
    } Dimesion;

    typedef struct Variable
    {
        std::string name;
        std::vector<Dimension> dimensions;
        double a;
        double b;
        std::string Units;
        std::vector<std::vector<double>> data2;
    } Variable;

    typedef struct Vectorial_Variable
    {
        std::vector<int> components_Id;
        bool east_west_oriented;
    } Vectorial_Variable;

    typedef struct Grid
    {
        InterpolationType interpolation_method;
		int interp_type;
        std::string dirname;
        std::string filename;

        Variable latitude;
        Variable longitude;

        Dimension dimension_x;
        Dimension dimension_y;

        std::string mpp_file;
		bool interpolation_in_latlon;

        bool loaded;

		bool masking;
		Variable masking_variable;
		std::vector<int> reduced_nodes_ind;

        int* pfindex;
        int pfnels;
        std::vector<double> gridX;
        std::vector<double> gridY;

        std::vector<double> gridLAT;
        std::vector<double> gridLON;
    } Grid;

    typedef struct Dataset
    {
        std::string dirname;
        std::string prefix;
        std::string postfix;
        std::string reference_date;

        std::vector<Variable> variables;
        std::vector<Vectorial_Variable> vectorial_variables;
        int target_size;
        Grid *grid;

        bool reloaded;

        int nb_timestep_day;
        Variable time;
        Dimension dimension_time;

        std::vector<double> ftime_range;
    } Dataset;

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
