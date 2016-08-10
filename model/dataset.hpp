/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   dataset.hpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#ifndef __DataSet_H
#define __DataSet_H 1

#include <environment.hpp>
#include <boost/ptr_container/ptr_vector.hpp>


/**
 * @class DataSet
 * @brief Describe the dataset
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

class DataSet
{

public:

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

};

} // Nextsim
#endif // __DataSet_H
