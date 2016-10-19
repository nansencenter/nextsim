/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gridoutput.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Thu Aug  4 09:47:27 CEST 2016
 */

#ifndef __GridOutput_H
#define __GridOutput_H 1

#include <environment.hpp>
#include <assert.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/format.hpp>
#include <gmshmesh.hpp>
#include <Bamgx.h>
#include <InterpFromMeshToGridx.h>
#include <BamgTriangulatex.h>
#include <netcdf>
#include <dataset.hpp>

/**
 * @class GridOutput
 * @brief Manage output of gridded data
 *
 * @see
 *
 */

namespace Nextsim
{
    class GridOutput: public DataSet
    {
        public:

            enum variableID
            {
                // Prognostic variables
                conc        =  0,
                thick       =  1,
                damage      =  2,
                snow_thick  =  3,
                VT_x        =  4,
                VT_y        =  5,
                tsurf       =  6,
                sst         =  7,
                sss         =  8

                // Diagnostic variables
            };

            enum variableKind
            {
                nodal      =  0,
                elemental  =  1
            };

            typedef struct Variable
            {
                std::string name;
                std::string longName;
                std::string stdName;
                std::vector<Dimension> dimensions;
                std::string Units;
                std::vector<double> data_mesh;
                std::vector<double> data_grid;
                int variableID;
            } Variable;

            GridOutput();

            GridOutput(int ncols, int nrows, double mooring_spacing, std::vector<Variable> variables, int kind);
            GridOutput(Grid grid, std::vector<Variable> variables, int kind);

            GridOutput(int ncols, int nrows, double mooring_spacing, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables);
            GridOutput(Grid grid, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables);

            GridOutput(int ncols, int nrows, double mooring_spacing, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables);
            GridOutput(Grid grid, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables);

            ~GridOutput();

            void updateGridMean(GmshMesh const &mesh);
            void resetGridMean();
            void resetMeshMean(GmshMesh const &mesh);
            void exportGridMeans(std::string filename, double time_step, double mooring_output_time_step);

            int M_ncols;
            int M_nrows;
            double M_mooring_spacing;
            int M_grid_size;
            Grid M_grid;
            std::vector<Variable> M_nodal_variables;
            std::vector<Variable> M_elemental_variables;
            std::vector<Vectorial_Variable> M_vectorial_variables;

            double M_miss_val = -1e24;

        private:

            GridOutput(std::vector<Variable> variables, int kind);

            GridOutput(std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables);

            GridOutput(std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables);

            void initRegularGrid(int ncols, int nrows, double mooring_spacing);

            void initArbitraryGrid(Grid grid);

            void updateGridMeanWorker(GmshMesh const &mesh, int mesh_size, std::vector<Variable> &variables);

            void rotateVectors(GmshMesh const &mesh, Vectorial_Variable const &vectorial_variable, std::vector<Variable> &variables);
    };
}
#endif // __GridOutput_H
