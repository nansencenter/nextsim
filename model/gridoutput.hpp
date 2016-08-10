/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

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
#include <gmshmesh.hpp>
#include <Bamgx.h>
#include <InterpFromMeshToGridx.h>
#include <BamgTriangulatex.h>
#include <netcdf>
#include <externaldata.hpp>

/**
 * @class GridOutput
 * @brief Manage output of gridded data
 *
 * @see
 *
 */

namespace Nextsim
{
    enum variableID
    {
        // Prognostic variables
        conc = 0,
        thick = 1,
        damage = 2,
        snow_thick = 3,
        VT_x = 4,
        VT_y = 5,

        // Diagnostic variables
    };

    class GridOutput: private ExternalData
    {
        public:

            //using ExternalData::Variable;
            using ExternalData::Dimension;
            using ExternalData::Grid;
            //using ExternalData::Dataset;
            using ExternalData::loadGrid;
    
            typedef struct Variable
            {
                std::string name;
                std::string longName;
                std::string stdName;
                std::vector<Dimension> dimensions;
                double sclFac;
                double addOff;
                std::string Units;
                std::vector<double> data_mesh;
                std::vector<double> data_grid;
                int variableID;
            } Variable;

            GridOutput();

            GridOutput(int ncols, int nrows, double mooring_spacing, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables);

            GridOutput(Grid grid, double mooring_spacing, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables);

            ~GridOutput();

            void updateGridMean(GmshMesh const &mesh);
            void resetGridMean();
            void resetMeshMean(GmshMesh const &mesh);
            void exportGridMeans(std::string filename, int time_step, int mooring_output_time_step);

            int M_ncols;
            int M_nrows;
            double M_mooring_spacing;
            std::vector<Variable> M_nodal_variables;
            std::vector<Variable> M_elemental_variables;

        private:
            double M_miss_val = -1e24;

            void updateGridMeanWorker(GmshMesh const &mesh, int mesh_size, std::vector<Variable> &variables);
    };
}
#endif // __GridOutput_H
