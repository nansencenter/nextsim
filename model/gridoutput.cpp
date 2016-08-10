/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   gridoutput.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Thu Aug  4 09:46:28 CEST 2016
 */

#include <gridoutput.hpp>

/**
 * @class GridOutput
 * @brief Manage output of gridded data
 *
 * @see
 *
*/

namespace Nextsim
{
    GridOutput::GridOutput()
    {}

    GridOutput::GridOutput(int ncols, int nrows, double mooring_spacing, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables)
        :
            M_ncols(ncols),
            M_nrows(nrows),
            M_mooring_spacing(mooring_spacing),
            M_nodal_variables(nodal_variables),
            M_elemental_variables(elemental_variables)
    {}

    GridOutput::GridOutput(Grid grid, double mooring_spacing, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables)
        :
            M_mooring_spacing(mooring_spacing),
            M_nodal_variables(nodal_variables),
            M_elemental_variables(elemental_variables)
    {}

    GridOutput::~GridOutput()
    {}

    // Interpolate from the mesh values to the grid
    void GridOutput::updateGridMean(GmshMesh const &mesh)
    {
        // Just call the worker routine once for each dataset
        updateGridMeanWorker(mesh, mesh.numTriangles(), M_elemental_variables);
        updateGridMeanWorker(mesh, mesh.numNodes(), M_nodal_variables);
    }

    // Interpolate from the mesh to the grid - updateing the gridded mean
    void GridOutput::updateGridMeanWorker(GmshMesh const &mesh, int mesh_size, std::vector<Variable> &variables)
    {
        int nb_var=variables.size();

        // TODO: We should check and make sure all variables in the vector are the same size
        // Input vector and output pointer
        std::vector<double> interp_in(nb_var*mesh_size);
        double* interp_out;

        // Stuff the input vector
        int tmp_nb_var=0;
        for (int i=0; i<mesh_size; ++i)
            for (int j=0; j<nb_var; j++)
                interp_in[nb_var*i+j] = variables[j].data_mesh[i];

        // interpolation from mesh to grid
        auto RX = mesh.coordX();
        auto RY = mesh.coordY();
        double xmin = *std::min_element( RX.begin(), RX.end() );
        double ymax = *std::max_element( RY.begin(), RY.end() );
        double dx = 1e3 * M_mooring_spacing;

        InterpFromMeshToGridx(interp_out,
                              &mesh.indexTr()[0],&mesh.coordX()[0],&mesh.coordY()[0],
                              mesh.numNodes(),mesh.numTriangles(),
                              &interp_in[0],
                              mesh_size, nb_var,
                              xmin,ymax,
                              dx,dx,
                              M_nrows, M_ncols,
                              M_miss_val);


        // Add the output pointer value to the grid vectors
        for (int i=0; i<M_nrows*M_ncols; ++i)
            for (int j=0; j<nb_var; j++)
                variables[j].data_grid[i] += interp_out[nb_var*i+j];

        xDelete<double>(interp_out);
    }

    // Set the _grid values back to zero
    void GridOutput::resetGridMean()
    {
        for (int i=0; i<M_nrows*M_ncols; ++i)
        {
            for (int j=0; j<M_nodal_variables.size(); j++)
                M_nodal_variables[j].data_grid[i] = 0.;

            for (int j=0; j<M_elemental_variables.size(); j++)
                M_elemental_variables[j].data_grid[i] = 0.;
        }
    }

    // Set the _mesh values back to zero
    void GridOutput::resetMeshMean(GmshMesh const &mesh)
    {
        for (int i=0; i<M_elemental_variables.size(); i++)
            M_elemental_variables[i].data_mesh.assign(mesh.numTriangles(), 0.);

        for (int i=0; i<M_nodal_variables.size(); i++)
            M_nodal_variables[i].data_mesh.assign(mesh.numNodes(), 0.);
    }

    // Save the _grid values to file ... this needs a lot of work, NetCDF etc.
    void GridOutput::exportGridMeans(std::string filename, int time_step, int mooring_output_time_step)
    {
        double time_factor = time_step/mooring_output_time_step;

        // Divide the element values with the time factor
        for (int i=0; i<M_nrows*M_ncols; ++i)
        {
            for (int j=0; j<M_nodal_variables.size(); j++)
                M_nodal_variables[j].data_grid[i] *= time_factor;

            for (int j=0; j<M_elemental_variables.size(); j++)
                M_elemental_variables[j].data_grid[i] *= time_factor;
        }

        // Save to file - this is still just an ascii dump of one parameter!
        std::ofstream myfile;
        myfile.open(filename);
        std::copy(M_elemental_variables[0].data_grid.begin(), M_elemental_variables[0].data_grid.end(), ostream_iterator<float>(myfile, " "));
        myfile.close();

        // Reset the grid
        this->resetGridMean();
    }
}



