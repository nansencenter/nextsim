/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gridoutput.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Thu Aug  4 09:46:28 CEST 2016
 */

#include <gridoutput.hpp>
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
    GridOutput::GridOutput()
    {}

    ////////////////////////////////////////////////////////////////////////////////
    // Constructor for only one set of variables
    ////////////////////////////////////////////////////////////////////////////////
    GridOutput::GridOutput(std::vector<Variable> variables, int kind)
    {
        M_vectorial_variables.resize(0);

        switch (kind)
        {
            case variableKind::nodal:
                M_nodal_variables = variables;
                M_elemental_variables.resize(0);
                break;

            case variableKind::elemental:
                M_elemental_variables = variables;
                M_nodal_variables.resize(0);

            default:
                std::logic_error("Incorrect varible kind in GridOutput::GridOutput");
        }
    }

    // Constructor for only one set of variables - regular grid
    GridOutput::GridOutput(int ncols, int nrows, double mooring_spacing, std::vector<Variable> variables, int kind)
        :
            GridOutput(variables, kind)
    {
        initRegularGrid(ncols, nrows, mooring_spacing);
    }

    // Constructor for only one set of variables - arbitrary grid
    GridOutput::GridOutput(Grid grid, std::vector<Variable> variables, int kind)
        :
            GridOutput(variables, kind)
    {
        initArbitraryGrid(grid);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // constructor for nodal and elemental variables only (no vectors)
    ////////////////////////////////////////////////////////////////////////////////
    GridOutput::GridOutput(std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables)
        :
            M_nodal_variables(nodal_variables),
            M_elemental_variables(elemental_variables)
    {
        M_vectorial_variables.resize(0);
    }

    // constructor for nodal and elemental variables only (no vectors) - regular grid
    GridOutput::GridOutput(int ncols, int nrows, double mooring_spacing, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables)
        :
            GridOutput(nodal_variables, elemental_variables)

    {
        initRegularGrid(ncols, nrows, mooring_spacing);
    }

    // constructor for nodal and elemental variables only (no vectors) - arbitrary grid
    GridOutput::GridOutput(Grid grid, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables)
        :
            GridOutput(nodal_variables, elemental_variables)
    {
        initArbitraryGrid(grid);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // constructor for nodal, elemental and vectorial variables
    ////////////////////////////////////////////////////////////////////////////////
    GridOutput::GridOutput(std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables)
        :
            M_nodal_variables(nodal_variables),
            M_elemental_variables(elemental_variables),
            M_vectorial_variables(vectorial_variables)
    {}

    // constructor for nodal, elemental and vectorial variables - regular grid
    GridOutput::GridOutput(int ncols, int nrows, double mooring_spacing, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables)
        :
            GridOutput(nodal_variables, elemental_variables, vectorial_variables)

    {
        initRegularGrid(ncols, nrows, mooring_spacing);
    }

    // constructor for nodal, elemental and vectorial variables - arbitrary grid
    GridOutput::GridOutput(Grid grid, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables)
        :
            GridOutput(nodal_variables, elemental_variables, vectorial_variables)
    {
        initArbitraryGrid(grid);
    }

    GridOutput::~GridOutput()
    {}

    ////////////////////////////////////////////////////////////////////////////////
    // Initialisation routines for the two kinds of grids
    ////////////////////////////////////////////////////////////////////////////////
    void GridOutput::initRegularGrid(int ncols, int nrows, double mooring_spacing)
    {
        M_ncols = ncols;
        M_nrows = nrows;
        M_mooring_spacing = mooring_spacing;
        M_grid_size = M_ncols*M_nrows;

        M_grid = Grid();
        M_grid.loaded = false;

        this->resetGridMean();
    }

    void GridOutput::initArbitraryGrid(Grid grid)
    {
        M_grid = grid;
        loadGrid(&M_grid, -1); // We don't expect there to be a time dependence for which grid to load

        M_ncols = -1;
        M_nrows = -1;
        M_mooring_spacing = -1.;
        M_grid_size = M_grid.gridX.size();

        this->resetGridMean();
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Functions other than construction and initialisation
    ////////////////////////////////////////////////////////////////////////////////

    // Interpolate from the mesh values to the grid
    void GridOutput::updateGridMean(GmshMesh const &mesh)
    {
        // Call the worker routine for the elements
        updateGridMeanWorker(mesh, mesh.numTriangles(), M_elemental_variables);

        // Rotate vectors if needed (these are assumed to be on the nodes)
        for ( auto it=M_vectorial_variables.begin(); it!=M_vectorial_variables.end(); it++ )
            if ( it->east_west_oriented | M_grid.loaded )
                rotateVectors(mesh, *it, M_nodal_variables);

        // Call the worker routine for the nodes
        updateGridMeanWorker(mesh, mesh.numNodes(), M_nodal_variables);
    }

    // Interpolate from the mesh to the grid - updateing the gridded mean
    void GridOutput::updateGridMeanWorker(GmshMesh const &mesh, int source_size, std::vector<Variable> &variables)
    {
        int nb_var=variables.size();

        // TODO: We should check and make sure all variables in the vector are the same size
        // Input vector and output pointer
        std::vector<double> interp_in(nb_var*source_size);
        double* interp_out;

        // Stuff the input vector
        for (int i=0; i<source_size; ++i)
            for (int j=0; j<nb_var; j++)
                interp_in[nb_var*i+j] = variables[j].data_mesh[i];

        // interpolation from mesh to grid
        auto RX = mesh.coordX();
        auto RY = mesh.coordY();
        double xmin = *std::min_element( RX.begin(), RX.end() );
        double ymax = *std::max_element( RY.begin(), RY.end() );

        // At the moment a non-regular grid, loaded into M_grid is handled by InterpFromMeshToMesh2dx.
        // Regular grids, based on the polar stereographic coordinate system and a regular spacing, are handled by InterpFromMeshToGridx.
        // This will likely change in the future with improved interpolation schemes for non-regular grids.
        if ( M_grid.loaded )
            InterpFromMeshToMesh2dx(&interp_out,
                                  &mesh.indexTr()[0],&mesh.coordX()[0],&mesh.coordY()[0],
                                  mesh.numNodes(),mesh.numTriangles(),
                                  &interp_in[0],
                                  source_size,nb_var,
                                  &M_grid.gridX[0],&M_grid.gridY[0],M_grid_size,
                                  true, M_miss_val);
        else if ( M_ncols>0 & M_nrows>0 & M_mooring_spacing>0 )
            InterpFromMeshToGridx(interp_out,
                                  &mesh.indexTr()[0],&mesh.coordX()[0],&mesh.coordY()[0],
                                  mesh.numNodes(),mesh.numTriangles(),
                                  &interp_in[0],
                                  source_size, nb_var,
                                  xmin,ymax,
                                  M_mooring_spacing,M_mooring_spacing,
                                  M_nrows, M_ncols,
                                  M_miss_val);
        else
            std::logic_error("GridOutput::updateGridMeanWorker: No grid loaded and one of M_ncols, M_nrows, or M_mooring_spacing not set properly.");

        // Add the output pointer value to the grid vectors
        for (int i=0; i<nb_var; i++)
            for (int j=0; j<M_grid_size; ++j)
                variables[i].data_grid[j] += interp_out[nb_var*j+i];

        xDelete<double>(interp_out);
    }

    // Rotate the vectors as needed
    void GridOutput::rotateVectors(GmshMesh const &mesh, Vectorial_Variable const &vectorial_variable, std::vector<Variable> &variables)
    {
        // First we decide the rotation angle
        // Get the rotation of the neXtSIM grid
        mapx_class *mapNextsim;
        std::string configfileNextsim = (boost::format( "%1%/%2%/%3%" )
                                  % Environment::nextsimDir().string()
                                  % "data"
                                  % "NpsNextsim.mpp"
                                  ).str();

        std::vector<char> strNextsim(configfileNextsim.begin(), configfileNextsim.end());
        strNextsim.push_back('\0');
        mapNextsim = init_mapx(&strNextsim[0]);

        // Try to get the rotation of the data set
        mapx_class *map;
        double rotation_angle;
        if(!vectorial_variable.east_west_oriented & M_grid.mpp_file!="")
        {
            std::string configfile = (boost::format( "%1%/%2%/%3%" )
                                  % Environment::nextsimDir().string()
                                  % M_grid.dirname
                                  % M_grid.mpp_file
                                  ).str();

            std::vector<char> str(configfile.begin(), configfile.end());
            str.push_back('\0');
            map = init_mapx(&str[0]);
            rotation_angle = (mapNextsim->rotation-map->rotation)*PI/180.;
        }
        else if (vectorial_variable.east_west_oriented)
        {
            // or rotate to zonal/meridional
            rotation_angle = mapNextsim->rotation*PI/180.;
        }
        else
        {
            // or do nothing
            rotation_angle=0.;
        }

        // Rotate!
        if ( rotation_angle!=0. | vectorial_variable.east_west_oriented )
        {
            double cosang = std::cos(rotation_angle);
            double sinang = std::sin(rotation_angle);
            std::vector<double> lon = mesh.lon();

            for (int i=0; i<mesh.numNodes(); ++i)
            {
                if (vectorial_variable.east_west_oriented)
                {
                    cosang = std::cos(-lon[i]*PI/180+rotation_angle);
                    sinang = std::sin(-lon[i]*PI/180+rotation_angle);
                }

                // Calculate u and v
                double ull = cosang*variables[vectorial_variable.components_Id[0]].data_mesh[i] - sinang*variables[vectorial_variable.components_Id[1]].data_mesh[i];
                double vll = sinang*variables[vectorial_variable.components_Id[0]].data_mesh[i] + cosang*variables[vectorial_variable.components_Id[1]].data_mesh[i];

                // Overwrite x and y with u and v
                variables[vectorial_variable.components_Id[0]].data_mesh[i] = ull;
                variables[vectorial_variable.components_Id[1]].data_mesh[i] = vll;
            }
        }

    }

    // Set the _grid values back to zero
    void GridOutput::resetGridMean()
    {
        for (int i=0; i<M_nodal_variables.size(); i++)
            M_nodal_variables[i].data_grid.assign(M_grid_size, 0.);

        for (int i=0; i<M_elemental_variables.size(); i++)
            M_elemental_variables[i].data_grid.assign(M_grid_size, 0.);
    }

    // Set the _mesh values back to zero
    void GridOutput::resetMeshMean(GmshMesh const &mesh)
    {
        for (int i=0; i<M_nodal_variables.size(); i++)
            M_nodal_variables[i].data_mesh.assign(mesh.numNodes(), 0.);

        for (int i=0; i<M_elemental_variables.size(); i++)
            M_elemental_variables[i].data_mesh.assign(mesh.numTriangles(), 0.);
    }

    // Save the _grid values to file ... this needs a lot of work, NetCDF etc.
    void GridOutput::exportGridMeans(std::string postfix, double time_step, double mooring_output_time_step)
    {
        double time_factor = time_step/mooring_output_time_step;
        std::ofstream myfile;
        std::string filename;

        // Multiply all the grid values with 'time_factor'
        // Nodes
        for ( auto it=M_nodal_variables.begin(); it!=M_nodal_variables.end(); ++it )
        {
            for ( auto jt=it->data_grid.begin(); jt!=it->data_grid.end(); jt++ )
                // We need to filter out the missing values so they remain the same for all output files
                if ( *jt > M_miss_val )
                    *jt *= time_factor;
                else
                    *jt = M_miss_val;

            // Save to file - this is still just an ascii dump of one parameter!
            filename = it->name + postfix;
            std::cout << "Writing " << it->name << " to " << filename << std::endl;
            myfile.open(filename);
            std::copy(it->data_grid.begin(), it->data_grid.end(), ostream_iterator<float>(myfile, " "));
            myfile.close();
        }

        // Elements
        for ( auto it=M_elemental_variables.begin(); it!=M_elemental_variables.end(); ++it )
        {
            for ( auto jt=it->data_grid.begin(); jt!=it->data_grid.end(); jt++ )
                // We need to filter out the missing values so they remain the same for all output files
                if ( *jt > M_miss_val )
                    *jt *= time_factor;
                else
                    *jt = M_miss_val;

            // Save to file - this is still just an ascii dump of one parameter!
            filename = it->name + postfix;
            std::cout << "Writing " << it->name << " to " << filename << std::endl;
            myfile.open(filename);
            std::copy(it->data_grid.begin(), it->data_grid.end(), ostream_iterator<float>(myfile, " "));
            myfile.close();
        }

    }
}
