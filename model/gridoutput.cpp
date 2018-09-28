/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gridoutput.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Thu Aug  4 09:46:28 CEST 2016
 */

#include <gridoutput.hpp>
#include <date.hpp>

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
GridOutput::GridOutput(std::vector<Variable> variables, variableKind kind)
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
GridOutput::GridOutput(GmshMesh const& mesh, int ncols, int nrows, double mooring_spacing, double xmin, double ymin, std::vector<Variable> variables, variableKind kind)
    :
    GridOutput(variables, kind)
{
    this->initRegularGrid(ncols, nrows, mooring_spacing, xmin, ymin);
    this->resetMeshMean(mesh);
    this->initMask();
}

// Constructor for only one set of variables - arbitrary grid
GridOutput::GridOutput(GmshMesh const& mesh, Grid grid, std::vector<Variable> variables, variableKind kind)
    :
    GridOutput(variables, kind)
{
    this->initArbitraryGrid(grid);
    this->resetMeshMean(mesh);
    this->initMask();
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
GridOutput::GridOutput(GmshMesh const& mesh, int ncols, int nrows, double mooring_spacing,
                       double xmin, double ymin,
                       std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables)
    :
    GridOutput(nodal_variables, elemental_variables)
{
    this->initRegularGrid(ncols, nrows, mooring_spacing, xmin, ymin);
    this->resetMeshMean(mesh);
    this->initMask();
}

// constructor for nodal and elemental variables only (no vectors) - arbitrary grid
GridOutput::GridOutput(GmshMesh const& mesh, Grid grid, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables)
    :
    GridOutput(nodal_variables, elemental_variables)
{
    this->initArbitraryGrid(grid);
    this->resetMeshMean(mesh);
    this->initMask();
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
GridOutput::GridOutput(GmshMesh const& mesh, int ncols, int nrows, double mooring_spacing,
                       double xmin, double ymin,
                       std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables)
    :
    GridOutput(nodal_variables, elemental_variables, vectorial_variables)
{
    this->initRegularGrid(ncols, nrows, mooring_spacing, xmin, ymin);
    this->resetMeshMean(mesh);
    this->initMask();
}

// constructor for nodal, elemental and vectorial variables - arbitrary grid
GridOutput::GridOutput(GmshMesh const& mesh, Grid grid, std::vector<Variable> nodal_variables,
                       std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables)
    :
    GridOutput(nodal_variables, elemental_variables, vectorial_variables)
{
    this->initArbitraryGrid(grid);
    this->resetMeshMean(mesh);
    this->initMask();
}

GridOutput::~GridOutput()
{}

////////////////////////////////////////////////////////////////////////////////
// Initialisation routines
////////////////////////////////////////////////////////////////////////////////
void
GridOutput::initRegularGrid(int ncols, int nrows, double mooring_spacing, double xmin, double ymin)
{
    // Set the grid size
    M_ncols = ncols;
    M_nrows = nrows;
    M_mooring_spacing = mooring_spacing;
    M_grid_size = M_ncols*M_nrows;

    M_xmin = xmin;

    M_grid = Grid();
    M_grid.loaded = false;

    // Calculate lat and lon
    M_grid.gridLAT.assign(M_grid_size, 0.);
    M_grid.gridLON.assign(M_grid_size, 0.);

    mapx_class *map;
    std::string filename = Environment::nextsimMppfile();
    std::vector<char> str(filename.begin(), filename.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    double lat;
    double lon;
    int i=0;
    double X = xmin;
    for (int ncols=0; ncols<M_ncols; ncols++)
    {
        double Y = ymin;
        for (int nrows=0; nrows<M_nrows; nrows++)
        {
            int status = inverse_mapx(map,X,Y,&lat,&lon);
            M_grid.gridLAT[i] = lat;
            M_grid.gridLON[i] = lon;
            Y += mooring_spacing;
            i++;
        }
        X += mooring_spacing;
        M_ymax = Y - mooring_spacing;
    }

    close_mapx(map);

    M_lsm_nodes.assign(M_grid_size, 1.);
    M_lsm_elements.assign(M_grid_size, 1.);
    this->resetGridMean();
}

void
GridOutput::initArbitraryGrid(Grid grid)
{
    M_grid = grid;

    // Load the grid from file
    // Check file
    std::string filename = (boost::format( "%1%/%2%" )
                            % Environment::nextsimDataDir().string()
                            % M_grid.gridFile
                            ).str();
    if ( ! boost::filesystem::exists(filename) )
        throw std::runtime_error("File not found: " + filename);

    // Open file
	netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

    // Read the dimension of the grid
    netCDF::NcDim dim;
    dim = dataFile.getDim(M_grid.dimNameX);
    M_ncols = dim.getSize();

    dim = dataFile.getDim(M_grid.dimNameY);
    M_nrows = dim.getSize();
    M_mooring_spacing = -1.;
    M_grid_size = M_ncols*M_nrows;

    // Read the lat and lon
    netCDF::NcVar VLAT = dataFile.getVar(M_grid.latName);
    netCDF::NcVar VLON = dataFile.getVar(M_grid.lonName);

    // Read data into M_grid.gridLON & M_grid.gridLAT
    M_grid.gridLAT.resize(M_grid_size);
    M_grid.gridLON.resize(M_grid_size);
    VLAT.getVar(&M_grid.gridLAT[0]);
    VLON.getVar(&M_grid.gridLON[0]);

    // Calculate x and y
    M_grid.gridX.resize(M_grid_size);
    M_grid.gridY.resize(M_grid_size);

    mapx_class *map;
    filename = Environment::nextsimMppfile();
    std::vector<char> str(filename.begin(), filename.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    for (int i=0; i<M_grid_size; ++i)
        forward_mapx(map, M_grid.gridLAT[i], M_grid.gridLON[i], &M_grid.gridX[i], &M_grid.gridY[i]);

    close_mapx(map);

    M_grid.loaded = true;

    M_lsm_nodes.assign(M_grid_size, 1.);
    M_lsm_elements.assign(M_grid_size, 1.);
    this->resetGridMean();
}

void
GridOutput::initMask()
{
    M_proc_mask_indx = -1;
    M_ice_mask_indx = -1;

    for ( int i=0; i<M_elemental_variables.size(); ++i )
    {
        if (M_elemental_variables[i].varID==variableID::proc_mask) M_proc_mask_indx = i;
        if (M_elemental_variables[i].varID==variableID::ice_mask)  M_ice_mask_indx = i;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Functions other than construction and initialisation
////////////////////////////////////////////////////////////////////////////////

// Interpolate from the mesh values to the grid
void
GridOutput::updateGridMean(GmshMesh const& mesh)
{
    // Reset proc_mask
    if (M_proc_mask_indx != -1)
        M_elemental_variables[M_proc_mask_indx].data_grid.assign(M_grid_size, 0.);

    // Call the worker routine for the elements
    this->updateGridMeanWorker(&mesh.indexTr()[0], &mesh.coordX()[0], &mesh.coordY()[0], mesh.numNodes(), mesh.numTriangles(),
            mesh.numTriangles(), M_elemental_variables);

    // Rotate vectors if needed (these are assumed to be on the nodes)
    for ( auto it=M_vectorial_variables.begin(); it!=M_vectorial_variables.end(); it++ )
        if ( (it->east_west_oriented) || (M_grid.loaded) )
            this->rotateVectors(mesh, *it, M_nodal_variables);

    // Call the worker routine for the nodes
    if (M_proc_mask_indx != -1)
        this->updateGridMeanWorker(&mesh.indexTr()[0], &mesh.coordX()[0], &mesh.coordY()[0], mesh.numNodes(), mesh.numTriangles(),
                mesh.numNodes(), M_nodal_variables, true, M_elemental_variables[M_proc_mask_indx]);
    else
        if ( M_nodal_variables.size() > 0 )
            throw std::logic_error("GridOutput::updateGridMean: There are nodal variables to be interpolated but no proc_mask set.");

    // Mask if that's requested
    for ( auto it=M_nodal_variables.begin(); it!=M_nodal_variables.end(); it++ )
        if ( it->mask )
            for ( int i=0; i<M_grid_size; ++i )
                if ( M_elemental_variables[M_ice_mask_indx].data_grid[i] <= 0. && it->data_grid[i] != M_miss_val )
                    it->data_grid[i] = 0.;

    for ( auto it=M_elemental_variables.begin(); it!=M_elemental_variables.end(); it++ )
        if ( it->mask )
            for ( int i=0; i<M_grid_size; ++i )
                if ( M_elemental_variables[M_ice_mask_indx].data_grid[i] <= 0. && it->data_grid[i] != M_miss_val )
                    it->data_grid[i] = 0.;
}

// Interpolate from the mesh to the grid - updateing the gridded mean
void
GridOutput::updateGridMeanWorker(int* indexTr, double* coordX, double* coordY, int numNodes, int numTriangles,
        int source_size, std::vector<Variable>& variables)
{
    bool apply_mask = false;
    Variable mask(variableID::proc_mask);

    updateGridMeanWorker(indexTr, coordX, coordY, numNodes, numTriangles, source_size, variables, apply_mask, mask);
}

void
GridOutput::updateGridMeanWorker(int* indexTr, double* coordX, double* coordY, int numNodes, int numTriangles,
        int source_size, std::vector<Variable>& variables, bool apply_mask, Variable mask)
{
    int nb_var = variables.size();
    if (nb_var==0) exit;

    // TODO: We should check and make sure all variables in the vector are the same size
    // Input vector and output pointer
    std::vector<double> interp_in(nb_var*source_size);
    double* interp_out;

    // Stuff the input vector
    for (int i=0; i<source_size; ++i)
        for (int j=0; j<nb_var; j++)
            interp_in[nb_var*i+j] = variables[j].data_mesh[i];

    // At the moment a non-regular grid, loaded into M_grid is handled by InterpFromMeshToMesh2dx.
    // Regular grids, based on the polar stereographic coordinate system and a regular spacing, are handled by InterpFromMeshToGridx.
    // This will likely change in the future with improved interpolation schemes for non-regular grids.
    if ( M_grid.loaded )
    {
        InterpFromMeshToMesh2dx(&interp_out,
                                indexTr,coordX,coordY,
                                numNodes,numTriangles,
                                &interp_in[0],
                                source_size,nb_var,
                                &M_grid.gridX[0],&M_grid.gridY[0],M_grid_size,
                                true, 0.);
    }
    else if ( (M_ncols>0) && (M_nrows>0) && (M_mooring_spacing>0) )
    {
        InterpFromMeshToGridx(interp_out,
                              indexTr,coordX,coordY,
                              numNodes,numTriangles,
                              &interp_in[0],
                              source_size, nb_var,
                              M_xmin,M_ymax,
                              M_mooring_spacing,M_mooring_spacing,
                              M_ncols, M_nrows,
                              0.);
    }
    else
    {
        std::logic_error("GridOutput::updateGridMeanWorker: No grid loaded from file and one of M_ncols, M_nrows, or M_mooring_spacing not set properly.");
    }

    // Add the output pointer value to the grid vectors
    for (int i=0; i<nb_var; i++)
        for (int j=0; j<M_grid_size; ++j)
                variables[i].data_grid[j] += apply_mask ? (interp_out[nb_var*j+i]*mask.data_grid[j]) : interp_out[nb_var*j+i];

    xDelete<double>(interp_out);
}

// Set the land-sea mask
void
GridOutput::setLSM(GmshMeshSeq const& mesh)
{
    M_lsm_nodes = getMask(mesh, variableKind::nodal);
    M_lsm_elements = getMask(mesh, variableKind::elemental);
    this->resetGridMean();
}

// Return a mask
std::vector<double>
GridOutput::getMask(GmshMeshSeq const& mesh, variableKind kind)
{
    double source_size;
    switch (kind)
    {
        case variableKind::nodal:
            source_size = mesh.numNodes();
            break;

        case variableKind::elemental:
            source_size = mesh.numTriangles();
            break;

        default:
            std::logic_error("Incorrect variable kind in GridOutput::getMask");
    }

    // Call the worker routine using a vector of ones and give zero for missing values (land mask)
    std::vector<double> data_mesh(source_size, 1.);
    std::vector<double> data_grid(M_grid_size);

    Variable lsm(variableID::lsm);
    lsm.data_mesh = data_mesh;
    lsm.data_grid = data_grid;

    std::vector<Variable> variables(1);
    variables[0] = lsm;
    this->updateGridMeanWorker(&mesh.indexTr()[0], &mesh.coordX()[0], &mesh.coordY()[0], mesh.numNodes(), mesh.numTriangles(),
            source_size, variables);

    return variables[0].data_grid;
}

// Rotate the vectors as needed
void
GridOutput::rotateVectors(GmshMesh const& mesh, Vectorial_Variable const& vectorial_variable, std::vector<Variable>& variables)
{
    // First we decide the rotation angle
    // Get the rotation of the neXtSIM grid
    mapx_class *mapNextsim;
    std::string configfileNextsim = Environment::nextsimMppfile();
    std::vector<char> strNextsim(configfileNextsim.begin(), configfileNextsim.end());
    strNextsim.push_back('\0');
    mapNextsim = init_mapx(&strNextsim[0]);

    // Try to get the rotation of the data set
    double rotation_angle;
    if((!vectorial_variable.east_west_oriented) && (M_grid.mpp_file!=""))
    {
        mapx_class *map;
        std::string configfile = (boost::format( "%1%/%2%" )
            % Environment::nextsimMeshDir().string()
            % M_grid.mpp_file
            ).str();

        std::vector<char> str(configfile.begin(), configfile.end());
        str.push_back('\0');
        map = init_mapx(&str[0]);
        rotation_angle = (mapNextsim->rotation-map->rotation)*PI/180.;
        close_mapx(map);
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
    close_mapx(mapNextsim);

    // Rotate!
    if ( (rotation_angle!=0.) || (vectorial_variable.east_west_oriented) )
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

// Set the _grid values back to zero - with land mask
void
GridOutput::resetGridMean()
{
    for (int i=0; i<M_nodal_variables.size(); i++)
    {
        M_nodal_variables[i].data_grid.resize(M_grid_size);
        for (int j=0; j<M_grid_size; j++)
            if (M_lsm_nodes[j] == 0.)
                M_nodal_variables[i].data_grid[j] = M_miss_val;
            else
                M_nodal_variables[i].data_grid[j] = 0.;
    }

    for (int i=0; i<M_elemental_variables.size(); i++)
    {
        M_elemental_variables[i].data_grid.resize(M_grid_size);
        for (int j=0; j<M_grid_size; j++)
            if (M_lsm_elements[j] == 0.)
                M_elemental_variables[i].data_grid[j] = M_miss_val;
            else
                M_elemental_variables[i].data_grid[j] = 0.;
    }
}

// Set the _mesh values back to zero
void
GridOutput::resetMeshMean(GmshMesh const& mesh)
{
    for (int i=0; i<M_nodal_variables.size(); i++)
        M_nodal_variables[i].data_mesh.assign(mesh.numNodes(), 0.);

    for (int i=0; i<M_elemental_variables.size(); i++)
        M_elemental_variables[i].data_mesh.assign(mesh.numTriangles(), 0.);
}

// Initialise a netCDF file and return the file name in an std::string
std::string
GridOutput::initNetCDF(std::string file_prefix, fileLength file_length, double current_time)
{
    // Choose the right file name, depending on how much data goes in there
    boost::gregorian::date now = Nextsim::parse_date(current_time);
    std::stringstream filename;

    filename << file_prefix;
    switch (file_length)
    {
       case fileLength::daily:
           filename << "_" << now.year() << "d" << setw(3) << setfill('0') << now.day_of_year();
           break;
       case fileLength::weekly:
           // The last week of the year is troublesome!
           // The boost library will (sometimes) give this the number 1, even though week 1 should be in January
           if ( (now.month().as_number() == 12) && (now.week_number() == 1) )
           {
               filename << "_" << now.year() + 1 << "w" << setw(2) << setfill('0') << now.week_number();
           }
           else
           {
               filename << "_" << now.year() << "w" << setw(2) << setfill('0') << now.week_number();
           }
           break;
       case fileLength::monthly:
           filename << "_" << now.year() << "m" << setw(2) << setfill('0') << now.month().as_number();
           break;
       case fileLength::yearly:
           filename << "_" << now.year();
           break;
       case fileLength::inf:
           break;
       default:
           throw std::logic_error("invalid file length");
    }
    filename << ".nc";

    // Create the netCDF file.
    //LOG(DEBUG) <<"Initialise mooring file named " << filename.str() << "\n";
    netCDF::NcFile dataFile(filename.str(), netCDF::NcFile::replace);

    // Create the time dimension
    netCDF::NcDim tDim = dataFile.addDim("time"); // unlimited

    // Create the time variable
    netCDF::NcVar time = dataFile.addVar("time", netCDF::ncFloat, tDim);
    time.putAtt("standard_name","time");
    time.putAtt("long_name","simulation time");
    time.putAtt("units","days since 1900-01-01 00:00:00");
    time.putAtt("calendar","standard");

    M_nc_step=0;

    // Create the two spatial dimensions.
    netCDF::NcDim xDim = dataFile.addDim("x", M_ncols);
    netCDF::NcDim yDim = dataFile.addDim("y", M_nrows);

    std::vector<netCDF::NcDim> dims2(2);
    dims2[0] = xDim;
    dims2[1] = yDim;

    //cell methods - combine time method with hard-coded area method defined for each variable
    std::string cell_methods_time = "time: point ";//for snapshot
    if (Environment::vm()["moorings.snapshot"].as<bool>())
    {
        double averaging_period = Environment::vm()["moorings.output_timestep"].as<double>()*24.;//hours
        cell_methods_time = (boost::format( "time: mean (interval: %1% hours) " )
                               % averaging_period
                               ).str();
    }

    // Create the longitude and latitude variables
    // Longitude
    netCDF::NcVar lon = dataFile.addVar("longitude", netCDF::ncFloat, dims2);
    lon.putAtt("standard_name","longitude");
    lon.putAtt("long_name","longitude");
    lon.putAtt("units","degrees_north");
    lon.putAtt("_CoordinateAxisType","Lon");
    lon.putVar(&M_grid.gridLON[0]);

    // Latitude
    netCDF::NcVar lat = dataFile.addVar("latitude", netCDF::ncFloat, dims2);
    lat.putAtt("standard_name","latitude");
    lat.putAtt("long_name","latitude");
    lat.putAtt("units","degrees_north");
    lat.putAtt("_CoordinateAxisType","Lat");
    lat.putVar(&M_grid.gridLAT[0]);

    // Create the output variables
    netCDF::NcVar data;
    std::vector<netCDF::NcDim> dims(3);
    dims[0] = tDim;
    dims[1] = xDim;
    dims[2] = yDim;
    for (auto it=M_nodal_variables.begin(); it!=M_nodal_variables.end(); ++it)
    {
        if ( it->varID < 0 ) // Skip non-outputting variables
            continue;
        data = dataFile.addVar(it->name, netCDF::ncFloat, dims);
        data.putAtt("standard_name",it->stdName);
        data.putAtt("long_name",it->longName);
        data.putAtt("units",it->Units);
        data.putAtt("cell_methods", cell_methods_time + it->cell_methods);
        data.putAtt("_FillValue", netCDF::ncFloat, M_miss_val);
    }
    for (auto it=M_elemental_variables.begin(); it!=M_elemental_variables.end(); ++it)
    {
        if ( it->varID < 0 ) // Skip non-outputing variables
            continue;
        data = dataFile.addVar(it->name, netCDF::ncFloat, dims);
        data.putAtt("standard_name",it->stdName);
        data.putAtt("long_name",it->longName);
        data.putAtt("units",it->Units);
        data.putAtt("cell_methods", cell_methods_time + it->cell_methods);
        data.putAtt("_FillValue", netCDF::ncFloat, M_miss_val);
    }

    // Create the projection variable
    // - get the projection
    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    std::vector<char> str(mppfile.begin(), mppfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    // - determine if false easting (if there are no vectorial variables it doesn't matter)
    std::string false_easting = "1";
    for (auto it=M_vectorial_variables.begin(); it!=M_vectorial_variables.end(); ++it)
        if(it->east_west_oriented)
            false_easting = "0";
    
    // - add the projection variable
    std::vector<netCDF::NcDim> dims0(0);
    netCDF::NcVar proj = dataFile.addVar("Polar_Stereographic_Grid", netCDF::ncInt, dims0);
    proj.putAtt("grid_mapping_name", "polar_stereographic");
    proj.putAtt("false_easting", false_easting);
    proj.putAtt("false_northing", false_easting);
    proj.putAtt("semi_major_axis", std::to_string(map->equatorial_radius));
    proj.putAtt("semi_minor_axis", std::to_string(map->polar_radius));
    proj.putAtt("straight_vertical_longitude_from_pole", std::to_string(map->rotation));
	proj.putAtt("latitude_of_projection_origin", std::to_string(map->lat0));
    proj.putAtt("standard_parallel", std::to_string(map->lat1));
    close_mapx(map);


    dataFile.putAtt("Conventions", "CF-1.6");
    dataFile.putAtt("institution", "NERSC, Thormoehlens gate 47, N-5006 Bergen, Norway");
    dataFile.putAtt("source", "neXtSIM model fields");

    return filename.str();
}

// Write data to the netCDF file
void
GridOutput::appendNetCDF(std::string filename, double timestamp)
{
    // Open the netCDF file
    netCDF::NcFile dataFile(filename, netCDF::NcFile::write);

    // Append to time
    std::vector<size_t> start;
    start.push_back(M_nc_step);

    std::vector<size_t> count;
    count.push_back(1);

    M_nc_step++;

    netCDF::NcVar time = dataFile.getVar("time");
    time.putVar(start, count, &timestamp);

    // Append to the output variables
    start.push_back(0);
    start.push_back(0);

    count.push_back(M_ncols);
    count.push_back(M_nrows);

    // Save to file
    netCDF::NcVar data;
    for (auto it=M_nodal_variables.begin(); it!=M_nodal_variables.end(); ++it)
    {
        if ( it->varID < 0 ) // Skip non-outputing variables
            continue;
        data = dataFile.getVar(it->name);
        data.putVar(start, count, &it->data_grid[0]);
    }
    for (auto it=M_elemental_variables.begin(); it!=M_elemental_variables.end(); ++it)
    {
        if ( it->varID < 0 ) // Skip non-outputing variables
            continue;
        data = dataFile.getVar(it->name);
        data.putVar(start, count, &it->data_grid[0]);
    }
}
} // Nextsim
