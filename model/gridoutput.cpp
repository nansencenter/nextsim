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
GridOutput::GridOutput(BamgMesh* bamgmesh, int ncols, int nrows, double mooring_spacing, double xmin, double ymin, std::vector<Variable> variables, variableKind kind)
    :
    GridOutput(variables, kind)
{
    this->initRegularGrid(bamgmesh, ncols, nrows, mooring_spacing, xmin, ymin);
}

// Constructor for only one set of variables - arbitrary grid
GridOutput::GridOutput(BamgMesh* bamgmesh, Grid grid, std::vector<Variable> variables, variableKind kind)
    :
    GridOutput(variables, kind)
{
    this->initArbitraryGrid(bamgmesh, grid);
}

////////////////////////////////////////////////////////////////////////////////
// Constructor for one set of variables plus vectors
////////////////////////////////////////////////////////////////////////////////

// Constructor for one set of variables plus vectors - regular grid
GridOutput::GridOutput(BamgMesh* bamgmesh, int ncols, int nrows, double mooring_spacing, double xmin, double ymin, std::vector<Variable> variables, variableKind kind, std::vector<Vectorial_Variable> vectorial_variables)
    :
    GridOutput(variables, kind)
{
    M_vectorial_variables = vectorial_variables;
    this->initRegularGrid(bamgmesh, ncols, nrows, mooring_spacing, xmin, ymin);
}

// Constructor for only one set of variables - arbitrary grid
GridOutput::GridOutput(BamgMesh* bamgmesh, Grid grid, std::vector<Variable> variables, variableKind kind, std::vector<Vectorial_Variable> vectorial_variables)
    :
    GridOutput(variables, kind)
{
    M_vectorial_variables = vectorial_variables;
    this->initArbitraryGrid(bamgmesh, grid);
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
GridOutput::GridOutput(BamgMesh* bamgmesh, int ncols, int nrows, double mooring_spacing,
                       double xmin, double ymin,
                       std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables)
    :
    GridOutput(nodal_variables, elemental_variables)
{
    this->initRegularGrid(bamgmesh, ncols, nrows, mooring_spacing, xmin, ymin);
}

// constructor for nodal and elemental variables only (no vectors) - arbitrary grid
GridOutput::GridOutput(BamgMesh* bamgmesh, Grid grid, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables)
    :
    GridOutput(nodal_variables, elemental_variables)
{
    this->initArbitraryGrid(bamgmesh, grid);
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
GridOutput::GridOutput(BamgMesh* bamgmesh, int ncols, int nrows, double mooring_spacing,
                       double xmin, double ymin,
                       std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables)
    :
    GridOutput(nodal_variables, elemental_variables, vectorial_variables)
{
    this->initRegularGrid(bamgmesh, ncols, nrows, mooring_spacing, xmin, ymin);
}

// constructor for nodal, elemental and vectorial variables - arbitrary grid
GridOutput::GridOutput(BamgMesh* bamgmesh, Grid grid, std::vector<Variable> nodal_variables,
                       std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables)
    :
    GridOutput(nodal_variables, elemental_variables, vectorial_variables)
{
    this->initArbitraryGrid(bamgmesh, grid);
}

GridOutput::~GridOutput()
{}

////////////////////////////////////////////////////////////////////////////////
// Initialisation routines
////////////////////////////////////////////////////////////////////////////////
void
GridOutput::initRegularGrid(BamgMesh* bamgmesh, int ncols, int nrows, double mooring_spacing, double xmin, double ymin)
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
    std::string filename = Environment::nextsimDir().string() + "/data/" + Environment::vm()["mesh.mppfile"].as<std::string>();
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

    M_lsm.assign(M_grid_size, 1);
    this->resetGridMean();
    this->resetMeshMean(bamgmesh);
    this->initMask();
}

void
GridOutput::initArbitraryGrid(BamgMesh* bamgmesh, Grid grid)
{
    M_grid = grid;

    // Load the grid from file
    // Check file
    std::string filename = (boost::format( "%1%/%2%/%3%" )
                            % Environment::nextsimDir().string()
                            % M_grid.dirname
                            % M_grid.gridFile
                            ).str();
    if ( ! boost::filesystem::exists(filename) )
        throw std::runtime_error("File not found: " + filename);

    // Open file
    netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

    // If we don't know the dimNameX and dimNameY we assume the first dimension
    // is x and the second is y ... unless M_grid.transpose is true, then we do
    // it the other way around.
    netCDF::NcVar VLAT = dataFile.getVar(M_grid.latName);
    netCDF::NcVar VLON = dataFile.getVar(M_grid.lonName);

    // Read the dimension of the grid
    netCDF::NcDim dim;
    if ( M_grid.dimNameX == "" )
        M_grid.dimNameX = VLAT.getDim((int)M_grid.transpose).getName();

    dim = dataFile.getDim(M_grid.dimNameX);
    M_ncols = dim.getSize();

    if ( M_grid.dimNameY == "" )
        M_grid.dimNameY = VLAT.getDim((int)!M_grid.transpose).getName();

    dim = dataFile.getDim(M_grid.dimNameY);
    M_nrows = dim.getSize();

    M_mooring_spacing = -1.;
    M_grid_size = M_ncols*M_nrows;

    // Read the lat and lon, and theta and corners if requested
    // Read data into M_grid.gridLON & M_grid.gridLAT
    M_grid.gridLAT.resize(M_grid_size);
    M_grid.gridLON.resize(M_grid_size);
    VLAT.getVar(&M_grid.gridLAT[0]);
    VLON.getVar(&M_grid.gridLON[0]);
    if ( M_grid.thetaName != "" )
    {
        netCDF::NcVar Theta = dataFile.getVar(M_grid.thetaName);
        M_grid.gridTheta.resize(M_grid_size);
        Theta.getVar(&M_grid.gridTheta[0]);
    }
    if ( M_grid.cornerLatName!= "" && M_grid.cornerLonName!= "")
    {
        netCDF::NcVar cornerLat = dataFile.getVar(M_grid.cornerLatName);
        netCDF::NcVar cornerLon = dataFile.getVar(M_grid.cornerLonName);
        M_grid.gridCornerLat.resize(M_grid_size*4);
        M_grid.gridCornerLon.resize(M_grid_size*4);
        cornerLat.getVar(&M_grid.gridCornerLat[0]);
        cornerLon.getVar(&M_grid.gridCornerLon[0]);
    }

    // Calculate x and y
    M_grid.gridX.resize(M_grid_size);
    M_grid.gridY.resize(M_grid_size);

    mapx_class *map;
    filename = Environment::nextsimDir().string() + "/data/" + Environment::vm()["mesh.mppfile"].as<std::string>();
    std::vector<char> str(filename.begin(), filename.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    for (int i=0; i<M_grid_size; ++i)
        forward_mapx(map, M_grid.gridLAT[i], M_grid.gridLON[i], &M_grid.gridX[i], &M_grid.gridY[i]);

    if ( M_grid.cornerLatName!= "" && M_grid.cornerLonName!= "")
    {
        M_grid.gridCornerX.resize(4*M_grid_size);
        M_grid.gridCornerY.resize(4*M_grid_size);
        for (int i=0; i<4*M_grid_size; ++i)
            forward_mapx(map, M_grid.gridCornerLat[i], M_grid.gridCornerLon[i], &M_grid.gridCornerX[i], &M_grid.gridCornerY[i]);
    }

    close_mapx(map);

    M_grid.loaded = true;

    M_lsm.assign(M_grid_size, 1);
    this->resetGridMean();
    this->resetMeshMean(bamgmesh, true);
    this->initMask();
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
GridOutput::updateGridMean(BamgMesh* bamgmesh)
{
    // Reset proc_mask
    if (M_proc_mask_indx != -1)
        M_elemental_variables[M_proc_mask_indx].data_grid.assign(M_grid_size, 0.);

    // Call the worker routine for the elements
    this->updateGridMeanWorker(bamgmesh, variableKind::elemental, M_elemental_variables, M_miss_val);

    // Call the worker routine for the nodes
    if (M_proc_mask_indx != -1)
        this->updateGridMeanWorker(bamgmesh, variableKind::nodal, M_nodal_variables, M_miss_val,
                true, M_elemental_variables[M_proc_mask_indx]);
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

// Interpolate from the mesh to the grid - updating the gridded mean
void
GridOutput::updateGridMeanWorker(BamgMesh* bamgmesh, variableKind kind, std::vector<Variable>& variables, double miss_val)
{
    bool apply_mask = false;
    Variable mask(variableID::proc_mask);

    this->updateGridMeanWorker(bamgmesh, kind, variables, miss_val, apply_mask, mask);
}

void
GridOutput::updateGridMeanWorker(BamgMesh* bamgmesh, variableKind kind, std::vector<Variable>& variables, double miss_val,
        bool apply_mask, Variable mask)
{
    int nb_var = variables.size();
    if ( nb_var == 0 )
        return;

    // Copy the triangle information
    int numTriangles = bamgmesh->TrianglesSize[0]; // mesh.numTriangles();
    std::vector<int> indexTr(3*numTriangles); // = mesh.indexTr();
    for (int tr=0; tr<numTriangles; ++tr)
    {
        indexTr[3*tr  ] = bamgmesh->Triangles[4*tr];
        indexTr[3*tr+1] = bamgmesh->Triangles[4*tr+1];
        indexTr[3*tr+2] = bamgmesh->Triangles[4*tr+2];
    }

    // Copy the node information
    int numNodes     = bamgmesh->VerticesSize[0]; //mesh.numNodes();
    std::vector<double> coordX(numNodes); // = mesh.coordX();
    std::vector<double> coordY(numNodes); // = mesh.coordY();
    for (int id=0; id<numNodes; ++id)
    {
        coordX[id] = bamgmesh->Vertices[3*id];
        coordY[id] = bamgmesh->Vertices[3*id+1];
    }

    int source_size = (kind==variableKind::nodal) ? numNodes : numTriangles;

    // Input vector and output pointer
    std::vector<double> interp_in(nb_var*source_size);
    double* interp_out;


    // Stuff the input vector
    for (int j=0; j<nb_var; j++)
    {
        assert( variables[j].data_mesh.size() == source_size );
        for (int i=0; i<source_size; ++i)
            interp_in[nb_var*i+j] = variables[j].data_mesh[i];
    }

    // Non-regular grids, loaded into M_grid is handled by InterpFromMeshToMesh2dx or ConservativeRemappingMeshToGrid.
    // Regular grids, based on the polar stereographic coordinate system and a regular spacing, are handled by InterpFromMeshToGridx.
    // TODO: Permit regular grids to use the conservative remapping.
    if ( M_grid.loaded )
    {
        if (false) // kind==variableKind::elemental && M_grid.cornerLatName!="" && M_grid.cornerLonName!="" )
            ConservativeRemappingMeshToGrid(interp_out, interp_in,
                                    nb_var, M_grid_size, miss_val,
                                    M_gridP, M_triangles, M_weights);
        else
            InterpFromMeshToMesh2dx(&interp_out,
                                    &indexTr[0],&coordX[0],&coordY[0],
                                    numNodes,numTriangles,
                                    &interp_in[0],
                                    source_size,nb_var,
                                    &M_grid.gridX[0],&M_grid.gridY[0],M_grid_size,
                                    true, miss_val);
    }
    else if ( (M_ncols>0) && (M_nrows>0) && (M_mooring_spacing>0) )
    {
        InterpFromMeshToGridx(interp_out,
                              &indexTr[0],&coordX[0],&coordY[0],
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

    // Rotate vectors if needed (these may only be on the nodes)
    if ( kind == variableKind::nodal )
        for ( auto it=M_vectorial_variables.begin(); it!=M_vectorial_variables.end(); it++ )
            if ( it->orientation != vectorOrientation::neXtSIM )
                this->rotateVectors(*it, nb_var, interp_out, miss_val);


    // Add the output pointer value to the grid vectors
    for (int i=0; i<nb_var; i++)
        for (int j=0; j<M_grid_size; ++j)
                variables[i].data_grid[j] += apply_mask ? (interp_out[nb_var*j+i]*mask.data_grid[j]) : interp_out[nb_var*j+i];

    xDelete<double>(interp_out);
}

// Set the land-sea mask
void
GridOutput::setLSM(BamgMesh* bamgmesh)
{
    M_lsm = getMask(bamgmesh, variableKind::elemental);
    this->resetGridMean();
}

// Return a mask
std::vector<int>
GridOutput::getMask(BamgMesh* bamgmesh, variableKind kind)
{
    // Call the worker routine using a vector of ones and give zero for missing values (land mask)
    std::vector<Variable> variables(1);
    variables[0] = Variable(variableID::lsm);

    variables[0].data_grid.assign(M_grid_size,0);
    switch (kind)
    {
        case variableKind::nodal:
            throw std::logic_error("variableKind::nodal not supported in GridOutput::getMask");
            variables[0].data_mesh.assign(bamgmesh->VerticesSize[0], 1.);
            break;

        case variableKind::elemental:
            variables[0].data_mesh.assign(bamgmesh->TrianglesSize[0], 1.);
            break;

        default:
            throw std::logic_error("Incorrect variable kind in GridOutput::getMask");
    }

    this->updateGridMeanWorker(bamgmesh, kind, variables, 0.);

    return std::vector<int>(variables[0].data_grid.begin(), variables[0].data_grid.end());
}

// Rotate the vectors as needed
void
GridOutput::rotateVectors(Vectorial_Variable const& vectorial_variable, int nb_var, double* &interp_out, double miss_val)
{
    // First we decide the rotation angle
    // Get the rotation of the neXtSIM grid
    mapx_class *mapNextsim;
    std::string configfileNextsim = (boost::format( "%1%/%2%/%3%" )
                                     % Environment::nextsimDir().string()
                                     % "data"
                                     % Environment::vm()["mesh.mppfile"].as<std::string>()
                                     ).str();

    std::vector<char> strNextsim(configfileNextsim.begin(), configfileNextsim.end());
    strNextsim.push_back('\0');
    mapNextsim = init_mapx(&strNextsim[0]);

    double rotation_angle = mapNextsim->rotation*PI/180.;

    close_mapx(mapNextsim);

    // One last check before we start
    if ( (vectorial_variable.orientation == vectorOrientation::grid) && (M_grid.thetaName == "") )
        throw std::logic_error("GridOutput::rotateVectors: You asked for vectors aligned with the output grid, but did not provide the rotation angle. Check your inputs when constructing your GridOutput::Grid object, or make sure you use simul.moorings_grid_file if you're using moorins.");

    // Rotate!
    for (int i=0; i<M_grid_size; ++i)
    {
        int first  = i*nb_var + vectorial_variable.components_Id.first;
        int second = i*nb_var + vectorial_variable.components_Id.second;
        if ( interp_out[first] == miss_val )
            continue;

        double cosang;
        double sinang;

        switch (vectorial_variable.orientation)
        {
            case vectorOrientation::east_west:
                cosang = std::cos(rotation_angle - M_grid.gridLON[i]*PI/180);
                sinang = std::sin(rotation_angle - M_grid.gridLON[i]*PI/180);
                break;
            case vectorOrientation::grid:
                cosang = std::cos(rotation_angle - M_grid.gridTheta[i]);
                sinang = std::sin(rotation_angle - M_grid.gridTheta[i]);
                break;
        }

        // Calculate u and v
        double u = cosang*interp_out[first] - sinang*interp_out[second];
        double v = sinang*interp_out[first] + cosang*interp_out[second];
        // Overwrite x and y with u and v
        interp_out[first]  = u;
        interp_out[second] = v;
    }
}

// Set the _grid values back to zero - with land mask
void
GridOutput::resetGridMean()
{
    // TODO: This is ugly - we should be able to loop over all the variables at once.
    for (int i=0; i<M_nodal_variables.size(); i++)
    {
        M_nodal_variables[i].data_grid.resize(M_grid_size);
        for (int j=0; j<M_grid_size; j++)
            if (M_lsm[j] == 0)
                M_nodal_variables[i].data_grid[j] = M_miss_val;
            else
                M_nodal_variables[i].data_grid[j] = 0.;
    }

    for (int i=0; i<M_elemental_variables.size(); i++)
    {
        M_elemental_variables[i].data_grid.resize(M_grid_size);
        for (int j=0; j<M_grid_size; j++)
            if (M_lsm[j] == 0)
                M_elemental_variables[i].data_grid[j] = M_miss_val;
            else
                M_elemental_variables[i].data_grid[j] = 0.;
    }
}

// Set the _mesh values back to zero and recalculate weights if needed (and if they're being used)
void
GridOutput::resetMeshMean(BamgMesh* bamgmesh, bool regrid)
{
    for (int i=0; i<M_nodal_variables.size(); i++)
        M_nodal_variables[i].data_mesh.assign(bamgmesh->VerticesSize[0], 0.);

    for (int i=0; i<M_elemental_variables.size(); i++)
        M_elemental_variables[i].data_mesh.assign(bamgmesh->TrianglesSize[0], 0.);

    if ( regrid && (M_grid.cornerLatName != "" && M_grid.cornerLonName != "") )
        ConservativeRemappingWeights(bamgmesh,
                                M_grid.gridX,M_grid.gridY,
                                M_grid.gridCornerX,M_grid.gridCornerY,
                                M_gridP, M_triangles, M_weights);
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
        if ( it->varID < 0 ) // Skip non-outputing variables
            continue;
        data = dataFile.addVar(it->name, netCDF::ncFloat, dims);
        data.putAtt("standard_name",it->stdName);
        data.putAtt("long_name",it->longName);
        data.putAtt("units",it->Units);
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
        data.putAtt("_FillValue", netCDF::ncFloat, M_miss_val);
    }

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
