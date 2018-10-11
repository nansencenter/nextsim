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
GridOutput::GridOutput(BamgMesh* bamgmesh, Grid grid, std::vector<Variable> variables, variableKind kind,
        BamgMesh* bamgmesh_root, bimap_type const & transfer_map, Communicator const & comm)
    :
    GridOutput(variables, kind)
{
    this->initArbitraryGrid(bamgmesh, grid, comm, bamgmesh_root, transfer_map);
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
GridOutput::GridOutput(BamgMesh* bamgmesh, Grid grid, std::vector<Variable> variables, variableKind kind, std::vector<Vectorial_Variable> vectorial_variables,
        BamgMesh* bamgmesh_root, bimap_type const & transfer_map, Communicator const & comm)
    :
    GridOutput(variables, kind)
{
    M_vectorial_variables = vectorial_variables;
    this->initArbitraryGrid(bamgmesh, grid, comm, bamgmesh_root, transfer_map);
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
GridOutput::GridOutput(BamgMesh* bamgmesh, Grid grid, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables,
        BamgMesh* bamgmesh_root, bimap_type const & transfer_map, Communicator const & comm)
    :
    GridOutput(nodal_variables, elemental_variables)
{
    this->initArbitraryGrid(bamgmesh, grid, comm, bamgmesh_root, transfer_map);
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
                       std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables,
        BamgMesh* bamgmesh_root, bimap_type const & transfer_map, Communicator const & comm)
    :
    GridOutput(nodal_variables, elemental_variables, vectorial_variables)
{
    this->initArbitraryGrid(bamgmesh, grid, comm, bamgmesh_root, transfer_map);
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
    std::string filename = (boost::format( "%1%/%2%" )
            % Environment::nextsimMeshDir().string()
            % Environment::vm()["mesh.mppfile"].as<std::string>()
            ).str();
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

    M_use_lsm = false;
    M_lsm.assign(M_grid_size, 1);
    this->initMask();
    this->resetGridMean();
    this->resetMeshMean(bamgmesh, true);
}

void
GridOutput::initArbitraryGrid(BamgMesh* bamgmesh, Grid& grid, Communicator const & comm,
        BamgMesh* bamgmesh_root, bimap_type const & transfer_map)
{
    M_grid = grid;
    M_comm = comm;

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

    // Do a conservative remapping if we find both cornerLatName and cornerLonName
    if ( M_grid.cornerLatName!= "" && M_grid.cornerLonName!= "")
    {
        netCDF::NcVar corner_x = dataFile.getVar("x_corners");
        netCDF::NcVar corner_y = dataFile.getVar("y_corners");
        netCDF::NcVar corner_z = dataFile.getVar("z_corners");

        std::vector<double> x(M_grid_size*4);
        std::vector<double> y(M_grid_size*4);
        std::vector<double> z(M_grid_size*4);

        corner_x.getVar(&x[0]);
        corner_y.getVar(&y[0]);
        corner_z.getVar(&z[0]);

        M_grid.gridCornerX.resize(4*M_grid_size);
        M_grid.gridCornerY.resize(4*M_grid_size);
        this->stereographicProjection(x, y, z, M_grid.gridCornerX, M_grid.gridCornerY);

        M_interp_method = interpMethod::conservative;
    }

    // Calculate x and y
    M_grid.gridX.resize(M_grid_size);
    M_grid.gridY.resize(M_grid_size);

    mapx_class *map;
    filename = (boost::format( "%1%/%2%" )
            % Environment::nextsimMeshDir().string()
            % Environment::vm()["mesh.mppfile"].as<std::string>()
            ).str();
    std::vector<char> str(filename.begin(), filename.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    for (int i=0; i<M_grid_size; ++i)
        forward_mapx(map, M_grid.gridLAT[i], M_grid.gridLON[i], &M_grid.gridX[i], &M_grid.gridY[i]);

    close_mapx(map);

    M_grid.loaded = true;

    M_use_lsm = false;
    M_lsm.assign(M_grid_size, 1);
    this->initMask();
    this->resetGridMean();
    if ( M_comm.rank() == 0 )
        this->resetMeshMean(bamgmesh, true, transfer_map, bamgmesh_root);
    else
        this->resetMeshMean(bamgmesh, true, transfer_map);
}

void
GridOutput::stereographicProjection(std::vector<double> const & x, std::vector<double> const & y, std::vector<double> const & z,
    std::vector<double>& ps_x, std::vector<double>& ps_y)
{
    // polar stereographic projection
    mapx_class *map;
    std::string filename = (boost::format( "%1%/%2%" )
            % Environment::nextsimMeshDir().string()
            % Environment::vm()["mesh.mppfile"].as<std::string>()
            ).str();
    std::vector<char> str(filename.begin(), filename.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    for ( int i=0; i<x.size(); i++ )
    {
        //compute latitude and longitude from cartesian coordinates
        double _x = x[i];
        double _y = y[i];
        double _z = z[i];

        // compute radius
        double radius = std::sqrt(std::pow(_x,2.)+std::pow(_y,2.)+std::pow(_z,2.));

        double latitude = std::asin(_z/radius)*(180./PI);
        double longitude = std::atan2(_y,_x);

        longitude = longitude-2*PI*std::floor(longitude/(2*PI));
        longitude = longitude*(180./PI);

        double x_, y_;
        int status = forward_mapx(map,latitude,longitude,&x_,&y_);

        ps_x[i] = x_;
        ps_y[i] = y_;
    }

    close_mapx(map);
}
void
GridOutput::initMask()
{
    M_ice_mask_indx = -1;

    for ( int i=0; i<M_elemental_variables.size(); ++i )
        if (M_elemental_variables[i].varID==variableID::ice_mask)
            M_ice_mask_indx = i;
}

void
GridOutput::setProcMask(BamgMesh* bamgmesh)
{
    // Call the worker routine using a vector of ones and give zero for missing values and gohsts
    std::vector<Variable> variables(1);
    variables[0] = Variable(variableID::proc_mask);

    variables[0].data_grid.assign(M_grid_size,0);
    variables[0].data_mesh.resize(bamgmesh->TrianglesSize[0]);

    std::fill( variables[0].data_mesh.begin(), variables[0].data_mesh.begin() + bamgmesh->TrianglesSize[2], 1. );
    std::fill( variables[0].data_mesh.begin() + bamgmesh->TrianglesSize[2],variables[0].data_mesh.end(),  0. );

    this->updateGridMeanWorker(bamgmesh, variableKind::elemental, interpMethod::meshToMesh, variables, 0.);

    M_proc_mask = variables[0].data_grid;
}

////////////////////////////////////////////////////////////////////////////////
// Functions other than construction and initialisation
////////////////////////////////////////////////////////////////////////////////

// Interpolate from the mesh values to the grid
void
GridOutput::updateGridMean(BamgMesh* bamgmesh)
{
    // Call the worker routine for the elements
    this->updateGridMeanWorker(bamgmesh, variableKind::elemental, M_interp_method, M_elemental_variables, M_miss_val);

    // Call the worker routine for the nodes
    this->updateGridMeanWorker(bamgmesh, variableKind::nodal, M_interp_method, M_nodal_variables, M_miss_val);

    // Ice mask if that's requested
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
GridOutput::updateGridMeanWorker(BamgMesh* bamgmesh, variableKind kind, interpMethod method, std::vector<Variable>& variables, double miss_val)
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

    /* Non-regular grids, loaded into M_grid is handled by InterpFromMeshToMesh2dx or
     * ConservativeRemappingMeshToGrid.
     * Regular grids, based on the polar stereographic coordinate system and a regular spacing, are
     * handled by InterpFromMeshToGridx.
     * Conservative remapping only really works when we can guarantee that the coastlines of mesh
     * and grid match. */
    if ( M_grid.loaded )
    {
        if ( kind==variableKind::elemental && method==interpMethod::conservative )
            ConservativeRemappingMeshToGrid(interp_out, interp_in,
                                    nb_var, M_grid_size, 0.,
                                    M_gridP, M_grid.gridCornerX, M_grid.gridCornerY,
                                    M_triangles, M_weights);
        else
            InterpFromMeshToMesh2dx(&interp_out,
                                    &indexTr[0],&coordX[0],&coordY[0],
                                    numNodes,numTriangles,
                                    &interp_in[0],
                                    source_size,nb_var,
                                    &M_grid.gridX[0],&M_grid.gridY[0],M_grid_size,
                                    true, 0.);
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
            if ( kind == variableKind::nodal )
                variables[i].data_grid[j] += interp_out[nb_var*j+i]*M_proc_mask[j];
            else
                variables[i].data_grid[j] += interp_out[nb_var*j+i];

    xDelete<double>(interp_out);
}

// Set the land-sea mask
/* This function should _only_ be called by the root and _only_ with bamgmesh_root as an argument.
 * You only need to call it if you want an LSM in your output or you need to get at it in the code
 * via GridOutput::getLSM() */
void
GridOutput::setLSM(BamgMesh* bamgmesh_root)
{
    M_use_lsm = true;

    // Call the worker routine using a vector of ones and give zero for missing values (land mask)
    std::vector<Variable> variables(1);
    variables[0] = Variable(variableID::lsm);
    variables[0].data_grid.assign(M_grid_size,0);
    variables[0].data_mesh.assign(bamgmesh_root->TrianglesSize[0], 1.);

    this->updateGridMeanWorker(bamgmesh_root, variableKind::elemental, interpMethod::meshToMesh, variables, 0.);

    M_lsm = std::vector<int>(variables[0].data_grid.begin(), variables[0].data_grid.end());
}

// Rotate the vectors as needed
void
GridOutput::rotateVectors(Vectorial_Variable const& vectorial_variable, int nb_var, double* &interp_out, double miss_val)
{
    // First we decide the rotation angle
    // Get the rotation of the neXtSIM grid
    mapx_class *mapNextsim;
    std::string configfileNextsim = (boost::format( "%1%/%2%" )
            % Environment::nextsimMeshDir().string()
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
    for ( auto &var : M_nodal_variables )
        var.data_grid.assign(M_grid_size, 0.);

    for ( auto &var : M_elemental_variables )
        var.data_grid.assign(M_grid_size, 0.);
}

// Apply the land-sea mask
void
GridOutput::applyLSM()
{
    // TODO: This is ugly - we should be able to loop over all the variables at once.
    for ( auto &var : M_nodal_variables )
        for (int i=0; i<M_grid_size; i++)
            if (M_lsm[i] == 0)
                var.data_grid[i] = M_miss_val;

    for ( auto &var : M_elemental_variables )
        for (int i=0; i<M_grid_size; i++)
            if (M_lsm[i] == 0)
                var.data_grid[i] = M_miss_val;
}

// Set the _mesh values back to zero, and recalculate weights and set proc_mask if needed
void
GridOutput::resetMeshMean(BamgMesh* bamgmesh,
        bool regrid, bimap_type const & transfer_map, BamgMesh* bamgmesh_root)
{
    for (int i=0; i<M_nodal_variables.size(); i++)
        M_nodal_variables[i].data_mesh.assign(bamgmesh->VerticesSize[0], 0.);

    for (int i=0; i<M_elemental_variables.size(); i++)
        M_elemental_variables[i].data_mesh.assign(bamgmesh->TrianglesSize[0], 0.);

    if ( regrid )
    {
        if ( M_nodal_variables.size() > 0 )
            this->setProcMask(bamgmesh);

        /* Calculate the weights on the root, broadcast them to ohers, and map from global to local
         * element id */
        if ( M_interp_method==interpMethod::conservative )
        {
            std::vector<int> gridP;
            std::vector<std::vector<int>> triangles;
            std::vector<std::vector<double>> weights;
            if ( M_comm.rank() == 0 )
                ConservativeRemappingWeights(bamgmesh_root,
                                    M_grid.gridX,M_grid.gridY,
                                    M_grid.gridCornerX,M_grid.gridCornerY,
                                    gridP, triangles, weights);

            this->broadcastWeights(gridP, triangles, weights);

            for ( int i=0; i<gridP.size(); i++ )
            {
                std::vector<int> local_triangles;
                std::vector<double> local_weights;
                // We'll almost always need either all the space or none of it, so reserving makes sense
                local_triangles.reserve(triangles[i].size());
                local_weights.reserve(weights[i].size());

                for ( int j=0; j<triangles[i].size(); j++ )
                {
                    int tr = triangles[i][j] + 1; // We need Bamg numbering
                    if (transfer_map.left.find(tr) != transfer_map.left.end())
                    {
                        local_triangles.push_back(transfer_map.left.find(tr)->second-1); // C++ numbering
                        local_weights.push_back(weights[i][j]);
                    }
                }

                if ( local_triangles.size() > 0 )
                {
                    M_gridP.push_back(gridP[i]);
                    M_triangles.push_back(local_triangles);
                    M_weights.push_back(local_weights);
                }
            }
        }
    }
}

// Broadcast gridP, triangles, and weights - only used when M_interp_method == interpMethod::conservative
void
GridOutput::broadcastWeights(std::vector<int>& gridP,
        std::vector<std::vector<int>>& triangles,
        std::vector<std::vector<double>>& weights )
{
    // Pack the triangles and weights
    std::vector<int> nb_weights;
    std::vector<int> packed_tr;
    std::vector<double> packed_w;
    if ( M_comm.rank() == 0 )
    {
        for ( int i=0; i<gridP.size(); i++ )
        {
            nb_weights.push_back(triangles[i].size());
            packed_tr.insert( packed_tr.end(), triangles[i].begin(), triangles[i].end() );
            packed_w.insert( packed_w.end(), weights[i].begin(), weights[i].end() );
        }
    }

    // Broadcast the sizes
    std::vector<int> tmp(2);
    tmp[0] = gridP.size();
    tmp[1] = packed_tr.size();

    boost::mpi::broadcast(M_comm, &tmp[0], 2, 0);

    int reduced_size = tmp[0];
    int packed_size = tmp[1];

    // Resize what's not on root and broadcast everything
    if ( M_comm.rank() != 0 )
    {
        gridP.resize(reduced_size);
        triangles.resize(reduced_size);
        weights.resize(reduced_size);

        nb_weights.resize(reduced_size);
        packed_tr.resize(packed_size);
        packed_w.resize(packed_size);
    }

    // TODO: Combine the first three broadcasts into one (they're all ints)
    boost::mpi::broadcast(M_comm, &gridP[0], reduced_size, 0);
    boost::mpi::broadcast(M_comm, &nb_weights[0], reduced_size, 0);
    boost::mpi::broadcast(M_comm, &packed_tr[0], packed_size, 0);
    boost::mpi::broadcast(M_comm, &packed_w[0], packed_size, 0);

    // Unpack
    if ( M_comm.rank() != 0 )
    {
        int first=0;
        int last=0;
        for ( int i=0; i<gridP.size(); i++ )
        {
            last += nb_weights[i];
            triangles[i] = std::vector<int>(packed_tr.begin()+first, packed_tr.begin()+last);
            weights[i] = std::vector<double>(packed_w.begin()+first, packed_w.begin()+last);
            first += nb_weights[i];
        }
    }
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

    // Apply mask - if needed
    if (M_use_lsm)
        this->applyLSM();

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
