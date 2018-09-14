/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:31:24 CET 2017
 */

#include <drifters.hpp>

/**
 * @class Drifters
 * @brief Manage drifters (base class)
 *
 * @see
 *
*/

namespace Nextsim
{
//! Constructors and destructors
//! * We have one constructor for regularly spaced drifters and one for
//!   drifter positions read in from file. In addition we have a private
//!   function (maskXY), called only by the constructors, to mask out drifters
//!   placed outside the ice cover.
//! * There is also a default, constructor which initialises things to zero and false.
 

Drifters::Drifters()
{
    M_is_initialised = false;
    M_no_drifters    = 0;
    M_X.resize(0);
    M_Y.resize(0);
    M_i.resize(0);
}
    
    
// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
Drifters::Drifters(double spacing, GmshMeshSeq const& mesh, std::vector<double>& conc, double climit)
{
    //! - 1) Calculates the grid spacing assuming a regular grid
    std::vector<double> RX = mesh.coordX();
    std::vector<double> RY = mesh.coordY();
    auto xcoords = std::minmax_element( RX.begin(), RX.end() );
    auto ycoords = std::minmax_element( RY.begin(), RY.end() );

    int ncols = (int) ( 0.5 + ( *xcoords.second - *xcoords.first )/spacing );
    int nrows = (int) ( 0.5 + ( *ycoords.second - *ycoords.first )/spacing );

    std::vector<double> X(ncols*nrows);
    std::vector<double> Y(ncols*nrows);

    int i=0;
    int y = *ycoords.first;
    for (int row=0; row<nrows; ++row)
    {
        int x = *xcoords.first;
        for (int col=0; col<ncols; ++col)
        {
            X[i] = x;
            Y[i] = y;
            ++i;
            x += spacing;
        }
        y += spacing;
    }

    //! -2) Applies the mask using conc and climit, and save to M_X and M_Y
    maskXY(mesh, X, Y, conc, climit);

    M_is_initialised = true;
}

    
// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
Drifters::Drifters(std::string filename, GmshMeshSeq const& mesh, std::vector<double>& conc, double climit, double current_time)
{
    //! - 1) Load the nodes from file
    // Check file
    std::fstream M_rgps_file;   // The file we read the IABP buoy data from

    M_rgps_file.open(filename, std::fstream::in);
    if ( ! M_rgps_file.good() )
        throw std::runtime_error("File not found: " + filename);

    //! -2) Read the current buoys from file
    int pos;    // To be able to rewind one line
    double time = current_time;
    std::vector<int> keepers;
    int gridSize=0;

    std::vector<double> LAT(0);
    std::vector<double> LON(0);

    int year, month, day, hour, number;
    double lat, lon;

    // Read the next line
    while ( M_rgps_file >> year >> month >> day >> hour >> number >> lat >> lon )
    {
        std::string date = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);
        time = from_date_string(date) + hour/24.;

        if(time== current_time )
        {
            LAT.push_back(lat);
            LON.push_back(lon);
            gridSize++;
        }
    }

    std::cout<<"gridSize: "<< gridSize <<"\n";

    //! - 3) Calculates cartesian (x,y) position from latitude and longitude
    std::vector<double> X(gridSize);
    std::vector<double> Y(gridSize);

    mapx_class *map;
    std::string mppfile = (boost::format( "%1%/%2%" )
            % Environment::nextsimMeshDir().string()
            % Environment::vm()["mesh.mppfile"].as<std::string>()
            ).str();
    std::vector<char> str(mppfile.begin(), mppfile.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    for (int i=0; i<gridSize; ++i)
        forward_mapx(map, LAT[i], LON[i], &X[i], &Y[i]);

    close_mapx(map);

    //! - 4) Applies mask using conc and climit, and save to M_X and M_Y
    maskXY(mesh, X, Y, conc, -1.);//climit);

    M_is_initialised = true;

    M_rgps_file.close();
}


// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
Drifters::Drifters(std::string gridFile,
                   std::string dimNameX, std::string dimNameY,
                   std::string latName, std::string lonName,
                   GmshMeshSeq const& mesh, std::vector<double>& conc, double climit)
{

    //! - 1) Loads the grid from file
    // Check file
    // TODO dirname not used
    std::string filename = (boost::format( "%1%/%2%" )
                            % Environment::nextsimDataDir().string()
                            % gridFile
                            ).str();

    if ( ! boost::filesystem::exists(filename) )
        throw std::runtime_error("Drifters::Drifters: File not found: " + filename);

    // Open file
    netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

    //! - 2) Reads the dimension of the grid
    netCDF::NcDim dim;
    dim = dataFile.getDim(dimNameX);
    if ( dim.isNull() )
        throw std::runtime_error("Drifters::Drifters: Empty dimension: " + dimNameX + " in file: " + filename);
    int ncols = dim.getSize();

    dim = dataFile.getDim(dimNameY);
    if ( dim.isNull() )
        throw std::runtime_error("Drifters::Drifters: Empty dimension: " + dimNameY + " in file: " + filename);
    int nrows = dim.getSize();

    int gridSize = ncols*nrows;

    //! - 3) Reads the latitude and longitude
    netCDF::NcVar VLAT = dataFile.getVar(latName);
    netCDF::NcVar VLON = dataFile.getVar(lonName);

    //! - 4) Reads the data in LON & LAT
    std::vector<double> LAT(gridSize);
    std::vector<double> LON(gridSize);
    VLAT.getVar(&LAT[0]);
    VLON.getVar(&LON[0]);

    //! - 5) Calculates cartesian (x,y) coordinates from the longitude and latitude
    std::vector<double> X(gridSize);
    std::vector<double> Y(gridSize);

    mapx_class *map;
    std::string mppfile = (boost::format( "%1%/%2%" )
            % Environment::nextsimMeshDir().string()
            % Environment::vm()["mesh.mppfile"].as<std::string>()
            ).str();
    std::vector<char> str(mppfile.begin(), mppfile.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    for (int i=0; i<gridSize; ++i)
        forward_mapx(map, LAT[i], LON[i], &X[i], &Y[i]);

    close_mapx(map);

    //! - 6) Applies mask using conc and climit, and save to M_X and M_Y
    maskXY(mesh, X, Y, conc, climit);

    M_is_initialised = true;
}

    
// --------------------------------------------------------------------------------------
//! Masks out the initial X and Y values so we only have drifters where there is ice.
void
Drifters::maskXY(GmshMeshSeq const& mesh, std::vector<double>& X, std::vector<double>& Y, std::vector<double>& conc, double clim)
{
    //! - 1) Interpolates the concentration onto the drifter positions
    int gridSize = X.size();
    double* interp_drifter_out;
    InterpFromMeshToMesh2dx(&interp_drifter_out,
                            &mesh.indexTr()[0],&mesh.coordX()[0],&mesh.coordY()[0],
                            mesh.numNodes(),mesh.numTriangles(),
                            &conc[0],
                            mesh.numTriangles(),1,
                            &X[0],&Y[0],gridSize,
                            true, 0.);

    //! - 2) Adds drifter positions where conc > conc_lim
    M_X.resize(0); // this shouldn't be necessary!
    M_Y.resize(0);
    M_i.resize(0);

    for ( int i=0; i<gridSize; ++i )
    {
        if ( interp_drifter_out[i] > clim )
        {
            M_X.push_back(X[i]);
            M_Y.push_back(Y[i]);
            M_i.push_back(i);
        }
    }

    M_no_drifters = M_X.size();

    xDelete<double>(interp_drifter_out);
}


// --------------------------------------------------------------------------------------
//! Checks to see if the drifters are correctly initialized.
bool
Drifters::isInitialised()
{
    return M_is_initialised;
}

//! * There are two types of moving drifters. The "normal" case where you overwrite
//!   the old coordinates with the new ones and the case where you
//!   want to store the new coordinates in a seperate variable. This latter
//!   case is useful when writing to file.

//! - 1) Moves drifters, with outputs into cartesian coordinates (X and Y).
void
Drifters::move(GmshMeshSeq const& mesh, std::vector<double> const& UT, std::vector<double>& X, std::vector<double>& Y)
{
    // Do nothing if we don't have to
    int no_drifters = X.size();
    assert(no_drifters==Y.size());
    if ( no_drifters == 0 )
        return;

    // Interpolate the total displacement onto the drifter positions
    int nb_var=2;
    int numNodes = mesh.numNodes();
    std::vector<double> interp_drifter_in(nb_var*numNodes);

    // Interpolate the velocity
    for (int i=0; i<numNodes; ++i)
    {
        interp_drifter_in[nb_var*i]   = UT[i];
        interp_drifter_in[nb_var*i+1] = UT[i+numNodes];
    }

    double* interp_drifter_out;
    InterpFromMeshToMesh2dx(&interp_drifter_out,
                            &mesh.indexTr()[0],&mesh.coordX()[0],&mesh.coordY()[0],
                            numNodes,mesh.numTriangles(),
                            &interp_drifter_in[0],
                            numNodes,nb_var,
                            &X[0],&Y[0],no_drifters,
                            true, 0.);

    for ( int i=0; i<no_drifters; ++i )
    {
        X[i] += interp_drifter_out[nb_var*i];
        Y[i] += interp_drifter_out[nb_var*i+1];
    }

    xDelete<double>(interp_drifter_out);
}

//! - 2) Moves drifters and replaces the old coordinates with the new ones.
void
Drifters::move(GmshMeshSeq const& mesh, std::vector<double> const& UT)
{
    this->move(mesh, UT, M_X, M_Y);
}


// File operations
// ----------------------------------------------------------------------------------------------
//! Initializes the netCDF file
//! * We use netCDF files to save a fixed number of drifters. The file
//! contains two vectors of latitude and longitude positions and can be
//! appended to through time. File name contains time stamp for the creation
//! date (format: yyyymmdd, we assume we won't be creating multiple files
//! with the same prefix per day.
 
void
Drifters::initNetCDF(std::string file_prefix, double current_time)
{
    // Construct the filename
    boost::gregorian::date now = Nextsim::parse_date(current_time);
    std::stringstream filename;

    filename << file_prefix;
    filename << now.year() << setw(2) << setfill('0') << now.month().as_number() << setw(2) << setfill('0') << now.day();
    filename << ".nc";
    M_filename = filename.str();

    /* Not sure if I want to check if the file exists
    // Throw an error if the file exists
    if ( ! boost::filesystem::exists(filename) )
    throw std::runtime_error("Drifters::initNetCDF: Trying to initialise a file that already exists: " + filename);
    */

    // Create the netCDF file.
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

    // Create the vector dimension
    netCDF::NcDim vecDim = dataFile.addDim("x", M_no_drifters);
    netCDF::NcDim emptyDim = dataFile.addDim("y", 1);
    std::vector<netCDF::NcDim> dims2(3);
    dims2[0] = tDim;
    dims2[1] = vecDim;
    dims2[2] = emptyDim;

    // Create the output variables
    // Longitude
    netCDF::NcVar lon = dataFile.addVar("longitude", netCDF::ncFloat, dims2);
    lon.putAtt("standard_name","longitude");
    lon.putAtt("long_name","longitude");
    lon.putAtt("units","degrees_north");

    // Latitude
    netCDF::NcVar lat = dataFile.addVar("latitude", netCDF::ncFloat, dims2);
    lat.putAtt("standard_name","latitude");
    lat.putAtt("long_name","latitude");
    lat.putAtt("units","degrees_north");

    // Reference indices if comming from a grid or a reference dataset
    netCDF::NcVar ind = dataFile.addVar("index", netCDF::ncInt, dims2);
    ind.putAtt("standard_name","index");
    ind.putAtt("long_name","reference index");
    ind.putAtt("units","");

    dataFile.putAtt("Conventions", "CF-1.6");
    dataFile.putAtt("institution", "NERSC, Thormoehlens gate 47, N-5006 Bergen, Norway");
    dataFile.putAtt("source", "neXtSIM model fields");
}


// -------------------------------------------------------------------------------------
//! Writes data to the netCDF file.
void
Drifters::appendNetCDF(double current_time, GmshMeshSeq const& mesh, std::vector<double> const& UT)
{
    // Move the drifters before export, but save the result in a different variable
    std::vector<double> X(M_no_drifters);
    std::vector<double> Y(M_no_drifters);
    this->move(mesh, UT, X, Y);

    // Calculate lat and lon
    mapx_class *map;
    std::string mppfile = (boost::format( "%1%/%2%" )
            % Environment::nextsimMeshDir().string()
            % Environment::vm()["mesh.mppfile"].as<std::string>()
            ).str();
    std::vector<char> str(mppfile.begin(), mppfile.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    std::vector<double> lat(M_no_drifters);
    std::vector<double> lon(M_no_drifters);
    for (int i=0; i<M_no_drifters; ++i)
        inverse_mapx(map, X[i], Y[i], &lat[i], &lon[i]);

    close_mapx(map);

    // Open the netCDF file
    netCDF::NcFile dataFile(M_filename, netCDF::NcFile::write);

    // Append to time
    std::vector<size_t> start;
    start.push_back(M_nc_step);

    std::vector<size_t> count;
    count.push_back(1);

    M_nc_step++;

    netCDF::NcVar time = dataFile.getVar("time");
    time.putVar(start, count, &current_time);

    // Append to the output variables
    start.push_back(0);
    count.push_back(M_no_drifters);

    start.push_back(0);
    count.push_back(1);

    netCDF::NcVar data;

    data = dataFile.getVar("longitude");
    data.putVar(start, count, &lon[0]);

    data = dataFile.getVar("latitude");
    data.putVar(start, count, &lat[0]);

    data = dataFile.getVar("index");
    data.putVar(start, count, &M_i[0]);
}
} // Nextsim
