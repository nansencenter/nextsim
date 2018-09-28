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
    M_num_drifters   = 0;
    M_X.resize(0);
    M_Y.resize(0);
    M_i.resize(0);
}

// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! * defines a uniform grid with resolution defined by spacing
//! Called by FiniteElement::initEquallySpacedDrifters()
Drifters::Drifters(double const& spacing,
        GmshMeshSeq const& movedmesh,
        std::vector<double> & conc, double const& climit,
        double const& current_time, double const& output_freq)
{
    //! -1) Set the time and output freq
    M_time_init = current_time;
    M_output_freq = output_freq;

    //! - 2) Calculates the grid spacing assuming a regular grid
    std::vector<double> RX = movedmesh.coordX();
    std::vector<double> RY = movedmesh.coordY();
    auto xcoords = std::minmax_element( RX.begin(), RX.end() );
    auto ycoords = std::minmax_element( RY.begin(), RY.end() );

    int ncols = (int) ( 0.5 + ( *xcoords.second - *xcoords.first )/spacing );
    int nrows = (int) ( 0.5 + ( *ycoords.second - *ycoords.first )/spacing );

    std::vector<double> X(ncols*nrows);
    std::vector<double> Y(ncols*nrows);
    std::vector<int> INDS(ncols*nrows);

    int i=0;
    int y = *ycoords.first;
    for (int row=0; row<nrows; ++row)
    {
        int x = *xcoords.first;
        for (int col=0; col<ncols; ++col)
        {
            X[i] = x;
            Y[i] = y;
            INDS[i] = i;
            ++i;
            x += spacing;
        }
        y += spacing;
    }

    //! -3) Applies the mask using conc and climit, and save to M_X and M_Y
    this->maskXY(movedmesh, X, Y, INDS, conc, climit);

    M_is_initialised = true;
}


// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! Called by FiniteElement::initSidfexDrifters() and FiniteElement::initRGPSDrifters()
Drifters::Drifters(std::string const& filename,
        GmshMeshSeq const& movedmesh,
        std::vector<double> & conc, double const& climit,
        double const& current_time, double const& output_freq)
{
    // interface for RGPS, SIDFEX
    // - reads a text file
    // - NB usually doesn't work for true SIDFEX buoy IDs as they are too large
    //   - therefore use an index file to map the true IDs to smaller integers
    //     (this is what is done in the forecast system)

    //! -1) Set the time and output freq
    M_time_init = current_time;
    M_output_freq = output_freq;

    //! - 2) Load the nodes from file

    // Check file
    std::fstream drifter_text_file;   // The file we read the buoy data from

    drifter_text_file.open(filename, std::fstream::in);
    if ( ! drifter_text_file.good() )
        throw std::runtime_error("File not found: " + filename);

    //! -3) Read the current buoys from file
    double time = current_time;
    std::vector<int> keepers;

    std::vector<double> LAT(0);
    std::vector<double> LON(0);
    std::vector<int> INDS(0);

    // skip header line
    std::string header;
    std::getline(drifter_text_file, header);

    int year, month, day, hour;//, number;
    int number;
    double lat, lon;

    // Read the next line
    // NB need to use eof() and old_number instead of commented code since that version stopped
    // after the first line for one example of input file
    int gridSize = 0;
    int old_number = 0;
    //while ( M_rgps_file >> year >> month >> day >> hour >> number >> lat >> lon )
    while (!drifter_text_file.eof())
    {

        if (gridSize>0)
            old_number = number;

        drifter_text_file >> year >> month >> day >> hour >> number >> lat >> lon;
        std::string date = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);
        time = from_date_string(date) + hour/24.;

        // may not be eof if \n at end of file
        // - this can lead to repetition of the last line
        if (gridSize>0)
            if (number == old_number)
                break;

        if(time== current_time )
        {
            LAT.push_back(lat);
            LON.push_back(lon);
            INDS.push_back(number);
            gridSize++;

            //std::cout << year << ", "<< month << ", "<< day << ", "
            //    << hour << ", "  << number << ", "<< lat << ", "<< lon << "\n";
        }
        else if(time>current_time)
            break;
    }

    //! - 3) Calculates cartesian (x,y) position from latitude and longitude
    std::vector<double> X(gridSize);
    std::vector<double> Y(gridSize);

    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    std::vector<char> str(mppfile.begin(), mppfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    for (int i=0; i<gridSize; ++i)
        forward_mapx(map, LAT[i], LON[i], &X[i], &Y[i]);

    close_mapx(map);

    //! - 4) Applies mask using conc and climit, and save to M_X and M_Y
    // - also save indices (order drifters are read) to M_i 
    this->maskXY(movedmesh, X, Y, INDS, conc, climit);

    M_is_initialised = true;

    drifter_text_file.close();
}


// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! * reads a netcdf file
//! Called by FiniteElement::initOsisafDrifters()
Drifters::Drifters(std::string const& gridFile,
                   std::string const& dimNameX, std::string const& dimNameY,
                   std::string const& latName, std::string const& lonName,
                   GmshMeshSeq const& movedmesh,
                   std::vector<double> & conc, double const& climit,
                   double const& current_time, double const& output_freq)
{

    //! -1) Set the time and output freq
    M_time_init = current_time;
    M_output_freq = output_freq;

    //! - 2) Loads the grid from file
    // Check file
    std::string filename = (boost::format( "%1%/%2%" )
                            % Environment::nextsimDataDir().string()
                            % gridFile
                            ).str();

    if ( ! boost::filesystem::exists(filename) )
        throw std::runtime_error("Drifters::Drifters: File not found: " + filename);

    // Open file
    netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

    //! - 3) Reads the dimension of the grid
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

    //! - 4) Reads the latitude and longitude
    netCDF::NcVar VLAT = dataFile.getVar(latName);
    netCDF::NcVar VLON = dataFile.getVar(lonName);

    //! - 5) Reads the data in LON & LAT
    std::vector<double> LAT(gridSize);
    std::vector<double> LON(gridSize);
    VLAT.getVar(&LAT[0]);
    VLON.getVar(&LON[0]);

    //! - 6) Calculates cartesian (x,y) coordinates from the longitude and latitude
    std::vector<double> X(gridSize);
    std::vector<double> Y(gridSize);
    std::vector<int> INDS(gridSize);

    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    std::vector<char> str(mppfile.begin(), mppfile.end());
    str.push_back('\0');

    map = init_mapx(&str[0]);

    for (int i=0; i<gridSize; ++i)
    {
        forward_mapx(map, LAT[i], LON[i], &X[i], &Y[i]);
        INDS[i] = i;
    }

    close_mapx(map);

    //! - 7) Applies mask using conc and climit, and save to M_X and M_Y
    this->maskXY(movedmesh, X, Y, INDS, conc, climit);

    M_is_initialised = true;
}


// --------------------------------------------------------------------------------------
//! Determine if we need to output a drifter
//! Called by outputtingDrifters()
//  TODO can now make a loop over all the drifters
bool
Drifters::isOutputTime(double const& current_time)
{
    // can only output if it's initialised
    if(!M_is_initialised)
        return false;

    bool do_output = false;
    if(current_time>M_time_init)
        // output is already done at init time
        do_output = std::fmod(current_time - M_time_init, M_output_freq) == 0;
    return do_output;
}


// --------------------------------------------------------------------------------------
//! Masks out the initial X and Y values so we only have drifters where there is ice.
//! Also fills M_i with the indices that are kept
void
Drifters::maskXY(GmshMeshSeq const& movedmesh,
        std::vector<double> & X, std::vector<double> & Y,
        std::vector<int> const& INDS,
        std::vector<double> & conc, double const& clim)
{

    //! - 1) Interpolates the concentration onto the drifter positions
    int gridSize = X.size();
    double* interp_drifter_out;
    InterpFromMeshToMesh2dx(&interp_drifter_out,
                            &movedmesh.indexTr()[0],
                            &movedmesh.coordX()[0],
                            &movedmesh.coordY()[0],
                            movedmesh.numNodes(),
                            movedmesh.numTriangles(),
                            &conc[0],
                            movedmesh.numTriangles(),1,
                            &X[0],&Y[0],gridSize,
                            true, 0.);

    //! - 2) Adds drifter positions where conc > conc_lim
    M_X.resize(0); // this shouldn't be necessary!
    M_Y.resize(0);
    M_i.resize(0);
    M_conc.resize(0);

    for ( int i=0; i<gridSize; ++i )
    {
        if ( interp_drifter_out[i] > clim )
        {
            M_X.push_back(X[i]);
            M_Y.push_back(Y[i]);
            M_i.push_back(INDS[i]);
            M_conc.push_back(interp_drifter_out[i]);
        }
    }

    M_num_drifters = M_X.size();

    xDelete<double>(interp_drifter_out);
}//maskXY()


// --------------------------------------------------------------------------------------
//! Checks to see if the drifters are correctly initialized.
bool
Drifters::isInitialised()
{
    return M_is_initialised;
}//isInitialised


// --------------------------------------------------------------------------------------
//! Move drifters and replace the old coordinates with the new ones
//! called by FiniteElement::checkDrifters()
void
Drifters::move(GmshMeshSeq const& mesh,
        std::vector<double> const& UT)
{
    // Do nothing if we don't have to
    if ( M_num_drifters == 0 )
        return;

    // Interpolate the total displacement onto the drifter positions
    int nb_var=2;
    int numNodes = mesh.numNodes();
    std::vector<double> interp_drifter_in(nb_var*numNodes);
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
                            &M_X[0],&M_Y[0],M_num_drifters,
                            true, 0.);

    // Add the displacement to the current position
    for ( int i=0; i<M_num_drifters; ++i )
    {
        M_X[i] += interp_drifter_out[nb_var*i];
        M_Y[i] += interp_drifter_out[nb_var*i+1];
    }

    xDelete<double>(interp_drifter_out);
}


// --------------------------------------------------------------------------------------
//! interp conc onto drifter positions
//! called by FiniteElement::checkDrifters()
void
Drifters::updateConc(GmshMeshSeq const& movedmesh,
        std::vector<double> &conc)
{
    // Do nothing if we don't have to
    if ( M_num_drifters == 0 )
        return;

    // move the mesh before interpolating
    int numNodes = movedmesh.numNodes();
    int numElements = movedmesh.numTriangles();

    // Interpolate the concentration onto the drifter positions
    int nb_var=1;
    double* interp_drifter_out;
    InterpFromMeshToMesh2dx(&interp_drifter_out,
                            &movedmesh.indexTr()[0], &movedmesh.coordX()[0], &movedmesh.coordY()[0],
                            numNodes, numElements,
                            &conc[0],
                            numElements, nb_var,
                            &M_X[0], &M_Y[0], M_num_drifters,
                            true, 0.);

    for ( int i=0; i<M_num_drifters; ++i )
        M_conc[i] = interp_drifter_out[i];

    xDelete<double>(interp_drifter_out);
}//updateConc


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
    time.putAtt("standard_name", "time");
    time.putAtt("long_name", "simulation time");
    time.putAtt("units", "days since 1900-01-01 00:00:00");
    time.putAtt("calendar", "standard");

    M_nc_step=0;

    // Create the vector dimension
    netCDF::NcDim vecDim = dataFile.addDim("x", M_num_drifters);
    netCDF::NcDim emptyDim = dataFile.addDim("y", 1);
    std::vector<netCDF::NcDim> dims2(3);
    dims2[0] = tDim;
    dims2[1] = vecDim;
    dims2[2] = emptyDim;

    // Create the output variables
    // Longitude
    netCDF::NcVar lon = dataFile.addVar("longitude", netCDF::ncFloat, dims2);
    lon.putAtt("standard_name", "longitude");
    lon.putAtt("long_name", "longitude");
    lon.putAtt("units", "degrees_north");

    // Latitude
    netCDF::NcVar lat = dataFile.addVar("latitude", netCDF::ncFloat, dims2);
    lat.putAtt("standard_name", "latitude");
    lat.putAtt("long_name", "latitude");
    lat.putAtt("units", "degrees_north");

    // Reference indices if coming from a grid or a reference dataset
    netCDF::NcVar ind = dataFile.addVar("index", netCDF::ncInt, dims2);
    ind.putAtt("standard_name", "index");
    ind.putAtt("long_name", "reference index");
    ind.putAtt("units","");

    // Concentration - to tell if we have gone out of the ice or not
    netCDF::NcVar sic = dataFile.addVar("sic", netCDF::ncFloat, dims2);
    sic.putAtt("standard_name", "sea_ice_area_fraction");
    sic.putAtt("long_name", "Sea Ice Concentration");
    sic.putAtt("units", "");

    // Global attributes
    dataFile.putAtt("Conventions", "CF-1.6");
    dataFile.putAtt("institution", "NERSC, Thormoehlens gate 47, N-5006 Bergen, Norway");
    dataFile.putAtt("source", "neXtSIM model fields");
}//initNetCDF


// -------------------------------------------------------------------------------------
//! Writes data to the netCDF file.
void
Drifters::appendNetCDF(double current_time)
{

    // Calculate lat and lon
    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    std::vector<char> str(mppfile.begin(), mppfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    std::vector<double> lat(M_num_drifters);
    std::vector<double> lon(M_num_drifters);
    for (int i=0; i<M_num_drifters; ++i)
        inverse_mapx(map, M_X[i], M_Y[i], &lat[i], &lon[i]);

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
    count.push_back(M_num_drifters);

    start.push_back(0);
    count.push_back(1);

    netCDF::NcVar data;

    data = dataFile.getVar("longitude");
    data.putVar(start, count, &lon[0]);

    data = dataFile.getVar("latitude");
    data.putVar(start, count, &lat[0]);
    
    data = dataFile.getVar("index");
    data.putVar(start, count, &M_i[0]);

    data = dataFile.getVar("sic");
    data.putVar(start, count, &M_conc[0]);
}//appendNetCDF

} // Nextsim
