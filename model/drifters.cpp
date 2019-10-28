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


//init from vectors (eg from restart)
Drifters::Drifters(
        std::string const& tag, std::string const& outfile_prefix,
        boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
        boost::unordered_map<std::string, std::vector<double>> & field_map_dbl,
        double const& output_freq, double const& conc_lim)
{
    M_tag = tag;
    M_is_initialised = true;
    M_output_freq = output_freq;
    M_conc_lim = conc_lim;
    this -> readFromRestart(field_map_int, field_map_dbl);
    this->initNetCDF(outfile_prefix, false);
}
 

// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! * defines a uniform grid with resolution defined by spacing
//! Called by FiniteElement::initEquallySpacedDrifters()
Drifters::Drifters(
        std::string const& tag, std::string const& outfile_prefix,
        double const& spacing, GmshMeshSeq const& moved_mesh,
        std::vector<double> & conc, double const& climit,
        double const& current_time, double const& output_freq)
{
    //! -1) Set the time and output freq
    M_tag = tag;
    M_time_init = current_time;
    M_output_freq = output_freq;
    M_conc_lim = climit;
    M_is_initialised = true;

    //! - 2) Calculates the grid spacing assuming a regular grid
    std::vector<double> RX = moved_mesh.coordX();
    std::vector<double> RY = moved_mesh.coordY();
    auto xcoords = std::minmax_element( RX.begin(), RX.end() );
    auto ycoords = std::minmax_element( RY.begin(), RY.end() );

    int ncols = (int) ( 0.5 + ( *xcoords.second - *xcoords.first )/spacing );
    int nrows = (int) ( 0.5 + ( *ycoords.second - *ycoords.first )/spacing );

    M_X.resize(ncols*nrows);
    M_Y.resize(ncols*nrows);
    M_i.resize(ncols*nrows);

    int i=0;
    int y = *ycoords.first;
    for (int row=0; row<nrows; ++row)
    {
        int x = *xcoords.first;
        for (int col=0; col<ncols; ++col)
        {
            M_X[i] = x;
            M_Y[i] = y;
            M_i[i] = i;
            ++i;
            x += spacing;
        }
        y += spacing;
    }

    //! -3) Set M_conc at all the drifters
    this->updateConc(moved_mesh, conc);

    //! -4) Applies the mask using M_conc and climit, and save to M_X and M_Y
    this->maskXY();

    //! -5) Init the netcdf output file
    this->initNetCDF(outfile_prefix, true);
}


// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! Called by FiniteElement::initSidfexDrifters() and FiniteElement::initRGPSDrifters()
Drifters::Drifters(
        std::string const& tag, std::string const& outfile_prefix,
        std::string const& infile,
        GmshMeshSeq const& moved_mesh,
        std::vector<double> & conc, double const& climit,
        double const& current_time, double const& output_freq)
{
    // interface for RGPS, SIDFEX
    // - reads a text file
    // - NB usually doesn't work for true SIDFEX buoy IDs as they are too large
    //   - therefore use an index file to map the true IDs to smaller integers
    //     (this is what is done in the forecast system)

    //! -1) Set the time and output freq
    M_tag = tag;
    M_time_init = current_time;
    M_output_freq = output_freq;
    M_i.resize(0);
    M_X.resize(0);
    M_Y.resize(0);
    M_conc_lim = climit;
    M_is_initialised = true;

    //! - 2) Load the nodes from file

    // Check file
    M_infile = infile;
    std::fstream drifter_text_file;   // The file we read the buoy data from
    drifter_text_file.open(M_infile, std::fstream::in);
    if ( !drifter_text_file.good() )
        throw std::runtime_error("File not found: " + M_infile);

    //skip header, save position, and close
    std::string header;
    std::getline(drifter_text_file, header);
    M_infile_position = drifter_text_file.tellg();
    drifter_text_file.close();

    //get the buoys
    M_X.resize(0);
    M_Y.resize(0);
    M_i.resize(0);
    this->grabBuoysFromInputFile(current_time);

    //! -3) Set M_conc at all the drifters
    this->updateConc(moved_mesh, conc);

    //! -4) Applies the mask using M_conc and climit, and save to M_X and M_Y
    this->maskXY();

    //! -5) Init the netcdf output file
    this->initNetCDF(outfile_prefix, true);
}


// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! * reads a netcdf file
//! Called by FiniteElement::initOsisafDrifters()
Drifters::Drifters(std::string const& tag, std::string const& outfile_prefix,
                   std::string const& gridFile,
                   std::string const& dimNameX, std::string const& dimNameY,
                   std::string const& latName, std::string const& lonName,
                   GmshMeshSeq const& moved_mesh,
                   std::vector<double> & conc, double const& climit,
                   double const& current_time, double const& output_freq)
{

    //! -1) Set the time and output freq
    M_tag = tag;
    M_time_init = current_time;
    M_output_freq = output_freq;
    M_conc_lim = climit;
    M_is_initialised = true;

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
    M_X.resize(gridSize);
    M_Y.resize(gridSize);
    M_i.resize(gridSize);

    //! - 4) Reads the latitude and longitude
    netCDF::NcVar VLAT = dataFile.getVar(latName);
    netCDF::NcVar VLON = dataFile.getVar(lonName);

    //! - 5) Reads the data in LON & LAT
    std::vector<double> LAT(gridSize);
    std::vector<double> LON(gridSize);
    VLAT.getVar(&LAT[0]);
    VLON.getVar(&LON[0]);

    //! - 6) Calculates cartesian (x,y) coordinates from the longitude and latitude
    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    map = init_mapx( const_cast<char *>(mppfile.c_str()) );

    for (int i=0; i<gridSize; ++i)
    {
        forward_mapx(map, LAT[i], LON[i], &M_X[i], &M_Y[i]);
        M_i[i] = i;
    }
    close_mapx(map);

    //! - 7) Set M_conc at all the drifters
    this->updateConc(moved_mesh, conc);

    //! - 8) Applies mask using conc and climit, and save to M_X and M_Y
    this->maskXY();

    //! -9) Init the netcdf output file
    this->initNetCDF(outfile_prefix, true);
}//init from text file


// File operations
// ----------------------------------------------------------------------------------------------
//! Initializes the netCDF file
//! * We use netCDF files to save a fixed number of drifters. The file
//! contains two vectors of latitude and longitude positions and can be
//! appended to through time. File name contains time stamp for the creation
//! date (format: yyyymmdd, we assume we won't be creating multiple files
//! with the same prefix per day.
 
void
Drifters::initNetCDF(std::string const& outfile_prefix, bool const& overwrite)
{
    M_outfile = outfile_prefix + datenumToString(M_time_init, "%Y%m%d.nc");
    if ( boost::filesystem::exists(M_outfile) && !overwrite )
        //don't want to overwrite if starting from restart
        return;

    // Create the netCDF file.
    netCDF::NcFile dataFile(M_outfile, netCDF::NcFile::replace);

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
Drifters::appendNetCDF(double const& current_time)
{

    // Calculate lat and lon
    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    map = init_mapx( const_cast<char *>(mppfile.c_str()) );

    std::vector<double> lat(M_num_drifters);
    std::vector<double> lon(M_num_drifters);
    for (int i=0; i<M_num_drifters; ++i)
        inverse_mapx(map, M_X[i], M_Y[i], &lat[i], &lon[i]);

    close_mapx(map);

    // Open the netCDF file
    netCDF::NcFile dataFile(M_outfile, netCDF::NcFile::write);

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


// --------------------------------------------------------------------------------------
//! Determine if we need to input a drifter
//! Called by outputtingDrifters()
void
Drifters::checkAndDoOutput(GmshMeshSeq const& movedmesh_root, std::vector<double> & conc_root, double const& current_time)
{
    if ( this->isOutputTime(current_time) )
    {
        this->updateConc(movedmesh_root, conc_root);
        this->appendNetCDF(current_time);
    }
}

} // Nextsim
