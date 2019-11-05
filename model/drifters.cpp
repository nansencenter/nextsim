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
void
Drifters::initFromRestart(
        boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
        boost::unordered_map<std::string, std::vector<double>> & field_map_dbl
        )
{
    double const restart_time = field_map_dbl["Time"][0];
    bool const in_restart = readFromRestart(field_map_int, field_map_dbl);
    if( !in_restart )
        //drifters are not in restart file - check init time and init output file
        this->fixInitTimeAtRestart(restart_time);

    M_outfile = M_output_prefix
        + datenumToString(M_time_init, "%Y%m%d.nc");

    if(!in_restart)
    {
        this->initNetCDF();
        return;
    }
    
    //if drifters are in the restart file
    M_is_initialised = true;
    std::string const backup = M_outfile + ".bak"; 
    bool const present = fs::exists(M_outfile);
    if (present)
        this->backupOutputFile(backup);
    this->initNetCDF();
    if (present)
        this->selectRecordsFromBackupNetCDF(backup, restart_time);
}
 

// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! * defines a uniform grid with resolution defined by spacing
//! Called by FiniteElement::initEquallySpacedDrifters()
Drifters::Drifters(
        std::string const& tag, std::string const& output_prefix,
        double const& spacing, double const& climit,
        TimingInfo const& timing_info)
{
    //! -1) Set the time and output freq
    M_tag = tag;
    M_output_prefix = output_prefix;
    M_conc_lim = climit;
    M_spacing = spacing;
    M_init_type = Drifters::initType::SPACING;
    this->setTimingInfo(timing_info);
}


// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! Called by FiniteElement::initSidfexDrifters() and FiniteElement::initRGPSDrifters()
Drifters::Drifters(
        std::string const& tag, std::string const& output_prefix,
        std::string const& infile, double const& climit,
        TimingInfo const& timing_info)
{
    // interface for RGPS, SIDFEX
    // - reads a text file
    // - NB usually doesn't work for true SIDFEX buoy IDs as they are too large
    //   - therefore use an index file to map the true IDs to smaller integers
    //     (this is what is done in the forecast system)

    //! -1) Set the time and output freq
    M_tag = tag;
    M_output_prefix = output_prefix;
    M_infile = infile;
    M_conc_lim = climit;
    M_init_type = Drifters::initType::TEXT_FILE;
    this->setTimingInfo(timing_info);
}


// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! * reads a netcdf file
//! Called by FiniteElement::initOsisafDrifters()
Drifters::Drifters(std::string const& tag, std::string const& output_prefix,
                   NetCDFInputInfo const& netcdf_input_info,
                   double const& climit,
                   TimingInfo const& timing_info)
{
    M_tag = tag;
    M_output_prefix = output_prefix;
    M_conc_lim = climit;
    M_netcdf_input_info = netcdf_input_info;
    M_init_type = initType::NETCDF;
    this->setTimingInfo(timing_info);
}//init from netcdf file


void
Drifters::initialise(GmshMeshSeq const& moved_mesh, std::vector<double> & conc)
{
    M_is_initialised = true;

    if (M_init_type == Drifters::initType::TEXT_FILE)
        this->initFromTextFile();
    else if (M_init_type == Drifters::initType::SPACING)
        this->initFromSpacing(moved_mesh);
    else if (M_init_type == Drifters::initType::NETCDF)
        this->initFromNetCDF();
    M_num_drifters = M_i.size();

    //! -3) Set M_conc at all the drifters
    this->updateConc(moved_mesh, conc);

    //! -4) Applies the mask using M_conc and climit, and save to M_X and M_Y
    this->maskXY();

    //! -5) Init the netcdf output file
    this->initNetCDF();
}


void
Drifters::initFromSpacing(GmshMeshSeq const& moved_mesh)
{
    //! - 2) Calculates the grid spacing assuming a regular grid
    std::vector<double> RX = moved_mesh.coordX();
    std::vector<double> RY = moved_mesh.coordY();
    auto xcoords = std::minmax_element( RX.begin(), RX.end() );
    auto ycoords = std::minmax_element( RY.begin(), RY.end() );

    int ncols = (int) ( 0.5 + ( *xcoords.second - *xcoords.first )/M_spacing );
    int nrows = (int) ( 0.5 + ( *ycoords.second - *ycoords.first )/M_spacing );

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
            x += M_spacing;
        }
        y += M_spacing;
    }
}


void
Drifters::selectRecordsFromBackupNetCDF(
        std::string const& backup, double const& current_time)
{
    // Open file
    netCDF::NcFile dataFile1(backup, netCDF::NcFile::read);
    netCDF::NcFile dataFile2(backup, netCDF::NcFile::write);

    //! Read x dimension of the grid
    netCDF::NcDim dim = dataFile1.getDim("x");
    if ( dim.isNull() )
    {
        std::stringstream msg;
        msg << "Drifters::selectRecordsFromBackupNetCDF: Empty dimension x in "
            << backup;
        throw std::runtime_error(msg.str());
    }
    M_num_drifters = dim.getSize();

    // time dimension
    dim = dataFile1.getDim("time");
    if ( dim.isNull() )
    {
        M_nc_step = 0;
        std::cout << "Drifters::selectRecordsFromBackupNetCDF: Empty dimension time: nothing to do\n";
        return;
    }
    size_t ntime = dim.getSize();
    std::vector<double> time1(ntime), time2(0);
    netCDF::NcVar vtime1 = dataFile1.getVar("time");
    vtime1.getVar(&time1[0]);
    for (auto t: time1)
        while(t<current_time)
            time2.push_back(t);
    netCDF::NcVar vtime2 = dataFile2.getVar("time");
    vtime2.putVar(&time2[0]);
    M_nc_step = time2.size();

    // Read the latitude and longitude and adds it to the new file
    std::vector<std::string> vars = {
        "latitude", "longitude", "index", "sic"};
    std::vector<double> tmp(M_num_drifters);
    std::vector<size_t> count = {1, (size_t) M_num_drifters, 1};//time, x, y
    for(int n=0; n<M_nc_step; n++)
    {
        //copy 1 time record at a time to save memory
        std::vector<size_t> start = {(size_t) n, 0, 0};//time, x, y
        for (auto vname : vars)
        {
            netCDF::NcVar v1 = dataFile1.getVar(vname);
            netCDF::NcVar v2 = dataFile2.getVar(vname);
            v1.getVar(start, count, &tmp[0]);
            v2.putVar(start, count, &tmp[0]);
        }
    }
}//selectRecordsFromBackupNetCDF


void
Drifters::initFromNetCDF()
{
    std::string const infile = (M_netcdf_input_info.gridFile);
    if ( ! boost::filesystem::exists(M_netcdf_input_info.gridFile) )
        throw std::runtime_error("Drifters::Drifters: File not found: " + infile);

    // Open file
    netCDF::NcFile dataFile(infile, netCDF::NcFile::read);

    //! - 3) Reads the dimension of the grid
    netCDF::NcDim dim;
    dim = dataFile.getDim(M_netcdf_input_info.dimNameX);
    if ( dim.isNull() )
        throw std::runtime_error("Drifters::Drifters: Empty dimension: " + M_netcdf_input_info.dimNameX + " in file: " + infile);
    int ncols = dim.getSize();

    dim = dataFile.getDim(M_netcdf_input_info.dimNameY);
    if ( dim.isNull() )
        throw std::runtime_error("Drifters::Drifters: Empty dimension: " + M_netcdf_input_info.dimNameY + " in file: " + infile);
    int nrows = dim.getSize();

    int gridSize = ncols*nrows;
    M_X.resize(gridSize);
    M_Y.resize(gridSize);
    M_i.resize(gridSize);

    //! - 4) Reads the latitude and longitude
    netCDF::NcVar VLAT = dataFile.getVar(M_netcdf_input_info.latName);
    netCDF::NcVar VLON = dataFile.getVar(M_netcdf_input_info.lonName);

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
}//init from netcdf file


void
Drifters::initFromTextFile()
{
    //! - 2) Load the nodes from file

    // Check file
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
    this->grabBuoysFromInputTextFile(M_time_init);
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
Drifters::initNetCDF()
{
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
Drifters::doIO(GmshMeshSeq const& movedmesh_root, std::vector<double> & conc_root, double const& current_time)
{
    if ( !this->isOutputTime(current_time) )
        return;
    this->updateConc(movedmesh_root, conc_root);
    this->appendNetCDF(current_time);
}

} // Nextsim
