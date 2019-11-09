/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:31:24 CET 2017
 */

#include <drifters_base.hpp>

/**
 * @class DriftersBase
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


// ---------------------------------------------
//! Initialise the drifter positions
//! Called from FiniteElement::updateDrifters()
void
DriftersBase::initialise(GmshMeshSeq const& moved_mesh, std::vector<double> & conc)
{
    M_is_initialised = true;

    if (M_init_type == DriftersBase::initType::TEXT_FILE)
        this->initFromTextFile();
    else if (M_init_type == DriftersBase::initType::SPACING)
        this->initFromSpacing(moved_mesh);
    else if (M_init_type == DriftersBase::initType::NETCDF)
        this->initFromNetCDF();

    //! -3) Set M_conc at all the drifters
    this->updateConc(moved_mesh, conc);

    //! -4) Applies the mask using M_conc and climit, and save to M_X and M_Y
    this->maskXY();

    //! -5) Init the netcdf output file
    this->setOutputFilename();
    this->initOutputFile();
}//initialise()


// ---------------------------------------------
//! Initialise the drifter positions on a regular grid
//! Called from DriftersBase::initialise()
void
DriftersBase::initFromSpacing(GmshMeshSeq const& moved_mesh)
{
    //! - 1) Calculates the grid spacing assuming a regular grid
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
}//initFromSpacing()


// ---------------------------------------------
//! Initialise the drifter positions from the grid of an input netcdf file
//! Called from DriftersBase::initialise()
void
DriftersBase::initFromNetCDF()
{
    std::string const infile = (M_netcdf_input_info.gridFile);
    if ( ! boost::filesystem::exists(M_netcdf_input_info.gridFile) )
        throw std::runtime_error("DriftersBase::DriftersBase: File not found: " + infile);

    // Open file
    netCDF::NcFile dataFile(infile, netCDF::NcFile::read);

    //! - 3) Reads the dimension of the grid
    netCDF::NcDim dim;
    dim = dataFile.getDim(M_netcdf_input_info.dimNameX);
    if ( dim.isNull() )
        throw std::runtime_error("DriftersBase::DriftersBase: Empty dimension: " + M_netcdf_input_info.dimNameX + " in file: " + infile);
    int ncols = dim.getSize();

    dim = dataFile.getDim(M_netcdf_input_info.dimNameY);
    if ( dim.isNull() )
        throw std::runtime_error("DriftersBase::DriftersBase: Empty dimension: " + M_netcdf_input_info.dimNameY + " in file: " + infile);
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
}//initFromNetCDF


// ---------------------------------------------
//! Initialise the drifter positions from a text file like the IABP, RGPS and SIDFEx drifters
//! Called from DriftersBase::initialise()
void
DriftersBase::initFromTextFile()
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
}//initFromTextFile()
 

// --------------------------------------------------------------------------------------
//! Add drifter info to restart
//! Called by DriftersBase::initFromTextFile()
std::vector<int>
DriftersBase::grabBuoysFromInputTextFile(double const& current_time)
{
    std::fstream fin(M_infile, std::fstream::in);
    fin.seekg(M_infile_position);

    int id_count;
    int year, month, day, hour, number;
    double lat, lon, x, y, time;
    std::vector<int> current_buoys(0);

    // NB need to use eof() and old_number instead of commented code since that version stopped
    // after the first line for one example of input file
    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    map = init_mapx( const_cast<char *>(mppfile.c_str()) );
    std::string line;
    while ( std::getline(fin, line) )
    {
        std::istringstream iss(line);
        iss >> year >> month >> day >> hour >> number >> lat >> lon;
        time = Nextsim::getDatenum(year, month, day, hour);
        if(time== current_time )
        {
            // Project and add the buoy to the map if it's missing
            // - don't try to add the conc yet
            id_count = std::count(M_i.begin(), M_i.end(), number);
            if ( id_count == 0 )
            {
                forward_mapx(map, lat, lon, &x, &y);
                M_i.push_back(number);
                M_X.push_back(x);
                M_Y.push_back(y);
                //LOG(DEBUG)<<"new transient buoy: time, index, lon, lat, x, y: "
                //    <<datenumToString(time)<<", "<<number<<", "
                //    <<lon<<", "<<lat<<", "<<x<<", "<<y<<"\n";
            }
            current_buoys.push_back(number);
        }
        else if(time>current_time)
            break;
        M_infile_position = fin.tellg();
    }
    close_mapx(map);
    fin.close();

    //sorting the drifter numbers helps with testing (both restarting and output)
    this->sortDrifterNumbers();
    return current_buoys;
}


//init from vectors (eg from restart)
void
DriftersBase::initFromRestart(
        boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
        boost::unordered_map<std::string, std::vector<double>> & field_map_dbl
        )
{
    double const restart_time = field_map_dbl["Time"][0];
    bool in_restart = false;
    if(!M_ignore_restart)
        in_restart = readFromRestart(field_map_int, field_map_dbl);
    if( !in_restart )
        //drifters are not in restart file - check init time and init output file
        this->fixInitTimeAtRestart(restart_time);
    this->setOutputFilename();

    if(!in_restart)
    {
        this->initOutputFile();
        return;
    }
    
    //if drifters are in the restart file
    M_is_initialised = true;
    std::string const backup = M_outfile + ".bak"; 
    bool const present = fs::exists(M_outfile);
    if (present)
        this->backupOutputFile(backup);
    this->initOutputFile();
    if (present)
        this->selectRecordsFromBackup(backup, restart_time);
}//initFromRestart()


// ---------------------------------------------
//! sort the drifters by ID number - helps testing of restart and output
//! called by DriftersBase::grabBuoysFromInputTextFile()
void
DriftersBase::sortDrifterNumbers()
{
    auto drifter_num = M_i;
    auto drifter_x = M_X;
    auto drifter_y = M_Y;
    std::vector<int> idx(M_i.size());
    int j = 0;
    std::iota(idx.begin(), idx.end(), j++); // {0, 1, 2, ..., num_drifters-1}
    sort(idx.begin(), idx.end(), [&](int i, int j){ return M_i[i]<M_i[j];} );
    for ( int j=0; j<M_i.size(); j++ )
    {
        M_i[j] = drifter_num[idx[j]];
        M_X[j] = drifter_x[idx[j]];
        M_Y[j] = drifter_y[idx[j]];
    }
}


// --------------------------------------------------------------------------------------
//! Add drifter info to restart
//! Called by FiniteElement::writeRestart()
void
DriftersBase::addToRestart(Exporter &exporter, std::fstream &outbin)
{
    // Do nothing if we don't have to
    if (!M_is_initialised)
        return;
    if (M_i.size() == 0)
        return;

    // write the fields to file
    std::vector<double> const t = {M_time_init};
    exporter.writeField(outbin, M_i, "Drifter_ID_"        + M_tag);
    exporter.writeField(outbin, M_X, "Drifter_x_"         + M_tag);
    exporter.writeField(outbin, M_Y, "Drifter_y_"         + M_tag);
    exporter.writeField(outbin, t,   "Drifter_time_init_" + M_tag);
}


// --------------------------------------------------------------------------------------
//! Read drifter info from restart
//! Called by DriftersBase::initFromRestart()
bool
DriftersBase::readFromRestart(
    boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
    boost::unordered_map<std::string, std::vector<double>> & field_map_dbl
    )
{
    std::string const key = "Drifter_ID_" + M_tag;
    if(field_map_int.count(key) == 0)
    {
        std::cout << "Warning: Couldn't read " << M_tag << " drifter positions from restart file."
            << " Drifter positions will be initialised as if there was no restart.\n";
        return false;
    }
    M_i         = field_map_int[key];
    M_X         = field_map_dbl["Drifter_x_"         + M_tag];
    M_Y         = field_map_dbl["Drifter_y_"         + M_tag];
    M_time_init = field_map_dbl["Drifter_time_init_" + M_tag][0];
    return true;
}


// --------------------------------------------------------------------------------------
//! Check init time is consistent with restart time
//! - try to fix; otherwise throw error
//! Called by DriftersBase::initFromRestart()
void
DriftersBase::fixInitTimeAtRestart(double const& restart_time)
{
    // if we are restarting before sceduled init time, there is no problem
    // - they will be initialised at that time
    if(restart_time<=M_time_init)
        return;

    //otherwise, if we can only start on a particular date, we raise an error...
    if(M_fixed_time_init)
    {
        std::stringstream msg;
        msg << "Error: " << M_tag<< " drifters have been requested but they are not in restart file,\n"
            << "and their initial time, which is fixed at " << datenumToString(M_time_init) << ",\n"
            << "is before the restart time ("<< datenumToString(restart_time) << ").\n"
            << "Either switch off the drifters or make their initialisation date later.";
        throw std::runtime_error(msg.str());
    }
    
    //...but if init time is not fixed we can just start at the same time on the next day
    double const shift = M_time_init - std::floor(M_time_init);//usually 0, but .5 for OSISAF
    M_time_init = std::ceil(restart_time) + shift;
}


// -------------------------------------------------------------------------------------
//! reset drifters
//! - so far only used by OSISAF drifters (reset them after 2 days)
//! Called by FiniteElement::updateDrifters()
void
DriftersBase::reset()
{
    M_is_initialised = false;
    M_i.resize(0);
    M_X.resize(0);
    M_Y.resize(0);
    M_conc.resize(0);
    M_time_init += M_lifetime;//new init time
    this->setOutputFilename();
}//reset()


// --------------------------------------------------------------------------------------
//! Move drifters and replace the old coordinates with the new ones
//! called by FiniteElement::updateDrifters()
void
DriftersBase::move(GmshMeshSeq const& mesh,
        std::vector<double> const& UT)
{
    // Do nothing if we don't have to
    if ( !M_is_initialised )
        return;
    int num_drifters = M_i.size();
    if ( num_drifters == 0 )
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
                            &mesh.indexTr()[0], &mesh.coordX()[0], &mesh.coordY()[0],
                            numNodes, mesh.numTriangles(),
                            &interp_drifter_in[0],
                            numNodes, nb_var,
                            &M_X[0], &M_Y[0], num_drifters,
                            true, 0.);

    // Add the displacement to the current position
    for ( int i=0; i<num_drifters; ++i )
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
DriftersBase::updateConc(GmshMeshSeq const& movedmesh,
        std::vector<double> & conc)
{
    // Do nothing if we don't have to
    int num_drifters = M_i.size();
    if ( num_drifters == 0 )
        return;
    M_conc.resize(num_drifters);

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
                            &M_X[0], &M_Y[0], num_drifters,
                            true, 0.);

    for ( int i=0; i<num_drifters; ++i )
        M_conc[i] = interp_drifter_out[i];

    xDelete<double>(interp_drifter_out);
}//updateConc


// --------------------------------------------------------------------------------------
//! Masks out X and Y values where there is no ice
//! Also fills M_i with the indices that are kept
void
DriftersBase::maskXY(std::vector<int> const& keepers)
{

    // Do nothing if we don't have to
    if ( M_i.size() == 0 )
        return;

    auto X = M_X;
    auto Y = M_Y;
    auto idx = M_i;
    auto conc = M_conc;

    //! - 2) Adds drifter positions where conc > conc_lim
    M_X.resize(0); // this shouldn't be necessary!
    M_Y.resize(0);
    M_i.resize(0);
    M_conc.resize(0);

    for ( int i=0; i<idx.size(); ++i )
    {
        int const id_count = std::count(keepers.begin(),
                    keepers.end(), idx[i]);
        if ( conc[i] > M_conc_lim && id_count>0 )
        {
            M_X.push_back(X[i]);
            M_Y.push_back(Y[i]);
            M_i.push_back(idx[i]);
            M_conc.push_back(conc[i]);
        }
    }
}//maskXY()


// --------------------------------------------------------------------------------------
//! Check the drifter output time step is consistent with model time step
//! Called by FiniteElement::instantiateDrifters()
void
DriftersBase::checkOutputTimeStep(int time_step)
{
    if( fmod(M_output_interval*24*3600, time_step) != 0 )
    {
        std::string msg = M_tag +" drifters' timestep not a multiple of model time step";
        throw std::runtime_error(msg);
    }
}


// --------------------------------------------------------------------------------------
//! Set timing info: init time, output/input intervals, lifetime,
//! output file name.
//! Also set switches: M_has_lifetime, M_fixed_time_init
//! Called by DriftersBase::DriftersBase()
void
DriftersBase::setTimingInfo(TimingInfo const& timing_info)
{
    //for DriftersBase
    M_time_init       = timing_info.time_init;
    M_output_interval = timing_info.output_interval;
    M_has_lifetime    = timing_info.has_lifetime;
    M_lifetime        = timing_info.lifetime;
    M_fixed_time_init = timing_info.fixed_time_init;
    this->setOutputFilename();

    if(M_has_lifetime)
    {
        if( M_output_interval > M_lifetime )
        {
            // output timestep should be <= lifetime
            std::stringstream msg;
            msg << M_tag << " drifters output timestep (" << M_output_interval
                << ") should be <= their lifetime (" << M_lifetime << ")";
            throw std::runtime_error(msg.str());
        }
        else if( fmod(M_lifetime, M_output_interval) != 0 )
        {
            // output timestep should fit into lifetime
            std::stringstream msg;
            msg << M_tag << " drifters lifetime (" << M_lifetime
                << ") should be a multiple of their output timestep (" << M_output_interval << ")";
            throw std::runtime_error(msg.str());
        }
    }

    //for transient drifters
    M_transient = timing_info.transient;
    if (!M_transient)
        return;

    M_input_interval = timing_info.input_interval;
    if( M_output_interval > M_input_interval )
    {
        std::string msg = M_tag + " drifters output timestep";
        msg += " should be <= "+M_tag+" input timestep";
        throw std::runtime_error(msg);
    }
    else if ( std::fmod(M_input_interval, M_output_interval) != 0 )
    {
        std::string const msg = M_tag + " drifter input timestep should be a multiple of the "
            + M_tag + " output timestep";
        throw std::runtime_error(msg);
    }
}


// --------------------------------------------------------------------------------------
//! Backup outfile at restart time
//! Called by DriftersBase::initFromRestart(),
void
DriftersBase::backupOutputFile(std::string const& backup)
{
    fs::path path1(M_outfile);
    if ( fs::exists(path1) )
    {
        fs::path path2(backup);
        fs::copy_file(path1, path2, fs::copy_option::overwrite_if_exists);
    }
}


// ----------------------------------------------------------------------------------------------
//! Initializes the netCDF file
//! * We use netCDF files to save a fixed number of drifters. The file
//! contains two vectors of latitude and longitude positions and can be
//! appended to through time. File name contains time stamp for the creation
//! date (format: yyyymmdd, we assume we won't be creating multiple files
//! with the same prefix per day.
void
DriftersBase::initOutputNetCDF()
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
    netCDF::NcDim vecDim = dataFile.addDim("x", M_i.size());
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
}//initOutputNetCDF


// ---------------------------------------------------------------------------------------
//! Initializes input/output text files
//! \note if starting from restart, we backup the output file, and only include the records
//! prior to the restart time in the new output file
//! Called by DriftersBase::initialise() and DriftersBase::initFromRestart()
void
DriftersBase::initOutputTextFile()
{
    std::fstream fout(M_outfile, std::fstream::out);
    if ( !fout.good() )
        throw std::runtime_error("Cannot write to file: " + M_outfile);
    fout << "Year Month Day Hour BuoyID Lat Lon Concentration\n";
    fout.close();
}//initTextFiles


// -----------------------------------------------------------------
//! if outfile exists already, we only keep the records for times
//! prior to the restart time
//! this is to prevent duplication of time records inthe netcdf file
//! Called by DriftersBase::initFromRestart()
void
DriftersBase::selectRecordsFromBackupNetCDF(
        std::string const& backup, double const& current_time)
{
    // Open file
    netCDF::NcFile dataFile1(backup, netCDF::NcFile::read);
    netCDF::NcFile dataFile2(M_outfile, netCDF::NcFile::write);

    //! Read x dimension of the grid
    netCDF::NcDim dim = dataFile1.getDim("x");
    if ( dim.isNull() )
    {
        std::stringstream msg;
        msg << "DriftersBase::selectRecordsFromBackupNetCDF: Empty dimension x in "
            << backup;
        throw std::runtime_error(msg.str());
    }
    int num_drifters = dim.getSize();

    // time dimension
    dim = dataFile1.getDim("time");
    if ( dim.isNull() )
    {
        M_nc_step = 0;
        std::cout << "DriftersBase::selectRecordsFromBackupNetCDF: Empty dimension time: nothing to do\n";
        return;
    }
    size_t ntime = dim.getSize();
    std::vector<double> time1(ntime), time2(0);
    netCDF::NcVar vtime1 = dataFile1.getVar("time");
    vtime1.getVar(&time1[0]);
    for (auto t: time1)
    {
        if (t>=current_time)
            break;
        time2.push_back(t);
    }
    M_nc_step = time2.size();
    netCDF::NcVar vtime2 = dataFile2.getVar("time");
    std::vector<size_t> start = {0};
    std::vector<size_t> count = {(size_t) M_nc_step};
    vtime2.putVar(start, count, &time2[0]);

    // Read the latitude and longitude and adds it to the new file
    std::vector<std::string> vars = {
        "latitude", "longitude", "index", "sic"};
    std::vector<double> tmp(num_drifters);
    count = {1, (size_t) num_drifters, 1};//time, x, y
    for(int n=0; n<M_nc_step; n++)
    {
        //copy 1 time record at a time to save memory
        start = {(size_t) n, 0, 0};//time, x, y
        for (auto vname : vars)
        {
            netCDF::NcVar v1 = dataFile1.getVar(vname);
            netCDF::NcVar v2 = dataFile2.getVar(vname);
            v1.getVar(start, count, &tmp[0]);
            v2.putVar(start, count, &tmp[0]);
        }
    }
}//selectRecordsFromBackupNetCDF


// ---------------------------------------------------------------------------------------
//! Backup output text file (if restarting), and only include the records
//! prior to the restart time in the new output file
//! Called by DriftersBase::selectRecordsFromBackup()
void
DriftersBase::selectRecordsFromBackupTextFile(
        std::string const& backup, double const& current_time)
{
    std::fstream fin(backup, std::fstream::in);
    std::fstream fout(M_outfile, std::fstream::out);
    if ( !fout.good() )
        throw std::runtime_error("Cannot write to file: " + M_outfile);

    int year, month, day, hour, number;
    double lat, lon, x, y, time;
    std::string line;
    int count = -1;
    while ( std::getline(fin, line) )
    {
        count++;
        if (count == 0)
        {
            //1st line is header - copy and continue
            fout << line << "\n";
            continue;
        }

        // check time
        std::istringstream iss(line);
        iss >> year >> month >> day >> hour >> number >> lat >> lon;
        time = Nextsim::getDatenum(year, month, day, hour);
        if(time>=current_time)
            break;
        fout << line << "\n";
    }
    fin.close();
    fout.close();
}//selectRecordsFromBackupTextFile()


// -------------------------------------------------------------------------------------
//! Writes data to the netCDF file.
void
DriftersBase::appendNetCDF(double const& current_time)
{

    // Calculate lat and lon
    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    map = init_mapx( const_cast<char *>(mppfile.c_str()) );

    int num_drifters = M_i.size();
    std::vector<double> lat(num_drifters);
    std::vector<double> lon(num_drifters);
    for (int i=0; i<num_drifters; ++i)
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
    count.push_back(num_drifters);

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


// -----------------------------------------------------------------------------------------------------------
//! Outputs the IABP drifter positions and conc
//! Called by the readRestart(), initDrifters() and checkDrifters() functions.
void
DriftersBase::appendTextFile(double const& current_time)
{
    // Initialize the map
    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    map = init_mapx( const_cast<char *>(mppfile.c_str()) );

    //open output file for appending
    std::fstream fout(M_outfile, std::fstream::out | std::fstream::app );
    if ( ! fout.good() )
        throw std::runtime_error("Cannot write to file: " + M_outfile);

    // Loop over the map and output
    auto ptime = datenumToPosixTime(current_time);
    for ( int j=0; j<M_i.size(); j++)
    {
        double lat, lon;
        inverse_mapx(map, M_X[j], M_Y[j], &lat, &lon);
        fout << std::setw(4) << posixTimeToString(ptime, "%Y")
             << " " << std::setw( 2) << posixTimeToString(ptime, "%m")
             << " " << std::setw( 2) << posixTimeToString(ptime, "%d")
             << " " << std::setw( 2) << posixTimeToString(ptime, "%H")
             << " " << std::setw(16) << M_i[j]
             << std::fixed << std::setprecision(5)
             << " " << std::setw( 8) << lat
             << " " << std::setw(10) << lon
             << " " << M_conc[j]
             << "\n";
    }

    //close output file and map
    fout.close();
    close_mapx(map);
}//appendTextFile


// --------------------------------------------------------------------------------------
//! Determine if we need to input a drifter
//! Called by outputtingDrifters()
bool
DriftersBase::isInputTime(double const& current_time)
{
    // can only output if it's initialised
    if(!M_is_initialised)
        return false;
    if(!M_transient)
        return false;

    bool do_input = false;
    if(current_time>M_time_init)
        // input is already done at init time
        do_input = std::fmod(current_time - M_time_init, M_input_interval) == 0;
    return do_input;
}//isInputTime()


// --------------------------------------------------------------------------------------
//! Determine if we need to input a drifter
//! Called by outputtingDrifters()
void
DriftersBase::doIO(GmshMeshSeq const& movedmesh_root, std::vector<double> & conc_root, double const& current_time)
{
    bool need_conc = true;
    if (this->isInputTime(current_time))
    {
        // check if we need to add new IABP drifters
        // NB do this after moving
        // NB this updates M_conc
        this->addRemoveDrifters(movedmesh_root, conc_root, current_time);
        need_conc = false;
    }

    if (this->isOutputTime(current_time))
    {
        if(need_conc)
            this->updateConc(movedmesh_root, conc_root);
        this->outputDrifters(current_time);
    }
}


// ---------------------------------------------------------------------------------------
//! update the transient drifters
//! - check if we should remove any (not in the input file anymore or too low SIC)
//! - check if we should add any (some in the input file that are not in the drifters)
//! Called by DriftersBase::doIO()
void
DriftersBase::addRemoveDrifters(GmshMeshSeq const& movedmesh_root, std::vector<double>& conc_root,
        double const& current_time)
{
    //add current buoys if not already there
    //(output is used for masking later)
    auto current_buoys = this->grabBuoysFromInputTextFile(current_time);

    // update M_iabp_drifters.conc (get model conc at all drifters)
    this->updateConc(movedmesh_root, conc_root);

    // Check the drifters map and throw out:
    // (i) the ones which IABP doesn't report as being in the ice anymore (not in keepers)
    // (ii) the ones which have a low conc according to the model
    this->maskXY(current_buoys);
}//addRemoveDrifters


void
DriftersBase::updateDrifters(
        GmshMeshSeq const& mesh_root,
        GmshMeshSeq const& movedmesh_root,
        std::vector<double> & conc_root,
        std::vector<double> const& UT_root,
        double const& current_time)
{
    //! 1) Move the drifters (if needed)
    // NB M_UT is relative to the fixed mesh, not the moved mesh
    this->move(mesh_root, UT_root);

    //! 2) Reset any temporary drifters if needed (eg OSISAF)
    if(this->resetting(current_time))
    {
        this->doIO(movedmesh_root, conc_root, current_time);
        this->reset();
    }

    //! 3) Initialize if needed
    //! - need conc on the moved mesh
    if(this->initialising(current_time))
        this->initialise(movedmesh_root, conc_root);

    //! 4) Do output if needed (for transient drifters: also do input if needed)
    this->doIO(movedmesh_root, conc_root, current_time);
}//updateDrifters

} // Nextsim
