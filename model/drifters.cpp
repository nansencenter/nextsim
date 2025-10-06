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
// ------------------------------------------
//! Main drifter interface to FiniteElement
//! Called from FiniteElement::checkUpdateDrifters()
//! \note call after moving
void
Drifters::updateDrifters(GmshMesh const& movedmesh,
                         std::vector<double> & conc, double const& current_time)
{
    std::vector<double> conc_drifters(0);

    //! 1) Reset temporary (i.e. OSISAF) drifters if needed
    //! \note this does outputting so needs conc_root
    if(this->resetting(current_time))
        this->reset(movedmesh, conc, current_time);

    //! 2) Initialize if needed
    //! - need conc on the moved mesh
    //! - \note updates conc_drifters
    if(this->initialising(current_time))
        this->initialise(movedmesh, conc, conc_drifters);

    //! 3) Add/remove drifters if needed
    //!    \note do this after moving
    if (this->isInputTime(current_time))
    {
        //! - update conc at drifter positions (conc_drifters)
        this->updateConc(movedmesh, conc, conc_drifters);

        //! - add current buoys if not already there
        //!   (output is used for masking later)
        if (M_comm.rank() == 0)
        {
            auto current_buoys = this->grabBuoysFromInputTextFile(current_time);

            //! - Check the drifters map and throw out:
            //!   i) the ones which IABP doesn't report as being in the ice anymore
            //!   (not in current_buoys)
            //!   ii) the ones which have a low conc according to the model
            this->maskXY(conc_drifters, current_buoys);
        }
    }

    //! 4) Add/remove drifters if needed
    if (this->isOutputTime(current_time))
    {
        int num_drifters = conc_drifters.size();
        boost::mpi::broadcast(M_comm, num_drifters, 0);
        if(num_drifters == 0)
            // get conc if needed
            // (haven't added new buoys or initialised this timestep)
            this->updateConc(movedmesh, conc, conc_drifters);
        if (M_comm.rank() == 0) this->outputDrifters(current_time, conc_drifters);
    }
}

// ---------------------------------------------
//! Initialise the drifter positions
//! Called from Drifters::updateDrifters()
void
Drifters::initialise(GmshMesh const& moved_mesh, std::vector<double> & conc,
        std::vector<double> &conc_drifters)
{
    //! -1) Initialise drifter positions
    if (M_comm.rank() == 0 && M_init_type == Drifters::initType::TEXT_FILE)
        this->initFromTextFile();
    else if (M_init_type == Drifters::initType::SPACING)
        this->initFromSpacing(moved_mesh);
    else if (M_comm.rank() == 0 && M_init_type == Drifters::initType::NETCDF)
        this->initFromNetCDF();

    //! -2) Update conc_drifters
    this->updateConc(moved_mesh, conc, conc_drifters);

    //! -3) Applies the mask using conc_drifters and climit, and save to M_X and M_Y
    if (M_comm.rank() == 0) this->maskXY(conc_drifters);

    //! -4) Init the netcdf output file
    if (M_comm.rank() == 0)
    {
        this->setOutputFilename();
        this->initOutputFile();
    }

    //! -5) Success 
    M_is_initialised = true;
}//initialise()


// ---------------------------------------------
//! Initialise the drifter positions on a regular grid
//! Called from Drifters::initialise()
void
Drifters::initFromSpacing(GmshMesh const& moved_mesh)
{
    //! Calculates the grid spacing assuming a regular grid
    std::vector<double> RX = moved_mesh.coordX();
    std::vector<double> RY = moved_mesh.coordY();
    auto xcoords_local = std::minmax_element( RX.begin(), RX.end() );
    auto ycoords_local = std::minmax_element( RY.begin(), RY.end() );

    double xmin = boost::mpi::all_reduce(M_comm, *xcoords_local.first,  boost::mpi::minimum<double>());
    double xmax = boost::mpi::all_reduce(M_comm, *xcoords_local.second, boost::mpi::maximum<double>());
    double ymin = boost::mpi::all_reduce(M_comm, *ycoords_local.first,  boost::mpi::minimum<double>());
    double ymax = boost::mpi::all_reduce(M_comm, *ycoords_local.second, boost::mpi::maximum<double>());

    if (M_comm.rank() != 0) return;

    int ncols = (int) ( 0.5 + ( xmax - xmin )/M_spacing );
    int nrows = (int) ( 0.5 + ( ymax - ymin )/M_spacing );

    M_X.resize(ncols*nrows);
    M_Y.resize(ncols*nrows);
    M_i.resize(ncols*nrows);

    int i=0;
    int y = ymin;
    for (int row=0; row<nrows; ++row)
    {
        int x = xmin;
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
//! Called from Drifters::initialise()
void
Drifters::initFromNetCDF()
{
    std::string const infile = (M_netcdf_input_info.gridFile);
    if ( ! boost::filesystem::exists(M_netcdf_input_info.gridFile) )
        throw std::runtime_error("Drifters::initFromNetCDF: File not found: " + infile);

    // Open file
    netCDF::NcFile dataFile(infile, netCDF::NcFile::read);

    //! - 3) Reads the dimension of the grid
    netCDF::NcDim dim;
    dim = dataFile.getDim(M_netcdf_input_info.dimNameX);
    if ( dim.isNull() )
        throw std::runtime_error("Drifters::initFromNetCDF: Empty dimension: " + M_netcdf_input_info.dimNameX + " in file: " + infile);
    int ncols = dim.getSize();

    dim = dataFile.getDim(M_netcdf_input_info.dimNameY);
    if ( dim.isNull() )
        throw std::runtime_error("Drifters::initFromNetCDF: Empty dimension: " + M_netcdf_input_info.dimNameY + " in file: " + infile);
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
//! Called from Drifters::initialise()
void
Drifters::initFromTextFile()
{
    //! Load the nodes from file

    // Check M_infile, set M_infile_position
    this->initInputTextFile();

    //get the buoys
    M_X.resize(0);
    M_Y.resize(0);
    M_i.resize(0);
    this->grabBuoysFromInputTextFile(M_time_init);
}//initFromTextFile()


// ---------------------------------------------
//!init drifters from fields in restart files
//! Called from FiniteElement::readRestart()
void
Drifters::initFromRestart(
        boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
        boost::unordered_map<std::string, std::vector<double>> & field_map_dbl
        )
{
    double const restart_time = field_map_dbl["Time"][0];
    bool in_restart = false;
    LOG(DEBUG) << M_tag << " drifters: restart time = " << restart_time << " = " << datenumToString(restart_time) << "\n";
    if(M_ignore_restart)
    {
        LOG(DEBUG)<< M_tag<<" drifters: ignoring restart and initialising from scratch\n";
    }
    else
    {
        in_restart = readFromRestart(field_map_int, field_map_dbl);
    }
    LOG(DEBUG) << M_tag << " drifters: init time = " << M_time_init << " = " << datenumToString(M_time_init) << "\n";

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

    //init input file if needed
    if (M_init_type == Drifters::initType::TEXT_FILE)
        this->initInputTextFile();
    
    //init output file if needed
    std::string const backup = M_outfile + ".bak"; 
    bool const present = fs::exists(M_outfile);
    if (present)
        this->backupOutputFile(backup);
    this->initOutputFile();
    if (present)
    {
        this->selectRecordsFromBackup(backup, restart_time);
        std::remove(backup.c_str());
    }
}//initFromRestart()


// ---------------------------------------------
//! Initialise the input text file (M_infile)
//! - check it exists and set M_infile_position
//! Called from Drifters::initFromTextFile and Drifters::initFromRestart()
void
Drifters::initInputTextFile()
{
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
}
 

// --------------------------------------------------------------------------------------
//! Add drifter info to restart
//! Called by Drifters::initFromTextFile()
std::vector<int>
Drifters::grabBuoysFromInputTextFile(double const& current_time)
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
}//grabBuoysFromInputTextFile()


// ---------------------------------------------
//! sort the drifters by ID number - helps testing of restart and output
//! called by Drifters::grabBuoysFromInputTextFile()
void
Drifters::sortDrifterNumbers()
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
}//sortDrifterNumbers()


// --------------------------------------------------------------------------------------
//! Add drifter info to restart
//! Called by FiniteElement::writeRestart()
void
Drifters::addToRestart(Exporter &exporter, MPI_File& outbin, Communicator M_comm, MPI_Offset& base_offset)
{
    // Do nothing if we don't have to
    int keep_going = M_i.size() + M_is_initialised;
    boost::mpi::broadcast(M_comm, &keep_going, 1, 0);
    if (!keep_going) return;

    // write the fields to file
    std::vector<double> const t = {M_time_init};
    std::vector<int> unused;
    exporter.writeField(outbin, M_i, "Drifter_ID_"        + M_tag, "int",  M_comm, base_offset, unused, 0, 0, 1);
    exporter.writeField(outbin, M_X, "Drifter_x_"         + M_tag, "double", M_comm, base_offset, unused, 0, 0, 1);
    exporter.writeField(outbin, M_Y, "Drifter_y_"         + M_tag, "double", M_comm, base_offset, unused, 0, 0, 1);
    exporter.writeField(outbin, t,   "Drifter_time_init_" + M_tag, "double", M_comm, base_offset, unused, 0, 0, 1);
}//addToRestart()


// --------------------------------------------------------------------------------------
//! Read drifter info from restart
//! Called by Drifters::initFromRestart()
bool
Drifters::readFromRestart(
    boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
    boost::unordered_map<std::string, std::vector<double>> & field_map_dbl
    )
{
    std::string const key = "Drifter_ID_" + M_tag;
    if(field_map_int.count(key) == 0)
    {
        LOG(DEBUG) << "Warning: Couldn't read " << M_tag << " drifter positions from restart file."
            << " Drifter positions will be initialised as if there was no restart.\n";
        return false;
    }
    M_i         = field_map_int[key];
    M_X         = field_map_dbl["Drifter_x_"         + M_tag];
    M_Y         = field_map_dbl["Drifter_y_"         + M_tag];
    M_time_init = field_map_dbl["Drifter_time_init_" + M_tag][0];
    return true;
}//readFromRestart()


// --------------------------------------------------------------------------------------
//! Check init time is consistent with restart time
//! - try to fix; otherwise throw error
//! Called by Drifters::initFromRestart()
void
Drifters::fixInitTimeAtRestart(double const& restart_time)
{
    // if we are restarting before scheduled init time, there is no problem
    // - they will be initialised at that time
    if(restart_time<=M_time_init)
        return;

    //otherwise, if we can only start on a particular date, we raise an error...
    if(M_fixed_time_init)
    {
        std::stringstream msg;
        msg << "Error: " << M_tag<< " drifters have been requested but they are not in restart file,\n"
            << "or M_ignore_restart has been set to true,\n"
            << "and their initial time, which is fixed at " << datenumToString(M_time_init) << ",\n"
            << "is before the restart time ("<< datenumToString(restart_time) << ").\n"
            << "Either switch off the drifters or make their initialisation date later.";
        throw std::runtime_error(msg.str());
    }
    
    //...but if init time is not fixed we can just start at the same time on the next day
    double const shift = M_time_init - std::floor(M_time_init);//usually 0, but .5 for OSISAF
    M_time_init = std::ceil(restart_time) + shift;
}//fixInitTimeAtRestart()


// -------------------------------------------------------------------------------------
//! reset drifters
//! - so far only used by OSISAF drifters (reset them after 2 days)
//! Called by FiniteElement::checkUpdateDrifters()
void
Drifters::reset(GmshMesh const& movedmesh, std::vector<double> & conc,
                double const& current_time)
{
    //! 1) Output final positions
    std::vector<double> conc_drifters;
    this->updateConc(movedmesh, conc, conc_drifters);
    if (M_comm.rank() == 0) this->outputDrifters(current_time, conc_drifters);

    //! 2) Set back to being uninitialised
    M_is_initialised = false;
    M_i.resize(0);
    M_X.resize(0);
    M_Y.resize(0);
    M_time_init += M_lifetime;//new init time
    if (M_comm.rank() == 0) this->setOutputFilename();
}//reset()


// --------------------------------------------------------------------------------------
//! Check whether a point is inside a triangle
//! called by Drifters::find_partition()
int
Drifters::inside(std::vector<double> const& points, double xp, double yp)
{
    std::vector<double> cross(3);
    for (int i = 0; i < 3; i++)
    {
        cross[i] = (points[(2+2*i)%6]  - points[2*i]) * (yp - points[2*i+1]) - (points[(2*i+3)%6] - points[2*i+1]) * (xp - points[2*i]);
        if (i > 0 && cross[i] * cross[i-1] < 0.) return 0;
    }

    return 1;
}

// --------------------------------------------------------------------------------------
//! Dispatch the drifters over the partitions and find the triangle containing it
//! called by Drifters::move() and Drifters::update_conc()
void
Drifters::find_partition(GmshMesh const& mesh, std::vector<double>& M_local_drifter_X, std::vector<double>& M_local_drifter_Y, std::vector<int>& M_triangle, std::vector<int>& M_nb_drifter)
{
    int size1 = M_X.size();
    boost::mpi::broadcast(M_comm, size1, 0);

    if (M_comm.rank() != 0)
    {
       M_X.resize(size1);
       M_Y.resize(size1);
    }

    boost::mpi::broadcast(M_comm, &M_X[0], size1, 0);
    boost::mpi::broadcast(M_comm, &M_Y[0], size1, 0);

    int GRID_SIZE = sqrt(mesh.numTriangles()/20);
    double xmax = -1.e20;
    double xmin = 1.e20;
    double ymax = -1.e20;
    double ymin = 1.e20;
    double EPSILON = 1.e-6;

    std::vector<double> coordX = mesh.coordX();
    std::vector<double> coordY = mesh.coordY();
    std::vector<int> triangles = mesh.indexTr();
    for (int k = 0; k <= coordX.size(); k++)
    {
        if (coordX[k] < xmin) xmin = coordX[k];
        if (coordX[k] > xmax) xmax = coordX[k];
        if (coordY[k] < ymin) ymin = coordY[k];
        if (coordY[k] > ymax) ymax = coordY[k];
    }

    double cellWidth = (xmax - xmin) / GRID_SIZE;
    double cellHeight = (ymax - ymin) / GRID_SIZE;

    // Organizing the triangles in a Cartesian grid of the local domain
    std::vector<std::vector<std::vector<int>>> list_triangles(GRID_SIZE, std::vector<std::vector<int>>(GRID_SIZE, std::vector<int>(0,0)));

    std::vector<int> list_idx(6), list_idy(6);
    int minx, maxx, miny, maxy;
    int index;

    for (int k = 0; k < mesh.numTrianglesWithoutGhost(); k++)
    {
        for (int i = 0; i < 3; i++) {
            index = triangles[3*k+i]-1;
            list_idx[i] = (int)((coordX[index] - xmin) / cellWidth + EPSILON);
            list_idy[i] = (int)((coordY[index] - ymin) / cellHeight + EPSILON);
            list_idx[3+i] = (int)((coordX[index] - xmin) / cellWidth - EPSILON);
            list_idy[3+i] = (int)((coordY[index] - ymin) / cellHeight - EPSILON);
        }

        minx = *std::min_element(list_idx.begin(), list_idx.end());
        if (minx < 0) minx = 0;
        maxx = *std::max_element(list_idx.begin(), list_idx.end());
        if (maxx >= GRID_SIZE) maxx = GRID_SIZE-1;
        miny = *std::min_element(list_idy.begin(), list_idy.end());
        if (miny < 0) miny = 0;
        maxy = *std::max_element(list_idy.begin(), list_idy.end());
        if (maxy >= GRID_SIZE) maxy = GRID_SIZE-1;

        for (int idx = minx; idx <= maxx; idx++) {
            for (int idy = miny; idy <= maxy; idy++) {
                list_triangles[idx][idy].push_back(k);
            }
        }
    }

    // Check whether the drifters are in a triangle of the local mesh
    for (int k = 0; k < M_X.size(); k++)
    {
        for (int i = 0; i < 3; i++) {
            list_idx[i] = (int)((M_X[k] - xmin) / cellWidth + EPSILON);
            list_idy[i] = (int)((M_Y[k] - ymin) / cellHeight + EPSILON);
            list_idx[3+i] = (int)((M_X[k] - xmin) / cellWidth - EPSILON);
            list_idy[3+i] = (int)((M_Y[k] - ymin) / cellHeight - EPSILON);
        }

        minx = *std::min_element(list_idx.begin(), list_idx.end());
        if (minx < 0) minx = 0;
        maxx = *std::max_element(list_idx.begin(), list_idx.end());
        if (maxx >= GRID_SIZE) maxx = GRID_SIZE-1;
        miny = *std::min_element(list_idy.begin(), list_idy.end());
        if (miny < 0) miny = 0;
        maxy = *std::max_element(list_idy.begin(), list_idy.end());
        if (maxy >= GRID_SIZE) maxy = GRID_SIZE-1;

        int should_break = 0;
        for (int idx = minx; idx <= maxx; idx++) {
            for (int idy = miny; idy <= maxy; idy++) {
                if (list_triangles[idx][idy].empty()) continue;
                for (int j = 0; j < list_triangles[idx][idy].size(); j++) {
                    int n = list_triangles[idx][idy][j];
                    if (!inside(std::vector<double> {coordX[triangles[3*n]-1], coordY[triangles[3*n]-1],
                                                     coordX[triangles[3*n+1]-1], coordY[triangles[3*n+1]-1],
                                                     coordX[triangles[3*n+2]-1], coordY[triangles[3*n+2]-1]},
                                                     M_X[k], M_Y[k])) continue;
                    M_local_drifter_X.push_back(M_X[k]);
                    M_local_drifter_Y.push_back(M_Y[k]);
                    M_nb_drifter.push_back(k);
                    M_triangle.push_back(n);
                    should_break = 1;
                    break;
                }
                if (should_break) break;
            }
            if (should_break) break;
        }
    }
}

// --------------------------------------------------------------------------------------
//! Move drifters and replace the old coordinates with the new ones
//! called by FiniteElement::updateDrifters()
void
Drifters::move(GmshMesh const& mesh, std::vector<double> const& UT)
{
    boost::mpi::broadcast(M_comm, M_is_initialised, 0);

    // Do nothing if we don't have to
    if ( !M_is_initialised ) return;

    std::vector<int> M_triangle, M_nb_drifter;
    std::vector<double> M_local_drifter_X, M_local_drifter_Y;

    // Send the drifters on the local corresponding partitions
    this->find_partition(mesh, M_local_drifter_X, M_local_drifter_Y, M_triangle, M_nb_drifter);

    // Interpolate the total displacement onto the drifter positions that are within the local domain
    int const numNodes = mesh.numNodes();
    std::vector<int> triangles = mesh.indexTr();
    std::vector<double> coordX = mesh.coordX();
    std::vector<double> coordY = mesh.coordY();
    for ( int i = 0; i < M_local_drifter_X.size(); ++i )
    {
        int k = M_triangle[i];
        double sumnx = 0.;
        double sumny = 0.;
        double sumd = 0.;
        for (int j = 0; j < 3; j++)
        {
            int index = triangles[3*k+j]-1;
            double distance = sqrt((M_local_drifter_X[i] - coordX[index]) * (M_local_drifter_X[i] - coordX[index]) +
                                   (M_local_drifter_Y[i] - coordY[index]) * (M_local_drifter_Y[i] - coordY[index]));

            if (distance < 1.e-15)
            {
                sumnx = UT[index];
                sumny = UT[index + numNodes];
                sumd = 1;
                break;
            }

            sumnx += UT[index] / distance;
            sumny += UT[index + numNodes] / distance;
            sumd  += 1./distance;
        }

        // Add the displacement to the current position 
        M_local_drifter_X[i] += sumnx / sumd;
        M_local_drifter_Y[i] += sumny / sumd;
    }

    // Reconstruction of the drifters positions on root
    std::vector<int> M_drifter_i;
    std::vector<double> M_drifter_X, M_drifter_Y;

    int size = M_nb_drifter.size();

    if (M_comm.rank() == 0)
    {
        M_drifter_i.resize(M_i.size());
        for (int i = 0; i < M_i.size(); i++) M_drifter_i[i] = -1;
        M_drifter_X.resize(M_i.size());
        M_drifter_Y.resize(M_i.size());

        for (int i = 0; i < M_nb_drifter.size(); i++)
        {
            M_drifter_i[i] = M_nb_drifter[i];
            M_drifter_X[i] = M_local_drifter_X[i];
            M_drifter_Y[i] = M_local_drifter_Y[i];
        }

        for (int i = 1; i < M_comm.size(); i++)
        {
            M_comm.recv(i,i,size);
            if (size == 0) continue;

            std::vector<int> local_drifter_i;
            std::vector<double> local_drifter_X, local_drifter_Y;
            M_comm.recv(i,i,local_drifter_i);
            M_comm.recv(i,i,local_drifter_X);
            M_comm.recv(i,i,local_drifter_Y);

            M_drifter_i.insert(M_drifter_i.end(), local_drifter_i.begin(), local_drifter_i.end());
            M_drifter_X.insert(M_drifter_X.end(), local_drifter_X.begin(), local_drifter_X.end());
            M_drifter_Y.insert(M_drifter_Y.end(), local_drifter_Y.begin(), local_drifter_Y.end());
        }

    }
    else
    {
        M_comm.send(0, M_comm.rank(), size);
        if (size != 0) 
        {
            M_comm.send(0, M_comm.rank(), M_nb_drifter);
            M_comm.send(0, M_comm.rank(), M_local_drifter_X);
            M_comm.send(0, M_comm.rank(), M_local_drifter_Y);
        }
    }

    if (M_comm.rank() != 0) {
        M_X.resize(0);
        M_Y.resize(0);
        return;
    }

    for (int i = 0; i < M_drifter_i.size(); i++)
    {
        if (M_drifter_i[i] == -1) continue; // Do not touch outside drifters
        M_X[M_drifter_i[i]] = M_drifter_X[i];
        M_Y[M_drifter_i[i]] = M_drifter_Y[i];
    }

} // move

// --------------------------------------------------------------------------------------
//! interp conc onto drifter positions
//! called by updateDrifters(), reset() and initialise()
void
Drifters::updateConc(GmshMesh const& moved_mesh,
                     std::vector<double> & conc, std::vector<double> &conc_drifters)
{
    // Do nothing if we don't have to
    int num_drifters = M_i.size();
    boost::mpi::broadcast(M_comm, num_drifters, 0);

    if ( num_drifters == 0 ) return;
    conc_drifters.resize(num_drifters);

    // move the mesh before interpolating
    int const numNodes = moved_mesh.numNodes();
    int const numElements = moved_mesh.numTriangles();

    std::vector<int> M_triangle, M_nb_drifter;
    std::vector<double> M_local_drifter_X, M_local_drifter_Y;

    // Send the drifters on the local corresponding partitions
    this->find_partition(moved_mesh, M_local_drifter_X, M_local_drifter_Y, M_triangle, M_nb_drifter);

    std::vector<double> conc_local_drifters(M_local_drifter_X.size());
    for ( int i = 0; i < M_local_drifter_X.size(); ++i ) conc_local_drifters[i] = conc[M_triangle[i]];

    // Reconstruction of the drifters positions on root
    std::vector<int> M_drifter_i;
    std::vector<double> global_conc;

    int size = M_nb_drifter.size();

    if (M_comm.rank() == 0)
    {
        M_drifter_i.resize(M_i.size());
        global_conc.resize(M_i.size());

        std::vector<int> rcounts(M_comm.size());
        boost::mpi::gather(M_comm, size, rcounts, 0);
        boost::mpi::gatherv(M_comm, M_nb_drifter, &M_drifter_i[0], rcounts, 0);
        boost::mpi::gatherv(M_comm, conc_local_drifters, &global_conc[0], rcounts, 0);
    }
    else
    {
        boost::mpi::gather(M_comm, size, 0);
        boost::mpi::gatherv(M_comm, M_nb_drifter, 0);
        boost::mpi::gatherv(M_comm, conc_local_drifters, 0);
    }

    if (M_comm.rank() != 0) return;

    for (int i = 0; i < M_drifter_i.size(); i++)
    {
        if (M_drifter_i[i] == -1)
        {
            conc_drifters[M_drifter_i[i]] = 0.;
            continue;
        }
        conc_drifters[M_drifter_i[i]] = global_conc[i];
    }

}//updateConc

// --------------------------------------------------------------------------------------
//! Masks out X and Y values where there is no ice
//! Also fills M_i with the indices that are kept
void
Drifters::maskXY(std::vector<double> & conc_drifters, std::vector<int> const& keepers)
{

    // Do nothing if we don't have to
    if ( M_i.size() == 0 )
        return;

    auto X = M_X;
    auto Y = M_Y;
    auto idx = M_i;
    auto conc = conc_drifters;

    //! Adds drifter positions where conc > conc_lim
    M_X.resize(0); // this shouldn't be necessary!
    M_Y.resize(0);
    M_i.resize(0);
    conc_drifters.resize(0);

    int identical = 0;
    if (idx == keepers) identical = 1;

    if (identical)
    {
        for ( int i=0; i<idx.size(); ++i )
        {
            if ( conc[i] > M_conc_lim )
            {
                M_X.push_back(X[i]);
                M_Y.push_back(Y[i]);
                M_i.push_back(idx[i]);
                conc_drifters.push_back(conc[i]);
            }
        }
    }
    else
    {
        int max_id = *std::max_element(idx.begin(), idx.end());
        std::vector<int> indices(max_id+1,0);
        for ( int i=0; i<keepers.size(); ++i) indices[keepers[i]] = 1;

        for ( int i=0; i<idx.size(); ++i )
        {
            if ( conc[i] > M_conc_lim && indices[idx[i]] )
            {
                M_X.push_back(X[i]);
                M_Y.push_back(Y[i]);
                M_i.push_back(idx[i]);
                conc_drifters.push_back(conc[i]);
            }
        }
    }
}//maskXY()


// --------------------------------------------------------------------------------------
//! Check the drifter output time step is consistent with model time step
//! Called by FiniteElement::instantiateDrifters()
void
Drifters::checkOutputTimeStep(int time_step)
{
    if( fmod(M_output_interval*24*3600, time_step) != 0 )
    {
        std::string msg = M_tag +" drifters' timestep not a multiple of model time step";
        throw std::runtime_error(msg);
    }
}//checkOutputTimeStep()


// --------------------------------------------------------------------------------------
//! Set timing info: init time, output/input intervals, lifetime,
//! output file name.
//! Also set switches: M_has_lifetime, M_fixed_time_init
//! Called by Drifters::Drifters()
void
Drifters::setTimingInfo(TimingInfo const& timing_info)
{
    //for Drifters
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
}//setTimingInfo()


// --------------------------------------------------------------------------------------
//! Backup outfile at restart time
//! Called by Drifters::initFromRestart(),
void
Drifters::backupOutputFile(std::string const& backup)
{
    fs::path path1(M_outfile);
    if ( fs::exists(path1) )
    {
        fs::path path2(backup);
        fs::copy_file(path1, path2, fs::copy_option::overwrite_if_exists);
    }
}//backupOutputFile()


// ----------------------------------------------------------------------------------------------
//! Initializes the netCDF file
//! * We use netCDF files to save a fixed number of drifters. The file
//! contains two vectors of latitude and longitude positions and can be
//! appended to through time. File name contains time stamp for the creation
//! date (format: yyyymmdd, we assume we won't be creating multiple files
//! with the same prefix per day.
void
Drifters::initOutputNetCDF()
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
//! Called by Drifters::initialise() and Drifters::initFromRestart()
void
Drifters::initOutputTextFile()
{
    std::fstream fout(M_outfile, std::fstream::out);
    if ( !fout.good() )
        throw std::runtime_error("Cannot write to file: " + M_outfile);
    fout << "Year Month Day Hour BuoyID Lat Lon Concentration\n";
    fout.close();
}//initOutputTextFile


// -----------------------------------------------------------------
//! if outfile exists already, we only keep the records for times
//! prior to the restart time
//! this is to prevent duplication of time records inthe netcdf file
//! Called by Drifters::initFromRestart()
void
Drifters::selectRecordsFromBackupNetCDF(
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
        msg << "Drifters::selectRecordsFromBackupNetCDF: Empty dimension x in "
            << backup;
        throw std::runtime_error(msg.str());
    }
    int num_drifters = dim.getSize();

    // time dimension
    dim = dataFile1.getDim("time");
    if ( dim.isNull() )
    {
        M_nc_step = 0;
        LOG(DEBUG) << "Drifters::selectRecordsFromBackupNetCDF: Empty dimension time: nothing to do\n";
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
//! Called by Drifters::selectRecordsFromBackup()
void
Drifters::selectRecordsFromBackupTextFile(
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
Drifters::appendNetCDF(double const& current_time,
        std::vector<double> const& conc_drifters)
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
    data.putVar(start, count, &conc_drifters[0]);
}//appendNetCDF


// -----------------------------------------------------------------------------------------------------------
//! Outputs the IABP drifter positions and conc
//! Called by the readRestart(), initDrifters() and checkDrifters() functions.
void
Drifters::appendTextFile(double const& current_time,
        std::vector<double> const& conc_drifters)
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
             << " " << conc_drifters[j]
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
Drifters::isInputTime(double const& current_time)
{
    // can only input if it's initialised
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


} // Nextsim
