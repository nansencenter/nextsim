/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:31:24 CET 2017
 */

#include <drifters.hpp>
#include <transient_drifters.hpp>

/**
 * @class TransientDrifters
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
 

TransientDrifters::TransientDrifters()
{
    M_is_initialised = false;
    M_num_drifters   = 0;
    M_X.resize(0);
    M_Y.resize(0);
    M_i.resize(0);
}


// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! Called by FiniteElement::initSidfexDrifters() and FiniteElement::initRGPSDrifters()
TransientDrifters::TransientDrifters(std::string const& infile,
        std::string const& outfile,
        GmshMeshSeq const& movedmesh,
        std::vector<double> & conc, double const& climit,
        double const& current_time, double const& output_freq,
        double const& input_freq)
{
    // interface for RGPS, SIDFEX
    // - reads a text file
    // - NB usually doesn't work for true SIDFEX buoy IDs as they are too large
    //   - therefore use an index file to map the true IDs to smaller integers
    //     (this is what is done in the forecast system)

    //! -1) Set the time and output freq
    M_time_init = current_time;
    M_output_freq = output_freq;
    M_input_freq = input_freq;
    M_filename = infile;
    M_outfile = outfile;
    M_is_initialised = true;
    M_conc_lim = climit;

    //! - 2) Load the current buoys from file
    this->initFiles();

#if 0
    //! -3) Read the current buoys from file
    double time = current_time;
    std::vector<int> keepers;

    int year, month, day, hour;
    int number;
    double lat, lon, x, y;

    // Read the next line
    // NB need to use eof() and old_number instead of commented code since that version stopped
    // after the first line for one example of input file
    int gridSize = 0;
    int old_number = 0;

    //while ( M_rgps_file >> year >> month >> day >> hour >> number >> lat >> lon )
    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    map = init_mapx( const_cast<char *>(mppfile.c_str()) );
    while (!M_fstream.eof())
    {

        if (gridSize>0)
            old_number = number;

        M_fstream >> year >> month >> day >> hour >> number >> lat >> lon;
        time = Nextsim::getDatenum(year, month, day, hour);

        // may not be eof if \n at end of file
        // - this can lead to repetition of the last line
        if (gridSize>0)
            if (number == old_number)
                break;

        if(time== current_time )
        {
            forward_mapx(map, lat, lon, &x, &y);
            M_X.push_back(x);
            M_Y.push_back(y);
            M_i.push_back(number);
            gridSize++;

            //std::cout << year << ", "<< month << ", "<< day << ", "
            //    << hour << ", "  << number << ", "<< lat << ", "<< lon << "\n";
        }
        else if(time>current_time)
            break;
    }
    close_mapx(map);
#endif

    //! - 3) Calculate conc for all the drifters
    this->updateConc(movedmesh, conc);

    //! - 4) Applies mask using conc and climit
    this->maskXY();
}


// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! Called by FiniteElement::initSidfexDrifters() and FiniteElement::initRGPSDrifters()
void
TransientDrifters::initFiles()
{
    // OUTPUT:
    // if M_outfile exists, we just append to it when it is time;
    // else we just add the header and close
    fs::path path1(M_outfile);
    if ( !fs::exists(path1) )
    {
        std::fstream fout(M_outfile, std::fstream::out);
        if ( !fout.good() )
            throw std::runtime_error("Cannot write to file: " + M_outfile);
        fout << "Year Month Day Hour BuoyID Lat Lon Concentration\n";
        fout.close();
    }

    // INPUT:
    //new buoy file has a header
    M_fstream.open(M_filename, std::fstream::in);
    if ( !M_fstream.good() )
        throw std::runtime_error("File not found: " + M_filename);

    //skip header
    std::string header;
    std::getline(M_fstream, header);
    LOG(DEBUG)<<"open transient drifter file: "<<M_filename<<"\n";
    LOG(DEBUG)<<"header: "<<header<<"\n";

    int pos;    // To be able to rewind one line
    double time = getDatenum(1979, 1, 1);
    int year, month, day, hour, number;
    double lat, lon, x, y;
    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    map = init_mapx( const_cast<char *>(mppfile.c_str()) );
    M_X.resize(0);
    M_Y.resize(0);
    M_i.resize(0);
    while ( time <= M_time_init )
    {
        // Remember where we were
        pos = M_fstream.tellg();

        // Read the next line
        M_fstream >> year >> month >> day >> hour >> number >> lat >> lon;
        time = getDatenum(year, month, day, hour);

        if(time=M_time_init)
        {
            //loads the current buoys
            forward_mapx(map, lat, lon, &x, &y);
            M_X.push_back(x);
            M_Y.push_back(y);
            M_i.push_back(number);
        }
    }
    close_mapx(map);

    // We must rewind one line so that updateIabpDrifters works correctly
    M_fstream.seekg(pos);
}//initFiles


void
TransientDrifters::updateDrifters(mesh_type_root const& movedmesh_root, std::vector<double>& conc_root,
        double const& current_time)
{

    // Initialize the map
    mapx_class *map;
    std::string mppfile = Environment::nextsimMppfile();
    map = init_mapx( const_cast<char *>(mppfile.c_str()) );

    // Read the current buoys from file
    double time = current_time;
    std::vector<int> keepers;
    int id_count;
    double x, y;
    while ( time == current_time )
    {
        // Read the next line
        int year, month, day, hour, number;
        double lat, lon;
        M_fstream >> year >> month >> day >> hour >> number >> lat >> lon;
        time = getDatenum(year, month, day, hour);

        // Remember which buoys are in the ice according to IABP
        keepers.push_back(number);

        // Project and add the buoy to the map if it's missing
        // - don't try to add the conc yet
        id_count = std::count(M_i.begin(), M_i.end(), number);
        if ( id_count == 0 )
        {
            forward_mapx(map, lat, lon, &x, &y);
            M_i.push_back(number);
            M_X.push_back(x);
            M_Y.push_back(y);
            LOG(DEBUG)<<"new transient buoy: time, index, lon, lat, x, y: "
                <<datenumToString(time)<<", "<<number<<", "
                <<lon<<", "<<lat<<", "<<x<<", "<<y<<"\n";
        }
    }
    close_mapx(map);

    // update M_iabp_drifters.conc (get model conc at all drifters)
    this->updateConc(movedmesh_root, conc_root);

    // Check the drifters map and throw out:
    // (i) the ones which IABP doesn't report as being in the ice anymore (not in keepers)
    // (ii) the ones which have a low conc according to the model
    this->maskXY(keepers);
}//updateDrifters


// -----------------------------------------------------------------------------------------------------------
//! Outputs the IABP drifter positions and conc
//! Called by the readRestart(), initDrifters() and checkDrifters() functions.
void
TransientDrifters::outputDrifters(double const& current_time)
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
    for ( int j=0; j<M_num_drifters; j++)
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
}//outputDrifters


// --------------------------------------------------------------------------------------
//! Add drifter info to restart
//! Called by FiniteElement::writeRestart()
//  TODO move to Drifters class and give each object a suffix string 
void
TransientDrifters::addToRestart(Exporter &exporter, std::fstream &outbin)
{
    // Sort the drifters so the restart files are identical
    // (this is just to make testing of restart procedure easier)
    std::vector<int> drifter_no = M_i;
    std::vector<double> drifter_x(M_num_drifters), drifter_y(M_num_drifters), drifter_conc(M_num_drifters);
    std::vector<int> idx(M_num_drifters);
    int j = 0;
    std::iota(idx.begin(), idx.end(), j++); // {0, 1, 2, ..., M_num_drifters-1}
    sort(idx.begin(), idx.end(), [&](int i, int j){ return drifter_no[i]<drifter_no[j];} );
    for ( int j=0; j<M_num_drifters; j++ )
    {
        drifter_no[j] = M_i[idx[j]];
        drifter_x[j] = M_X[idx[j]];
        drifter_y[j] = M_Y[idx[j]];
        drifter_conc[j] = M_conc[idx[j]];
    }

    // write the fields to file
    std::vector<double> t = {M_time_init};
    exporter.writeField(outbin, drifter_no, "Drifter_no");
    exporter.writeField(outbin, drifter_x, "Drifter_x");
    exporter.writeField(outbin, drifter_y, "Drifter_y");
    exporter.writeField(outbin, drifter_conc, "Drifter_conc");
    exporter.writeField(outbin, t, "Drifter_time_init");
}


// --------------------------------------------------------------------------------------
//! Determine if we need to input a drifter
//! Called by outputtingDrifters()
bool
TransientDrifters::isInputTime(double const& current_time)
{
    // can only output if it's initialised
    if(!M_is_initialised)
        return false;

    bool do_input = false;
    if(current_time>M_time_init)
        // output is already done at init time
        do_input = std::fmod(current_time - M_time_init, M_input_freq) == 0;
    return do_input;
}


} // Nextsim
