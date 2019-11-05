/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:31:24 CET 2017
 */

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


// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! Called by FiniteElement::initSidfexDrifters() and FiniteElement::initRGPSDrifters()
TransientDrifters::TransientDrifters(std::string const& tag, std::string const& output_prefix,
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
    M_infile = infile;
    M_output_prefix = output_prefix;
    M_conc_lim = climit;
    this->setTimingInfo(timing_info);
}


void
TransientDrifters::initialise(GmshMeshSeq const& movedmesh, std::vector<double> & conc)
{
    M_is_initialised = true;

    //! - 2) Prepare input and output files
    this->initTextFiles(true, M_time_init);

    //! - 3) Load the current buoys from file
    M_X.resize(0);
    M_Y.resize(0);
    M_i.resize(0);
    this->grabBuoysFromInputTextFile(M_time_init);
    M_num_drifters = M_i.size();

    //! - 4) Calculate conc for all the drifters
    this->updateConc(movedmesh, conc);

    //! - 5) Applies mask using conc and climit
    this->maskXY();
}


void
TransientDrifters::initFromRestart(
                boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
                boost::unordered_map<std::string, std::vector<double>> & field_map_dbl)
{
    double const restart_time = field_map_dbl["Time"][0];
    if( this->readFromRestart(field_map_int, field_map_dbl) )
    {
        M_is_initialised = true;
        this->initTextFiles(false, restart_time);
    }
    else
    {
        this->fixInitTimeAtRestart(restart_time);
    }
}


// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! Called by FiniteElement::initSidfexDrifters() and FiniteElement::initRGPSDrifters()
void
TransientDrifters::initTextFiles(bool const& overwrite, double const& current_time)
{
    // INPUT:
    //new buoy file has a header
    std::fstream fin(M_infile, std::fstream::in);
    if ( !fin.good() )
        throw std::runtime_error("File not found: " + M_infile);

    //skip header
    std::string header;
    std::getline(fin, header);
    M_infile_position = fin.tellg();
    fin.close();
    //LOG(DEBUG)<<"open transient drifter file: "<<M_infile<<"\n";
    //LOG(DEBUG)<<"header: "<<header<<"\n";

    // OUTPUT:
    M_outfile = M_output_prefix + datenumToString(M_time_init, "%Y%m%d.txt");
    fs::path path1(M_outfile);
    if ( fs::exists(path1) && !overwrite )
        this->backupOutputTextFile(current_time);
        // don't overwrite if starting from restart
        // - just overwrite any times>=current_time
        return;

    std::fstream fout(M_outfile, std::fstream::out);
    if ( !fout.good() )
        throw std::runtime_error("Cannot write to file: " + M_outfile);
    fout << "Year Month Day Hour BuoyID Lat Lon Concentration\n";
    fout.close();
}//initTextFiles


void
TransientDrifters::backupOutputTextFile(double const& current_time)
{
    std::string const backup = M_outfile + ".bak";
    this->backupOutputFile(backup);

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
            fout << line;
            continue;
        }

        // check time
        std::istringstream iss(line);
        iss >> year >> month >> day >> hour >> number >> lat >> lon;
        time = Nextsim::getDatenum(year, month, day, hour);
        if(time < current_time )
            fout << line;
        else if(time>=current_time)
            break;
    }
    fin.close();
    fout.close();
}

void
TransientDrifters::updateDrifters(GmshMeshSeq const& movedmesh_root, std::vector<double>& conc_root,
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
        do_input = std::fmod(current_time - M_time_init, M_input_interval) == 0;
    return do_input;
}


// --------------------------------------------------------------------------------------
//! Determine if we need to input a drifter
//! Called by outputtingDrifters()
void
TransientDrifters::doIO(GmshMeshSeq const& movedmesh_root, std::vector<double> & conc_root, double const& current_time)
{
    bool const inputting = this->isInputTime(current_time);
    bool const outputting = this->isOutputTime(current_time);
    if (inputting)
        // check if we need to add new IABP drifters
        // NB do this after moving
        // NB this updates M_iabp_conc
        this->updateDrifters(movedmesh_root, conc_root, current_time);

    if (outputting)
    {
        // output IABP drifters
        // NB do this after moving
        if(!inputting)
            //still need to update its conc
            this->updateConc(movedmesh_root, conc_root);
        this->outputDrifters(current_time);
    }
}


} // Nextsim
