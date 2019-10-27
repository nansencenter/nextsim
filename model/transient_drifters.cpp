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

        //init from vectors (eg from restart)
TransientDrifters::TransientDrifters(std::vector<int> const& buoy_id_in, std::vector<double> const& x_in,
                std::vector<double> const& y_in, std::vector<double> const& conc_in,
                std::string const& infile, std::string const& outfile,
                double const& init_time, double const& output_freq,
                double const& input_freq, double const& conc_lim)
{
                M_i = buoy_id_in;
                M_X = x_in;
                M_Y = y_in;
                M_conc = conc_in;
                M_is_initialised = true;
                M_time_init = init_time;
                M_output_freq = output_freq;
                M_input_freq = input_freq;
                M_conc_lim = conc_lim;
}
 

// ---------------------------------------------------------------------------------------
//! Initializes drifters : seeds and destroys drifters.
//! Called by FiniteElement::initSidfexDrifters() and FiniteElement::initRGPSDrifters()
TransientDrifters::TransientDrifters(std::string const& infile,
        std::string const& outfile,
        GmshMeshSeq const& movedmesh,
        std::vector<double> conc, double const& climit,
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

    //! - 2) Prepare input and output files
    this->initFiles();

    //! - 3) Load the current buoys from file
    M_X.resize(0);
    M_Y.resize(0);
    M_i.resize(0);
    this->grabBuoysFromInputFile(current_time);

    //! - 4) Calculate conc for all the drifters
    this->updateConc(movedmesh, conc);

    //! - 5) Applies mask using conc and climit
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
    std::fstream fin(M_filename, std::fstream::in);
    if ( !fin.good() )
        throw std::runtime_error("File not found: " + M_filename);

    //skip header
    std::string header;
    std::getline(fin, header);
    M_infile_position = fin.tellg();
    fin.close();
    //LOG(DEBUG)<<"open transient drifter file: "<<M_filename<<"\n";
    //LOG(DEBUG)<<"header: "<<header<<"\n";
}//initFiles


void
TransientDrifters::updateDrifters(GmshMeshSeq const& movedmesh_root, std::vector<double>& conc_root,
        double const& current_time)
{
    //add current buoys if not already there
    //(output is used for masking later)
    auto current_buoys = this->grabBuoysFromInputFile(current_time);

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
        do_input = std::fmod(current_time - M_time_init, M_input_freq) == 0;
    return do_input;
}


// --------------------------------------------------------------------------------------
//! Determine if we need to input a drifter
//! Called by outputtingDrifters()
void
TransientDrifters::checkAndDoIO(GmshMeshSeq movedmesh_root, std::vector<double> & conc_root, double const& current_time)
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
