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
    return current_buoys;
}


// --------------------------------------------------------------------------------------
//! Add drifter info to restart
//! Called by FiniteElement::writeRestart()
//  TODO give each object a suffix string so we can put all drifters in the restart file
void
DriftersBase::addToRestart(Exporter &exporter, std::fstream &outbin)
{
    // Do nothing if we don't have to
    if (!M_is_initialised)
        return;
    if (M_num_drifters == 0)
        return;

    // Sort the drifters so the restart files are identical
    // (this is just to make testing of restart procedure easier)
    std::vector<int> drifter_no = M_i;
    std::vector<double> drifter_x(M_num_drifters), drifter_y(M_num_drifters);
    std::vector<int> idx(M_num_drifters);
    int j = 0;
    std::iota(idx.begin(), idx.end(), j++); // {0, 1, 2, ..., M_num_drifters-1}
    sort(idx.begin(), idx.end(), [&](int i, int j){ return drifter_no[i]<drifter_no[j];} );
    for ( int j=0; j<M_num_drifters; j++ )
    {
        drifter_no[j] = M_i[idx[j]];
        drifter_x[j] = M_X[idx[j]];
        drifter_y[j] = M_Y[idx[j]];
    }

    // write the fields to file
    std::vector<double> const t = {M_time_init};
    exporter.writeField(outbin, drifter_no,   "Drifter_ID_"        + M_tag);
    exporter.writeField(outbin, drifter_x,    "Drifter_x_"         + M_tag);
    exporter.writeField(outbin, drifter_y,    "Drifter_y_"         + M_tag);
    exporter.writeField(outbin, t,            "Drifter_time_init_" + M_tag);
}


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
    M_num_drifters = M_i.size();
    return true;
}


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
            << "Either switch off the drifters or initialise at a later date.";
        throw std::runtime_error(msg.str());
    }
    
    //...but if init time is not fixed we can just start at the same time on the next day
    double const shift = M_time_init - std::floor(M_time_init);//usually 0, but .5 for OSISAF
    M_time_init = std::ceil(restart_time) + shift;
}


void
DriftersBase::reset()
{
    M_is_initialised = false;
    M_num_drifters = 0;
    M_i.resize(0);
    M_X.resize(0);
    M_Y.resize(0);
    M_conc.resize(0);
    M_time_init += M_lifetime;//new init time
    std::string const ext = M_outfile.substr(
            M_outfile.find_last_of("."));
    M_outfile = M_output_prefix + datenumToString(M_time_init, "%Y%m%d")
        + ext;
}

// --------------------------------------------------------------------------------------
//! Move drifters and replace the old coordinates with the new ones
//! called by FiniteElement::checkDrifters()
void
DriftersBase::move(GmshMeshSeq const& mesh,
        std::vector<double> const& UT)
{
    // Do nothing if we don't have to
    if ( !M_is_initialised )
        return;
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
DriftersBase::updateConc(GmshMeshSeq const& movedmesh,
        std::vector<double> & conc)
{
    // Do nothing if we don't have to
    if ( M_num_drifters == 0 )
        return;
    M_conc.resize(M_num_drifters);

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


// --------------------------------------------------------------------------------------
//! Masks out X and Y values where there is no ice
//! Also fills M_i with the indices that are kept
void
DriftersBase::maskXY(std::vector<int> const& keepers)
{

    // Do nothing if we don't have to
    if ( M_num_drifters == 0 )
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

    M_num_drifters = M_X.size();
}//maskXY()


void
DriftersBase::checkOutputTimeStep(int time_step)
{
    if( fmod(M_output_interval*24*3600, time_step) != 0 )
    {
        std::string msg = M_tag +" drifters' timestep not a multiple of model time step";
        throw std::runtime_error(msg);
    }
}


void
DriftersBase::setTimingInfo(TimingInfo const& timing_info)
{
    //for Drifters & TransientDrifters
    M_time_init       = timing_info.time_init;
    M_output_interval = timing_info.output_interval;
    M_has_lifetime    = timing_info.has_lifetime;
    M_lifetime        = timing_info.lifetime;
    M_fixed_time_init = timing_info.fixed_time_init;
    M_outfile = M_output_prefix
        + datenumToString(M_time_init, "%Y%m%d.nc");
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

    //for TransientDrifters
    if (timing_info.input_interval <0 )
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

} // Nextsim
