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
DriftersBase::grabBuoysFromInputFile(double const& current_time)
{
    std::fstream fin(M_filename, std::fstream::in);
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
//  TODO move to Drifters class and give each object a suffix string 
void
DriftersBase::addToRestart(Exporter &exporter, std::fstream &outbin)
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


} // Nextsim
