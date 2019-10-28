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
    std::vector<double> const t = {M_time_init};
    exporter.writeField(outbin, drifter_no,   "Drifter_no_"        + M_tag);
    exporter.writeField(outbin, drifter_x,    "Drifter_x_"         + M_tag);
    exporter.writeField(outbin, drifter_y,    "Drifter_y_"         + M_tag);
    exporter.writeField(outbin, drifter_conc, "Drifter_conc_"      + M_tag);
    exporter.writeField(outbin, t,            "Drifter_time_init_" + M_tag);
}


void
DriftersBase::readFromRestart(
    boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
    boost::unordered_map<std::string, std::vector<double>> & field_map_dbl
    )
{
    M_i         = field_map_int["Drifter_no_"        + M_tag];
    M_X         = field_map_dbl["Drifter_x_"         + M_tag];
    M_Y         = field_map_dbl["Drifter_y_"         + M_tag];
    M_conc      = field_map_dbl["Drifter_conc_"      + M_tag];
    M_time_init = field_map_dbl["Drifter_time_init_" + M_tag][0];
    M_num_drifters = M_i.size();
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

    auto X = M_X;
    auto Y = M_Y;
    auto INDS = M_i;
    auto conc = M_conc;

    //! - 2) Adds drifter positions where conc > conc_lim
    M_X.resize(0); // this shouldn't be necessary!
    M_Y.resize(0);
    M_i.resize(0);
    M_conc.resize(0);
    int id_count;

    for ( int i=0; i<INDS.size(); ++i )
    {
        id_count = std::count(keepers.begin(),
                    keepers.end(), INDS[i]);
        if ( conc[i] > M_conc_lim
                && id_count>0 )
        {
            M_X.push_back(X[i]);
            M_Y.push_back(Y[i]);
            M_i.push_back(INDS[i]);
            M_conc.push_back(conc[i]);
        }
    }

    M_num_drifters = M_X.size();
}//maskXY()


} // Nextsim
