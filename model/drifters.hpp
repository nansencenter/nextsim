/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:49:20 CET 2017
 */

#ifndef __Drifters_H
#define __Drifters_H 1

#include <gmshmeshseq.hpp>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromMeshToGridx.h>
#include <netcdf>
#include <date.hpp>

/**
 * @class GridOutput
 * @brief Manage output of gridded data
 *
 * @see
 *
 */

namespace Nextsim
{
    class Drifters
    {
public:

        typedef struct Grid
        {
            std::string gridFile;
            std::string dirname;
            std::string mpp_file;
            std::string dimNameX;
            std::string dimNameY;
            std::string latName;
            std::string lonName;

            bool loaded;
            std::vector<double> gridLAT;
            std::vector<double> gridLON;
            std::vector<double> gridX;
            std::vector<double> gridY;
        } Grid;

        Drifters() {}

        //init from vectors (eg from restart)
        Drifters(std::vector<int> const& buoy_id_in, std::vector<double> const& x_in,
                std::vector<double> const& y_in, std::vector<double> const& conc_in,
                double const& init_time, double const& output_freq,
                double const& conc_lim) :
                M_i(buoy_id_in), M_X(x_in), M_Y(y_in),
                M_conc(conc_in), M_is_initialised(true),
                M_time_init(init_time), M_output_freq(output_freq),
                M_conc_lim(conc_lim)
        {}

        //! init equally-spaced drifters
        Drifters(double const& spacing, GmshMeshSeq const& movedmesh,
                std::vector<double> & conc, double const& climit,
                double const& current_time, double const& output_freq);

        //! init drifters from netcdf file
        Drifters(std::string const& gridFile,
                 std::string const& dimNameX, std::string const& dimNameY,
                 std::string const& latName, std::string const& lonName,
                 GmshMeshSeq const& movedmesh,
                 std::vector<double> & conc, double const& climit,
                 double const& current_time, double const& output_freq);

        //! init drifters from text file
        Drifters(std::string const& filename,
                GmshMeshSeq const& movedmesh,
                std::vector<double> & conc, double const& climit,
                double const& current_time, double const& output_freq);

        void move(GmshMeshSeq const& mesh, std::vector<double> const& UT);
        void updateConc( GmshMeshSeq const& movedmesh,
                std::vector<double> & conc);

        void initNetCDF(std::string file_prefix, double current_time);

        void appendNetCDF(double current_time);

        bool isInitialised();
        bool isOutputTime(double const& current_time);
        void checkAndDoOutput(GmshMeshSeq const& movedmesh_root, std::vector<double> & conc_root,
                double const& current_time);


protected://only derived classes have access to these
        bool M_is_initialised = false;
        int M_num_drifters = 0;
        double M_time_init;
        double M_output_freq;
        double M_conc_lim;


        size_t M_nc_step;

        std::string M_filename;

        std::vector<double> M_X;
        std::vector<double> M_Y;
        std::vector<int> M_i;
        std::vector<double> M_conc;

        void maskXY(std::vector<int> const& keepers); //check if buoy IDs are in a given list and remove if they are not
        void maskXY() { this->maskXY(M_i); } //don't remove any buoys unless the conc is too low
    };
} // Nextsim

#endif // __Drifters_H
