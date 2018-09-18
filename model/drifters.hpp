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

        Drifters();

        //! init equally-spaced drifters
        Drifters(double spacing, GmshMeshSeq const& movedmesh,
                std::vector<double>& conc, double climit);

        //! init drifters from netcdf file
        Drifters(std::string gridFile,
                 std::string dimNameX, std::string dimNameY,
                 std::string latName, std::string lonName,
                 GmshMeshSeq const& movedmesh,
                 std::vector<double>& conc,
                 double climit);

        //! init drifters from text file
        Drifters(std::string filename,
                GmshMeshSeq const& movedmesh,
                std::vector<double>& conc, double climit, double time);

        void move(GmshMeshSeq const& mesh, std::vector<double> const& UT);
        void updateConc( GmshMeshSeq const& movedmesh,
                std::vector<double> & conc);

        void initNetCDF(std::string file_prefix, double current_time);

        void appendNetCDF(double current_time);

        bool isInitialised();
        bool isOutputTime(double const& current_time);

        // --------------------------------------------------------------------------------------
        //! Set the time to initialise the drifters
        //! Called by initDrifterOpts()
        //  TODO can now make a loop over all the drifters
        void setInitTimeAndOutputFreq(double const& current_time, double const& output_freq)
        {
            M_time_init = current_time;
            M_output_freq = output_freq;
        }


private:
        bool M_is_initialised;
        double M_time_init;
        double M_output_freq;

        int M_num_drifters;

        size_t M_nc_step;

        std::string M_filename;

        std::vector<double> M_X;
        std::vector<double> M_Y;
        std::vector<long int> M_i;
        std::vector<double> M_conc;

        void maskXY(GmshMeshSeq const& mesh,
                std::vector<double>& X, std::vector<double>& Y,
                std::vector<long int>& INDS,
                std::vector<double>& conc, double clim);
    };
} // Nextsim

#endif // __Drifters_H
