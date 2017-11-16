/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:49:20 CET 2017
 */

#ifndef __Drifters_H
#define __Drifters_H 1

#include <gmshmesh.hpp>
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
            Drifters(double spacing, GmshMesh const &mesh, std::vector<double> &conc, double climit);
            Drifters(std::string dirname, std::string gridFile, std::string dimNameX, std::string dimNameY, std::string latName, std::string lonName, GmshMesh const &mesh, std::vector<double> &conc, double climit);
            Drifters(std::string filename, GmshMesh const &mesh, std::vector<double> &conc, double climit, double time);

            void move(GmshMesh const &mesh, std::vector<double> const &UT);

            void initNetCDF(std::string file_prefix, double current_time);
            void appendNetCDF(double current_time, GmshMesh const &mesh, std::vector<double> const &UT);

            bool isInitialised();

        private:
            bool M_is_initialised;

            int M_no_drifters;

            size_t M_nc_step;

            std::string M_filename;

            std::vector<double> M_X;
            std::vector<double> M_Y;
            std::vector<int> M_i;

            void move(GmshMesh const &mesh, std::vector<double> const &UT, std::vector<double> &X, std::vector<double> &Y);

            void maskXY(GmshMesh const &mesh, std::vector<double> &X, std::vector<double> &Y, std::vector<double> &conc, double clim);
    };
} // Nextsim

#endif // __Drifters_H
