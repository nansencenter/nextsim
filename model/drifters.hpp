/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:49:20 CET 2017
 */

#ifndef __Drifters_H
#define __Drifters_H 1

#include <drifters_base.hpp>
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
    class Drifters: public DriftersBase
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
        Drifters(std::string const& tag, std::string const& outfile_prefix,
                boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
                boost::unordered_map<std::string, std::vector<double>> & field_map_dbl,
                double const& output_freq, double const& conc_lim);

        //! init equally-spaced drifters
        Drifters(std::string const& tag, std::string const& outfile_prefix,
                double const& spacing, GmshMeshSeq const& movedmesh,
                std::vector<double> & conc, double const& climit,
                double const& current_time, double const& output_freq);

        //! init drifters from netcdf file
        Drifters(std::string const& tag, std::string const& outfile_prefix,
                 std::string const& gridFile,
                 std::string const& dimNameX, std::string const& dimNameY,
                 std::string const& latName, std::string const& lonName,
                 GmshMeshSeq const& movedmesh,
                 std::vector<double> & conc, double const& climit,
                 double const& current_time, double const& output_freq);

        //! init drifters from text file
        Drifters(std::string const& tag, std::string const& outfile_prefix,
                std::string const& filename,
                GmshMeshSeq const& movedmesh,
                std::vector<double> & conc, double const& climit,
                double const& current_time, double const& output_freq);

        void appendNetCDF(double const& current_time);
        void checkAndDoOutput(GmshMeshSeq const& movedmesh_root, std::vector<double> & conc_root,
                double const& current_time);

private:
        void initNetCDF(std::string const& outfile_prefix, bool const& overwrite);
        size_t M_nc_step;
    };
} // Nextsim

#endif // __Drifters_H
