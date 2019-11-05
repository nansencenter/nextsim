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

        enum initType
        {
            TEXT_FILE = 1,
            SPACING   = 2,
            NETCDF    = 3
        };

        typedef struct NetCDFInputInfo
        {
            NetCDFInputInfo() {}
            NetCDFInputInfo(std::string const& infile,
                    std::string const& xname, std::string const& yname,
                    std::string const& lat_name, std::string const& lon_name) :
                gridFile(infile),
                dimNameX(xname), dimNameY(yname),
                latName(lat_name), lonName(lon_name)
            {}
            std::string gridFile;
            std::string dimNameX;
            std::string dimNameY;
            std::string latName;
            std::string lonName;
        } NetCDFInputInfo;

        Drifters() {}

        //! init equally-spaced drifters
        Drifters(std::string const& tag, std::string const& output_prefix,
                double const& spacing, double const& climit,
                DriftersBase::TimingInfo const& timing_info,
                bool const& no_start_from_restart):
            DriftersBase(tag, output_prefix, climit, timing_info,
                    no_start_from_restart)
        {
            M_init_type = Drifters::initType::SPACING;
            M_spacing = spacing;
        }

        //! init drifters from netcdf file
        Drifters(std::string const& tag, std::string const& output_prefix,
                 NetCDFInputInfo const& netcdf_input_info,
                 double const& climit,
                 DriftersBase::TimingInfo const& timing_info,
                 bool const& no_start_from_restart):
            DriftersBase(tag, output_prefix, climit, timing_info,
                    no_start_from_restart)
        {
            M_init_type = Drifters::initType::NETCDF;
            M_netcdf_input_info = netcdf_input_info;
        }

        //! init drifters from text file
        Drifters(std::string const& tag, std::string const& output_prefix,
                std::string const& filename, double const& climit,
                DriftersBase::TimingInfo const& timing_info,
                bool const& no_start_from_restart):
            DriftersBase(tag, output_prefix, climit, timing_info,
                    no_start_from_restart)
        {
            M_init_type = Drifters::initType::TEXT_FILE;
            M_infile = filename;
        }
        

        void initFromRestart(
                boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
                boost::unordered_map<std::string, std::vector<double>> & field_map_dbl);
        void initialise(GmshMeshSeq const& moved_mesh, std::vector<double> & conc);
        void appendNetCDF(double const& current_time);
        void doIO(GmshMeshSeq const& movedmesh_root, std::vector<double> & conc_root,
                double const& current_time);

private:
        void initFromSpacing(GmshMeshSeq const& moved_mesh);
        void initFromTextFile();
        void initFromNetCDF();
        void initNetCDF();
        void backupOutputNetCDF(std::string const& backup);
        void selectRecordsFromBackupNetCDF(
            std::string const& backup, double const& current_time);

        size_t M_nc_step;
        initType M_init_type;
        NetCDFInputInfo M_netcdf_input_info;
        double M_spacing;
    };
} // Nextsim

#endif // __Drifters_H
