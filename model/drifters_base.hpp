/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:49:20 CET 2017
 */

#ifndef __DriftersBase_H
#define __DriftersBase_H 1

#include <gmshmeshseq.hpp>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromMeshToGridx.h>
#include <netcdf>
#include <date.hpp>
#include <exporter.hpp>
#include <environment.hpp>
#include "debug.hpp"


/**
 * @class GridOutput
 * @brief Manage output of gridded data
 *
 * @see
 *
 */

namespace Nextsim
{
    class DriftersBase
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

        typedef struct TimingInfo
        {
            TimingInfo() {}

            TimingInfo(
                    double const& ti,
                    double const& oint,
                    bool const& hlt,
                    double const& lt,
                    bool const& fti
                    )
            {
                time_init = ti;
                output_interval = oint;
                has_lifetime = hlt;
                lifetime = lt;
                fixed_time_init = fti;
            }

            TimingInfo(
                    double const& ti,
                    double const& oint,
                    double const& iint,
                    bool const& hlt,
                    double const& lt,
                    bool const& fti
                    ) : TimingInfo(ti, oint, hlt, lt, fti)
            {
                input_interval = iint;
                transient = true;
            }

            double time_init;
            double output_interval;
            double input_interval = -1;
            bool has_lifetime;
            double lifetime;
            bool fixed_time_init;
            bool transient = false;
        } TimingInfo;

        DriftersBase() {}
        DriftersBase(std::string const& tag, std::string const& output_prefix,
                double const& climit,
                DriftersBase::TimingInfo const& timing_info,
                bool const& ignore_restart) : 
            M_tag(tag), M_output_prefix(output_prefix),
            M_conc_lim(climit), M_ignore_restart(ignore_restart)
        { this->setTimingInfo(timing_info); }

        //! init equally-spaced drifters
        DriftersBase(std::string const& tag, std::string const& output_prefix,
                double const& spacing, double const& climit,
                DriftersBase::TimingInfo const& timing_info,
                bool const& no_start_from_restart):
            DriftersBase(tag, output_prefix, climit, timing_info,
                    no_start_from_restart)
        {
            M_init_type = DriftersBase::initType::SPACING;
            M_spacing = spacing;
        }

        //! init drifters from netcdf file
        DriftersBase(std::string const& tag, std::string const& output_prefix,
                 NetCDFInputInfo const& netcdf_input_info,
                 double const& climit,
                 DriftersBase::TimingInfo const& timing_info,
                 bool const& no_start_from_restart):
            DriftersBase(tag, output_prefix, climit, timing_info,
                    no_start_from_restart)
        {
            M_init_type = DriftersBase::initType::NETCDF;
            M_netcdf_input_info = netcdf_input_info;
        }

        //! init drifters from text file
        DriftersBase(std::string const& tag, std::string const& output_prefix,
                std::string const& filename, double const& climit,
                DriftersBase::TimingInfo const& timing_info,
                bool const& no_start_from_restart):
            DriftersBase(tag, output_prefix, climit, timing_info,
                    no_start_from_restart)
        {
            M_init_type = DriftersBase::initType::TEXT_FILE;
            M_infile = filename;
        }
        bool initialising(double const &current_time)
        {
            if (M_is_initialised)
                return false;
            return current_time == M_time_init;
        }
        void initialise(GmshMeshSeq const& moved_mesh, std::vector<double> & conc);
        void initFromRestart(
                boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
                boost::unordered_map<std::string, std::vector<double>> & field_map_dbl
                );

        void addToRestart(Exporter &exporter, std::fstream &outbin);
        void move(GmshMeshSeq const& mesh, std::vector<double> const& UT);
        void reset();
        void updateConc( GmshMeshSeq const& movedmesh,
                std::vector<double> & conc);
        bool isOutputTime(double const& current_time)
        {
            if(!M_is_initialised)
                return false;
            else
                return std::fmod(current_time - M_time_init, M_output_interval) == 0;
        }
        void checkOutputTimeStep(int time_step);
        double getInitTime() { return M_time_init; }
        void setInitTime(double const& t) { M_time_init = t; }
        bool isInitialised() { return M_is_initialised; }
        bool resetting(double const& current_time)
        {
            if(!M_is_initialised)
                return false;
            if(!M_has_lifetime)
                return false;
            return (current_time == M_time_init + M_lifetime);
        }
        void updateDrifters(GmshMeshSeq const& mesh_root,
                GmshMeshSeq const& movedmesh_root,
                std::vector<double> & conc_root,
                std::vector<double> const& UT_root,
                double const& current_time);

private:
        void initFromSpacing(GmshMeshSeq const& moved_mesh);
        void initFromTextFile();
        void initFromNetCDF();
        bool readFromRestart(
                boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
                boost::unordered_map<std::string, std::vector<double>> & field_map_dbl);
        void fixInitTimeAtRestart(double const& restart_time);
        std::vector<int> grabBuoysFromInputTextFile(double const& current_time);
        void sortDrifterNumbers();
        void maskXY(std::vector<int> const& current_buoys);
            // - check if buoy IDs are in a given list and remove if they are not
            // - also check for if conc > a threshold
        void maskXY()
        {
            auto keepers = M_i;
            // - don't remove any buoys unless the conc is too low
            // - take a copy as we resize M_i inside maskXY
            this->maskXY(keepers);
        }
        void setTimingInfo(TimingInfo const& timing_info);

        void setOutputFilename()
        {
            M_outfile = M_output_prefix + datenumToString(M_time_init, "%Y%m%d");
            if(M_transient)
                M_outfile += ".txt";
            else
                M_outfile += ".nc";
        }
        void initOutputFile()
        {
            if(M_transient)
                this->initOutputTextFile();
            else
                this->initOutputNetCDF();
        }
        void initOutputTextFile();
        void initOutputNetCDF();

        void backupOutputFile(std::string const& backup);
        void selectRecordsFromBackup(
                std::string const& backup, double const& current_time)
        {
            if(M_transient)
                this->selectRecordsFromBackupTextFile(
                        backup, current_time);
            else
                this->selectRecordsFromBackupNetCDF(
                        backup, current_time);
        }
        void selectRecordsFromBackupTextFile(
                std::string const& backup, double const& current_time);
        void selectRecordsFromBackupNetCDF(
                std::string const& backup, double const& current_time);

        bool isInputTime(double const& current_time);
        void doIO(GmshMeshSeq const& movedmesh_root,
                std::vector<double> & conc_root, double const& current_time);
        void addRemoveDrifters(GmshMeshSeq const& movedmesh_root,
                std::vector<double>& conc_root, double const& current_time);
        void outputDrifters(double const& current_time)
        {
            if(M_transient)
                this->appendTextFile(current_time);
            else
                this->appendNetCDF(current_time);
        }
        void appendTextFile(double const& current_time);
        void appendNetCDF(double const& current_time);

        bool M_is_initialised = false;
        double M_time_init;//init time
        double M_output_interval;//time between outputs
        double M_input_interval;//time between inputs (transient drifters)
        double M_conc_lim;//threshold concentration for removing buoys
        std::string M_output_prefix;//prefix of output file - init time and extension are added at instantiation/restart time
        double M_lifetime;//lifetime of drifters in days - eg 2. for OSISAF drifters, after which they are reset
        bool M_has_lifetime = false;//do they have a limited lifetime, like the OSISAF drifters
        bool M_fixed_time_init = false;// some drifters have a fixed init time (eg depend on a specific input file);
                                       // true: error is raised if we try to restart after the fixed init time
                                       // false: init time is flexible
        bool M_ignore_restart = false;//don't use restart file, even if drifter is present (eg for SIDFEx in forecast scenario)
        bool M_transient = true;

        size_t M_nc_step;
        initType M_init_type;
        NetCDFInputInfo M_netcdf_input_info;
        double M_spacing;
        int M_infile_position;//position in input text file (if used)
        
        std::string M_infile;//input filename (if used)
        std::string M_outfile;//output filename
        std::string M_tag;//drifter tag, mainly for name in restart files and name of output files

        std::vector<double> M_X;//x coordinates of drifters
        std::vector<double> M_Y;//y coordinates of drifters
        std::vector<int> M_i;//buoy ID of drifters
        std::vector<double> M_conc;//conc at drifter positions TODO needed?
    };
} // Nextsim

#endif /* __DriftersBase_H */
