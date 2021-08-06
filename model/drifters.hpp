/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:49:20 CET 2017
 */

#ifndef __Drifters_H
#define __Drifters_H 1

#include <environment.hpp>
#include <gmshmeshseq.hpp>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromMeshToGridx.h>
#include <netcdf>
#include <date.hpp>
#include <exporter.hpp>
#include <environment.hpp>
#include "debug.hpp"
#include <numeric>

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

        //how are we initialising?
        enum initType
        {
            TEXT_FILE = 1,
            SPACING   = 2,
            NETCDF    = 3
        };

        //info for netcdf file needed when
        //initialising from restart
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

        // info for drifters about "timing"
        // - input/output intervals, lifetime, init time...
        typedef struct TimingInfo
        {
            TimingInfo() {}

            //constructor for "persistent" drifters
            // - drifters are only added at initialisation
            // (unlike eg IABP)
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

            //constructor for "transient" drifters
            // - drifters can be added and removed periodically
            //  (eg IABP)
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

            double time_init;//init time
            double output_interval;//time between outputs
            double input_interval = -1;//time between inputs (transient drifters)
            bool has_lifetime;//do they have a limited lifetime, like the OSISAF drifters
            double lifetime;//lifetime of drifters in days - eg 2. for OSISAF drifters, after which they are reset
            bool fixed_time_init;// some drifters have a fixed init time (eg depend on a specific input file);
            bool transient = false;// can drifters be added/removed after init?
        } TimingInfo;

        //! Constructors and destructors
        //! * We have one constructor for regularly spaced drifters and two for
        //!   drifter positions read in from file
        //!   (one uses a text file and one uses a netcdf file).
        //! * There is also an empty constructor which has M_is_initialised=false (default value).
        Drifters() {}

        //! Construct a Drifters object with equally-spaced drifters
        Drifters(std::string const& tag, std::string const& output_prefix,
                double const& spacing, double const& climit,
                Drifters::TimingInfo const& timing_info,
                bool const& ignore_restart):
            Drifters(tag, output_prefix, climit, timing_info,
                    ignore_restart)
        {
            M_init_type = Drifters::initType::SPACING;
            M_spacing = spacing;
        }

        //! Construct a Drifters object with drifter positions read from a netcdf file
        Drifters(std::string const& tag, std::string const& output_prefix,
                 NetCDFInputInfo const& netcdf_input_info,
                 double const& climit,
                 Drifters::TimingInfo const& timing_info,
                 bool const& ignore_restart):
            Drifters(tag, output_prefix, climit, timing_info,
                    ignore_restart)
        {
            M_init_type = Drifters::initType::NETCDF;
            M_netcdf_input_info = netcdf_input_info;
        }

        //! Construct a Drifters object with drifter positions from a text file
        Drifters(std::string const& tag, std::string const& output_prefix,
                std::string const& filename, double const& climit,
                Drifters::TimingInfo const& timing_info,
                bool const& ignore_restart):
            Drifters(tag, output_prefix, climit, timing_info,
                    ignore_restart)
        {
            M_init_type = Drifters::initType::TEXT_FILE;
            M_infile = filename;
        }
        bool initialising(double const &current_time)
        {
            if (M_is_initialised)
                return false;
            return current_time == M_time_init;
        }
        void initFromRestart(
                boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
                boost::unordered_map<std::string, std::vector<double>> & field_map_dbl
                );

        void addToRestart(Exporter &exporter, std::fstream &outbin);
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

        //main interface to FiniteElement
        void updateDrifters(GmshMeshSeq const& movedmesh_root,
                std::vector<double> & conc_root, double const& current_time);
        void move(GmshMeshSeq const& mesh, std::vector<double> const& UT);

private:
        //initialising
        Drifters(std::string const& tag, std::string const& output_prefix,
                double const& climit,
                Drifters::TimingInfo const& timing_info,
                bool const& ignore_restart) : 
            M_tag(tag),
            M_output_prefix(output_prefix),
            M_conc_lim(climit),
            M_ignore_restart(ignore_restart),
            M_log_level(Environment::logLevel()),
            M_log_all(Environment::logAll()),
            M_comm(Environment::comm())
        {
            this->setTimingInfo(timing_info);
        }

        void initialise(GmshMeshSeq const& moved_mesh, std::vector<double> & conc,
                std::vector<double> & conc_drifters);
        void initFromSpacing(GmshMeshSeq const& moved_mesh);
        void initFromTextFile();
        void initFromNetCDF();
        bool readFromRestart(
                boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
                boost::unordered_map<std::string, std::vector<double>> & field_map_dbl);
        void fixInitTimeAtRestart(double const& restart_time);
        std::vector<int> grabBuoysFromInputTextFile(double const& current_time);
        void sortDrifterNumbers();

        //main ops
        void reset(GmshMeshSeq const& moved_mesh, std::vector<double> & conc_root,
                double const& current_time);
        void updateConc( GmshMeshSeq const& moved_mesh,
                std::vector<double> & conc, std::vector<double> & conc_drifters);
        bool resetting(double const& current_time)
        {
            if(!M_is_initialised)
                return false;
            if(!M_has_lifetime)
                return false;
            return (current_time == M_time_init + M_lifetime);
        }
        void maskXY(std::vector<double> & conc_drifters,
                std::vector<int> const& current_buoys);
            // - check if buoy IDs are in a given list and remove if they are not
            // - also check for if conc > a threshold
        void maskXY(std::vector<double> & conc_drifters)
        {
            auto keepers = M_i;
            // - don't remove any buoys unless the conc is too low
            // - take a copy as we resize M_i inside maskXY
            this->maskXY(conc_drifters, keepers);
        }
        void setTimingInfo(TimingInfo const& timing_info);

        //IO
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
        void initInputTextFile();
        void initOutputTextFile();
        void initOutputNetCDF();
        bool isInputTime(double const& current_time);
        void outputDrifters(double const& current_time, std::vector<double> const& conc_drifters)
        {
            if(M_transient)
                this->appendTextFile(current_time, conc_drifters);
            else
                this->appendNetCDF(current_time, conc_drifters);
        }
        void appendTextFile(double const& current_time, std::vector<double> const& conc_drifters);
        void appendNetCDF(double const& current_time, std::vector<double> const& conc_drifters);

        // used at restart time - append to output files,
        // but overwrite records after restart time
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

        //variables
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

        //for LOG
        LogLevel M_log_level;
        bool M_log_all;
        Communicator M_comm;
    };
} // Nextsim

#endif /* __Drifters_H */
