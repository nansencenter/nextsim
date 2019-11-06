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
            { input_interval = iint; }

            double time_init;
            double output_interval;
            double input_interval = -1;
            bool has_lifetime;
            double lifetime;
            bool fixed_time_init;
        } TimingInfo;

        DriftersBase() {}
        DriftersBase(std::string const& tag, std::string const& output_prefix,
                double const& climit,
                DriftersBase::TimingInfo const& timing_info,
                bool const& ignore_restart) : 
            M_tag(tag), M_output_prefix(output_prefix),
            M_conc_lim(climit), M_ignore_restart(ignore_restart)
        { this->setTimingInfo(timing_info); }

        void addToRestart(Exporter &exporter, std::fstream &outbin);
        void move(GmshMeshSeq const& mesh, std::vector<double> const& UT);
        void reset();
        void updateConc( GmshMeshSeq const& movedmesh,
                std::vector<double> & conc);
        bool initialising(double const &current_time)
        {
            if (M_is_initialised)
                return false;
            return current_time == M_time_init;
        }
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

protected:
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
        void backupOutputFile(std::string const& backup);

        bool M_is_initialised = false;
        double M_time_init;//init time
        double M_output_interval;//time between outputs
        double M_input_interval;//time between inputs (TransientDrifters)
        double M_conc_lim;//threshold concentration for removing buoys
        std::string M_output_prefix;//prefix of output file - init time and extension are added at instantiation/restart time
        double M_lifetime;//lifetime of drifters in days - eg 2. for OSISAF drifters, after which they are reset
        bool M_has_lifetime = false;//do they have a limited lifetime, like the OSISAF drifters
        bool M_fixed_time_init = false;// some drifters have a fixed init time (eg depend on a specific input file);
                                       // true: error is raised if we try to restart after the fixed init time
                                       // false: init time is flexible
        bool M_ignore_restart = false;//
        
        std::string M_infile;
        std::string M_outfile;
        std::string M_tag;
        int M_infile_position;

        std::vector<double> M_X;
        std::vector<double> M_Y;
        std::vector<int> M_i;
        std::vector<double> M_conc;
    };
} // Nextsim

#endif /* __TransientDrifters_H */
