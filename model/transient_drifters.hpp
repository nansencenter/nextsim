/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:49:20 CET 2017
 */

#ifndef __TransientDrifters_H
#define __TransientDrifters_H 1

#include <drifters_base.hpp>
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
    class TransientDrifters: public DriftersBase
    {
public:
        TransientDrifters() {}

        //! init drifters from text file
        TransientDrifters(std::string const& tag, std::string const& output_prefix,
                std::string const& infile, double const& climit,
                DriftersBase::TimingInfo const& timing_info,
                bool const& no_start_from_restart):
            DriftersBase(tag, output_prefix, climit, timing_info,
                    no_start_from_restart)
            { M_infile = infile; }

        void initFromRestart(
                boost::unordered_map<std::string, std::vector<int>>    & field_map_int,
                boost::unordered_map<std::string, std::vector<double>> & field_map_dbl);
        void initialise(GmshMeshSeq const& movedmesh, std::vector<double> & conc);
        void initTextFiles(bool const& overwrite, double const& current_time);
        void backupOutputTextFile(double const& current_time);
        void updateDrifters(GmshMeshSeq const& movedmesh_root, std::vector<double>& conc_root,
            double const& current_time);
        void outputDrifters(double const& current_time);
        bool isInputTime(double const& current_time);
        void doIO(GmshMeshSeq const& movedmesh_root, std::vector<double> & conc_root,
                double const& current_time);
    };
} // Nextsim

#endif /* __TransientDrifters_H */
