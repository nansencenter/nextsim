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

        //init from vectors (eg from restart)
        TransientDrifters(std::vector<int> const& buoy_id_in, std::vector<double> const& x_in,
                std::vector<double> const& y_in, std::vector<double> const& conc_in,
                std::string const& infile, std::string const& outfile,
                double const& init_time, double const& output_freq,
                double const& input_freq, double const& conc_lim);

        //! init drifters from text file
        TransientDrifters(std::string const& infile,
                std::string const& outfile,
                GmshMeshSeq const& movedmesh,
                std::vector<double> conc, double const& climit,
                double const& current_time, double const& output_freq,
                double const& input_freq);

        void initFiles();
        void updateDrifters(GmshMeshSeq const& movedmesh_root, std::vector<double>& conc_root,
            double const& current_time);
        void outputDrifters(double const& current_time);
        bool isInputTime(double const& current_time);
        void checkAndDoIO(GmshMeshSeq movedmesh_root, std::vector<double> & conc_root,
                double const& current_time);
    };
} // Nextsim

#endif /* __TransientDrifters_H */
