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

        DriftersBase() {}
        void addToRestart(Exporter &exporter, std::fstream &outbin);
        bool isOutputTime(double const& current_time);
        void move(GmshMeshSeq const& mesh, std::vector<double> const& UT);
        void updateConc( GmshMeshSeq const& movedmesh,
                std::vector<double> & conc);
        bool isInitialised() { return M_is_initialised; }

protected:
        std::vector<int> grabBuoysFromInputFile(double const& current_time);
        void maskXY(std::vector<int> const& current_buoys); //check if buoy IDs are in a given list and remove if they are not
        void maskXY() { this->maskXY(M_i); } //don't remove any buoys unless the conc is too low

        bool M_is_initialised = false;
        int M_num_drifters = 0;
        double M_time_init;
        double M_output_freq;
        double M_input_freq;
        double M_conc_lim;

        std::string M_infile;
        std::string M_outfile;
        int M_infile_position;

        std::vector<double> M_X;
        std::vector<double> M_Y;
        std::vector<int> M_i;
        std::vector<double> M_conc;
    };
} // Nextsim

#endif /* __TransientDrifters_H */
