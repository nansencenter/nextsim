/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:49:20 CET 2017
 */

#ifndef __TransientDrifters_H
#define __TransientDrifters_H 1

#include <gmshmeshseq.hpp>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromMeshToGridx.h>
#include <date.hpp>
#include <exporter.hpp>

/**
 * @class GridOutput
 * @brief Manage output of gridded data
 *
 * @see
 *
 */

namespace Nextsim
{
    class TransientDrifters: public Drifters
    {
public:

        typedef GmshMeshSeq mesh_type_root;

        TransientDrifters();

        //init from vectors (eg from restart)
        TransientDrifters(std::vector<int> const& buoy_id_in, std::vector<double> const& x_in,
                std::vector<double> const& y_in, std::vector<double> const& conc_in,
                std::string const& infile, std::string const& outfile,
                double const& init_time, double const& output_freq,
                double const& input_freq, double const& conc_lim) :
                M_i(buoy_id_in), M_X(x_in), M_Y(y_in),
                M_conc(conc_in), M_is_initialised(true),
                M_time_init(init_time), M_output_freq(output_freq),
                M_input_freq(input_freq), M_conc_lim(conc_lim)
        {}

        //! init drifters from text file
        TransientDrifters(std::string const& infile,
                std::string const& outfile,
                GmshMeshSeq const& movedmesh,
                std::vector<double> & conc, double const& climit,
                double const& current_time, double const& output_freq,
                double const& input_freq);

        void initFiles();
        void updateDrifters(mesh_type_root const& movedmesh_root, std::vector<double>& conc_root,
            double const& current_time);
        void outputDrifters(double const& current_time);
        void addToRestart(Exporter &exporter, std::fstream &outbin);
        bool isInputTime(double const& current_time);

private:
        bool M_is_initialised;
        double M_time_init;
        double M_output_freq;
        double M_input_freq;
        double M_conc_lim;
        int M_num_drifters;

        std::string M_filename; // The file we read the IABP buoy data from
        std::fstream M_fstream;
        std::string M_outfile; // The (text) file we output to

        std::vector<double> M_X;
        std::vector<double> M_Y;
        std::vector<int> M_i;
        std::vector<double> M_conc;
    };
} // Nextsim

#endif /* __TransientDrifters_H */
