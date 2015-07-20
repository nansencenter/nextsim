/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   options.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Tue Jul 14 13:23:45 2015
 */


#include "boost/program_options.hpp"

namespace po = boost::program_options;

namespace Nextsim
{
    po::options_description
    descrOptions()
    {
        po::options_description desc("Options");

        desc.add_options()
            ("help,h", "Print help messages")
            ("info", "Print help messages")
            ("log_summary", "Print help messages")
            ("config-file", po::value<std::string>(), "specify a configuration file")
            ("config-files", po::value<std::vector<std::string> >()->multitoken(), "specify a list of configuration files" )
            ("nx", po::value<int>()->default_value( 150 ), "Record length in x direction")
            ("ny", po::value<int>()->default_value( 4 ), "Record length in y direction")
            ("dx", po::value<double>()->default_value( 4e+3 ), "Resolution in x direction")
            ("dy", po::value<double>()->default_value( 4e+4 ), "Resolution in y direction")
            ;
        return desc;
    }

} // WIMOPT
