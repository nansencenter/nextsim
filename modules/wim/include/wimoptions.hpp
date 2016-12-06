/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

/**
 * @file   options.hpp
 * @author Abdoulaye Samake <abdama@beijing.ad.nersc.no>
 * @date   Fri Aug  5 16:07:13 2016
 */

#ifndef __WIMOPTIONS_H
#define __WIMOPTIONS_H 1

#include "boost/program_options.hpp"

namespace Wim
{
namespace po = boost::program_options;

    po::options_description
    descrWimOptions();

} // WIMOPT

#endif
