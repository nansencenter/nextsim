/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   clock.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 13 Mar 2019 09:56:02 CET
 */

#ifndef __Clock_HPP
#define __Clock_HPP 1

#include <boost/mpi/timer.hpp>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>

class Clock
{
    public:
        ///////////////////////////////////////////////////////////////////////
        // Constructors (and destructor)
        ///////////////////////////////////////////////////////////////////////
        
        Clock();

        Clock(const std::vector<std::string> & names);

        ~Clock();

        ///////////////////////////////////////////////////////////////////////
        // Public functions
        ///////////////////////////////////////////////////////////////////////
        
        void tick(const std::string & name);
        const double tock(const std::string & name);
        const double elapsed(const std::string & name);
        const std::string printAll();

    private:
        typedef struct Works
        {
            boost::mpi::timer timer;
            double elapsed;
            std::string parent;
        } Works;

        std::map<std::string,Works> M_clock;
        std::vector<std::string> M_names;
        std::string M_global_clock = "main";
};

#endif
