/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   timer.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 13 Mar 2019 09:56:02 CET
 */

#ifndef __Timer_HPP
#define __Timer_HPP 1

#include <boost/mpi/timer.hpp>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>

class Timer
{
    public:
        ///////////////////////////////////////////////////////////////////////
        // Constructors (and destructor)
        ///////////////////////////////////////////////////////////////////////
        
        Timer();

        ~Timer();

        ///////////////////////////////////////////////////////////////////////
        // Public functions
        ///////////////////////////////////////////////////////////////////////
        
        void tick(const std::string & name);
        void tock(const std::string & name);
        const double lap(const std::string & name);
        const double elapsed(const std::string & name);
        const std::string printAll();

        typedef Timer timer;

    private:
        std::string printChildren(const std::string & parent, const double wall_time);
        void printTimer(const std::string & name, const std::string & parent, const double wall_time,
                std::stringstream & return_string, double & not_counted);

        typedef struct Works
        {
            boost::mpi::timer timer;
            double elapsed; // Total elapsed time since first "tick"
            double lap; // Elapsed time since last "tick" or of last "tick"-"tock" interval
            std::string parent;
            int generation;
            bool running;
        } Works;

        std::map<std::string,Works> M_timer; // The central timer object
        std::vector<std::string> M_names;  // A vector of all timer names
        std::vector<std::string> M_lineage; // A vector of timer lineage {... ,grandparent, parent, child}
        std::string M_global_timer = "Total"; // Row label / name for the global timer
        std::string M_not_counted = "Unaccounted for"; // Row label for the "not counted time"

        std::size_t M_width = M_not_counted.size(); // Width of the name column (dynamically adjusted)
        static const int M_padding = 2; // Multiplier for the indent for each generation
};

#endif
