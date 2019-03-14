/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   clock.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 13 Mar 2019 09:56:02 CET
 */

#include "clock.hpp"

///////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////

Clock::Clock()
{}

//! Construct a "Clock" instance. This constructor takes a std::vector of
//! std::strings that list the clocks that we want to keep track of over the
//! entire run. We do this so we can provide run-wide statistics at the end.
//! Clocks not in this list can still be tick-tock'ed, but the return value from
//! "elapsed" is unpredictable for those.
Clock::Clock(const std::vector<std::string> & names)
    : M_names(names)
{
    for ( const std::string & name : M_names )
        M_clock[name].second = 0.;

    // Start the global timer
    M_global_clock.restart();
}

Clock::~Clock()
{}

///////////////////////////////////////////////////////////////////////
// Public functions
///////////////////////////////////////////////////////////////////////

//! Start a clock named "name"
void Clock::tick(const std::string & name)
{
    M_clock[name].first.restart();
}

//! Stop a clock named "name" and returned elapsed time since the last tick
const double Clock::tock(const std::string & name)
{
    double tock = M_clock[name].first.elapsed();

    M_clock[name].second += tock;

    return tock;
}

//! Return the total elapsed time of a clock named "name" since construction
const double Clock::elapsed(const std::string & name)
{
    return M_clock[name].second;
}

//! Pretty-print all the timers initialised by the constructor
const std::string Clock::print_all()
{

    // Don't write out anything if M_names is empty
    if ( M_names.size() == 0 )
        return std::string();

    double wall_time = M_global_clock.elapsed();
    double not_counted = wall_time;

    std::stringstream return_string;
    return_string << "   =====   Timer results =====   " << std::endl;
    return_string << " clock name          time spent % of total" << std::endl;

    for ( const std::string & name : M_names )
    {
        // Integer division!
        int hours   = M_clock[name].second/3600;
        int minutes = (M_clock[name].second - hours*3600)/60;
        int seconds = (M_clock[name].second - hours*3600 - minutes*60);
        double fraction = M_clock[name].second/wall_time * 100.;
        not_counted -= M_clock[name].second;
        return_string << " " 
            << std::setfill(' ') << std::setw(20) << std::left << name
            << std::setfill('0') << std::setw(2) << hours << ":"
            << std::setfill('0') << std::setw(2) << minutes << ":"
            << std::setfill('0') << std::setw(2) << seconds << "   " << fraction << std::endl;
    }

    // Not counted
    int hours   = not_counted/3600;
    int minutes = (not_counted - hours*3600)/60;
    int seconds = (not_counted - hours*3600 - minutes*60);
    double fraction = not_counted/wall_time * 100.;
    return_string << " not counted         "
        << std::setfill('0') << std::setw(2) << hours << ":"
        << std::setfill('0') << std::setw(2) << minutes << ":"
        << std::setfill('0') << std::setw(2) << seconds << "   " << fraction << std::endl;

    return return_string.str();
}




