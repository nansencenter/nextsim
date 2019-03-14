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
//! "elapsed" is unpredictable for those and they are not included when calling
//! printAll. Clocks can be nested once using a "parent.child" notation.
//! Children report the fraction of time spent within the parent's scope.
Clock::Clock(const std::vector<std::string> & names)
    : M_names(names)
{
    // Check for use of reserved name set by M_global_clock (currently "main")
    for ( const std::string & name : M_names )
        if ( name == M_global_clock )
            throw std::runtime_error("Clock::Clock: Illegal use of reserved clock name '" + M_global_clock + "'\n");

    // Sort the names so that children get put after parents
    std::sort(M_names.begin(), M_names.end());

    // The global clock is the first one on the list
    M_names.push_back(M_global_clock);
    std::rotate(M_names.rbegin(), M_names.rbegin() + 1, M_names.rend());

    // Initialise the elapsed time and set parents
    for ( const std::string & name : M_names )
    {
        M_clock[name].elapsed = 0.;
        M_clock[name].parent = name.substr(0, name.find("."));
        if ( M_clock[name].parent == name )
            M_clock[name].parent = M_global_clock;
    }

    // Start the global timer
    this->tick(M_global_clock);
}

Clock::~Clock()
{}

///////////////////////////////////////////////////////////////////////
// Public functions
///////////////////////////////////////////////////////////////////////

//! Start a clock named "name"
void Clock::tick(const std::string & name)
{
    M_clock[name].timer.restart();
}

//! Stop a clock named "name" and returned elapsed time since the last tick.
//! Add to the total elapsed time.
const double Clock::tock(const std::string & name)
{
    double tock = M_clock[name].timer.elapsed();

    M_clock[name].elapsed += tock;

    return tock;
}

//! Return the total elapsed time of a clock named "name"
const double Clock::elapsed(const std::string & name)
{
    return M_clock[name].elapsed;
}

//! Pretty-print all the clocks initialised by the constructor
const std::string Clock::printAll()
{

    // Don't write out anything if M_names is empty
    if ( M_names.size() == 0 )
        return std::string();

    double wall_time = this->tock(M_global_clock);
    double not_counted = wall_time;

    std::stringstream return_string;
    return_string << "   =====   Timer results =====   " << std::endl;
    return_string << " clock name                    | time spent | % of parent" << std::endl;

    for ( const std::string & name : M_names )
    {
        // Integer division!
        double elapsed = M_clock[name].elapsed;
        std::string parent = M_clock[name].parent;
        int hours   = elapsed/3600;
        int minutes = (elapsed - hours*3600)/60;
        int seconds = (elapsed - hours*3600 - minutes*60);
        double fraction = elapsed / M_clock[parent].elapsed * 100.;
        if ( parent == M_global_clock && name != M_global_clock )
            not_counted -= elapsed;
        return_string << " " 
            << std::setfill(' ') << std::setw(30) << std::left << name << "|  "
            << std::setfill('0') << std::setw(2) << hours << ":"
            << std::setfill('0') << std::setw(2) << minutes << ":"
            << std::setfill('0') << std::setw(2) << seconds << "  | "
            << std::setfill(' ') << std::setw(10) << std::setprecision(2) << std::fixed << std::right << fraction << std::endl;
    }

    // Not counted
    int hours   = not_counted/3600;
    int minutes = (not_counted - hours*3600)/60;
    int seconds = (not_counted - hours*3600 - minutes*60);
    double fraction = not_counted/wall_time * 100.;
    return_string << " not counted                   |  "
        << std::setfill('0') << std::setw(2) << hours << ":"
        << std::setfill('0') << std::setw(2) << minutes << ":"
        << std::setfill('0') << std::setw(2) << seconds << "  | "
        << std::setfill(' ') << std::setw(10) << std::setprecision(2) << std::fixed << std::right << fraction << std::endl;

    return return_string.str();
}




