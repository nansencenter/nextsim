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

// Construct a "Clock" instance and start the global clock
Clock::Clock()
{
    // Just start the global clock - we can't use this->tick though
    M_names.push_back(M_global_clock);
    M_clock[M_global_clock].elapsed = 0.;
    M_clock[M_global_clock].timer.restart();

}

Clock::~Clock()
{}

///////////////////////////////////////////////////////////////////////
// Public functions
///////////////////////////////////////////////////////////////////////

//! Start a clock named "name". The "name" is added to a global list if it does
//! not already exist, and the total elapsed time is initialised.  Clocks can be
//! nested once using a "parent.child" notation.  Children report the fraction
//! of time spent within the parent's scope. There's no checking whether a child
//! is called within a parent's scope or not.
void Clock::tick(const std::string & name)
{
    // Add the name to M_names if needed and initialise the elapsed counter
    if ( M_clock.count(name) == 0 )
    {
        if ( name == M_global_clock )
            throw std::runtime_error("Clock::tick: Illegal use of reserved clock name '" + M_global_clock + "'\n");

        M_names.push_back(name);
        M_clock[name].elapsed = 0.;
    }

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

    // Sort out the familial relations and field width
    std::size_t width = 0;
    for ( const std::string & name : M_names )
    {
        M_clock[name].parent = name.substr(0, name.find("."));
        if ( M_clock[name].parent == name )
            M_clock[name].parent = M_global_clock;

        if ( M_clock.count(M_clock[name].parent) == 0 )
            throw std::logic_error("Clock::printAll: No parent found for " + name + " (assuming 'parent.child' naming convention).\n");

        width = std::max(name.length(), width);
    }

    std::stringstream return_string;
    return_string << "   =====   Timer results =====   " << std::endl;
    return_string << " " << std::setfill(' ') << std::setw(width) << std::left << "Clock name"
        << " | time spent | % of parent" << std::endl;

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
        return_string << " " << std::setfill(' ') << std::setw(width) << std::left << name << " |  "
            << std::right
            << std::setfill(' ') << std::setw(3) << hours << ":"
            << std::setfill('0') << std::setw(2) << minutes << ":"
            << std::setfill('0') << std::setw(2) << seconds << " | "
            << std::setfill(' ') << std::setw(10) << std::setprecision(2) << std::fixed << std::right << fraction << std::endl;
    }

    // Not counted
    int hours   = not_counted/3600;
    int minutes = (not_counted - hours*3600)/60;
    int seconds = (not_counted - hours*3600 - minutes*60);
    double fraction = not_counted/wall_time * 100.;
    return_string << std::setfill(' ') << std::setw(width) << " unaccounted for " << " |  "
        << std::right
        << std::setfill(' ') << std::setw(3) << hours << ":"
        << std::setfill('0') << std::setw(2) << minutes << ":"
        << std::setfill('0') << std::setw(2) << seconds << " | "
        << std::setfill(' ') << std::setw(10) << std::setprecision(2) << std::fixed << std::right << fraction << std::endl;

    return return_string.str();
}

