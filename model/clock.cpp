/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   clock.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 13 Mar 2019 09:56:02 CET
 */

#include "clock.hpp"

///////////////////////////////////////////////////////////////////////
// Constructor and destructor
///////////////////////////////////////////////////////////////////////

// Construct a "Clock" instance and start the global clock
Clock::Clock()
{
    // Just start the global clock
    // We can't use this->tick though, because there's a check preventing a
    // call to tick with M_global_clock as an argument
    M_names.push_back(M_global_clock);
    M_lineage.push_back(M_global_clock);

    M_clock[M_global_clock].elapsed = 0.;
    M_clock[M_global_clock].parent = M_global_clock;
    M_clock[M_global_clock].generation = 0;
    M_clock[M_global_clock].running = true;
    M_clock[M_global_clock].timer.restart();

}

Clock::~Clock()
{}

///////////////////////////////////////////////////////////////////////
// Public functions
///////////////////////////////////////////////////////////////////////

//! (Re)start a clock named "name". The "name" is added to a global list if it
//! does not already exist, the total elapsed time is initialised, and the
//! parent is registered. Children report the fraction of time spent within the
//! parent's scope. There is no checking of correct use; in particular, the
//! "name" needs to be globally unique, otherwise an existing clock (with that
//! name) is used, regardless of lineage.
void Clock::tick(const std::string & name)
{
    // Add the name to M_names if needed and initialise the elapsed counter
    if ( M_clock.count(name) == 0 )
    {
        if ( name == M_global_clock )
            throw std::runtime_error("Clock::tick: Illegal use of reserved clock name '" + M_global_clock + "'\n");

        M_names.push_back(name);
        M_clock[name].elapsed = 0.;
        M_clock[name].parent = M_lineage.back();
        M_clock[name].generation = M_lineage.size();
        M_clock[name].running = false;
    }

    // Check for consistency
    if ( M_clock[name].running )
        throw std::logic_error("Clock:tick: Clock "+name+" is already running.");

    // Keep track of lineage
    M_lineage.push_back(name);

    M_clock[name].running = true;

    M_clock[name].timer.restart();
}

//! Stop a clock named "name" and add to the total elapsed time
void Clock::tock(const std::string & name)
{
    // Check for consistency
    if ( !M_clock[name].running )
        throw std::runtime_error("Clock:tock: Clock "+name+" is not running.");

    M_clock[name].lap = M_clock[name].timer.elapsed();

    M_clock[name].elapsed += M_clock[name].lap;

    M_lineage.pop_back();

    // Now repeated calls to tock won't mess everything up (but who would anyone do this?)
    M_clock[name].timer.restart();

    M_clock[name].running = false;
}

//! For a clock named "name" returned elapsed time since the last tick of a
//! running clock or the length of the last tick-tock interval of a stopped clock
const double Clock::lap(const std::string & name)
{
    if ( M_clock[name].running )
        return M_clock[name].timer.elapsed();
    else
        return M_clock[name].lap;
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

    // Wall clock
    this->tock(M_global_clock);
    double wall_time = this->elapsed(M_global_clock);
    double not_counted = wall_time;

    // Column width
    std::size_t width = 15; // Minimum width of the name column
    int padding = 2; // Multiplier for the indent for each generation
    for ( const std::string & name : M_names )
        width = std::max(name.length() + padding*M_clock[name].generation, width);

    // Actual output
    std::stringstream return_string;
    return_string << "   =====   Timer results =====   " << std::endl;
    return_string << " " << std::setfill(' ') << std::setw(width) << std::left << "Clock name"
        << " | time spent | % of parent | % of total" << std::endl;

    for ( const std::string & name : M_names )
    {
        // Just here for readability
        double elapsed = this->elapsed(name);
        std::string parent = M_clock[name].parent;

        // Integer division!
        int hours   = elapsed/3600;
        int minutes = (elapsed - hours*3600)/60;
        int seconds = (elapsed - hours*3600 - minutes*60);
        double fraction_parent = elapsed / this->elapsed(parent) * 100.;
        double fraction_total  = elapsed / wall_time * 100.;

        // We don't subtract children's time from the total
        if ( parent == M_global_clock &&  name != M_global_clock )
                not_counted -= elapsed;

        int indent = padding*M_clock[name].generation;
        return_string
            << std::setfill(' ') << std::setw(indent+1) << " "
            << std::setfill(' ') << std::setw(width-indent) << std::left << name << " | "
            << std::right
            << std::setfill(' ') << std::setw(4) << hours << ":"
            << std::setfill('0') << std::setw(2) << minutes << ":"
            << std::setfill('0') << std::setw(2) << seconds << " | ";

        if ( parent == M_global_clock &&  name != M_global_clock )
            return_string << std::setfill(' ') << std::setw(11) << " ";
        else
            return_string << std::setfill(' ') << std::setw(11) << std::setprecision(2) << std::fixed
                << fraction_parent;

        return_string << " | "
            << std::setfill(' ') << std::setw(10) << std::setprecision(2) << std::fixed
            << fraction_total << std::endl;
    }

    // Not counted
    int hours   = not_counted/3600;
    int minutes = (not_counted - hours*3600)/60;
    int seconds = (not_counted - hours*3600 - minutes*60);
    double fraction = not_counted/wall_time * 100.;
    return_string << " " << std::setfill(' ') << std::setw(width) << std::left << "Unaccounted for" << " | "
        << std::right
        << std::setfill(' ') << std::setw(4) << hours << ":"
        << std::setfill('0') << std::setw(2) << minutes << ":"
        << std::setfill('0') << std::setw(2) << seconds << " | "
        << std::setfill(' ') << std::setw(11) << " " << " | "
        << std::setfill(' ') << std::setw(10) << std::setprecision(2) << std::fixed << fraction << std::endl;

    return return_string.str();
}

