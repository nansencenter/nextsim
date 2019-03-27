/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   timer.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 13 Mar 2019 09:56:02 CET
 */

#include "timer.hpp"

///////////////////////////////////////////////////////////////////////
// Constructor and destructor
///////////////////////////////////////////////////////////////////////

// Construct a "Timer" instance and start the global timer
Timer::Timer()
{
    // Just start the global timer
    // We can't use this->tick though, because there's a check preventing a
    // call to tick with M_global_timer as an argument
    M_names.push_back(M_global_timer);
    M_lineage.push_back(M_global_timer);

    M_timer[M_global_timer].elapsed = 0.;
    M_timer[M_global_timer].parent = M_global_timer;
    M_timer[M_global_timer].generation = 0;
    M_timer[M_global_timer].running = true;
    M_timer[M_global_timer].timer.restart();

}

Timer::~Timer()
{}

///////////////////////////////////////////////////////////////////////
// Public functions
///////////////////////////////////////////////////////////////////////

//! (Re)start a timer named "name". The "name" is added to a global list if it
//! does not already exist, the total elapsed time is initialised, and the
//! parent is registered. Children report the fraction of time spent within the
//! parent's scope. There is no checking of correct use; in particular, the
//! "name" needs to be globally unique, otherwise an existing timer (with that
//! name) is used, regardless of lineage.
void Timer::tick(const std::string & name)
{
    // Add the name to M_names if needed and initialise the elapsed counter
    if ( M_timer.count(name) == 0 )
    {
        if ( name == M_global_timer )
            throw std::runtime_error("Timer::tick: Illegal use of reserved timer name '" + M_global_timer + "'\n");

        M_names.push_back(name);
        M_timer[name].elapsed = 0.;
        M_timer[name].parent = M_lineage.back();
        M_timer[name].generation = M_lineage.size();
        M_timer[name].running = false;
    }

    // Check for consistency
    if ( M_timer[name].running )
        throw std::logic_error("Timer:tick: Timer "+name+" is already running.");

    // Keep track of lineage
    M_lineage.push_back(name);

    M_timer[name].running = true;

    M_timer[name].timer.restart();
}

//! Stop a timer named "name" and add to the total elapsed time
void Timer::tock(const std::string & name)
{
    // Check for consistency
    if ( !M_timer[name].running )
        throw std::runtime_error("Timer:tock: Timer "+name+" is not running.");

    M_timer[name].lap = M_timer[name].timer.elapsed();

    M_timer[name].elapsed += M_timer[name].lap;

    M_lineage.pop_back();

    // Now repeated calls to tock won't mess everything up (but who would anyone do this?)
    M_timer[name].timer.restart();

    M_timer[name].running = false;
}

//! For a timer named "name" returned elapsed time since the last tick of a
//! running timer or the length of the last tick-tock interval of a stopped timer
const double Timer::lap(const std::string & name)
{
    if ( M_timer[name].running )
        return M_timer[name].timer.elapsed();
    else
        return M_timer[name].lap;
}

//! Return the total elapsed time of a timer named "name"
const double Timer::elapsed(const std::string & name)
{
    return M_timer[name].elapsed;
}

//! Pretty-print all the timers initialised by the constructor
const std::string Timer::printAll()
{
    // Don't write out anything if M_names is empty
    if ( M_names.size() == 0 )
        return std::string();

    // Wall timer
    this->tock(M_global_timer);
    double wall_time = this->elapsed(M_global_timer);
    double not_counted = wall_time;

    // Column width
    std::size_t width = 15; // Minimum width of the name column
    int padding = 2; // Multiplier for the indent for each generation
    for ( const std::string & name : M_names )
        width = std::max(name.length() + padding*M_timer[name].generation, width);

    // Actual output
    std::stringstream return_string;
    return_string << "   =====   Timer results =====   " << std::endl;
    return_string << " " << std::setfill(' ') << std::setw(width) << std::left << "Timer name"
        << " | time spent | % of parent | % of total" << std::endl;

    for ( const std::string & name : M_names )
    {
        // Just here for readability
        double elapsed = this->elapsed(name);
        std::string parent = M_timer[name].parent;

        // Integer division!
        int hours   = elapsed/3600;
        int minutes = (elapsed - hours*3600)/60;
        int seconds = (elapsed - hours*3600 - minutes*60);
        double fraction_parent = elapsed / this->elapsed(parent) * 100.;
        double fraction_total  = elapsed / wall_time * 100.;

        // We don't subtract children's time from the total
        if ( parent == M_global_timer &&  name != M_global_timer )
                not_counted -= elapsed;

        int indent = padding*M_timer[name].generation;
        return_string
            << std::setfill(' ') << std::setw(indent+1) << " "
            << std::setfill(' ') << std::setw(width-indent) << std::left << name << " | "
            << std::right
            << std::setfill(' ') << std::setw(4) << hours << ":"
            << std::setfill('0') << std::setw(2) << minutes << ":"
            << std::setfill('0') << std::setw(2) << seconds << " | ";

        if ( parent == M_global_timer &&  name != M_global_timer )
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

