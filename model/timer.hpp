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

/**
 * @class Timer
 *
 * @brief A class to provide named timers.
 *
 * @details This class should can be used as a singleton (although it doesn't
 * 			have to be) to provide named timers throughout an application. It
 * 			not only provides start (\ref tick()) and stop (\ref tock())
 * 			actions but also a stopwatch-like \ref lap() action.
 *
 */
class Timer
{
    public:
        ///////////////////////////////////////////////////////////////////////
        // Constructors (and destructor)
        ///////////////////////////////////////////////////////////////////////
        
	    //! Construct a "Timer" instance and start the global timer
	    Timer();

        ~Timer();

        ///////////////////////////////////////////////////////////////////////
        // Public functions
        ///////////////////////////////////////////////////////////////////////
        
        //! @brief Starts or restarts a named timer.
        //!
        //! @details (Re)starts a timer named "name". The "name" is added to a
        //! 		 global list if it does not already exist, the total
        //!          elapsed time is initialised, and the parent is registered.
        //!          Children report the fraction of time spent within the
        //!          parent's scope. There is no checking of correct use; in
        //!          particular, the "name" needs to be globally unique,
        //!          otherwise an existing timer (with that name) is used,
        //!          regardless of lineage.
        //! @param  name of the timer to be (re)started.
        void tick(const std::string & name);
        //! @brief Stops a timer named "name" and adds to the total elapsed time.
        //! @param  name of the timer to be stopped.
        void tock(const std::string & name);
        //! @brief Returns the current time of a named timer.
        //!
        //! @details For a timer named "name" returned elapsed time since the
        //!			 last tick of a running timer or the length of the last
        //! 		 tick-tock interval of a stopped timer.
        //! @param  name of the timer to be interrogated.
        //! @returns elapsed time since the last tick or tock.
        const double lap(const std::string & name);
        //! @brief Return the total elapsed time of a named timer.
        //! @param  name of the timer to be interrogated.
        //! @returns the total elapsed time.
        const double elapsed(const std::string & name);
        //! @brief Pretty-print all the timers.
        //! @returns pretty printed information as a string.
        const std::string printAll();

    private:
        // Pretty-print all children of a given timer - and recursively its children as well
        std::string printChildren(const std::string & parent, const double wall_time);
        // Pretty-print a given timer
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
