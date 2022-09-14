/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   logger.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 17 Aug 2022 07:59:06 CEST
 */

#ifndef __Logger_HPP
#define __Logger_HPP 1

#include <environment.hpp>

namespace Nextsim
{
    class Logger
    {
        public:
            Logger()
            {
                //! Get the MPI communicator
                M_comm = Environment::comm();

                //! Sets the characteristics of the output log (INFOR, WARNING, DEBUG, ERROR),
                M_log_level = Environment::logLevel();

                //! Do we output the log on all processors?
                M_log_all = Environment::logAll();
            }

            LogLevel M_log_level;
            bool M_log_all;
            Communicator M_comm;

    };
}

#endif
