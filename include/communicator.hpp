/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   communicator.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Thu Jul  2 10:37:49 2015
 */

#ifndef __Communicator_H
#define __Communicator_H 1

#include <boost/mpi.hpp>
#include <boost/timer.hpp>

namespace Nextsim
{
class Communicator : public boost::mpi::communicator
{
	typedef boost::mpi::communicator super;

public:

	typedef Communicator communicator_type;

	Communicator()
		:
		super()
	{}

    static communicator_type commSelf() { return communicator_type();}

};

} // Nextsim
#endif // __Communicator_H
