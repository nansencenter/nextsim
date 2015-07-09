/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   parameter.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Tue Jul  7 14:57:53 2015
 */

#ifndef __Parameter_H
#define __Parameter_H 1

#if !defined(BOOST_PARAMETER_MAX_ARITY)
#define BOOST_PARAMETER_MAX_ARITY 10
#endif

#include <boost/parameter.hpp>
#include <boost/type_traits.hpp>

namespace Nextsim
{
	namespace parameter = boost::parameter;

	BOOST_PARAMETER_NAME(rtolerance)
	BOOST_PARAMETER_NAME(dtolerance)
	BOOST_PARAMETER_NAME(atolerance)
	BOOST_PARAMETER_NAME(maxit)

} // Nextsim

#endif // __Parameter_H
