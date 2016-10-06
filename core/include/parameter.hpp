/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   parameter.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Tue Jul  7 14:57:53 2015
 */

#ifndef __Parameter_H
#define __Parameter_H 1

#if !defined(BOOST_PARAMETER_MAX_ARITY)
#define BOOST_PARAMETER_MAX_ARITY 12
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
    BOOST_PARAMETER_NAME(reuse_prec)
    BOOST_PARAMETER_NAME(rebuild)

    BOOST_PARAMETER_NAME(matrix)
    BOOST_PARAMETER_NAME(solution)
    BOOST_PARAMETER_NAME(rhs)
    BOOST_PARAMETER_NAME(ksp)
    BOOST_PARAMETER_NAME(pc)
    BOOST_PARAMETER_NAME(pcfactormatsolverpackage)

} // Nextsim

#endif // __Parameter_H
