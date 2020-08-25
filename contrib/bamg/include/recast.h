/*
 * recast.h
 *
 *  Created on: Jun 26, 2012
 *      Author: utke
 */

#ifndef _RECAST_H_
#define _RECAST_H_

#if !defined(_HAVE_ADOLC_) || defined(_WRAPPERS_)

template<class To, class From> To reCast(const From& from) {
    return (To)from;
}

#else
#include "./types.h"

template<class To, class From> struct ForPartialSpecialization {
    static  To reCast(const From& from ) { return (To) from;}
};

template<class To, class From> To reCast(const From& from) {
    return ForPartialSpecialization<To,From>::reCast(from);
}
/**
 * partial specialization
 */
template<class To> struct ForPartialSpecialization<To,adouble> {
    static  To reCast(const adouble& from ) { return (To) (from.getValue());}
};

#endif

/*Morlighem's change: we do not want dynamic_casts because of performance
 * issue, so for now, we just use C-like cast*/
template<class To,class From> To xDynamicCast(const From& from) {

	/*C-like cast (fast but not safe)*/
	return (To) from;

	/*C++ dynamic_cast, poor performance but safer*/
	//return dynamic_cast<To>(from);
}

#endif
