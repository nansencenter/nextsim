/*!\file:  isnan.h
 * \brief: only used for intel compiler.
 */

#ifndef _XISNAN_H_
#define _XISNAN_H_

/*If include cmath instead of math, isnan on linux64 murdo does not work: */
#include <math.h>

template <class T> int xIsNan(const T& X) {
#ifdef _INTEL_WIN_
		return (X!=X)?1:0;
#else
		return isnan(X);
#endif
}

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
#include "./types.h"
template <> int xIsNan<adouble> (const adouble& X);
#endif

#endif
