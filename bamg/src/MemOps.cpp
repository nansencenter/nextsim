/*
 * MemOps.cpp
 *
 *  Created on: Sep 10, 2013
 *      Author: utke
 */

#include "MemOps.h"

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
template <> adouble*  xNew(unsigned int size) {
	ensureContiguousLocations(size);
	adouble* aT_p=new adouble[size];
	assert(aT_p);
	return aT_p;
}
#endif
