#ifndef _BINARYRAND_H_
#define _BINARYRAND_H_

#include <cstdlib>

/*Return 1 or 0 randomly*/
inline int BinaryRand(){
	#ifdef RAND_MAX
		/*RAND_MAX is defined by stdlib.h and is usually 32767*/
		const long HalfRandMax = RAND_MAX/2;
		return rand() < HalfRandMax;
	#else
		/*For sun machines, RAND_MAX is not defined, use 2^24*/
		return rand() & 16384;
	#endif
} 

#endif
