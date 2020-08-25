#ifndef _BAMGTYPEDEFS_H
#define _BAMGTYPEDEFS_H

#include "./R2.h"

namespace bamg {

	/*Integer coordinates types*/
	typedef int  Icoor1; 
	typedef long long Icoor2;

	/*I2 and R2*/
	typedef P2<Icoor1,Icoor2>  I2;
	typedef P2<double,double>  R2;
}

#endif
