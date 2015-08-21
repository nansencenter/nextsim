#ifndef _BAMGDET_H_
#define _BAMGDET_H_

#include "./typedefs.h"

namespace bamg {

	Icoor2 inline det(const I2 &a,const I2 & b,const I2 &c){
		Icoor2 bax = b.x - a.x ,bay = b.y - a.y; 
		Icoor2 cax = c.x - a.x ,cay = c.y - a.y; 
		return  bax*cay - bay*cax;
	}

}
#endif
