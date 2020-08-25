#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "./Direction.h"
#include "./macros.h"

namespace bamg {

	/*Constructors/Destructors*/
	Direction::Direction():/*{{{*/
		dir(MaxICoor){

	}/*}}}*/
	Direction::Direction(Icoor1 i,Icoor1 j) {/*{{{*/
		Icoor2 n2 = 2*(Abs(i)+Abs(j));  
		Icoor2 r  = MaxICoor* (Icoor2) i;
		Icoor1 r1 = (Icoor1) (2*(r/ n2)); // odd number 
		dir = (j>0) ? r1 : r1+1;          // odd-> j>0 even-> j<0
	}/*}}}*/

	/*Methods*/
	int Direction::direction(Icoor1 i,Icoor1 j) {/*{{{*/
		int r =1; 
		if (dir!= MaxICoor) {
			Icoor2 x(dir/2),y1(MaxICoor/2-Abs(x)),y((dir%2)?-y1:y1);
			r = (x*i + y*j) >=0;
		}
		return r;
	}
	/*}}}*/

} 
