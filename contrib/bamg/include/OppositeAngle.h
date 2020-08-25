#ifndef _OPPOSITEANGLE_H_
#define _OPPOSITEANGLE_H_

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

/*Return the opposite angle modulo 2 Pi*/
namespace bamg {
	inline float  OppositeAngle(float  a){return a<0 ? PI+a:a-PI;}
	inline double OppositeAngle(double a){return a<0 ?  PI+a:a- PI;}
}

#endif
