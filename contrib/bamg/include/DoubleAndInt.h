#ifndef _DOUBLEANDINT_H_
#define _DOUBLEANDINT_H_

#include "./include.h"

namespace bamg {

	class DoubleAndInt {
		//class used by Mesh::MakeQuadrangles

		public:
			double q;
			long i3j;

			//Operators
			int operator<(DoubleAndInt a){return q > a.q;}
	};
}
#endif
