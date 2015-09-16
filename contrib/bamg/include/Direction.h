#ifndef _DIRECTIONO_H_
#define _DIRECTIONO_H_

#include "./typedefs.h"

namespace bamg {

	class Direction {
		private:
			Icoor1 dir;

		public:
			//Methods
			Direction();
			Direction(Icoor1 i,Icoor1 j);
			int direction(Icoor1 i,Icoor1 j);
	};
}
#endif
