/*
 * Physics1D.hpp
 *
 *  Created on: 6 Aug 2021
 *      Author: Tim Spain, <timothy.spain@nersc.no>
 */

#ifndef MODEL_PHYSICS1D_HPP_
#define MODEL_PHYSICS1D_HPP_

#include "finiteelement.hpp"

namespace Nextsim {
class Physics1D {
public:
	Physics1D();

	void setFE(FiniteElement & fe);

private:
	FiniteElement& fe;
}; // class Physics1D
} // namespace Nextsim
#endif /* MODEL_PHYSICS1D_HPP_ */
