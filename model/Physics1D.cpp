/*
 * Physics1D.cpp
 *
 *  Created on: 6 Aug 2021
 *      Author: Tim Spain, <timothy.spain@nersc.no>
 */

#include "Physics1D.hpp"

namespace Nextsim {
Physics1D::Physics1D():
fe(0)
{
} // Physics1D::Physics1D()

void Physics1D::setFE(FiniteElement& fe) {
	this->fe = fe;
} // void Physics1D::setFE(FiniteElement&)


} // namespace Nextsim
