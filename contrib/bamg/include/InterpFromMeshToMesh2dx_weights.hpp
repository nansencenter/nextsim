/*!\file:  InterpFromMeshToMesh2dx_weights.h
 * \brief header file for Bamg module
 */

#ifndef _INTERPFROMMESHTOMESH2DX_WEIGHTS_H
#define _INTERPFROMMESHTOMESH2DX_WEIGHTS_H

//#include "../../classes/classes.h"
#ifndef _INTERPFROMMESHTOMESH2DX_H
#include "./Options.h"
#endif
#include <omp.h>

int InterpFromMeshToMesh2dx_weights(
      std::vector<std::vector<double>>& areacoord, std::vector<std::vector<int>>& vertex,
      std::vector<int> it,
      int* index_data,double* x_data,double* y_data,
      int nods_data,int nels_data,
      int M_data,
      const double* x_interp, const double* y_interp,int N_interp);

#endif
