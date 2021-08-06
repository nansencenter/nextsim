/*!\file:  InterpFromMeshToMesh2dx_apply.h
 * \brief header file for Bamg module
 */

#ifndef _INTERPFROMMESHTOMESH2DX_APPLY_H
#define _INTERPFROMMESHTOMESH2DX_APPLY_H

//#include "../../classes/classes.h"
#ifndef _INTERPFROMMESHTOMESH2DX_H
#include "./Options.h"
#endif
#include <omp.h>

int InterpFromMeshToMesh2dx_apply(double** pdata_interp,
      std::vector<std::vector<double>>& areacoord, std::vector<std::vector<int>>& vertex,
      std::vector<int> it,
      int nods_data,int nels_data,
      double* data,int M_data,int N_data,
      double* x_interp,double* y_interp);

#endif
