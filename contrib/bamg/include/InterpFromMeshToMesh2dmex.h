/*!\file:  InterpFromMeshToMesh2dx.h
 * \brief header file for Bamg module
 */

#ifndef _INTERPFROMMESHTOMESH2DMEX_H
#define _INTERPFROMMESHTOMESH2DMEX_H

//#include "../../classes/classes.h"
#include "./Options.h"
//#include <omp.h>

int InterpFromMeshToMesh2dmex(double** pdata_interp,int* index_data,double* x_data,double* y_data,int nods_data,int nels_data,
                            double* data,int M_data,int N_data,double* x_interp,double* y_interp,int N_interp,bool isdefault=false, double defaultvalue=1e-24 /*Options* options*/);

#endif
