/*!\file InterpFromGridToMeshx.h
 * \brief: header file for Data interpolation routines.
 */

#ifndef _INTERPFROMGRIDTOMESHX_H
#define _INTERPFROMGRIDTOMESHX_H

#include "./EnumDefinitions.h"
#include "issm_threads.h"
#include "./IssmSeqVec.h"

#define _NUMTHREADS_ 4

template <class doubletype> class IssmSeqVec;

/*threading: */
typedef struct{
	double*             x;
	int                 x_rows;
	double*             y;
	int                 y_rows;
	double*             data;
	double              default_value;
	int                 interp;
	int                 M;
	int                 N;
	int                 N_data;
	int                 nods;
	double*             x_mesh;
	double*             y_mesh;
	double*             data_mesh;
	bool                row_major;
} InterpFromGridToMeshxThreadStruct;

int    InterpFromGridToMeshx(double* &data_mesh,double* x, int x_rows, double* y, int y_rows, double* data, int M, int N, int N_data, double* x_mesh, double* y_mesh, int nods, double default_value, int interpenum=BilinearInterpEnum, bool row_major=false);

int InterpFromGridToMeshxt(InterpFromGridToMeshxThreadStruct gate, double* data_mesh);
void*  InterpFromGridToMeshxt(void* vInterpFromGridToMeshxThreadStruct);
bool   findindices(int* pn,int* pm,double* x,int x_rows, double* y,int y_rows, double xgrid,double ygrid);
double triangleinterp(double x1,double x2,double y1,double y2,double Q11,double Q12,double Q21,double Q22,double x,double y);
double bilinearinterp(double x1,double x2,double y1,double y2,double Q11,double Q12,double Q21,double Q22,double x,double y);
double nearestinterp(double x1,double x2,double y1,double y2,double Q11,double Q12,double Q21,double Q22,double x,double y);

#endif /* _INTERPFROMGRIDTOMESHX_H */
