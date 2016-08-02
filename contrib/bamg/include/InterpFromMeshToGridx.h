/*!\file InterpFromMeshToGridx.h
 * \brief: header file for Data interpolation routines.
 */

#ifndef _INTERPFROMMESHTOGRIDX_H
#define _INTERPFROMMESHTOGRIDX_H

void InterpFromMeshToGridx(double* &pgriddata, int* index_mesh, double* x_mesh, double* y_mesh, int nods, int nels, double* data_mesh, int data_length, int N_data, double xmin, double ymax, double xposting, double yposting, int nlines, int ncols, double default_value);

#endif /* _INTERPFROMMESHTOGRIDX_H */
