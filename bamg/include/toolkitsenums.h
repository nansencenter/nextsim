/* \file toolkitsenums.h
 * \brief: enums that encompasse all of the toolkit APIs. This mainly maps into PETSC enums.
 */

#ifndef _TOOLKITSENUMS_H_
#define _TOOLKITSENUMS_H_

typedef enum {INS_VAL, ADD_VAL} InsMode;
typedef enum {NORM_INF,NORM_TWO,NORM_FROB} NormMode;
typedef enum {DENSE_SEQUENTIAL,SPARSE_SEQUENTIAL} MatrixType;

#endif
