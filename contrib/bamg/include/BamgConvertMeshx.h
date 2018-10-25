/*!\file:  BamgConvertMeshx.h
 * \brief header file for Bamg module
 */

#ifndef _BAMGCONVERTMESHX_H
#define _BAMGCONVERTMESHX_H

#include "./Mesh.h"
#include "./Geometry.h"
#include "./Metric.h"

class BamgMesh;
class BamgGeom;
class BamgOpts;

/* local prototypes: */
int BamgConvertMeshx(BamgMesh* bamgmesh,BamgGeom* bamggeom,int* index,double* x,double* y,int nods,int nels);
int BamgConvertMeshx(BamgMesh* bamgmesh,BamgGeom* bamggeom,int* index,double* x,double* y,int nods,int nels, int nods_local, int nels_local);

#endif
