/*!\file:  Bamgx.h
 * \brief header file for Bamg module
 */

#ifndef _BAMGX_H
#define _BAMGX_H

#include "./Mesh.h"
#include "./Geometry.h"
#include "./Metric.h"

class BamgMesh;
class BamgGeom;
class BamgOpts;

int Bamgx(BamgMesh* bamgmesh_out,BamgGeom* bamggeom_out,BamgMesh* bamgmesh_in,BamgGeom* bamggeom_in,BamgOpts* bamgopts);

#endif  /* _BAMGX_H */
