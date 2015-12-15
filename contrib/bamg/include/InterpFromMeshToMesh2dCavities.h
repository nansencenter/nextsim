/*!\file:  InterpFromMeshToMesh2dCavities.h
 * \brief header file for Bamg module
 */

#ifndef _INTERPFROMMESHTOMESH2DCAVITIES_H
#define _INTERPFROMMESHTOMESH2DCAVITIES_H

#include "./Mesh.h"
//#include "../../classes/classes.h"
#include "./IssmSeqVec.h"

/*threading: */
typedef struct{
	int* 	PreviousNumbering;
	int 	nb_new_nodes;
	int* 	old_elements;
	int* 	new_elements;
	int* 	size_born_cavity;
	int* 	size_dead_cavity;
	int 	nb_cavities;
	int* 	element_born_cavity;
	int* 	element_dead_cavity;
} InterpFromMeshToMesh2dCavitiesThreadStruct;

int InterpFromMeshToMesh2dCavities(double** pdata_interp,double* data,int N_data,
			double* surface_old, double* surface_new, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new);

int DetectCavities(InterpFromMeshToMesh2dCavitiesThreadStruct* gate, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new);

#endif