/*!\file:  InterpFromMeshToMesh2dCavities.h
 * \brief header file for Bamg module
 */

#ifndef _INTERPFROMMESHTOMESH2DCAVITIES_H
#define _INTERPFROMMESHTOMESH2DCAVITIES_H

#include "./Mesh.h"
//#include "../../classes/classes.h"
#include "./IssmSeqVec.h"

typedef unsigned long ulong;

/*threading: */
typedef struct{
	ulong* 	PreviousNumbering;
	ulong* 	old_elements;
	ulong* 	new_elements;
	ulong* 	size_born_cavity;
	ulong* 	size_dead_cavity;
	ulong* 	element_born_cavity;
	ulong* 	element_dead_cavity;
	ulong 	new_bamg_mesh_Nn;
	ulong 	nb_matching_elements;
	ulong   nb_cavities;
	ulong	max_size_born_cavity;
    ulong	max_size_dead_cavity;
} InterpFromMeshToMesh2dCavitiesThreadStruct;

#define next(a) ( ( (a) == 0 )  ?  1   : ( ( (a) == 1 ) ?  2   : 0 ) )

int InterpFromMeshToMesh2dCavities(double** pdata_interp,double* data,int N_data,
			double* surface_old, double* surface_new, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new);

int DetectCavities(InterpFromMeshToMesh2dCavitiesThreadStruct* gate, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new);

#endif