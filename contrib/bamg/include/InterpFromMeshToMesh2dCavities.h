/**
 * @file   InterpFromMeshToMesh2dCavities.h
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Wed Dec 23 15:00:00 2015
 */

#ifndef _INTERPFROMMESHTOMESH2DCAVITIES_H
#define _INTERPFROMMESHTOMESH2DCAVITIES_H

#include "./Mesh.h"
//#include "../../classes/classes.h"
#include "./IssmSeqVec.h"
#include <omp.h>

typedef unsigned long ulong;

/*threading: */
typedef struct{
	double* PreviousNumbering;
	double* old_elements;
	double* new_elements;
	int* 	size_born_cavity;
	int* 	size_dead_cavity;
	double* element_born_cavity;
	double* element_dead_cavity;
	int 	new_bamg_mesh_Nn;
	int 	nb_matching_elements;
	int     nb_cavities;
	int  	max_size_born_cavity;
    int  	max_size_dead_cavity;
} InterpFromMeshToMesh2dCavitiesThreadStruct;

#define next_node(a) ( ( (a) == 0 )  ?  1   : ( ( (a) == 1 ) ?  2   : 0 ) )
#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )

int InterpFromMeshToMesh2dCavities(double** pdata_interp,double* data, int* method,int N_data,
			double* surface_old, double* surface_new, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new);

int DetectCavities(InterpFromMeshToMesh2dCavitiesThreadStruct* gate, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new);

int InterpCavity(double* tmp_mean_variables, double* tmp_integrated_area,
			int nb_dead_elements, int nb_born_elements, int nb_variables,
			double* dead_elements, double* born_elements, double* PreviousNumbering,
			double* IntMatrix_in, int* method_in, BamgMesh* bamgmesh_old, BamgMesh* bamgmesh_new, int debug_born_elements_i);

#endif
