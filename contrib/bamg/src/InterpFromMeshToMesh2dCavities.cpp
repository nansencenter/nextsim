/*!\file InterpFromMeshToMesh2dCavities
 */

#include "./InterpFromMeshToMesh2dCavities.h"

#include "./Mesh.h"
#include "./shared.h"
#include "./typedefs.h"
//#include "../../toolkits/toolkits.h"
//#include "../../classes/classes.h"
//#include "../../bamg/bamgobjects.h"

using namespace bamg;
using namespace std;

int InterpFromMeshToMesh2dCavities(double** pdata_interp,double* data,int N_data,
		double* surface_old, double* surface_new, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new){

	/*Output*/
	double* data_interp=NULL;

    /*Intermediary*/
    int new_nb_elements=bamgmesh_new->TrianglesSize[0];


    /*---- beginning of the detection of the cavities------*/
    /*initialize thread parameters: */
    InterpFromMeshToMesh2dCavitiesThreadStruct gate;
    gate.PreviousNumbering      =NULL;
    gate.old_elements           =NULL;
    gate.new_elements           =NULL;
    gate.size_born_cavity       =NULL;
    gate.size_dead_cavity       =NULL;
    gate.element_born_cavity    =NULL;
    gate.element_dead_cavity    =NULL;
    gate.nb_cavities            =0;
    gate.nb_new_nodes           =0;

    DetectCavities(&gate, bamgmesh_old, bamgmesh_new);

    /*---- end of the detection of the cavities------*/


	/*Initialize output*/
	data_interp=xNew<double>(new_nb_elements*N_data);

	*pdata_interp=data_interp;
	return 1;

}

int DetectCavities(InterpFromMeshToMesh2dCavitiesThreadStruct* gate, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new){

	/*---------- Input  ----------*/

	int *bamg_mesh_Triangles                 = (int *)bamgmesh_old->Triangles;
    int *bamg_mesh_NodalElementConnectivity  = (int *)bamgmesh_old->NodalElementConnectivity;
    int *bamg_mesh_NodalConnectivity         = (int *)bamgmesh_old->NodalConnectivity ;
    int *bamg_mesh_ElementConnectivity       = (int *)bamgmesh_old->ElementConnectivity;
    
    int bamg_mesh_Nn                         = (int)bamgmesh_old->VerticesSize[0];
    int bamg_mesh_NNodalConnectivity         = (int)bamgmesh_old->NodalConnectivitySize[0];
    int bamg_mesh_NNodalElementConnectivity  = (int)bamgmesh_old->NodalElementConnectivitySize[0];
    int bamg_mesh_Ne                         = (int)bamgmesh_old->TrianglesSize[0];
    
    int *new_bamg_mesh_Triangles             = (int *)bamgmesh_new->Triangles;
    int *new_bamg_mesh_ElementConnectivity   = (int *)bamgmesh_new->ElementConnectivity;
    int *new_bamg_mesh_PreviousNumbering     = (int *)bamgmesh_new->PreviousNumbering;
    
    int new_bamg_mesh_Nn                     = (int)bamgmesh_new->VerticesSize[0];
    int new_bamg_mesh_NVerticesOnGeomVertex  = (int)bamgmesh_new->VerticesOnGeomVertexSize[0];
    int new_bamg_mesh_Ne                     = (int)bamgmesh_new->TrianglesSize[0];


    #if 0
    /*---------- Local variables ----------*/
    int i, j, k, l, l0, l1, l2, j_next;
    
    int i_intersect, number_of_cavities, max_size_born_cavity, max_size_dead_cavity;
    
    int i_born, j_born, i_dead, j_dead;
    int edge_node_1, edge_node_2;
    int Previous_edge_node_1, Previous_edge_node_2;
    int Previous_connected_element, tmp_born, tmp_dead;
    int Previous_num_node[3], Previous_connected_element_num_node[3];
    int max_candidate_element_to_treat=10000, number_candidate_element_to_treat, number_candidate_element_treated;    
    int candidate_element_to_treat[max_candidate_element_to_treat];
    int candidate_dead_element[2], nb_candidate_dead_element;
    int nb_dead, nb_born, new_cavity_number, neighbour, nb_connected_nodes;
    int i_dead_element_node_1, dead_element_node_1;
    int i_dead_element_node_2, dead_element_node_2;
    int i_dead_element, max_cavity_size;
       
    int  *is_born_element_to_treat = xNew<int>(new_bamg_mesh_Ne); /* boolean to know if the born element is to be treated */      
    int  *is_dead_element_to_treat = xNew<int>(bamg_mesh_Ne);	/* boolean to know if the dead element is to be treated */ 
    
    int found_intersect;
    #endif

    int  *PreviousNumbering = xNew<int>(new_bamg_mesh_Nn);

    int i_intersect=1;

    /* reducing the old_elements and new_elements arrays */
    int  *old_elements = xNew<int>(i_intersect);
    int  *new_elements = xNew<int>(i_intersect);

    int number_of_cavities=1;

    int  *size_born_cavity = xNew<int>(number_of_cavities);
    int  *size_dead_cavity = xNew<int>(number_of_cavities);

    int max_size_born_cavity=1;
    int max_size_dead_cavity=1;

    int  *element_born_cavity = xNew<int>(number_of_cavities*max_size_born_cavity);
    int  *element_dead_cavity = xNew<int>(number_of_cavities*max_size_dead_cavity);


    gate->PreviousNumbering     =PreviousNumbering;
    gate->nb_new_nodes          =new_bamg_mesh_Nn;
    gate->old_elements          =old_elements;
    gate->new_elements          =new_elements;
    gate->size_born_cavity      =size_born_cavity;
    gate->size_dead_cavity      =size_dead_cavity;
    gate->nb_cavities           =number_of_cavities;
    gate->element_born_cavity   =element_born_cavity;
    gate->element_dead_cavity   =element_dead_cavity;

return 1;	
}
