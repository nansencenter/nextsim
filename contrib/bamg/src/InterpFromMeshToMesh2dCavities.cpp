/**
 * @file   InterpFromMeshToMesh2dCavities.cpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Wed Dec 23 15:00:00 2015
 */

#include "./InterpFromMeshToMesh2dCavities.h"

#include "./Mesh.h"
#include "./shared.h"
#include "./typedefs.h"
#include <cmath>
//#include "../../toolkits/toolkits.h"
//#include "../../classes/classes.h"
//#include "../../bamg/bamgobjects.h"

using namespace bamg;
using namespace std;



int InterpFromMeshToMesh2dCavities(double** pdata_interp,double* IntMatrix_in,
      int* method_in, int nb_variables, double* surface_old,
      double* surface_new, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new)
{

    /*  Conservative interpolation method
        Use the cavities detected by detect_cavity to transfer physical
        quantities from the old element to the new ones*/

    int verbosity;

    /*Bamg options*/
    verbosity=0;//bamgopts->verbose;

    /* Which conservative method to be used:
      method 0 is simply averaging over the whole cavity (conservative but highly diffusive)
      method 1 computes the intersection between born and dead triangles (conservative and non-diffusive)*/
    int method=1; /* method 1 should be used, method 0 is kept for comparison and debugging*/

    /*Output*/
    double* IntMatrix_out=NULL;
    double* AgeVector_out=NULL;

    /*Intermediary*/
    int new_nb_elements=bamgmesh_new->TrianglesSize[0];

    /*Initialize output*/
    IntMatrix_out=xNew<double>(new_nb_elements*nb_variables);
    AgeVector_out=xNewZeroInit<double>(new_nb_elements);

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
    gate.new_bamg_mesh_Nn       =0;
    gate.nb_matching_elements   =0;
    gate.max_size_born_cavity   =0;
    gate.max_size_dead_cavity   =0;

    DetectCavities(&gate, bamgmesh_old, bamgmesh_new);
    /*---- end of the detection of the cavities        ------*/

    /*---- beginning of the interpolation cavity by cavity------*/
    if (verbosity>1) _printf_("nb_cavities" << gate.nb_cavities   << "\n");
    if (verbosity>1) _printf_("gate.new_bamg_mesh_Nn" << gate.new_bamg_mesh_Nn << "\n");
    if (verbosity>1) _printf_("gate.nb_matching_elements" << gate.nb_matching_elements << "\n");
    if (verbosity>1) _printf_("gate.max_size_born_cavity" << gate.max_size_born_cavity << "\n");
    if (verbosity>1) _printf_("gate.max_size_dead_cavity" << gate.max_size_dead_cavity  << "\n");

    if (verbosity>1) _printf_("   Interp Cavities: copy...\n");

    /* Copy the values of the old to the corresponding new elements*/
    double new_element_i;
    double old_element_i;
    for(int i=0; i<gate.nb_matching_elements; i++)
    {
        new_element_i=gate.new_elements[i];
        old_element_i=gate.old_elements[i];

        for(int j=0; j<nb_variables; j++)
            IntMatrix_out[((int)new_element_i-1)*nb_variables+j]=IntMatrix_in[((int)old_element_i-1)*nb_variables+j];
    }

    if (verbosity>1) _printf_("   Interp Cavities: initialize intermediary...\n");
    double  *dead_elements = xNew<double>(gate.max_size_dead_cavity);
    double  *born_elements = xNew<double>(gate.max_size_born_cavity);
    if (verbosity>1) _printf_("nb_dead_elements, nb_born_elements, tmp_element" << "\n");
    int nb_dead_elements, nb_born_elements;
    double tmp_element;
    double integrated_area, new_integrated_area, correction_factor;

    // integrated_variables is used with method 0
    double *integrated_variables = xNew<double>(nb_variables);

    // tmp_mean_variables and tmp_integrated_area are used with method 1
    double *tmp_mean_variables  = xNew<double>(gate.max_size_born_cavity*nb_variables);
    double *tmp_integrated_area = xNew<double>(gate.max_size_born_cavity);

    double area_error_tolerance=1000;
    double area_error;

    if (verbosity>1) _printf_("   Interp Cavities: loop1...\n");

#ifndef BAMG_NO_OMP
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
#pragma omp parallel for num_threads(max_threads) private(thread_id)
#endif
    for (int cavity_id=0; cavity_id < gate.nb_cavities; ++cavity_id)
    {
        if (verbosity>1) _printf_("   Interp Cavities: enter loop1...\n");
        nb_dead_elements= gate.size_dead_cavity[cavity_id];
        nb_born_elements= gate.size_born_cavity[cavity_id];

        if (nb_dead_elements<1)  _error_("nb_dead_elements should not be larger than 1. Here it is: " << nb_dead_elements<<" ");
        if (nb_born_elements<1)  _error_("nb_born_elements should not be larger than 1. Here it is: " << nb_born_elements<<" ");

        if (verbosity>1) _printf_("   Interp Cavities: init...\n");
        /* Initialization */
        integrated_area=0.;
        new_integrated_area=0.;
        for(int i=0; i<nb_variables; i++)
            integrated_variables[i]=0.;

        for(int i=0; i<gate.max_size_born_cavity; i++)
        {
            tmp_integrated_area[i] = 0.;
            for(int j=0; j<nb_variables; j++)
                tmp_mean_variables[i*nb_variables+j]=0.;
        }

        if (verbosity>1) _printf_("   Interp Cavities: integrated area, nb_dead_elements:" << nb_dead_elements << "\n");

        for(int i=0; i<nb_dead_elements; i++)
        {
            tmp_element=gate.element_dead_cavity[cavity_id*gate.max_size_dead_cavity+i];
            if (verbosity>1) _printf_("tmp_element: " << tmp_element << "\n");

            dead_elements[i]=tmp_element;
            integrated_area+=surface_old[(int)tmp_element-1];
        }

        if (verbosity>1) _printf_("   Interp Cavities: integrated area...\n");
        for(int i=0; i<nb_born_elements; i++)
        {
            tmp_element=gate.element_born_cavity[cavity_id*gate.max_size_born_cavity+i];
            born_elements[i]=tmp_element;
            new_integrated_area+=surface_new[(int)tmp_element-1];
        }

        if (new_integrated_area<=0.)  _error_("new_integrated_area should be larger than 0. Here it is: " << new_integrated_area<<" ");
        if (integrated_area<=0.)  _error_("integrated_area should be larger than 0. Here it is: " << integrated_area<<" ");

        if (verbosity>1) _printf_("   Interp Cavities: correction factor...\n");
        /* A correction factor is applied because the area of the cavity may
          slightly differs when computed with the dead or the born elements
          This is due to slight difference in the node position due to numerical precision used in ISSM*/
        correction_factor=integrated_area/new_integrated_area;
        if (correction_factor<=0.)  _error_("correction_factor should not be larger than 0. Here it is: " << correction_factor<<" ");

        if(method==0)
        {
            if (verbosity>1) _printf_("   Interp Cavities: method1...\n");
            /* First method (fully conservative but too diffusive)
              Compute the weighted average over the dead element of the cavity
              and give this value to the new elements of the cavity    */

            if (verbosity>1) _printf_("   Interp Cavities: integrated_variables..\n");


            for(int i=0; i<nb_dead_elements; i++)
            {
                tmp_element=dead_elements[i];

                for(int j=0; j<nb_variables; j++)
                    integrated_variables[j]+=surface_old[(int)tmp_element-1]*IntMatrix_in[((int)tmp_element-1)*nb_variables+j];
            }

            if (verbosity>1) _printf_("   Interp Cavities: InrMatrix_out..\n");
            for(int i=0; i<nb_born_elements; i++)
            {
                tmp_element=born_elements[i];

                for(int j=0; j<nb_variables; j++)
                    IntMatrix_out[((int)tmp_element-1)*nb_variables+j]=integrated_variables[j]/integrated_area*correction_factor;
            }
        }
        else
        {
            /* Second method (fully conservative)*/

            /* interpolation over the cavity*/
            InterpCavity(tmp_mean_variables, tmp_integrated_area,
                  nb_dead_elements, nb_born_elements, nb_variables,
                  dead_elements, born_elements, gate.PreviousNumbering,
                  IntMatrix_in, method_in, bamgmesh_old, bamgmesh_new, 0);

            /* Sanity check
              does the integrated area corresponds to the area of the born element*/
            for(int i=0; i<nb_born_elements; i++)
            {
                area_error=fabs(tmp_integrated_area[i]-surface_new[(int)born_elements[i]-1]);
                if(area_error>area_error_tolerance)
                {
                    _printf_("tmp_integrated_area[i]=" << tmp_integrated_area[i] << "\n");
                    _printf_("(int)born_elements[i]-1=" << (int)born_elements[i]-1 << "\n");
                    _printf_("surface_new[born_elements[i]-1]=" << surface_new[(int)born_elements[i]-1] << "\n");
                    _printf_("area_error=" << area_error << "\n");
                    _printf_("area_error_tolerance=" << area_error_tolerance << "\n");
                    _error_("Conservation issue!!! The area of the new element does not match with the sum of the triangle intersections...");
                }
            }

            for(int i=0; i<nb_born_elements; i++)
            {
                tmp_element=born_elements[i];

                for(int j=0; j<nb_variables; j++)
                    IntMatrix_out[((int)tmp_element-1)*nb_variables+j]=tmp_mean_variables[i*nb_variables+j]*correction_factor;
            }
        }
    }

    if (verbosity>1) _printf_("   Interp Cavities: end!\n");



    // destroy temporary arrays
    xDelete<double>(dead_elements);
    xDelete<double>(born_elements);

    xDelete<double>(integrated_variables);
    xDelete<double>(tmp_mean_variables);
    xDelete<double>(tmp_integrated_area);

    xDelete<double>(gate.PreviousNumbering);
    xDelete<double>(gate.old_elements);
    xDelete<double>(gate.new_elements);
    xDelete<int>(gate.size_born_cavity);
    xDelete<int>(gate.size_dead_cavity);
    xDelete<double>(gate.element_born_cavity);
    xDelete<double>(gate.element_dead_cavity);

    /*---- end of the interpolation cavity by cavity      ------*/
    *pdata_interp=IntMatrix_out;

    return 1;
}

int DetectCavities(InterpFromMeshToMesh2dCavitiesThreadStruct* gate, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new){

    int verbosity;

    /*Bamg options*/
    verbosity=0;//bamgopts->verbose;

    if (verbosity>1) _printf_("   Detect_cavities: Read input...\n");
    /*---------- Input  ----------*/

    double *bamg_mesh_Triangles                 = bamgmesh_old->Triangles;
    double *bamg_mesh_NodalElementConnectivity  = bamgmesh_old->NodalElementConnectivity;
    double *bamg_mesh_NodalConnectivity         = bamgmesh_old->NodalConnectivity ;
    double *bamg_mesh_ElementConnectivity       = bamgmesh_old->ElementConnectivity;

    int bamg_mesh_Nn                         = bamgmesh_old->VerticesSize[0];
    int bamg_mesh_NNodalConnectivity         = bamgmesh_old->NodalConnectivitySize[1];
    int bamg_mesh_NNodalElementConnectivity  = bamgmesh_old->NodalElementConnectivitySize[1];
    int bamg_mesh_Ne                         = bamgmesh_old->TrianglesSize[0];
    int bamg_mesh_NVerticesOnGeomVertex      = bamgmesh_old->VerticesOnGeomVertexSize[0];

    double *new_bamg_mesh_Triangles             = bamgmesh_new->Triangles;
    double *new_bamg_mesh_ElementConnectivity   = bamgmesh_new->ElementConnectivity;
    double *new_bamg_mesh_PreviousNumbering     = bamgmesh_new->PreviousNumbering;

    int new_bamg_mesh_Nn                     = bamgmesh_new->VerticesSize[0];
    int new_bamg_mesh_NVerticesOnGeomVertex  = bamgmesh_new->VerticesOnGeomVertexSize[0];
    int new_bamg_mesh_Ne                     = bamgmesh_new->TrianglesSize[0];


    if (verbosity>1) _printf_("   bamg_mesh_Nn = " << bamg_mesh_Nn << "\n");
    if (verbosity>1) _printf_("   bamg_mesh_NNodalConnectivity = " << bamg_mesh_NNodalConnectivity << "\n");
    if (verbosity>1) _printf_("   bamg_mesh_NNodalElementConnectivity = " << bamg_mesh_NNodalElementConnectivity << "\n");
    if (verbosity>1) _printf_("   bamg_mesh_Ne = " << bamg_mesh_Ne << "\n");
    if (verbosity>1) _printf_("   bamg_mesh_NVerticesOnGeomVertex = " << bamg_mesh_NVerticesOnGeomVertex << "\n");

    if (verbosity>1) _printf_("   new_bamg_mesh_Nn = " << new_bamg_mesh_Nn << "\n");
    if (verbosity>1) _printf_("   new_bamg_mesh_NVerticesOnGeomVertex = " << new_bamg_mesh_NVerticesOnGeomVertex << "\n");
    if (verbosity>1) _printf_("   new_bamg_mesh_Ne = " << new_bamg_mesh_Ne << "\n");

    if (verbosity>1) _printf_("   Detect_cavities: Initialize arrays...\n");

    if(bamg_mesh_NVerticesOnGeomVertex!=new_bamg_mesh_NVerticesOnGeomVertex)
        _error_("not the same number of vertices on the Geom in the old and new mesh");


    /*---------- Output arrays      ----------*/
    double  *PreviousNumbering = xNew<double>(new_bamg_mesh_Nn);

    /*---------- Intermediary arrays---------- */
    double *tmp_old_elements = xNew<double>(bamg_mesh_Ne);
    double *tmp_new_elements = xNew<double>(new_bamg_mesh_Ne);
    double  *cavity_number_dead = xNew<double>(bamg_mesh_Ne);
    int  *cavity_number_born = xNew<int>(new_bamg_mesh_Ne);
    int  *is_born_element_to_treat = xNew<int>(new_bamg_mesh_Ne); /* boolean to know if the born element is to be treated */
    int  *is_dead_element_to_treat = xNew<int>(bamg_mesh_Ne);   /* boolean to know if the dead element is to be treated */

    if (verbosity>1) _printf_("   Detect_cavities: Initialize local variables...\n");
    /*---------- Local variables    ----------*/
    int l, l0, l1, l2, j_next;

    int nb_matching_elements, nb_cavities;
    int max_size_born_cavity, max_size_dead_cavity;

    int edge_node_1, edge_node_2;
    int Previous_edge_node_1, Previous_edge_node_2;
    int tmp_node_1, tmp_node_2;
    double Previous_connected_element;
    int tmp_born, tmp_dead;
    double Previous_num_node[3], Previous_connected_element_num_node[3];
    int max_candidate_element_to_treat=10000, number_candidate_element_to_treat, number_candidate_element_treated;

    int  *candidate_element_to_treat = xNew<int>(max_candidate_element_to_treat);   /* boolean to know if the dead element is to be treated */

    int candidate_dead_element[2], nb_candidate_dead_element;
    int nb_dead, nb_born, nb_connected_nodes;
    int neighbour;
    double neighbour_double;
    int new_cavity_number;
    int i_dead_element_node_1;
    int i_dead_element_node_2;
    double dead_element_node_1, dead_element_node_2;
    int i_dead_element, max_cavity_size;

    int found_intersect;

    if (verbosity>1) _printf_("   Detect_cavities: Begin code...\n");
    /* ---------------------------------------------------------------
     * BEGIN CODE
     * --------------------------------------------------------------- */

    /*  The previous numbering of the nodes on the boundaries given by bamg is wrong */
    /*  Bamg adds these nodes at the beginning of the list of nodes. */
    /*  In our case, these nodes are always the same as the boundary is not */
    /*  moving nor adapted. */
    for (int i=0; i<new_bamg_mesh_NVerticesOnGeomVertex; i++)
    {
        if (verbosity>1) _printf_(" i " << i << "(ulong) new_bamg_mesh_PreviousNumbering[i] " << new_bamg_mesh_PreviousNumbering[i] << "\n");
        PreviousNumbering[i]=double(i)+1;
    }
    if (verbosity>1) _printf_("   Detect_cavities: loop1...\n");
    for (int i=new_bamg_mesh_NVerticesOnGeomVertex; i<new_bamg_mesh_Nn; i++)
    {
        if (verbosity>1) _printf_(" i " << i << "(ulong) new_bamg_mesh_PreviousNumbering[i] " << new_bamg_mesh_PreviousNumbering[i] << "\n");
        PreviousNumbering[i]=new_bamg_mesh_PreviousNumbering[i];
    }


    if (verbosity>1) _printf_("   Detect_cavities: find the correspondance...\n");
    /*  To find the correspondance between the new and old elements */
    nb_matching_elements=0;

    for (int i=0; i<new_bamg_mesh_Ne; i++)
    {
        for (int j=0; j<3; j++)
        {
            Previous_num_node[j]=PreviousNumbering[(int)new_bamg_mesh_Triangles[4*i+j]-1];
        }
        found_intersect=0;
        if((Previous_num_node[0]>0) && (Previous_num_node[1]>0) && (Previous_num_node[2]>0))
        {
            // if (verbosity>1) _printf_("   Detect_cavities: found intersect...\n");
            for (int j=0; (j<3) && (found_intersect==0); j++)
            {
                for (int k=0; (k<bamg_mesh_NNodalElementConnectivity)  && (found_intersect==0); k++)
                {
                    // if (verbosity>1) _printf_("   Detect_cavities: loop k...\n");
                    Previous_connected_element=bamg_mesh_NodalElementConnectivity[((int)Previous_num_node[j]-1)*bamg_mesh_NNodalElementConnectivity+k];

                    /* If NaN in matlab (0 after calling uint32)  */
                    if (std::isnan(Previous_connected_element))
                        continue;

                    for (l=0; l<3; l++)
                    {
                        Previous_connected_element_num_node[l]=bamg_mesh_Triangles[((int)Previous_connected_element-1)*4+l];
                    }

                    for (l=0; (l<3)  && (found_intersect==0); l++)
                    {
                        /* l0, l1, l2 */
                        l0=l;
                        l1=next_node(l0);
                        l2=next_node(l1);

                        if((Previous_connected_element_num_node[0]==Previous_num_node[l0])
                                && (Previous_connected_element_num_node[1]==Previous_num_node[l1])
                                && (Previous_connected_element_num_node[2]==Previous_num_node[l2]))
                        {
                            tmp_old_elements[nb_matching_elements]=Previous_connected_element;
                            tmp_new_elements[nb_matching_elements]=double(i)+1;
                            nb_matching_elements=nb_matching_elements+1;
                            found_intersect=1;
                        }
                    }
                }
            }
        }
    }

    if (verbosity>1) _printf_("   Detect_cavities: reducing the old_elements and new_elements arrays...\n");
    /* reducing the old_elements and new_elements arrays */
    /*---------- Output arrays      ---------- */
    double *old_elements = xNew<double>(nb_matching_elements);
    double *new_elements = xNew<double>(nb_matching_elements);

    /*---------- Intermediary arrays---------- */
    std::vector<double> is_not_dead(bamg_mesh_Ne,0.);
    std::vector<double> is_not_born(new_bamg_mesh_Ne,0.);

    std::vector<double> ind_dead(bamg_mesh_Ne,0.);
    std::vector<double> ind_born(new_bamg_mesh_Ne,0.);

    //double *boundary_between_two_cavities = xNew<double>(bamg_mesh_Ne*3);
    std::vector<bool> boundary_between_two_cavities(bamg_mesh_Ne*3,false);

    if (verbosity>1) _printf_("   nb_matching_elements = " << nb_matching_elements << "\n");

    if (verbosity>1) _printf_("   Detect_cavities: loop2...\n");
    for (int i=0; i<nb_matching_elements; i++)
    {
        old_elements[i]=tmp_old_elements[i];
        new_elements[i]=tmp_new_elements[i];

        is_not_dead[(int)tmp_old_elements[i]-1]=1.;
        is_not_born[(int)tmp_new_elements[i]-1]=1.;
    }

    if (verbosity>1) _printf_("   Detect_cavities: destroy...\n");
    /* destroy the temporary matrix */
    xDelete<double>(tmp_old_elements);
    xDelete<double>(tmp_new_elements);

    if (verbosity>1) _printf_("   Detect_cavities: loop3...\n");
    /*  The old elements without correspondance are called "dead element" */
    /*  The new elements without correspondance are called "born element" */
    nb_dead=0;
    for (int i=0; i<bamg_mesh_Ne; i++)
    {
        cavity_number_dead[i]=0.;
        if(is_not_dead[i]!=1.)
        {
            ind_dead[nb_dead]=double(i)+1;
            cavity_number_dead[i]=-1.;
            nb_dead++;
        }
    }

    if (verbosity>1) _printf_("   Detect_cavities: loop4...\n");
    nb_born=0;
    for (int i=0; i<new_bamg_mesh_Ne; i++)
    {
        cavity_number_born[i]=0;
        if(is_not_born[i]!=1.)
        {
            ind_born[nb_born]=double(i)+1;
            cavity_number_born[i]=-1;
            nb_born++;
        }
    }
    if (verbosity>1) _printf_("   nb_dead = " << nb_dead << "\n");
    if (verbosity>1) _printf_("   nb_born = " << nb_born << "\n");

    if (verbosity>1) _printf_("   Detect_cavities: destroy2...\n");

    if (verbosity>1) _printf_("   Detect_cavities: detect the cavities...\n");
    /*  Loop that detect the cavities */
    new_cavity_number=0;
    for (int i_born=0; i_born<nb_born; i_born++)
    {
        /*  first element to be checked */
        tmp_born=(int)ind_born[i_born];

        /*  if the born element is already treated, do nothing */
        if(cavity_number_born[tmp_born-1]!=-1)
            continue;

        /*  We start the detection of a new cavity */
        new_cavity_number++;

        /*  we will store the candidate element to treat in an expandable vector, whose */
        /*  initial size is max_candidate_element_to_treat */
        number_candidate_element_to_treat=0;
        for (int i=0; i<max_candidate_element_to_treat; i++)
            candidate_element_to_treat[i]=0;

        for (int i=0; i<new_bamg_mesh_Ne; i++)
            is_born_element_to_treat[i]=0;

        /*  Initialization with this first element */
        candidate_element_to_treat[number_candidate_element_to_treat]=tmp_born;
        number_candidate_element_to_treat++;

        is_born_element_to_treat[tmp_born-1]=1;

        /*  Starting the treating loop */
        number_candidate_element_treated=0;
        while(number_candidate_element_treated<number_candidate_element_to_treat)
        {
            tmp_born=candidate_element_to_treat[number_candidate_element_treated];
            number_candidate_element_treated++;

            cavity_number_born[tmp_born-1]=new_cavity_number;

            /*  Loop over the neighbours, */
            for (int j_born=0; j_born<3; j_born++)
            {
                neighbour_double=new_bamg_mesh_ElementConnectivity[(tmp_born-1)*3+j_born];
                neighbour=(int)neighbour_double;

                /*  if neighbour exists and is already detected in this cavity, do nothing */
                if (!std::isnan(neighbour_double))
                {
                    if(cavity_number_born[neighbour-1]==new_cavity_number)
                        continue;
                }

                /*  if the edge exists in the old mesh */
                /*  the edge j goes from node j to node jnext */
                edge_node_1 = next_node(j_born);
                edge_node_2 =next_node(edge_node_1);

                Previous_edge_node_1=PreviousNumbering[(int)new_bamg_mesh_Triangles[(tmp_born-1)*4+edge_node_1]-1];
                Previous_edge_node_2=PreviousNumbering[(int)new_bamg_mesh_Triangles[(tmp_born-1)*4+edge_node_2]-1];

                // _printf_("   Previous_edge_node_1 = " << Previous_edge_node_1 << "\n");
                // _printf_("   Previous_edge_node_2 = " << Previous_edge_node_2 << "\n");

                if((Previous_edge_node_1>0) && (Previous_edge_node_2>0))
                {
                    /*  the last indices in the connectivity vector is the number */
                    /*  of connected nodes */

                    nb_connected_nodes=(int)bamg_mesh_NodalConnectivity[(Previous_edge_node_1-1)*bamg_mesh_NNodalConnectivity+bamg_mesh_NNodalConnectivity-1];

                    found_intersect=0;
                    double _edge_node_2 = (double)Previous_edge_node_2;
                    for (int i=0; i<nb_connected_nodes; i++)
                        if(bamg_mesh_NodalConnectivity[(Previous_edge_node_1-1)*bamg_mesh_NNodalConnectivity+i]==_edge_node_2)
                        {
                            found_intersect=1;
                            break;
                        }

                    /*  if this edge exists in the old mesh, we are at the */
                    /*  boundary of the cavity. */
                    if(found_intersect==1)
                    {
                        /*  we look for the dead element that share the same */
                        /*  edge, to include it in the cavity */

                        /* resetting candidate_dead_elements */
                        candidate_dead_element[0]=0;
                        candidate_dead_element[1]=0;
                        nb_candidate_dead_element=0;

                        /* searching for candidates */
                        for (i_dead_element_node_1=0; i_dead_element_node_1<bamg_mesh_NNodalElementConnectivity; i_dead_element_node_1++)
                        {
                            dead_element_node_1=bamg_mesh_NodalElementConnectivity[(Previous_edge_node_1-1)*bamg_mesh_NNodalElementConnectivity+i_dead_element_node_1];

                            if(dead_element_node_1==0)
                                break;

                            for (i_dead_element_node_2=0; i_dead_element_node_2<bamg_mesh_NNodalElementConnectivity; i_dead_element_node_2++)
                            {
                                dead_element_node_2=bamg_mesh_NodalElementConnectivity[(Previous_edge_node_2-1)*bamg_mesh_NNodalElementConnectivity+i_dead_element_node_2];

                                if(dead_element_node_2==0.)
                                    break;

                                if(dead_element_node_1==dead_element_node_2)
                                {
                                    candidate_dead_element[nb_candidate_dead_element]=(int)dead_element_node_1;
                                    nb_candidate_dead_element++;
                                    // _printf_("   nb_candidate_dead_element = " << nb_candidate_dead_element << "\n");
                                    // _printf_("   dead_element_node_1 = " << dead_element_node_1 << "\n");
                                    // _printf_("   dead_element_node_2 = " << dead_element_node_2 << "\n");
                                }
                            }
                        }

                        if(nb_candidate_dead_element>2)
                            _error_("nb_candidate_dead_element should not be larger than 2. Here it is: " << nb_candidate_dead_element <<" ");

                        /*  register if the edge is the boundary between two cavities */
                        for (i_dead_element=0; i_dead_element<nb_candidate_dead_element; i_dead_element++)
                        {
                            tmp_dead=candidate_dead_element[i_dead_element];

                            int j;
                            for (j=0; j<3; j++)
                            {
                                /*  if the edge exists in the old mesh */
                                /*  the edge j goes from node j to node jnext */
                                edge_node_1=next_node(j);
                                edge_node_2=next_node(edge_node_1);

                                tmp_node_1=(int)bamg_mesh_Triangles[(tmp_dead-1)*4+edge_node_1];
                                tmp_node_2=(int)bamg_mesh_Triangles[(tmp_dead-1)*4+edge_node_2];

                                if( (tmp_node_1==Previous_edge_node_1) && (tmp_node_2==Previous_edge_node_2))
                                    break;
                            }

                            if(j<3)
                            {
                                if(cavity_number_dead[tmp_dead-1]==-1.)
                                    cavity_number_dead[tmp_dead-1]=(double)new_cavity_number;

                                boundary_between_two_cavities[(tmp_dead-1)*3+j]=true;
                                break;
                            }
                        }
                        continue;
                    } /* end of the if on found_intersect */
                } /* end of the if */

                /*  if no neighbour, do nothing */
                if (std::isnan(neighbour_double))
                    continue;

                /*  if neighbour is not born, do nothing */
                if(cavity_number_born[neighbour-1]==0.)
                    continue;

                /*  if neighbour is already selected to be treated */
                if(is_born_element_to_treat[neighbour-1]==1)
                    continue;

                if(number_candidate_element_to_treat==max_candidate_element_to_treat)
                    _error_("You have to increase the value of max_candidate_element_to_treat in detect_cavities_mex.c. Sorry for that...");

                candidate_element_to_treat[number_candidate_element_to_treat]=neighbour;
                number_candidate_element_to_treat++;
                is_born_element_to_treat[neighbour-1]=1;

            } /*  end of the loop over the neighbours */
        } /*  end of the while loop */
    } /*  end of the loop of the born elements */


    nb_cavities= new_cavity_number;

    if (verbosity>1) _printf_("   Detect_cavities: define their cavity number...\n");
    /*  New loop on the dead element to define their cavity number. */
    for (int i_dead=0; i_dead<nb_dead; i_dead++)
    {
        /*  first element to be checked */
        tmp_dead=(int)ind_dead[i_dead];

        /*  We want to start from a dead element that has already be assigned to */
        /*  a cavity. */

        /*  if the dead element is not already treated, do nothing. It will be */
        /*  treated from another element */
        if(cavity_number_dead[tmp_dead-1]==-1.)
            continue;

        new_cavity_number=(int)cavity_number_dead[tmp_dead-1];

        /*  we will store the candidate element to treat in an expandable vector, whose */
        /*  initial size is max_candidate_element_to_treat */
        number_candidate_element_to_treat=0;
        for (int i=0; i<max_candidate_element_to_treat; i++)
            candidate_element_to_treat[i]=0;

        for (int i=0; i<bamg_mesh_Ne; i++)
            is_dead_element_to_treat[i]=0;

        /*  Initialization with this first element */
        candidate_element_to_treat[number_candidate_element_to_treat]=tmp_dead;
        number_candidate_element_to_treat++;
        is_dead_element_to_treat[tmp_dead-1]=1;

        /*  Starting the treating loop */
        number_candidate_element_treated=0;

        while(number_candidate_element_treated<number_candidate_element_to_treat)
        {
            tmp_dead=candidate_element_to_treat[number_candidate_element_treated];
            number_candidate_element_treated++;

            cavity_number_dead[tmp_dead-1]=(double)new_cavity_number;

            for (int j_dead=0; j_dead<3; j_dead++)
            {
                neighbour_double = bamg_mesh_ElementConnectivity[(tmp_dead-1)*3+j_dead];
                neighbour=(int)neighbour_double;

                /*  if no neighbour, do nothing */
                if (std::isnan(neighbour_double))
                    continue;

                /*  if neighbour is already treated, do nothing */
                if(cavity_number_dead[neighbour-1]!=-1.)
                    continue;

                /*  if the edge corresponds to the boudary between two cavities */
                if( boundary_between_two_cavities[(tmp_dead-1)*3+j_dead] )
                    continue;

                /*  if neighbour is already selected to be treated */
                if(is_dead_element_to_treat[neighbour-1]==1)
                    continue;

                if(number_candidate_element_to_treat==max_candidate_element_to_treat)
                    _error_("You have to increase the value of max_candidate_element_to_treat in detect_cavities_mex.c. Sorry for that...");

                candidate_element_to_treat[number_candidate_element_to_treat]=neighbour;
                number_candidate_element_to_treat++;
                is_dead_element_to_treat[neighbour-1]=1;
            }
        } /*  end of the while loop */
    } /*  end of the loop of the dead elements */

    int *size_born_cavity = xNew<int>(nb_cavities);
    int *size_dead_cavity = xNew<int>(nb_cavities);

    for (int i_born=0; i_born<nb_cavities; i_born++)
        size_born_cavity[i_born]=0;

    for (int i_dead=0; i_dead<nb_cavities; i_dead++)
        size_dead_cavity[i_dead]=0;


    if (verbosity>1) _printf_("   Detect_cavities: loop6...\n");
    for (int i_born=0; i_born<nb_born; i_born++)
    {
        new_cavity_number=cavity_number_born[(int)ind_born[i_born]-1];

        if(new_cavity_number<0)
            _error_("A born element has not been assigned to a cavity");

        if(new_cavity_number>0)
            size_born_cavity[new_cavity_number-1]++;
    }

    bool not_all_dead_assigned=0;
    if (verbosity>1) _printf_("   Detect_cavities: loop7...\n");
    for (int i_dead=0; i_dead<nb_dead; i_dead++)
    {
        new_cavity_number=cavity_number_dead[(int)ind_dead[i_dead]-1];
        if (verbosity>1) _printf_("i_dead:" << i_dead << ", ind_dead[i_dead]:" << ind_dead[i_dead] << ", new_cavity_number: " << new_cavity_number <<"\n");

        if(new_cavity_number<0)
        {
            if (verbosity>1) _printf_("not assigned \n");
            not_all_dead_assigned=1;
        }

        if(new_cavity_number>0)
            size_dead_cavity[new_cavity_number-1]++;
    }
    if(not_all_dead_assigned)
        _error_("A dead element has not been assigned to a cavity");

    if (verbosity>1) _printf_("   Detect_cavities: loop8...\n");
    max_size_born_cavity=0;
    max_size_dead_cavity=0;
    for (int i=0; i<nb_cavities; i++)
    {
        if(max_size_born_cavity<size_born_cavity[i])
            max_size_born_cavity= size_born_cavity[i];

        size_born_cavity[i]=0;

        if(max_size_dead_cavity<size_dead_cavity[i])
            max_size_dead_cavity= size_dead_cavity[i];

        size_dead_cavity[i]=0;
    }

    if (verbosity>1) _printf_("  max_size_born_cavity" << max_size_born_cavity << "\n");
    if (verbosity>1) _printf_("  max_size_dead_cavity" << max_size_dead_cavity << "\n");
    if (verbosity>1) _printf_("  nb_cavities*max_size_born_cavity" << nb_cavities*max_size_born_cavity << "\n");
    if (verbosity>1) _printf_("  nb_cavities*max_size_dead_cavity" << nb_cavities*max_size_dead_cavity << "\n");

    double *element_born_cavity = xNew<double>(nb_cavities*max_size_born_cavity);
    double *element_dead_cavity = xNew<double>(nb_cavities*max_size_dead_cavity);

    if (verbosity>1) _printf_("   Detect_cavities: loop9...\n");
    for (int i_born=0; i_born<nb_born; i_born++)
    {
        new_cavity_number=cavity_number_born[(int)ind_born[i_born]-1];

        if(new_cavity_number>0)
        {
            element_born_cavity[(new_cavity_number-1)*max_size_born_cavity+size_born_cavity[new_cavity_number-1]]=ind_born[i_born];
            //element_born_cavity[nb_cavities*size_born_cavity[new_cavity_number-1]+new_cavity_number-1]=ind_born[i_born];
            size_born_cavity[new_cavity_number-1]++;
        }
    }

    if (verbosity>1) _printf_("   Detect_cavities: loop10...\n");
    for (int i_dead=0; i_dead<nb_dead; i_dead++)
    {
        new_cavity_number=cavity_number_dead[(int)ind_dead[i_dead]-1];

        if(new_cavity_number>0)
        {
            element_dead_cavity[(new_cavity_number-1)*max_size_dead_cavity+size_dead_cavity[new_cavity_number-1]]=ind_dead[i_dead];
            //element_dead_cavity[nb_cavities*size_dead_cavity[new_cavity_number-1]+new_cavity_number-1]=ind_dead[i_dead];
            size_dead_cavity[new_cavity_number-1]++;
        }
    }

    if (verbosity>1) _printf_("   Detect_cavities: destroy final...\n");
    /* destroy the temporary matrix */
    // xDelete<double>(boundary_between_two_cavities);
    xDelete<double>(cavity_number_dead);

    xDelete<int>(cavity_number_born);
    xDelete<int>(is_born_element_to_treat);
    xDelete<int>(is_dead_element_to_treat);
    xDelete<int>(candidate_element_to_treat);

    // destroy temporary arrays



    if (verbosity>1) _printf_("   Detect_cavities: output...\n");
    /*Output*/
    gate->PreviousNumbering     =PreviousNumbering;
    gate->old_elements          =old_elements;
    gate->new_elements          =new_elements;
    gate->size_born_cavity      =size_born_cavity;
    gate->size_dead_cavity      =size_dead_cavity;
    gate->element_born_cavity   =element_born_cavity;
    gate->element_dead_cavity   =element_dead_cavity;
    gate->new_bamg_mesh_Nn      =new_bamg_mesh_Nn;
    gate->nb_matching_elements  =nb_matching_elements;
    gate->nb_cavities           =nb_cavities;
    gate->max_size_born_cavity  =max_size_born_cavity;
    gate->max_size_dead_cavity  =max_size_dead_cavity;

    return 1;
}

int InterpCavity(double* tmp_mean_variables, double* tmp_integrated_area,
            int nb_dead_elements, int nb_born_elements, int nb_variables,
            double* dead_elements, double* born_elements, double* PreviousNumbering,
            double* IntMatrix_in, int* method_in, BamgMesh* bamgmesh_old, BamgMesh* bamgmesh_new, int debug_born_elements_i){

    /*---------- Input  ----------*/
    int *dead_num_node      = xNew<int>(3*nb_dead_elements);
    int *born_num_node      = xNew<int>(3*nb_born_elements);

    double *x_dead_elements = xNew<double>(3*nb_dead_elements);
    double *x_born_elements = xNew<double>(3*nb_born_elements);
    double *y_dead_elements = xNew<double>(3*nb_dead_elements);
    double *y_born_elements = xNew<double>(3*nb_born_elements);

    double *dead_variables  = xNew<double>(nb_variables*nb_dead_elements);

    /*---------- Local variables ----------*/
    int i_born, i_dead, j_born, j_dead, jnext, jnext_dead, j_edge_polygon, k;

    int tmp_sign_cross, max_nb_intersections, sign_polygon;

    double *integrated_variables = xNew<double>(nb_variables);
    double tmp_distance, min_distance;

    int born_node_on_dead_node[3], born_node_in_dead_element[3], dead_node_in_born_element[3], born_node_on_dead_edge[3], dead_node_on_born_edge[3];
    double born_edge_intersect_dead_edge[3*3];
    double alpha_tmp[3];
    double alpha_tol=1e-12;

    double test_integrated_area, area_polygon, area_intersection;

    double dist_max, dist, tmp_point_x, tmp_point_y;
    int ind_max;

    double x_vect_born_dead_node_j, y_vect_born_dead_node_j, x_vect_dead_node_j_jnext, y_vect_dead_node_j_jnext;
    double x_vect_dead_born_node_j, y_vect_dead_born_node_j, x_vect_born_node_j_jnext, y_vect_born_node_j_jnext;
    double cross_prod;

    double P0x, P0y, P1x, P1y;
    double Q0x, Q0y, Q1x, Q1y;
    double Xx, Xy;

    double alpha_denominator, alpha_numerator, alpha, beta;

    double intersection_points_x[18], intersection_points_y[18];
    int nb_intersection_points, i_intersection_points;

    double polygon_points_x[6], polygon_points_y[6];
    int next_polygon_point[6];
    int nb_polygon_points, i_polygon_points, inext_polygon_points;

    double x_target, y_target, x_vect_target, y_vect_target, x_edge_polygon, y_edge_polygon;
    int flag_print;
    div_t divresult;
    int tmp_ind, tmp_element;

    /* ---------------------------------------------------------------
     * BEGIN CODE
     * --------------------------------------------------------------- */
    for (i_born=0; i_born<nb_born_elements; i_born++)
    {
        tmp_element=born_elements[i_born];
        for (j_born=0; j_born<3; j_born++)
        {
            tmp_ind=bamgmesh_new->Triangles[4*(tmp_element-1)+j_born];
            born_num_node[i_born*3+j_born]=tmp_ind;
            x_born_elements[i_born*3+j_born]=bamgmesh_new->Vertices[3*(tmp_ind-1)];
            y_born_elements[i_born*3+j_born]=bamgmesh_new->Vertices[3*(tmp_ind-1)+1];
        }
    }

    for (i_dead=0; i_dead<nb_dead_elements; i_dead++)
    {
        tmp_element=dead_elements[i_dead];
        for(int k=0; k<nb_variables; k++)
            dead_variables[i_dead*nb_variables+k]=IntMatrix_in[(tmp_element-1)*nb_variables+k];

        for (j_dead=0; j_dead<3; j_dead++)
        {
            tmp_ind=bamgmesh_old->Triangles[4*(tmp_element-1)+j_dead];
            dead_num_node[i_dead*3+j_dead]=tmp_ind;
            x_dead_elements[i_dead*3+j_dead]=bamgmesh_old->Vertices[3*(tmp_ind-1)];
            y_dead_elements[i_dead*3+j_dead]=bamgmesh_old->Vertices[3*(tmp_ind-1)+1];
        }

    }

    /* Loop over the born elements of the cavity*/
    for (i_born=0; i_born<nb_born_elements; i_born++)
    {
        flag_print=0;
        if(born_elements[i_born]==debug_born_elements_i)
            flag_print=1;

        /*-------- Initializing the integrated quantities --------- -*/
        for (k=0; k<nb_variables; k++)
            integrated_variables[k]=0.;

        min_distance=-1.;

        test_integrated_area=0.;

        /*-------- Lets starting the loop on the dead elements of the cavity ---------- */
        /* Loop over the dead elements of the cavity */
        for (i_dead=0; i_dead<nb_dead_elements; i_dead++)
        {
            /*-------- Initialization ---------- */


            /* 1) Shared nodes: */
            for (j_born=0; j_born<3; j_born++)
                born_node_on_dead_node[j_born]=-1;          /* -1 if node is not on a node of dead element, local indice of the node of the dead element otherwise */

            /* 2) Born nodes included in dead element: */
            for (j_born=0; j_born<3; j_born++)
                born_node_in_dead_element[j_born]=0;       /* 0 if node from born element is not in dead element, 1 otherwise */

            /* 3) Dead nodes included in born element: */
            for (j_dead=0; j_dead<3; j_dead++)
                dead_node_in_born_element[j_dead]=0;       /* 0 if node from dead element is not in born element, 1 otherwise */

            /* 4) Edge to edge intersections: */
            /* Each line   corresponds to a born edge */
            /* Each column corresponds to a dead edge */
            /* 0<=alpha<=1 defines the position of the intersection relatively to born edge */
            for (j_born=0; j_born<3; j_born++)
                for (j_dead=0; j_dead<3; j_dead++)
                    born_edge_intersect_dead_edge[j_born*3+j_dead]=-1.;   /* -1 if edge i does no intersect edge j, alpha otherwise */

            /* 4 bis) Born nodes on dead edge: */
            for (j_born=0; j_born<3; j_born++)
                born_node_on_dead_edge[j_born]=0;       /* 0 if node from born element is not in dead element, 1 otherwise */

            /* 4 tris) Born nodes on dead edge: */
            for (j_born=0; j_born<3; j_born++)
                dead_node_on_born_edge[j_born]=0;       /* 0 if node from born element is not in dead element, 1 otherwise */

            /* 4 tris) Born nodes on dead edge: */
            dist_max=0.;
            for (j_born=0; j_born<3; j_born++)
            {
                if(j_born<2)
                    dist=pow(pow(x_born_elements[i_born*3+j_born]-x_born_elements[i_born*3+j_born+1],2.)+pow(y_born_elements[i_born*3+j_born]-y_born_elements[i_born*3+j_born+1],2.),0.5);
                else
                    dist=pow(pow(x_born_elements[i_born*3+j_born]-x_born_elements[i_born*3+0  ],2.)+pow(y_born_elements[i_born*3+j_born]-y_born_elements[i_born*3+0  ],2.),0.5);

                if(dist>dist_max)
                    dist_max=dist;
            }

            /*-------- Detection ---------- */

            /* Loop over the born nodes */
            for (j_born=0; j_born<3; j_born++)
            {

                /* 1) Shared nodes: */

                /* Loop over the dead nodes to detect if the node correspond to a node of the dead element */
                for (j_dead=0; j_dead<3; j_dead++)
                    if(PreviousNumbering[born_num_node[i_born*3+j_born]-1]==dead_num_node[i_dead*3+j_dead])
                        born_node_on_dead_node[j_born]=j_dead;
                    else /* Bamg sometimes says that the nodes are different, whereas they are almost at the same position */
                    {

                        dist=pow(pow(x_dead_elements[i_dead*3+j_dead]-x_born_elements[i_born*3+j_born],2.)+pow(y_dead_elements[i_dead*3+j_dead]-y_born_elements[i_born*3+j_born],2.),0.5);
                        if(flag_print==1)
                            _printf_("dist/dist_max=" << dist/dist_max <<"\n");

                        if( dist/dist_max<alpha_tol )
                        {
                            born_node_on_dead_node[j_born]=j_dead;
                            if(flag_print==1)
                                _printf_("set to same node as dist/dist_max<alpha_tol smaller than =" << alpha_tol <<"\n");
                        }
                    }


                /* 2) Born nodes included in dead element: */

                /* detect if the node is inside the dead element */
                /* only if the node is not shared */
                if(born_node_on_dead_node[j_born]==-1)
                {
                    tmp_sign_cross=0;
                    /* Loop over the dead nodes */
                    for (j_dead=0; j_dead<3; j_dead++)
                    {
                        /* vector from dead node j_dead to born node */
                        x_vect_born_dead_node_j=x_born_elements[i_born*3+j_born]-x_dead_elements[i_dead*3+j_dead];
                        y_vect_born_dead_node_j=y_born_elements[i_born*3+j_born]-y_dead_elements[i_dead*3+j_dead];

                        /* vector from dead node j_dead to j_dead+1 */
                        divresult = div(j_dead+1,3);
                        jnext_dead=divresult.rem;

                        x_vect_dead_node_j_jnext=x_dead_elements[i_dead*3+jnext_dead]-x_dead_elements[i_dead*3+j_dead];
                        y_vect_dead_node_j_jnext=y_dead_elements[i_dead*3+jnext_dead]-y_dead_elements[i_dead*3+j_dead];

                        /* add the sign of the cross product */
                        cross_prod=x_vect_born_dead_node_j*y_vect_dead_node_j_jnext-x_vect_dead_node_j_jnext*y_vect_born_dead_node_j;
                        tmp_sign_cross+=sign(cross_prod);
                    }
                    if(flag_print==1)
                        _printf_("in dead, i_born=" << i_born << ", j_born=" << j_born << ",tmp_sign_cross=" << tmp_sign_cross <<"\n");

                    /* check if the 3 cross products have the same sign */
                    if((tmp_sign_cross==(3))||(tmp_sign_cross==-(3)))
                    {
                        born_node_in_dead_element[j_born]=1;
                    }
                }
            } /* end of a loop on the born nodes */

            /* 3) Dead nodes included in born element: */

            /* Loop over the dead nodes */
            for (j_dead=0; j_dead<3; j_dead++)
            {
                /* detect if the node is inside the dead element */
                /* only if the node is not shared */
                if(born_node_on_dead_node[0]!=j_dead && born_node_on_dead_node[1]!=j_dead && born_node_on_dead_node[2]!=j_dead)
                {
                    tmp_sign_cross=0;
                    /* Loop over the born nodes */
                    for (j_born=0; j_born<3; j_born++)
                    {
                        /* vector from born node j to born node */
                        x_vect_dead_born_node_j=x_dead_elements[i_dead*3+j_dead]-x_born_elements[i_born*3+j_born];
                        y_vect_dead_born_node_j=y_dead_elements[i_dead*3+j_dead]-y_born_elements[i_born*3+j_born];

                        /* vector from born node j to j+1 */
                        divresult = div(j_born+1,3);
                        jnext=divresult.rem;

                        x_vect_born_node_j_jnext=x_born_elements[i_born*3+jnext]-x_born_elements[i_born*3+j_born];
                        y_vect_born_node_j_jnext=y_born_elements[i_born*3+jnext]-y_born_elements[i_born*3+j_born];

                        /* add the sign of the cross product */
                        cross_prod=x_vect_dead_born_node_j*y_vect_born_node_j_jnext-x_vect_born_node_j_jnext*y_vect_dead_born_node_j;
                        tmp_sign_cross+=sign(cross_prod);
                    }
                    /* check if the 3 cross products have the same sign */
                    if((tmp_sign_cross==3)||(tmp_sign_cross==-3))
                    {
                        dead_node_in_born_element[j_dead]=1;
                    }
                }
            } /* end of a loop on the dead nodes */

            /* 4) Edge to edge intersections: */
            /* Each line   corresponds to a born edge */
            /* Each column corresponds to a dead edge */
            /* 0<=alpha<=1 defines the position of the intersection relatively to born edge */

            /* Loop over the born edges */
            for (j_born=0; j_born<3; j_born++)
            {
                /* the edge j goes from node j to node jnext */
                divresult = div(j_born+1,3);
                jnext=divresult.rem;

                /* Extremities of the born edge */
                P0x=x_born_elements[i_born*3+j_born];
                P0y=y_born_elements[i_born*3+j_born];
                P1x=x_born_elements[i_born*3+jnext];
                P1y=y_born_elements[i_born*3+jnext];

                /* define the maximum number of intersections to search for */
                max_nb_intersections=2;

                if(born_node_in_dead_element[j_born]==1)
                    max_nb_intersections--;

                if(born_node_in_dead_element[jnext]==1)
                    max_nb_intersections--;

                if(born_node_on_dead_node[j_born]>=0)
                    max_nb_intersections--;

                if(born_node_on_dead_node[jnext]>=0)
                    max_nb_intersections--;

                /* the nodes of the born edge are either in the dead element or on one of its nodes (no intersection to compute) */
                if(max_nb_intersections<=0)
                    continue;

                /* Loop over the dead edges */
                for (j_dead=0; j_dead<3; j_dead++)
                    alpha_tmp[j_dead]=-1.;

                /* Loop over the dead edges */
                for (j_dead=0; j_dead<3; j_dead++)
                {
                    /* the edge j_dead goes from node j_dead to node jnext_dead */
                    divresult = div(j_dead+1,3);
                    jnext_dead=divresult.rem;

                    /* ---------- Compute the intersection of the two edges -------- */

                    /* Extremities of the dead edge */
                    Q0x=x_dead_elements[i_dead*3+j_dead];
                    Q0y=y_dead_elements[i_dead*3+j_dead];
                    Q1x=x_dead_elements[i_dead*3+jnext_dead];
                    Q1y=y_dead_elements[i_dead*3+jnext_dead];

                    /* Denuminator and numerator for alpha */
                    alpha_denominator= (P1x-P0x)*(Q1y-Q0y)-(P1y-P0y)*(Q1x-Q0x);
                    alpha_numerator =-((P0x-Q0x)*(Q1y-Q0y)-(P0y-Q0y)*(Q1x-Q0x));

                    if(flag_print==1)
                        _printf_("alpha_denominator=" << alpha_denominator << "\n");

                    /* the edges are parallel, do nothing (no intersection to compute) */
                    if(alpha_denominator==0.)
                    {
                        if(flag_print==1)
                            _printf_("parrallel edge");

                        /* Check if edge j in on the dead edge */
                        if((born_node_on_dead_node[j_born]!=j_dead) && (born_node_on_dead_node[j_born]!=jnext_dead))
                        {
                            alpha=-1.;
                            beta=-1.;
                            /* alpha in case of parrallel edge */
                            if(Q1x==Q0x)
                            {
                                if(P0x==Q0x)
                                    alpha=-2.;
                                else
                                    continue;
                            }
                            else
                                alpha=(P0x-Q0x)/(Q1x-Q0x);

                            /* beta in case of parrallel edge */
                            if(Q1y==Q0y)
                            {
                                if(P0y==Q0y)
                                    beta=-2.;
                                else
                                    continue;
                            }
                            else
                                beta=(P0y-Q0y)/(Q1y-Q0y);

                            if(beta==-2.)
                                if (alpha>0. && alpha<1.)
                                    born_node_on_dead_edge[j_born]=1;

                            if(alpha==-2.)
                                if(beta>0. && beta<1.)
                                    born_node_on_dead_edge[j_born]=1;

                            if(alpha!=-2. && beta!=-2. && alpha==beta && alpha>0. && alpha<1. && beta>0. && beta<1.)
                                born_node_on_dead_edge[j_born]=1;

                            if(flag_print==1)
                                _printf_("alpha=" << alpha << ", beta=" << beta << ",\n");
                        }

                        /* is the node jnext on the edge */
                        /* Check if edge jnext in on the dead edge */
                        if((born_node_on_dead_node[jnext]!=j_dead) && (born_node_on_dead_node[jnext]!=jnext_dead))
                        {
                            alpha=-1.;
                            beta=-1.;

                            /* alpha in case of parrallel edge */
                            if((Q1x-Q0x)==0.)
                            {
                                if(P1x==Q0x)
                                    alpha=-2.;
                                else
                                    continue;
                            }
                            else
                                alpha=(P1x-Q0x)/(Q1x-Q0x);

                            /* beta in case of parrallel edge */
                            if((Q1y-Q0y)==0.)
                            {
                                if(P1y==Q0y)
                                    beta=-2.;
                                else
                                    continue;
                            }
                            else
                                beta=(P1y-Q0y)/(Q1y-Q0y);

                            /* is the node on the edge */
                            if(beta==-2.)
                                if (alpha>0. && alpha<1.)
                                    born_node_on_dead_edge[jnext]=1;

                            if(alpha==-2.)
                                if(beta>0. && beta<1.)
                                    born_node_on_dead_edge[jnext]=1;

                            if(alpha!=-2. && beta!=-2. && alpha==beta  && alpha>0. && alpha<1. && beta>0. && beta<1.)
                                born_node_on_dead_edge[jnext]=1;

                            if(flag_print==1)
                                _printf_("alpha=" << alpha << ", beta=" << beta << ",\n");
                        }

                        /* Check if dead edge j_dead is on the born edge */
                        if((born_node_on_dead_node[j_born]!=j_dead) && (born_node_on_dead_node[jnext]!=j_dead))
                        {
                            alpha=-1.;
                            beta=-1.;
                            /* alpha in case of parrallel edge */
                            if(P1x==P0x)
                            {
                                if(Q0x==P0x)
                                    alpha=-2.;
                                else
                                    continue;
                            }
                            else
                                alpha=(Q0x-P0x)/(P1x-P0x);

                            /* beta in case of parrallel edge */
                            if(P1y==P0y)
                            {
                                if(Q0y==P0y)
                                    beta=-2.;
                                else
                                    continue;
                            }
                            else
                                beta=(Q0y-P0y)/(P1y-P0y);

                            if(beta==-2.)
                                if (alpha>0. && alpha<1.)
                                    dead_node_on_born_edge[j_dead]=1;

                            if(alpha==-2.)
                                if(beta>0. && beta<1.)
                                    dead_node_on_born_edge[j_dead]=1;

                            if(alpha!=-2. && beta!=-2. && alpha==beta && alpha>0. && alpha<1. && beta>0. && beta<1.)
                                dead_node_on_born_edge[j_dead]=1;

                            if(flag_print==1)
                                _printf_("alpha=" << alpha << ", beta=" << beta << ",\n");
                        }

                        /* Check if node jnext_dead in on the born edge */
                        if((born_node_on_dead_node[j_born]!=jnext_dead) && (born_node_on_dead_node[jnext]!=jnext_dead))
                        {
                            alpha=-1.;
                            beta=-1.;

                            /* alpha in case of parrallel edge */
                            if((P1x-P0x)==0.)
                            {
                                if(Q1x==P0x)
                                    alpha=-2.;
                                else
                                    continue;
                            }
                            else
                                alpha=(Q1x-P0x)/(P1x-P0x);

                            /* beta in case of parrallel edge */
                            if((P1y-P0y)==0.)
                            {
                                if(Q1y==P0y)
                                    beta=-2.;
                                else
                                    continue;
                            }
                            else
                                beta=(Q1y-P0y)/(P1y-P0y);

                            /* is the node on the edge */
                            if(beta==-2.)
                                if (alpha>0. && alpha<1.)
                                    dead_node_on_born_edge[jnext_dead]=1;

                            if(alpha==-2.)
                                if(beta>0. && beta<1.)
                                    dead_node_on_born_edge[jnext_dead]=1;

                            if(alpha!=-2. && beta!=-2. && alpha==beta  && alpha>0. && alpha<1. && beta>0. && beta<1.)
                                dead_node_on_born_edge[jnext_dead]=1;

                            if(flag_print==1)
                                _printf_("alpha=" << alpha << ", beta=" << beta << ",\n");
                        }

                        /* end of this check in case of parallel edge */
                        continue;
                    }

                     /* check if the edges are connected (no other intersection to compute) */
                    if(born_node_on_dead_node[j_born]==j_dead    )
                        continue;
                    if(born_node_on_dead_node[j_born]==jnext_dead)
                        continue;
                    if(born_node_on_dead_node[jnext]==j_dead    )
                        continue;
                    if(born_node_on_dead_node[jnext]==jnext_dead)
                        continue;

                    /* alpha and beta defines the position of the intersection */
                    /* relatively to the born and dead edge, respectively. */


                    alpha=alpha_numerator/alpha_denominator;


                    if(fabs(Q1x-Q0x)<fabs(Q1y-Q0y))
                        beta=((P0y-Q0y)+alpha*(P1y-P0y))/(Q1y-Q0y);
                    else
                        beta=((P0x-Q0x)+alpha*(P1x-P0x))/(Q1x-Q0x);

                    /*if(flag_print==1)
                            _printf_("interesct: alpha=" << alpha << ", beta=" << beta << ", beta-1.=" << (beta-1.0)*1e12 << "\n");
                        */

                    /* Check on the alpha and beta */
                    if((alpha>=0.-alpha_tol) && (alpha<=1.+alpha_tol) && (beta>=0.-alpha_tol) && (beta<=1.+alpha_tol))
                    {
                        /*if(flag_print==1)
                            _printf_("select interesct: alpha=" << alpha << ", beta=" << beta << ", beta-1.=" << beta-1.0) << "\n");
                        */
                        if((alpha<=alpha_tol) && (beta>alpha_tol) && (beta<(1.-alpha_tol))) /* add intersection if not on a dead node */
                            if(born_node_on_dead_node[j_born]<0) /* thick check may be not useful */
                                born_node_on_dead_edge[j_born]=1;

                        if((alpha>=(1.-alpha_tol)) && (beta>alpha_tol) && (beta<(1.-alpha_tol))) /* add intersection if not on a dead node */
                            if(born_node_on_dead_node[jnext]<0) /* thick check may be not useful */
                                born_node_on_dead_edge[jnext]=1;

                        if((beta<=alpha_tol) && (alpha>alpha_tol) && (alpha<(1.-alpha_tol))) /* add intersection if not on a dead node */
                            dead_node_on_born_edge[j_dead]=1;

                        if((alpha>alpha_tol) && (alpha<(1.-alpha_tol)) && (beta>=(1.-alpha_tol))) /* add intersection if not on a dead node */
                            dead_node_on_born_edge[jnext_dead]=1;

                        if((beta<(1.-alpha_tol)) && (alpha>alpha_tol) && (alpha<(1.-alpha_tol)) && (beta>alpha_tol)) /* add intersection */
                        {
                            if(flag_print==1)
                                _printf_("alpha=" <<  alpha << "");

                            alpha_tmp[j_dead]=alpha;
                        }

                    } /* else, do nothing, no intersection */

                } /* End loop over the dead edges */

                /* storing information in born_edge_intersect_dead_edge */
                for (j_dead=0; j_dead<3; j_dead++)
                    born_edge_intersect_dead_edge[j_born*3+j_dead]=alpha_tmp[j_dead];

            } /* End loop over the born edges */

            /*-------- Building the intersection (a convex polygon) ---------- */

            /* List of interesection points (shared node, inside node, or edge to edge intersection) */
            nb_intersection_points=0;

            /* 1) Shared nodes: */
            for (j_born=0; j_born<3; j_born++)
            {
                if(flag_print==1)
                    _printf_("j_born=" << j_born << ",born_node_on_dead_node[j]=" <<  born_node_on_dead_node[j_born] << "\n");

                if(born_node_on_dead_node[j_born]>=0)
                {

                    intersection_points_x[nb_intersection_points]=x_born_elements[i_born*3+j_born];
                    intersection_points_y[nb_intersection_points]=y_born_elements[i_born*3+j_born];
                    nb_intersection_points++;
                }
            }

            /* 2) Born nodes included in dead element: */
            for (j_born=0; j_born<3; j_born++)
            {
                if(flag_print==1)
                    _printf_("j_born=" << j_born << ",born_node_in_dead_element[j_born]=" << born_node_in_dead_element[j_born] << "\n");

                if(born_node_in_dead_element[j_born]>0)
                {
                    intersection_points_x[nb_intersection_points]=x_born_elements[i_born*3+j_born];
                    intersection_points_y[nb_intersection_points]=y_born_elements[i_born*3+j_born];
                    nb_intersection_points++;
                }
            }

            /* 3) Dead nodes included in born element: */
            for (j_dead=0; j_dead<3; j_dead++)
            {
                if(flag_print==1)
                    _printf_("j_dead=" << j_dead << ",dead_node_in_born_element[j_dead]=" << dead_node_in_born_element[j_dead] << "\n");

                if(dead_node_in_born_element[j_dead]>0)
                {
                    intersection_points_x[nb_intersection_points]=x_dead_elements[i_dead*3+j_dead];
                    intersection_points_y[nb_intersection_points]=y_dead_elements[i_dead*3+j_dead];
                    nb_intersection_points++;
                }
            }

            /* 4) Edge to edge intersections: */
            for (j_born=0; j_born<3; j_born++)
            {
                for (j_dead=0; j_dead<3; j_dead++)
                {
                    if(flag_print==1)
                        _printf_("j_born=" << j_born << ",j_dead=" << j_dead
                                << ",born_edge_intersect_dead_edge[j_born*3+j_dead]="
                                << born_edge_intersect_dead_edge[j_born*3+j_dead] << "\n");

                    if(born_edge_intersect_dead_edge[j_born*3+j_dead]!=-1.)
                    {
                        alpha=born_edge_intersect_dead_edge[j_born*3+j_dead];

                        /* the edge j goes from node j to node jnext */
                        divresult = div(j_born+1,3);
                        jnext=divresult.rem;

                        P0x=x_born_elements[i_born*3+j_born];
                        P0y=y_born_elements[i_born*3+j_born];
                        P1x=x_born_elements[i_born*3+jnext];
                        P1y=y_born_elements[i_born*3+jnext];

                        Xx=P0x+alpha*(P1x-P0x);
                        Xy=P0y+alpha*(P1y-P0y);

                        intersection_points_x[nb_intersection_points]=Xx;
                        intersection_points_y[nb_intersection_points]=Xy;
                        nb_intersection_points++;
                    }
                }
            }

            /* 4 bis) Born nodes on dead edge: */
            for (j_born=0; j_born<3; j_born++)
            {
                if(born_node_on_dead_edge[j_born]>0)
                {
                    intersection_points_x[nb_intersection_points]=x_born_elements[i_born*3+j_born];
                    intersection_points_y[nb_intersection_points]=y_born_elements[i_born*3+j_born];
                    nb_intersection_points++;
                }
            }

            /* 4 tris) Dead nodes on born edge: */
            for (j_dead=0; j_dead<3; j_dead++)
            {
                if(dead_node_on_born_edge[j_dead]>0)
                {
                    intersection_points_x[nb_intersection_points]=x_dead_elements[i_dead*3+j_dead];
                    intersection_points_y[nb_intersection_points]=y_dead_elements[i_dead*3+j_dead];
                    nb_intersection_points++;
                }
            }

            if(flag_print==1)
                for (i_intersection_points=0; i_intersection_points<nb_intersection_points; i_intersection_points++)
                    _printf_("intersection_points_x[i_intersection_points]="
                            << intersection_points_x[i_intersection_points] << ",intersection_points_y[i_intersection_points]="
                            << intersection_points_y[i_intersection_points] << "\n");

            /* Building the polygon as in https://hal.archives-ouvertes.fr/inria-00354509/document */
            /* "The convex polygon is meshed by primarily constructing an oriented triangle with three points chosen randomly. */
            /* Notice that the three points cannot be aligned by construction. */
            /* Then, a new point is selected and the unique triangle */
            /* edge which is viewing this point is looked for, i.e., the only edge for which the barycentric coordinate is negative. */
            /* A new oriented triangle is built by connecting the point to this edge. */
            /* The process is iterated until all points are inserted by only checking edges forming the boundary of the current polygon */
            /* (i.e., edges which are not shared by two triangles). */
            /* Notice that, at each iteration, as the current polygon is convex, there is only one boundary edge that views the selected point. */
            /* A mesh of the polygon is then obtained." */

            /* Resetting the list of polygon points */
            nb_polygon_points=0;
            for (i_polygon_points=0; i_polygon_points<6; i_polygon_points++)
                next_polygon_point[i_polygon_points]=0;

            if(flag_print==1)
                _printf_("i_born=" << i_born << ",i_dead=" << i_dead << "\n");

            if(flag_print==1)
                _printf_("nb_intersection_points=" << nb_intersection_points << "\n");

            /* Defining a polygon with the intersection points */
            if(nb_intersection_points>=3)
            {

                /* we select 3 nodes that are far away one from each other */

                /* select the node that is the farest from the first node */
                dist_max=0.;
                ind_max=0;
                for (i_intersection_points=1; i_intersection_points<nb_intersection_points; i_intersection_points++)
                {
                    dist=pow(pow(intersection_points_x[i_intersection_points]-intersection_points_x[0],2.)
                            +pow(intersection_points_y[i_intersection_points]-intersection_points_y[0],2.),0.5);
                    if(dist>dist_max)
                    {
                        dist_max=dist;
                        ind_max=i_intersection_points;
                    }
                }
                tmp_point_x=intersection_points_x[1];
                tmp_point_y=intersection_points_y[1];

                intersection_points_x[1]=intersection_points_x[ind_max];
                intersection_points_y[1]=intersection_points_y[ind_max];

                intersection_points_x[ind_max]=tmp_point_x;
                intersection_points_y[ind_max]=tmp_point_y;

                /* select the node whose the summed distance to the first and the select node is the largest */
                dist_max=0.;
                ind_max=0;
                for (i_intersection_points=2; i_intersection_points<nb_intersection_points; i_intersection_points++)
                {
                    dist= pow(pow(intersection_points_x[i_intersection_points]-intersection_points_x[0],2.)
                            +pow(intersection_points_y[i_intersection_points]-intersection_points_y[0],2.),0.5);
                    dist+=pow(pow(intersection_points_x[i_intersection_points]-intersection_points_x[1],2.)
                            +pow(intersection_points_y[i_intersection_points]-intersection_points_y[1],2.),0.5);
                    if(dist>dist_max)
                    {
                        dist_max=dist;
                        ind_max=i_intersection_points;
                    }
                }
                tmp_point_x=intersection_points_x[2];
                tmp_point_y=intersection_points_y[2];

                intersection_points_x[2]=intersection_points_x[ind_max];
                intersection_points_y[2]=intersection_points_y[ind_max];

                intersection_points_x[ind_max]=tmp_point_x;
                intersection_points_y[ind_max]=tmp_point_y;

                /* Select the first three nodes of the intersection */
                for (i_intersection_points=0; i_intersection_points<3; i_intersection_points++)
                {
                    polygon_points_x[i_intersection_points]=intersection_points_x[i_intersection_points];
                    polygon_points_y[i_intersection_points]=intersection_points_y[i_intersection_points];
                    nb_polygon_points++;
                    if(i_intersection_points>0)
                        next_polygon_point[i_intersection_points-1]=i_intersection_points;
                }

                /* compute the sign of the polygon */
                area_polygon=0.;
                area_polygon+= (polygon_points_x[0]*polygon_points_y[1]-polygon_points_y[0]*polygon_points_x[1]);
                area_polygon+= (polygon_points_x[1]*polygon_points_y[2]-polygon_points_y[1]*polygon_points_x[2]);
                area_polygon+= (polygon_points_x[2]*polygon_points_y[0]-polygon_points_y[2]*polygon_points_x[0]);
                sign_polygon=sign(area_polygon);

                /* While loop to treat the remaining intersection nodes */
                while(i_intersection_points<nb_intersection_points)
                {
                    x_target=intersection_points_x[i_intersection_points];
                    y_target=intersection_points_y[i_intersection_points];

                    /* Loop over the edge of the current polygon */
                    i_polygon_points=0;
                    for (j_edge_polygon=0; j_edge_polygon<nb_polygon_points; j_edge_polygon++)
                    {
                        /* Defining the next polygon point */
                        inext_polygon_points=next_polygon_point[i_polygon_points];

                        /* vector from born node j to born node */
                        x_vect_target=x_target-polygon_points_x[i_polygon_points];
                        y_vect_target=y_target-polygon_points_y[i_polygon_points];

                        /* vector from j_edge_polygon to j_edge_polygon+1 */
                        x_edge_polygon=polygon_points_x[inext_polygon_points]-polygon_points_x[i_polygon_points];
                        y_edge_polygon=polygon_points_y[inext_polygon_points]-polygon_points_y[i_polygon_points];

                        /* add the sign of the cross product */
                        cross_prod=x_vect_target*y_edge_polygon-x_edge_polygon*y_vect_target;
                        /*if(flag_print==1)
                            _printf_("cross_prod=" << << ",sign(cross_prod)=" << << ", sign_polygon=" << << "\n", cross_prod, sign(cross_prod), sign_polygon);
                        */
                        if(cross_prod!=0. && sign(cross_prod)==sign_polygon)
                        {
                            polygon_points_x[nb_polygon_points]=x_target;
                            polygon_points_y[nb_polygon_points]=y_target;
                            next_polygon_point[nb_polygon_points]=inext_polygon_points;
                            next_polygon_point[i_polygon_points]=nb_polygon_points;
                            nb_polygon_points++;
                            break;
                        }

                        i_polygon_points=inext_polygon_points;
                    }
                    i_intersection_points++;
                }

                /*if(flag_print==1)
                    _printf_("nb_polygon_points=" << << "\n", nb_polygon_points);*/

                /* Compute the area of this convex polygon by */
                /* first meshing it with triangles and then summing */
                /* the area of the triangles */
                area_intersection=0.;

                i_polygon_points=0;
                for (j_edge_polygon=0; j_edge_polygon<nb_polygon_points; j_edge_polygon++)
                {
                    /*if(flag_print==1)
                        _printf_("polygon_points_x[i_polygon_points]=" << << ",polygon_points_y[i_polygon_points]=" << << "\n",polygon_points_x[i_polygon_points],polygon_points_y[i_polygon_points]);
                    */
                    /* Defining the next polygon point */
                    inext_polygon_points=next_polygon_point[i_polygon_points];

                    /*if(flag_print==1)
                        _printf_("i_polygon_points=" << << ", inext_polygon_points=" << << "\n", i_polygon_points, inext_polygon_points);
                    if(flag_print==1)
                        _printf_("area_intersection=" << << ", darea=" << << "\n", area_intersection, (polygon_points_x[i_polygon_points]*polygon_points_y[inext_polygon_points]-polygon_points_y[i_polygon_points]*polygon_points_x[inext_polygon_points]));
                    */
                    area_intersection+=(polygon_points_x[i_polygon_points]*polygon_points_y[inext_polygon_points]-polygon_points_y[i_polygon_points]*polygon_points_x[inext_polygon_points]);

                    i_polygon_points=inext_polygon_points;
                }
                area_intersection=fabs(area_intersection/2.);
            }
            else /* no triangle interesection */
                area_intersection=0.;


            /*if(flag_print==1)
                _printf_("area_intersection=" << << "\n \n", area_intersection);*/

            /*-------- Computing the integrated quantities ---------- */


            tmp_distance= pow(
                pow(
                     (x_dead_elements[i_dead*3+0]+x_dead_elements[i_dead*3+1]+x_dead_elements[i_dead*3+2])/3.
                    -(x_born_elements[i_born*3+0]+x_born_elements[i_born*3+1]+x_born_elements[i_born*3+2])/3.
                    ,2.)
                +pow(
                     (y_dead_elements[i_dead*3+0]+y_dead_elements[i_dead*3+1]+y_dead_elements[i_dead*3+2])/3.
                    -(y_born_elements[i_born*3+0]+y_born_elements[i_born*3+1]+y_born_elements[i_born*3+2])/3.
                    ,2.)
                              ,0.5);

            if(min_distance<0. || tmp_distance<min_distance)
                min_distance = tmp_distance;

            for (k=0; k<nb_variables; k++)
            {
                if(method_in[k]==1)
                    integrated_variables[k]+=area_intersection*dead_variables[i_dead*nb_variables+k];
                else
                    if(min_distance == tmp_distance)
                        integrated_variables[k]=dead_variables[i_dead*nb_variables+k];
            }

            test_integrated_area+=area_intersection;
        } /* end of loop on i_dead */

        /* Weighted average is computed as the */
        /* integrated variables divided by integrated area */
        for (k=0; k<nb_variables; k++)
            if(method_in[k]==1)
                tmp_mean_variables[i_born*nb_variables+k]=integrated_variables[k]/test_integrated_area;
            else
                tmp_mean_variables[i_born*nb_variables+k]=integrated_variables[k];

        tmp_integrated_area[i_born]=test_integrated_area;

        if(test_integrated_area<=0.)
            _printf_("born_elements[i_born]=" << born_elements[i_born] << "\n");

    } /* end of loop on i */


    // destroy temporary arrays
    xDelete<int>(dead_num_node);
    xDelete<int>(born_num_node);

    xDelete<double>(x_dead_elements);
    xDelete<double>(x_born_elements);
    xDelete<double>(y_dead_elements);
    xDelete<double>(y_born_elements);

    xDelete<double>(dead_variables);
    xDelete<double>(integrated_variables);


    return 1;
}
