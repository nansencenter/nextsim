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

int InterpFromMeshToMesh2dCavities(double** pdata_interp,double* IntMatrix_in,int nb_variables,
		double* surface_old, double* surface_new, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new){
    /*  Conservative interpolation method
        Use the cavities detected by detect_cavity to transfer physical
        quantities from the old element to the new ones*/

    int verbosity;
    
    /*Bamg options*/
    verbosity=2;//bamgopts->verbose;

    /* Which conservative method to be used:
      method 0 is simply averaging over the whole cavity (conservative but highly diffusive)
      method 1 computes the intersection between born and dead triangles (conservative and non-diffusive)*/
    bool method=0; /* method 1 should be used, method 0 is kept for comparison and debugging*/

	/*Output*/
	double* IntMatrix_out=NULL;

    /*Intermediary*/
    int new_nb_elements=bamgmesh_new->TrianglesSize[0];

    /*Initialize output*/
    IntMatrix_out=xNew<double>(new_nb_elements*nb_variables);

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
    ulong new_element_i;
    ulong old_element_i;
    for(ulong i=0; i<gate.nb_matching_elements; i++)
    {
        new_element_i=gate.new_elements[i];
        old_element_i=gate.old_elements[i];
        for(int j=0; j<nb_variables; j++)
            IntMatrix_out[(new_element_i-1)*nb_variables+j]=IntMatrix_in[(old_element_i-1)*nb_variables+j];
    }
    
    if (verbosity>1) _printf_("   Interp Cavities: initialize intermediary...\n");
    ulong  *dead_elements = xNew<ulong>(gate.max_size_dead_cavity);
    ulong  *born_elements = xNew<ulong>(gate.max_size_born_cavity);
    if (verbosity>1) _printf_("nb_dead_elements, nb_born_elements, tmp_element" << "\n");
    ulong nb_dead_elements, nb_born_elements;
    ulong tmp_element;
    double integrated_area, new_integrated_area, correction_factor;
    if (verbosity>1) _printf_("nb_variables" << nb_variables << "\n");
    double *integrated_variables = xNew<double>(nb_variables);


    if (verbosity>1) _printf_("   Interp Cavities: loop1...\n");

    for(ulong cavity_id=0; cavity_id< gate.nb_cavities; cavity_id++) 
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

        if (verbosity>1) _printf_("   Interp Cavities: integrated area, nb_dead_elements:" << nb_dead_elements << "\n");

        for(ulong i=0; i<nb_dead_elements; i++)
        {
            tmp_element=gate.element_dead_cavity[cavity_id*gate.max_size_dead_cavity+i];
            if (verbosity>1) _printf_("tmp_element: " << tmp_element << "\n");

            dead_elements[i]=tmp_element;
            integrated_area+=surface_old[tmp_element-1];
        }

        if (verbosity>1) _printf_("   Interp Cavities: integrated area...\n");
        for(ulong i=0; i<nb_born_elements; i++)
        {
            tmp_element=gate.element_born_cavity[cavity_id*gate.max_size_born_cavity+i];
            born_elements[i]=tmp_element;
            new_integrated_area+=surface_new[tmp_element-1];
        }
        
        if (new_integrated_area<=0.)  _error_("new_integrated_area should not be larger than 0. Here it is: " << new_integrated_area<<" ");
        if (integrated_area<=0.)  _error_("integrated_area should not be larger than 0. Here it is: " << integrated_area<<" ");
        
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


            for(ulong i=0; i<nb_dead_elements; i++)
            {
                tmp_element=dead_elements[i];
                for(int j=0; j<nb_variables; j++)
                    integrated_variables[j]+=surface_old[tmp_element-1]*IntMatrix_in[(tmp_element-1)*nb_variables+j];
            }
            
            if (verbosity>1) _printf_("   Interp Cavities: InrMatrix_out..\n");
            for(ulong i=0; i<nb_born_elements; i++)
            {
                tmp_element=born_elements[i];
                for(int j=0; j<nb_variables; j++)
                    IntMatrix_out[(tmp_element-1)*nb_variables+j]=integrated_variables[j]/integrated_area*correction_factor;
            }
        }
        #if 0           
        else
            /* Second method (fully conservative)*/
            
            /* preparation*/
            dead_num_node=element_old.num_node(dead_elements,:);
            born_num_node=element.num_node(born_elements,:);
            
            x_dead_elements=mesh_old.node.x(dead_num_node);
            y_dead_elements=mesh_old.node.y(dead_num_node);
            
            x_born_elements=mesh.node.x(born_num_node);
            y_born_elements=mesh.node.y(born_num_node);
            
            dead_variables=IntMatrix_in(dead_elements,:);
            
            /* interpolation over the cavity*/
            [tmp_mean_variables,tmp_integrated_area]=interp_cavity_mex(...
                nb_dead_elements,...
                nb_born_elements,...
                nb_variables,...
                gate.new_bamg_mesh_Nn,...
                uint32(dead_elements'),...
                uint32(born_elements'),...
                uint32(dead_num_node'),...
                uint32(born_num_node'),...
                uint32(PreviousNumbering'),...
                x_dead_elements',...
                x_born_elements',...
                y_dead_elements',...
                y_born_elements',...
                dead_variables',...
                uint32(0));
            
            /* reshaping */
            tmp_mean_variables=tmp_mean_variables;
            
            /* Sanity check
              does the integrated area corresponds to the area of the born element*/
            area_error_tolerance=1000; 
            area_error=abs(tmp_integrated_area-new_element_area(born_elements)');
            if(max(area_error)>area_error_tolerance)
                ind_error=find(area_error>area_error_tolerance);
                
                disp('Indices of the problematic born elements')
                born_elements(ind_error)
                
                % Rerunning interpolation over the cavity with debugging
                % activated
                [tmp_mean_variables,tmp_integrated_area]=interp_cavity_mex(...
                    nb_dead_elements,...
                    nb_born_elements,...
                    nb_variables,...
                    gate.new_bamg_mesh_Nn,...
                    uint32(dead_elements'),...
                    uint32(born_elements'),...
                    uint32(dead_num_node'),...
                    uint32(born_num_node'),...
                    uint32(PreviousNumbering'),...
                    x_dead_elements',...
                    x_born_elements',...
                    y_dead_elements',...
                    y_born_elements',...
                    dead_variables',...
                    uint32(born_elements(ind_error(1))));
                
                tmp_integrated_area(ind_error)
                new_element_area(born_elements(ind_error))
                
                figure
                patch(x_dead_elements',y_dead_elements',0*y_dead_elements','EdgeColor','green','LineWidth',16,'FaceAlpha',0.0)
                hold on;
                patch(x_born_elements',y_born_elements',0*y_born_elements','EdgeColor','yellow','LineWidth',8,'FaceAlpha',0.0)
                
                text(x_dead_elements(:),y_dead_elements(:)+1,num2str(dead_num_node(:)),'Color','black')
                
                for i=1:length(dead_elements)
                    text(mean(x_dead_elements(i,:)),mean(y_dead_elements(i,:)),num2str(dead_elements(i)),'Color','red')
                end
                    
                wrong_x=x_born_elements(ind_error',:);
                wrong_y=y_born_elements(ind_error',:);
                wrong_num=born_num_node(ind_error',:);
                text(wrong_x(:),wrong_y(:)+0.1,num2str(wrong_num(:)),'Color','blue');
                
                error(' Conservation issue!!! The area of the new element does not match with the sum of the triangle intersections...')
            end
            %-------- Parrallel stuff (may be not useful in case parfor is not used) ----------
            
            % store information in the shared vectors (Not useful if parfor id not used)
            par_born_elements{cavity_id,1}=born_elements;
            par_mean_variables{cavity_id,1}=tmp_mean_variables*correction_factor; 
        end
        #endif
    }

    if (verbosity>1) _printf_("   Interp Cavities: end!\n");

    /*---- end of the interpolation cavity by cavity      ------*/
	*pdata_interp=IntMatrix_out;
	return 1;
}

int DetectCavities(InterpFromMeshToMesh2dCavitiesThreadStruct* gate, BamgMesh* bamgmesh_old,BamgMesh* bamgmesh_new){

    int verbosity;
    
    /*Bamg options*/
    verbosity=2;//bamgopts->verbose;

    if (verbosity>1) _printf_("   Detect_cavities: Read input...\n");
	/*---------- Input  ----------*/

	double *bamg_mesh_Triangles                 = bamgmesh_old->Triangles;
    double *bamg_mesh_NodalElementConnectivity  = bamgmesh_old->NodalElementConnectivity;
    double *bamg_mesh_NodalConnectivity         = bamgmesh_old->NodalConnectivity ;
    double *bamg_mesh_ElementConnectivity       = bamgmesh_old->ElementConnectivity;
    
    ulong bamg_mesh_Nn                         = (ulong)bamgmesh_old->VerticesSize[0];
    ulong bamg_mesh_NNodalConnectivity         = (ulong)bamgmesh_old->NodalConnectivitySize[1];
    ulong bamg_mesh_NNodalElementConnectivity  = (ulong)bamgmesh_old->NodalElementConnectivitySize[1];
    ulong bamg_mesh_Ne                         = (ulong)bamgmesh_old->TrianglesSize[0];
    ulong bamg_mesh_NVerticesOnGeomVertex      = (ulong)bamgmesh_old->VerticesOnGeomVertexSize[0];
    
    double *new_bamg_mesh_Triangles             = bamgmesh_new->Triangles;
    double *new_bamg_mesh_ElementConnectivity   = bamgmesh_new->ElementConnectivity;
    double *new_bamg_mesh_PreviousNumbering     = bamgmesh_new->PreviousNumbering;
    
    ulong new_bamg_mesh_Nn                     = (ulong)bamgmesh_new->VerticesSize[0];
    ulong new_bamg_mesh_NVerticesOnGeomVertex  = (ulong)bamgmesh_new->VerticesOnGeomVertexSize[0];
    ulong new_bamg_mesh_Ne                     = (ulong)bamgmesh_new->TrianglesSize[0];


    if (verbosity>1) _printf_("   bamg_mesh_Nn = " << bamg_mesh_Nn << "\n");
    if (verbosity>1) _printf_("   bamg_mesh_NNodalConnectivity = " << bamg_mesh_NNodalConnectivity << "\n");
    if (verbosity>1) _printf_("   bamg_mesh_NNodalElementConnectivity = " << bamg_mesh_NNodalElementConnectivity << "\n");
    if (verbosity>1) _printf_("   bamg_mesh_Ne = " << bamg_mesh_Ne << "\n");
    if (verbosity>1) _printf_("   bamg_mesh_NVerticesOnGeomVertex = " << bamg_mesh_NVerticesOnGeomVertex << "\n");

    if (verbosity>1) _printf_("   new_bamg_mesh_Nn = " << new_bamg_mesh_Nn << "\n");
    if (verbosity>1) _printf_("   new_bamg_mesh_NVerticesOnGeomVertex = " << new_bamg_mesh_NVerticesOnGeomVertex << "\n");
    if (verbosity>1) _printf_("   new_bamg_mesh_Ne = " << new_bamg_mesh_Ne << "\n");

    if (verbosity>1) _printf_("   Detect_cavities: Initialize arrays...\n");

    /*if(bamg_mesh_NVerticesOnGeomVertex!=new_bamg_mesh_NVerticesOnGeomVertex)
        _error_("not the same number of vertices on the Geom in the old and new mesh");*/


    /*---------- Output arrays      ----------*/
    ulong  *PreviousNumbering = xNew<ulong>(new_bamg_mesh_Nn);

    /*---------- Intermediary arrays---------- */
    ulong *tmp_old_elements = xNew<ulong>(bamg_mesh_Ne);
    ulong *tmp_new_elements = xNew<ulong>(new_bamg_mesh_Ne);
    long  *cavity_number_dead = xNew<long>(bamg_mesh_Ne);
    long  *cavity_number_born = xNew<long>(new_bamg_mesh_Ne);
    long  *is_born_element_to_treat = xNew<long>(new_bamg_mesh_Ne); /* boolean to know if the born element is to be treated */      
    long  *is_dead_element_to_treat = xNew<long>(bamg_mesh_Ne);   /* boolean to know if the dead element is to be treated */ 

    if (verbosity>1) _printf_("   Detect_cavities: Initialize local variables...\n");
    /*---------- Local variables    ----------*/
    int l, l0, l1, l2, j_next;
    
    ulong nb_matching_elements, nb_cavities;
    ulong max_size_born_cavity, max_size_dead_cavity;
    
    int edge_node_1, edge_node_2;
    ulong Previous_edge_node_1, Previous_edge_node_2;
    ulong tmp_node_1, tmp_node_2;
    int Previous_connected_element;
    ulong tmp_born, tmp_dead;
    ulong Previous_num_node[3], Previous_connected_element_num_node[3];
    ulong max_candidate_element_to_treat=10000, number_candidate_element_to_treat, number_candidate_element_treated;    
    ulong  *candidate_element_to_treat = xNew<ulong>(max_candidate_element_to_treat);   /* boolean to know if the dead element is to be treated */ 

    int candidate_dead_element[2], nb_candidate_dead_element;
    int nb_dead, nb_born, nb_connected_nodes;
    ulong neighbour;
    long new_cavity_number;
    int i_dead_element_node_1, dead_element_node_1;
    int i_dead_element_node_2, dead_element_node_2;
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
    for (ulong i=0; i<new_bamg_mesh_NVerticesOnGeomVertex; i++)
    {
        if (verbosity>1) _printf_(" i " << i << "(ulong) new_bamg_mesh_PreviousNumbering[i] " << (ulong) new_bamg_mesh_PreviousNumbering[i] << "\n");
        PreviousNumbering[i]=i+1;
    }
    if (verbosity>1) _printf_("   Detect_cavities: loop1...\n");
    for (ulong i=new_bamg_mesh_NVerticesOnGeomVertex; i<new_bamg_mesh_Nn; i++)
    {
        if (verbosity>1) _printf_(" i " << i << "(ulong) new_bamg_mesh_PreviousNumbering[i] " << (ulong) new_bamg_mesh_PreviousNumbering[i] << "\n");
        PreviousNumbering[i]=(ulong) new_bamg_mesh_PreviousNumbering[i];
    }
    
          
    if (verbosity>1) _printf_("   Detect_cavities: find the correspondance...\n");
    /*  To find the correspondance between the new and old elements */
    nb_matching_elements=0;

    for (ulong i=0; i<new_bamg_mesh_Ne; i++)
    {
 /*       if (verbosity>1) _printf_("   Detect_cavities: Previous_num_node...\n");


            _printf_("   new_bamg_mesh_Triangles = " << new_bamg_mesh_Triangles[4*i+0] << "\n");
            _printf_("   new_bamg_mesh_Triangles = " << new_bamg_mesh_Triangles[4*i+1] << "\n");
            _printf_("   new_bamg_mesh_Triangles = " << new_bamg_mesh_Triangles[4*i+2] << "\n");
            _printf_("   new_bamg_mesh_Triangles = " << new_bamg_mesh_Triangles[4*i+3] << "\n");
            _printf_("   new_bamg_mesh_Triangles = " << bamgmesh_old->Triangles[4*i+0] << "\n");
            _printf_("   new_bamg_mesh_Triangles = " << bamgmesh_old->Triangles[4*i+1] << "\n");
            _printf_("   new_bamg_mesh_Triangles = " << bamgmesh_old->Triangles[4*i+2] << "\n");
            _printf_("   new_bamg_mesh_Triangles = " << bamgmesh_old->Triangles[4*i+3] << "\n");*/

        for (int j=0; j<3; j++)
        {   
           /* _printf_("   i = " << i << "\n");
            _printf_("   j = " << j << "\n");
            _printf_("   new_bamg_mesh_Triangles = " << new_bamg_mesh_Triangles[4*i+j] << "\n");
            _printf_("   new_bamg_mesh_Triangles = " << (ulong) new_bamg_mesh_Triangles[4*i+j] << "\n");
            _printf_("   Previous = " << PreviousNumbering[(ulong) new_bamg_mesh_Triangles[4*i+j]-1] << "\n");*/
            Previous_num_node[j]=PreviousNumbering[(ulong) new_bamg_mesh_Triangles[4*i+j]-1];
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
                    Previous_connected_element=(ulong) bamg_mesh_NodalElementConnectivity[(Previous_num_node[j]-1)*bamg_mesh_NNodalElementConnectivity+k];
                    
                    /* If NaN in matlab (0 after calling uint32)  */
                    if(Previous_connected_element==0)
                        break;
                    
                    for (l=0; l<3; l++)
                        Previous_connected_element_num_node[l]=(ulong) bamg_mesh_Triangles[(Previous_connected_element-1)*4+l];
                    
                    for (l=0; (l<3)  && (found_intersect==0); l++)
                    {
                        /* l0, l1, l2 */
                        l0=l;
                        l1=next(l0);
                        l2=next(l1);
                       
                        if((Previous_connected_element_num_node[0]==Previous_num_node[l0]) && (Previous_connected_element_num_node[1]==Previous_num_node[l1]) && (Previous_connected_element_num_node[2]==Previous_num_node[l2]))
                        {
                            tmp_old_elements[nb_matching_elements]=Previous_connected_element;
                            tmp_new_elements[nb_matching_elements]=i+1;
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
    ulong *old_elements = xNew<ulong>(nb_matching_elements);
    ulong *new_elements = xNew<ulong>(nb_matching_elements);

    /*---------- Intermediary arrays---------- */
    ulong *is_not_dead = xNew<ulong>(bamg_mesh_Ne);
    ulong *is_not_born = xNew<ulong>(new_bamg_mesh_Ne);
    
    ulong *ind_dead = xNew<ulong>(bamg_mesh_Ne);
    ulong *ind_born = xNew<ulong>(new_bamg_mesh_Ne);
    
    ulong *boundary_between_two_cavities = xNew<ulong>(bamg_mesh_Ne*3);

    if (verbosity>1) _printf_("   nb_matching_elements = " << nb_matching_elements << "\n");

    if (verbosity>1) _printf_("   Detect_cavities: loop2...\n");
    for (int i=0; i<nb_matching_elements; i++)
    {
        old_elements[i]=tmp_old_elements[i];
        new_elements[i]=tmp_new_elements[i]; 
        
        is_not_dead[tmp_old_elements[i]-1]=1;
        is_not_born[tmp_new_elements[i]-1]=1;
    }
    
    if (verbosity>1) _printf_("   Detect_cavities: destroy...\n");
    /* destroy the temporary matrix */
    delete tmp_old_elements;
    delete tmp_new_elements;

    if (verbosity>1) _printf_("   Detect_cavities: loop3...\n");
    /*  The old elements without correspondance are called "dead element" */
    /*  The new elements without correspondance are called "born element" */
    nb_dead=0;
    for (ulong i=0; i<bamg_mesh_Ne; i++)
    {
        cavity_number_dead[i]=0;
        if(is_not_dead[i]!=1)
        {
            ind_dead[nb_dead]=i+1;
            cavity_number_dead[i]=-1;
            nb_dead++;
        }
    }

    if (verbosity>1) _printf_("   Detect_cavities: loop4...\n");
    nb_born=0;
    for (ulong i=0; i<new_bamg_mesh_Ne; i++)
    {
        cavity_number_born[i]=0;
        if(is_not_born[i]!=1)
        {
            ind_born[nb_born]=i+1;
            cavity_number_born[i]=-1;
            nb_born++;
        }
    }
    if (verbosity>1) _printf_("   nb_dead = " << nb_dead << "\n");
    if (verbosity>1) _printf_("   nb_born = " << nb_born << "\n");

    if (verbosity>1) _printf_("   Detect_cavities: destroy2...\n");
    /* destroy the temporary matrix */
    delete is_not_dead;
    delete is_not_born;

    if (verbosity>1) _printf_("   Detect_cavities: detect the cavities...\n");
    /*  Loop that detect the cavities */
    new_cavity_number=0; 
    for (ulong i_born=0; i_born<nb_born; i_born++)
        {
        /*  first element to be checked */
        tmp_born=ind_born[i_born];
        
        /*  if the born element is already treated, do nothing */
        if(cavity_number_born[tmp_born-1]!=-1)
            continue;
        
        /*  We start the detection of a new cavity */
         new_cavity_number++;
        
        /*  we will store the candidate element to treat in an expandable vector, whose */
        /*  initial size is max_candidate_element_to_treat */
        number_candidate_element_to_treat=0;
        for (ulong i=0; i<max_candidate_element_to_treat; i++)
            candidate_element_to_treat[i]=0;
        
        for (ulong i=0; i<new_bamg_mesh_Ne; i++)
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
                neighbour=(ulong) new_bamg_mesh_ElementConnectivity[(tmp_born-1)*3+j_born];
                
                /*  if neighbour exists and is already detected in this cavity, do nothing */
                if(neighbour>0)
                {
                    if(cavity_number_born[neighbour-1]==new_cavity_number)
                        continue;
                }
                /*  if the edge exists in the old mesh */
                /*  the edge j goes from node j to node jnext */
                edge_node_1 = next(j_born);
                edge_node_2 =next(edge_node_1);      
                
                Previous_edge_node_1=PreviousNumbering[(ulong) new_bamg_mesh_Triangles[(tmp_born-1)*4+edge_node_1]-1];
                Previous_edge_node_2=PreviousNumbering[(ulong) new_bamg_mesh_Triangles[(tmp_born-1)*4+edge_node_2]-1];
                
                // _printf_("   Previous_edge_node_1 = " << Previous_edge_node_1 << "\n");
                // _printf_("   Previous_edge_node_2 = " << Previous_edge_node_2 << "\n");

                 if((Previous_edge_node_1>0) && (Previous_edge_node_2>0))
                    {        
                     /*  the last indices in the connectivity vector is the number */
                     /*  of connected nodes */
     
                    nb_connected_nodes=(ulong) bamg_mesh_NodalConnectivity[(Previous_edge_node_1-1)*bamg_mesh_NNodalConnectivity+bamg_mesh_NNodalConnectivity-1];
                     
                    found_intersect=0;
                    for (int i=0; i<nb_connected_nodes; i++)
                        if((ulong) bamg_mesh_NodalConnectivity[(Previous_edge_node_1-1)*bamg_mesh_NNodalConnectivity+i]==Previous_edge_node_2)
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
                            dead_element_node_1=(ulong) bamg_mesh_NodalElementConnectivity[(Previous_edge_node_1-1)*bamg_mesh_NNodalElementConnectivity+i_dead_element_node_1];
                            
                            if(dead_element_node_1==0)
                                break;
                            
                            for (i_dead_element_node_2=0; i_dead_element_node_2<bamg_mesh_NNodalElementConnectivity; i_dead_element_node_2++)
                            {
                                dead_element_node_2=(ulong) bamg_mesh_NodalElementConnectivity[(Previous_edge_node_2-1)*bamg_mesh_NNodalElementConnectivity+i_dead_element_node_2];
                                
                                if(dead_element_node_2==0)
                                    break;
                                
                                if(dead_element_node_1==dead_element_node_2)
                                {
                                    candidate_dead_element[nb_candidate_dead_element]=dead_element_node_1;
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
                               edge_node_1=next(j);
                               edge_node_2=next(edge_node_1);

                               tmp_node_1=(ulong) bamg_mesh_Triangles[(tmp_dead-1)*4+edge_node_1];
                               tmp_node_2=(ulong) bamg_mesh_Triangles[(tmp_dead-1)*4+edge_node_2];

                               if( (tmp_node_1==Previous_edge_node_1) && (tmp_node_2==Previous_edge_node_2))
                                break;
                            }
                            
                            if(j<3)
                            {
                                if(cavity_number_dead[tmp_dead-1]==-1)
                                    cavity_number_dead[tmp_dead-1]=new_cavity_number;
                                
                                boundary_between_two_cavities[(tmp_dead-1)*3+j]=1;
                                break;
                            }
                        }
                        continue;
                     } /* end of the if on found_intersect */
                 } /* end of the if */

                /*  if no neighbour, do nothing */
                if(neighbour==0) 
                    continue; 
                
                /*  if neighbour is not born, do nothing */
                if(cavity_number_born[neighbour-1]==0) 
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
    for (ulong i_dead=0; i_dead<nb_dead; i_dead++)
    {
        /*  first element to be checked */
        tmp_dead=ind_dead[i_dead];
        
        /*  We want to start from a dead element that has already be assigned to */
        /*  a cavity. */
        
        /*  if the dead element is not already treated, do nothing. It will be */
        /*  treated from another element */
        if(cavity_number_dead[tmp_dead-1]==-1)
            continue;
        
        new_cavity_number=cavity_number_dead[tmp_dead-1];
        
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
            
            cavity_number_dead[tmp_dead-1]=new_cavity_number;
            
            for (int j_dead=0; j_dead<3; j_dead++)
            {
                neighbour=(ulong) bamg_mesh_ElementConnectivity[(tmp_dead-1)*3+j_dead];
                
                /*  if no neighbour, do nothing */
                if(neighbour==0)
                    continue;
                
                /*  if neighbour is already treated, do nothing */
                if(cavity_number_dead[neighbour-1]!=-1)
                    continue;
                
                /*  if the edge corresponds to the boudary between two cavities */
                if(boundary_between_two_cavities[(tmp_dead-1)*3+j_dead]==1)
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

    ulong *size_born_cavity = xNew<ulong>(nb_cavities);
    ulong *size_dead_cavity = xNew<ulong>(nb_cavities);

    for (ulong i_born=0; i_born<nb_cavities; i_born++)
        size_born_cavity[i_born]=0;

    for (ulong i_dead=0; i_dead<nb_cavities; i_dead++)
        size_dead_cavity[i_dead]=0;


    if (verbosity>1) _printf_("   Detect_cavities: loop6...\n");
    for (ulong i_born=0; i_born<nb_born; i_born++)
    {
        new_cavity_number=cavity_number_born[ind_born[i_born]-1];
        
        if(new_cavity_number<0)
            _error_("A born element has not been assigned to a cavity");
        
        if(new_cavity_number>0)
            size_born_cavity[new_cavity_number-1]++;
    }

    bool not_all_dead_assigned=0;
    if (verbosity>1) _printf_("   Detect_cavities: loop7...\n");
    for (ulong i_dead=0; i_dead<nb_dead; i_dead++)
    {
        new_cavity_number=cavity_number_dead[ind_dead[i_dead]-1];
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
    for (ulong i=0; i<nb_cavities; i++)
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

    ulong *element_born_cavity = xNew<ulong>(nb_cavities*max_size_born_cavity);
    ulong *element_dead_cavity = xNew<ulong>(nb_cavities*max_size_dead_cavity);

    if (verbosity>1) _printf_("   Detect_cavities: loop9...\n");
    for (ulong i_born=0; i_born<nb_born; i_born++)
    {
        new_cavity_number=cavity_number_born[ind_born[i_born]-1];
        
        if(new_cavity_number>0)
        {
            element_born_cavity[(new_cavity_number-1)*max_size_born_cavity+size_born_cavity[new_cavity_number-1]]=ind_born[i_born];
            //element_born_cavity[nb_cavities*size_born_cavity[new_cavity_number-1]+new_cavity_number-1]=ind_born[i_born];
            size_born_cavity[new_cavity_number-1]++;
        }
    }

    if (verbosity>1) _printf_("   Detect_cavities: loop10...\n");
    for (ulong i_dead=0; i_dead<nb_dead; i_dead++)
    {
        new_cavity_number=cavity_number_dead[ind_dead[i_dead]-1];
        
        if(new_cavity_number>0)
        {
            element_dead_cavity[(new_cavity_number-1)*max_size_dead_cavity+size_dead_cavity[new_cavity_number-1]]=ind_dead[i_dead];
            //element_dead_cavity[nb_cavities*size_dead_cavity[new_cavity_number-1]+new_cavity_number-1]=ind_dead[i_dead];
            size_dead_cavity[new_cavity_number-1]++;
        }
    }

    if (verbosity>1) _printf_("   Detect_cavities: destroy final...\n");
    /* destroy the temporary matrix */
    delete ind_dead; 
    delete ind_born;
    delete boundary_between_two_cavities;
    delete is_born_element_to_treat;
    delete is_dead_element_to_treat;
    delete cavity_number_dead;
    delete cavity_number_born;

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
