/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   finiteelement.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Mon Aug 24 11:02:45 2015
 */

#include <meshtools.hpp>


namespace MeshTools
{


#if 0
std::vector<double>
rotatedWimElementsX(double const& rotangle)
{
    //get x coord of WIM centers (rotated)
    std::vector<double> x(num_elements_wim_grid);
    double cos_rotangle=std::cos(rotangle);
    double sin_rotangle=std::sin(rotangle);
    for (int i=0; i<num_elements_wim_grid; ++i)
    {
        x[i] = cos_rotangle*wim_grid.X[i] + sin_rotangle*wim_grid.Y[i];
    }

    return x;
}

std::vector<double>
rotatedWimElementsY(double const& rotangle)
{
    //get x coord of WIM centers (rotated)
    std::vector<double> y(num_elements_wim_grid);
    double cos_rotangle=std::cos(rotangle);
    double sin_rotangle=std::sin(rotangle);
    for (int i=0; i<num_elements_wim_grid; ++i)
    {
        y[i] = -sin_rotangle*wim_grid.X[i] + cos_rotangle*wim_grid.Y[i];
    }

    return y;
}
#endif


double
measure(double const x0,double const y0,double const x1,double const y1,double const x2,double const y2)
{
    double jac = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
    return .5*std::abs(jac);
}
double
jacobian(double const x0,double const y0,double const x1,double const y1,double const x2,double const y2)
{
    return (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
}

double
jacobian(element_type const& element, mesh_type const& mesh)
{
    std::vector<double> vertex_0 = mesh.nodes()[element.indices[0]-1].coords;//{x0,y0}
    std::vector<double> vertex_1 = mesh.nodes()[element.indices[1]-1].coords;//{x1,y1}
    std::vector<double> vertex_2 = mesh.nodes()[element.indices[2]-1].coords;//{x2,y2}

    //jacobian = determinant([x1-x0,y1-y0;x2-x0,y2-y0]);
    double jac = (vertex_1[0]-vertex_0[0])*(vertex_2[1]-vertex_0[1]);//(x1-x0)*(y2-y0)
    jac -= (vertex_2[0]-vertex_0[0])*(vertex_1[1]-vertex_0[1]);//(x2-x0)*(y1-y0)

    return  jac;
}

double
jacobian(element_type const& element, mesh_type const& mesh,
                        std::vector<double> const& um, double factor)
{
    std::vector<double> vertex_0 = mesh.nodes()[element.indices[0]-1].coords;
    std::vector<double> vertex_1 = mesh.nodes()[element.indices[1]-1].coords;
    std::vector<double> vertex_2 = mesh.nodes()[element.indices[2]-1].coords;
    int M_num_nodes = mesh.numNodes();

    for (int i=0; i<2; ++i)
    {
        vertex_0[i] += factor*um[element.indices[0]-1+i*(M_num_nodes)];
        vertex_1[i] += factor*um[element.indices[1]-1+i*(M_num_nodes)];
        vertex_2[i] += factor*um[element.indices[2]-1+i*(M_num_nodes)];
    }

    double jac = (vertex_1[0]-vertex_0[0])*(vertex_2[1]-vertex_0[1]);
    jac -= (vertex_2[0]-vertex_0[0])*(vertex_1[1]-vertex_0[1]);

    return  jac;
}

std::vector<double>
sides(element_type const& element, mesh_type const& mesh)
{
    std::vector<double> vertex_0 = mesh.nodes()[element.indices[0]-1].coords;
    std::vector<double> vertex_1 = mesh.nodes()[element.indices[1]-1].coords;
    std::vector<double> vertex_2 = mesh.nodes()[element.indices[2]-1].coords;

    std::vector<double> side(3);

    side[0] = std::hypot(vertex_1[0]-vertex_0[0], vertex_1[1]-vertex_0[1]);
    side[1] = std::hypot(vertex_2[0]-vertex_1[0], vertex_2[1]-vertex_1[1]);
    side[2] = std::hypot(vertex_2[0]-vertex_0[0], vertex_2[1]-vertex_0[1]);

    return side;
}

std::vector<double>
minMaxSide(mesh_type const& mesh)
{
    std::vector<double> minmax(2);
    std::vector<double> all_min_side(mesh.numTriangles());
    std::vector<double> all_max_side(mesh.numTriangles());

    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        auto side = MeshTools::sides(*it,mesh);
        all_min_side[cpt] = *std::min_element(side.begin(),side.end());
        all_max_side[cpt] = *std::max_element(side.begin(),side.end());
        ++cpt;
    }

    // minmax[0] = *std::min_element(all_min_side.begin(),all_min_side.end());
    // minmax[1] = *std::max_element(all_max_side.begin(),all_max_side.end());

    minmax[0] = std::accumulate(all_min_side.begin(),all_min_side.end(),0.)/(all_min_side.size());
    minmax[1] = std::accumulate(all_max_side.begin(),all_max_side.end(),0.)/(all_max_side.size());

    return minmax;
}

double
minAngles(element_type const& element, mesh_type const& mesh)
{
    std::vector<double> side = MeshTools::sides(element,mesh);
    //std::for_each(side.begin(), side.end(), [&](double& f){ f = 1000.*f; });
    std::sort(side.begin(),side.end());
    double minang = std::acos( (std::pow(side[1],2.) + std::pow(side[2],2.) - std::pow(side[0],2.) )/(2*side[1]*side[2]) );
    minang = minang*45.0/std::atan(1.0);

    return minang;
}

double
minAngle(mesh_type const& mesh)
{
    int M_num_elements = mesh.numTriangles();
    auto M_elements = mesh.triangles();
    std::vector<double> all_min_angle(M_num_elements);
    double min_angle;

#if 0
    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        all_min_angle[cpt] = MeshTools::minAngles(*it,mesh);
        ++cpt;
    }
#endif

#if 1
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        all_min_angle[cpt] = MeshTools::minAngles(M_elements[cpt],mesh);
    }
#endif

    min_angle = *std::min_element(all_min_angle.begin(),all_min_angle.end());
    return min_angle;
}//minAngle

double
minAngle(mesh_type const& mesh, std::vector<double> const& um, double factor)
{
    auto movedmesh = mesh;
    movedmesh.move(um,factor);

    int M_num_elements = movedmesh.numTriangles();
    std::vector<double> all_min_angle(movedmesh.numTriangles());

#if 0
    // int cpt = 0;
    // for (auto it=movedmesh.triangles().begin(), end=movedmesh.triangles().end(); it!=end; ++it)
    // {
    //     all_min_angle[cpt] = MeshTools::minAngles(*it,movedmesh);
    //     ++cpt;
    // }
#endif

#if 1
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        all_min_angle[cpt] = MeshTools::minAngles(movedmesh.triangles()[cpt],movedmesh);
    }
#endif

    return *std::min_element(all_min_angle.begin(),all_min_angle.end());;
}//minAngle

bool
flip(mesh_type const& mesh, std::vector<double> const& um, double factor)
{
    auto movedmesh = mesh;
    movedmesh.move(um,factor);

    std::vector<double> area(movedmesh.numTriangles());
    double area_init;

 #if 0
    int cpt = 0;
    for (auto it=movedmesh.triangles().begin(), end=movedmesh.triangles().end(); it!=end; ++it)
    {
        area[cpt] = MeshTools::jacobian(*it,movedmesh);
        ++cpt;
    }
#endif

#if 1
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
    int M_num_elements = mesh.numTriangles();

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        area_init = MeshTools::jacobian(mesh.triangles()[cpt],mesh);
        area[cpt] = MeshTools::jacobian(movedmesh.triangles()[cpt],movedmesh);

        if(area_init*area[cpt]<=0.)
        {
            std::cout <<"FLIP DETECTED element:"<< cpt <<"\n";
        }
    }
#endif

    double minarea = *std::min_element(area.begin(),area.end());
    double maxarea = *std::max_element(area.begin(),area.end());

    return ((minarea <= 0.) && (maxarea >= 0.));
}//flip

double
resolution(mesh_type const& mesh)
{
    std::vector<double> all_min_measure(mesh.numTriangles());

    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        all_min_measure[cpt] = MeshTools::measure(*it,mesh);
        ++cpt;
    }

    double resol = std::accumulate(all_min_measure.begin(),all_min_measure.end(),0.)/(all_min_measure.size());
    resol = std::pow(resol,0.5);//sqrt of the mean area of the triangles

    return resol;
}//resolution


std::vector<double>
hminVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh)
{
    std::vector<double> hmin(bamg_mesh->NodalElementConnectivitySize[0]);

    for (int i=0; i<bamg_mesh->NodalElementConnectivitySize[0]; ++i)
    {
        std::vector<double> measure(bamg_mesh->NodalElementConnectivitySize[1]);
        int j = 0;
        for (j=0; j<bamg_mesh->NodalElementConnectivitySize[1]; ++j)
        {
            int elt_num = bamg_mesh->NodalElementConnectivity[bamg_mesh->NodalElementConnectivitySize[1]*i+j]-1;

            if ((0 <= elt_num) && (elt_num < mesh.numTriangles()) && (elt_num != NAN))
            {
                measure[j] = MeshTools::measure(mesh.triangles()[elt_num],mesh);
            }
            else
            {
                break;
            }
        }

        measure.resize(j);
        hmin[i] = std::sqrt(2.)*std::sqrt(*std::min_element(measure.begin(),measure.end()))*0.8;
    }

    return hmin;
}//hminVertices

std::vector<double>
hmaxVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh)
{
    std::vector<double> hmax(bamg_mesh->NodalElementConnectivitySize[0]);

    for (int i=0; i<bamg_mesh->NodalElementConnectivitySize[0]; ++i)
    {
        std::vector<double> measure(bamg_mesh->NodalElementConnectivitySize[1]);
        int j = 0;
        for (j=0; j<bamg_mesh->NodalElementConnectivitySize[1]; ++j)
        {
            int elt_num = bamg_mesh->NodalElementConnectivity[bamg_mesh->NodalElementConnectivitySize[1]*i+j]-1;

            if ((0 <= elt_num) && (elt_num < mesh.numTriangles()) && (elt_num != NAN))
            {
                measure[j] = MeshTools::measure(mesh.triangles()[elt_num],mesh);
            }
            else
            {
                break;
            }
        }

        measure.resize(j);
        hmax[i] = std::sqrt(2.)*std::sqrt(*std::max_element(measure.begin(),measure.end()))*1.2;
    }

    return hmax;
}//hmaxVertices

std::vector<double>
AllMinAngle(mesh_type const& mesh, std::vector<double> const& um, double factor)
{
    auto movedmesh = mesh;
    movedmesh.move(um,factor);

    std::vector<double> all_min_angle(movedmesh.numTriangles());

#if 0
    // int cpt = 0;
    // for (auto it=movedmesh.triangles().begin(), end=movedmesh.triangles().end(); it!=end; ++it)
    // {
    //     all_min_angle[cpt] = MeshTools::minAngles(*it,movedmesh);
    //     ++cpt;
    // }
#endif

#if 1
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
    int M_num_elements = mesh.numTriangles();

//#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        all_min_angle[cpt] = MeshTools::minAngles(movedmesh.triangles()[cpt],movedmesh);
    }
#endif

    return all_min_angle;
}//AllMinAngle

double
measure(element_type const& element, mesh_type const& mesh)
{
    return (1./2)*std::abs(jacobian(element,mesh));
}//measure

double
measure(element_type const& element, mesh_type const& mesh,
                       std::vector<double> const& um, double factor)
{
    return (1./2)*std::abs(jacobian(element,mesh,um,factor));
}//measure


std::vector<double>
shapeCoeff(element_type const& element, mesh_type const& mesh)
{
    std::vector<double> x(3);
    std::vector<double> y(3);

    for (int i=0; i<3; ++i)
    {
        x[i] = mesh.nodes()[element.indices[i]-1].coords[0];
        y[i] = mesh.nodes()[element.indices[i]-1].coords[1];
    }

    std::vector<double> coeff(6);
    double jac = jacobian(element,mesh);

    for (int k=0; k<6; ++k)
    {
        int kp1 = (k+1)%3;
        int kp2 = (k+2)%3;

        if (k<3)
        {
            coeff[k] = (y[kp1]-y[kp2])/jac;
        }
        else
        {
            coeff[k] = (x[kp2]-x[kp1])/jac;
        }
    }

    return coeff;
}//shapeCoeff


#if 1
void
advect(double** interp_elt_out_ptr,     // pointer to pointer to output data
        double* interp_elt_in,          // pointer to input data
        mesh_info_type_dbl* mesh_info,  // pointer to structure with mesh info: positions of nodes and elements,
                                        //  index (maps elements to nodes), element connectivity
        double* VC_in,                  // pointer to convective velocities (len = 2*num_nodes)
        int* interp_method,             // pointer to interp methods for each variable
        int nb_var,                     // number of variables
        double time_step)               // time step (s)
{

	/*Initialize output*/
    int Nels = mesh_info->num_elements;
    int Nnod = mesh_info->num_nodes;
	double* interp_elt_out=NULL;

    interp_elt_out=xNew<double>(nb_var*Nels);

    int thread_id;
    int total_threads;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < Nels; ++cpt)
    {
        /* some variables used for the advection*/
        double surface  = mesh_info->surface[cpt];
        double integrated_variable;

        /* some variables used for the advection*/
        double x[3],y[3],x_new[3],y_new[3];
        int x_ind, y_ind, neighbour_int, vertex_1, vertex_2;
        double outer_fluxes_area[3], vector_edge[2], outer_vector[2], VC_middle[2], VC_x[3], VC_y[3];
        int fluxes_source_id[3];

        double neighbour_double;
        int other_vertex[3*2]={1,2 , 2,0 , 0,1};


        /* convective velocity */
        for(int i=0;i<3;i++)
        {
            x_ind   = mesh_info->index[3*cpt+i]-1;//NB bamg index starts at 1
            y_ind   = x_ind+Nnod;

            //positions of nodes
            x[i] = mesh_info->nodes_x[x_ind];
            y[i] = mesh_info->nodes_y[x_ind];

            //convective velocity
            VC_x[i] = VC_in[x_ind];
            VC_y[i] = VC_in[y_ind];
        }

        //fluxes
        for(int i=0;i<3;i++)
        {
            outer_fluxes_area[i]=0;

            vertex_1=other_vertex[2*i  ];
            vertex_2=other_vertex[2*i+1];

            vector_edge[0]=x[vertex_2]-x[vertex_1];
            vector_edge[1]=y[vertex_2]-y[vertex_1];

            outer_vector[0]= vector_edge[1];
            outer_vector[1]=-vector_edge[0];

            VC_middle[0] = (VC_x[vertex_2]+VC_x[vertex_1])/2.;
            VC_middle[1] = (VC_y[vertex_2]+VC_y[vertex_1])/2.;

            outer_fluxes_area[i]=outer_vector[0]*VC_middle[0]+outer_vector[1]*VC_middle[1];


            if(outer_fluxes_area[i]>0)
            {
                outer_fluxes_area[i]=std::min(surface/time_step/3.,outer_fluxes_area[i]);
                fluxes_source_id[i]=cpt;
            }
            else
            {
			    neighbour_double = mesh_info->element_connectivity[cpt*3+i] -1;//NB bamg indices start at 1
                neighbour_int=(int) neighbour_double;
			    if (!std::isnan(neighbour_double) && neighbour_int>0)
                {
                    double surface = mesh_info->surface[neighbour_int];//NB don't want to reset "surface" outside this scope
                    outer_fluxes_area[i]=-std::min(surface/time_step/3.,-outer_fluxes_area[i]);
                    fluxes_source_id[i]=neighbour_int;
                }
                else // open boundary with incoming fluxes
                    fluxes_source_id[i]=cpt;
            }
        }


        for(int j=0; j<nb_var; j++)
        {
            if(interp_method[j]==1)
            {
                double tmp = 0.;
                for (int k=0;k<3;k++)
                    tmp += interp_elt_in[fluxes_source_id[k]*nb_var+j]*outer_fluxes_area[k];
                interp_elt_out[cpt*nb_var+j] = interp_elt_in[cpt*nb_var+j] - (tmp/surface)*time_step;
            }
            else
                interp_elt_out[cpt*nb_var+j] = interp_elt_in[cpt*nb_var+j];
        }
    }
	*interp_elt_out_ptr=interp_elt_out;
}//advect
#endif

#if 0
void
diffuse(double* variable_elt, double diffusivity_parameters, double dx)
{
    if(diffusivity_parameters<=0.)
    {
        LOG(DEBUG) <<"diffusivity parameter lower or equal to 0 \n";
        LOG(DEBUG) <<"nothing to do\n";
        return;
    }   

    double factor=diffusivity_parameters*time_step/std::pow(dx,2.);
    double* old_variable_elt=xNew<double>(M_num_elements);
    
    for (int cpt=0; cpt < M_num_elements; ++cpt)
        old_variable_elt[cpt]=variable_elt[cpt];

    int thread_id;
    int total_threads;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    int Nd = bamgmesh->NodalConnectivitySize[1];

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        /* some variables used for the advection*/
        int neighbour_int;
        double fluxes_source[3];
        int fluxes_source_id;

        double neighbour_double;
      
        for(int i=0;i<3;i++)
        {
            neighbour_double=bamgmesh->ElementConnectivity[cpt*3+i];
            neighbour_int=(int) bamgmesh->ElementConnectivity[cpt*3+i];
                    
            if (!std::isnan(neighbour_double) && neighbour_int>0)
            {
                fluxes_source_id=neighbour_int-1;
                fluxes_source[i]=factor*(old_variable_elt[fluxes_source_id]-old_variable_elt[cpt]);
            }
            else // no diffusion crosses open nor closed boundaries
                 fluxes_source[i]=0.;        
        }
        variable_elt[cpt] += fluxes_source[0] + fluxes_source[1] + fluxes_source[2];                
    }
    // Cleaning
    xDelete<double>(old_variable_elt);
}//diffuse
#endif


#if 0
void
adaptMesh()
{
    delete bamgopt_previous;
    delete bamggeom_previous;
    delete bamgmesh_previous;

    bamgopt_previous = new BamgOpts();
    bamggeom_previous = new BamgGeom();
    bamgmesh_previous = new BamgMesh();

    *bamgmesh_previous = *bamgmesh;
    *bamggeom_previous = *bamggeom;
    *bamgopt_previous = *bamgopt;

    int fnd = 0;

    for (int edg=0; edg<bamgmesh_previous->EdgesSize[0]; ++edg)
    {
        fnd = bamgmesh_previous->Edges[3*edg]-1;

        if ((std::binary_search(M_dirichlet_flags.begin(),M_dirichlet_flags.end(),fnd)))
        {
            bamggeom_previous->Edges[3*edg+2] = M_flag_fix;
            bamgmesh_previous->Edges[3*edg+2] = M_flag_fix;
        }
        else
        {
            bamggeom_previous->Edges[3*edg+2] = M_flag_fix+1; // we just want it to be different than M_flag_fix 
            bamgmesh_previous->Edges[3*edg+2] = M_flag_fix+1; // we just want it to be different than M_flag_fix 
        }
        
    }

    //Environment::logMemoryUsage("before adaptMesh");
    Bamgx(bamgmesh,bamggeom,bamgmesh_previous,bamggeom_previous,bamgopt_previous);
    //Environment::logMemoryUsage("after adaptMesh");

    // Save the old id_node before redefining it
    std::vector<int> old_node_id=M_mesh.id();

    // Import the mesh from bamg
    MeshTools::importBamg(bamgmesh);

    // We mask out the boundary nodes
    M_mask.assign(bamgmesh->VerticesSize[0],false) ;
    M_mask_dirichlet.assign(bamgmesh->VerticesSize[0],false) ;

    for (int vert=0; vert<bamgmesh->VerticesOnGeomVertexSize[0]; ++vert)
        M_mask[bamgmesh->VerticesOnGeomVertex[2*vert]-1]=true; // The factor 2 is because VerticesOnGeomVertex has 2 dimensions in bamg

        // Recompute the node ids
    if(bamgopt->KeepVertices)
    {
        std::vector<int> new_nodes_id=M_mesh.id();

        int Boundary_id=0;
        int nb_new_nodes=0;

        // The new id will have values higher than the previous one
        int first_new_node=*std::max_element(old_node_id.begin(),old_node_id.end())+1;

        for (int vert=0; vert<bamgmesh->VerticesSize[0]; ++vert)
        {
            if(M_mask[vert])
            {
                Boundary_id++;
                new_nodes_id[vert]=Boundary_id;
            }
            else
            {
                if(bamgmesh->PreviousNumbering[vert]==0)
                {
                    nb_new_nodes++;
                    new_nodes_id[vert]=first_new_node+nb_new_nodes-1;
                }
                else
                {
                    new_nodes_id[vert]=old_node_id[bamgmesh->PreviousNumbering[vert]-1];
                }
            }
        }
        M_mesh.setId(new_nodes_id);
    }


    // update dirichlet nodes
    M_boundary_flags.resize(0);
    M_dirichlet_flags.resize(0);
    for (int edg=0; edg<bamgmesh->EdgesSize[0]; ++edg)
    {
        M_boundary_flags.push_back(bamgmesh->Edges[3*edg]-1);
        //M_boundary_flags.push_back(bamgmesh->Edges[3*edg+1]-1);

        if (bamgmesh->Edges[3*edg+2] == M_flag_fix)
        {
            M_dirichlet_flags.push_back(bamgmesh->Edges[3*edg]-1);
            //M_dirichlet_flags.push_back(bamgmesh->Edges[3*edg+1]-1);

            //std::cout<<"NODES["<< edg <<"]= "<< bamgmesh->Edges[3*edg]-1 << " and "<< bamgmesh->Edges[3*edg+1]-1 <<"\n";
        }
    }

    std::sort(M_dirichlet_flags.begin(), M_dirichlet_flags.end());
    //M_dirichlet_flags.erase( std::unique(M_dirichlet_flags.begin(), M_dirichlet_flags.end() ), M_dirichlet_flags.end());

    std::sort(M_boundary_flags.begin(), M_boundary_flags.end());
    //M_boundary_flags.erase(std::unique(M_boundary_flags.begin(), M_boundary_flags.end() ), M_boundary_flags.end());

    M_neumann_flags.resize(0);
    std::set_difference(M_boundary_flags.begin(), M_boundary_flags.end(),
                        M_dirichlet_flags.begin(), M_dirichlet_flags.end(),
                        std::back_inserter(M_neumann_flags));

    // for (const int& edg : M_dirichlet_flags)
    // {
    //     std::cout<<"AFTER["<< edg << "]= ("<< M_mesh.coordX()[edg] <<","<< M_mesh.coordY()[edg] <<")\n";
    // }

    M_dirichlet_nodes.resize(2*(M_dirichlet_flags.size()));
    for (int i=0; i<M_dirichlet_flags.size(); ++i)
    {
        M_dirichlet_nodes[2*i] = M_dirichlet_flags[i];
        M_dirichlet_nodes[2*i+1] = M_dirichlet_flags[i]+M_num_nodes;
        M_mask_dirichlet[M_dirichlet_flags[i]]=true;
    }


    M_neumann_nodes.resize(2*(M_neumann_flags.size()));
    for (int i=0; i<M_neumann_flags.size(); ++i)
    {
        M_neumann_nodes[2*i] = M_neumann_flags[i];
        M_neumann_nodes[2*i+1] = M_neumann_flags[i]+M_num_nodes;
    }

    M_surface.assign(M_num_elements,0.);

    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        M_surface[cpt] = MeshTools::measure(*it,M_mesh);
        ++cpt;
    }
}//adaptMesh
#endif


#if 0
void
coriolis()
{
    // Interpolation of the latitude
    std::vector<double> lat = M_mesh.meanLat();

    for (int i=0; i<M_fcor.size(); ++i)
    {
        if (vm["simul.use_coriolis"].as<bool>())
        {
            M_fcor[i] = 2*(physical::omega)*std::sin(lat[i]*PI/180.);
        }
        else
        {
            M_fcor[i] = 0.;
        }
    }
}//coriolis
#endif


#if 0
void
nodesToElements(double const* depth, std::vector<double>& v)
{
    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        double sum = 0;
        for (int j=0; j<3; ++j)
        {
            sum += depth[it->indices[j]-1];
        }

        v[cpt] = sum/3.0;
        ++cpt;
    }
}//nodesToElements
#endif


std::string
gitRevision()
{
    //std::string command = "git rev-parse HEAD";
    return MeshTools::system("git rev-parse HEAD");
}//gitRevision

std::string
system(std::string const& command)
{
    char buffer[128];
    //std::string command = "git rev-parse HEAD";
    std::string result = "";
    std::shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe)
    {
        throw std::runtime_error("popen() failed!");
    }

    while (!feof(pipe.get()))
    {
        if (fgets(buffer, 128, pipe.get()) != NULL)
        {
            // remove newline from the buffer
            int len = strlen(buffer);
            if( buffer[len-1] == '\n' )
                buffer[len-1] = 0;

            result += buffer;
        }
    }

    // return the result of the command
    return result;
}//system

std::string
getEnv(std::string const& envname)
{
    const char* senv = ::getenv(envname.c_str());
    if ( senv == NULL )
        senv = "NULL";
    return std::string(senv);
}//getEnv

} // namespace MeshTools
