/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   finiteelement.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Mon Aug 24 11:02:45 2015
 */

#include <meshtools.hpp>
#include <BamgTriangulatex.h>

namespace Wim
{

// ===============================================================
// MeshInfo class

template<typename T>
MeshInfo<T>::MeshInfo(T_gmsh const &movedmesh)
{
    M_mesh_type = E_mesh_type::simple;

    //update mesh with moved mesh
    this->setFields(movedmesh);

    // set some other things (eg resolution)
    this->initSimple();
}//MeshInfo constructor (simple)


template<typename T>
MeshInfo<T>::MeshInfo(std::vector<int> const& index,T_val_vec const &nodes_x,T_val_vec const &nodes_y)
{
    //set some scalars
    M_mesh_type     = E_mesh_type::simple;
    M_num_nodes     = nodes_x.size();
    M_num_elements  = (index.size())/3;

    //set some vectors
    M_index   = index;
    M_nodes_x = nodes_x;
    M_nodes_y = nodes_y;

    // set some other things (eg resolution)
    this->initSimple();

}//MeshInfo constructor (simple)


template<typename T>
MeshInfo<T>::MeshInfo(T_val_vec const &nodes_x,T_val_vec const &nodes_y)
{
    //takes nodes and does triangulation

    //set some scalars
    M_mesh_type     = E_mesh_type::simple;
    M_num_nodes     = nodes_x.size();

    //set some vectors
    M_nodes_x = nodes_x;
    M_nodes_y = nodes_y;

    //do triangulation
    std::cout<<"#nodes for triangulation = "<<M_num_nodes<<"\n";

    int* index;
    BamgTriangulatex(&index,          //pointer to index              (output)
                     &M_num_elements, //pointer to num elements       (output)
                     &(M_nodes_x)[0], //pointer to x-coord of nodes   (input)
                     &(M_nodes_y)[0], //pointer to y-coord of nodes   (input)
                     M_num_nodes);    //num nodes                     (input)

    std::cout<<"#elements from triangulation = "<<M_num_elements<<"\n";

    //get M_index from the pointer index
    M_index.resize(3*M_num_elements);
    for (int i=0; i<3*M_num_elements; i++)
        M_index[i]   = index[i];

    xDelete<int>(index);

    // set some other things (eg resolution)
    this->initSimple();

}//MeshInfo constructor (simple)


template<typename T>
MeshInfo<T>::MeshInfo(std::vector<int> const& index,T_val_vec const &nodes_x,T_val_vec const &nodes_y,
        T_val_vec const &elements_x,T_val_vec const &elements_y)
{
    //set some scalars
    M_mesh_type     = E_mesh_type::simple;
    M_num_nodes     = nodes_x.size();
    M_num_elements  = (index.size())/3;

    //set some vectors
    M_index      = index;
    M_nodes_x    = nodes_x;
    M_nodes_y    = nodes_y;
    M_elements_x = elements_x;
    M_elements_y = elements_y;

    // set some other things (eg resolution)
    this->initSimple();

}//MeshInfo constructor (simple)


template<typename T>
MeshInfo<T>::MeshInfo(T_gmsh const &movedmesh,BamgMesh* bamgmesh,int const& flag_fix)
{
    //interface for M_wim_on_mesh
    M_mesh_type = E_mesh_type::full;

    //set the basic fields
    this->setFields(movedmesh);


    //calculate the surface areas,
    //get the element connectivity from bamgmesh
    int Nels = M_num_elements;
    M_surface.assign(Nels,0);
    M_element_connectivity.assign(3*Nels,0);
    for (int i=0;i<Nels;i++)
    {
        T_val_vec xnods(3);
        T_val_vec ynods(3);
        for (int k=0;k<3;k++)
        {
            int ind  = M_index[3*i+k]-1;//NB bamg indices go from 1 to Nels
            xnods[k] = M_nodes_x[ind];
            ynods[k] = M_nodes_y[ind];

            M_element_connectivity[3*i+k]
                = bamgmesh->ElementConnectivity[3*i+k];//NB stick to bamg convention (indices go from 1 to Nels)
        }

        T_val area = .5*MeshTools::jacobian(
                xnods[0],ynods[0],xnods[1],ynods[1],xnods[2],ynods[2]);
        if(area>=0.)
            M_surface[i] = area;
        else
        {
            std::cout<<"Area of triangle "<<i<<" <0 : "<<area<<"\n";
            throw std::runtime_error("MeshInfo<T>::MeshInfo (full): negative area found\n");
        }
    }

    // ================================================================
    //get the Dirichlet mask
    int Nn  = M_num_nodes;
    std::vector<int> dirichlet_flags(0);
    for (int edg=0; edg<bamgmesh->EdgesSize[0]; ++edg)
        if (bamgmesh->Edges[3*edg+2] == flag_fix)
            dirichlet_flags.push_back(bamgmesh->Edges[3*edg]-1);
    
    M_mask_dirichlet.assign(Nn,false);
    for (int i=0; i<dirichlet_flags.size(); ++i)
        M_mask_dirichlet[dirichlet_flags[i]] = true;
    // ================================================================

    int max_nec = bamgmesh->NodalElementConnectivitySize[1];
    M_max_node_el_conn = max_nec;
    M_node_element_connectivity.resize(Nn*max_nec);
    for (int i=0;i<Nn;i++)
        for (int j=0; j<max_nec; ++j)
        {
            M_node_element_connectivity[max_nec*i+j]
                = bamgmesh->NodalElementConnectivity[max_nec*i+j];
            // NB stick to bamg convention (element indices go from 1 to Nels)
            // To test if element is OK:
            // elt_num  = M_node_element_connectivity[max_nec*i+j]-1;
            // OK if ((0 <= elt_num) && (elt_num < mesh.numTriangles()) && (elt_num != NAN))
        }
    // ================================================================================

    M_resolution = MeshTools::resolution(movedmesh);
}//MeshInfo constructor (full)


template<typename T>
void MeshInfo<T>::initSimple()
{
    if(M_elements_x.size()==0)
        this->setElements();

    int Nels = M_num_elements;
    T_val area_tot = 0.;
    for (int i=0;i<Nels;i++)
    {
        T_val_vec xnods(3);
        T_val_vec ynods(3);
        for (int k=0;k<3;k++)
        {
            int ind  = M_index[3*i+k]-1;//NB bamg indices go from 1 to Nels
            xnods[k] = M_nodes_x[ind];
            ynods[k] = M_nodes_y[ind];
        }

        area_tot += .5*MeshTools::jacobian(
                xnods[0],ynods[0],xnods[1],ynods[1],xnods[2],ynods[2]);
    }

    M_resolution = std::sqrt(area_tot/Nels);//sqrt of the mean area of the triangles
}

template<typename T>
void MeshInfo<T>::setElements()
{
    M_elements_x.assign(M_num_elements,0.);
    M_elements_y.assign(M_num_elements,0.);

    //take the average of the nodal positions
    T_val wt = 1./3.;
    for (int i=0;i<M_num_elements;i++)
        for (int k=0;k<3;k++)
        {
            int ind  = M_index[3*i+k]-1;//NB bamg indices go from 1 to Nels
            M_elements_x[i] += wt*M_nodes_x[ind];
            M_elements_y[i] += wt*M_nodes_y[ind];
        }
}//setElements()


template<typename T>
void MeshInfo<T>::setFields(T_gmsh const &mesh_in)
{
    //sets the basic mesh fields
    M_num_nodes     = mesh_in.numNodes();
    M_num_elements  = mesh_in.numTriangles();
    M_index         = mesh_in.indexTr();
    M_id            = mesh_in.id();
    M_nodes_x       = mesh_in.coordX();
    M_nodes_y       = mesh_in.coordY();
    M_elements_x    = mesh_in.bcoordX();
    M_elements_y    = mesh_in.bcoordY();
}//setFields


template<typename T>
void MeshInfo<T>::advect(T_val** interp_elt_out_ptr, // pointer to pointer to output data
        T_val* interp_elt_in,                        // pointer to input data
        T_val* VC_in,                                // pointer to convective velocities (len = 2*num_nodes)
        int* interp_method,                           // pointer to interp methods for each variable
        int nb_var,                                   // number of variables
        T_val time_step,                             // time step (s)
        T_val* inc_values)                           // values incoming from open boundaries
                                                      // - currently just specify a constant (ie spatially invariant) value for each variable
{

    if(M_mesh_type!=E_mesh_type::full)
    {
        std::cout<<"M_mesh_type (E_mesh_type::full) = "<<(int) M_mesh_type<<" ("<< (int) E_mesh_type::full<<")\n";
        throw runtime_error("advect: Mesh object not fully initialised\n");
    }

    //general variables
    int Nels = M_num_elements;
    int Nnod = M_num_nodes;
    int thread_id;
    int total_threads;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

	/*Initialize output*/
	T_val* interp_elt_out   = NULL;
    //interp_elt_out          = new T_val(nb_var*Nels);
    interp_elt_out          = xNew<T_val>(nb_var*Nels);


#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < Nels; ++cpt)
    {
        //std::cout<<"advect: cpt,Nels = "<<cpt<<","<<Nels<<"\n";
        /* some variables used for the advection*/
        T_val surface  = M_surface[cpt];

        /* some variables used for the advection*/
        T_val x[3],y[3];
        T_val outer_fluxes_area[3], vector_edge[2], outer_vector[2], VC_middle[2], VC_x[3], VC_y[3];
        int fluxes_source_id[3], node_nums[3];
        int other_vertex[3*2]={1,2 , 2,0 , 0,1};


        /* convective velocity */
        for(int i=0;i<3;i++)
        {
            node_nums[i] = M_index[3*cpt+i]-1;//NB bamg index starts at 1

            //positions of nodes
            x[i] = M_nodes_x[node_nums[i]];
            y[i] = M_nodes_y[node_nums[i]];

            //convective velocity
            VC_x[i] = VC_in[node_nums[i]];
            VC_y[i] = VC_in[node_nums[i]+Nnod];
        }

        //fluxes
        //std::cout<<"fluxes\n";
        for(int i=0;i<3;i++)
        {
            //std::cout<<"node = "<<i<<"\n";
            outer_fluxes_area[i]=0;

            int vertex_1 = other_vertex[2*i  ];
            int vertex_2 = other_vertex[2*i+1];

            int node_num1   = node_nums[vertex_1];
            int node_num2   = node_nums[vertex_2];

            vector_edge[0]=x[vertex_2]-x[vertex_1];
            vector_edge[1]=y[vertex_2]-y[vertex_1];

            outer_vector[0]= vector_edge[1];
            outer_vector[1]=-vector_edge[0];

            VC_middle[0] = (VC_x[vertex_2]+VC_x[vertex_1])/2.;
            VC_middle[1] = (VC_y[vertex_2]+VC_y[vertex_1])/2.;

            outer_fluxes_area[i]=outer_vector[0]*VC_middle[0]+outer_vector[1]*VC_middle[1];

            if(outer_fluxes_area[i]>0)
            {
                // std::cout<<"outward flux\n";
                // fluxes are going out of the cell
                // - this should let the quantity leave the cell without problems, even on boundaries
                // - we treat open and closed boundaries the same
                outer_fluxes_area[i]=std::min(surface/time_step/3.,outer_fluxes_area[i]);//limit the flux
                fluxes_source_id[i]=cpt;
            }//outgoing fluxes
            else
            {
                //std::cout<<"inward flux\n";
                // fluxes are coming into the cell
                // - find which element it is coming from
			    T_val neighbour_double = M_element_connectivity[cpt*3+i] -1;//NB bamg indices start at 1
                int neighbour_int=(int) neighbour_double;
                //std::cout<<"neighbour = "<<neighbour_double<<","<<neighbour_int<<"\n";
			    if (!std::isnan(neighbour_double) && neighbour_int>=0 && neighbour_int<Nels)
                {
                    // not on a boundary
                    //std::cout<<"neighbour = "<<neighbour_double<<","<<neighbour_int<<"\n";
                    T_val surface2 = M_surface[neighbour_int];
                    outer_fluxes_area[i]=-std::min(surface2/time_step/3.,-outer_fluxes_area[i]);//limit the flux
                    fluxes_source_id[i]=neighbour_int;
                }//not on boundary
                else
                {
                    // on a boundary with incoming fluxes
                    if (M_mask_dirichlet[node_num1]&&M_mask_dirichlet[node_num2])
                    {
                        // closed boundary (on coast) - zero flux coming in
                        fluxes_source_id[i] = cpt;
                        outer_fluxes_area[i] = 0.;
                    }
                    else
                    {
                        // open boundary - zero flux coming in
                        fluxes_source_id[i] = -1;//flag to use inc_values
                    }
                }//on boundary
            }//incoming fluxes
        }// loop over nodes (fluxes)


        //loop over variables
        //std::cout<<"variables\n";
        for(int j=0; j<nb_var; j++)
        {
            if(interp_method[j]==1)
            {
                T_val tmp = 0.;
                for (int i=0;i<3;i++)
                    if(fluxes_source_id[i]>=0)
                        tmp += interp_elt_in[fluxes_source_id[i]*nb_var+j]*outer_fluxes_area[i];
                    else
                        tmp += inc_values[j]*outer_fluxes_area[i];
                interp_elt_out[cpt*nb_var+j] = interp_elt_in[cpt*nb_var+j] - (tmp/surface)*time_step;
            }
            else
                interp_elt_out[cpt*nb_var+j] = interp_elt_in[cpt*nb_var+j];
        }
    }
	*interp_elt_out_ptr=interp_elt_out;
}//advect


template<typename T>
void MeshInfo<T>::interpToPoints(
        T_val_vec_ptrs &output_data,         //output data
        T_val_vec_ptrs const &input_data,    //input data
        T_val_vec &Rx,                       //location of output data (x-coord)
        T_val_vec &Ry,                       //location of output data (y-coord)
        std::vector<int> const& wet_indices) //optional (say which cells in input_data are wet)
{

    if(M_mesh_type==E_mesh_type::uninitialised)
        throw runtime_error("interpToPoints: Mesh object uninitialised\n");

    int nb_var      = input_data.size();
    if (nb_var==0)
        //nothing to do
        return;

    int Ninterp = (*(input_data[0])).size();//get pointer, then get size
    T_val_vec interp_in(Ninterp*nb_var);    //input to interp routine

    if (wet_indices.size()==0)
    {
        std::cout<<"interpToPoints: Ninterp,nb_var = "<<Ninterp<<","<<nb_var<<"\n";
        for (int i=0;i<Ninterp;i++)
            for (int p=0;p<nb_var;p++)
                interp_in[nb_var*i+p]   = (*(input_data[p]))[i];
    }
    else
    {
        int Ninterp = wet_indices.size();//get pointer, then get size
        std::cout<<"interpToPoints: Ninterp (wet),nb_var = "<<Ninterp<<","<<nb_var<<"\n";

        interp_in.resize(Ninterp*nb_var);   //input to interp routine
        for (int i=0;i<Ninterp;i++)
            for (int p=0;p<nb_var;p++)
                interp_in[nb_var*i+p]   = (*(input_data[p]))[wet_indices[i]];
    }

    T_val* interp_out;
    int target_size = Rx.size();

    std::cout<<"interpToPoints: InterpFromMeshToMesh2dx\n";
    InterpFromMeshToMesh2dx(&interp_out,  // output data
                          &M_index[0],    // index 
                          &M_nodes_x[0],  // location of input nodes (x-coord)
                          &M_nodes_y[0],  // location of input nodes (y-coord)
                          M_num_nodes,    // num nodes
                          M_num_elements, // num elements
                          &interp_in[0],  // input data
                          Ninterp,        // num input locations
                          nb_var,         // num input variables
                          &Rx[0],         // output location (x-coord)
                          &Ry[0],         // output location (y-coord)
                          target_size,    // num output locations
                          false);         // use default value if outside mesh (use nearest)
    

    //output
    for (int p=0;p<nb_var;p++)
    {
        output_data[p]->assign(target_size,0.);
        for (int i=0;i<target_size;i++)
            (*(output_data[p]))[i]  = interp_out[nb_var*i+p];
    }

    xDelete<T_val>(interp_out);

}//interpToPoints


template<typename T>
void MeshInfo<T>::elementsToNodes(
        T_val_vec_ptrs &output_data,       //output data
        T_val_vec_ptrs const &input_data)  //input data
{

    if(M_mesh_type==E_mesh_type::uninitialised)
        throw runtime_error("elementsToNodes: Mesh object uninitialised\n");

    if(0)
    {
        // just use meshToPoints()
        // - currently crashes
        this->interpToPoints(output_data,input_data,M_nodes_x,M_nodes_y);
        return;
    }

    // meshToPoints method currently crashes
    // - here we take nodal value to be the average of connected elements
    // TODO just do this for the boundary nodes, and interp to interior?
    // TODO better way?
    int Nn      = M_num_nodes;
    for (auto it=output_data.begin();it!=output_data.end();it++)
        (*it)->assign(Nn,0.);

    int max_nec = M_max_node_el_conn;
    int nb_var  = input_data.size();
    for (int i=0;i<Nn;i++)
    {
        int nec = 0;

        // check which elements are connected to the nodes
        // and accumulate the values
        // (take average once we know the number of good values)
        for (int j=0; j<max_nec; ++j)
        {
            int elt_num  = M_node_element_connectivity[max_nec*i+j]-1;//NB bamg indices start at 1
            if ((0 <= elt_num) && (elt_num < M_num_elements) && (elt_num != NAN))
            {
                nec++;//it's a good element
                for(int tmp_nb_var=0;tmp_nb_var<nb_var;tmp_nb_var++)
                    (*(output_data[tmp_nb_var]))[i] += (*(input_data[tmp_nb_var]))[elt_num];
            }
        }

        //divide by number of connected elements (nec) to get the average
        if(nec>0)
            for(int tmp_nb_var=0;tmp_nb_var<nb_var;tmp_nb_var++)
                (*(output_data[tmp_nb_var]))[i] /= nec;
    }

}//elementsToNodes()


template<typename T>
void MeshInfo<T>::interpToGrid(
        T_val_vec_ptrs &output_data,            //output data
        T_val_vec_ptrs const &input_data,       //input data
        T_val const &xmin, T_val const &ymax,   //grid parameters
        int const&nx, int const&ny,
        T_val const &dx, T_val const &dy,
        T_val const &default_value)

{

    if(M_mesh_type==E_mesh_type::uninitialised)
        throw runtime_error("interpToGrid: Mesh object uninitialised\n");

    int nb_var      = input_data.size();
    int Ninterp     = (*(input_data[0])).size();//get pointer, then get size

    T_val_vec interp_in(Ninterp*nb_var);   //input to interp routine
    for (int i=0;i<Ninterp;i++)
        for (int p=0;p<nb_var;p++)
            interp_in[nb_var*i+p]   = (*(input_data[p]))[i];

    T_val* interp_out;

    std::cout<<"before interp mesh2grid\n";
    InterpFromMeshToGridx(interp_out,     //output data (pointer)
                          &M_index[0],    //mesh index (element->node map)
                          &M_nodes_x[0],  //node positions (x-coords)
                          &M_nodes_y[0],  //node positions (y-coords)
                          M_num_nodes,    //num nodes
                          M_num_elements, //num elements
                          &interp_in[0],  //input data
                          Ninterp,        //data length
                          nb_var,         //no of variables
                          xmin,           //xmin (of grid elements' positions)
                          ymax,           //ymax (of grid elements' positions)
                          dx,             //grid dx
                          dy,             //grid dy
                          nx,             //grid nx
                          ny,             //grid ny
                          default_value); //default value (given to points outside mesh)


    //output
    for (int p=0;p<nb_var;p++)
    {
        output_data[p]->assign(nx*ny,0.);
        for (int i=0;i<nx*ny;i++)
            (*(output_data[p]))[i]  = interp_out[nb_var*i+p];
    }

    xDelete<T_val>(interp_out);

}//interpToGrid
// ===============================================================


// ===============================================================
// general functions
namespace MeshTools
{


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
jacobian(T_gmsh_el const& element, T_gmsh const& mesh)
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
jacobian(T_gmsh_el const& element, T_gmsh const& mesh,
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
sides(T_gmsh_el const& element, T_gmsh const& mesh)
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
minMaxSide(T_gmsh const& mesh)
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
minAngles(T_gmsh_el const& element, T_gmsh const& mesh)
{
    std::vector<double> side = MeshTools::sides(element,mesh);
    //std::for_each(side.begin(), side.end(), [&](double& f){ f = 1000.*f; });
    std::sort(side.begin(),side.end());
    double minang = std::acos( (std::pow(side[1],2.) + std::pow(side[2],2.) - std::pow(side[0],2.) )/(2*side[1]*side[2]) );
    minang = minang*45.0/std::atan(1.0);

    return minang;
}


double
minAngle(T_gmsh const& mesh)
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
minAngle(T_gmsh const& mesh, std::vector<double> const& um, double factor)
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
flip(T_gmsh const& mesh, std::vector<double> const& um, double factor)
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
resolution(T_gmsh const& mesh)
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
hminVertices(T_gmsh const& mesh, BamgMesh const* bamg_mesh)
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
hmaxVertices(T_gmsh const& mesh, BamgMesh const* bamg_mesh)
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
AllMinAngle(T_gmsh const& mesh, std::vector<double> const& um, double factor)
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
measure(T_gmsh_el const& element, T_gmsh const& mesh)
{
    return (1./2)*std::abs(jacobian(element,mesh));
}//measure


double
measure(T_gmsh_el const& element, T_gmsh const& mesh,
                       std::vector<double> const& um, double factor)
{
    return (1./2)*std::abs(jacobian(element,mesh,um,factor));
}//measure


std::vector<double>
shapeCoeff(T_gmsh_el const& element, T_gmsh const& mesh)
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

// instantiate wim class for type double
//template class MeshInfo<float>;

// instantiate wim class for type double
template class MeshInfo<double>;

} // namespace Wim
