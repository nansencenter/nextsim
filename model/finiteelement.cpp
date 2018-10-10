/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   finiteelement.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Mon Aug 24 11:02:45 2015
 */

#include <finiteelement.hpp>
#include <constants.hpp>
#include <date.hpp>
#include <redistribute.hpp>
#include <exporter.hpp>
#include <numeric>

#define GMSH_EXECUTABLE gmsh

namespace Nextsim
{

//------------------------------------------------------------------------------------------------------
//!Despite its name, this is the main model file. All functions pertaining to neXtSIM are defined here.
FiniteElement::FiniteElement()
    :
    vm(Environment::vm()),
    M_mesh(),
    M_solver(),
    M_matrix(),
    M_vector(),
    timer()
{}

    
//------------------------------------------------------------------------------------------------------
//! Initialisation of the mesh.
//! Called by the init() function.
void
FiniteElement::initMesh()
{
    this->initBamg();

    M_comm = M_mesh.comm();
    M_rank = M_comm.rank();

    this->rootMeshProcessing();

    if (!M_use_restart)
    {
        this->distributedMeshProcessing(true);
    }
}//initMesh

    
//------------------------------------------------------------------------------------------------------
//! Distribution of mesh processing for parallel computing.
//! Called by the interpFields(), initMesh() and distributedMeshProcessing() functions.
void
FiniteElement::distributedMeshProcessing(bool start)
{
    M_comm.barrier();

    if (!start)
    {
        M_mesh = mesh_type();
    }

    M_mesh.setOrdering("gmsh");

    if (M_rank == 0)
        LOG(INFO) <<"filename= "<< M_partitioned_mesh_filename <<"\n";

    timer["meshread"].first.restart();
    M_mesh.readFromFile(M_partitioned_mesh_filename, M_mesh_fileformat);
    if (M_rank == 0)
        std::cout<<"-------------------MESHREAD done in "<< timer["meshread"].first.elapsed() <<"s\n";

    if (!start)
    {
        delete bamggeom;
        delete bamgmesh;

        bamggeom = new BamgGeom();
        bamgmesh = new BamgMesh();
    }

    timer["createbamg"].first.restart();
    BamgConvertMeshx(
                     bamgmesh,bamggeom,
                     &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                     M_mesh.numNodes(), M_mesh.numTriangles()
                     );

    if (M_rank == 0)
        std::cout<<"-------------------CREATEBAMG done in "<< timer["createbamg"].first.elapsed() <<"s\n";

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_ndof = M_mesh.numGlobalNodes();

    M_local_ndof = M_mesh.numLocalNodesWithoutGhost();
    M_local_ndof_ghost = M_mesh.numLocalNodesWithGhost();

    M_local_nelements = M_mesh.numTrianglesWithoutGhost();
    M_num_nodes = M_local_ndof_ghost;

    timer["bcmarker"].first.restart();
    this->bcMarkedNodes();
    if (M_rank == 0)
        std::cout<<"-------------------BCMARKER done in "<< timer["bcmarker"].first.elapsed() <<"s\n";

    timer["creategraph"].first.restart();
    this->createGraph();
    if (M_rank == 0)
        std::cout<<"-------------------CREATEGRAPH done in "<< timer["creategraph"].first.elapsed() <<"s\n";

    timer["gathersize"].first.restart();
    this->gatherSizes();
    if (M_rank == 0)
        std::cout<<"-------------------GATHERSIZE done in "<< timer["gathersize"].first.elapsed() <<"s\n";

#if 1
    timer["scattercvt"].first.restart();
    this->scatterElementConnectivity();
    //if (M_rank == 0)
    std::cout<<"-------------------CONNECTIVITY done in "<< timer["gathersize"].first.elapsed() <<"s\n";
#endif

#if 0
    // LOG(INFO) <<"["<< M_rank << "] NODES   = "<< M_mesh.numGlobalNodes() << " --- "<< M_local_ndof <<"\n";
    // LOG(INFO) <<"["<< M_rank << "] ELEMENTS= "<< M_mesh.numGlobalElements() << " --- "<< M_local_nelements <<"\n";

    std::cout << "[INFO]: " <<"["<< M_rank << "] NODES   = "<< M_mesh.numGlobalNodes() << " --- "<< M_local_ndof <<"\n";
    std::cout << "[INFO]: " <<"["<< M_rank << "] ELEMENTS= "<< M_mesh.numGlobalElements() << " --- "<< M_local_nelements <<"\n";

    int num_nodes = boost::mpi::all_reduce(M_comm, M_local_ndof, std::plus<int>());
    int num_elements = boost::mpi::all_reduce(M_comm, M_local_nelements, std::plus<int>());

    std::cout<<"NODE COMPARE: "<< M_mesh.numGlobalNodesFromSarialMesh() << " and "<< num_nodes <<"\n";

    if(M_mesh.numGlobalNodesFromSarialMesh() != num_nodes)
    {
        throw std::logic_error("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@INCONSISTANT NODAL PARTITIONS");
    }

    std::cout<<"ELEMENT COMPARE: "<< M_mesh.numGlobalElementsFromSarialMesh() << " and "<< num_elements <<"\n";

    if(M_mesh.numGlobalElementsFromSarialMesh() != num_elements)
    {
        throw std::logic_error("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@INCONSISTANT ELEMENT PARTITIONS");
    }
    //M_comm.barrier();
    //std::abort();
#endif

}//distributedMeshProcessing
    
    
//------------------------------------------------------------------------------------------------------
//! Marks the mesh nodes where boundary conditions should apply.
//! Called by the distributedMeshProcessing() function.
void
FiniteElement::bcMarkedNodes()
{
#if 0
    M_dirichlet_flags_root.resize(5);
    M_neumann_flags_root.resize(11);

    M_dirichlet_flags_root = {2,3,11,12,13};
    M_neumann_flags_root = {4,5,6,7,1,8,9,10,14,15,16};

    // M_dirichlet_flags_root.resize(16);
    // std::iota(M_dirichlet_flags_root.begin(), M_dirichlet_flags_root.end(), 1);
    // M_neumann_flags_root.resize(0);

#endif

    //std::cout<<"["<< M_rank << "] NDOFS= "<< M_num_nodes << " --- "<< M_local_ndof <<"\n";

    std::vector<int> flags_size_root(2);
    if (M_rank == 0)
    {
        flags_size_root = {(int)M_dirichlet_flags_root.size(), (int)M_neumann_flags_root.size()};
    }

    boost::mpi::broadcast(M_comm, &flags_size_root[0], 2, 0);

    int dir_size = flags_size_root[0];
    int nmn_size = flags_size_root[1];

    std::vector<int> flags_root;
    if (M_rank == 0)
    {
        std::copy(M_dirichlet_flags_root.begin(), M_dirichlet_flags_root.end(), std::back_inserter(flags_root));
        std::copy(M_neumann_flags_root.begin(), M_neumann_flags_root.end(), std::back_inserter(flags_root));

#if 0
        for (int i=0; i<M_dirichlet_flags_root.size(); ++i)
        {
            std::cout<<"M_dirichlet_flags_root["<< i <<"]= "<< M_dirichlet_flags_root[i] <<"\n";
        }

        for (int i=0; i<M_neumann_flags_root.size(); ++i)
        {
            std::cout<<"M_neumann_flags_root["<< i <<"]= "<< M_neumann_flags_root[i] <<"\n";
        }

        for (int i=0; i<flags_root.size(); ++i)
        {
            std::cout<<"flags_root["<< i <<"]= "<< flags_root[i] <<"\n";
        }
#endif
    }

    if (M_rank != 0)
    {
        flags_root.resize(dir_size+nmn_size);
    }

    // broadcast the dirichlet and neumann nodes from master process to other processes
    boost::mpi::broadcast(M_comm, &flags_root[0], dir_size+nmn_size, 0);

    auto transfer_bimap = M_mesh.transferMapInit();

    M_dirichlet_flags.resize(0);
    M_neumann_flags.resize(0);

    // We mask out the boundary nodes
    M_mask.assign(M_num_nodes,false);
    M_mask_dirichlet.assign(M_num_nodes,false);

    for (int i=0; i<flags_root.size(); ++i)
    {
        if (transfer_bimap.left.find(flags_root[i]) != transfer_bimap.left.end())
        {
            int lindex = transfer_bimap.left.find(flags_root[i])->second-1;

            if (i < dir_size)
            {
                // exclude ghost nodes
                if (lindex < M_local_ndof)
                {
                    M_dirichlet_flags.push_back(lindex);
                    //std::cout<<"["<< M_comm.rank() <<"] " << "-----------------here  = "<< flags_root[i]  <<"\n";
                    // add mask for dirichlet nodes
                    M_mask_dirichlet[lindex] = true;
                }
            }
            else
            {
                // including ghost nodes
                M_neumann_flags.push_back(lindex);
                //std::cout<<"["<< M_comm.rank() <<"] " << "-----------------here  = "<< flags_root[i]  <<"\n";
            }

            // add mask for boundary nodes
            M_mask[lindex] = true;
        }
    }

    std::sort(M_dirichlet_flags.begin(), M_dirichlet_flags.end());
    M_dirichlet_flags.erase(std::unique(M_dirichlet_flags.begin(), M_dirichlet_flags.end() ), M_dirichlet_flags.end());

    std::sort(M_neumann_flags.begin(), M_neumann_flags.end());
    M_neumann_flags.erase(std::unique(M_neumann_flags.begin(), M_neumann_flags.end() ), M_neumann_flags.end());

    LOG(DEBUG) <<"["<< M_comm.rank() <<"] " << "Dirichlet flags= "<< M_dirichlet_flags.size() <<"\n";
    LOG(DEBUG) <<"["<< M_comm.rank() <<"] " << "Neumann flags  = "<< M_neumann_flags.size() <<"\n";


    M_dirichlet_nodes.resize(2*(M_dirichlet_flags.size()));
    for (int i=0; i<M_dirichlet_flags.size(); ++i)
    {
        M_dirichlet_nodes[2*i] = M_dirichlet_flags[i];
        M_dirichlet_nodes[2*i+1] = M_dirichlet_flags[i]+M_num_nodes;
    }

    M_neumann_nodes.resize(2*(M_neumann_flags.size()));
    for (int i=0; i<M_neumann_flags.size(); ++i)
    {
        M_neumann_nodes[2*i] = M_neumann_flags[i];
        M_neumann_nodes[2*i+1] = M_neumann_flags[i]+M_num_nodes;
    }
}//bcMarkedNodes

    
//------------------------------------------------------------------------------------------------------
//! Reads, converts and applies boundary conditions to the mesh.
//! Called by the initMesh() function.
void
FiniteElement::rootMeshProcessing()
{
    if (M_rank == 0)
    {

        // read the original input mesh
        M_mesh_root.setOrdering("gmsh");
        LOG(DEBUG) <<"Reading root mesh starts\n";
        chrono.restart();
        M_mesh_root.readFromFile(M_mesh_filename);
        LOG(DEBUG) <<"Reading root mesh done in "<< chrono.elapsed() <<"s\n";

        chrono.restart();
        M_mesh_root.stereographicProjection();
        LOG(DEBUG) <<"Projection root mesh done in "<< chrono.elapsed() <<"s\n";

        M_mesh_init_root = M_mesh_root;

        LOG(DEBUG) <<"Convert mesh starts\n";
        BamgConvertMeshx(
                         bamgmesh_root,bamggeom_root,
                         &M_mesh_root.indexTr()[0],&M_mesh_root.coordX()[0],&M_mesh_root.coordY()[0],
                         M_mesh_root.numNodes(), M_mesh_root.numTriangles()
                         );

        // ------ Boundary conditions -----------


        // for (auto it=M_mesh_root.markerNames().begin(), end=M_mesh_root.markerNames().end(); it!=end; ++it)
        // {
        //     std::cout<<"*******************************markerNames["<< it->first <<"]= ("<< it->second[0] << "," << it->second[1] <<")\n";
        // }

        // set M_flag_fix to its correct value when PhysicalNames section is present in the msh file (version 2.2)
        if (!(M_mesh_root.markerNames().empty()))
        {
            // get the id associated to the physical name "coast" and assign it to M_flag_fix
            M_flag_fix = M_mesh_root.markerNames().find("coast")->second[0];
        }
        else
        {
            throw std::runtime_error(
                    "No \"coast\" marker in mesh file. Check your input file (" + M_mesh_filename + ")");
        }

        for (auto it=M_mesh_root.edges().begin(), end=M_mesh_root.edges().end(); it!=end; ++it)
        {
            if (it->physical==M_flag_fix)
            {
                M_dirichlet_flags_root.push_back(it->indices[0]/*-1*/);
                M_dirichlet_flags_root.push_back(it->indices[1]/*-1*/);
            }
        }

        std::sort(M_dirichlet_flags_root.begin(), M_dirichlet_flags_root.end());
        M_dirichlet_flags_root.erase(std::unique(M_dirichlet_flags_root.begin(), M_dirichlet_flags_root.end() ), M_dirichlet_flags_root.end());

        // for (int i=0; i<M_dirichlet_flags_root.size(); ++i)
        // {
        //     LOG(DEBUG) <<"M_dirichlet_flags_root["<< i <<"]= "<< M_dirichlet_flags_root[i] <<"\n";
        // }

        LOG(DEBUG) <<"Convert mesh done\n";

        // Definition of the hmin, hmax, hminVertices or hmaxVertices
        auto h = this->minMaxSide(M_mesh_root);
        M_res_root_mesh = this->resolution(M_mesh_root);

        std::cout <<"MESH: HMIN= "<< h[0] <<"\n";
        std::cout <<"MESH: HMAX= "<< h[1] <<"\n";
        std::cout <<"MESH: RESOLUTION= "<< M_res_root_mesh <<"\n";

        switch (M_mesh_type)
        {
        case setup::MeshType::FROM_UNREF:
            // For the other meshes, we use a constant hmin and hmax
            bamgopt->hmin = h[0];
            bamgopt->hmax = h[1];
            break;
        case setup::MeshType::FROM_SPLIT:
            bamgopt->hmin = h[0];
            bamgopt->hmax = h[1];

            M_hminVertices = this->hminVertices(M_mesh_init_root, bamgmesh_root);
            M_hmaxVertices = this->hmaxVertices(M_mesh_init_root, bamgmesh_root);

            // LOG(DEBUG) <<"HMIN MIN= "<< *std::min_element(M_hminVertices.begin(), M_hminVertices.end()) <<"\n";
            // LOG(DEBUG) <<"HMIN MAX= "<< *std::max_element(M_hminVertices.begin(), M_hminVertices.end()) <<"\n";
            // LOG(DEBUG) <<"HMAX MIN= "<< *std::min_element(M_hmaxVertices.begin(), M_hmaxVertices.end()) <<"\n";
            // LOG(DEBUG) <<"HMAX MAX= "<< *std::max_element(M_hmaxVertices.begin(), M_hmaxVertices.end()) <<"\n";

            bamgopt->hminVertices = new double[M_mesh_init_root.numNodes()];
            bamgopt->hmaxVertices = new double[M_mesh_init_root.numNodes()];
            for (int i=0; i<M_mesh_init_root.numNodes(); ++i)
            {
                bamgopt->hminVertices[i] = M_hminVertices[i];
                bamgopt->hmaxVertices[i] = M_hmaxVertices[i];
            }
            break;
        default:
            LOG(DEBUG)  << "invalid mesh type"<<"\n";
            throw std::logic_error("invalid mesh type");
        }

        if(M_mesh_type==setup::MeshType::FROM_SPLIT)
        {

            chrono.restart();
            LOG(DEBUG) <<"First adaptation starts\n";
            // step 1 (only for the first time step): Start by having bamg 'clean' the mesh with KeepVertices=0
            bamgopt->KeepVertices=0;
            bamgopt->splitcorners=1;
            this->adaptMesh();
            bamgopt->KeepVertices=1;
            bamgopt->splitcorners=0;
            LOG(DEBUG) <<"First adaptation done in "<< chrono.elapsed() <<"s\n";

            // Interpolate hminVertices and hmaxVertices onto the current mesh
            this->interpVertices();
        }

        if (!M_use_restart)
        {
            chrono.restart();
            LOG(DEBUG) <<"AdaptMesh starts\n";
            this->adaptMesh();
            LOG(DEBUG) <<"AdaptMesh done in "<< chrono.elapsed() <<"s\n";

            // Add information on the number of partition to mesh filename
            LOG(DEBUG) <<"["<< M_rank <<"] " <<"filename= "<< M_partitioned_mesh_filename <<"\n";

            std::cout<<"------------------------------version       = "<< M_mesh_root.version() <<"\n";
            std::cout<<"------------------------------ordering      = "<< M_mesh_root.ordering() <<"\n";
            std::cout<<"------------------------------format        = "<< M_mesh_fileformat <<"\n";
            std::cout<<"------------------------------space         = "<< vm["mesh.partitioner-space"].as<std::string>() <<"\n";
            std::cout<<"------------------------------partitioner   = "<< vm["mesh.partitioner"].as<std::string>() <<"\n";


            // save mesh (only root process)
            chrono.restart();
            if (M_partition_space == mesh::PartitionSpace::MEMORY)
            {
                // Environment::logMemoryUsage("before gmodel...");
                M_mesh_root.initGModel();
                M_mesh_root.writeToGModel();
                // Environment::logMemoryUsage("before after...");
            }
            else if (M_partition_space == mesh::PartitionSpace::DISK)
                M_mesh_root.writeToFile(M_partitioned_mesh_filename);
            //LOG(DEBUG) <<"Saving mesh done in "<< chrono.elapsed() <<"s\n";
            std::cout <<"Writing mesh done in "<< chrono.elapsed() <<"s\n";

            // partition the mesh on root process (rank 0)
            chrono.restart();
            M_mesh_root.partition(M_partitioned_mesh_filename,
                    M_partitioner, M_partition_space, M_mesh_fileformat);
            //LOG(DEBUG) <<"Partitioning mesh done in "<< chrono.elapsed() <<"s\n";
            std::cout <<"Partitioning mesh done in "<< chrono.elapsed() <<"s\n";
        }
    }

    boost::mpi::broadcast(M_comm, M_flag_fix, 0);

}//rootMeshProcessing

    
//------------------------------------------------------------------------------------------------------
//! Performs a re-numbering of the mesh nodes and elements
//! !Does not seem to be used!
void
FiniteElement::rootMeshRenumbering()
{
#if 0
    if (M_rank == 0)
    {
        M_mesh_root.reorder(M_mesh.mapNodes(),M_mesh.mapElements());

        delete bamggeom_root;
        delete bamgmesh_root;

        bamggeom_root = new BamgGeom();
        bamgmesh_root = new BamgMesh();

        std::cout<<"Re-numbering: Convert mesh starts\n";

        chrono.restart();

        BamgConvertMeshx(
                         bamgmesh_root,bamggeom_root,
                         &M_mesh_root.indexTr()[0],&M_mesh_root.coordX()[0],&M_mesh_root.coordY()[0],
                         M_mesh_root.numNodes(), M_mesh_root.numTriangles()
                         );
        std::cout<<"Re-numbering: Convert mesh done in "<< chrono.elapsed() <<"s\n";
    }
#endif
}//rootMeshRenumbering

    
//------------------------------------------------------------------------------------------------------
//! Initializes the size of all physical variables with values set to zero.
//! Called by the init() and readRestart() functions.
//! \note
//! - The prefix D_ is used for diagnostic variables (outputs),
//! - The prefix M_ is used for global variables of the finite element class, accessible for all functions defined in finiteelement.cpp.
//! - The suffix M is used for a quantity at the previous (nth) time step (e.g., VTM)
//! - The suffix M is used for a quantity at the second-previous (n-1 th) time step (e.g., VTMM)

void
FiniteElement::initVariables()
{
    chrono_tot.restart();

    // Global variables are assigned the prefix M_
    M_nb_regrid = 0; //! \param M_nb_regrid (int) Number of times remeshing has been called since the beginning of the run

    M_solver = solver_ptrtype(new solver_type());
    M_matrix = matrix_ptrtype(new matrix_type());
    M_vector = vector_ptrtype(new vector_type());
    M_solution = vector_ptrtype(new vector_type());

    M_reuse_prec = true;

    M_VT.resize(2*M_num_nodes,0.); //! \param M_VT (double) Instantaneous velocity vector at the (n+1)th (current) t-step [m/s]
    M_VTM.resize(2*M_num_nodes,0.); //! \param M_VTM (double) Instantaneous velocity vector at the nth t-step [m/s]
    M_VTMM.resize(2*M_num_nodes,0.); //! \param M_VTMM (double) Instantaneous velocity vector at the (n-1)th t-step [m/s]
    
    M_h_thin.assign(M_num_elements,0.); //! \param M_h_thin (double) Thickness of thin ice [m]
    M_conc_thin.assign(M_num_elements,0.); //! \param M_conc_thin (double) Concentration thin ice
    M_hs_thin.assign(M_num_elements,0.); //! \param M_hs_thin (double) Thickness of snow on top of thin ice [m]
    M_tsurf_thin.assign(M_num_elements,0.); //! \param M_tsurf_thin (double) Temperature at the surface of thin ice [C]
    
    // stresses
    M_sigma.assign(3*M_num_elements,0.); //! \param M_sigma (double) Internal stress tensor [N/m2]
    
    // random numbers
    //M_random_number.resize(M_num_elements);

    M_random_number_root.resize(M_mesh.numGlobalElements());

    if (M_rank == 0)
    {
        boost::minstd_rand intgen;
        boost::uniform_01<boost::minstd_rand> gen(intgen);

        for (int i=0; i<M_random_number_root.size(); ++i)
        {
            M_random_number_root[i] = gen();
            //M_random_number_root[i] = static_cast <double> (std::rand()) / static_cast <double> (RAND_MAX);
        }
    }

#if 0
    //std::vector<int> sizes_elements = M_sizes_elements_with_ghost;
    M_comm.barrier();

    if (M_rank == 0)
    {
        //std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });
        boost::mpi::scatterv(M_comm, M_random_number_root, M_sizes_elements_with_ghost, &M_random_number[0], 0);
    }
    else
    {
        boost::mpi::scatterv(M_comm, &M_random_number[0], M_num_elements, 0);
    }

    std::cout<<"TERMINATED........................................\n";
#endif

#if 1
    // if (M_rank == 0)
    // {
    //     std::cout<<"Random MIN= "<< *std::min_element(M_random_number_root.begin(), M_random_number_root.end()) <<"\n";
    //     std::cout<<"Random MAX= "<< *std::max_element(M_random_number_root.begin(), M_random_number_root.end()) <<"\n";
    // }

    boost::mpi::broadcast(M_comm, &M_random_number_root[0], M_mesh.numGlobalElements(), 0);

    M_random_number.resize(M_num_elements);
    auto id_elements = M_mesh.trianglesIdWithGhost();

    for (int i=0; i<M_random_number.size(); ++i)
        M_random_number[i] = M_random_number_root[id_elements[i]-1];
#endif

    M_conc.resize(M_num_elements); //! \param M_conc (double) Concentration of thick ice
    M_thick.resize(M_num_elements); //! \param M_thick (double) Thickness of thick ice [m]
    M_damage.resize(M_num_elements); //! \param M_damage (double) Level of damage
    M_ridge_ratio.assign(M_num_elements,0.); //! \param M_ridge_ratio (double) Ratio of ridged vs unridged ice
    M_snow_thick.resize(M_num_elements); //! \param M_snow_thick (double) Snow thickness (on top of thick ice) [m]
    
    M_sst.resize(M_num_elements); //! \param M_sst (double) Sea surface temperature [C]
    M_sss.resize(M_num_elements); //! \param M_sss (double) Sea surface salinity [C]
    
    switch (M_thermo_type)
    {
        case (setup::ThermoType::ZERO_LAYER):
            M_tice.resize(1);   //! \param M_tice (double) Ice temperature [C]
            break;
        case (setup::ThermoType::WINTON):
            M_tice.resize(3);
            break;
        default:
            std::cout << "thermo_type= " << (int)M_thermo_type << "\n";
            throw std::logic_error("Wrong thermo_type");
    }

    for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
        it->assign(M_num_elements,0.);

    for (int i=0; i<M_num_elements; ++i)
    {
        M_sigma[3*i]=0.;
        M_sigma[3*i+1]=0.;
        M_sigma[3*i+2]=0.;

        if ((M_conc[i] <= 0.) || (M_thick[i] <= 0.) )
        {
            M_conc[i] = 0.;
            M_thick[i] = 0.;
        }
    }

    // Diagnostic variables are assigned the prefix D_
    D_Qa.resize(M_num_elements); //! \param D_Qa (double) Total heat flux to the atmosphere
    D_Qsh.resize(M_num_elements); //! \param D_Qsh (double) Sensible heat flux to the atmosphere
    D_Qlh.resize(M_num_elements); //! \param D_Qlh (double) Latent heat flux to the atmosphere
    D_Qlw.resize(M_num_elements); //! \param D_Qlw (double) Long wave heat flux to the atmosphere
    D_Qsw.resize(M_num_elements); //! \param D_Qsw (double) Short wave heat flux to the atmosphere
    D_Qo.resize(M_num_elements); //! \param D_Qo (double) Total heat lost by the ocean
    D_delS.resize(M_num_elements); //! \param D_delS (double) Salt release to the ocean [kg/day]
    
    M_UT.assign(2*M_num_nodes,0.); //! \param M_UT (double) Total ice displacement (M_UT[] = time_step*M_VT[]) [m]
    
    
    if (M_rank == 0)
    {
        M_surface_root.assign(M_mesh_root.numTriangles(),0.);

        int cpt = 0;
        for (auto it=M_mesh_root.triangles().begin(), end=M_mesh_root.triangles().end(); it!=end; ++it)
        {
            M_surface_root[cpt] = this->measure(*it,M_mesh_root);
            ++cpt;
        }
    }

    if ((M_rank == 0) && (M_use_drifters))
    {
        M_UT_root.assign(2*M_ndof,0.);
        M_conc_root.resize(M_mesh_root.numTriangles());
    }

    this->assignVariables();
}//initVariables

    
//------------------------------------------------------------------------------------------------------
//! Calls the functions for assimilation of ocean and ice data: assimilateSlabOcean() and assimilateIce().
//! Called by the init() function.
void
FiniteElement::DataAssimilation()
{

    LOG(DEBUG) << "assimilateSlabOcean\n";
    this->assimilateSlabOcean();

    LOG(DEBUG) << "assimilateIce\n";
    this->assimilateIce();
}//DataAssimilation

    
//------------------------------------------------------------------------------------------------------
//! Assigns variables in the context of remeshing : the size of variables needs to be update when remeshing because the nb of elements/nodes has changed.
//! Called by the regrid() and initVariables() functions.
void
FiniteElement::assignVariables()
{
    M_surface.assign(M_num_elements,0.);

    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        M_surface[cpt] = this->measure(*it,M_mesh);
        ++cpt;
    }

    M_matrix->init(2*M_ndof,2*M_ndof,
                   2*M_local_ndof,2*M_local_ndof,
                   M_graphmpi);

    M_vector->init(2*M_ndof,2*M_local_ndof,M_graphmpi);
    M_solution->init(2*M_ndof,2*M_local_ndof,M_graphmpi);

    M_UM.assign(2*M_num_nodes,0.);
    //M_UT.assign(2*M_num_nodes,0.);

    M_fcor.assign(M_num_elements, 0.);

#if 1
    // reload the dataset
    M_atmosphere_nodes_dataset.loaded=false;
    M_atmosphere_elements_dataset.loaded=false;
    M_atmosphere_bis_elements_dataset.loaded=false;
    M_ocean_nodes_dataset.loaded=false;
    M_ocean_elements_dataset.loaded=false;

    M_ice_topaz_elements_dataset.loaded=false;
    M_ice_piomas_elements_dataset.loaded=false;
    M_ice_amsre_elements_dataset.loaded=false;
    M_ice_osisaf_elements_dataset.loaded=false;
    M_ice_osisaf_type_elements_dataset.loaded=false;
    M_ice_amsr2_elements_dataset.loaded=false;
    M_ice_cs2_smos_elements_dataset.loaded=false;
    M_ice_smos_elements_dataset.loaded=false;
    M_bathymetry_elements_dataset.loaded=false;


    // reload the grid
    M_atmosphere_nodes_dataset.grid.loaded=false;
    M_atmosphere_elements_dataset.grid.loaded=false;
    M_atmosphere_bis_elements_dataset.grid.loaded=false;
    M_ocean_nodes_dataset.grid.loaded=false;
    M_ocean_elements_dataset.grid.loaded=false;

    M_ice_topaz_elements_dataset.grid.loaded=false;
    M_ice_piomas_elements_dataset.grid.loaded=false;
    M_ice_amsre_elements_dataset.grid.loaded=false;
    M_ice_osisaf_elements_dataset.grid.loaded=false;
    M_ice_osisaf_type_elements_dataset.grid.loaded=false;
    M_ice_amsr2_elements_dataset.grid.loaded=false;
    M_ice_cs2_smos_elements_dataset.grid.loaded=false;
    M_ice_smos_elements_dataset.grid.loaded=false;
    M_bathymetry_elements_dataset.grid.loaded=false;

    // --------------------------------------------------------------
    // interpolation of the dataset
    M_atmosphere_nodes_dataset.interpolated=false;
    M_atmosphere_elements_dataset.interpolated=false;
    M_atmosphere_bis_elements_dataset.interpolated=false;
    M_ocean_nodes_dataset.interpolated=false;
    M_ocean_elements_dataset.interpolated=false;

    M_ice_topaz_elements_dataset.interpolated=false;
    M_ice_piomas_elements_dataset.interpolated=false;
    M_ice_amsre_elements_dataset.interpolated=false;
    M_ice_osisaf_elements_dataset.interpolated=false;
    M_ice_osisaf_type_elements_dataset.interpolated=false;
    M_ice_amsr2_elements_dataset.interpolated=false;
    M_ice_cs2_smos_elements_dataset.interpolated=false;
    M_ice_smos_elements_dataset.interpolated=false;
    M_bathymetry_elements_dataset.interpolated=false;
    // --------------------------------------------------------------
#endif

    // //loop over vector of pointers to datasets defined in initForcings()
    // for (auto it=M_datasets_regrid.begin(), end=M_datasets_regrid.end(); it!=end; ++it)
    // {
    //     if (M_rank == 0)
    //         std::cout<<"REGRIDDING: need to re-interpolate dataset "<<(*it)->name<<"\n";

    //     (*it)->interpolated=false;

    //     // for the parallel code, it will be necessary to add these lines
    //     // as the domain covered by the partitions changes at each remeshing/partitioning
    //     // (*it)->grid.interpolated=false;
    //     (*it)->grid.loaded=false;
    // }

    M_Cohesion.resize(M_num_elements); // \param M_Cohesion (double) Ice cohesive strength [N/m2]
    M_Compressive_strength.resize(M_num_elements); // \param M_Compressive_strength (double) Ice maximum compressive strength [N/m2]
    M_time_relaxation_damage.resize(M_num_elements,time_relaxation_damage); // \param M_time_relaxation_damage (double) Characteristic time for healing [s]
    
#if 1
    // root
    // M_UM_root.assign(2*M_mesh.numGlobalNodes(),0.);

    if (M_rank == 0)
    {
        M_surface_root.assign(M_mesh_root.numTriangles(),0.);

        cpt = 0;
        for (auto it=M_mesh_root.triangles().begin(), end=M_mesh_root.triangles().end(); it!=end; ++it)
        {
            M_surface_root[cpt] = this->measure(*it,M_mesh_root);
            ++cpt;
        }
    }
#endif

    // number of variables to interpolate
    M_nb_var_element = /*11*/15 + M_tice.size();
}//assignVariables
    

//------------------------------------------------------------------------------------------------------
//! Initializes the physical state of the model - ocean and ice - by calling the initSlabOcean() and initIce() functions.
//! Called by the init() function.
void
FiniteElement::initModelState()
{

    LOG(DEBUG) << "initSlabOcean\n";
    this->initSlabOcean();

    LOG(DEBUG) << "initIce\n";
    this->initIce();
}//initModelState

    
//------------------------------------------------------------------------------------------------------
//! Initializes the datasets for forcing. Replaced the initDatasets() function of the serial code.
//! Called by the init() function.
void
FiniteElement::initForcings()
{
    //! - 1) Initializes the dataset definitions,
    switch(M_atmosphere_type){
        case setup::AtmosphereType::CONSTANT:
            break;

        case setup::AtmosphereType::ASR:
            M_atmosphere_nodes_dataset=DataSet("asr_nodes");
            M_atmosphere_elements_dataset=DataSet("asr_elements");
            break;

        case setup::AtmosphereType::ERAi:
            M_atmosphere_nodes_dataset=DataSet("ERAi_nodes");
            M_atmosphere_elements_dataset=DataSet("ERAi_elements");
            break;

        case setup::AtmosphereType::EC:
            M_atmosphere_nodes_dataset=DataSet("ec_nodes");
            M_atmosphere_elements_dataset=DataSet("ec_elements");
            break;

        case setup::AtmosphereType::EC2:
            M_atmosphere_nodes_dataset=DataSet("ec2_nodes");
            M_atmosphere_elements_dataset=DataSet("ec2_elements");
            break;

        case setup::AtmosphereType::EC_ERAi:
            M_atmosphere_nodes_dataset=DataSet("ec_nodes");
            M_atmosphere_elements_dataset=DataSet("ec_elements");
            M_atmosphere_bis_elements_dataset=DataSet("ERAi_elements");
            break;

        case setup::AtmosphereType::CFSR:
            M_atmosphere_nodes_dataset=DataSet("cfsr_nodes");
            M_atmosphere_elements_dataset=DataSet("cfsr_elements");
            break;

        case setup::AtmosphereType::CFSR_HI:
            M_atmosphere_nodes_dataset=DataSet("cfsr_nodes_hi");
            M_atmosphere_elements_dataset=DataSet("cfsr_elements");
            break;

        default:
            std::cout << "invalid wind forcing"<<"\n";throw std::logic_error("invalid wind forcing");
    }

    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
        break;

        case setup::OceanType::TOPAZR:
        case setup::OceanType::TOPAZR_atrest:
            M_ocean_nodes_dataset=DataSet("topaz_nodes");
            M_ocean_elements_dataset=DataSet("topaz_elements");
            break;

        case setup::OceanType::TOPAZR_ALTIMETER:
            M_ocean_nodes_dataset=DataSet("ocean_currents_nodes");
            M_ocean_elements_dataset=DataSet("topaz_elements");
            break;

        case setup::OceanType::TOPAZF:
            M_ocean_nodes_dataset=DataSet("topaz_forecast_nodes");
            M_ocean_elements_dataset=DataSet("topaz_forecast_elements");
            break;

        default:
            std::cout << "invalid ocean forcing"<<"\n";throw std::logic_error("invalid ocean forcing");
    }
    if (M_use_nesting)
    {
        M_nesting_nodes_dataset=DataSet("nesting_nodes");
        M_nesting_ocean_elements_dataset=DataSet("nesting_ocean_elements");
        M_nesting_ice_elements_dataset=DataSet("nesting_ice_elements");
        M_nesting_dynamics_elements_dataset=DataSet("nesting_dynamics_elements");
        M_nesting_distance_elements_dataset=DataSet("nesting_distance_elements");
        M_nesting_distance_nodes_dataset=DataSet("nesting_distance_nodes");
    }

    M_ice_topaz_elements_dataset=DataSet("ice_topaz_elements");

    M_ice_icesat_elements_dataset=DataSet("ice_icesat_elements");

    M_ice_piomas_elements_dataset=DataSet("ice_piomas_elements");

    M_ice_amsre_elements_dataset=DataSet("ice_amsre_elements");

    M_ice_osisaf_elements_dataset=DataSet("ice_osisaf_elements");

    M_ice_osisaf_type_elements_dataset=DataSet("ice_osisaf_type_elements");

    M_ice_amsr2_elements_dataset=DataSet("ice_amsr2_elements");

    M_ice_nic_elements_dataset=DataSet("ice_nic_elements");

    M_ice_nic_weekly_elements_dataset=DataSet("ice_nic_weekly_elements");

    M_ice_cs2_smos_elements_dataset=DataSet("ice_cs2_smos_elements");

    M_ice_smos_elements_dataset=DataSet("ice_smos_elements");

    // datasets that need to be re-interpolated after regridding
    // - not needed if only used at initialisation, or if not interpolated onto
    // mesh (eg wave datasets are interpolated onto a fixed grid)
    M_datasets_regrid.push_back(&M_atmosphere_nodes_dataset);
    M_datasets_regrid.push_back(&M_atmosphere_elements_dataset);
    M_datasets_regrid.push_back(&M_atmosphere_bis_elements_dataset);
    M_datasets_regrid.push_back(&M_ocean_nodes_dataset);
    M_datasets_regrid.push_back(&M_ocean_elements_dataset);

    //! - 2) populates the forcing variables.
    LOG(DEBUG) <<"Initialize forcingAtmosphere\n";
    this->forcingAtmosphere();

    LOG(DEBUG) <<"Initialize forcingOcean\n";
    this->forcingOcean();
}//initForcings


//------------------------------------------------------------------------------------------------------
//! Loads and checks on the loading of various datasets.
//! * for external data to be interpolated onto the mesh elements use
//!   RX = M_mesh.bCoordX(), RY = M_mesh.bCoordY()
//! * for external data to be interpolated onto the mesh nodes use
//!   RX = M_mesh.coordX(), RY = M_mesh.coordY()
//! * for external data to be interpolated onto the WIM elements use
//!   RX = M_wim.getX(), RY = M_wim.getY()
//! * NB we don't rotate RX, RY yet since rotation angle is not always defined at this point
//! Called by checkReloadMainDatasets(), and all the ice initialisation and assimilation routines.
void
FiniteElement::checkReloadDatasets(external_data_vec const& ext_data_vec,
        double const& CRtime, std::vector<double> &RX, std::vector<double> &RY)
{
    if ( ext_data_vec.size()==0 )
    {
        LOG(DEBUG) <<"checkReloadDatasets - nothing to do\n";
        return;
    }

    //loop over ext_data_vec and call check and reload for each:
    int i = 0;
    for ( auto it = ext_data_vec.begin(); it != ext_data_vec.end(); ++it, ++i )
    {
        ASSERT((*it)->isInitialized(),
                "checkReloadDatasets: ExternalData object "
                +std::to_string(i) + " is not initialised yet");
        (*it)->check_and_reload(RX, RY, CRtime);
    }
}//checkReloadDatasets


//------------------------------------------------------------------------------------------------------
//! Loads and checks on the loading of the time-dependant datasets.
//! * In practice this is done by looping over the ExternalData objects in
//!   M_external_data_elements and M_external_data_nodes and checking if the corresponding Dataset
//!   needs to be reloaded and/or reinterpolated
//! Called by init() and step()
void
FiniteElement::checkReloadMainDatasets(double const& CRtime)
{
    // check the time-dependant ExternalData objects to see if they need to be reloaded
    // - mesh elements
    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG) <<"checkReloadDatasets (time-dependant elements)\n";
    this->checkReloadDatasets(M_external_data_elements, CRtime, RX, RY);

    // - mesh nodes
    RX = M_mesh.coordX();
    RY = M_mesh.coordY();
    if(M_rank==0)
        LOG(DEBUG) <<"checkReloadDatasets (time-dependant nodes)\n";
    this->checkReloadDatasets(M_external_data_nodes, CRtime, RX, RY);
}//checkReloadMainDatasets

    
//------------------------------------------------------------------------------------------------------
//! Initializes a Bamg mesh grid.
//! Called by the initMesh() function.
void
FiniteElement::initBamg()
{
    bamgopt = new BamgOpts();

    bamgopt->Crack             = 0;
    bamgopt->anisomax          = 1e30;
    bamgopt->coeff             = 1;
    bamgopt->cutoff            = 1e-5;
    //bamgopt->err               = 0.01;
    bamgopt->errg              = 0.1;
    bamgopt->field             = NULL;
    bamgopt->gradation         = 1.5;
    bamgopt->Hessiantype       = 0;
    bamgopt->hmin              = 1e-100;
    bamgopt->hmax              = 1e100;
    bamgopt->hminVertices      = NULL;
    bamgopt->hmaxVertices      = NULL;
    bamgopt->hVertices         = NULL;
    bamgopt->KeepVertices      = 1;
    bamgopt->MaxCornerAngle    = 10;
    bamgopt->maxnbv            = 1e7;
    bamgopt->maxsubdiv         = 10;
    bamgopt->metric            = NULL;
    bamgopt->Metrictype        = 0;
    bamgopt->nbjacobi          = 1;
    bamgopt->nbsmooth          = 3;
    bamgopt->omega             = 1.8;
    bamgopt->power             = 1.;
    bamgopt->splitcorners      = 0; //the Devil!  Changed to 0, original 1 Phil
    bamgopt->geometricalmetric = 0;
    bamgopt->random            = true;
    bamgopt->verbose           = vm["debugging.bamg_verbose"].as<int>();

    bamggeom = new BamgGeom();
    bamgmesh = new BamgMesh();

    if (M_mesh.comm().rank() == 0)
    {
        bamggeom_root = new BamgGeom();
        bamgmesh_root = new BamgMesh();

        bamgopt_previous = new BamgOpts();
        bamggeom_previous = new BamgGeom();
        bamgmesh_previous = new BamgMesh();
    }

    bamgopt->Check();
}//initBamg

    
//------------------------------------------------------------------------------------------------------
//! Defines output options and parameters such as the different time steps (output, thermo, mooring), etc.
//! Called by the init() function.
//! \note These options and parameters are defined in the options.cpp file.
void
FiniteElement::initOptAndParam()
{
    //! Sets the characteristics of the output log (INFOR, WARNING, DEBUG, ERROR),
    const boost::unordered_map<const std::string, LogLevel> str2log = boost::assign::map_list_of
        ("info", INFO)
        ("warning", WARNING)
        ("debug", DEBUG)
        ("error", ERROR);

    M_log_level = str2log.find(vm["debugging.log-level"].as<std::string>())->second;


    //! Defines the export (output) path.
    M_export_path = vm["output.exporter_path"].as<std::string>(); //! \param M_export_path (string) Path of the export files
    // Changes directory for outputs if the option "output.exporter_path" is not empty
    fs::path output_path(M_export_path);

    // Creates the output directory if it does not exist
    if ( (!fs::exists(output_path)) && (M_comm.rank()==0) )
    fs::create_directories(output_path);


    //! Sets Poisson's ratio
    nu0 = vm["dynamics.nu0"].as<double>(); //! \param nu0 (double) Poisson's ratio
    
    //! Sets the Young's modulus, ice density and snow density
    young = vm["dynamics.young"].as<double>(); //! \param young (double) Young modulus of undamaged ice
    rhoi = physical::rhoi; //! \param rhoi (double) Ice density [kg/m3]
    rhos = physical::rhos; //! \param rhos (double) Snow density [kg/m3]
    
    
    //! Sets various time steps (init, thermo, output, mooring, restart) and options on data assimilation and restarts
    days_in_sec = 24.0*3600.0; // Conversion factor from days to seconds
    if (vm["simul.time_init"].as<std::string>() == "")
        throw std::runtime_error("Please provide simul.time_init option (start time)\n");
    else
    time_init = Nextsim::from_date_time_string(vm["simul.time_init"].as<std::string>()); //! \param time_init (string) Time at which the simulation is started
    ptime_step =  days_in_sec/vm["debugging.ptime_per_day"].as<int>(); //! \param ptime_step (int) Debugging time step?
    
    time_step = vm["simul.timestep"].as<double>(); //! \param time_step (double) Model time step [s]
    thermo_timestep = vm["simul.thermo_timestep"].as<double>(); //! \param thermo_timestep (double) Thermodynamic time step [s]
    if ( fmod(thermo_timestep,time_step) != 0)
    {
        std::cout << thermo_timestep << " " << time_step << "\n";
        throw std::runtime_error("thermo_timestep is not an integer multiple of time_step");
    }
    // Temporarly disabling super-stepping of the thermodynamcis. The model hangs randomly when it's enabled
    thermo_timestep = time_step;

    output_time_step =  (vm["output.output_per_day"].as<int>()<0) ? time_step : time_step * floor(days_in_sec/vm["output.output_per_day"].as<int>()/time_step); //! \param output_time_step (int) Time step of model outputs
    mooring_output_time_step =  vm["moorings.output_timestep"].as<double>()*days_in_sec; //! \param mooring_output_time_step (double) Time step for mooring outputs [s]
    mooring_time_factor = time_step/mooring_output_time_step; 
    if ( fmod(mooring_output_time_step,time_step) != 0)
    {
        std::cout << mooring_output_time_step << " " << time_step << "\n";
        throw std::runtime_error("mooring_output_time_step is not an integer multiple of time_step");
    }

    duration = (vm["simul.duration"].as<double>())*days_in_sec; //! \param duration (double) Duration of the simulation [s]
    restart_time_step =  vm["restart.output_time_step"].as<double>()*days_in_sec; //! \param restart_time_step (double) Time step for outputting restart files [s]
    M_use_assimilation   = vm["setup.use_assimilation"].as<bool>(); //! \param M_use_assimilation (boolean) Option on using data assimilation
    M_use_restart   = vm["restart.start_from_restart"].as<bool>(); //! \param M_write_restart (boolean) Option on using starting simulation from a restart file
    M_write_restart = vm["restart.write_restart"].as<bool>(); //! \param M_write_restart (double) Option on writing restart files
    if ( fmod(restart_time_step,time_step) != 0)
    {
        std::cout << restart_time_step << " " << time_step << "\n";
        throw std::runtime_error("restart_time_step not an integer multiple of time_step");
    }

    
    //! Sets the value of some parameters relevant for ocean forcing (turning angle, surface drag coef, basal drag )
    ocean_turning_angle_rad = 0.; //! \param ocean_turning_angle_rad (double) Ocean turning angle [rad]
    if (vm["dynamics.use_coriolis"].as<bool>())
        ocean_turning_angle_rad = (PI/180.)*vm["dynamics.oceanic_turning_angle"].as<double>();
    ridging_exponent = vm["dynamics.ridging_exponent"].as<double>(); //! \param ridging_exponent (double) Ridging exponent
    
    quad_drag_coef_water = vm["dynamics.quad_drag_coef_water"].as<double>(); //! \param quad_drag_coef_water (double) Quadratic ocean drag coefficient
    
    basal_k2 = vm["dynamics.Lemieux_basal_k2"].as<double>(); //! \param basal_k2 (double) Free parameter that determines the maximum basal stress (ice keels scheme of Lemieux et al., 2016)
    basal_u_0 = vm["dynamics.Lemieux_basal_u_0"].as<double>(); //! \param basal_u_0 (double) "Small velocity" parameter (ice keels scheme of Lemieux et al., 2016)
    basal_Cb = vm["dynamics.Lemieux_basal_Cb"].as<double>(); //! \param basal_Cb (double) Basal stress coefficient (ice keels scheme of Lemieux et al., 2016)
    
    
    //! Sets the values of parameters related to healing
    time_relaxation_damage = vm["dynamics.time_relaxation_damage"].as<double>()*days_in_sec; //! \param time_relaxation_damage (double) Characteristic healing time [s]
    deltaT_relaxation_damage = vm["dynamics.deltaT_relaxation_damage"].as<double>(); //! \param deltaT_relaxation_damage (double) Difference between the air and ocean temperature considered to set the characteristic time of damage [C]
    

    //! Sets the minimum and maximum thickness of thin ice
    h_thin_max = vm["thermo.h_thin_max"].as<double>(); //! \param h_thin_max (double) Maximum thickness of thin ice [m]
    h_thin_min = vm["thermo.h_thin_min"].as<double>(); //! \param h_thin_min (double) Minimum thickness of thin ice [m]
    
    
    //! Sets mechanical parameter values
    compr_strength = vm["dynamics.compr_strength"].as<double>(); //! \param compr_strength (double) Maximum compressive strength [N/m2]
    tract_coef = vm["dynamics.tract_coef"].as<double>(); //! \param tract_coef (double) Coefficient to set the maximum tensile strenght as a function of the cohesive strength
    // scale_coef is now set after initialising the mesh
    // scale_coef = vm["dynamics.scale_coef"].as<double>();
    alea_factor = vm["dynamics.alea_factor"].as<double>(); //! \param alea_factor (double) Sets the width of the distribution of cohesion
    cfix = vm["dynamics.cfix"].as<double>(); //! \param cfix (double) Fixed part of the cohesion [Pa]
    // C_fix    = cfix*scale_coef;          // C_fix;...  : cohesion (mohr-coulomb) in MPa (40000 Pa)
    // C_alea   = alea_factor*C_fix;        // C_alea;... : alea sur la cohesion (Pa)
    tan_phi = vm["dynamics.tan_phi"].as<double>(); //! \param tan_phi (double) Internal friction coefficient (mu)
    

    //! Sets options on the thermodynamics scheme
    if ( vm["thermo.newice_type"].as<int>() == 4 )
        M_ice_cat_type = setup::IceCategoryType::THIN_ICE; //! \param M_ice_cat_type (int) Option on using ice categories (thin ice or "classic")
    else
        M_ice_cat_type = setup::IceCategoryType::CLASSIC;

    const boost::unordered_map<const std::string, setup::ThermoType> str2thermo = boost::assign::map_list_of
        ("zero-layer", setup::ThermoType::ZERO_LAYER)
        ("winton", setup::ThermoType::WINTON);
    M_thermo_type = str2thermo.find(vm["setup.thermo-type"].as<std::string>())->second; //! \param M_thermo_type (string) Option on the thermodynamic scheme (Winton or zero-layer model)
    LOG(DEBUG)<<"ThermoType= "<< (int)M_thermo_type <<"\n";

    
    //! Sets options on the atmospheric and ocean forcing, initialization of ice, type of dynamics, bathymetry and on the use of nested meshes
    const boost::unordered_map<const std::string, setup::AtmosphereType> str2atmosphere = boost::assign::map_list_of
        ("constant", setup::AtmosphereType::CONSTANT)
        ("asr", setup::AtmosphereType::ASR)
        ("erai", setup::AtmosphereType::ERAi)
        ("ec", setup::AtmosphereType::EC)
        ("ec2", setup::AtmosphereType::EC2)
        ("ec_erai", setup::AtmosphereType::EC_ERAi)
        ("cfsr", setup::AtmosphereType::CFSR)
        ("cfsr_hi", setup::AtmosphereType::CFSR_HI);
    M_atmosphere_type = str2atmosphere.find(vm["setup.atmosphere-type"].as<std::string>())->second; //! \param M_atmosphere_type (string) Option on the type of atm. forcing (constant or reanalyses)
        switch(M_atmosphere_type){
        case setup::AtmosphereType::CONSTANT:   quad_drag_coef_air = vm["dynamics.ASR_quad_drag_coef_air"].as<double>(); break;
        case setup::AtmosphereType::ASR:        quad_drag_coef_air = vm["dynamics.ASR_quad_drag_coef_air"].as<double>(); break;
        case setup::AtmosphereType::CFSR_HI:
        case setup::AtmosphereType::CFSR:       quad_drag_coef_air = vm["dynamics.CFSR_quad_drag_coef_air"].as<double>(); break;
        case setup::AtmosphereType::ERAi:       quad_drag_coef_air = vm["dynamics.ERAi_quad_drag_coef_air"].as<double>(); break;
        case setup::AtmosphereType::EC:
        case setup::AtmosphereType::EC2:
        case setup::AtmosphereType::EC_ERAi:
                    quad_drag_coef_air = vm["dynamics.ECMWF_quad_drag_coef_air"].as<double>(); break;
        default:        std::cout << "invalid wind forcing"<<"\n";throw std::logic_error("invalid wind forcing");
    }
    LOG(DEBUG)<<"AtmosphereType= "<< (int)M_atmosphere_type <<"\n";

    M_use_nesting= vm["nesting.use_nesting"].as<bool>(); //! \param M_use_nesting (boolean) Option on the use of nested model meshes
    
    if (M_use_nesting)
    {
        M_use_ocean_nesting = vm["nesting.use_ocean_nesting"].as<bool>();
        M_nest_outer_mesh   = vm["nesting.outer_mesh"].as<std::string>();
        M_nest_inner_mesh   = vm["nesting.inner_mesh"].as<std::string>();
        M_nest_method       = vm["nesting.method"].as<std::string>();
        M_nudge_function    = vm["nesting.nudge_function"].as<std::string>();
        M_nudge_timescale   = vm["nesting.nudge_timescale"].as<double>();
        M_nudge_lengthscale = vm["nesting.nudge_lengthscale"].as<double>();
        M_nest_dynamic_vars = vm["nesting.nest_dynamic_vars"].as<bool>();
    }

    const boost::unordered_map<const std::string, setup::OceanType> str2ocean = boost::assign::map_list_of
        ("constant", setup::OceanType::CONSTANT)
        ("topaz", setup::OceanType::TOPAZR)
        ("topaz_atrest", setup::OceanType::TOPAZR_atrest)
        ("topaz_forecast", setup::OceanType::TOPAZF)
        ("topaz_altimeter", setup::OceanType::TOPAZR_ALTIMETER);
    M_ocean_type = str2ocean.find(vm["setup.ocean-type"].as<std::string>())->second; //! \param M_ocean_type (string) Option on the type of ocean forcing (constant or Topaz options)
    LOG(DEBUG) <<"OCEANTYPE= "<< (int)M_ocean_type <<"\n";

    const boost::unordered_map<const std::string, setup::IceType> str2conc = boost::assign::map_list_of
        ("constant", setup::IceType::CONSTANT)
        ("constant_partial", setup::IceType::CONSTANT_PARTIAL)
        ("target", setup::IceType::TARGET)
        ("topaz", setup::IceType::TOPAZ4)
        ("topaz_forecast", setup::IceType::TOPAZ4F)
        ("topaz_forecast_amsr2", setup::IceType::TOPAZ4FAMSR2)
        ("topaz_forecast_amsr2_osisaf", setup::IceType::TOPAZ4FAMSR2OSISAF)
        ("topaz_forecast_amsr2_osisaf_nic", setup::IceType::TOPAZ4FAMSR2OSISAFNIC)
        ("topaz_forecast_amsr2_osisaf_nic_weekly", setup::IceType::TOPAZ4FAMSR2OSISAFNICWEEKLY)
        ("amsre", setup::IceType::AMSRE)
        ("amsr2", setup::IceType::AMSR2)
        ("osisaf", setup::IceType::OSISAF)
        ("piomas", setup::IceType::PIOMAS)
        ("cs2_smos", setup::IceType::CS2_SMOS)
        ("cs2_smos_amsr2", setup::IceType::CS2_SMOS_AMSR2)
        ("smos", setup::IceType::SMOS)
        ("topaz_osisaf_icesat", setup::IceType::TOPAZ4OSISAFICESAT);
    M_ice_type = str2conc.find(vm["setup.ice-type"].as<std::string>())->second;
    LOG(DEBUG) <<"ICETYPE= "<< (int)M_ice_type <<"\n";

    const boost::unordered_map<const std::string, setup::DynamicsType> str2dynamics = boost::assign::map_list_of
        ("default", setup::DynamicsType::DEFAULT)
        ("no_motion", setup::DynamicsType::NO_MOTION)
        ("free_drift", setup::DynamicsType::FREE_DRIFT);
    M_dynamics_type = str2dynamics.find(vm["setup.dynamics-type"].as<std::string>())->second; //! \param M_dynamics_type (string) Option on the type of dynamics (default, no motion or freedrift)
    LOG(DEBUG) <<"DYNAMICSTYPE= "<< (int)M_dynamics_type <<"\n";

    const boost::unordered_map<const std::string, setup::BathymetryType> str2bathymetry = boost::assign::map_list_of
        ("constant", setup::BathymetryType::CONSTANT)
        ("etopo", setup::BathymetryType::ETOPO);
    M_bathymetry_type = str2bathymetry.find(vm["setup.bathymetry-type"].as<std::string>())->second; //! \param M_bathymetry_type (string) Option on the type of bathymetry (constant or ETOPO)
    LOG(DEBUG) <<"BATHYMETRYTYPE= "<< (int) M_bathymetry_type <<"\n";

    const boost::unordered_map<const std::string, setup::BasalStressType> str2basal_stress= boost::assign::map_list_of
        ("none", setup::BasalStressType::NONE)
        ("lemieux", setup::BasalStressType::LEMIEUX)
        ("bouillon", setup::BasalStressType::BOUILLON);
    M_basal_stress_type = str2basal_stress.find(vm["setup.basal_stress-type"].as<std::string>())->second; //! \param M_basal_stress_type (string) Option on the type of basal stress (none, from Lemieux et al., 2016 or from Bouillon)
    LOG(DEBUG) <<"BASALSTRESTYPE= "<< (int) M_basal_stress_type <<"\n";



    //! Sets the type and format of the mesh and the mesh filename
    const boost::unordered_map<const std::string, setup::MeshType> str2mesh = boost::assign::map_list_of
        ("from_unref", setup::MeshType::FROM_UNREF)
        ("from_split", setup::MeshType::FROM_SPLIT);
    M_mesh_type = str2mesh.find(vm["mesh.type"].as<std::string>())->second; //! \param M_mesh_type (string) Mesh type (unref or split)
    LOG(DEBUG) <<"MESHTYPE= "<< (int) M_mesh_type <<"\n";

    M_mesh_basename = vm["mesh.filename"].as<std::string>(); //! \param M_mesh_basename (string) Mesh filename
    M_mesh_filename = (boost::format( "%1%/%2%" ) // \param M_mesh_filename (string) Mesh filename (with path)
                       % Environment::nextsimMeshDir().string()
            % M_mesh_basename
            ).str();
    M_partitioned_mesh_filename = (boost::format( "%1%/par%2%%3%" )
            % M_export_path
            % M_comm.size()
            % M_mesh_basename
            ).str();
    M_mesh_fileformat = vm["mesh.partitioner-fileformat"].as<std::string>(); //! \param M_mesh_fileformat (string) Format of the partitioned mesh file (used if mesh.partitioner-space=="disk")
    M_mesh.setOrdering("bamg");
    
    
    //! Sets options on the use of moorings
    M_use_moorings =  vm["moorings.use_moorings"].as<bool>(); //! \param M_use_moorings (boolean) Option on the use of moorings
    M_moorings_snapshot =  vm["moorings.snapshot"].as<bool>(); //! \param M_moorings_snapshot (boolean) Option on outputing snapshots of mooring records
    M_moorings_parallel_output =  vm["moorings.parallel_output"].as<bool>(); //! \param M_moorings_parallel_output (boolean) Option on parallel outputs
    const boost::unordered_map<const std::string, GridOutput::fileLength> str2mooringsfl = boost::assign::map_list_of
        ("inf", GridOutput::fileLength::inf)
        ("daily", GridOutput::fileLength::daily)
        ("weekly", GridOutput::fileLength::weekly)
        ("monthly", GridOutput::fileLength::monthly)
        ("yearly", GridOutput::fileLength::yearly);
    M_moorings_file_length = str2mooringsfl.find(vm["moorings.file_length"].as<std::string>())->second;
        //! \param M_moorings_file_length (string) Length (in time) of the mooring output file
        //! (set according to daily, weekly, monthly or yearly outputs or to the "unlimited" option.)
    M_moorings_false_easting = vm["moorings.false_easting"].as<bool>();
        //! \param M_moorings_false_easting (boolean) Orientation of output vectors (true: components relative to output grid; false: or north/east components)
    M_moorings_averaging_period = 0.;//! \param M_moorings_averaging_period (double) averaging period in days. Zero if outputting snapshots. Used in netcdf metadata
    if(!M_moorings_snapshot)
        M_moorings_averaging_period = vm["moorings.output_timestep"].as<double>();
    
    //! Sets the type of partitioner and partition space
    const boost::unordered_map<const std::string, mesh::Partitioner> str2partitioner = boost::assign::map_list_of
        ("chaco", mesh::Partitioner::CHACO)
        ("metis", mesh::Partitioner::METIS);
    M_partitioner = str2partitioner.find(vm["mesh.partitioner"].as<std::string>())->second; //! \param M_partitioner (string) Sets the type of partioner (CHACO or METIS)
    
    const boost::unordered_map<const std::string, mesh::PartitionSpace> str2partitionspace = boost::assign::map_list_of
        ("memory", mesh::PartitionSpace::MEMORY)
        ("disk", mesh::PartitionSpace::DISK);

    M_partition_space = str2partitionspace.find(vm["mesh.partitioner-space"].as<std::string>())->second; //! \param M_partition_space (string) Sets the space for partitions (memory or disk)
    //! - Set the drifter options
    //  NB needs to be done before readRestart()
    this->initDrifterOpts();

}//initOptAndParam

    
//------------------------------------------------------------------------------------------------------
//! Creates a GMSH mesh grid.
//! !Does not seem to be used!
void
FiniteElement::createGMSHMesh(std::string const& geofilename)
{
    std::string gmshgeofile = (boost::format( "%1%/%2%" )
            % Environment::nextsimMeshDir().string()
            % geofilename
            ).str();

    if (fs::exists(gmshgeofile))
    {
        //std::cout<<"NOT FOUND " << fs::absolute( gmshgeofile ).string() <<"\n";
        std::ostringstream gmshstr;
        gmshstr << BOOST_PP_STRINGIZE( GMSH_EXECUTABLE )
                << " -" << 2 << " -part " << 1 << " -clmax " << vm["mesh.hsize"].as<double>() << " " << gmshgeofile;

        std::cout << "[Gmsh::generate] execute '" <<  gmshstr.str() << "'\n";
        auto err = ::system( gmshstr.str().c_str() );
    }
    else
    {
        std::cout << "Cannot found " << gmshgeofile <<"\n";
    }
}//createGMSHMesh

    
//------------------------------------------------------------------------------------------------------
//! Calculates the Jacobian Matrix Determinate:  measure of the normals of the element faces relative to each other.
//! Called by the flip(), measure() and shapeCoeff() functions.
//! \note
//! * This is used to calculate the finite element shape coefficient.
//! * The Jacobian an indicator of the distortion of the current mesh with respect to an undistorted mesh.

double
FiniteElement::jacobian(element_type const& element, mesh_type const& mesh) const
{
    std::vector<double> vertex_0 = mesh.nodes().find(element.indices[0])->second.coords;
    std::vector<double> vertex_1 = mesh.nodes().find(element.indices[1])->second.coords;
    std::vector<double> vertex_2 = mesh.nodes().find(element.indices[2])->second.coords;

    double jac = (vertex_1[0]-vertex_0[0])*(vertex_2[1]-vertex_0[1]);
    jac -= (vertex_2[0]-vertex_0[0])*(vertex_1[1]-vertex_0[1]);

    return  jac;
}//jacobian


double
FiniteElement::jacobian(element_type const& element, mesh_type const& mesh,
                        std::vector<double> const& um, double factor) const
{
    std::vector<double> vertex_0 = mesh.nodes().find(element.indices[0])->second.coords;
    std::vector<double> vertex_1 = mesh.nodes().find(element.indices[1])->second.coords;
    std::vector<double> vertex_2 = mesh.nodes().find(element.indices[2])->second.coords;

    for (int i=0; i<2; ++i)
    {
        vertex_0[i] += factor*um[element.indices[0]-1+i*(M_num_nodes)];
        vertex_1[i] += factor*um[element.indices[1]-1+i*(M_num_nodes)];
        vertex_2[i] += factor*um[element.indices[2]-1+i*(M_num_nodes)];
    }

    double jac = (vertex_1[0]-vertex_0[0])*(vertex_2[1]-vertex_0[1]);
    jac -= (vertex_2[0]-vertex_0[0])*(vertex_1[1]-vertex_0[1]);

    return  jac;
}//jacobian


//------------------------------------------------------------------------------------------------------
//! Calculates the Jacobian Matrix Determinate:  measure of the normals of the element faces relative to each other.
//! * This is used to calculate the finite element shape coefficient.
//! * The Jacobian an indicator of the distortion of the current mesh with respect to an undistorted mesh.
//! Called by the flip(), measure() and shapeCoeff() functions.
double
FiniteElement::jacobian(element_type const& element, mesh_type_root const& mesh) const
{
    std::vector<double> vertex_0 = mesh.nodes()[element.indices[0]-1].coords;
    std::vector<double> vertex_1 = mesh.nodes()[element.indices[1]-1].coords;
    std::vector<double> vertex_2 = mesh.nodes()[element.indices[2]-1].coords;

    double jac = (vertex_1[0]-vertex_0[0])*(vertex_2[1]-vertex_0[1]);
    jac -= (vertex_2[0]-vertex_0[0])*(vertex_1[1]-vertex_0[1]);

    return  jac;
}//jacobian


double
FiniteElement::jacobian(element_type const& element, mesh_type_root const& mesh,
                        std::vector<double> const& um, double factor) const
{
    std::vector<double> vertex_0 = mesh.nodes()[element.indices[0]-1].coords;
    std::vector<double> vertex_1 = mesh.nodes()[element.indices[1]-1].coords;
    std::vector<double> vertex_2 = mesh.nodes()[element.indices[2]-1].coords;

    int num_nodes = M_mesh_root.numNodes();

    for (int i=0; i<2; ++i)
    {
        vertex_0[i] += factor*um[element.indices[0]-1+i*(num_nodes)];
        vertex_1[i] += factor*um[element.indices[1]-1+i*(num_nodes)];
        vertex_2[i] += factor*um[element.indices[2]-1+i*(num_nodes)];
    }

    double jac = (vertex_1[0]-vertex_0[0])*(vertex_2[1]-vertex_0[1]);
    jac -= (vertex_2[0]-vertex_0[0])*(vertex_1[1]-vertex_0[1]);

    return  jac;
}//jacobian


//------------------------------------------------------------------------------------------------------
//! Calculates the length of the vertices of the triangular mesh elements.
//! Called by the minMaxSides() and minAngles() functions.
std::vector<double>
FiniteElement::sides(element_type const& element, mesh_type const& mesh) const
{
    std::vector<double> vertex_0 = mesh.nodes().find(element.indices[0])->second.coords;
    std::vector<double> vertex_1 = mesh.nodes().find(element.indices[1])->second.coords;
    std::vector<double> vertex_2 = mesh.nodes().find(element.indices[2])->second.coords;

    std::vector<double> side(3);

    side[0] = std::hypot(vertex_1[0]-vertex_0[0], vertex_1[1]-vertex_0[1]);
    side[1] = std::hypot(vertex_2[0]-vertex_1[0], vertex_2[1]-vertex_1[1]);
    side[2] = std::hypot(vertex_2[0]-vertex_0[0], vertex_2[1]-vertex_0[1]);

    return side;
}//sides


std::vector<double>
FiniteElement::sides(element_type const& element, mesh_type const& mesh,
                     std::vector<double> const& um, double factor) const
{
    std::vector<double> vertex_0 = mesh.nodes().find(element.indices[0])->second.coords;
    std::vector<double> vertex_1 = mesh.nodes().find(element.indices[1])->second.coords;
    std::vector<double> vertex_2 = mesh.nodes().find(element.indices[2])->second.coords;

    for (int i=0; i<2; ++i)
    {
        vertex_0[i] += factor*um[element.indices[0]-1+i*(M_num_nodes)];
        vertex_1[i] += factor*um[element.indices[1]-1+i*(M_num_nodes)];
        vertex_2[i] += factor*um[element.indices[2]-1+i*(M_num_nodes)];
    }

    std::vector<double> side(3);

    side[0] = std::hypot(vertex_1[0]-vertex_0[0], vertex_1[1]-vertex_0[1]);
    side[1] = std::hypot(vertex_2[0]-vertex_1[0], vertex_2[1]-vertex_1[1]);
    side[2] = std::hypot(vertex_2[0]-vertex_0[0], vertex_2[1]-vertex_0[1]);

    return side;
}//sides


std::vector<double>
FiniteElement::sides(element_type const& element, mesh_type_root const& mesh) const
{
    std::vector<double> vertex_0 = mesh.nodes()[element.indices[0]-1].coords;
    std::vector<double> vertex_1 = mesh.nodes()[element.indices[1]-1].coords;
    std::vector<double> vertex_2 = mesh.nodes()[element.indices[2]-1].coords;

    std::vector<double> side(3);

    side[0] = std::hypot(vertex_1[0]-vertex_0[0], vertex_1[1]-vertex_0[1]);
    side[1] = std::hypot(vertex_2[0]-vertex_1[0], vertex_2[1]-vertex_1[1]);
    side[2] = std::hypot(vertex_2[0]-vertex_0[0], vertex_2[1]-vertex_0[1]);

    return side;
}//sides


std::vector<double>
FiniteElement::sides(element_type const& element, mesh_type_root const& mesh,
                     std::vector<double> const& um, double factor) const
{
    std::vector<double> vertex_0 = mesh.nodes()[element.indices[0]-1].coords;
    std::vector<double> vertex_1 = mesh.nodes()[element.indices[1]-1].coords;
    std::vector<double> vertex_2 = mesh.nodes()[element.indices[2]-1].coords;

    int num_nodes = M_mesh_root.numNodes();

    for (int i=0; i<2; ++i)
    {
        vertex_0[i] += factor*um[element.indices[0]-1+i*(num_nodes)];
        vertex_1[i] += factor*um[element.indices[1]-1+i*(num_nodes)];
        vertex_2[i] += factor*um[element.indices[2]-1+i*(num_nodes)];
    }

    std::vector<double> side(3);

    side[0] = std::hypot(vertex_1[0]-vertex_0[0], vertex_1[1]-vertex_0[1]);
    side[1] = std::hypot(vertex_2[0]-vertex_1[0], vertex_2[1]-vertex_1[1]);
    side[2] = std::hypot(vertex_2[0]-vertex_0[0], vertex_2[1]-vertex_0[1]);

    return side;
}//sides

    
//------------------------------------------------------------------------------------------------------
//! Calculates the maximum and minimum sizes of the vertices of the triangular mesh elements
//! Called by the rootMeshProcessing function.
std::vector<double>
FiniteElement::minMaxSide(mesh_type_root const& mesh) const
{
    std::vector<double> minmax(2);
    std::vector<double> all_min_side(mesh.numTriangles());
    std::vector<double> all_max_side(mesh.numTriangles());

    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        auto side = this->sides(*it,mesh);
        all_min_side[cpt] = *std::min_element(side.begin(),side.end());
        all_max_side[cpt] = *std::max_element(side.begin(),side.end());
        ++cpt;
    }

    // minmax[0] = *std::min_element(all_min_side.begin(),all_min_side.end());
    // minmax[1] = *std::max_element(all_max_side.begin(),all_max_side.end());

    minmax[0] = std::accumulate(all_min_side.begin(),all_min_side.end(),0.)/(all_min_side.size());
    minmax[1] = std::accumulate(all_max_side.begin(),all_max_side.end(),0.)/(all_max_side.size());

    return minmax;
}//minMaxSide

    
//------------------------------------------------------------------------------------------------------
//! Calculates the minimum angle of each triangular mesh element.
//! Called by the minAngle() function.
template<typename FEMeshType>
double
FiniteElement::minAngles(element_type const& element, FEMeshType const& mesh) const
{
    std::vector<double> side = this->sides(element,mesh);
    //std::for_each(side.begin(), side.end(), [&](double& f){ f = 1000.*f; });
    std::sort(side.begin(),side.end());
    double minang = std::acos( (std::pow(side[1],2.) + std::pow(side[2],2.) - std::pow(side[0],2.) )/(2*side[1]*side[2]) );
    minang = minang*45.0/std::atan(1.0);

    return minang;
}//minAngles


template<typename FEMeshType>
double
FiniteElement::minAngles(element_type const& element, FEMeshType const& mesh,
                         std::vector<double> const& um, double factor) const
{
    std::vector<double> side = this->sides(element,mesh,um,factor);
    //std::for_each(side.begin(), side.end(), [&](double& f){ f = 1000.*f; });
    std::sort(side.begin(),side.end());
    double minang = std::acos( (std::pow(side[1],2.) + std::pow(side[2],2.) - std::pow(side[0],2.) )/(2*side[1]*side[2]) );
    minang = minang*45.0/std::atan(1.0);

    return minang;
}//minAngles

    
//------------------------------------------------------------------------------------------------------
//! Finds the minimum angle within all mesh elements: condition for the remeshing.
//! Called by the init() function.
template<typename FEMeshType>
double
FiniteElement::minAngle(FEMeshType const& mesh) const
{
    std::vector<double> all_min_angle(mesh.numTriangles());

    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        all_min_angle[cpt] = this->minAngles(*it,mesh);
        ++cpt;
    }

    double min_angle = *std::min_element(all_min_angle.begin(),all_min_angle.end());
    //return min_angle;
    return boost::mpi::all_reduce(M_comm, min_angle, boost::mpi::minimum<double>());
}//minAngle


template<typename FEMeshType>
double
FiniteElement::minAngle(FEMeshType const& mesh, std::vector<double> const& um, double factor, bool root) const
{
    std::vector<double> all_min_angle(mesh.numTriangles());

    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        all_min_angle[cpt] = this->minAngles(*it,mesh,um,factor);
        ++cpt;
    }

    double min_angle = *std::min_element(all_min_angle.begin(),all_min_angle.end());

    double res_value = min_angle;

    if (!root)
    {
        res_value = boost::mpi::all_reduce(M_comm, min_angle, boost::mpi::minimum<double>());
    }

    return res_value;
}//minAngle


//------------------------------------------------------------------------------------------------------
//! Detects a flipped element.
//! Called by the regrid() function.
template<typename FEMeshType>
bool
FiniteElement::flip(FEMeshType const& mesh, std::vector<double> const& um, double factor) const
{
    std::vector<double> area(mesh.numTriangles());

    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        area[cpt] = this->jacobian(*it,mesh,um,factor);
        ++cpt;
    }

    double minarea = *std::min_element(area.begin(),area.end());
    double maxarea = *std::max_element(area.begin(),area.end());

    return ((minarea <= 0.) && (maxarea >= 0.));
}//flip


//------------------------------------------------------------------------------------------------------
//! Calculates the mean resolution of the non-structured triangular elements mesh grid.
//! Called by the rootMeshProcessing() function.
double
FiniteElement::resolution(mesh_type_root const& mesh) const
{
    std::vector<double> all_min_measure(mesh.numTriangles());

    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        all_min_measure[cpt] = this->measure(*it,mesh);
        ++cpt;
    }

    double resol = std::accumulate(all_min_measure.begin(),all_min_measure.end(),0.)/(all_min_measure.size());
    resol = std::pow(resol,0.5);

    return resol;
}//resolution


//------------------------------------------------------------------------------------------------------
//! Calculates the minimum height (i.e., dimension) of each triangular mesh element.
//! Called by the rootMeshProcessing() function.
std::vector<double>
FiniteElement::hminVertices(mesh_type_root const& mesh, BamgMesh const* bamg_mesh) const
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
                measure[j] = this->measure(mesh.triangles()[elt_num],mesh);
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

    
//------------------------------------------------------------------------------------------------------
//! Returns the maximum height (i.e., dimension) of all mesh elements
//! Called by the rootMeshProcessing() function.
std::vector<double>
FiniteElement::hmaxVertices(mesh_type_root const& mesh, BamgMesh const* bamg_mesh) const
{
    std::vector<double> hmax = this->hminVertices(mesh,bamg_mesh);

    std::for_each(hmax.begin(), hmax.end(), [&](double& f){ f = 1.2*f; });

    return hmax;
}//hmaxVertices

    
//------------------------------------------------------------------------------------------------------
//! Calculates the area of triangular mesh elements.
//! Called by the advect() and other functions.
template<typename FEMeshType>
double
FiniteElement::measure(element_type const& element, FEMeshType const& mesh) const
{
    return (1./2)*std::abs(jacobian(element,mesh));
}//measure


//------------------------------------------------------------------------------------------------------
//! Calculates the area of triangular mesh elements.
//! Called by the advect() and other functions.
template<typename FEMeshType>
double
FiniteElement::measure(element_type const& element, FEMeshType const& mesh,
                       std::vector<double> const& um, double factor) const
{
    return (1./2)*std::abs(jacobian(element,mesh,um,factor));
}//measure

    
//------------------------------------------------------------------------------------------------------
//! Calculates finite element shape coefficients.
//! Called by the FETensors() function.
std::vector<double>
FiniteElement::shapeCoeff(element_type const& element, mesh_type const& mesh) const
{
    std::vector<double> x(3);
    std::vector<double> y(3);

    for (int i=0; i<3; ++i)
    {
        x[i] = mesh.nodes().find(element.indices[i])->second.coords[0];
        y[i] = mesh.nodes().find(element.indices[i])->second.coords[1];
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

    
//------------------------------------------------------------------------------------------------------
//! Calculates finite element shape coefficients.
//! Called by the FETensors() function.
std::vector<double>
FiniteElement::shapeCoeff(element_type const& element, mesh_type_root const& mesh) const
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

    
//------------------------------------------------------------------------------------------------------
//! Interpolates hminVertices and hmaxVertices onto the current mesh.
//! Called by the rootMeshProcessing() and regrid() functions.
void
FiniteElement::interpVertices()
{
    chrono.restart();
    LOG(DEBUG) <<"Interpolate hminVertices starts\n";
    
    // NODAL INTERPOLATION
    int init_num_nodes = M_mesh_init_root.numNodes();

    std::vector<double> interp_Vertices_in(2*init_num_nodes);

    double* interp_Vertices_out;

    for (int i=0; i<init_num_nodes; ++i)
    {
        interp_Vertices_in[2*i]   = M_hminVertices[i];
        interp_Vertices_in[2*i+1] = M_hmaxVertices[i];
    }

    InterpFromMeshToMesh2dx(&interp_Vertices_out,
                            &M_mesh_init_root.indexTr()[0],&M_mesh_init_root.coordX()[0],&M_mesh_init_root.coordY()[0],
                            M_mesh_init_root.numNodes(),M_mesh_init_root.numTriangles(),
                            &interp_Vertices_in[0],
                            M_mesh_init_root.numNodes(),2,
                            &M_mesh_root.coordX()[0],&M_mesh_root.coordY()[0],M_mesh_root.numNodes(),
                            false);

    bamgopt->hminVertices = new double[M_mesh_root.numNodes()];
    bamgopt->hmaxVertices = new double[M_mesh_root.numNodes()];

    for (int i=0; i<M_mesh_root.numNodes(); ++i)
    {
        bamgopt->hminVertices[i] = interp_Vertices_out[2*i];
        bamgopt->hmaxVertices[i] = interp_Vertices_out[2*i+1];
    }

    xDelete<double>(interp_Vertices_out);
    LOG(DEBUG) <<"Interpolate hmin done in "<< chrono.elapsed() <<"s\n";
}//interpVertices

    
//------------------------------------------------------------------------------------------------------
//! Interpolates hminVertices and hmaxVertices onto the current mesh.
//! Called by the distributedMeshProcessing() function.
void
FiniteElement::gatherSizes()
{

    std::vector<int> fesizes_local(4);
    fesizes_local = {M_local_ndof, M_num_nodes, M_local_nelements, M_num_elements};

    std::vector<int> fesizes;

    boost::mpi::all_gather(M_comm, &fesizes_local[0], 4, fesizes);

    M_sizes_nodes.resize(M_comm.size());
    M_sizes_nodes_with_ghost.resize(M_comm.size());
    M_sizes_elements.resize(M_comm.size());
    M_sizes_elements_with_ghost.resize(M_comm.size());

    for (int i=0; i<M_comm.size(); ++i)
    {
        // nodes
        M_sizes_nodes[i] = fesizes[4*i];
        M_sizes_nodes_with_ghost[i] = fesizes[4*i+1];

        // elements
        M_sizes_elements[i] = fesizes[4*i+2];
        M_sizes_elements_with_ghost[i] = fesizes[4*i+3];
    }

    std::vector<int> sizes_nodes = M_sizes_nodes_with_ghost;

    if (M_rank == 0)
    {
        int out_size_nodes = std::accumulate(sizes_nodes.begin(),sizes_nodes.end(),0);
        M_id_nodes.resize(out_size_nodes);

        boost::mpi::gatherv(M_comm, M_mesh.localDofWithGhostInit(), &M_id_nodes[0], sizes_nodes, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, M_mesh.localDofWithGhostInit(), 0);
    }

    std::vector<int> sizes_elements = M_sizes_elements_with_ghost;

    if (M_rank == 0)
    {
        int out_size_elements = std::accumulate(sizes_elements.begin(),sizes_elements.end(),0);
        M_id_elements.resize(out_size_elements);

        boost::mpi::gatherv(M_comm, M_mesh.trianglesIdWithGhost(), &M_id_elements[0], sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, M_mesh.trianglesIdWithGhost(), 0);
    }

    // -------------------------------------------------------------
    if (M_rank == 0)
    {
        M_rmap_nodes = M_mesh.mapNodes();
        M_rmap_elements = M_mesh.mapElements();
    }
    // -------------------------------------------------------------
}//gatherSizes

    
//------------------------------------------------------------------------------------------------------
//! Collects model variables and stores them into a single vector, interp_elt_in_local: called by the update() function,
//! before updating all variables after solving.
//! Called by the gatherFieldsElement() function.
void
FiniteElement::collectVariables(std::vector<double>& interp_elt_in_local, bool ghosts)
{
    // ELEMENT INTERPOLATION With Cavities
    int nb_var_element = M_nb_var_element;

    int num_elements = M_local_nelements;
    if (ghosts)
    {
        num_elements = M_num_elements;
    }

    interp_elt_in_local.resize(nb_var_element*num_elements);
    M_interp_method.resize(nb_var_element); // 0 for non conservative method, 1 for conservative method (for variables defined in terms of blabla/per unit area)
    M_diffusivity_parameters.resize(nb_var_element); // 0 for non added diffusion, positive value for active diffusion in [m^2/s] (only non conservative implementation available)

    int tmp_nb_var=0;
    for (int i=0; i<num_elements; ++i)
    {
        tmp_nb_var=0;

        // concentration
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_conc[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // thickness
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_thick[i];
        M_interp_method[tmp_nb_var] = 1; // crash if equal 1
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // snow thickness
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_snow_thick[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // integrated_stress1
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // integrated_stress2
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i+1];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // integrated_stress3
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i+2];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // damage
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_damage[i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // ridge_ratio
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_ridge_ratio[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // random_number
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_random_number[i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // sea surface salinity
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sss[i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=vm["thermo.diffusivity_sss"].as<double>();
        tmp_nb_var++;

		// sea surface temperature
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sst[i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=vm["thermo.diffusivity_sst"].as<double>();
        tmp_nb_var++;

        // Ice temperature
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_tice[0][i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = ( M_tice[1][i] - physical::mu*physical::si*physical::Lf/(physical::C*M_tice[1][i]) ) * M_thick[i]; // (39) times volume with f1=1
            M_interp_method[tmp_nb_var] = 1;
            M_diffusivity_parameters[tmp_nb_var]=0.;
            tmp_nb_var++;

            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = ( M_tice[2][i] ) * M_thick[i]; // (39) times volume with f1=0
            M_interp_method[tmp_nb_var] = 1;
            M_diffusivity_parameters[tmp_nb_var]=0.;
            tmp_nb_var++;
        }

        // thin ice thickness
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_h_thin[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // thin ice thickness
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_conc_thin[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // snow on thin ice
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_hs_thin[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // Ice surface temperature for thin ice
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_tsurf_thin[i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        if(tmp_nb_var>nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }
    }
}//collectVariables

    
//------------------------------------------------------------------------------------------------------
//! Collects model variables and stores them into a single vector, interp_elt_in_local, for outputing.
//! Called by the gatherFieldsElementIO() function.
    void
FiniteElement::collectVariablesIO(std::vector<double>& interp_elt_in_local, bool ghosts, bool thin_ice)
{
    int nb_var_element = M_nb_var_element;
    if (!thin_ice)
    {
        nb_var_element -= 4;
    }

    if (vm["output.save_diagnostics"].as<bool>())
    {
        nb_var_element += 7;
    }

    int num_elements = M_local_nelements;
    if (ghosts)
    {
        num_elements = M_num_elements;
    }

    interp_elt_in_local.resize(nb_var_element*num_elements);

    int tmp_nb_var=0;
    for (int i=0; i<num_elements; ++i)
    {
        tmp_nb_var=0;

        // concentration
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_conc[i];
        tmp_nb_var++;

        // thickness
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_thick[i];
        tmp_nb_var++;

        // snow thickness
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_snow_thick[i];
        tmp_nb_var++;

        // integrated_stress1
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i];
        tmp_nb_var++;

        // integrated_stress2
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i+1];
        tmp_nb_var++;

        // integrated_stress3
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i+2];
        tmp_nb_var++;

        // damage
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_damage[i];
        tmp_nb_var++;

        // ridge_ratio
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_ridge_ratio[i]/**M_thick[i]*/;
        tmp_nb_var++;

        // random_number
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_random_number[i];
        tmp_nb_var++;

        // diffusivity_sss
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sss[i];
        tmp_nb_var++;

        // diffusivity_sst
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sst[i];
        tmp_nb_var++;

        // Ice temperature
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_tice[0][i];
        tmp_nb_var++;

        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_tice[1][i];
            tmp_nb_var++;

            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_tice[2][i];
            tmp_nb_var++;
        }

        if (thin_ice)
        {
            // thin ice thickness
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_h_thin[i];
            tmp_nb_var++;

            // thin ice thickness
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_conc_thin[i];
            tmp_nb_var++;

            // snow on thin ice
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_hs_thin[i];
            tmp_nb_var++;

            // Ice surface temperature for thin ice
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_tsurf_thin[i];
            tmp_nb_var++;
        }

        if (vm["output.save_diagnostics"].as<bool>())
        {
            // Qatm
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = D_Qa[i];
            tmp_nb_var++;

            // Qsw
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = D_Qsw[i];
            tmp_nb_var++;

            // Qlw
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = D_Qlw[i];
            tmp_nb_var++;

            // Qsh
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = D_Qsh[i];
            tmp_nb_var++;

            // Qlh
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = D_Qlh[i];
            tmp_nb_var++;

            // Qocean
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = D_Qo[i];
            tmp_nb_var++;

            // Saltflux
            interp_elt_in_local[nb_var_element*i+tmp_nb_var] = D_delS[i];
            tmp_nb_var++;
        }

        if(tmp_nb_var>nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }
    }
}//collectVariablesIO


//------------------------------------------------------------------------------------------------------
//! Interpolates variables (other than velocities and displacements) onto the mesh grid once updated.
//! Called by the updated() function, after the advect() function.
void
FiniteElement::redistributeVariables(std::vector<double> const& out_elt_values, bool check_conc)
{
    int nb_var_element = M_nb_var_element;

    int tmp_nb_var=0;

    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_nb_var=0;

        // concentration
        M_conc[i] = std::max(0., out_elt_values[nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        // thickness
        M_thick[i] = std::max(0., out_elt_values[nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        // snow thickness
        M_snow_thick[i] = std::max(0., out_elt_values[nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        if (M_thick[i] != 0.)
        {
            // integrated_stress1
            M_sigma[3*i] = out_elt_values[nb_var_element*i+tmp_nb_var]/*/M_thick[i]*/;
            tmp_nb_var++;

            // integrated_stress2
            M_sigma[3*i+1] = out_elt_values[nb_var_element*i+tmp_nb_var]/*/M_thick[i]*/;
            tmp_nb_var++;

            // integrated_stress3
            M_sigma[3*i+2] = out_elt_values[nb_var_element*i+tmp_nb_var]/*/M_thick[i]*/;
            tmp_nb_var++;

            // damage
            M_damage[i] = std::max(0., std::min(1.,out_elt_values[nb_var_element*i+tmp_nb_var]));
            tmp_nb_var++;

            // ridge_ratio
            M_ridge_ratio[i] = std::max(0., std::min(1.,out_elt_values[nb_var_element*i+tmp_nb_var]/*/M_thick[i]*/));
            tmp_nb_var++;
        }
        else
        {
            // tmp_nb_var += 5;
            // integrated_stress1
            M_sigma[3*i] = 0.;
            tmp_nb_var++;

            // integrated_stress2
            M_sigma[3*i+1] = 0.;
            tmp_nb_var++;

            // integrated_stress3
            M_sigma[3*i+2] = 0.;
            tmp_nb_var++;

            // damage
            M_damage[i] = 0.;
            tmp_nb_var++;

            // damage
            M_ridge_ratio[i] = 0.;
            tmp_nb_var++;
        }

        // random_number
        M_random_number[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // SSS
        M_sss[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // SST
        M_sst[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // Ice temperature
        M_tice[0][i] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            if(M_thick[i]>0.)
            {
                double tmp = out_elt_values[nb_var_element*i+tmp_nb_var]/M_thick[i];
                M_tice[1][i] = 0.5*( tmp - std::sqrt(tmp*tmp + 4*physical::mu*physical::si*physical::Lf/physical::C) ); // (38) divided with volume with f1=1
                tmp_nb_var++;

                M_tice[2][i] = out_elt_values[nb_var_element*i+tmp_nb_var]/M_thick[i]; // (40) divided with volume with f1=0
                tmp_nb_var++;
            }
            else
            {
                M_tice[1][i] = -physical::mu*M_sss[i];
                tmp_nb_var++;

                M_tice[2][i] = -physical::mu*M_sss[i];
                tmp_nb_var++;
            }
        }

        // thin ice thickness
        M_h_thin[i] = std::max(0., out_elt_values[nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        // thin ice concentration
        M_conc_thin[i] = std::max(0., out_elt_values[nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        // snow on thin ice
        M_hs_thin[i] = std::max(0., out_elt_values[nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        // Ice surface temperature for thin ice
        M_tsurf_thin[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        if(tmp_nb_var!=nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }

        if (check_conc)
        {
            M_conc[i]=(M_conc[i]>1.) ? 1. : M_conc[i];
            double conc_thin_tmp = ( (M_conc[i]+M_conc_thin[i])>1.) ? 1.-M_conc[i] : M_conc_thin[i];
            double h_thin_tmp = 0.;

            if(M_conc_thin[i]>0.)
            {
                h_thin_tmp = M_h_thin[i]*conc_thin_tmp/M_conc_thin[i];
            }

            M_thick[i] += M_h_thin[i] - h_thin_tmp;
            M_h_thin[i] = h_thin_tmp;
            M_conc_thin[i] = conc_thin_tmp;
        }
    }
}//redistributeVariables

    
//------------------------------------------------------------------------------------------------------
//! Gets the names of the variables in the run restart file, takes a list of names and for each name,
//! adds a pointer to the appropriate vector, adds the number of components in the variable to another vector
//! (this is usually 1, but can be 3 eg M_sigma).
//! These outputs are then used in loops in collectVariablesIO and scatterFieldsElementIO.
//! Called from the readRestart() function.
void
FiniteElement::getVariablesIO(
        std::vector<std::vector<double>*> &data,
        std::vector<int> &num_components,
        std::vector<std::string> const &names)
{


    //!1st set pointers to the data requested in "names"
    for(auto it=names.begin(); it!=names.end(); it++)
    {
        LOG(DEBUG)<<"collectVariablesIO: adding "<<*it<<"\n";
        int num_comp = 1;
        if(*it=="M_conc")
        {
            // concentration
            data.push_back(&M_conc);
        }
        else if(*it=="M_thick")
        {
            // thickness
            data.push_back(&M_thick);
        }
        else if(*it=="M_snow_thick")
        {
            // snow thickness
            data.push_back(&M_snow_thick);
        }
        else if(*it=="M_sigma")
        {
            // stress
            data.push_back(&M_sigma);
            num_comp = 3;
        }
        else if(*it=="M_damage")
        {
            // damage
            data.push_back(&M_damage);
        }
        else if(*it=="M_ridge_ratio")
        {
            // damage
            data.push_back(&M_ridge_ratio);
        }
        else if(*it == "M_random_number")
        {
            // random_number
            data.push_back(&M_random_number);
        }
        else if(*it == "M_sss")
        {
            // SSS
            data.push_back(&M_sss);
        }
        else if(*it == "M_sst")
        {
            // SST
            data.push_back(&M_sst);
        }
        else if(*it == "M_tice_0")
        {
            // M_tice[0] - Ice temperature
            data.push_back(&(M_tice[0]));
        }
        else if(*it == "M_tice_1")
        {
            // M_tice[1] - Ice temperature
            data.push_back(&(M_tice[1]));
        }
        else if(*it == "M_tice_2")
        {
            // M_tice[2] - Ice temperature
            data.push_back(&(M_tice[2]));
        }
        else if(*it == "M_h_thin")
        {
            // thin ice thickness
            data.push_back(&M_h_thin);
        }
        else if(*it == "M_conc_thin")
        {
            // thin ice concentration
            data.push_back(&M_conc_thin);
        }
        else if(*it == "M_hs_thin")
        {
            // snow thickness on thin ice
            data.push_back(&M_hs_thin);
        }
        else if(*it == "M_tsurf_thin")
        {
            // surface temperature over thin ice
            data.push_back(&M_tsurf_thin);
        }
        else
            throw std::runtime_error("Unimplemented name: "+*it);

        //! 2nd, sets the number of components to loop over and resize the variables
        num_components.push_back(num_comp);
    }
}//getVariableIO


//------------------------------------------------------------------------------------------------------
//! Redistributes variables (parallel computing).
//! Called by function scatterFieldsElementIO().
//! * out_elt_values is vector containing all the variables to be redistributed (eg after scattering from root) into the individual variables (eg M_conc, M_thick,...)
//! * data is a vector of pointers to the variables to be assigned values from out_elt_values
//! * num_components is a vector with the number of components in each variable (usually 1, but can be 3 eg for M_sigma)
void
FiniteElement::redistributeVariablesIO(std::vector<double> const& out_elt_values,
        std::vector<std::vector<double>*> &data,
        std::vector<int> const& num_components)
{

    //! 1st, initializes the data
    int nb_var_element = 0;
    for(int j=0; j<data.size(); j++)
    {
        int num_comp = num_components[j];
        data[j]->assign(num_comp*M_num_elements, 0.);
        nb_var_element += num_comp;
    }

    //! 2nd, loops over the data and get their values from out_elt_values
    for (int i=0; i<M_num_elements; ++i)
    {
        int tmp_nb_var=0;
        for(int j=0; j<data.size(); j++)
        {
            int num_comp = num_components[j];
            for (int k=0; k<num_comp; k++)
            {
                (*(data[j]))[num_comp*i+k] = out_elt_values[nb_var_element*i+tmp_nb_var];
                tmp_nb_var++;
            }//loop over each component of variables
        }//loop over variables
        if(tmp_nb_var!=nb_var_element)
            throw std::logic_error("tmp_nb_var not equal to nb_var_element");
    }//loop over elements
}//redistributeVariablesIO

    
// Hotfix for issue #53 - we only have pure Lagrangian now.
//------------------------------------------------------------------------------------------------------
//! Performs the advection, using a Eulerian or ALE scheme
//! Called by the update() function.
#if 0
void
FiniteElement::advect(std::vector<double> const& interp_elt_in, std::vector<double>& interp_elt_out)
{
    M_comm.barrier();


    std::vector<double> UM_P = M_UM;

    interp_elt_out.resize(M_nb_var_element*M_num_elements);

    bool use_eulerian = false;
    bool use_lagrangian = false;
    bool use_ALE = false;
    int ALE_smoothing_step_nb = vm["numerics.ALE_smoothing_step_nb"].as<int>();
    if (vm["numerics.advection_scheme"].as<string>()=="Eulerian")
        // Diffusive eulerian case where M_UM is not changed and then =0.
        use_eulerian = true;
    else if (vm["numerics.advection_scheme"].as<string>()=="Lagrangian")
        // Purely Lagrangian case where M_UM is updated with M_VT
        use_lagrangian = true;
    else if (vm["numerics.advection_scheme"].as<string>()=="ALE")
    {
        // ALE case where M_UM is updated with a smoothed version of M_VT
        use_ALE = true;
        if (ALE_smoothing_step_nb<=0)
            throw std::runtime_error("numerics.ALE_smoothing_step_nb option should be >0");
    }
    else
        throw std::runtime_error("numerics.advection_scheme option should be Eulerian, Lagrangian or ALE");

    //change M_UM if not in Eulerian mode
    if(!use_eulerian)
    {
        if(use_lagrangian)
        {
            // pure Lagrangian case
            for (int nd=0; nd<M_UM.size(); ++nd)
            {
                M_UM[nd] += time_step*M_VT[nd];
            }
        }
        else
        {
            // ALE case - need to calculate the smoothed velocity
            std::vector<double> vt_root;
            std::vector<double> M_VT_smoothed;
            std::vector<double> M_VT_smoothed_root;

            this->gatherNodalField(M_VT,vt_root);

            if(M_rank == 0)
            {
                M_VT_smoothed_root = vt_root;
                std::vector<double> M_VT_tmp = M_VT_smoothed_root;

                // get the global number of nodes
                int num_nodes = M_mesh_root.numNodes();
                int Nd = bamgmesh_root->NodalConnectivitySize[1];

                for (int k=0; k<ALE_smoothing_step_nb; ++k)
                {
                    M_VT_tmp = M_VT_smoothed_root;

                    //for (int i=0; i<M_ndof; ++i)
                    for (int i=0; i<num_nodes; ++i)
                    {
                        int Nc;
                        double UM_x, UM_y;

                        if(M_mask_dirichlet_root[i]==false)
                        {
                            Nc = bamgmesh_root->NodalConnectivity[Nd*(i+1)-1];

                            UM_x = 0.;
                            UM_y = 0.;
                            for (int j=0; j<Nc; ++j)
                            {
                                UM_x += M_VT_tmp[bamgmesh_root->NodalConnectivity[Nd*i+j]-1];
                                UM_y += M_VT_tmp[bamgmesh_root->NodalConnectivity[Nd*i+j]-1+num_nodes];
                            }

                            M_VT_smoothed_root[i          ] = UM_x/Nc;
                            M_VT_smoothed_root[i+num_nodes] = UM_y/Nc;
                        }
                    }
                }
            }

            // M_comm.barrier();

            this->scatterNodalField(M_VT_smoothed_root,M_VT_smoothed);

            // M_comm.barrier();

            for (int nd=0; nd<M_UM.size(); ++nd)
            {
                M_UM[nd] += time_step*M_VT_smoothed[nd];
            }
        }//using ALE

        // set back the neumann nodes (open boundaries) at their position, the fluxes will be computed thanks to the convective velocity
        for (const int& nd : M_neumann_nodes)
        {
            M_UM[nd] = UM_P[nd];
        }
    }//not Eulerian

    LOG(DEBUG) <<"VT MIN= "<< *std::min_element(M_VT.begin(),M_VT.end()) <<"\n";
    LOG(DEBUG) <<"VT MAX= "<< *std::max_element(M_VT.begin(),M_VT.end()) <<"\n";

    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        /* some variables used for the advection*/
        double surface, surface_new;
        double integrated_variable;

        /* some variables used for the advection*/
        double x[3],y[3],x_new[3],y_new[3];
        int x_ind, y_ind, neighbour_int, vertex_1, vertex_2;
        double outer_fluxes_area[3], vector_edge[2], outer_vector[2], VC_middle[2], VC_x[3], VC_y[3];
        int fluxes_source_id[3];

        double neighbour_double;
        int other_vertex[3*2]={1,2 , 2,0 , 0,1};

        /*======================================================================
         * Update:
         * Ice and snow thickness, and concentration using a Lagrangian or an Eulerian scheme
         *======================================================================
         */

        /* convective velocity */
        for(int i=0;i<3;i++)
        {
            x_ind = (M_elements[cpt]).indices[i]-1;
            y_ind = (M_elements[cpt]).indices[i]-1+M_num_nodes;

            x[i] = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[0];
            y[i] = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[1];

            /* old and new positions of the mesh */
            x_new[i] = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[0]+M_UM[x_ind];
            y_new[i] = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[1]+M_UM[y_ind];
            x[i]     = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[0]+UM_P[x_ind];
            y[i]     = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[1]+UM_P[y_ind];

            // convective velocity
            // - this is zero for pure Lagrangian, except for at open boundaries,
            // where it is M_VT
            VC_x[i] = M_VT[x_ind]-(M_UM[x_ind]-UM_P[x_ind])/time_step;
            VC_y[i] = M_VT[y_ind]-(M_UM[y_ind]-UM_P[y_ind])/time_step;
        }

        for(int i=0;i<3;i++)
        {
            outer_fluxes_area[i] = 0;

            vertex_1 = other_vertex[2*i  ];
            vertex_2 = other_vertex[2*i+1];

            vector_edge[0] = x[vertex_2]-x[vertex_1];
            vector_edge[1] = y[vertex_2]-y[vertex_1];

            outer_vector[0] =  vector_edge[1];
            outer_vector[1] = -vector_edge[0];

            VC_middle[0] = (VC_x[vertex_2]+VC_x[vertex_1])/2.;
            VC_middle[1] = (VC_y[vertex_2]+VC_y[vertex_1])/2.;

            outer_fluxes_area[i] = outer_vector[0]*VC_middle[0]+outer_vector[1]*VC_middle[1];

            if(outer_fluxes_area[i]>=0)
            {
                surface = this->measure(M_elements[cpt],M_mesh, UM_P);
                outer_fluxes_area[i] = std::min(surface/time_step/3.,outer_fluxes_area[i]);
                fluxes_source_id[i]  = cpt;
            }
            else
            {
		        neighbour_double = bamgmesh->ElementConnectivity[cpt*3+i];
                neighbour_int = (int)bamgmesh->ElementConnectivity[cpt*3+i];

                // neighbour_double = M_element_connectivity[cpt*3+i];
                // neighbour_int = (int)M_element_connectivity[cpt*3+i];

		        if (!std::isnan(neighbour_double) && neighbour_int>0)
                {
                    surface = this->measure(M_elements[neighbour_int-1],M_mesh, UM_P);
                    outer_fluxes_area[i] = -std::min(surface/time_step/3.,-outer_fluxes_area[i]);
                    fluxes_source_id[i]  = neighbour_int-1;
                }
                else // open boundary with incoming fluxes
                    fluxes_source_id[i] = cpt;
            }
        }

        surface = this->measure(M_elements[cpt],M_mesh, UM_P);
        surface_new = this->measure(M_elements[cpt],M_mesh,M_UM);
        M_surface[cpt] = surface_new;

        for(int j=0; j<M_nb_var_element; j++)
        {
            if(M_interp_method[j]==1)
            {
                integrated_variable = interp_elt_in[cpt*M_nb_var_element+j]*surface
                    - (
                       interp_elt_in[fluxes_source_id[0]*M_nb_var_element+j]*outer_fluxes_area[0]
                       + interp_elt_in[fluxes_source_id[1]*M_nb_var_element+j]*outer_fluxes_area[1]
                       + interp_elt_in[fluxes_source_id[2]*M_nb_var_element+j]*outer_fluxes_area[2]
                       )*time_step;

                interp_elt_out[cpt*M_nb_var_element+j] = integrated_variable/surface_new;
            }
            else
            {
                interp_elt_out[cpt*M_nb_var_element+j] = interp_elt_in[cpt*M_nb_var_element+j];
            }
        }
    }
}//advect
#endif


#if 0//advectRoot not used - looks weird too
//------------------------------------------------------------------------------------------------------
//! Other advection function: unused
void
FiniteElement::advectRoot(std::vector<double> const& interp_elt_in, std::vector<double>& interp_elt_out)
{
    M_comm.barrier();

    bool use_eulerian = false;
    bool use_lagrangian = false;
    bool use_ALE = false;
    int ALE_smoothing_step_nb = vm["numerics.ALE_smoothing_step_nb"].as<int>();
    if (vm["numerics.advection_scheme"].as<string>()=="Eulerian")
        // Diffusive eulerian case where M_UM is not changed and then =0.
        use_eulerian = true;
    else if (vm["numerics.advection_scheme"].as<string>()=="Lagrangian")
        // ALE_smoothing_step_nb=0 is the purely Lagrangian case where M_UM is updated with M_VT
        use_lagrangian = true;
    else if (vm["numerics.advection_scheme"].as<string>()=="ALE")
    {
        // ALE case where M_UM is updated with a smoothed version of M_VT
        use_ALE = true;
        if (ALE_smoothing_step_nb<=0)
            throw std::runtime_error("numerics.ALE_smoothing_step_nb option should be >0");
    }
    else
        throw std::runtime_error("numerics.advection_scheme option should be Eulerian, Lagrangian or ALE");

    interp_elt_out.resize(M_nb_var_element*M_num_elements);
    std::vector<double> interp_elt_out_root;

    std::vector<double> M_VT_root;
    std::vector<double> M_VT_smoothed;
    std::vector<double> M_VT_smoothed_root;
    std::vector<double> M_UM_root;

    this->gatherNodalField(M_VT, M_UM, M_VT_root, M_UM_root);

    std::vector<double> interp_elt_in_root;
    this->gatherElementField(interp_elt_in, interp_elt_in_root, M_nb_var_element);

    //std::vector<double> M_surface_root;

    if (M_rank == 0)
    {
        std::vector<double> UM_P_root = M_UM_root;

        M_VT_smoothed_root = M_VT_root;
        std::vector<double> M_VT_tmp = M_VT_smoothed_root;

        // get the global number of nodes
        int num_nodes = M_mesh_root.numNodes();
        int Nd = bamgmesh_root->NodalConnectivitySize[1];

        for (int k=0; k<ALE_smoothing_step_nb; ++k)
        {
            M_VT_tmp = M_VT_smoothed_root;

            for (int i=0; i<num_nodes; ++i)
            {
                int Nc;
                double UM_x, UM_y;

                if(M_mask_dirichlet_root[i]==false)
                {
                    Nc = bamgmesh_root->NodalConnectivity[Nd*(i+1)-1];

                    UM_x = 0.;
                    UM_y = 0.;

                    for (int j=0; j<Nc; ++j)
                    {
                        UM_x += M_VT_tmp[bamgmesh_root->NodalConnectivity[Nd*i+j]-1];
                        UM_y += M_VT_tmp[bamgmesh_root->NodalConnectivity[Nd*i+j]-1+num_nodes];
                    }

                    M_VT_smoothed_root[i          ] = UM_x/Nc;
                    M_VT_smoothed_root[i+num_nodes] = UM_y/Nc;
                }
            }
        }

        for (int nd=0; nd<M_UM_root.size(); ++nd)
        {
            M_UM_root[nd] += time_step*M_VT_smoothed_root[nd];
        }

        for (const int& nd : M_neumann_nodes_root)
        {
            M_UM_root[nd] = UM_P_root[nd];
        }


        int num_elements = M_mesh_root.numTriangles();
        auto elements_root = M_mesh_root.triangles();
        auto nodes_root = M_mesh_root.nodes();

        interp_elt_out_root.resize(M_nb_var_element*num_elements);

        // std::cout<<"interp_elt_in_root size = "<< interp_elt_in_root.size() <<"\n";
        // std::cout<<"interp_elt_out_root size= "<< interp_elt_out_root.size() <<"\n";
        // M_surface_root.resize(num_elements);
#if 1
        for (int cpt=0; cpt < num_elements; ++cpt)
        {
            /* some variables used for the advection*/
            double surface, surface_new;
            double integrated_variable;

            /* some variables used for the advection*/
            double x[3],y[3],x_new[3],y_new[3];
            int x_ind, y_ind, neighbour_int, vertex_1, vertex_2;
            double outer_fluxes_area[3], vector_edge[2], outer_vector[2], VC_middle[2], VC_x[3], VC_y[3];
            int fluxes_source_id[3];

            double neighbour_double;
            int other_vertex[3*2]={1,2 , 2,0 , 0,1};

            /*======================================================================
             * Update:
             * Ice and snow thickness, and concentration using a Lagrangian or an Eulerian scheme
             *======================================================================
             */

            for(int i=0;i<3;i++)
            {
                x_ind = (elements_root[cpt]).indices[i]-1;
                y_ind = (elements_root[cpt]).indices[i]-1+num_nodes;

                x[i] = nodes_root[(elements_root[cpt]).indices[i]-1].coords[0];
                y[i] = nodes_root[(elements_root[cpt]).indices[i]-1].coords[1];

                /* old and new positions of the mesh */
                x_new[i] = x[i]+M_UM_root[x_ind];
                y_new[i] = y[i]+M_UM_root[y_ind];
                x[i]     = x[i]+UM_P_root[x_ind];
                y[i]     = y[i]+UM_P_root[y_ind];

                // VC_x,y are 0 in Lagrangian case
                VC_x[i] = M_VT_root[x_ind]-(M_UM_root[x_ind]-UM_P_root[x_ind])/time_step;
                VC_y[i] = M_VT_root[y_ind]-(M_UM_root[y_ind]-UM_P_root[y_ind])/time_step;
            }

            surface = this->measure(elements_root[cpt],M_mesh_root, UM_P_root);
            surface_new = this->measure(elements_root[cpt],M_mesh_root, M_UM_root);


            for(int i=0;i<3;i++)
            {
                outer_fluxes_area[i] = 0;

                vertex_1 = other_vertex[2*i  ];
                vertex_2 = other_vertex[2*i+1];

                vector_edge[0] = x[vertex_2]-x[vertex_1];
                vector_edge[1] = y[vertex_2]-y[vertex_1];

                outer_vector[0] = vector_edge[1];
                outer_vector[1] = -vector_edge[0];

                VC_middle[0] = (VC_x[vertex_2]+VC_x[vertex_1])/2.;
                VC_middle[1] = (VC_y[vertex_2]+VC_y[vertex_1])/2.;

                outer_fluxes_area[i] = outer_vector[0]*VC_middle[0]+outer_vector[1]*VC_middle[1];

                if(outer_fluxes_area[i]>0)
                {
                    // surface = this->measure(elements_root[cpt],M_mesh_root, UM_P_root);
                    outer_fluxes_area[i] = std::min(surface/time_step/3.,outer_fluxes_area[i]);
                    fluxes_source_id[i] = cpt;
                }
                else
                {
                    neighbour_double = bamgmesh_root->ElementConnectivity[cpt*3+i];
                    neighbour_int = (int)bamgmesh_root->ElementConnectivity[cpt*3+i];

                    // neighbour_double = M_element_connectivity[cpt*3+i];
                    // neighbour_int = (int)M_element_connectivity[cpt*3+i];

                    if (!std::isnan(neighbour_double) && neighbour_int>0)
                    {
                        double surface_local = this->measure(elements_root[neighbour_int-1],M_mesh_root, UM_P_root);
                        outer_fluxes_area[i] = -std::min(surface_local/time_step/3.,-outer_fluxes_area[i]);
                        fluxes_source_id[i] = neighbour_int-1;
                    }
                    else // open boundary with incoming fluxes
                    {
                        fluxes_source_id[i] = cpt;
                    }
                }
            }

            // surface = this->measure(elements_root[cpt],M_mesh_root, UM_P_root);
            // surface_new = this->measure(elements_root[cpt],M_mesh_root, M_UM_root);


#if 1
            for(int j=0; j<M_nb_var_element; j++)
            {
                if(M_interp_method[j]==1)
                {
                    integrated_variable = interp_elt_in_root[cpt*M_nb_var_element+j]*surface
                        - (
                           interp_elt_in_root[fluxes_source_id[0]*M_nb_var_element+j]*outer_fluxes_area[0]
                           + interp_elt_in_root[fluxes_source_id[1]*M_nb_var_element+j]*outer_fluxes_area[1]
                           + interp_elt_in_root[fluxes_source_id[2]*M_nb_var_element+j]*outer_fluxes_area[2]
                           )*time_step;

                    interp_elt_out_root[cpt*M_nb_var_element+j] = integrated_variable/surface_new;
                }
                else
                {
                    interp_elt_out_root[cpt*M_nb_var_element+j] = interp_elt_in_root[cpt*M_nb_var_element+j];
                }
            }
#endif
        }
#endif
    }//advection on root

#if 1
    this->scatterElementField(interp_elt_out_root, interp_elt_out, M_nb_var_element);

    this->scatterNodalField(M_VT_smoothed_root,M_VT_smoothed);

    std::vector<double> UM_P = M_UM;

    for (int nd=0; nd<M_UM.size(); ++nd)
    {
        M_UM[nd] += time_step*M_VT_smoothed[nd];
    }

    for (const int& nd : M_neumann_nodes)
    {
        M_UM[nd] = UM_P[nd];
    }

    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        M_surface[cpt] = this->measure(M_elements[cpt],M_mesh,M_UM);
    }
#endif
}//advectRoot
#endif


//------------------------------------------------------------------------------------------------------
//! Adds diffusion in slab ocean (for sst and sss).
//! Called by the update() function.
void
FiniteElement::diffuse(std::vector<double>& variable_elt, double diffusivity_parameters, double dx)
{
    if(diffusivity_parameters<=0.)
    {
        LOG(DEBUG) <<"diffusivity parameter lower or equal to 0 \n";
        LOG(DEBUG) <<"nothing to do\n";
        return;
    }

    // new addition
    std::vector<double> variable_elt_root;
    this->gatherElementField(variable_elt, variable_elt_root);

    if (M_rank == 0)
    {
        double factor = diffusivity_parameters*time_step/std::pow(dx,2.);
        std::vector<double> old_variable_elt = variable_elt_root;

        // get the global number of nodes
        int num_elements = M_mesh_root.numTriangles();

        for (int cpt=0; cpt < num_elements; ++cpt)
        {
            /* some variables used for the advection*/
            double fluxes_source[3];
            int fluxes_source_id;

            int neighbour_int;
            double neighbour_double;

            for(int i=0;i<3;i++)
            {
                neighbour_double = bamgmesh_root->ElementConnectivity[cpt*3+i];
                neighbour_int    = (int)bamgmesh_root->ElementConnectivity[cpt*3+i];

                // neighbour_double = M_element_connectivity[cpt*3+i];
                // neighbour_int = (int)M_element_connectivity[cpt*3+i];

                if (!std::isnan(neighbour_double) && neighbour_int>0)
                {
                    fluxes_source_id=neighbour_int-1;
                    fluxes_source[i]=factor*(old_variable_elt[fluxes_source_id]-old_variable_elt[cpt]);
                }
                else // no diffusion crosses open nor closed boundaries
                {
                    fluxes_source[i]=0.;
                }
            }

            variable_elt_root[cpt] += fluxes_source[0] + fluxes_source[1] + fluxes_source[2];
        }
    }

    // scatter back verctor from root to all processes
    this->scatterElementField(variable_elt_root, variable_elt);
}//diffuse

    
//------------------------------------------------------------------------------------------------------
//! ?? Has to do with the parallel computing.
//! Called by distributedMeshProcessing(), initMesh and  functions.
void
FiniteElement::scatterElementConnectivity()
{
    auto transfer_map_local = M_mesh.transferMapElt();
    std::vector<int> sizes_elements = M_sizes_elements_with_ghost;

    int nb_var_element = 3;
    std::vector<double> connectivity_root;

    if (M_rank == 0)
    {

        // std::cout<<"bamgmesh_root->ElementConnectivitySize[0]= "<< bamgmesh_root->ElementConnectivitySize[0] <<"\n";
        // std::cout<<"bamgmesh_root->ElementConnectivitySize[1]= "<< bamgmesh_root->ElementConnectivitySize[1] <<"\n";

        connectivity_root.resize(3*nb_var_element*M_id_elements.size());

        for (int i=0; i<M_id_elements.size(); ++i)
        {
            int ri = M_id_elements[i]-1;

            for (int j=0; j<nb_var_element; ++j)
            {
                double neighbour_id_db = bamgmesh_root->ElementConnectivity[nb_var_element*ri+j];
                int neighbour_id_int = (int)neighbour_id_db;

                if (!std::isnan(neighbour_id_db) && neighbour_id_int>0)
                {
                    connectivity_root[nb_var_element*i+j] = neighbour_id_db;
                }
                else
                {
                    connectivity_root[nb_var_element*i+j] = -100.;
                }
            }
        }
    }

    M_element_connectivity.resize(nb_var_element*M_num_elements);

    if (M_rank == 0)
    {
        std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });
        boost::mpi::scatterv(M_comm, connectivity_root, sizes_elements, &M_element_connectivity[0], 0);
    }
    else
    {
        boost::mpi::scatterv(M_comm, &M_element_connectivity[0], nb_var_element*M_num_elements, 0);
    }

    M_comm.barrier();

    auto element_connectivity_gid = M_element_connectivity;

    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        for (int j=0; j<nb_var_element; ++j)
        {
            double neighbour_id_global_db = element_connectivity_gid[nb_var_element*cpt+j];
            int neighbour_id_global_int = (int)element_connectivity_gid[nb_var_element*cpt+j];

            if (neighbour_id_global_int>0)
            {
                if (transfer_map_local.left.find(neighbour_id_global_int) != transfer_map_local.left.end())
                {
                    int neighbour_id_local = transfer_map_local.left.find(neighbour_id_global_int)->second;
                    M_element_connectivity[nb_var_element*cpt+j] = neighbour_id_local;
                }
                else
                {
                    M_element_connectivity[nb_var_element*cpt+j] = -100.;
                }
            }
        }
    }
}//scatterElementConnectivity


//------------------------------------------------------------------------------------------------------
//! Gathers information about the fields for interpolation onto the mesh grid.
//! Called by interpFields() function.
void
FiniteElement::gatherFieldsElement(std::vector<double>& interp_in_elements)
{
    //M_comm.barrier();

    int nb_var_element = M_nb_var_element;

    timer["gather"].first.restart();

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------GATHER ELEMENT starts\n";

    std::vector<int> sizes_elements = M_sizes_elements;
    //std::cout<<"------------------------------------------------------------------------------------M_nb_var_element= "<< M_nb_var_element <<"\n";
    std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });

    std::vector<double> interp_elt_in_local;
    this->collectVariables(interp_elt_in_local, false);

    if (M_rank == 0)
    {
        interp_in_elements.resize(nb_var_element*M_mesh_previous_root.numTriangles());
        boost::mpi::gatherv(M_comm, interp_elt_in_local, &interp_in_elements[0], sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, interp_elt_in_local, 0);
    }

    if (M_rank == 0)
    {
        auto interp_in_elements_nrd = interp_in_elements;

        for (int i=0; i<M_mesh_previous_root.numTriangles(); ++i)
        {
            int ri = M_rmap_elements[i];

            for (int j=0; j<nb_var_element; ++j)
            {
                interp_in_elements[nb_var_element*i+j] = interp_in_elements_nrd[nb_var_element*ri+j];
            }
        }
    }

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------GATHER ELEMENT done in "<< timer["gather"].first.elapsed() <<"s\n";
}//gatherFieldsElement


//------------------------------------------------------------------------------------------------------
//! Gathers information about the fields for outputing.
//! Called by the exportResults() function.
void
FiniteElement::gatherFieldsElementIO(std::vector<double>& interp_in_elements, bool thin_ice)
{
    int nb_var_element = M_nb_var_element;

    if (!thin_ice)
    {
        nb_var_element -= 4;
    }

    if (vm["output.save_diagnostics"].as<bool>())
    {
        nb_var_element += 7;
    }

    timer["gather"].first.restart();

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------IO: GATHER ELEMENT starts\n";

    std::vector<int> sizes_elements = M_sizes_elements;
    //std::cout<<"------------------------------------------------------------------------------------M_nb_var_element= "<< M_nb_var_element <<"\n";
    std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });

    std::vector<double> interp_elt_in_local;
    this->collectVariablesIO(interp_elt_in_local, false, thin_ice);

    if (M_rank == 0)
    {
        interp_in_elements.resize(nb_var_element*M_mesh_root.numTriangles());
        boost::mpi::gatherv(M_comm, interp_elt_in_local, &interp_in_elements[0], sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, interp_elt_in_local, 0);
    }

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------IO: GATHER ELEMENT done in "<< timer["gather"].first.elapsed() <<"s\n";
}//gatherFieldsElementsIO


//------------------------------------------------------------------------------------------------------
//! Scatters (redistributes) P0 (elemental) field values to the subdomains (parallel computing).
//! Called by the interpFields() function.
void
FiniteElement::scatterFieldsElement(double* interp_elt_out)
{
    timer["scatter"].first.restart();

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------SCATTER ELEMENT starts\n";

    int nb_var_element = M_nb_var_element;

    std::vector<int> sizes_elements = M_sizes_elements_with_ghost;
    std::vector<double> in_elt_values;

    if (M_rank == 0)
    {
        in_elt_values.resize(nb_var_element*M_id_elements.size());

        for (int i=0; i<M_id_elements.size(); ++i)
        {
            int ri = M_id_elements[i]-1;

            for (int j=0; j<nb_var_element; ++j)
            {
                in_elt_values[nb_var_element*i+j] = interp_elt_out[nb_var_element*ri+j];
            }
        }
    }

    std::vector<double> out_elt_values(nb_var_element*M_num_elements);

    if (M_rank == 0)
    {
        std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });
        boost::mpi::scatterv(M_comm, in_elt_values, sizes_elements, &out_elt_values[0], 0);
    }
    else
    {
        boost::mpi::scatterv(M_comm, &out_elt_values[0], nb_var_element*M_num_elements, 0);
    }

    // LOG(DEBUG) <<"["<< M_rank <<"]: " <<"Min val= "<< *std::min_element(out_elt_values.begin(), out_elt_values.end()) <<"\n";
    // LOG(DEBUG) <<"["<< M_rank <<"]: " <<"Max val= "<< *std::max_element(out_elt_values.begin(), out_elt_values.end()) <<"\n";


    M_conc.assign(M_num_elements,0.);
    M_thick.assign(M_num_elements,0.);
    M_snow_thick.assign(M_num_elements,0.);
    M_sigma.assign(3*M_num_elements,0.);
    M_damage.assign(M_num_elements,0.);
    M_ridge_ratio.assign(M_num_elements,0.);
    M_random_number.resize(M_num_elements);

    M_sst.resize(M_num_elements);
    M_sss.resize(M_num_elements);

    for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
        it->assign(M_num_elements,0.);

    M_h_thin.assign(M_num_elements,0.);
    M_conc_thin.assign(M_num_elements,0.);
    M_hs_thin.assign(M_num_elements,0.);
    M_tsurf_thin.assign(M_num_elements,0.);

    // Diagnostics
    D_Qa.assign(M_num_elements,0.);
    D_Qsh.assign(M_num_elements,0.);
    D_Qlh.assign(M_num_elements,0.);
    D_Qlw.assign(M_num_elements,0.);
    D_Qsw.assign(M_num_elements,0.);
    D_Qo.assign(M_num_elements,0.);
    D_delS.assign(M_num_elements,0.);

    this->redistributeVariables(out_elt_values,true);

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------SCATTER ELEMENT done in "<< timer["scatter"].first.elapsed() <<"s\n";
}//scatterFieldsElement

    
//------------------------------------------------------------------------------------------------------
//! Gets the names of the variables that need to be gathered and scattered when reading or saving restarts.
//! Called by the readRestart() function.
std::vector<std::string>
FiniteElement::getRestartVariableNames()
{

    std::vector<std::string> names = {
        "M_conc",
        "M_thick",
        "M_snow_thick",
        "M_sigma",//3 components
        "M_damage",
        "M_ridge_ratio",
        "M_random_number",
        "M_sss",
        "M_sst"};
    
    for(int i=0; i<M_tice.size(); i++)
        names.push_back("M_tice_" + std::to_string(i));
    if( M_ice_cat_type == setup::IceCategoryType::THIN_ICE)
    {
        names.push_back("M_h_thin");
        names.push_back("M_conc_thin");
        names.push_back("M_hs_thin");
        names.push_back("M_tsurf_thin");
    }
    return names;
}//getRestartVariableNames


//------------------------------------------------------------------------------------------------------
//! Redistributes all variables into the individual variables after scaterring from root.
//! Called by the readRestart() function.
void
FiniteElement::scatterFieldsElementIO(std::vector<double> const& interp_elt_out,
        std::vector<std::vector<double>*> &data,
        std::vector<int> const& num_components)
{
    //! * interp_elt_out is a vector containing all the variables to be
    //!   redistributed (eg after scattering from root) into the
    //!   individual variables (eg M_conc, M_thick,...)
    //!   - rearranged using M_id_elements and passed to
    //!     boost::mpi::scatterv
    //! * data is a vector of pointers to the variables to be assigned
    //!   values from out_elt_values
    //!   - passed to redistributeVariablesIO
    //! * num_components is a vector with the number of components in
    //!   each variable (usually 1, but can be 3 eg for M_sigma)
    //!   - passed to redistributeVariablesIO
    timer["scatter"].first.restart();

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------SCATTER ELEMENT starts\n";

    int const nb_var_element = std::accumulate(
            num_components.begin(), num_components.end(), 0);

    std::vector<int> sizes_elements = M_sizes_elements_with_ghost;
    std::vector<double> in_elt_values;

    if (M_rank == 0)
    {
        in_elt_values.resize(nb_var_element*M_id_elements.size());

        for (int i=0; i<M_id_elements.size(); ++i)
        {
            int ri = M_id_elements[i]-1;
            for (int j=0; j<nb_var_element; ++j)
            {
                in_elt_values[nb_var_element*i+j]
                    = interp_elt_out[nb_var_element*ri+j];
            }
        }
    }

    std::vector<double> out_elt_values(nb_var_element*M_num_elements);
    if (M_rank == 0)
    {
        std::for_each(sizes_elements.begin(), sizes_elements.end(),
                [&](int& f){ f = nb_var_element*f; });
        boost::mpi::scatterv(M_comm, in_elt_values, sizes_elements,
                &out_elt_values[0], 0);
    }
    else
    {
        boost::mpi::scatterv(M_comm, &out_elt_values[0],
                nb_var_element*M_num_elements, 0);
    }

    this->redistributeVariablesIO(out_elt_values, data, num_components);

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------SCATTER ELEMENT done in "<< timer["scatter"].first.elapsed() <<"s\n";
}//scatterFieldsElementIO


//------------------------------------------------------------------------------------------------------
//! Interpolates fields onto the mesh grid, e.g., after remeshing.
//! Called by the regrid() function.
void
FiniteElement::interpFields(std::vector<int> const& rmap_nodes, std::vector<int> sizes_nodes)
{
    std::vector<double> interp_in_elements;
    std::vector<double> interp_in_nodes;

    timer["gather"].first.restart();
    this->gatherFieldsElement(interp_in_elements);
    this->gatherFieldsNode(interp_in_nodes, rmap_nodes, sizes_nodes);
    if (M_rank == 0)
        std::cout<<"-------------------GATHER done in "<< timer["gather"].first.elapsed() <<"s\n";

    double* interp_elt_out;
    double* interp_nd_out;

    if (M_rank == 0)
    {
        // Interpolate elements

        std::vector<double> surface_previous(M_mesh_previous_root.numTriangles());
        std::vector<double> surface_root(M_mesh_root.numTriangles());
        //M_surface_root.resize(M_mesh_root.numTriangles());

        int cpt = 0;
        for (auto it=M_mesh_previous_root.triangles().begin(), end=M_mesh_previous_root.triangles().end(); it!=end; ++it)
        {
            surface_previous[cpt] = this->measure(*it,M_mesh_previous_root);
            ++cpt;
        }

        cpt = 0;
        for (auto it=M_mesh_root.triangles().begin(), end=M_mesh_root.triangles().end(); it!=end; ++it)
        {
            surface_root[cpt] = this->measure(*it,M_mesh_root);
            ++cpt;
        }

        //! The interpolation with the cavities still needs to be tested on a long run.
        //! By default, we then use the non-conservative MeshToMesh interpolation

        timer["cavities"].first.restart();
        InterpFromMeshToMesh2dCavities(&interp_elt_out,&interp_in_elements[0],
                                       &M_interp_method[0], M_nb_var_element,
                                       &surface_previous[0], &surface_root[0], bamgmesh_previous, bamgmesh_root);

        if (M_rank == 0)
            std::cout<<"-------------------CAVITIES done in "<< timer["cavities"].first.elapsed() <<"s\n";

#if 0
        // chrono.restart();
        // std::cout<<"InterpFromMeshToMesh2dx starts\n";

        // bamg::Mesh* Th;
        // Th = new bamg::Mesh(bamggeom_previous, bamgmesh_previous, bamgopt_previous);

        InterpFromMeshToMesh2dx(&interp_elt_out,
                                //Th,&M_mesh_previous_root.coordX()[0],&M_mesh_previous_root.coordY()[0],
                                &M_mesh_previous_root.indexTr()[0],&M_mesh_previous_root.coordX()[0],&M_mesh_previous_root.coordY()[0],
                                M_mesh_previous_root.numNodes(),M_mesh_previous_root.numTriangles(),
                                &interp_in_elements[0],
                                M_mesh_previous_root.numTriangles(),M_nb_var_element,
                                &M_mesh_root.bCoordX()[0],&M_mesh_root.bCoordY()[0],M_mesh_root.numTriangles(),
                                false);

        // std::cout<<"InterpFromMeshToMesh2dx done in "<< chrono.elapsed() <<"\n";
#endif

        // Interpolate nodes
        InterpFromMeshToMesh2dx(&interp_nd_out,
                                &M_mesh_previous_root.indexTr()[0],&M_mesh_previous_root.coordX()[0],&M_mesh_previous_root.coordY()[0],
                                M_mesh_previous_root.numNodes(),M_mesh_previous_root.numTriangles(),
                                &interp_in_nodes[0],
                                M_mesh_previous_root.numNodes(),M_nb_var_node,
                                &M_mesh_root.coordX()[0],&M_mesh_root.coordY()[0],M_mesh_root.numNodes(),
                                false);
    } // rank 0

    timer["distributed"].first.restart();
    this->distributedMeshProcessing();
    if (M_rank == 0)
        std::cout<<"-------------------DISTRIBUTED done in "<< timer["distributed"].first.elapsed() <<"s\n";

    timer["scatter"].first.restart();
    this->scatterFieldsElement(interp_elt_out);
    this->scatterFieldsNode(interp_nd_out);
    if (M_rank == 0)
        std::cout<<"-------------------SCATTER done in "<< timer["scatter"].first.elapsed() <<"s\n";

    if (M_rank == 0)
    {
        xDelete<double>(interp_elt_out);
        xDelete<double>(interp_nd_out);
    }
}//interpFields

    
//------------------------------------------------------------------------------------------------------
//! Gathers field values (velocities, displacements) at the nodes into a single structure, interp_node_in_local.
//! Called by the interpFields() function.
void
FiniteElement::gatherFieldsNode(std::vector<double>& interp_in_nodes, std::vector<int> const& rmap_nodes, std::vector<int> sizes_nodes)
{
    timer["gather.node"].first.restart();

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------GATHER NODE starts\n";

    M_nb_var_node = 10;
    std::vector<double> interp_node_in_local(M_nb_var_node*M_prv_local_ndof,0.);

    chrono.restart();
    //std::cout<<"Nodal Interp starts\n";
    //std::cout<<"NODAL: Interp starts\n";

    for (int i=0; i<M_prv_local_ndof; ++i)
    {
        // VT
        interp_node_in_local[M_nb_var_node*i] = M_VT[i];
        interp_node_in_local[M_nb_var_node*i+1] = M_VT[i+M_prv_num_nodes];

        // VTM
        interp_node_in_local[M_nb_var_node*i+2] = M_VTM[i];
        interp_node_in_local[M_nb_var_node*i+3] = M_VTM[i+M_prv_num_nodes];

        // VTMM
        interp_node_in_local[M_nb_var_node*i+4] = M_VTMM[i];
        interp_node_in_local[M_nb_var_node*i+5] = M_VTMM[i+M_prv_num_nodes];

        // UM
        interp_node_in_local[M_nb_var_node*i+6] = M_UM[i];
        interp_node_in_local[M_nb_var_node*i+7] = M_UM[i+M_prv_num_nodes];

        // UT
        interp_node_in_local[M_nb_var_node*i+8] = M_UT[i];
        interp_node_in_local[M_nb_var_node*i+9] = M_UT[i+M_prv_num_nodes];
    }

    std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = M_nb_var_node*f; });

    if (M_rank == 0)
    {
        interp_in_nodes.resize(M_nb_var_node*M_prv_global_num_nodes);
        boost::mpi::gatherv(M_comm, interp_node_in_local, &interp_in_nodes[0], sizes_nodes, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, interp_node_in_local, 0);
    }

    if (M_rank == 0)
    {
        auto interp_in_nodes_nrd = interp_in_nodes;

        for (int i=0; i<M_prv_global_num_nodes; ++i)
        {
            int ri =  rmap_nodes[i];

            for (int j=0; j<M_nb_var_node; ++j)
            {
                interp_in_nodes[M_nb_var_node*i+j] = interp_in_nodes_nrd[M_nb_var_node*ri+j];
            }
        }
    }

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------GATHER NODE done in "<< timer["gather.node"].first.elapsed() <<"s\n";
}//gatherFieldsNode


//------------------------------------------------------------------------------------------------------
//! Scatters field values (velocities, displacements) at the field nodes from the root
//! Called by the interpFields() and readRestart() functions.
void
FiniteElement::scatterFieldsNode(double* interp_nd_out)
{
    timer["scatter.node"].first.restart();

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------SCATTER NODE starts\n";

    std::vector<double> in_nd_values;

    if (M_rank == 0)
    {
        in_nd_values.resize(M_nb_var_node*M_id_nodes.size());

        for (int i=0; i<M_id_nodes.size(); ++i)
        {
            //int ri = rmap_nodes.right.find(id_nodes[i])->second-1;
            int ri = M_id_nodes[i]-1;

            for (int j=0; j<M_nb_var_node; ++j)
            {
                in_nd_values[M_nb_var_node*i+j] = interp_nd_out[M_nb_var_node*ri+j];
            }
        }
    }

    std::vector<double> out_nd_values(M_nb_var_node*M_num_nodes);
    std::vector<int> sizes_nodes = M_sizes_nodes_with_ghost;

    if (M_rank == 0)
    {
        std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = M_nb_var_node*f; });
        boost::mpi::scatterv(M_comm, in_nd_values, sizes_nodes, &out_nd_values[0], 0);
    }
    else
    {
        boost::mpi::scatterv(M_comm, &out_nd_values[0], M_nb_var_node*M_num_nodes, 0);
    }


    M_VT.assign(2*M_num_nodes,0.);
    M_VTM.assign(2*M_num_nodes,0.);
    M_VTMM.assign(2*M_num_nodes,0.);
    M_UM.assign(2*M_num_nodes,0.);
    M_UT.assign(2*M_num_nodes,0.);


    for (int i=0; i<M_num_nodes; ++i)
    {
        // VT
        M_VT[i] = out_nd_values[M_nb_var_node*i];
        M_VT[i+M_num_nodes] = out_nd_values[M_nb_var_node*i+1];

        // VTM
        M_VTM[i] = out_nd_values[M_nb_var_node*i+2];
        M_VTM[i+M_num_nodes] = out_nd_values[M_nb_var_node*i+3];

        // VTMM
        M_VTMM[i] = out_nd_values[M_nb_var_node*i+4];
        M_VTMM[i+M_num_nodes] = out_nd_values[M_nb_var_node*i+5];

        // UM
        M_UM[i] = out_nd_values[M_nb_var_node*i+6];
        M_UM[i+M_num_nodes] = out_nd_values[M_nb_var_node*i+7];

        // UT
        M_UT[i] = out_nd_values[M_nb_var_node*i+8];
        M_UT[i+M_num_nodes] = out_nd_values[M_nb_var_node*i+9];
    }


    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------SCATTER NODE done in "<< timer["scatter.node"].first.elapsed() <<"s\n";
}//scatterFieldsNode

    
//------------------------------------------------------------------------------------------------------
//! Sends displacement vector to the root process.
//! Called by the regrid() function.
void
FiniteElement::gatherNodalField(std::vector<double> const& field_local, std::vector<double>& field_root)
{
    std::vector<double> um_local(2*M_local_ndof,0.);
    for (int i=0; i<M_local_ndof; ++i)
    {
        um_local[2*i] = field_local[i]; //M_UM[i];
        um_local[2*i+1] = field_local[i+M_num_nodes]; //M_UM[i+M_num_nodes];
    }

    std::vector<int> sizes_nodes = M_sizes_nodes;
    std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = 2*f; });

    // send displacement vector to the root process (rank 0)

    if (M_rank == 0)
    {
        field_root.resize(2*M_ndof);
        boost::mpi::gatherv(M_comm, um_local, &field_root[0], sizes_nodes, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, um_local, 0);
    }

    if (M_rank == 0)
    {
        int global_num_nodes = M_mesh.numGlobalNodes();

        auto field_root_nrd = field_root;

        for (int i=0; i<global_num_nodes; ++i)
        {
            int ri =  M_rmap_nodes[i];

            field_root[i] = field_root_nrd[2*ri];
            field_root[i+global_num_nodes] = field_root_nrd[2*ri+1];
        }
    }
}//gatherNodalField

//------------------------------------------------------------------------------------------------------
//! Gathers nodal fields.
//! Only called by the advecRoot() function, which is not used anymore.
#if 0
void
FiniteElement::gatherNodalField(std::vector<double> const& field1_local, std::vector<double> const& field2_local,
                                std::vector<double>& field1_root, std::vector<double>& field2_root)
{
    std::vector<double> um_local(4*M_local_ndof,0.);
    for (int i=0; i<M_local_ndof; ++i)
    {
        // field1
        um_local[4*i] = field1_local[i];
        um_local[4*i+1] = field1_local[i+M_num_nodes];

        // field2
        um_local[4*i+2] = field2_local[i];
        um_local[4*i+3] = field2_local[i+M_num_nodes];
    }

    std::vector<int> sizes_nodes = M_sizes_nodes;
    std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = 4*f; });

    // send displacement vector to the root process (rank 0)

    std::vector<double> um_root;

    if (M_rank == 0)
    {
        um_root.resize(4*M_ndof);
        boost::mpi::gatherv(M_comm, um_local, &um_root[0], sizes_nodes, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, um_local, 0);
    }

    if (M_rank == 0)
    {
        int global_num_nodes = M_mesh.numGlobalNodes();

        field1_root.resize(2*M_ndof);
        field2_root.resize(2*M_ndof);

        auto um_root_nrd = um_root;

        for (int i=0; i<global_num_nodes; ++i)
        {
            int ri =  M_rmap_nodes[i];

            field1_root[i] = um_root_nrd[4*ri];
            field1_root[i+global_num_nodes] = um_root_nrd[4*ri+1];

            field2_root[i] = um_root_nrd[4*ri+2];
            field2_root[i+global_num_nodes] = um_root_nrd[4*ri+3];
        }
    }
}//gatherNodalField
#endif


//------------------------------------------------------------------------------------------------------
//! Scatter nodal fields.
//! Called by the advect() function.
void
FiniteElement::scatterNodalField(std::vector<double> const& field_root, std::vector<double>& field_local)
{
    std::vector<double> in_nd_values;

    if (M_rank == 0)
    {
        int global_num_nodes = M_mesh.numGlobalNodes();

        in_nd_values.resize(2*M_id_nodes.size());

        for (int i=0; i<M_id_nodes.size(); ++i)
        {
            int ri = M_id_nodes[i]-1;

            in_nd_values[2*i]   = field_root[ri];
            in_nd_values[2*i+1] = field_root[ri+global_num_nodes];

            // for (int j=0; j<2; ++j)
            // {
            //     in_nd_values[2*i+j] = field_root[2*ri+j];
            // }
        }
    }

    field_local.resize(2*M_num_nodes);
    std::vector<int> sizes_nodes = M_sizes_nodes_with_ghost;

    if (M_rank == 0)
    {
        std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = 2*f; });
        boost::mpi::scatterv(M_comm, in_nd_values, sizes_nodes, &field_local[0], 0);
    }
    else
    {
        boost::mpi::scatterv(M_comm, &field_local[0], 2*M_num_nodes, 0);
    }

    std::vector<double> field_local_copy = field_local;

    for (int i=0; i<M_num_nodes; ++i)
    {
        // U component
        field_local[i] = field_local_copy[2*i];
        // V component
        field_local[i+M_num_nodes] = field_local_copy[2*i+1];
    }
}//scatterNodalField


//------------------------------------------------------------------------------------------------------
//! Gather field values over elements.
//! Called by the advect(), diffuse() and updateDrifterPosition() functions.
void
FiniteElement::gatherElementField(std::vector<double> const& field_local, std::vector<double>& field_root, int nb_fields)
{
    std::vector<double> field_local_copy(nb_fields*M_local_nelements);

    for (int i=0; i<M_local_nelements; ++i)
    {
        // copy values without ghosts
        for (int j=0; j<nb_fields; ++j)
        {
            field_local_copy[nb_fields*i+j] = field_local[nb_fields*i+j];
        }
    }

    std::vector<int> sizes_elements = M_sizes_elements;

    if (nb_fields != 1)
    {
        std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_fields*f; });
    }

    if (M_rank == 0)
    {
        field_root.resize(nb_fields*M_mesh_root.numTriangles());
        boost::mpi::gatherv(M_comm, field_local_copy, &field_root[0], sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, field_local_copy, 0);
    }

    if (M_rank == 0)
    {
        auto field_root_nrd = field_root;

        for (int i=0; i<M_mesh_root.numTriangles(); ++i)
        {
            int ri = M_rmap_elements[i];

            for (int j=0; j<nb_fields; ++j)
            {
                field_root[nb_fields*i+j] = field_root_nrd[nb_fields*ri+j];
            }
        }
    }

}//gatherElementField


//------------------------------------------------------------------------------------------------------
//! Scatters back vector of field values at the elements from root to all processes.
//! Called by the advectRoot() function.
void
FiniteElement::scatterElementField(std::vector<double> const& field_root, std::vector<double>& field_local, int nb_fields)
{
    std::vector<double> field_root_extended;

    if (M_rank == 0)
    {
        field_root_extended.resize(nb_fields*M_id_elements.size());

        for (int i=0; i<M_id_elements.size(); ++i)
        {
            int ri = M_id_elements[i]-1;

            for (int j=0; j<nb_fields; ++j)
            {
                field_root_extended[nb_fields*i+j] = field_root[nb_fields*ri+j];
            }
        }
    }

    field_local.resize(nb_fields*M_num_elements);

    std::vector<int> sizes_elements = M_sizes_elements_with_ghost;
    if (nb_fields != 1)
    {
        std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_fields*f; });
    }

    if (M_rank == 0)
    {
        boost::mpi::scatterv(M_comm, field_root_extended, sizes_elements, &field_local[0], 0);
    }
    else
    {
        boost::mpi::scatterv(M_comm, &field_local[0], nb_fields*M_num_elements, 0);
    }
}//scatterElementField

    
//------------------------------------------------------------------------------------------------------
//! Sends the displacement vector to the root process
//! !Does not seem to be used!
#if 0
void
FiniteElement::gatherUM(std::vector<double>& um)
{
    std::vector<double> um_local(2*M_local_ndof,0.);
    for (int i=0; i<M_local_ndof; ++i)
    {
        um_local[2*i] = M_UM[i];
        um_local[2*i+1] = M_UM[i+M_num_nodes];
    }

    std::vector<int> sizes_nodes = M_sizes_nodes;
    std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = 2*f; });

    // send displacement vector to the root process (rank 0)

    chrono.restart();
    if (M_rank == 0)
    {
        um.resize(2*M_ndof);
        boost::mpi::gatherv(M_comm, um_local, &um[0], sizes_nodes, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, um_local, 0);
    }

    if (M_rank == 0)
    {
        //auto rmap_nodes = M_mesh.mapNodes();
        int prv_global_num_nodes = M_mesh.numGlobalNodes();

        auto interp_in_nodes_nrd = um;

        for (int i=0; i<prv_global_num_nodes; ++i)
        {
            //int ri =  rmap_nodes.left.find(i+1)->second-1;
            int ri =  M_rmap_nodes[i];

            um[i] = interp_in_nodes_nrd[2*ri];
            um[i+prv_global_num_nodes] = interp_in_nodes_nrd[2*ri+1];
        }
    }
}//gatherUM
#endif


//------------------------------------------------------------------------------------------------------
//! Performs the re-gridding.
//! -Called by the step() function.
void
FiniteElement::regrid(bool step)
{
    M_comm.barrier();

    timer["regrid"].first.restart();

    double displacement_factor = 2.;
    int substep_nb=1;
    int step_order=-1;
    bool flip = true;
    int substep = 0;

    std::vector<double> um_root;

    this->gatherNodalField(M_UM,um_root);

    if (M_rank == 0)
    {
        chrono.restart();
        LOG(DEBUG) <<"Flip starts\n";

        while (flip /*|| (minang<(vm["numerics.regrid_angle"].as<double>())/10.)*/)
        {
            ++substep;
            displacement_factor /= 2.;
            step_order++;

            flip = this->flip(M_mesh_root,um_root,displacement_factor);

            if (substep > 1)
                LOG(DEBUG) <<"FLIP DETECTED "<< substep-1 <<"\n";
        }

        LOG(DEBUG) <<"displacement_factor= "<< displacement_factor <<"\n";

        // int step_order_max = step_order;
        // boost::mpi::reduce(M_comm, step_order, step_order_max, boost::mpi::maximum<int>(), 0);
        // step_order = step_order_max;

        LOG(DEBUG) <<"["<< M_rank <<"] " << "STEP ORDER= "<< step_order <<"\n";

        substep_nb = std::pow(2,step_order);

        if(substep_nb!=1)
        {
            std::cout << substep_nb << "substeps will be needed for the remeshing!" <<"\n";
            std::cout << "Warning: It is probably due to very high ice speed, check your fields!\n";
        }

        LOG(DEBUG) <<"Flip done in "<< chrono.elapsed() <<"s\n";

        substep_nb = 1;

        for (int substep_i = 0; substep_i < substep_nb; substep_i++ )
        {
            //LOG(DEBUG) <<"substep_nb= "<< substep_nb <<"\n";

            timer["movemesh"].first.restart();
            LOG(DEBUG) <<"Move starts\n";
            M_mesh_root.move(um_root,displacement_factor);
            LOG(DEBUG) <<"Move done in "<< timer["movemesh"].first.elapsed() <<"s\n";

            timer["movevertices"].first.restart();
            LOG(DEBUG) <<"Move bamgmesh->Vertices starts\n";
            auto RX = M_mesh_root.coordX();
            auto RY = M_mesh_root.coordY();

            for (int id=0; id<bamgmesh_root->VerticesSize[0]; ++id)
            {
                bamgmesh_root->Vertices[3*id] = RX[id];
                bamgmesh_root->Vertices[3*id+1] = RY[id] ;
            }

            LOG(DEBUG) <<"Move bamgmesh->Vertices done in "<< timer["movevertices"].first.elapsed() <<"s\n";

            if(M_mesh_type==setup::MeshType::FROM_SPLIT)
            {
                timer["interpvertices"].first.restart();
                LOG(DEBUG) <<"Interp vertices starts\n";
                this->interpVertices();
                std::cout <<"Interp vertices done in "<< timer["interpvertices"].first.elapsed() <<"\n";
            }

#if 0
            had_remeshed=true;
            if(step && (vm["numerics.regrid_output_flag"].as<bool>()))
            {

                std::string tmp_string1    = (boost::format( "before_adaptMesh_%1%_mesh_adapt_step_%2%_substep_%3%" )
                                                                                 % step
                                                                                 % mesh_adapt_step
                                              % substep ).str();

                this->exportResults(tmp_string1,true,true,false);
            }
#endif

            timer["adaptmesh"].first.restart();
            LOG(DEBUG) <<"---TRUE AdaptMesh starts\n";
            this->adaptMesh();
            std::cout <<"---TRUE AdaptMesh done in "<< timer["adaptmesh"].first.elapsed() <<"s\n";

#if 0
            if(step && (vm["numerics.regrid_output_flag"].as<bool>()))
            {
                std::string tmp_string2    = (boost::format( "after_adaptMesh_%1%_mesh_adapt_step_%2%_substep_%3%" )
                                              % step
                                              % mesh_adapt_step
                                              % substep ).str();

                this->exportResults(tmp_string2,true,false,false);
            }
#endif

            // save mesh (only root process)

            std::cout<<"------------------------------version       = "<< M_mesh_root.version() <<"\n";
            std::cout<<"------------------------------ordering      = "<< M_mesh_root.ordering() <<"\n";
            std::cout<<"------------------------------format        = "<< M_mesh_fileformat <<"\n";
            std::cout<<"------------------------------space         = "<< vm["mesh.partitioner-space"].as<std::string>() <<"\n";
            std::cout<<"------------------------------partitioner   = "<< vm["mesh.partitioner"].as<std::string>() <<"\n";

            // Environment::logMemoryUsage("before partitioning...");
            timer["savemesh"].first.restart();
            LOG(DEBUG) <<"Saving mesh starts\n";
            if (M_partition_space == mesh::PartitionSpace::MEMORY)
                M_mesh_root.writeToGModel();
            else if (M_partition_space == mesh::PartitionSpace::DISK)
                M_mesh_root.writeToFile(M_partitioned_mesh_filename);
            std::cout <<"Saving mesh done in "<< timer["savemesh"].first.elapsed() <<"s\n";

            // partition the mesh on root process (rank 0)
            timer["meshpartition"].first.restart();
            LOG(DEBUG) <<"Partitioning mesh starts\n";
            M_mesh_root.partition(M_partitioned_mesh_filename,
                    M_partitioner, M_partition_space, M_mesh_fileformat);
            std::cout <<"Partitioning mesh done in "<< timer["meshpartition"].first.elapsed() <<"s\n";

            // Environment::logMemoryUsage("after partitioning...");
        }
    } // rank 0

    // --------------------------------BEGINNING-------------------------

    M_prv_local_ndof = M_local_ndof;
    M_prv_num_nodes = M_num_nodes;
    M_prv_num_elements = M_local_nelements;
    //bimap_type prv_rmap_nodes = M_mesh.mapNodes();
    std::vector<int> prv_rmap_nodes = M_mesh.mapNodes();
    M_prv_global_num_nodes = M_mesh.numGlobalNodes();
    M_prv_global_num_elements = M_mesh.numGlobalElements();
    std::vector<int> sizes_nodes = M_sizes_nodes;

    timer["interpFields"].first.restart();
    this->interpFields(prv_rmap_nodes, sizes_nodes);
    if (M_rank == 0)
        std::cout <<"interpFields done in "<< timer["interpFields"].first.elapsed() <<"s\n";

    // --------------------------------END-------------------------------

    if (M_rank == 0)
        std::cout <<"TIMER REGRIDDING= "<< timer["regrid"].first.elapsed() <<"s\n";

    this->assignVariables();
}//regrid


//------------------------------------------------------------------------------------------------------
//! Adapts the mesh grid.
//! Called by the regrid() function.
void
FiniteElement::adaptMesh()
{
    delete bamgopt_previous;
    delete bamggeom_previous;
    delete bamgmesh_previous;

    bamgopt_previous = new BamgOpts();
    bamggeom_previous = new BamgGeom();
    bamgmesh_previous = new BamgMesh();

    *bamgmesh_previous = *bamgmesh_root;
    *bamggeom_previous = *bamggeom_root;
    *bamgopt_previous = *bamgopt;

    // set dirichlet flags
    for (int edg=0; edg<bamgmesh_previous->EdgesSize[0]; ++edg)
    {
        int fnd = bamgmesh_previous->Edges[3*edg]/*-1*/;

        if ((std::binary_search(M_dirichlet_flags_root.begin(),M_dirichlet_flags_root.end(),fnd)))
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

    timer["bamgmesh"].first.restart();
    Bamgx(bamgmesh_root,bamggeom_root,bamgmesh_previous,bamggeom_previous,bamgopt_previous);
    std::cout <<"---BAMGMESH done in "<< timer["bamgmesh"].first.elapsed() <<"s\n";

    //! Imports the mesh from bamg, updates the boundary flags and node ID's
    this->importBamg(bamgmesh_root);
    this->updateBoundaryFlags();
    if(bamgopt->KeepVertices)
        this->updateNodeIds();
}//adaptMesh


//------------------------------------------------------------------------------------------------------
//!  Updates the node ID's after regriding and mesh adaptation. Called by the adaptMesh() function.
void
FiniteElement::updateNodeIds()
{
    // Recompute the node ids
    // - for during adaptMesh() if bamgopt->KeepVertices is set to 1
    //M_mesh_previous_root is the mesh before importBamg was called
    std::vector<int> old_nodes_id = M_mesh_previous_root.id();

    //M_mesh_root is the mesh after importBamg was called
    std::vector<int> new_nodes_id = M_mesh_root.id();

    int Boundary_id  = 0;
    int nb_new_nodes = 0;

    // The new id's will have values higher than the previous ones
    int first_new_node = *std::max_element(old_nodes_id.begin(), old_nodes_id.end())+1;

    for (int vert=0; vert<bamgmesh_root->VerticesSize[0]; ++vert)
    {
        if(M_mask_root[vert])
        {
            Boundary_id++;
            new_nodes_id[vert] = Boundary_id;
        }
        else
        {
            if(bamgmesh_root->PreviousNumbering[vert]==0)
            {
                nb_new_nodes++;
                new_nodes_id[vert] = first_new_node+nb_new_nodes-1;
            }
            else
            {
                new_nodes_id[vert] = old_nodes_id[bamgmesh_root->PreviousNumbering[vert]-1];
            }
        }
    }
    M_mesh_root.setId(new_nodes_id);
}//updateNodeIds

    
//------------------------------------------------------------------------------------------------------
//! Updates the boundary flags (Neumann vs Dirichlet) after regriding and mesh adaptation.
//! Called by the adaptMesh() function.
void
FiniteElement::updateBoundaryFlags()
{
    LOG(DEBUG) <<"CLOSED: FLAGS SIZE BEFORE= "<< M_dirichlet_flags_root.size() <<"\n";
    LOG(DEBUG) <<"OPEN  : FLAGS SIZE BEFORE= "<< M_neumann_flags_root.size() <<"\n";

    //! 1) Updates Dirichlet nodes
    M_dirichlet_flags_root.resize(0);
    M_neumann_flags_root.resize(0);

    //! 2) Gets the global number of nodes
    int num_nodes = M_mesh_root.numNodes();

    //! 3) Masks out the boundary nodes
    M_mask_root.assign(bamgmesh_root->VerticesSize[0],false) ;
    M_mask_dirichlet_root.assign(bamgmesh_root->VerticesSize[0],false) ;
    for (int vert=0; vert<bamgmesh_root->VerticesOnGeomVertexSize[0]; ++vert)
        M_mask_root[bamgmesh_root->VerticesOnGeomVertex[2*vert]-1] = true; // The factor 2 is because VerticesOnGeomVertex has 2 dimensions in bamg

    //! 4) Updates Dirichlet and Neumann flags
    std::vector<int> boundary_flags_root;
    for (int edg=0; edg<bamgmesh_root->EdgesSize[0]; ++edg)
    {
        boundary_flags_root.push_back(bamgmesh_root->Edges[3*edg]);
        if (bamgmesh_root->Edges[3*edg+2] == M_flag_fix)
            M_dirichlet_flags_root.push_back(bamgmesh_root->Edges[3*edg]);
    }

    std::sort(M_dirichlet_flags_root.begin(), M_dirichlet_flags_root.end());
    std::sort(boundary_flags_root.begin(), boundary_flags_root.end());
    std::set_difference(boundary_flags_root.begin(), boundary_flags_root.end(),
                        M_dirichlet_flags_root.begin(), M_dirichlet_flags_root.end(),
                        std::back_inserter(M_neumann_flags_root));


    //! 5) Updates dirichlet nodes
    M_dirichlet_nodes_root.resize(2*(M_dirichlet_flags_root.size()));
    for (int i=0; i<M_dirichlet_flags_root.size(); ++i)
    {
        M_dirichlet_nodes_root[2*i] = M_dirichlet_flags_root[i];
        M_dirichlet_nodes_root[2*i+1] = M_dirichlet_flags_root[i]+num_nodes;
        M_mask_dirichlet_root[M_dirichlet_flags_root[i]] = true;
    }

    //! 6) Updates neumann nodes
    M_neumann_nodes_root.resize(2*(M_neumann_flags_root.size()));
    for (int i=0; i<M_neumann_flags_root.size(); ++i)
    {
        M_neumann_nodes_root[2*i] = M_neumann_flags_root[i];
        M_neumann_nodes_root[2*i+1] = M_neumann_flags_root[i]+num_nodes;
    }

    LOG(DEBUG) <<"CLOSED: FLAGS SIZE AFTER= "<< M_dirichlet_flags_root.size() <<"\n";
    LOG(DEBUG) <<"OPEN  : FLAGS SIZE AFTER= "<< M_neumann_flags_root.size() <<"\n";
}//updateBoundaryFlags
    

//------------------------------------------------------------------------------------------------------
//! Assembles matrices for solver: fvdata and data, with each term of the momentum equation.
//! Called by the step() function.
void
FiniteElement::assemble(int pcpt)
{
    M_comm.barrier();
    LOG(DEBUG) << "Reinitialize matrix and vector to zero starts\n";
    M_matrix->zero();
    M_vector->zero();
    LOG(DEBUG) << "Reinitialize matrix and vector to zero done\n";


    //std::vector<int> extended_dirichlet_nodes = M_dirichlet_nodes;

    // ---------- Identical values for all the elements -----------
    // coriolis term
    double beta0;
    double beta1;
    double beta2;
    if (pcpt > 1)
    {
        // Adams-Bashfort 3 (AB3)
        beta0 = 23./12;
        beta1 =-16./12;
        beta2 =  5./12;
    }
    else if (pcpt == 1)
    {
        // Adams-Bashfort 2 (AB2)
        beta0 = 3/2;
        beta1 =-1/2;
        beta2 = 0  ;
    }
    else if (pcpt == 0)
    {
        // Euler explicit (Fe)
        beta0 = 1 ;
        beta1 = 0 ;
        beta2 = 0 ;
    }

    double cos_ocean_turning_angle = std::cos(ocean_turning_angle_rad);
    double sin_ocean_turning_angle = std::sin(ocean_turning_angle_rad);

    // ---------- Assembling starts -----------
    LOG(DEBUG) <<"Assembling starts\n";
    timer["assembly"].first.restart();

    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {

        // total thickness and concentration
        double total_concentration=M_conc[cpt];
        double total_thickness=M_thick[cpt];
        double total_snow=M_snow_thick[cpt];

        // Add the thin ice concentration and thickness
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            total_concentration += M_conc_thin[cpt];
            total_thickness     += M_h_thin[cpt];
            total_snow          += M_hs_thin[cpt];
        }

        // Limits to avoid very small values
        total_thickness =       (vm["dynamics.min_h"].as<double>()>total_thickness)        ? vm["dynamics.min_h"].as<double>() : total_thickness;
        total_concentration =   (vm["dynamics.min_c"].as<double>()>total_concentration)    ? vm["dynamics.min_c"].as<double>() : total_concentration;

        double coef_min = 100.;

        // values used when no ice or when ice too thin
        double coef_drag    = 0.;  // coef_drag is a switch that set the external forcings to 0 (wind, ocean, bottom drag, waves stress) where there is too little ice
        double coef         = coef_min;
        double mass_e       = 0.;
        double coef_C       = 0.;
        double coef_V       = 0.;
        double coef_X       = 0.;
        double coef_Y       = 0.;
        double coef_sigma   = 0.;


        // temporary variables
        double mloc = 0.;
        double dloc = 0.;

        double b0tj_sigma_hu = 0.;
        double b0tj_sigma_hv = 0.;

        double coef_Vair = 0.;
        double coef_Voce = 0.;
        double coef_basal = 0.;

        double undamaged_time_relaxation_sigma = vm["dynamics.undamaged_time_relaxation_sigma"].as<double>();
        double exponent_relaxation_sigma = vm["dynamics.exponent_relaxation_sigma"].as<double>();

        double time_viscous = undamaged_time_relaxation_sigma*std::pow(1.-M_damage[cpt],exponent_relaxation_sigma-1.);
        double multiplicator = time_viscous/(time_viscous+time_step);

        double norm_Voce_ice = 0.;
        double norm_Voce_ice_min = 0.01; // minimum value to avoid 0 water drag term.

        double norm_Vair_ice = 0.;
        double norm_Vair_ice_min = 0.01; // minimum value to avoid 0 water drag term.

        double norm_Vice = 0.;

        double element_ssh = 0.;
        double critical_h = 0.;
        double max_keel_height=28; // [m] from "A comprehensive analysis of the morphology of first-year sea ice ridges"
        double ice_to_keel_factor=19.28; // from "A comprehensive analysis of the morphology of first-year sea ice ridges"
        double keel_height_estimate;
        double critical_h_mod = 0.;

        if( (M_conc[cpt] > vm["dynamics.min_c"].as<double>()) && (M_thick[cpt] > vm["dynamics.min_h"].as<double>()) )
        {
            /* Compute the value that only depends on the element */
            double welt_ice = 0.;
            double welt_ssh = 0.;
            int nind;

            for (int i=0; i<3; ++i)
            {
                nind = (M_elements[cpt]).indices[i]-1;
                welt_ssh += M_ssh[nind];
            }

            element_ssh = welt_ssh/3.;

            if(M_conc[cpt]>vm["dynamics.min_c"].as<double>())
            {
                switch ( M_basal_stress_type )
                {
                case setup::BasalStressType::NONE:
                    // No grounding
                    critical_h     = 0.;
                    critical_h_mod = 0.;
                case setup::BasalStressType::BOUILLON:
                    // Sylvain's grounding scheme
                    keel_height_estimate = ice_to_keel_factor*std::pow(M_thick[cpt]/M_conc[cpt],0.5);
                    keel_height_estimate = ( keel_height_estimate > max_keel_height ) ? max_keel_height : keel_height_estimate;
                    critical_h      = M_conc[cpt]*std::pow((M_element_depth[cpt]+element_ssh)/ice_to_keel_factor,2.);
                    critical_h_mod  = M_conc[cpt]*std::pow(keel_height_estimate/ice_to_keel_factor,2.);
                    break;
                case setup::BasalStressType::LEMIEUX:
                    // JF Lemieux's grounding
                    keel_height_estimate = vm["dynamics.Lemieux_basal_k1"].as<double>()*M_thick[cpt]/M_conc[cpt];
                    keel_height_estimate = ( keel_height_estimate > max_keel_height ) ? max_keel_height : keel_height_estimate;

                    critical_h      = M_conc[cpt]*(M_element_depth[cpt]+element_ssh)/(vm["dynamics.Lemieux_basal_k1"].as<double>());
                    critical_h_mod  = M_conc[cpt]*keel_height_estimate/(vm["dynamics.Lemieux_basal_k1"].as<double>());
                    break;
                }
            }

            if(young>0.) // EB rheology
            {
                coef = multiplicator*young*(1.-M_damage[cpt])*M_thick[cpt]*std::exp(ridging_exponent*(1.-M_conc[cpt]));
            }
            else // Linear viscous rheology where nominal viscosity is defined as -young*time_step
            {
                double norm_factor=vm["dynamics.cohesion_thickness_normalisation"].as<double>();
                double exponent=vm["dynamics.cohesion_thickness_exponent"].as<double>();
                double mult_factor = std::pow(M_thick[cpt]/norm_factor,exponent);
                coef = -young*M_thick[cpt]*mult_factor*std::exp(ridging_exponent*(1.-M_conc[cpt]));
            }

            coef = (coef<coef_min) ? coef_min : coef ;

            if (vm["dynamics.use_coriolis"].as<bool>())
                mass_e = (rhoi*total_thickness + rhos*total_snow)/total_concentration;
            else
                mass_e = 0.;

            /* compute the x and y derivative of g*ssh, for the sea surface tilt term */
            double g_ssh_e_x = 0.;
            double g_ssh_e_y = 0.;
            double g_ssh_e;

            for(int i=0; i<3; i++)
            {
                g_ssh_e = (physical::gravity)*M_ssh[(M_elements[cpt]).indices[i]-1] /*g_ssh*/;   /* g*ssh at the node k of the element e */
                g_ssh_e_x += M_shape_coeff[cpt][i]*g_ssh_e; /* x derivative of g*ssh */
                g_ssh_e_y += M_shape_coeff[cpt][i+3]*g_ssh_e; /* y derivative of g*ssh */
            }

            coef_drag  = 1.;
            coef_C     = mass_e*M_fcor[cpt];              /* for the Coriolis term */
            coef_V     = mass_e/time_step;                /* for the inertial term */
            coef_X     = - mass_e*g_ssh_e_x;              /* for the ocean slope */
            coef_Y     = - mass_e*g_ssh_e_y;              /* for the ocean slope */
            coef_sigma = M_thick[cpt]*multiplicator;      /* for the internal stress */
        }

        std::vector<int> rindices(6); //new
        std::vector<int> cindices(6);
        int vs = 0;

        for (int s=0; s<3; ++s)
        {
            int index_u = (M_elements[cpt]).indices[s]-1;

            if (!((M_elements[cpt]).ghostNodes[s]))
            {
                rindices[2*vs] = index_u;
                rindices[2*vs+1] = index_u+M_local_ndof;
                ++vs;
            }

            cindices[2*s] = index_u;
            cindices[2*s+1] = index_u+M_local_ndof_ghost;
        }

        rindices.resize(2*vs);

        /* Loop over the 6 by 6 components of the finite element integral
         * this is done smartly by looping over j=0:2 and i=0:2
         * col = (mwIndex)it[2*j]-1  , row = (mwIndex)it[2*i]-1;
         * col  , row   -> UU component
         * col  , row+1 -> VU component
         * col+1, row   -> VV component
         * col+1, row+1 -> UV component */

        double Vcor_index_v, Vcor_index_u;
        double surface_e = M_surface[cpt];
        double duu, dvu, duv, dvv;
        std::vector<double> data(36);
        std::vector<double> fvdata(6,0.);

        int l_j = -1; // node counter to skip gohsts
        for(int j=0; j<3; j++)
        {
            /* Column corresponding to indice j (we also assemble terms in col+1) */
            //col = (mwIndex)it[2*j]-1; /* -1 to use the indice convention of C */

            int index_u = (M_elements[cpt]).indices[j]-1;
            int index_v = (M_elements[cpt]).indices[j]-1+M_num_nodes;

            double vt_u = M_VT[index_u];
            double vt_v = M_VT[index_v];

            double ocean_u = M_ocean[index_u];
            double ocean_v = M_ocean[index_v];

            double wind_u = M_wind[index_u];
            double wind_v = M_wind[index_v];

            Vcor_index_v = beta0*vt_v + beta1*M_VTM[index_v] + beta2*M_VTMM[index_v];
            Vcor_index_u = beta0*vt_u + beta1*M_VTM[index_u] + beta2*M_VTMM[index_u];

            /* Skip ghost nodes */
            if (!((M_elements[cpt]).ghostNodes[j]))
            {
                l_j = l_j + 1;

                for(int i=0; i<3; i++)
                {
                    /* Row corresponding to indice i (we also assemble terms in row+1) */

                    /* Select the nodal weight values from M_loc */
                    mloc = M_Mass[3*j+i];
                    dloc = M_Diag[3*j+i];

                    b0tj_sigma_hu = 0.;
                    b0tj_sigma_hv = 0.;

                    for(int k=0; k<3; k++)
                    {
                        b0tj_sigma_hu += M_B0T[cpt][k*6+2*i]*(M_sigma[3*cpt+k]*coef_sigma/*+sigma_P[k]*/);
                        b0tj_sigma_hv += M_B0T[cpt][k*6+2*i+1]*(M_sigma[3*cpt+k]*coef_sigma/*+sigma_P[k]*/);
                    }

                    /* ---------- UU component */

                    norm_Voce_ice = std::hypot(M_VT[index_u]-M_ocean[index_u],M_VT[index_v]-M_ocean[index_v]);
                    norm_Voce_ice = (norm_Voce_ice > norm_Voce_ice_min) ? (norm_Voce_ice):norm_Voce_ice_min;

                    coef_Voce = (vm["dynamics.lin_drag_coef_water"].as<double>()+(quad_drag_coef_water*norm_Voce_ice));
                    coef_Voce *= coef_drag*physical::rhow;

                    norm_Vair_ice = std::hypot(M_VT[index_u]-M_wind [index_u],M_VT[index_v]-M_wind [index_v]);
                    norm_Vair_ice = (norm_Vair_ice > norm_Vair_ice_min) ? (norm_Vair_ice):norm_Vair_ice_min;

                    coef_Vair = (vm["dynamics.lin_drag_coef_air"].as<double>()+(quad_drag_coef_air*norm_Vair_ice));
                    coef_Vair *= coef_drag*(physical::rhoa);

                    norm_Vice = std::hypot(M_VT[index_u],M_VT[index_v]);
                    norm_Vice = (norm_Vice > basal_u_0) ? (norm_Vice):basal_u_0;

                    coef_basal = basal_k2/norm_Vice;
                    coef_basal *= coef_drag*std::max(0., critical_h_mod-critical_h)*std::exp(-basal_Cb*(1.-M_conc[cpt]));

                    duu = surface_e*( mloc*(coef_V)
                                      +dloc*(coef_Vair+coef_basal+coef_Voce*cos_ocean_turning_angle)
                                      +M_B0T_Dunit_B0T[cpt][(2*i)*6+2*j]*coef*time_step );

                    /* ---------- VU component */
                    dvu = surface_e*( +M_B0T_Dunit_B0T[cpt][(2*i+1)*6+2*j]*coef*time_step );

                    /* ---------- UV component */
                    duv = surface_e*( +M_B0T_Dunit_B0T[cpt][(2*i)*6+2*j+1]*coef*time_step );

                    /* ---------- VV component */
                    dvv = surface_e*( mloc*(coef_V)
                                      +dloc*(coef_Vair+coef_basal+coef_Voce*cos_ocean_turning_angle)
                                      +M_B0T_Dunit_B0T[cpt][(2*i+1)*6+2*j+1]*coef*time_step );

                    data[(2*l_j  )*6+2*i  ] = duu;
                    data[(2*l_j+1)*6+2*i  ] = dvu;
                    data[(2*l_j  )*6+2*i+1] = duv;
                    data[(2*l_j+1)*6+2*i+1] = dvv;

                    fvdata[2*i] += surface_e*( mloc*(
                                                     // TODO: +coef_drag*M_tau[index_u] <- waves are not here yet!!
                                                     +coef_X
                                                     +coef_V*vt_u
                                                     +coef_C*Vcor_index_v
                                                     )
                                               +dloc*(
                                                      +coef_Vair*wind_u
                                                      +coef_Voce*cos_ocean_turning_angle*ocean_u
                                                      -coef_Voce*sin_ocean_turning_angle*(ocean_v-vt_v)
                                                      )
                                               - b0tj_sigma_hu/3);

                    fvdata[2*i+1] += surface_e*( mloc*(
                                                       // TODO: +coef_drag*M_tau[index_v] <- waves are not here yet!!
                                                       +coef_Y
                                                       +coef_V*vt_v
                                                       -coef_C*Vcor_index_u
                                                       )
                                                 +dloc*(
                                                        +coef_Vair*wind_v
                                                        +coef_Voce*cos_ocean_turning_angle*ocean_v
                                                        +coef_Voce*sin_ocean_turning_angle*(ocean_u-vt_u)
                                                        )
                                                 - b0tj_sigma_hv/3);
                }
            }
        }

        // update matrix
        M_matrix->addMatrix(&rindices[0], rindices.size(), &cindices[0], cindices.size(), &data[0]);

        // update vector
        M_vector->addVector(&cindices[0], cindices.size(), &fvdata[0]);


    }

    // close petsc matrix
    LOG(DEBUG) <<"Closing matrix starts\n";
    M_matrix->close();
    LOG(DEBUG) <<"Closing matrix done\n";

    // close petsc vector
    LOG(DEBUG) <<"Closing vector starts\n";
    M_vector->close();
    LOG(DEBUG) <<"Closing vector done\n";

    M_matrix->on(M_dirichlet_nodes,*M_vector);

}//assemble


//------------------------------------------------------------------------------------------------------
//! Assembles the mass matrix.
//! Called by the step() function.
void
FiniteElement::FETensors()
{
    M_Dunit.assign(9,0);
    M_Mass.assign(9,0);
    M_Diag.assign(9,0);

    for (int k=0; k<6; k+=3)
    {
        for (int kk=0; kk<2; ++kk )
        {
            M_Dunit[k+kk] = (1-((k+kk)%2)*(1-nu0))/(1-std::pow(nu0,2.));
        }
    }
    M_Dunit[8] = (1-nu0)/(2.*(1-std::pow(nu0,2.)));

    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            M_Mass[3*i+j] = ((i == j) ? 2.0 : 1.0)/12.0;
            M_Diag[3*i+j] = ((i == j) ? 1.0 : 0.0)/3.0;
            //std::cout<< std::left << std::setw(12) << Mass[3*i+j] <<"  ";
        }
    }

    M_B0T.resize(M_num_elements);
    M_B0T_Dunit_B0T.resize(M_num_elements);
    M_shape_coeff.resize(M_num_elements);

    std::vector<double> B0T(18,0);
    std::vector<double> B0Tj_Dunit(6,0);
    std::vector<double> B0T_Dunit_B0T(36,0);

    double B0Tj_Dunit_tmp0, B0Tj_Dunit_tmp1;
    double B0Tj_Dunit_B0Ti_tmp0, B0Tj_Dunit_B0Ti_tmp1, B0Tj_Dunit_B0Ti_tmp2, B0Tj_Dunit_B0Ti_tmp3;

    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        std::vector<double> shapecoeff = this->shapeCoeff(*it,M_mesh);

        for (int i=0; i<18; ++i)
        {
            if (i < 3)
            {
                B0T[2*i] = shapecoeff[i];
                B0T[12+2*i] = shapecoeff[i+3];
                B0T[13+2*i] = shapecoeff[i];
            }
            else if (i < 6)
            {
                B0T[2*i+1] = shapecoeff[i];
            }
        }

        //std::cout<<"\n";
        // if (cpt == 0)
        //     for (int i=0; i<3; ++i)
        //     {
        //         for (int j=0; j<6; ++j)
        //         {
        //             std::cout<< std::left << std::setw(12) << B0T[6*i+j] <<"  ";
        //         }
        //         std::cout<<"\n";
        //     }

        for(int j=0; j<3; j++)
        {

            /* The rigidity matrix that will be multiplied by E and by the surface
             * is given by the product B0'*matrix.Dunit*B0
             * This product is computed for the indices:
             * 2*i  ,2*j   -> B0Tj_Dunit_B0Ti[0]
             * 2*i  ,2*j+1 -> B0Tj_Dunit_B0Ti[1]
             * 2*i+1,2*j   -> B0Tj_Dunit_B0Ti[2]
             * 2*i+1,2*j+1 -> B0Tj_Dunit_B0Ti[3] */

            /* new version without assembling zero */
            for(int i=0; i<3; i++)
            {
                /* product of the first line of B0T' and the matrix Dunit */
                B0Tj_Dunit_tmp0 = 0.;
                B0Tj_Dunit_tmp1 = 0.;

                for(int kk=0; kk<3; kk++)
                {
                    B0Tj_Dunit_tmp0 += B0T[kk*6+2*j]*M_Dunit[3*i+kk];
                    B0Tj_Dunit_tmp1 += B0T[kk*6+2*j+1]*M_Dunit[3*i+kk];
                }

                B0Tj_Dunit[2*i] = B0Tj_Dunit_tmp0;
                B0Tj_Dunit[2*i+1] = B0Tj_Dunit_tmp1;
            }

            for(int i=0; i<3; i++)
            {
                /* The rigidity matrix */
                /* scalar product of B0Ti_Dunit and the first column of B0T */
                B0Tj_Dunit_B0Ti_tmp0 = 0.;
                B0Tj_Dunit_B0Ti_tmp1 = 0.;
                B0Tj_Dunit_B0Ti_tmp2 = 0.;
                B0Tj_Dunit_B0Ti_tmp3 = 0.;

                for(int kk=0; kk<3; kk++)
                {
                    B0Tj_Dunit_B0Ti_tmp0 += B0Tj_Dunit[2*kk]*B0T[kk*6+2*i];
                    B0Tj_Dunit_B0Ti_tmp1 += B0Tj_Dunit[2*kk]*B0T[kk*6+2*i+1];
                    B0Tj_Dunit_B0Ti_tmp2 += B0Tj_Dunit[2*kk+1]*B0T[kk*6+2*i];
                    B0Tj_Dunit_B0Ti_tmp3 += B0Tj_Dunit[2*kk+1]*B0T[kk*6+2*i+1];
                }

                B0T_Dunit_B0T[(2*i)*6+2*j] = B0Tj_Dunit_B0Ti_tmp0;
                B0T_Dunit_B0T[(2*i+1)*6+2*j] = B0Tj_Dunit_B0Ti_tmp1;
                B0T_Dunit_B0T[(2*i)*6+2*j+1] = B0Tj_Dunit_B0Ti_tmp2;
                B0T_Dunit_B0T[(2*i+1)*6+2*j+1] = B0Tj_Dunit_B0Ti_tmp3;
            }
        }

        /*
         * B0T_Dunit_B0T should be symmetric but is not exactly after the calcultation here above
         * because the sequence of operation is not the same for the component i,j and j,i.
         * We force the matrix to be symmetric by copying the upper part onto the lower part
         */
        for(int i=1; i<6; i++)
            for(int j=0; j<i; j++)
                B0T_Dunit_B0T[i*6+j] = B0T_Dunit_B0T[j*6+i];

        M_shape_coeff[cpt]        = shapecoeff;
        M_B0T[cpt]                = B0T;
        M_B0T_Dunit_B0T[cpt]      = B0T_Dunit_B0T;

        ++cpt;
    }
}//FETensors

    
//------------------------------------------------------------------------------------------------------
//! Calculates the cohesion field (sum of a fixed value and a random component) and the maximum compressive strength of sea ice.
//! Called by the step() function.
void
FiniteElement::calcCohesion()
{
    for (int i=0; i<M_Cohesion.size(); ++i)
        M_Cohesion[i] = C_fix+C_alea*(M_random_number[i]-0.5);

    for (int i=0; i<M_Compressive_strength.size(); ++i)
        M_Compressive_strength[i] = compr_strength*scale_coef;
}//calcCohesion

    
//------------------------------------------------------------------------------------------------------
//! Update all relevant fields and physical variables after solving. Called by the step() function.
void
FiniteElement::update()
{
    // Hotfix for issue #53 - we only have pure Lagrangian now.
    // // collect the variables into a single structure
    // std::vector<double> interp_elt_in_local;
    // this->collectVariables(interp_elt_in_local, true);
    //
    // // advect
    // std::vector<double> interp_elt_out;
    // this->advect(interp_elt_in_local, interp_elt_out);
    //
    // // redistribute the interpolated values
    // this->redistributeVariables(interp_elt_out);

    std::vector<double> UM_P = M_UM;
    for (int nd=0; nd<M_UM.size(); ++nd)
        M_UM[nd] += time_step*M_VT[nd];

    for (const int& nd : M_neumann_nodes)
        M_UM[nd] = UM_P[nd];

    // Horizontal diffusion
    this->diffuse(M_sst,vm["thermo.diffusivity_sst"].as<double>(),M_res_root_mesh);
    this->diffuse(M_sss,vm["thermo.diffusivity_sss"].as<double>(),M_res_root_mesh);

    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        double old_damage;

        /* deformation, deformation rate and internal stress tensor and temporary variables */
        double epsilon_veloc_i;
        std::vector<double> epsilon_veloc(3);

        /* divergence */
        double divergence_rate;

        /* shear*/
        double shear_rate;

        // ridging scheme
        double delta_ridging;
        double G_star=0.15;
        double e_factor=2.;

        std::vector<double> sigma_pred(3);
        double sigma_dot_i;

        /* invariant of the internal stress tensor and some variables used for the damaging process*/
        double sigma_s, sigma_n, sigma_1, sigma_2;
        double tract_max, sigma_t, sigma_c, q;
        double tmp, sigma_target;

        // Temporary memory
        old_damage = M_damage[cpt];

        /*======================================================================
         * Diagnostic:
         * Elastic deformation and instantaneous deformation rate
         *======================================================================
         */

        //! - Computes the elastic deformation and the instantaneous deformation rate
        for(int i=0;i<3;i++)
        {
            epsilon_veloc_i = 0.0;
            for(int j=0;j<3;j++)
            {
                /* deformation */
                epsilon_veloc_i += M_B0T[cpt][i*6 + 2*j]*M_VT[(M_elements[cpt]).indices[j]-1];
                epsilon_veloc_i += M_B0T[cpt][i*6 + 2*j + 1]*M_VT[(M_elements[cpt]).indices[j]-1+M_num_nodes];
            }

            epsilon_veloc[i] = epsilon_veloc_i;
        }

        divergence_rate = (epsilon_veloc[0]+epsilon_veloc[1]);
        shear_rate= std::hypot(epsilon_veloc[0]-epsilon_veloc[1],epsilon_veloc[2]);
        delta_ridging= std::hypot(divergence_rate,shear_rate/e_factor);

        /*======================================================================
        //! - Updates the ice and snow thickness and ice concentration using a Lagrangian or an Eulerian advection scheme
         *======================================================================
         */

        // We update only elements which have deformed. Not strictly neccesary, but may improve performance.
        bool to_be_updated=false;
        if( divergence_rate!=0.)
            to_be_updated=true;

	/* Important: We don't update elements on the open boundary. This means
         * that ice will flow out as if there was no resistance and in as if the ice
         * state outside the boundary was the same as that inside it. */
        if(std::binary_search(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[0]-1) ||
           std::binary_search(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[1]-1) ||
           std::binary_search(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[2]-1))
            to_be_updated=false;

        // We update only elements where there's ice. Not strictly neccesary, but may improve performance.
        if((M_conc[cpt]>0.)  && (to_be_updated))
        {
            double surf_ratio = this->measure(M_elements[cpt],M_mesh, UM_P) / this->measure(M_elements[cpt],M_mesh,M_UM);

            M_conc[cpt] *= surf_ratio;
            M_thick[cpt] *= surf_ratio;
            M_snow_thick[cpt] *= surf_ratio;
            M_sigma[3*cpt] *= surf_ratio;
            M_sigma[3*cpt+1] *= surf_ratio;
            M_sigma[3*cpt+2] *= surf_ratio;
            M_ridge_ratio[cpt] *= surf_ratio;

            if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                M_h_thin[cpt] *= surf_ratio;
                M_conc_thin[cpt] *= surf_ratio;
                M_hs_thin[cpt] *= surf_ratio;
            }
        }

        /*======================================================================
        //! - Performs the mechanical redistribution (after the advection the concentration can be higher than 1, meaning that ridging should have occured)
         *======================================================================
         */
        double open_water_concentration=1.-M_conc[cpt];

        /* Thin ice category */
        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            open_water_concentration -= M_conc_thin[cpt];
        }

        // ridging scheme
        double opening_factor=((1.-M_conc[cpt])>G_star) ? 0. : std::pow(1.-(1.-M_conc[cpt])/G_star,2.);

        // limit open_water concentration to 0.
        open_water_concentration=(open_water_concentration<0.)?0.:open_water_concentration;

        // opening factor set to 0 for viscous case.
        opening_factor=(young>0.) ? opening_factor : 0.;

        //open_water_concentration += time_step*0.5*(delta_ridging-divergence_rate)*opening_factor;
        open_water_concentration += time_step*0.5*shear_rate/e_factor*opening_factor;

        // limit open_water concentration to 1.
        open_water_concentration=(open_water_concentration>1.)?1.:open_water_concentration;

        /* Thin ice category */
        double new_conc_thin=0.;
        double new_h_thin=0.;
        double new_hs_thin=0.;

        double newice = 0.;
        double del_c = 0.;
        double newsnow = 0.;

        double ridge_thin_ice_aspect_ratio=10.;

#if 1
        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            if(M_conc_thin[cpt]>0. )
            {
                new_conc_thin   = std::min(1.,std::max(1.-M_conc[cpt]-open_water_concentration,0.));

                // Ridging
                if( (M_conc[cpt] > vm["dynamics.min_c"].as<double>()) && (M_thick[cpt] > vm["dynamics.min_h"].as<double>()) && (new_conc_thin < M_conc_thin[cpt] ))
                {
                    new_h_thin      = new_conc_thin*M_h_thin[cpt]/M_conc_thin[cpt]; // so that we keep the same h0, no preferences for the ridging
                    new_hs_thin     = new_conc_thin*M_hs_thin[cpt]/M_conc_thin[cpt];

                    newice = M_h_thin[cpt]-new_h_thin;
                    del_c   = (M_conc_thin[cpt]-new_conc_thin)/ridge_thin_ice_aspect_ratio;
                    newsnow = M_hs_thin[cpt]-new_hs_thin;

                    M_h_thin[cpt]   = new_h_thin;
                    M_hs_thin[cpt]  = new_hs_thin;

                    M_thick[cpt]        += newice;
                    M_conc[cpt]         += del_c;
                    M_conc[cpt] = std::min(1.,std::max(M_conc[cpt],0.));

                    M_snow_thick[cpt]   += newsnow;

                    if( newice>0. )
                        M_ridge_ratio[cpt]=std::max(0.,std::min(1.,(M_ridge_ratio[cpt]*(M_thick[cpt]-newice)+newice)/M_thick[cpt]));
                }

                M_conc_thin[cpt] = new_conc_thin;
            }
            else
            {
                M_conc_thin[cpt]=0.;
                M_h_thin[cpt]=0.;
                M_hs_thin[cpt]=0.;
            }
        }
#endif

        double new_conc=std::min(1.,std::max(1.-M_conc_thin[cpt]-open_water_concentration+del_c,0.));

        if((new_conc+M_conc_thin[cpt])>1.)
            new_conc=1.-M_conc_thin[cpt];

        if(new_conc<M_conc[cpt])
        {
            M_ridge_ratio[cpt]=std::max(0.,std::min(1.,(M_ridge_ratio[cpt]+(1.-M_ridge_ratio[cpt])*(M_conc[cpt]-new_conc)/M_conc[cpt])));
        }
        M_conc[cpt]=new_conc;

        double max_true_thickness = 50.;
        if(M_conc[cpt]>0.)
        {
            double test_h_thick=M_thick[cpt]/M_conc[cpt];
            test_h_thick = (test_h_thick>max_true_thickness) ? max_true_thickness : test_h_thick ;
            M_conc[cpt]=std::min(1.-M_conc_thin[cpt],M_thick[cpt]/test_h_thick);
        }
        else
        {
            M_ridge_ratio[cpt]=0.;
            M_thick[cpt]=0.;
            M_snow_thick[cpt]=0.;
        }

        // END: Ridging scheme and mechanical redistribution

        /*======================================================================
         //! - Updates the internal stress
         *======================================================================
         */

        if( (M_conc[cpt] > vm["dynamics.min_c"].as<double>()) && (M_thick[cpt] > vm["dynamics.min_h"].as<double>()) && (young>0.))
        {

            double damaging_exponent = ridging_exponent;
            double undamaged_time_relaxation_sigma=vm["dynamics.undamaged_time_relaxation_sigma"].as<double>();
            double exponent_relaxation_sigma=vm["dynamics.exponent_relaxation_sigma"].as<double>();

            double time_viscous=undamaged_time_relaxation_sigma*std::pow(1.-old_damage,exponent_relaxation_sigma-1.);
            double multiplicator=time_viscous/(time_viscous+time_step);

        for(int i=0;i<3;i++)
        {
            sigma_dot_i = 0.0;
            for(int j=0;j<3;j++)
            {
                sigma_dot_i += std::exp(damaging_exponent*(1.-M_conc[cpt]))*young*(1.-old_damage)*M_Dunit[3*i + j]*epsilon_veloc[j];
            }

            sigma_pred[i] = (M_sigma[3*cpt+i]+4.*time_step*sigma_dot_i)*multiplicator;
            sigma_pred[i] = (M_conc[cpt] > vm["dynamics.min_c"].as<double>()) ? (sigma_pred[i]):0.;

            M_sigma[3*cpt+i] = (M_sigma[3*cpt+i]+time_step*sigma_dot_i)*multiplicator;
            M_sigma[3*cpt+i] = (M_conc[cpt] > vm["dynamics.min_c"].as<double>()) ? (M_sigma[3*cpt+i]):0.;
        }

        /*======================================================================
         //! - Estimates the level of damage from the updated internal stress and the local damage criterion
         *======================================================================
         */

        /* Compute the shear and normal stress, which are two invariants of the internal stress tensor */

        sigma_s = std::hypot((sigma_pred[0]-sigma_pred[1])/2.,sigma_pred[2]);
        sigma_n =-          (sigma_pred[0]+sigma_pred[1])/2.;

        sigma_1 = sigma_n+sigma_s; // max principal component following convention (positive sigma_n=pressure)
        sigma_2 = sigma_n-sigma_s; // max principal component following convention (positive sigma_n=pressure)

        double ridge_to_normal_cohesion_ratio=vm["dynamics.ridge_to_normal_cohesion_ratio"].as<double>();
        double norm_factor=vm["dynamics.cohesion_thickness_normalisation"].as<double>();
        double exponent=vm["dynamics.cohesion_thickness_exponent"].as<double>();

        double hi=0.;
        if(M_conc[cpt]>0.1)
            hi = M_thick[cpt]/M_conc[cpt];
        else
            hi = M_thick[cpt]/0.1;

        //* REMOVE THIS: SHOULD NOT NEED TO SCALE COHESION WITH THICKNESS OF ICE
        double mult_factor = std::pow(hi/norm_factor,exponent)*(1. + M_ridge_ratio[cpt]*(ridge_to_normal_cohesion_ratio-1.) );

        double effective_cohesion = mult_factor * M_Cohesion[cpt];
        double effective_compressive_strength = mult_factor * M_Compressive_strength[cpt];

        q = std::pow(std::pow(std::pow(tan_phi,2.)+1,.5)+tan_phi,2.);
        sigma_c=2.*effective_cohesion/(std::pow(std::pow(tan_phi,2.)+1,.5)-tan_phi);
        sigma_t=-sigma_c/q;
        tract_max=-tract_coef*effective_cohesion/tan_phi; /* minimum and maximum normal stress */
            
            
        /* Calculate the adjusted level of damage */
            //! \warning{sigma_target is actually not effective: critical states of stress are not projected back onto the damage envelope.}
        if(sigma_n>effective_compressive_strength)
        {
            sigma_target=effective_compressive_strength;

            tmp=1.0-sigma_target/sigma_n*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }

        if((sigma_1-q*sigma_2)>sigma_c)
        {
            sigma_target = sigma_c;

            tmp = 1.0-sigma_target/(sigma_1-q*sigma_2)*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt] = tmp;
            }
        }

        if(sigma_n<tract_max)
        {
            sigma_target = tract_max;

            tmp = 1.0-sigma_target/sigma_n*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt] = tmp;
            }
        }

        }
        else // if M_conc or M_thick too low, set sigma to 0.
        {
            for(int i=0;i<3;i++)
            {
                M_sigma[3*cpt+i] = 0.;
            }
        }

        /*======================================================================
         * Check:
         *======================================================================
         */

        /* lower bounds */
        M_conc[cpt] = ((M_conc[cpt]>0.)?(M_conc[cpt] ):(0.)) ;
        M_thick[cpt]        = ((M_thick[cpt]>0.)?(M_thick[cpt]     ):(0.)) ;
        M_snow_thick[cpt]   = ((M_snow_thick[cpt]>0.)?(M_snow_thick[cpt]):(0.)) ;

        /* Ice damage
         * We use now a constant healing rate defined as 1/time_recovery_damage
         * so that we are now able to reset the damage to 0.
         * otherwise, it will never heal completely.
         * time_recovery_damage still depends on the temperature when themodynamics is activated.
         */
        tmp=M_damage[cpt]-time_step/M_time_relaxation_damage[cpt];
        if(M_thick[cpt]==0.)
            tmp=0.;

        M_damage[cpt]=((tmp>0.)?(tmp):(0.));
    }//loop over elements
}//update


//------------------------------------------------------------------------------------------------------
//! Solves the momentum equation for the sea ice velocity. Called by step(), after the assemble() function.
void
FiniteElement::solve()
{
    M_comm.barrier();
    //SolverPetsc ksp;
    //ksp.solve(M_matrix, M_solution, M_vector);

    //timer["solution"].first.restart();

    M_solver->solve(_matrix=M_matrix,
                    _solution=M_solution,
                    _rhs=M_vector,
                    _rtolerance=1.e-09,
                    _ksp=vm["solver.ksp-type"].as<std::string>()/*"preonly"*/,
                    _pc=vm["solver.pc-type"].as<std::string>()/*"cholesky"*/,
                    _pcfactormatsolverpackage=vm["solver.mat-package-type"].as<std::string>()/*"cholmod"*/,
                    _reuse_prec=vm["solver.ksp-reuse-prec"].as<bool>()/*true*/,
                    _rebuild=M_regrid
                    );


    M_solution->close();

    // if (M_rank==0)
    //     std::cout<<"TIMER SOLUTION= " << timer["solution"].first.elapsed() <<"s\n";

    //M_solution->printMatlab("solution.m");
}//solve


//------------------------------------------------------------------------------------------------------
//! Nests the ice variables (thickness, concentration, snow thickness) from outer domain: used in a regional forecast context.
//! Called by the step() function.
void
FiniteElement::nestingIce()
{
    double const nudge_time  = M_nudge_timescale;
    double const nudge_scale = M_nudge_lengthscale*M_res_root_mesh;

    for (int i=0; i < M_num_elements; ++i)
    {
        double fNudge;
        if ( ( M_nudge_function != "linear" ) && ( M_nudge_function != "exponential" ) )
        {
            throw std::logic_error("M_nudge_function not known: use exponential or linear");
        }
        else
        {
            if ( M_nudge_function == "linear" )
                fNudge = 1. - std::min(1.,(M_nesting_dist_elements[i]/nudge_scale));
            if ( M_nudge_function == "exponential" )
                fNudge = std::exp(-M_nesting_dist_elements[i]/nudge_scale);

            if ( Environment::vm()["thermo.newice_type"].as<int>() == 4 ) {
                M_conc_thin[i]  += (fNudge*(time_step/nudge_time)*(M_ice_conc_thin[i]-M_conc_thin[i]));
                M_h_thin[i]     += (fNudge*(time_step/nudge_time)*(M_ice_h_thin[i]-M_h_thin[i]));
                M_hs_thin[i]    += (fNudge*(time_step/nudge_time)*(M_ice_hs_thin[i]-M_hs_thin[i]));
                M_conc[i]       += (fNudge*(time_step/nudge_time)*(M_ice_conc[i]-M_conc[i]));
                M_thick[i]      += (fNudge*(time_step/nudge_time)*(M_ice_thick[i]-M_thick[i]));
                M_snow_thick[i] +=(fNudge*(time_step/nudge_time)*(M_ice_snow_thick[i]-M_snow_thick[i]));
            }
            else {
                M_conc[i]       += (fNudge*(time_step/nudge_time)*(M_ice_conc[i]-M_conc[i]));
                M_thick[i]      += (fNudge*(time_step/nudge_time)*(M_ice_thick[i]-M_thick[i]));
                M_snow_thick[i] += (fNudge*(time_step/nudge_time)*(M_ice_snow_thick[i]-M_snow_thick[i]));
            }
        }
    }
}//nestingIce

    
//------------------------------------------------------------------------------------------------------
//! Nests the dynamical variable (velocity, damage, stress) from outer domain.
//! Called by the step() function.
void
FiniteElement::nestingDynamics()
{
    double const nudge_time  = M_nudge_timescale;
    double const nudge_scale = M_nudge_lengthscale*M_res_root_mesh;

    for (int i=0; i < M_num_elements; ++i)
    {
        double fNudge;
        if ( ( M_nudge_function != "linear" ) && ( M_nudge_function != "exponential" ) )
        {
            throw std::logic_error("M_nudge_function not known: use exponential or linear");
        }
        else
        {
            if ( M_nudge_function == "linear" )
                fNudge = 1. - std::min(1.,(M_nesting_dist_elements[i]/nudge_scale));
            if ( M_nudge_function == "exponential" )
                fNudge = std::exp(-M_nesting_dist_elements[i]/nudge_scale);
            M_sigma[3*i]     += (fNudge*(time_step/nudge_time)*(M_nesting_sigma1[i]-M_sigma[3*i]));
            M_sigma[3*i+1]   += (fNudge*(time_step/nudge_time)*(M_nesting_sigma2[i]-M_sigma[3*i+1]));
            M_sigma[3*i+2]   += (fNudge*(time_step/nudge_time)*(M_nesting_sigma3[i]-M_sigma[3*i+2]));
            M_damage[i]      += (fNudge*(time_step/nudge_time)*(M_nesting_damage[i]-M_damage[i]));
            M_ridge_ratio[i] += (fNudge*(time_step/nudge_time)*(M_nesting_ridge_ratio[i]-M_ridge_ratio[i]));
        }
    }

    for (int i=0; i < M_num_nodes; ++i)
    {
        double fNudge;
        if ( ( M_nudge_function != "linear" ) && ( M_nudge_function != "exponential" ) )
        {
            throw std::logic_error("M_nudge_function not known: use exponential or linear");
        }
        else
        {
            if ( M_nudge_function == "linear" )
                fNudge = 1. - std::min(1.,(M_nesting_dist_nodes[i]/nudge_scale));
            if ( M_nudge_function == "exponential" )
                fNudge = std::exp(-M_nesting_dist_nodes[i]/nudge_scale);
            M_VT[i]             += (fNudge*(time_step/nudge_time)*(M_nesting_VT1[i]-M_VT[i]));
            M_VT[i+M_num_nodes] += (fNudge*(time_step/nudge_time)*(M_nesting_VT2[i]-M_VT[i+M_num_nodes]));
        }
    }
}//nestingDynamics

    
//------------------------------------------------------------------------------------------------------
//! Performs thermodynamics calculation based on the 1D thermodynamical model.
    
//! \note
//! - Uses either the Winton et al. 2000 or a zero-layer scheme (thermoWinton(), thermoIce0()).
//! - No stability dependent atmospheric drag for now.
//! - There is now only one big loop for the thermodynamics, so that we can use multithreading.
void
FiniteElement::thermo(double dt)
{
    M_comm.barrier();


    // constant variables
    //! 1) Sets local variables to values defined in options.cpp
    double const timeT = vm["thermo.ocean_nudge_timeT"].as<double>(); //! \param timeT (double const) Nudging time for temperature
    double const timeS = vm["thermo.ocean_nudge_timeS"].as<double>(); //! \param timeS (double const) Nudging time for salinity
    double const Qdw_const = vm["ideal_simul.constant_Qdw"].as<double>(); //! \param Qdw_const (double const) Heat flux from ocean nudging
    double const Fdw_const = vm["ideal_simul.constant_Fdw"].as<double>(); //! \param Qdw_const (double const) Fresh water flux from ocean nudging
    double const ocean_albedo = vm["thermo.albedoW"].as<double>(); //! \param Qdw_const (double const) Ocean albedo
    double const drag_ocean_t = vm["thermo.drag_ocean_t"].as<double>(); //! \param drag_ocean_t (double const) Ocean drag parameter, to calculate sensible heat flux
    double const drag_ocean_q = vm["thermo.drag_ocean_q"].as<double>(); //! \param drag_ocean_q (double const) Ocean drag parameter, to calculate latent heat flux
    double const rh0   = 1./vm["thermo.hnull"].as<double>(); //! \param rh0 (double const)
    double const rPhiF = 1./vm["thermo.PhiF"].as<double>(); //! \param rPhiF (double const)
    
    double const qi = physical::Lf * physical::rhoi; //! \param qi (double const) Latent heat of fusion * ice density [J m^{-3}]
    double const qs = physical::Lf * physical::rhos; //! \param qi (double const) Latent heat of fusion * snow density [J m^{-3}]

    int const newice_type = vm["thermo.newice_type"].as<int>(); //! \param newice_type (int const) Type of new ice thermo scheme (4 diff. cases: Hibler 1979, Olason 2009, ...)
    int const melt_type = vm["thermo.melt_type"].as<int>(); //! \param melt_type (int const) Type of melting scheme (2 diff. cases : Hibler 1979, Mellor and Kantha 1989)
    double const PhiM = vm["thermo.PhiM"].as<double>(); //! \param PhiM (double const) Parameter for melting?
    double const PhiF = vm["thermo.PhiF"].as<double>(); //! \param PhiF (double const) Parameter for freezing?
    
    const double aw=6.1121e2, bw=18.729, cw=257.87, dw=227.3; //! \param aw, bw, cw, dw (double const) Constants for the calculation of specific humidity (atmosphere)
    const double Aw=7.2e-4, Bw=3.20e-6, Cw=5.9e-10; //! \param Aw, Bw, Cw (double const) Other set of constants for the calculation of specific humidity (atmosphere)

    const double alpha=0.62197, beta=0.37803; //! \param alpha, beta (double const) Constants for the calculation of specific humidity (at the ocean surface)

    for (int i=0; i < M_num_elements; ++i)
    {
        // -------------------------------------------------
        //! 1.1) Initializes temporary variables

        double  hi=0.;          //! \param hi (double) Ice thickness (slab) [m]
        double  hi_old=0.;      //! \param hi_old (double) Ice thickness at the start of the time step (slab) [m]
        double  hs=0.;          //! \param hs (double) Snow thickness (slab) [m]

        double  hi_thin=0.;     //! \param hi_thin (double) Thin ice thickness (slab) [m]
        double  hi_thin_old=0.; //! \param hi_thin_old (double) Thin ice thickness at the start of the time step (slab) [m]
        double  hs_thin=0.;     //! \param hs_thin (double) Snow thickness on thin ice (slab) [m]

        double  del_hi=0.;      //! \param del_hi (double) Rate of change in ice thickness (slab only) [m/s]
        double  del_hi_thin=0.; //! \param del_hi_thin (double) Rate of change in thin ice thickness (slab only) [m/s]

        double  evap=0.;        //! \param evap (double) Evaporation (rate or amount?)

        double  Qdw=0.;         //! \param Qdw (double) Heat flux from ocean nudging
        double  Fdw=0.;         //! \param Fdw (double) Fresh water flux from ocean nudging

        double  Qio=0.;         //! \param Qio (double) Ice-ocean heat flux
        double  Qio_thin=0.;    //! \param Qio_thin (double) Ice-ocean heat flux through thin ice
        double  Qai=0.;         //! \param Qai (double) Total atmosphere-ice heat flux
        double  Qai_thin=0.;    //! \param Qai_thin (double) Total atmosphere-ice heat flux over thin ice
        double  Qswi=0.;        //! \param Qswi (double) Short-wave atmosphere-ice heat flux
        double  Qsw_thin=0.;    //! \param Qsw_thin (double) Short-wave atmosphere-ice heat flux over thin ice
        double  Qlwi=0.;        //! \param Qlwi (double) Latent atmosphere-ice heat flux
        double  Qlw_thin=0.;    //! \param Qlw_thin (double) Latent atmosphere-ice heat flux over thin ice
        double  Qshi=0.;        //! \param Qshi (double) Sensible atmosphere-ice heat flux
        double  Qsh_thin=0.;    //! \param Qsh_thin (double) Sensible atmosphere-ice heat flux over thin ice
        double  Qlhi=0.;        //! \param Qlhi (double) Long-wave atmosphere-ice heat flux
        double  Qlh_thin=0.;    //! \param Qlh_thin (double) Long-wave atmosphere-ice heat flux over thin ice
        double  Qow=0.;         //! \param Qow (double) Open water heat flux

        //! 1.2) Saves old _volumes_ and concentrations and calculates wind speed
        double  old_vol=M_thick[i];
        double  old_snow_vol=M_snow_thick[i];
        double  old_conc=M_conc[i];

        double  old_h_thin = 0.;
        double  old_hs_thin = 0.;
        double  old_conc_thin=0.;
        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            old_h_thin  = M_h_thin[i];
            old_conc_thin  = M_conc_thin[i];
            old_hs_thin = M_hs_thin[i];
        }

        double sum_u=0.;
        double sum_v=0.;

        for (int j=0; j<3; ++j)
        {
            // calculate wind per element
            sum_u += M_wind[M_elements[i].indices[j]-1];
            sum_v += M_wind[M_elements[i].indices[j]-1+M_num_nodes];
        }
        double  wspeed = std::hypot(sum_u, sum_v)/3.;

        // definition of the snow fall in kg/m^2/s
        double tmp_snowfall;
        if(M_snowfr.isInitialized())
            tmp_snowfall=M_precip[i]*M_snowfr[i];
        else if (M_snowfall.isInitialized())
            tmp_snowfall=M_snowfall[i];
        else
        {
            if(M_tair[i]<0)
                tmp_snowfall=M_precip[i];
            else
                tmp_snowfall=0.;
        }

        double Qsw_in;
        if(M_Qsw_in.isInitialized())
        {
            Qsw_in=M_Qsw_in[i];
        }
        else
        {
            throw std::logic_error(
                    "The function approxSW is not yet implemented, you need to initialize M_Qsw_in");
            //Qsw_in=approxSW();
        }

        double mld=( M_mld[i] > vm["ideal_simul.constant_mld"].as<double>() ) ? M_mld[i] : vm["ideal_simul.constant_mld"].as<double>();

        // -------------------------------------------------
        //! 2) Calculates or sets the flux due to nudging
        if ( M_ocean_type == setup::OceanType::CONSTANT )
        {
            Qdw=Qdw_const;
            Fdw=Fdw_const;
        }
        else
        {
            // nudgeFlux
            if ( M_ocean_salt[i] > physical::si )
            {
                Qdw = -(M_sst[i]-M_ocean_temp[i]) * mld * physical::rhow * physical::cpw/timeT;

                double delS = M_sss[i] - M_ocean_salt[i];
                Fdw = delS * mld * physical::rhow /(timeS*M_sss[i] - dt*delS);
            }
            else
            {
                Qdw = Qdw_const;
                Fdw = Fdw_const;
            }
        }

        // SYL: Calculate the drag coefficients missing??

        // -------------------------------------------------
        //! 3) Calculates fluxes in the open water portion
        //(openWaterFlux in old c code, Qow_mex in matlab)
        double sphuma, sphumw;

        // Calculate atmospheric fluxes

        /* Out-going long-wave flux */
        double Qlw_out = physical::eps*physical::sigma_sb*std::pow(M_sst[i]+physical::tfrwK,4.);

        // -------------------------------------------------
        //! 3.1) Calculates specific humidity of the atmosphere.
        //! - There are two ways to calculate this. We decide which one by checking mixrat - the calling routine must set this to a negative value if the dewpoint should be used.
        if ( M_sphuma.isInitialized() )
        {
            sphuma = M_sphuma[i];
        }
        else if ( M_dair.isInitialized() )
        {
            double fa     = 1. + Aw + M_mslp[i]*1e-2*( Bw + Cw*M_dair[i]*M_dair[i] );
            double esta   = fa*aw*std::exp( (bw-M_dair[i]/dw)*M_dair[i]/(M_dair[i]+cw) );
            sphuma = alpha*fa*esta/(M_mslp[i]-beta*fa*esta) ;
        }
        else if ( M_mixrat.isInitialized() )
        {
            sphuma = M_mixrat[i]/(1.+M_mixrat[i]) ;
        }
        else
        {
            throw std::logic_error("Neither M_sphuma, M_mixrat nor M_dair have been initialized. I cannot calculate sphuma.");
        }

        // -------------------------------------------------
        //! 3.2) Calculates specific humidity at the ocean surface (calcSphumW in matlab)
        double fw     = 1. + Aw + M_mslp[i]*1e-2*( Bw + Cw*M_sst[i]*M_sst[i] );
        double estw   = aw*std::exp( (bw-M_sst[i]/dw)*M_sst[i]/(M_sst[i]+cw) )*(1-5.37e-4*M_sss[i]);
        sphumw = alpha*fw*estw/(M_mslp[i]-beta*fw*estw) ;

        // -------------------------------------------------
        /* Density of air */
        double rhoair = M_mslp[i]/(physical::Ra*(M_tair[i]+physical::tfrwK)) * (1.+sphuma)/(1.+1.609*sphuma);

        /* Sensible heat flux */
        double Qsh_ow = drag_ocean_t*rhoair*physical::cpa*wspeed*( M_sst[i] - M_tair[i] );

        /* Latent heat flux */
        double Lv  = physical::Lv0 - 2.36418e3*M_tair[i] + 1.58927*M_tair[i]*M_tair[i] - 6.14342e-2*std::pow(M_tair[i],3.);
        double Qlh_ow = drag_ocean_q*rhoair*Lv*wspeed*( sphumw - sphuma );

        /* Evaporation */
        evap = Qlh_ow/(physical::rhofw*Lv);

        // Sum them up
        double Qlw_in;
        if(M_Qlw_in.isInitialized())
        {
            Qlw_in=M_Qlw_in[i];
        }
        else
        {
            double tsa = M_tice[0][i] + physical::tfrwK;
            double taa = M_tair[i]  + physical::tfrwK;
            // s.b.idso & r.d.jackson, thermal radiation from the atmosphere, j. geophys. res. 74, 5397-5403, 1969
        	Qlw_in = sigma_sb*std::pow(taa,4) \
        			*( 1. - 0.261*std::exp(-7.77e-4*std::pow(taa-physical::tfrwK,2)) ) \
        			*( 1. + 0.275*M_tcc[i] );
        }

        // Qow>0 => flux out of ocean:
        // - subtract shortwave and longwave input;
        // add heat loss from longwave radiation, sensible heat loss (temp changes)
        // and evaporation (latent heat loss - temp doesn't change, but phase changes)
        double Qsw_ow = -Qsw_in*(1.-ocean_albedo);
        double Qlw_ow = -Qlw_in + Qlw_out;
        Qow = Qsw_ow + Qlw_ow + Qsh_ow + Qlh_ow;

        // -------------------------------------------------
        //! 4) Calculates the thickness change of the ice slab (thermoIce0 in matlab)

        switch ( M_thermo_type )
        {
            case setup::ThermoType::ZERO_LAYER:
                this->thermoIce0(i, dt, wspeed, sphuma, M_conc[i], M_thick[i], M_snow_thick[i],
                        Qlw_in, Qsw_in, mld, tmp_snowfall,//end of inputs - rest are outputs
                        hi, hs, hi_old, Qio, del_hi, M_tice[0][i],
                        Qai, Qswi, Qlwi, Qshi, Qlhi);
                break;
            case setup::ThermoType::WINTON:
                this->thermoWinton(i, dt, wspeed, sphuma, M_conc[i], M_thick[i], M_snow_thick[i],
                        Qlw_in, Qsw_in, mld, tmp_snowfall,//end of inputs - rest are outputs
                        hi, hs, hi_old, Qio, del_hi,
                        M_tice[0][i], M_tice[1][i], M_tice[2][i],
                        Qai, Qswi, Qlwi, Qshi, Qlhi);
                break;
        }

        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            this->thermoIce0(i, dt, wspeed, sphuma, old_conc_thin, M_h_thin[i], M_hs_thin[i],
                    Qlw_in, Qsw_in, mld, tmp_snowfall, hi_thin, hs_thin, hi_thin_old, Qio_thin, del_hi_thin, M_tsurf_thin[i],
                        Qai_thin, Qsw_thin, Qlw_thin, Qsh_thin, Qlh_thin);
            M_h_thin[i]  = hi_thin * old_conc_thin;
            M_hs_thin[i] = hs_thin * old_conc_thin;
        }

        // -------------------------------------------------
        //! 5) Calculates the ice growth over open water and lateral melt (thermoOW in matlab)

        /* Local variables */
        double tw_new, tfrw, newice, del_c, newsnow, h0;

        /* dT/dt due to heatflux ocean->atmosphere */
        tw_new = M_sst[i] - Qow*dt/(mld*physical::rhow*physical::cpw);
        tfrw   = -physical::mu*M_sss[i];

        /* Form new ice in case of super cooling, and reset Qow and evap */
        if ( tw_new < tfrw )
        {
            newice  = (1.-M_conc[i]-M_conc_thin[i])*(tfrw-tw_new)*mld*physical::rhow*physical::cpw/qi;// m
            Qow  = -(tfrw-M_sst[i])*mld*physical::rhow*physical::cpw/dt;
            // evap = 0.;
        }
        else
        {
            newice  = 0.;
        }

        /* Decide the change in ice fraction (del_c) */
        /* Initialise to be safe */
        del_c = 0.;
        newsnow = 0.;

        /* Freezing conditions */
        switch ( newice_type )
        {
            case 1:
                /* Hibler's (79) approach */
                del_c = newice*rh0;
                break;
            case 2:
                /* Mellor and Kantha (89) */
                if ( hi_old > 0. )
                {
                    del_c = newice*PhiF/hi_old;
                }
                else
                {
                    if ( newice > 0. )
                        del_c = 1.;
                    else
                        del_c = 0.;
                }
                break;
            case 3:
                /* Olason and Harms (09) */
                h0    = (1.+0.1*wspeed)/15.;
                del_c = newice/std::max(rPhiF*hi_old,h0);
                break;
            case 4:
                /* Thin ice category */
                M_h_thin[i]+=newice;
                M_conc_thin[i]=std::min(1.-M_conc[i],M_conc_thin[i]+newice/h_thin_min);
                newice  = 0.;
                newsnow = 0.;

                if(M_conc_thin[i]>0.)
                {
                    /* Two cases: Thin ice fills the cell or not */
                    if ( M_h_thin[i] < h_thin_min*M_conc_thin[i] )
                    {
                        M_conc_thin[i] = M_h_thin[i]/h_thin_min;
                    }
                    else
                    {
                        h0 = h_thin_min + 2.*(M_h_thin[i]-h_thin_min*M_conc_thin[i])/(M_conc_thin[i]);
                        if(h0>h_thin_max)
                        {
                            del_c = M_conc_thin[i]/(h0-h_thin_min) * (h0-h_thin_max);
                            double del_h_thin = del_c*(h0+h_thin_max)/2.;
                            double del_hs_thin = del_c*M_hs_thin[i]/M_conc_thin[i];

                            M_thick[i] += del_h_thin;
                            // M_conc[i]  += del_c; ; <- this is done properly below

                            newice  = del_h_thin; // Reset newice to use below
                            newsnow = del_hs_thin;
                            // M_snow_thick[i] += newsnow; <- this is done properly below

                            M_conc_thin[i] -= del_c;
                            M_h_thin[i]    -= del_h_thin;
                            M_hs_thin[i]   -= del_hs_thin;
                        }
                    }
                }
                else // we should not have thin ice, no space for it
                {
                    M_thick[i] += M_h_thin[i];

                    newice  = M_h_thin[i];
                    newsnow = M_hs_thin[i];

                    M_h_thin[i] = 0.;
                    M_hs_thin[i]= 0.;
                }
                break;
            default:
                std::cout << "newice_type = " << newice_type << "\n";
                throw std::logic_error("Wrong newice_type");
            }
            /* Check bounds on del_c */
            del_c = std::min( 1.-M_conc[i], del_c );

        if ( del_hi < 0. )
        {
            /* Melting conditions */
            switch ( melt_type )
            {
                case 1:
                    /* Hibler (79) using PhiM to tune. PhiM = 0.5 is
                     * equivalent to Hibler's (79) approach */
                    if ( M_conc[i] < 1. )
                        del_c += del_hi*M_conc[i]*PhiM/hi_old;
                    else
                        del_c += 0.;
                    break;
                case 2:
                    /* Mellor and Kantha (89) */
                    if ( hi > 0. )
                    {
                        /* Use the fraction PhiM of (1-c)*Qow to melt laterally */
                        del_c += PhiM*(1.-M_conc[i])*std::min(0.,Qow)*dt/( hi*qi+hs*qs );
                        /* Deliver the fraction (1-PhiM) of Qow to the ocean */
                        Qow = (1.-PhiM)*Qow;
                    }
                    else
                    {
                        del_c = -M_conc[i];
                    }
                    // This is handled below
                    // /* + Deliver excess energy to the ocean when there's no ice left */
                    //         + std::min(0., std::max(0.,M_conc[i]+del_c)*( hi*qi+hs*qs )/dt);
                    // /* Don't suffer negative c! */
                    // del_c = std::max(del_c, -M_conc[i]);
                    break;
                default :
                    std::cout << "newice_type = " << newice_type << "\n";
                    throw std::logic_error("Wrong newice_type");
            }
        }

        /* New concentration */
        M_conc[i] = M_conc[i] + del_c;

        /* New thickness */
        /* We conserve volume and energy */
        if ( M_conc[i] >= physical::cmin )
        {
            hi = ( hi*old_conc + newice )/M_conc[i];
            if ( del_c < 0. )
            {
                /* We conserve the snow height, but melt away snow as the concentration decreases */
                Qow = Qow + del_c*hs*qs/dt;
            }
            else
            {
                /* Snow volume is conserved as concentration increases */
                hs  = ( hs*old_conc + newsnow )/M_conc[i];
            }

            if ( M_thermo_type == setup::ThermoType::WINTON )
            {
                // Add newice evenly to both layers and recalculate temperature
                double f1    = M_thick[i]/(M_thick[i]+newice); // Fraction of old ice (as opposed to newice) in the upper layer
                double Tbar  = f1*( M_tice[1][i] - physical::Lf*physical::mu*physical::si/(physical::C*M_tice[1][i]) ) + (1-f1)*tfrw; // (39)
                M_tice[1][i] = ( Tbar - std::sqrt(Tbar*Tbar + 4*physical::mu*physical::si*physical::Lf/physical::C) )/2.; // (38)
                M_tice[2][i] = f1*M_tice[2][i] + (1-f1)*tfrw; // (26) slightly rewritten
            }
        }
#if 1
        /* Check limits */
        if ( (M_conc[i] < physical::cmin) || (hi < physical::hmin) )
        {
            // Extract heat from the ocean corresponding to the heat in the
            // remaining ice and snow
            Qow    = Qow + M_conc[i]*hi*qi/dt + M_conc[i]*hs*qs/dt;
            M_conc[i]  = 0.;

            for (int j=0; j<M_tice.size(); j++)
                M_tice[j][i] = tfrw;

            //M_tsurf_thin[i] = tfrw;
            hi     = 0.;
            hs     = 0.;
        }
#endif
        // -------------------------------------------------
        //! 6) Calculates effective ice and snow thickness
        M_thick[i] = hi*M_conc[i];
        M_snow_thick[i] = hs*M_conc[i];

        // -------------------------------------------------
        //! 7) Applies slab Ocean model
        // (slabOcean in matlab)

        // local variables
        double del_vi;      // Change in ice volume
        double del_vs;      // Change in snow volume
        double rain;        // Liquid precipitation
        double emp;         // Evaporation minus liquid precipitation
        double Qio_mean;    // Element mean ice-ocean heat flux
        double Qow_mean;    // Element mean open water heat flux

        //! * Calculates changes in ice and snow volumes to calculate salt rejection
        del_vi = M_thick[i] - old_vol + M_h_thin[i] - old_h_thin;
        del_vs = M_snow_thick[i] - old_snow_vol + M_hs_thin[i] - old_hs_thin;

        // Rain falling on ice falls straight through. We need to calculate the
        // bulk freshwater input into the entire cell, i.e. everything in the
        // open-water part plus rain in the ice-covered part.
        rain = (1.-old_conc-old_conc_thin)*M_precip[i] + (old_conc+old_conc_thin)*(M_precip[i]-tmp_snowfall);
        emp  = (evap*(1.-old_conc-old_conc_thin)-rain);

        Qio_mean = Qio*old_conc + Qio_thin*old_conc_thin;
        Qow_mean = Qow*(1.-old_conc-old_conc_thin);

        /* Heat-flux */
        M_sst[i] = M_sst[i] - dt*( Qio_mean + Qow_mean - Qdw )/(physical::rhow*physical::cpw*mld);

        /* Change in salinity */
        double denominator= ( mld*physical::rhow - del_vi*physical::rhoi - ( del_vs*physical::rhos + (emp-Fdw)*dt) );
        denominator = ( denominator > 1.*physical::rhow ) ? denominator : 1.*physical::rhow;

        double sss_old = M_sss[i];
        M_sss[i] = M_sss[i] + ( (M_sss[i]-physical::si)*physical::rhoi*del_vi + M_sss[i]*(del_vs*physical::rhos + (emp-Fdw)*dt) )
            / denominator;

        // -------------------------------------------------
        //! 8) Damage manipulation

        // local variables
        double deltaT;      // Temperature difference between ice bottom and the snow-ice interface

        //! * Newly formed ice is undamaged and unridged: Hence calculates damage and ridge ratio as a weighted average of the old damage - ridge ratio and 0, weighted with volume.
        if ( M_thick[i] > old_vol )
        {
            M_damage[i] = M_damage[i]*old_vol/M_thick[i];
            M_ridge_ratio[i] = M_ridge_ratio[i]*old_vol/M_thick[i];
        }

        if ( vm["dynamics.use_temperature_dependent_healing"].as<bool>() )
        {
            //! * Sets time_relaxation_damage to be inversely proportional to the temperature difference between bottom and snow-ice interface
            if ( M_thick[i] > 0. )
            {
                double Tbot = -physical::mu*M_sss[i];
                double C;
                switch (M_thermo_type)
                {
                    case (setup::ThermoType::ZERO_LAYER):
                        C = physical::ki*M_snow_thick[i]/(physical::ks*M_thick[i]);
                        deltaT = std::max(1e-36, Tbot - M_tice[0][i] ) / ( 1. + C );
                        break;
                    case (setup::ThermoType::WINTON):
                        C = physical::ki*M_snow_thick[i]/(physical::ks*M_thick[i]/4.);
                        deltaT = std::max(1e-36, Tbot + C*(Tbot-M_tice[1][i]) - M_tice[0][i] ) / ( 1. + C );
                        break;
                    default:
                        std::cout << "thermo_type= " << (int)M_thermo_type << "\n";
                        throw std::logic_error("Wrong thermo_type");
                }
                M_time_relaxation_damage[i] = std::max(time_relaxation_damage*deltaT_relaxation_damage/deltaT, dt);
            }
            else
            {
                M_time_relaxation_damage[i] = 1e36;
            }
        }
        // -------------------------------------------------

        //! 9) Computes diagnostics (open water fraction and heat fluxes to the atmosphere and ocean)
        double ow_fraction = 1. - old_conc - old_conc_thin;

        // Total heat flux to the atmosphere
        D_Qa[i] = Qai*old_conc + Qai_thin*old_conc_thin + (Qsw_ow+Qlw_ow+Qsh_ow+Qlh_ow)*ow_fraction;

        // Short wave flux to the atmosphere
        D_Qsw[i] = Qswi*old_conc + Qsw_thin*old_conc_thin + Qsw_ow*ow_fraction;

        // Long wave flux to the atmosphere
        D_Qlw[i] = Qlwi*old_conc + Qlw_thin*old_conc_thin + Qlw_ow*ow_fraction;

        // Sensible heat flux to the atmosphere
        D_Qsh[i] = Qshi*old_conc + Qsh_thin*old_conc_thin + Qsh_ow*ow_fraction;

        // Latent heat flux to the atmosphere
        D_Qlh[i] = Qlhi*old_conc + Qlh_thin*old_conc_thin + Qlh_ow*ow_fraction;

        // Total heat lost by ocean
        D_Qo[i] = Qio_mean + Qow_mean;

        // Salt release into the ocean - kg/day
        D_delS[i] = (M_sss[i] - sss_old)*physical::rhow*mld/dt;
    }// end for loop
}//thermo
    

//------------------------------------------------------------------------------------------------------
//! Calculates atmospheric fluxes through bulk formula.
//! Called by the thermoWinton() and thermoIce0() functions.
void
FiniteElement::atmFluxBulk(int i, double Tsurf, double sphuma, double drag_ice_t, double Qsw, double Qlw_in, double wspeed,
        double &Qai, double &dQaidT, double &subl,
        double &Qsh, double &Qlh, double &Qlw)
{
    // Constants
    const double ai=6.1115e2, bi=23.036, ci=279.82, di=333.7;
    const double Ai=2.2e-4, Bi=3.83e-6, Ci=6.4e-10;

    const double alpha=0.62197, beta=0.37803;

    // -------------------------------------------------
    // 4.1) BULK ICE
    /* Out-going long-wave flux and derivative */
    double Qlw_out =   physical::eps * physical::sigma_sb * std::pow(Tsurf+physical::tfrwK,4);
    double dQlwdT  = 4.*physical::eps * physical::sigma_sb * std::pow(Tsurf+physical::tfrwK,3);

    // -------------------------------------------------
    // calcSphumI
    /* Specific humidity - ice surface */
    double fi      = 1. + Ai + M_mslp[i]*1e-2*( Bi + Ci*Tsurf*Tsurf );
    double esti    = ai*std::exp( (bi-Tsurf/di)*Tsurf/(Tsurf+ci) );
    double sphumi = alpha*fi*esti/(M_mslp[i]-beta*fi*esti);

    /* We need the derivative of sphumi wrt. tsurf */
    double dsphumdesti = alpha/(M_mslp[i]-beta*fi*esti)*( 1. + beta*fi*esti/(M_mslp[i]-beta*fi*esti) );
    double destidT     = ( bi*ci*di-Tsurf*( 2.*ci+Tsurf) )/( di*std::pow(ci+Tsurf,2) )*esti;
    double dfidT       = 2.*Ci*Bi*Tsurf;
    double dsphumidT   = dsphumdesti*(fi*destidT+esti*dfidT);

    // -------------------------------------------------

    /* Density of air */
    double tairK  = M_tair[i] + physical::tfrwK;
    double rhoair = M_mslp[i]/(physical::Ra*tairK) * (1.+sphuma)/(1.+1.609*sphuma);

    /* Sensible heat flux and derivative */
    Qsh    = drag_ice_t * rhoair * physical::cpa * wspeed*( Tsurf - M_tair[i] );
    double dQshdT = drag_ice_t * rhoair * physical::cpa * wspeed;

    /* Latent heat flux and derivative */
    Qlh    = drag_ice_t*rhoair*(physical::Lf+physical::Lv0)*wspeed*( sphumi - sphuma );
    double dQlhdT = drag_ice_t*(physical::Lf+physical::Lv0)*rhoair*wspeed*dsphumidT;

    /* Sum them up */
    dQaidT = dQlwdT + dQshdT + dQlhdT;

    /* Sublimation */
    subl    = Qlh/(physical::Lf+physical::Lv0);

    /* Sum them up */
    Qlw = Qlw_out - Qlw_in;
    Qai = Qsw + Qlw + Qsh + Qlh;
}//atmFluxBulk
    

//------------------------------------------------------------------------------------------------------
//! Calculates ice-ocean heat fluxes.
//! Called by the thermoWinton() and thermoIce0() functions.
double
FiniteElement::iceOceanHeatflux(int cpt, double sst, double sss, double mld, double dt)
{
    //! - Uses all of the excess heat to melt or grow ice. This is not accurate, but it will have to do for now!
    double const Tbot = -physical::mu*sss; // Temperature at ice base (bottom), also freezing point of sea-water
    if ( vm["thermo.Qio-type"].as<std::string>() == "basic" )
    {
        return (sst-Tbot)*physical::rhow*physical::cpw*mld/dt;
    } else if ( vm["thermo.Qio-type"].as<std::string>() == "exchange" ) {
        double welt_oce_ice = 0.;
        for (int i=0; i<3; ++i)
        {
            int nind = (M_elements[cpt]).indices[i]-1;
            welt_oce_ice += std::hypot(M_VT[nind]-M_ocean[nind],M_VT[nind+M_num_nodes]-M_ocean[nind+M_num_nodes]);
        }
        double norm_Voce_ice = welt_oce_ice/3.;
        double Csens_io = 1e-3;
        return (sst-Tbot)*norm_Voce_ice*Csens_io*physical::rhow*physical::cpw;
    } else {
        std::cout << "Qio-type = " << vm["thermo.Qio-type"].as<std::string>() << "\n";
        throw std::logic_error("Wrong Qio-type");
    }
}//iceOceanHeatflux


//------------------------------------------------------------------------------------------------------
//! Calculates the surface albedo. Called by the thermoWinton() function.
//! - Different schemes can be implemented, e.g., Semtner 1976, Untersteiner 1971, CCSM3, ...
double
FiniteElement::albedo(int alb_scheme, double Tsurf, double hs, double alb_sn, double alb_ice, double I_0)
{
    double albedo;

    /* Calculate albedo - we can impliment different schemes if we want */
    switch ( alb_scheme )
    {
        case 1:
        case 2:
            /* This scheme mimics Semtner 76 and Maykut and Untersteiner 71 when
             * alb_ice = 0.64 and alb_sn = 0.85 */
            if ( hs > 0. )
            {
                /* decrease the albedo towards ice albedo for snow thinner than 0.2 m */
                if ( alb_scheme == 2 )
                    albedo = std::min(alb_sn, alb_ice + (alb_sn-alb_ice)*hs/0.2);
                else
                    albedo = alb_sn;
            } else {
                /* account for penetrating shortwave radiation */
                albedo = alb_ice + 0.4*( 1.-alb_ice )*I_0;
            }
            break;
        case 3:
            { // create a scope for abli, albs, and frac_sn
            /* Albedo scheme from ccsm 3 */
            /* The scheme is simplified using the assumption that visible solar
             * radiation accounts for 52% of the spectrum and the near-infrared for
             * 48% (the same split as used in cice when running in stand-alone
             * mode). */

            /* This is the ccsm3 scheme when alb_ice = 0.538 and alb_sn = 0.8256 */
            double albi, albs;
            if ( Tsurf > -1. )
            {
                albi = alb_ice - 0.075*(Tsurf+1.);
                albs = alb_sn  - 0.124*(Tsurf+1.);
            } else {
                albi = alb_ice;
                albs = alb_sn;
            }

            /* Snow cover fraction */
            double frac_sn = hs/(hs+0.02);

            /* Final albedo */
            albedo = frac_sn*albs + (1.-frac_sn)*albi;

            break;
            }
        default:
            std::cout << "alb_scheme = " << alb_scheme << "\n";
            throw std::logic_error("Wrong albedo_scheme");
    }

    return albedo;
}//albedo


//------------------------------------------------------------------------------------------------------
//! Caculates heat fluxes through the ice according to the Winton scheme (ice temperature, growth, and melt).
//! Called by the thermo() function.
void
FiniteElement::thermoWinton(int i, double dt, double wspeed, double sphuma, double conc, double voli, double vols,
        double Qlw_in, double Qsw_in, double mld, double snowfall,
        double &hi, double &hs, double &hi_old, double &Qio, double &del_hi, double &Tsurf, double &T1, double &T2,
        double &Qai, double &Qsw, double &Qlw, double &Qsh, double &Qlh)
{
    // Constants
    double const alb_ice = vm["thermo.alb_ice"].as<double>();
    double const alb_sn  = vm["thermo.alb_sn"].as<double>();
    double const I_0     = vm["thermo.I_0"].as<double>();
    int const alb_scheme = vm["thermo.alb_scheme"].as<int>();

    double const drag_ice_t = vm["thermo.drag_ice_t"].as<double>();

    bool const flooding = vm["thermo.flooding"].as<bool>();

    // Useful volumetric quantities
    double const qi   = physical::Lf * physical::rhoi;
    double const qs   = physical::Lf * physical::rhos;
    double const Crho = physical::C * physical::rhoi;

    double const Tbot     = -physical::mu*M_sss[i]; // Temperature at ice base (bottom), also freezing point of sea-water
    double const Tfr_ice  = -physical::mu*physical::si; // Freezing point of ice

    /* Local variables */
    double dQaidT, subl;

    /* Don't do anything if there's no ice */
    if ( conc <=0. || voli<=0.)
    {
        hi       = 0.;
        hs       = 0.;
        hi_old   = 0.;
        Qio      = 0.;
        Qai      = 0.;
        del_hi   = 0.;
        Tsurf    = Tbot;
        T1       = Tbot;
        T2       = Tbot;
        Qsw      = 0.;
        Qlw      = 0.;
        Qsh      = 0.;
        Qlh      = 0.;
    } else {
        /* Calculate the slab thickness */
        hi     = voli/conc;
        hi_old = hi;
        hs     = vols/conc;

        double const Tfr_surf = ( hs > 0 ) ? 0. : Tfr_ice; // Freezing point at surface (snow or ice)

        /*
         * Internal temperatures
         * Following Winton (2000) we start by calculating the temperature - first T1, then Tsurf, and finally T2
         * Numers in parentheses refer to equations in the paper
         */

        /* Calculate atmospheric fluxes */
        // Shortwave is modulated by the albedo
        Qsw = -Qsw_in*(1.-this->albedo(alb_scheme, Tsurf, hs, alb_sn, alb_ice, I_0))*(1.-I_0);
        // The rest is calculated by bulk formula
        this->atmFluxBulk(i, Tsurf, sphuma, drag_ice_t, Qsw, Qlw_in, wspeed, Qai, dQaidT,subl,
                Qsh, Qlh, Qlw);

        // First some coefficients based on temperatures from the previous time step
        double K12 = 4*physical::ki*physical::ks / ( physical::ks*hi + 4*physical::ki*hs ); // (5)
        double A   = Qai - Tsurf*dQaidT; // (7)
        double B   = dQaidT; // (8)
        double K32 = 2*physical::ki/hi; // (10)

        double A1 = hi*Crho/(2*dt) + K32*( 4*dt*K32 + hi*Crho ) / ( 6*dt*K32 + hi*Crho ) + K12*B/(K12+B); // (16)
        double B1 = -hi/(2*dt) * ( Crho*T1 + qi*Tfr_ice/T1 ) - I_0*Qsw
            - K32*( 4*dt*K32*Tbot + hi*Crho*T2 ) / ( 6*dt*K32 + hi*Crho ) + A*K12/(K12+B); // (17)
        double C1 = hi*qi*Tfr_ice/(2*dt); // (18)

        // Then calculate T1 and Tsurf
        T1    = - ( B1 + std::sqrt(B1*B1-4*A1*C1) ) / ( 2*A1 ); // (21)
        Tsurf = ( K12*T1 - A ) / ( K12 + B ); // (6)

        // Recalculate T1 and Tsurf if we're melting at the surface
        double Msurf = 0.;
        if ( Tsurf > Tfr_surf )
        {
            Tsurf = Tfr_surf;
            // A1 = hi*Crho/(2*dt) + K32*( 4*dt*K32 + hi*Crho ) / ( 6*dt*K32 + hi*Crho ) + K12; // (19)
            A1   += K12 - K12*B/(K12+B);
            // B1 = -hi/(2*dt) * ( Crho*T1 + qi*Tfr_ice/T1 ) - I_0*Qsw
            //     - K32 * ( 4*dt*K32*Tbot + hi*Crho*T2 ) / (6*dt*K32 + hi*Crho ) - K12*Tsurf ; // (20)
            B1   -= K12*Tsurf + A*K12/(K12+B);
            T1    = - ( B1 + std::sqrt(B1*B1-4*A1*C1) ) / ( 2*A1 ); // (21)

            // Surface melt
            Msurf = K12*(T1-Tsurf) - (A+B*Tsurf); // (22)
        }

        // Calculate T2 based on the new T1 and T2 from the previous time step
        T2 = ( 2*dt*K32*(T1+2*Tbot) + hi*Crho*T2 ) / ( 6*dt*K32 + hi*Crho ); // (15)

        /*
         * Thickness changes
         */

        double h1 = hi/2.;
        double h2 = hi/2.;
        double E1 = Crho*(T1 - Tfr_ice) - qi*( 1 - Tfr_ice/T1 ); // (1) - but I've multiplied with rhoi, because it's missing in the paper
        double E2 = Crho*(T2 - Tfr_ice) - qi; // (25) - but I've multiplied with rhoi, because it's missing in the paper

        // Snowfall [kg/m^2/s]
        hs += snowfall/physical::rhos*dt;

        // Sublimation at the surface
        if ( subl*dt <= hs*physical::rhos)
            hs -= subl*dt/physical::rhos;
        else if ( subl*dt - hs*physical::rhos <= h1*physical::rhoi )
        {
            h1 -= (subl*dt - hs*physical::rhos)/physical::rhoi;
            hs  = 0.;
        }
        else if ( subl*dt - h1*physical::rhoi - hs*physical::rhos <= h2*physical::rhoi )
        {
            h2 -= (subl*dt - h1*physical::rhoi - hs*physical::rhos)/physical::rhoi;
            h1  = 0.;
            hs  = 0.;
        }
        else
        {
            double ocn_evap_err = ( subl*dt - (h1+h2)*physical::rhoi - hs*physical::rhos )/physical::rhow;
			LOG(WARNING) << "All the ice has sublimated. This shouldn't happen and will result in lack of evaporation from the ocean of "
                << ocn_evap_err*1e3 << " mm over the current time step, in element " << i << ".\n";
            h2 = 0.;
            h1 = 0.;
            hs = 0.;
        }

        // Bottom melt/freezing
        Qio    = this->iceOceanHeatflux(i, M_sst[i], M_sss[i], mld, dt);
        double Mbot  = Qio - 4*physical::ki*(Tbot-T2)/hi; // (23)

        // Growth/melt at the ice-ocean interface
        if ( Mbot <= 0. )
        {
            // Growth
            double Ebot  = Crho*(Tbot - Tfr_ice) - qi; // (25) - but I've multiplied with rhoi, because it's missing in the paper
            double delh2 = Mbot*dt/Ebot; // (24)
            T2 = ( delh2*Tbot + h2*T2 ) / ( delh2 + h2 ); // (26)
            h2 += delh2;
        } else {
            // Melt
            double delh2 = -std::min(          -  Mbot*dt/E2,                        h2); // (31) - with added division with rhoi
            double delh1 = -std::min(std::max( -( Mbot*dt + E2*h2 )/E1,         0.), h1); // (32) - with added division with rhoi
            double delhs = -std::min(std::max(  ( Mbot*dt + E2*h2 + E1*h1 )/qs, 0.), hs); // (32) - with added division with rhoi and rhos

            // If everyting melts we need to give back to the ocean
            if ( h2+h1+hs -delh2-delh1-delhs <= 0. )
                Qio -= std::max(Mbot*dt - qs*hs + E1*h1 + E2*h2, 0.)/dt; // (34) - with added multiplication of rhoi and rhos and division with dt

            hs += delhs;
            h1 += delh1;
            h2 += delh2;
        }

        // Melting at the surface
        assert(Msurf >= 0); // Sanity check
        double delhs = -std::min(             Msurf*dt/qs,                          hs); // (27) - with division of rhos
        double delh1 = -std::min(std::max( -( Msurf*dt - qs*hs )/E1,           0.), h1); // (28) - with division of rhoi and rhos
        double delh2 = -std::min(std::max( -( Msurf*dt - qs*hs + E1*h1 ) / E2, 0.), h2); // (29) - with division of rhoi and rhos

        // If everyting melts we need to give back to the ocean
        if ( h2+h1+hs -delh2-delh1-delhs <= 0. )
            Qio -= std::max(Msurf*dt - qs*hs + E1*h1 + E2*h2, 0.)/dt; // (30) - with multiplication of rhoi and rhos and division with dt

        hs += delhs;
        h1 += delh1;
        h2 += delh2;
        hi  = h1 + h2;

        // Snow-to-ice conversion
        double freeboard = ( hi*(physical::rhow-physical::rhoi) - hs*physical::rhos) / physical::rhow;
        if ( flooding && freeboard < 0)
        {
            // double delhs = -std::max( ( hs - (physical::rhow-physical::rhoi)*hi/physical::rhos )*physical::rhoi/physical::rhow, 0.); // (35)
            hs += std::min( freeboard*physical::rhoi/physical::rhos, 0. );
            // double delh1 =  std::max( ( hs - (physical::rhow-physical::rhoi)*hi/physical::rhos )*physical::rhos/physical::rhow, 0.); // (36)
            double delh1 = std::max( -freeboard, 0. );

            double f1   = 1-delh1/(delh1+h1); // Fraction of layer 1 ice in the new upper layer
            double Tbar = f1*( T1 + qi*Tfr_ice/(Crho*T1) ) + (1-f1)*Tfr_ice; // (39)

            T1 = ( Tbar - std::sqrt(Tbar*Tbar - 4*Tfr_ice*qi/Crho) )/2.; // (38)
            h1 += delh1;
        }

        // Even out the layer structure and temperatures
        if ( h2 > h1 )
        {
            // Lower layer ice is added to the upper layer
            // T1 changes, but T2 not
            double f1   = h1/hi*2.; // Fraction of layer 1 ice found in the new layer 1
            double Tbar = f1*( T1 + qi*Tfr_ice/(Crho*T1) ) + (1-f1)*T2; // (39)
            T1 = ( Tbar - std::sqrt(Tbar*Tbar - 4*Tfr_ice*qi/Crho) )/2.; // (38)
        } else {
            // Upper layer ice is added to the lower layer
            // T2 changes, but T1 not
            double f1   = (2.*h1-hi)/hi; // Fraction of layer 1 ice found in new layer 2
            T2 = f1*( T1 + qi*Tfr_ice/(Crho*T1) ) + (1-f1)*T2; // (40)

            // Melt from top and bottom if T2 is too high
            if ( T2 > Tfr_ice )
            {
                // This is:
                // hi -= h2*C*(T2-Tfr_ice) / ( E1 + Ebot );
                // But h2 hasn't been updated, E1 may have changed and Ebot is not in this scope
                // so we just write it out:
                hi -= hi/2*Crho*(T2-Tfr_ice)*T1/( qi*T1 + (Crho*T1-qi)*(Tfr_ice-T1) );
                T2  = Tfr_ice;
            }
        }

        // Book keeping
        del_hi = hi-hi_old;

        /* Make sure we don't get too small hi_new */
        if ( hi < physical::hmin )
        {
            Qio   -= ( -qs*hs + (E1+E2)*hi/2. )/dt; // modified (30) - with multiplication of rhoi and rhos and division with dt

            del_hi = -hi_old;
            hi     = 0.;
            hs     = 0.;
            Tsurf  = Tbot;
            T1     = Tbot;
            T2     = Tbot;
        }
    }
}//thermoWinton

    
//------------------------------------------------------------------------------------------------------
//! Calculates heat fluxes through the ice according to the Semtner zero layer model.
//! Called by the thermo() function.
void
FiniteElement::thermoIce0(int i, double dt, double wspeed, double sphuma, double conc, double voli, double vols, double Qlw_in, double Qsw_in, double mld, double snowfall,
        double &hi, double &hs, double &hi_old, double &Qio, double &del_hi, double &Tsurf,
        double &Qai, double &Qsw, double &Qlw, double &Qsh, double &Qlh)
{

    // Constants
    double const alb_ice = vm["thermo.alb_ice"].as<double>();
    double const alb_sn  = vm["thermo.alb_sn"].as<double>();
    double const I_0     = vm["thermo.I_0"].as<double>();
    int const alb_scheme = vm["thermo.alb_scheme"].as<int>();

    double const drag_ice_t = vm["thermo.drag_ice_t"].as<double>();

    bool const flooding = vm["thermo.flooding"].as<bool>();

    double const qi = physical::Lf * physical::rhoi;
    double const qs = physical::Lf * physical::rhos;

    /* Don't do anything if there's no ice */
    if ( conc <=0. || voli<=0.)
    {
        hi      = 0.;
        hi_old  = 0.;
        hs      = 0.;
        Tsurf   = 0.;
        Qio     = 0.;
        Qai     = 0.;
        del_hi  = 0.;
        Qsw     = 0.;
        Qlw     = 0.;
        Qsh     = 0.;
        Qlh     = 0.;
    } else {
        /* Calculate the slab thickness */
        hi     = voli/conc;
        hi_old = hi;
        hs     = vols/conc;

        /* Local variables */
        double dQaidT, subl;
        double Qic, del_hs, del_ht, del_hb, draft;

        /* ---------------------------------------------------------------
        * Calculate the surface temperature within a while-loop
        * --------------------------------------------------------------- */
        double albedo = this->albedo(alb_scheme, Tsurf, hs, alb_sn, alb_ice, I_0);
        double dtsurf   = 1.;
        double Tbot     = -physical::mu*M_sss[i];
        int nb_iter_while=0;
        while ( dtsurf > 1e-4 )
        {
            nb_iter_while++;

            /* Calculate atmospheric fluxes */
            // Shortwave is modulated by the albedo
            Qsw = -Qsw_in*(1.-albedo)*(1.-I_0);
            // The rest is calculated by bulk formula
            this->atmFluxBulk(i, Tsurf, sphuma, drag_ice_t, Qsw, Qlw_in, wspeed, Qai, dQaidT,subl,
                    Qsh, Qlh, Qlw);

            // -------------------------------------------------
            /* Recalculate Tsurf */
            dtsurf = Tsurf;
            Qic    = physical::ks*( Tbot-Tsurf )/( hs + physical::ks*hi/physical::ki );
            Tsurf = Tsurf + ( Qic - Qai )/
                ( physical::ks/(hs+physical::ks*hi/physical::ki) + dQaidT );

            /* Set Tsurf to the freezing point of snow or ice */
            if ( hs > 0. )
                Tsurf = std::min(0., Tsurf);
            else
                Tsurf = std::min(-physical::mu*physical::si, Tsurf);

            /* Re-evaluate the exit condition */
            dtsurf = std::abs(dtsurf-Tsurf);
        }

        if(nb_iter_while>10)
        {
            LOG(DEBUG) << "nb_iter_while = " << nb_iter_while << "\n";
            //throw std::logic_error("nb_iter_while larger than 10");
        }

        /* Conductive flux through the ice */
        Qic = physical::ks*( Tbot-Tsurf )/( hs + physical::ks*hi/physical::ki );

        /* ---------------------------------------------------------------
         * Melt and growth
         * --------------------------------------------------------------- */

        /* Top melt */
        /* Snow melt and sublimation */
        del_hs = std::min(Qai-Qic,0.)*dt/qs - subl*dt/physical::rhos;
        /* Use the energy left over after snow melts to melt the ice */
        del_ht = std::min(hs+del_hs,0.)*qs/qi;
        /* Can't have negative hs! */
        del_hs = std::max(del_hs,-hs);
        // snowfall in kg/m^2/s
        hs  = hs + del_hs + snowfall/physical::rhos*dt;

        /* Heatflux from ocean */
        Qio = this->iceOceanHeatflux(i, M_sst[i], M_sss[i], mld, dt);
        /* Bottom melt/growth */
        del_hb = (Qic-Qio)*dt/qi;

        /* Combine top and bottom */
        del_hi = del_ht+del_hb;
        hi     = hi + del_hi;

        /* Snow-to-ice conversion */
        draft = ( hi*physical::rhoi + hs*physical::rhos ) / physical::rhow;
        if ( flooding && draft > hi )
        {
            /* Subtract the mass of snow converted to ice from hs_new */
            hs = hs - ( draft - hi )*physical::rhoi/physical::rhos;
            hi = draft;
        }

        /* Make sure we don't get too small hi_new */
        if ( hi < physical::hmin )
        {
            del_hi  = -hi_old; //del_hi-hi;
            Qio     = Qio + hi*qi/dt + hs*qs/dt;

            hi      = 0.;
            hs      = 0.;
            Tsurf   = Tbot;
        }
    }
}//thermoIce0

    
//------------------------------------------------------------------------------------------------------
//! Initializes constants, dataset descriptions, the time, mesh, variables, forcings, bathymetry, moorings and drifters.
//! * Also outputs initial mooring snapshots, drifter outputs, normal outputs and  restarts.
//! Called by the run() function.
void
FiniteElement::init()
{
    //! - 1) Initializes everything that doesn't depend on the mesh (constants, dataset descriptions and time) using the initOptAndParam() function,

    M_comm.barrier();

    pcpt = 0;
    mesh_adapt_step=0;
    had_remeshed=false;

    this->initOptAndParam();
    M_current_time = time_init;

    //! - 2) Initializes the mesh using the initMesh() function,
    this->initMesh();

    if (M_rank==0)
    {
        LOG(INFO) << "-----------------------Simulation started on "<< Nextsim::current_time_local() <<"\n";
        LOG(INFO) <<"TIMESTEP= "<< time_step <<"\n";
        LOG(INFO) <<"DURATION= "<< duration <<"\n";
    }

    // We need to set the scale_coeff et al after initialising the mesh - this was previously done in initConstants
    // The mean resolution of the small_arctic_10km mesh is 7446.71 m. Using 74.5 gives scale_coef = 0.100022, for that mesh
    boost::mpi::broadcast(M_comm, M_res_root_mesh, 0);
    scale_coef = std::sqrt(74.5/M_res_root_mesh);
    C_fix    = cfix*scale_coef;          // C_fix;...  : cohesion (mohr-coulomb) in MPa (40000 Pa)
    C_alea   = alea_factor*C_fix;        // C_alea;... : alea sur la cohesion (Pa)
    LOG(DEBUG) << "SCALE_COEF = " << scale_coef << "\n";

    if ( M_use_restart )
    {
        std::string resfil = vm["restart.filename"].as<std::string>();
        LOG(DEBUG) <<"Reading restart file: "<< resfil <<"\n";
        if ( resfil.empty() )
            throw std::runtime_error("Please provide restart.filename");
        else if ( resfil.length()<=6 )
            throw std::runtime_error("Please provide valid option for restart.filename (currently too short)");
        std::string res_str = resfil.substr(6, resfil.length());
        this->readRestart(res_str);

        //write fields from restart to file (needed?)
        if(M_rank==0)
            LOG(DEBUG) <<"export starts\n";
        this->exportResults("restart", true, true, true);
        if(M_rank==0)
            LOG(DEBUG) <<"export done in " << chrono.elapsed() <<"s\n";
    }
    else
    {
        // Do one regrid to get the mesh right
        //this->regrid(pcpt);

        //! - 3) Initializes variables using the initVariables() function,
        chrono.restart();
        LOG(DEBUG) <<"Initialize variables\n";
        this->initVariables();
    }

    // Check the minimum angle of the grid
    // - needs to be after readRestart, otherwise M_mesh is not initialised yet
    double minang = this->minAngle(M_mesh);
    if (minang < vm["numerics.regrid_angle"].as<double>())
    {
        LOG(INFO) <<"invalid regridding angle: should be smaller than the minimal angle in the initial grid\n";
        throw std::logic_error("invalid regridding angle: should be smaller than the minimal angle in the intial grid");
    }


    //! - 4) Initializes atmospheric and oceanic forcings using the initForcings() function,
    this->initForcings();

    //! - 5) Initializes the bathymetry using the initBathymetry() function,
    LOG(DEBUG) <<"Initialize bathymetry\n";
    this->initBathymetry();

    //! - 6) Loads the data from the datasets initialized in 1) using the checkReloadDatasets(),
    if(M_rank==0)
        LOG(DEBUG) << "init - time-dependant ExternalData objects\n";
    timer["reload"].first.restart();
    this->checkReloadMainDatasets(M_current_time);
    if (M_rank == 0)
        LOG(DEBUG) <<"check_and_reload in "<< timer["reload"].first.elapsed() <<"s\n";

    if ( !M_use_restart )
    {
        timer["state"].first.restart();
        this->initModelState();
        if (M_rank == 0)
            LOG(DEBUG) <<"initModelState done in "<< timer["state"].first.elapsed() <<"s\n";
    }

    if ( M_use_restart && M_use_assimilation )
    {
        timer["assimilation"].first.restart();
        this->DataAssimilation();
        LOG(DEBUG) <<"DataAssimilation done in "<< timer["assimilation"].first.elapsed() <<"s\n";
    }


    //! - 7) Initializes the moorings - if requested - using the initMoorings() function,
    LOG(DEBUG) << "initMoorings\n";
    if ( M_use_moorings )
        this->initMoorings();


    //! - 8) Checks if anything has to be output now using the checkOutputs() function.
    // 1. moorings:
    // - check if we are adding snapshot to netcdf file
    // 2. do we need to init any drifters (also save output at init time)
    // 3. check if writing outputs, and do it if it's time
    // 4. check if writing restart, and do it if it's time
    this->checkOutputs(true);
}//init

    
//------------------------------------------------------------------------------------------------------
//! Increments the model by one time step. Called by the run() function.
//!    * updates drifters,
//!    * remeshes and remaps prognostic variables,
//!    * performs the thermodynamics,
//!    * performs the dynamics.
void
FiniteElement::step()
{
    if (vm["debugging.check_fields"].as<bool>())
        // check fields for nans and if thickness is too big
        this->checkFields();

    //! 1) Remeshes and remaps the prognostic variables

    // The first time step we behave as if we just did a regrid
    M_regrid = (pcpt==0);

    if (vm["numerics.regrid"].as<std::string>() == "bamg")
    {
        double displacement_factor = 1.;
        double minang = this->minAngle(M_mesh,M_UM,displacement_factor);
        LOG(DEBUG) <<"REGRID ANGLE= "<< minang <<"\n";

        if (M_rank == 0)
        {
            if(fmod(pcpt*time_step, ptime_step) == 0)
                std::cout <<"NUMBER OF REGRIDDINGS = " << M_nb_regrid <<"\n";

            std::cout <<"REGRID ANGLE= "<< minang <<"\n";
        }

        if ( minang < vm["numerics.regrid_angle"].as<double>() )
        {
            M_regrid = true;

            if ( M_use_moorings && !M_moorings_snapshot )
                M_moorings.updateGridMean(M_mesh);

            LOG(DEBUG) <<"Regridding starts\n";
            chrono.restart();
            if ( M_use_restart && pcpt==0)
                this->regrid(1); // Special case where the restart conditions imply to remesh
            else
                this->regrid(pcpt);

            LOG(DEBUG) <<"Regridding done in "<< chrono.elapsed() <<"s\n";
            if ( M_use_moorings )
                M_moorings.resetMeshMean(M_mesh);

            ++M_nb_regrid;
        }//M_regrid
    }//bamg-regrid

    M_comm.barrier();

    if ( M_regrid || M_use_restart )
    {
        timer["FETensors"].first.restart();
        this->FETensors();
        if (M_rank == 0)
            std::cout <<"---timer FETensors:              "<< timer["FETensors"].first.elapsed() <<"\n";

        timer["calcCohesion"].first.restart();
        this->calcCohesion();
        if (M_rank == 0)
            std::cout <<"---timer calcCohesion:             "<< timer["calcCohesion"].first.elapsed() <<"\n";

        if (vm["dynamics.use_coriolis"].as<bool>())
        {
            timer["calcCoriolis"].first.restart();
            this->calcCoriolis();
            if (M_rank == 0)
                std::cout <<"---timer calcCoriolis:             "<< timer["calcCoriolis"].first.elapsed() <<"\n";
        }
    }

    if(M_rank==0)
        LOG(DEBUG) << "step - time-dependant ExternalData objects\n";
    timer["reload"].first.restart();
    this->checkReloadMainDatasets(M_current_time+time_step/(24*3600.0));
    if (M_rank == 0)
        std::cout <<"---timer check_and_reload:     "<< timer["reload"].first.elapsed() <<"s\n";

    //======================================================================
    //! 2) Performs the thermodynamics
    //======================================================================
    if ( vm["thermo.use_thermo_forcing"].as<bool>() && (fmod(pcpt*time_step,thermo_timestep) == 0) )
    {
        timer["thermo"].first.restart();
        this->thermo(thermo_timestep);
        if (M_rank == 0)
            std::cout <<"---timer thermo:               "<< timer["thermo"].first.elapsed() <<"s\n";
    }

    //======================================================================
    //! 3) Performs the nesting of the tracers
    //======================================================================
    if( M_use_nesting )
    {
        chrono.restart();
        LOG(DEBUG) <<"nestingIce starts\n";
        this->nestingIce();
        LOG(DEBUG) <<"nestingIce done in "<< chrono.elapsed() <<"s\n";

   //======================================================================
   //! 4) Performs the nesting of the dynamical variables
   //======================================================================

        if( M_nest_dynamic_vars )
        {
            chrono.restart();
            LOG(DEBUG) <<"nestingDynamics starts\n";
            this->nestingDynamics();
            LOG(DEBUG) <<"nestingDynamics done in "<< chrono.elapsed() <<"s\n";
        }
    }

    //======================================================================
    //! 5) Performs the dynamics
    //======================================================================

    if ( M_dynamics_type == setup::DynamicsType::DEFAULT )
    {
        //======================================================================
        //! - 5.1) Assembles the rigidity matrix by calling the assemble() function,
        //======================================================================
        timer["assemble"].first.restart();
        this->assemble(pcpt);
        if (M_rank == 0)
            std::cout <<"---timer assemble:             "<< timer["assemble"].first.elapsed() <<"s\n";

#if 0
        if(had_remeshed && (vm["numerics.regrid_output_flag"].as<bool>()))
        {
            std::string tmp_string3    = (boost::format( "after_assemble_%1%_mesh_adapt_step_%2%" )
                                   % pcpt
                                   % mesh_adapt_step ).str();

            this->exportResults(tmp_string3, true, true, true);

            had_remeshed=false;
        }
#endif

        //======================================================================
        //! - 5.2) Solves the linear problem by calling the solve() function
        //! - 5.3) Updates the velocities by calling the updateVelocity() function
        //! - 5.4) Uptates relevant variables by calling the update() function
        //======================================================================
        timer["solve"].first.restart();
        this->solve();
        if (M_rank == 0)
            std::cout <<"---timer solve:                "<< timer["solve"].first.elapsed() <<"s\n";

        timer["updatevelocity"].first.restart();
        this->updateVelocity();
        if (M_rank == 0)
            std::cout <<"---timer updateVelocity:       "<< timer["updatevelocity"].first.elapsed() <<"s\n";

        timer["update"].first.restart();
        this->update();
        if (M_rank == 0)
            std::cout <<"---timer update:               "<< timer["update"].first.elapsed() <<"s\n";
    }

    if ( M_dynamics_type == setup::DynamicsType::FREE_DRIFT )
    {
        this->updateFreeDriftVelocity();
    }


    //======================================================================
    //! 6) Updates the time
    //======================================================================

    ++pcpt;
    M_current_time = time_init + pcpt*time_step/(24*3600.0);


    //======================================================================
    //! 7) Does the post-processing, checks the output and updates moorings.
    //======================================================================
    
    // 1. moorings:
    // - update fields on grid if outputting mean fields
    // - check if we are adding records to netcdf file
    // 2. check if writing outputs, and do it if it's time
    // 3. check if writing restart, and do it if it's time
    // TODO also add drifter check here
    // - if we move at restart output time we can remove M_UT from
    //   restart files (then it would always be 0)
    this->checkOutputs(false);
 }//step


//------------------------------------------------------------------------------------------------------
//! checks if we are doing any outputs
//! called by <FiniteElement::init>() and <FiniteElement::step>()
void
FiniteElement::checkOutputs(bool const& at_init_time)
{
    //! 1) moorings:
    //! - update fields on grid if outputting mean fields
    //! - check if we are adding records to netcdf file
    //! 2) call checkDrifters
    //!    - do we need to init any, move any, input IABP, output any ... ?
    //! 3) check if writing outputs, and do it if it's time
    //! 4) check if writing restart, and do it if it's time
    // TODO also add drifter check here
    // - if we move at restart output time we can remove M_UT from
    //   restart files (then it would always be 0)

    
    if(M_use_moorings)
    {
        if(!at_init_time)
            this->updateMoorings();
        else if(    M_moorings_snapshot
                && fmod(pcpt*time_step, mooring_output_time_step) == 0 )
        {
            // write initial conditions to moorings file if using snapshot option
            // (only if at the right time though)

            // - set the fields on the mesh
            this->updateMeans(M_moorings, 1.);

            // - interpolate to the grid and write them to the netcdf file
            this->mooringsAppendNetcdf(M_current_time);
        }
    }

    if(M_use_drifters)
    {
        // 1. Gather the fields needed by the drifters
        // 2. Check if we need to:
        //  i.   move any drifters
        //  ii.  input &/or output drifters
        //  iii. init any drifters
        // NB do all this before we output or write restart to make sure
        //  IABP drifters written to restart file are correct
        // TODO also move the drifters if we are writing a restart
        // - then M_UT will always be zero at restart time and can be removed from restart files
        this->checkDrifters();
    }

    // check if we are outputting results file
    if(fmod(pcpt*time_step, output_time_step) == 0)
    {
        chrono.restart();
        LOG(DEBUG) <<"export starts\n";
        this->exportResults(true, true, true);
        LOG(DEBUG) <<"export done in " << chrono.elapsed() <<"s\n";
    }


    // check if writing restart
    if(this->writingRestart())
        this->writeRestart();
}//checkOutputs


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//------------------------------------------------------------------------------------------------------
//! This is the main working function. Called from the main() function. Calls:
//! - init(),
//! - step(),
//! - exportResults(),
//! - finalise().
//------------------------------------------------------------------------------------------------------
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void
FiniteElement::run()
{
    std::string current_time_system = current_time_local();
    // **********************************************************************
    // Initializing
    // **********************************************************************
    this->init();
    int maxiter = vm["debugging.maxiteration"].as<int>();
    int niter = 0;

    // write the logfile: assigned to the process master (rank 0)
    if (M_comm.rank() == 0)
    {
        this->writeLogFile();
    }

    M_current_time = time_init + pcpt*time_step/(24*3600.0);
    bool is_running = true;
    if(duration<0)
        throw std::runtime_error("Set simul.duration >= 0\n");
    else if(duration==0)
        is_running = false;

#if 1
    // main loop for nextsim program
    while (is_running)
    {
         M_comm.barrier();

        if (M_rank == 0)
        {
            std::cout <<"---------------------- TIME STEP "<< pcpt << " : "
                      << model_time_str(vm["simul.time_init"].as<std::string>(), pcpt*time_step);

            if(fmod(pcpt*time_step, ptime_step) == 0)
            {
                std::string time_spent_str = time_spent(current_time_system);
                std::cout <<" ---------- progression: ("<< 100.0*(pcpt*time_step/duration) <<"%) ---------- time spent: "<< time_spent_str <<"\n";
            }

            std::cout <<"\n";
        }

        is_running = ((pcpt+1)*time_step) < duration;

        // **********************************************************************
        // Time-stepping
        // **********************************************************************
        this->step();

        //stop early if debugging
        niter++;
        if (niter == maxiter)
            is_running = false;
    }

    // **********************************************************************
    // Exporting results
    // **********************************************************************
    this->exportResults("final", true, true, true);

    // **********************************************************************
    // Finalizing
    // **********************************************************************
    this->finalise();

#endif

    M_comm.barrier();

    if (M_rank==0)
    {
        std::cout <<"nb regrid total = " << M_nb_regrid <<"\n";

        std::cout << "[INFO]: " << "-----------------------Simulation done on "<< current_time_local() <<"\n";
        std::cout << "[INFO]: " << "-----------------------Total time spent:  "<< time_spent(current_time_system) <<"\n";
    }
}//run


//------------------------------------------------------------------------------------------------------
//! Updates moorings variables (both elemental and nodal) on the grid and multiply by a time factor (time_factor)
//! to weight their contribution to the mean fields.
//! Also updates the position of the grid nodes.
//! Called by the step() function.
void
FiniteElement::updateMeans(GridOutput& means, double time_factor)
{
    // Update elements and multiply with time_factor
    for ( auto it=means.M_elemental_variables.begin(); it!=means.M_elemental_variables.end(); ++it )
    {
        switch (it->varID)
        {
            // Prognostic variables
            case (GridOutput::variableID::conc):
                for (int i=0; i<M_local_nelements; i++)
                {
                    it->data_mesh[i] += M_conc[i]*time_factor;
                    if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
                        it->data_mesh[i] += M_conc_thin[i]*time_factor;
                }
                break;

            case (GridOutput::variableID::thick):
                for (int i=0; i<M_local_nelements; i++)
                {
                    it->data_mesh[i] += M_thick[i]*time_factor;
                    if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
                        it->data_mesh[i] += M_h_thin[i]*time_factor;
                }
                break;

            case (GridOutput::variableID::damage):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += M_damage[i]*time_factor;
                break;

            case (GridOutput::variableID::snow):
                for (int i=0; i<M_local_nelements; i++)
                {
                    it->data_mesh[i] += M_snow_thick[i]*time_factor;
                    if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
                        it->data_mesh[i] += M_hs_thin[i]*time_factor;
                }
                break;

            case (GridOutput::variableID::tsurf):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += ( M_conc[i]*M_tice[0][i] + M_conc_thin[i]*M_tsurf_thin[i] + (1-M_conc[i]-M_conc_thin[i])*M_sst[i] )*time_factor;
                break;

            case (GridOutput::variableID::sst):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += M_sst[i]*time_factor;
                break;

            case (GridOutput::variableID::sss):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += M_sss[i]*time_factor;
                break;

            case (GridOutput::variableID::tsurf_ice):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += M_tice[0][i]*time_factor;
                break;

            case (GridOutput::variableID::t1):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += M_tice[1][i]*time_factor;
                break;

            case (GridOutput::variableID::t2):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += M_tice[2][i]*time_factor;
                break;

            case (GridOutput::variableID::conc_thin):
                for (int i=0; i<M_local_nelements; i++)
                {
                    it->data_mesh[i] += M_conc_thin[i]*time_factor;
                }
                break;

            case (GridOutput::variableID::h_thin):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += M_h_thin[i]*time_factor;
                break;

            case (GridOutput::variableID::hs_thin):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += M_hs_thin[i]*time_factor;
                break;

            // Diagnostic variables
            case (GridOutput::variableID::Qa):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += D_Qa[i]*time_factor;
                break;
            case (GridOutput::variableID::Qsw):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += D_Qsw[i]*time_factor;
                break;
            case (GridOutput::variableID::Qlw):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += D_Qlw[i]*time_factor;
                break;
            case (GridOutput::variableID::Qsh):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += D_Qsh[i]*time_factor;
                break;
            case (GridOutput::variableID::Qlh):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += D_Qlh[i]*time_factor;
                break;
            case (GridOutput::variableID::Qo):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += D_Qo[i]*time_factor;
                break;
            case (GridOutput::variableID::delS):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += D_delS[i]*time_factor;
                break;

            // Non-output variables
            case (GridOutput::variableID::proc_mask):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] = 1;
                break;

            case (GridOutput::variableID::ice_mask):
                for (int i=0; i<M_local_nelements; i++)
                    it->data_mesh[i] += (M_thick[i]>0.) ? 1. : 0.;
                break;

            default: std::logic_error("Updating of given variableID not implimented (elements)");
        }
    }

    // Update nodes
    for ( auto it=means.M_nodal_variables.begin(); it!=means.M_nodal_variables.end(); ++it )
    {
        switch (it->varID)
        {
            case (GridOutput::variableID::VT_x):
                for (int i=0; i<M_num_nodes; i++)
                    it->data_mesh[i] += M_VT[i]*time_factor;
                break;

            case (GridOutput::variableID::VT_y):
                for (int i=0; i<M_num_nodes; i++)
                    it->data_mesh[i] += M_VT[i+M_num_nodes]*time_factor;
                break;

            default: std::logic_error("Updating of given variableID not implimented (nodes)");
        }
    }
}//updateMeans


//------------------------------------------------------------------------------------------------------
//! Initializes the moorings datasets and variables recorded by the moorings.
//! Called by the init() function.
void
FiniteElement::initMoorings()
{

    if (       (!M_moorings_snapshot) 
            && fmod(pcpt*time_step, mooring_output_time_step) != 0 )
    {
        std::string msg = "FE::initMoorings: Start time (or restart time) incompatible with\n";
        msg += "mooring output time step (from moorings.output_timestep option)\n";
        msg += "(when moorings.snapshot = false).";
        throw std::runtime_error(msg);
    }

    bool use_ice_mask = false;

    // Output variables - elements
    std::vector<GridOutput::Variable> elemental_variables(0);

    // Output variables - nodes
    std::vector<GridOutput::Variable> nodal_variables(0);

    // The vectorial variables are (always on the nodes) ...
    std::vector<GridOutput::Vectorial_Variable> vectorial_variables(0);

    // map config file names to the GridOutput enum's
    boost::unordered_map<std::string, GridOutput::variableID>
        mooring_name_map_elements = boost::assign::map_list_of
            ("conc", GridOutput::variableID::conc)
            ("thick", GridOutput::variableID::thick)
            ("snow", GridOutput::variableID::snow)
            ("tsurf", GridOutput::variableID::tsurf)
            ("Qa", GridOutput::variableID::Qa)
            ("Qo", GridOutput::variableID::Qo)
            ("Qsw", GridOutput::variableID::Qsw)
            ("Qlw", GridOutput::variableID::Qlw)
            ("Qsh", GridOutput::variableID::Qsh)
            ("Qlh", GridOutput::variableID::Qlh)
            ("delS", GridOutput::variableID::delS)
            ("conc_thin", GridOutput::variableID::conc_thin)
            ("h_thin", GridOutput::variableID::h_thin)
            ("hs_thin", GridOutput::variableID::hs_thin)
        ;

    std::vector<std::string> names = vm["moorings.variables"].as<std::vector<std::string>>();
    std::vector<std::string> names_thin = {"conc_thin", "h_thin", "hs_thin"};

    for ( auto it=names.begin(); it!=names.end(); ++it )
    {
        // error if trying to output thin ice variables if not using thin ice category
        if (std::count(names_thin.begin(), names_thin.end(), *it) > 0)
        {
            if(M_ice_cat_type!=setup::IceCategoryType::THIN_ICE)
            {
                LOG(ERROR)<<"initMoorings: skipping <<"<< *it<<">> as not running with thin ice\n";
                throw std::runtime_error("Invalid mooring name");
            }
        }

        // Nodal variables and vectors
        if (*it == "velocity")
        {
            use_ice_mask = true; // Needs to be set so that an ice_mask variable is added to elemental_variables below
            GridOutput::Variable siu(GridOutput::variableID::VT_x, use_ice_mask);
            GridOutput::Variable siv(GridOutput::variableID::VT_y, use_ice_mask);
            nodal_variables.push_back(siu);
            nodal_variables.push_back(siv);

            std::vector<int> siuv_id(2);
            siuv_id[0] = 0;
            siuv_id[1] = 1;

            GridOutput::Vectorial_Variable siuv;
            siuv.components_Id = siuv_id;
            vectorial_variables.push_back(siuv);
        }

        // Element variables
        else if (mooring_name_map_elements.count(*it)==0)
        {
            LOG(ERROR)<<"Unimplemented moorings name: "<<*it<<"\n\n";
            LOG(ERROR)<<"Available names are:\n";
            LOG(ERROR)<<"  velocity\n";
            for (auto ptr=mooring_name_map_elements.begin();
                    ptr!=mooring_name_map_elements.end(); ptr++)
                LOG(ERROR)<<"  "<< ptr->first <<"\n";
            throw std::runtime_error("Invalid mooring name");
        }
        else
        {
            GridOutput::Variable tmp(mooring_name_map_elements[*it]);
            elemental_variables.push_back(tmp);
        }
    }

    // Need a mask for the nodal variables
    if ( nodal_variables.size() > 0 )
    {
        GridOutput::Variable proc_mask(GridOutput::variableID::proc_mask);
        elemental_variables.push_back(proc_mask);
    }

    // A mask for velocity (if we want it)
    if ( use_ice_mask )
    {
        GridOutput::Variable ice_mask(GridOutput::variableID::ice_mask);
        elemental_variables.push_back(ice_mask);
    }

    if(vm["moorings.grid_type"].as<std::string>()=="regular")
    {
        // Calculate the grid spacing (assuming a regular grid for now)
        auto RX = M_mesh.coordX();
        auto RY = M_mesh.coordY();
        auto xcoords = std::minmax_element( RX.begin(), RX.end() );
        auto ycoords = std::minmax_element( RY.begin(), RY.end() );

        double xmin = boost::mpi::all_reduce(M_comm, *xcoords.first,  boost::mpi::minimum<double>());
        double xmax = boost::mpi::all_reduce(M_comm, *xcoords.second, boost::mpi::maximum<double>());
        double ymin = boost::mpi::all_reduce(M_comm, *ycoords.first,  boost::mpi::minimum<double>());
        double ymax = boost::mpi::all_reduce(M_comm, *ycoords.second, boost::mpi::maximum<double>());

        double mooring_spacing = 1e3 * vm["moorings.spacing"].as<double>();
        int ncols = (int) ( 0.5 + ( xmax - xmin )/mooring_spacing );
        int nrows = (int) ( 0.5 + ( ymax - ymin )/mooring_spacing );

        // Define the mooring dataset
        M_moorings = GridOutput(M_mesh, ncols, nrows, mooring_spacing, xmin, ymin, nodal_variables,
                elemental_variables, vectorial_variables, M_moorings_averaging_period, M_moorings_false_easting);
    }
    else if(vm["moorings.grid_type"].as<std::string>()=="from_file")
    {
        // Read the grid in from file
        // - add an .mpp file if the projection is different to the neXtSIM projection
        // - if not given, assume it is the same as neXtSIM
        std::string mpp_file = Environment::vm()["moorings.mppfile"].as<std::string>();
        if(!mpp_file.empty())
            mpp_file = (boost::format( "%1%/%2%" )
                    % Environment::nextsimMeshDir().string()
                    % mpp_file
                    ).str();
        else
            mpp_file = Environment::nextsimMppfile();

        GridOutput::Grid grid{
            gridFile: Environment::vm()["moorings.grid_file"].as<std::string>(),
            dirname: "data",
            mpp_file: mpp_file,
            dimNameX: "y",
            dimNameY: "x",
            latName: "latitude",
            lonName: "longitude"
        };

        // Define the mooring dataset
        M_moorings = GridOutput(M_mesh, grid, nodal_variables, elemental_variables, vectorial_variables,
                M_moorings_averaging_period, M_moorings_false_easting);
    }
    else
    {
        throw std::runtime_error("FiniteElement::initMoorings: invalid moorings.grid_type " + vm["moorings.grid_type"].as<std::string>()
                + ". It must be either 'regular' or 'from_file'.");
    }

    // As only the root processor knows the entire grid we set the land mask using it
    if ( M_rank == 0 )
        M_moorings.setLSM(M_mesh_root);

    // Initialise netCDF output
    if ( (M_rank==0) || M_moorings_parallel_output )
    {
        double output_time;
        if ( M_moorings_snapshot )
            output_time = M_current_time;
        else
            // shift the timestamp in the file to the centre of the output interval
            output_time = M_current_time - mooring_output_time_step/86400/2;

        std::string filename_root;
        if ( M_moorings_parallel_output )
            filename_root = M_export_path + "/Moorings_" + std::to_string(M_rank);
        else
            filename_root = M_export_path + "/Moorings";

        M_moorings_file = M_moorings.initNetCDF(filename_root, M_moorings_file_length, output_time);
    }

}//initMoorings


//------------------------------------------------------------------------------------------------------
//! Updates the data recorded by moorings, by calling the updateMeans() function.
//! Called by the step() function.
void
FiniteElement::updateMoorings()
{
    // If we're taking snapshots then we only call updateMeans before writing to file
    // - otherwise we update every time step
    if ( !M_moorings_snapshot )
        this->updateMeans(M_moorings, mooring_time_factor);

    //check if we are outputting
    if ( fmod(pcpt*time_step, mooring_output_time_step) == 0 )
    {
        double output_time = M_current_time;
        if ( M_moorings_snapshot )
        {
            // Update the snapshot
            this->updateMeans(M_moorings, 1.);
        }
        else
        {
            // shift the timestamp in the file to the centre of the output interval
            output_time = M_current_time - mooring_output_time_step/86400/2;
        }

        // If it's a new day we check if we need a new file
        double not_used;
        if (       (M_rank==0 || M_moorings_parallel_output)
                && (M_moorings_file_length != GridOutput::fileLength::inf)
                && (modf(output_time, &not_used) < time_step*86400) )
        {
            std::string filename_root;
            if ( M_moorings_parallel_output )
                filename_root = M_export_path + "/Moorings_" + std::to_string(M_rank);
            else
                filename_root = M_export_path + "/Moorings";

            boost::gregorian::date now = Nextsim::parse_date(output_time);
            switch (M_moorings_file_length)
            {
            case GridOutput::fileLength::daily:
                M_moorings_file = M_moorings.initNetCDF(filename_root, M_moorings_file_length, output_time);
                break;
            case GridOutput::fileLength::weekly:
                if ( now.day_of_week().as_number() == 1 )
                    M_moorings_file = M_moorings.initNetCDF(filename_root, M_moorings_file_length, output_time);
                break;
            case GridOutput::fileLength::monthly:
                if ( now.day().as_number() == 1 )
                    M_moorings_file = M_moorings.initNetCDF(filename_root, M_moorings_file_length, output_time);
                break;
            case GridOutput::fileLength::yearly:
                if ( now.day_of_year() == 1 )
                    M_moorings_file = M_moorings.initNetCDF(filename_root, M_moorings_file_length, output_time);
            }
        }

        // get data on grid and write to netcdf
        // (gathering to master if necessary)
        this->mooringsAppendNetcdf(output_time);

    }//outputting
}//updateMoorings


// -------------------------------------------------------------------------------------
//! update grid means, gather to root processor if requested, and append the fields to the output netcdf file
void
FiniteElement::mooringsAppendNetcdf(double const &output_time)
{
    // update data on grid
    M_moorings.updateGridMean(M_mesh);

    if ( ! M_moorings_parallel_output )
    {
        //gather fields to root processor if not using parallel output
        for (auto it=M_moorings.M_nodal_variables.begin(); it!=M_moorings.M_nodal_variables.end(); ++it)
        {
            std::vector<double> result;
            boost::mpi::reduce(M_comm, it->data_grid, result, std::plus<double>(), 0.);
            if (M_rank==0) it->data_grid = result;
        }
        for (auto it=M_moorings.M_elemental_variables.begin(); it!=M_moorings.M_elemental_variables.end(); ++it)
        {
            std::vector<double> result;
            boost::mpi::reduce(M_comm, it->data_grid, result, std::plus<double>(), 0.);
            if (M_rank==0) it->data_grid = result;
        }
    }

    //append to netcdf
    if ( (M_rank==0) || M_moorings_parallel_output )
        M_moorings.appendNetCDF(M_moorings_file, output_time);

    //reset means on mesh and grid
    M_moorings.resetMeshMean(M_mesh);
    M_moorings.resetGridMean();
}//mooringsAppendNetcdf


//------------------------------------------------------------------------------------------------------
//! do we write a restart file this time step?
//! called by checkOutputs()
bool
FiniteElement::writingRestart()
{
    //check if it's time to write a restart
    if(!vm["restart.write_restart"].as<bool>())
        return false;
    else if(vm["restart.debugging"].as<bool>())
        return true;
    else if ( fmod(pcpt*time_step, restart_time_step) == 0)
        return true;
    else
        return false;
}//writingRestart


//------------------------------------------------------------------------------------------------------
//! Writes restart files.
//! Called by the checkOutputs() function.
void
FiniteElement::writeRestart()
{
    //Determines the name to be passed to writeRestart
    std::string name_str;
    if (vm["output.datetime_in_filename"].as<bool>())
        name_str = to_date_time_string_for_filename(M_current_time);
    else if(vm["restart.debugging"].as<bool>())
        name_str = (boost::format( "%1%" ) % pcpt).str();
    else
    {
        int rstep = pcpt*time_step/restart_time_step;
        name_str = (boost::format( "%1%" ) % rstep).str();
    }
    this->writeRestart(name_str);
}//writeRestart

    
void
FiniteElement::writeRestart(std::string const& name_str)
{
    M_prv_local_ndof = M_local_ndof;
    M_prv_num_nodes = M_num_nodes;
    M_prv_num_elements = M_local_nelements;
    M_prv_global_num_nodes = M_mesh.numGlobalNodes();
    M_prv_global_num_elements = M_mesh.numGlobalElements();

    std::cout<< "["<< M_rank << "] "<<"M_prv_local_ndof         = "<< M_prv_local_ndof <<"\n";
    std::cout<< "["<< M_rank << "] "<<"M_prv_num_nodes          = "<< M_prv_num_nodes <<"\n";
    std::cout<< "["<< M_rank << "] "<< "M_prv_num_elements      = "<< M_prv_num_elements <<"\n";
    std::cout<< "["<< M_rank << "] "<<"M_prv_global_num_nodes   = "<< M_prv_global_num_nodes <<"\n";
    std::cout<< "["<< M_rank << "] "<<"M_prv_global_num_elements= "<< M_prv_global_num_elements <<"\n";
    std::cout<< "["<< M_rank << "] "<<"M_ndof                   = "<< M_ndof <<"\n";

    // fields defined on mesh nodes
    std::vector<double> interp_in_nodes;
    this->gatherFieldsNode(interp_in_nodes, M_rmap_nodes, M_sizes_nodes);

    std::vector<double> M_VT_root;
    std::vector<double> M_VTM_root;
    std::vector<double> M_VTMM_root;
    std::vector<double> M_UM_root;
    std::vector<double> M_UT_root;

    int tmp_nb_var=0;

    if (M_rank == 0)
    {
        M_VT_root.resize(2*M_ndof);
        M_VTM_root.resize(2*M_ndof);
        M_VTMM_root.resize(2*M_ndof);
        M_UM_root.resize(2*M_ndof);
        M_UT_root.resize(2*M_ndof);

        for (int i=0; i<M_ndof; ++i)
        {
            tmp_nb_var = 0;

            // VT_X
            M_VT_root[i] = interp_in_nodes[M_nb_var_node*i+tmp_nb_var];
            tmp_nb_var++;

            // VT_Y
            M_VT_root[i+M_ndof] = interp_in_nodes[M_nb_var_node*i+tmp_nb_var];
            tmp_nb_var++;

            // VTM_X
            M_VTM_root[i] = interp_in_nodes[M_nb_var_node*i+tmp_nb_var];
            tmp_nb_var++;

            // VTM_Y
            M_VTM_root[i+M_ndof] = interp_in_nodes[M_nb_var_node*i+tmp_nb_var];
            tmp_nb_var++;

            // VTMM_X
            M_VTMM_root[i] = interp_in_nodes[M_nb_var_node*i+tmp_nb_var];
            tmp_nb_var++;

            // VTMM_Y
            M_VTMM_root[i+M_ndof] = interp_in_nodes[M_nb_var_node*i+tmp_nb_var];
            tmp_nb_var++;

            // UM_X
            M_UM_root[i] = interp_in_nodes[M_nb_var_node*i+tmp_nb_var];
            tmp_nb_var++;

            // UM_Y
            M_UM_root[i+M_ndof] = interp_in_nodes[M_nb_var_node*i+tmp_nb_var];
            tmp_nb_var++;

            // UT_X
            M_UT_root[i] = interp_in_nodes[M_nb_var_node*i+tmp_nb_var];
            tmp_nb_var++;

            // UT_Y
            M_UT_root[i+M_ndof] = interp_in_nodes[M_nb_var_node*i+tmp_nb_var];
            tmp_nb_var++;

            if(tmp_nb_var>M_nb_var_node)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }
    }

    M_nb_var_element = 15 + M_tice.size();//15;
    int nb_var_element = M_nb_var_element;
    if (M_ice_cat_type!=setup::IceCategoryType::THIN_ICE)
    {
        nb_var_element -= 4;
    }

    std::vector<double> interp_in_elements;
    this->gatherFieldsElementIO(interp_in_elements,M_ice_cat_type==setup::IceCategoryType::THIN_ICE);

    M_comm.barrier();

    if (M_rank == 0)
    {
        int num_elements_root = M_mesh_root.numTriangles();
        int tice_size = M_tice.size();

        std::vector<double> M_conc_root(num_elements_root);
        std::vector<double> M_thick_root(num_elements_root);
        std::vector<double> M_snow_thick_root(num_elements_root);
        std::vector<double> M_sigma_root(3*num_elements_root);
        std::vector<double> M_damage_root(num_elements_root);
        std::vector<double> M_ridge_ratio_root(num_elements_root);
        std::vector<double> M_random_number_root(num_elements_root);
        std::vector<double> M_sss_root(num_elements_root);
        std::vector<double> M_sst_root(num_elements_root);
        std::vector<std::vector<double>> M_tice_root(M_tice.size());
        for(auto it=M_tice_root.begin(); it!=M_tice_root.end(); it++)
            it->resize(num_elements_root);

        std::vector<double> M_h_thin_root;
        std::vector<double> M_conc_thin_root;
        std::vector<double> M_hs_thin_root;
        std::vector<double> M_tsurf_thin_root;

        if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_h_thin_root.resize(num_elements_root);
            M_conc_thin_root.resize(num_elements_root);
            M_hs_thin_root.resize(num_elements_root);
            M_tsurf_thin_root.resize(num_elements_root);
        }

        for (int i=0; i<M_mesh_root.numTriangles(); ++i)
        {
            tmp_nb_var=0;
            int ri = M_rmap_elements[i];

            // concentration
            M_conc_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // thickness
            M_thick_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // snow thickness
            M_snow_thick_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress1
            M_sigma_root[3*i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress2
            M_sigma_root[3*i+1] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress3
            M_sigma_root[3*i+2] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // damage
            M_damage_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // ridge_ratio
            M_ridge_ratio_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // ridge_ratio
            M_random_number_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // sss
            M_sss_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // sst
            M_sst_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // M_tice
            for(auto it=M_tice_root.begin(); it!=M_tice_root.end(); it++)
            {
                (*it)[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;
            }

            if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                // h_thin
                M_h_thin_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;

                // conc_thin
                M_conc_thin_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;

                // hs_thin
                M_hs_thin_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;

                // tsurf_thin
                M_tsurf_thin_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;
            }

            if(tmp_nb_var>nb_var_element)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }


        Exporter exporter("double");
        std::string filename;

        // === Start with the mesh ===
        // First the data
        // NB directory is never empty, due to the default of output.exporter_path
        std::string directory = vm["output.exporter_path"].as<std::string>() + "/restart";

        // create the output directory if it does not exist
        fs::path output_path(directory);
        if ( !fs::exists(output_path) )
            fs::create_directories(output_path);

        filename = (boost::format( "%1%/mesh_%2%.bin" )
                    % directory
                    % name_str ).str();

        std::fstream meshbin(filename, std::ios::binary | std::ios::out | std::ios::trunc);
        if ( ! meshbin.good() )
            throw std::runtime_error("Cannot write to file: " + filename);
        exporter.writeMesh(meshbin, M_mesh_root);
        meshbin.close();

        // Then the record
        filename = (boost::format( "%1%/mesh_%2%.dat" )
                    % directory
                    % name_str ).str();

        std::fstream meshrecord(filename, std::ios::out | std::ios::trunc);
        if ( ! meshrecord.good() )
            throw std::runtime_error("Cannot write to file: " + filename);
        exporter.writeRecord(meshrecord,"mesh");
        meshrecord.close();

        // === Write the prognostic variables ===
        // First the data
        filename = (boost::format( "%1%/field_%2%.bin" )
                    % directory
                    % name_str ).str();
        std::fstream outbin(filename, std::ios::binary | std::ios::out | std::ios::trunc );
        if ( ! outbin.good() )
            throw std::runtime_error("Cannot write to file: " + filename);

        std::vector<int> misc_int(4);
        misc_int[0] = pcpt;
        misc_int[1] = M_flag_fix;
        misc_int[2] = mesh_adapt_step;
        misc_int[3] = M_nb_regrid;

        exporter.writeField(outbin, misc_int, "Misc_int");
        exporter.writeField(outbin, M_dirichlet_flags_root, "M_dirichlet_flags");

        std::vector<double> timevec(1);
        timevec[0] = M_current_time;
        exporter.writeField(outbin, timevec, "Time");
        exporter.writeField(outbin, M_conc_root, "M_conc");
        exporter.writeField(outbin, M_thick_root, "M_thick");
        exporter.writeField(outbin, M_snow_thick_root, "M_snow_thick");
        exporter.writeField(outbin, M_sigma_root, "M_sigma");
        exporter.writeField(outbin, M_damage_root, "M_damage");
        exporter.writeField(outbin, M_ridge_ratio_root, "M_ridge_ratio");
        exporter.writeField(outbin, M_random_number_root, "M_random_number");
        exporter.writeField(outbin, M_sst_root, "M_sst");
        exporter.writeField(outbin, M_sss_root, "M_sss");
        for (int i=0; i<M_tice.size(); i++)
            exporter.writeField(outbin, M_tice_root[i],
                    "M_tice_"+std::to_string(i));

        exporter.writeField(outbin, M_VT_root, "M_VT");
        exporter.writeField(outbin, M_VTM_root, "M_VTM");
        exporter.writeField(outbin, M_VTMM_root, "M_VTMM");
        exporter.writeField(outbin, M_UM_root, "M_UM");
        exporter.writeField(outbin, M_UT_root, "M_UT");

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            exporter.writeField(outbin, M_h_thin_root, "M_h_thin");
            exporter.writeField(outbin, M_conc_thin_root, "M_conc_thin");
            exporter.writeField(outbin, M_hs_thin_root, "M_hs_thin");
            exporter.writeField(outbin, M_tsurf_thin_root, "M_tsurf_thin");
        }

        if (       M_use_iabp_drifters
                && M_iabp_drifters.size()>0)
        {
            // if drifters not initialised yet, don't try to write them
            std::vector<int> drifter_no(M_iabp_drifters.size());
            std::vector<double> drifter_x(M_iabp_drifters.size());
            std::vector<double> drifter_y(M_iabp_drifters.size());

            // Sort the drifters so the restart files are identical (this is just to make testing easier)
            int j=0;
            for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); ++it )
            {
                drifter_no[j] = it->first;
                ++j;
            }
            std::sort(drifter_no.begin(), drifter_no.end());

            j=0;
            for ( auto it = drifter_no.begin(); it != drifter_no.end(); it++ )
            {
                drifter_x[j] = M_iabp_drifters[drifter_no[j]][0];
                drifter_y[j] = M_iabp_drifters[drifter_no[j]][1];
                j++;
            }

            exporter.writeField(outbin, drifter_no, "Drifter_no");
            exporter.writeField(outbin, drifter_x, "Drifter_x");
            exporter.writeField(outbin, drifter_y, "Drifter_y");
        }
        
        // Finally add the previous numbering to the restart file
        // used in adaptMesh (updateNodeIds)
        std::vector<double> PreviousNumbering(M_mesh_root.numNodes());
        for ( int i=0; i<M_mesh_root.numNodes(); ++i )
            PreviousNumbering[i] = bamgmesh_root->PreviousNumbering[i];
        exporter.writeField(outbin, PreviousNumbering, "PreviousNumbering");

        outbin.close();

        // Then the record
        filename = (boost::format( "%1%/field_%2%.dat" )
                    % directory
                    % name_str ).str();

        std::fstream outrecord(filename, std::ios::out | std::ios::trunc);
        if ( ! outrecord.good() )
            throw std::runtime_error("Cannot write to file: " + filename);

        exporter.writeRecord(outrecord);
        outrecord.close();
    }
}//writeRestart

    
//------------------------------------------------------------------------------------------------------
//! Reads restart files. Called by the init() function.
void
FiniteElement::readRestart(std::string const& name_str)
{
    Exporter exp_field("double"), exp_mesh("double");
    std::string filename;
    boost::unordered_map<std::string, std::vector<int>>    field_map_int;
    boost::unordered_map<std::string, std::vector<double>> field_map_dbl;
    std::vector<int> misc_int;
    std::vector<double> time_vec;

    if (M_rank == 0)
    {

        //! - Reads in the mesh restart files,
        std::string restart_path = vm["restart.input_path"].as<std::string>();
        if ( restart_path.empty() )
            throw std::runtime_error("need to define restart.input option if starting from restart");

        //! - Starts with the record,
        filename = (boost::format( "%1%/mesh_%2%.dat" )
                    % restart_path
                    % name_str ).str();
        LOG(DEBUG)<<"restart file = "<<filename<<"\n";

        std::ifstream meshrecord(filename);
        if ( ! meshrecord.good() )
            throw std::runtime_error("File not found: " + filename);

        exp_mesh.readRecord(meshrecord);
        meshrecord.close();

        //! - Reads in the data itself
        filename = (boost::format( "%1%/mesh_%2%.bin" )
                    % restart_path
                    % name_str ).str();

        std::fstream meshbin(filename, std::ios::binary | std::ios::in );
        if ( ! meshbin.good() )
            throw std::runtime_error("File not found: " + filename);

        exp_mesh.loadFile(meshbin, field_map_int, field_map_dbl);
        meshbin.close();

        std::vector<int>   indexTr = field_map_int["Elements"];
        std::vector<double> coordX = field_map_dbl["Nodes_x"];
        std::vector<double> coordY = field_map_dbl["Nodes_y"];
        std::vector<int>   nodeId = field_map_int["id"];

        // === Read in the prognostic variables ===
        // Start with the record
        filename = (boost::format( "%1%/field_%2%.dat" )
                    % restart_path
                    % name_str ).str();

        std::ifstream inrecord(filename);
        if ( ! inrecord.good() )
            throw std::runtime_error("File not found: " + filename);

        exp_field.readRecord(inrecord);
        inrecord.close();

        // Then onto the data itself
        filename = (boost::format( "%1%/field_%2%.bin" )
                    % restart_path
                    % name_str ).str();

        std::fstream inbin(filename, std::ios::binary | std::ios::in);

        if ( ! inbin.good() )
            throw std::runtime_error("File not found: " + filename);

        field_map_int.clear();
        field_map_dbl.clear();
        exp_field.loadFile(inbin, field_map_int, field_map_dbl);
        inbin.close();

        //! - Recreates the mesh grid
        // Create bamgmesh and bamggeom
        BamgConvertMeshx(
                         bamgmesh_root,bamggeom_root,
                         &indexTr[0],&coordX[0],&coordY[0],
                         coordX.size(), indexTr.size()/3.
                         );

        // read time and misc_int
        time_vec = field_map_dbl["Time"];
        misc_int = field_map_int["Misc_int"];

        // Fix boundaries
        M_flag_fix = misc_int[1];
        std::vector<int> dirichlet_flags = field_map_int["M_dirichlet_flags"];

        for (int edg=0; edg<bamgmesh_root->EdgesSize[0]; ++edg)
        {
            int fnd = bamgmesh_root->Edges[3*edg];

            if ((std::binary_search(dirichlet_flags.begin(),dirichlet_flags.end(),fnd)))
            {
                bamggeom_root->Edges[3*edg+2] = M_flag_fix;
                bamgmesh_root->Edges[3*edg+2] = M_flag_fix;
            }
            else
            {
                bamggeom_root->Edges[3*edg+2] = M_flag_fix+1; // we just want it to be different than M_flag_fix
                bamgmesh_root->Edges[3*edg+2] = M_flag_fix+1; // we just want it to be different than M_flag_fix
            }
        }

        //! - Imports the bamg structs
        this->importBamg(bamgmesh_root);
        this->updateBoundaryFlags();// update boundary flags
        M_mesh_root.setId(nodeId);  // set the node id's

        //! - Adds the previous numbering from the restart file used in adaptMesh (updateNodeIds)
        std::vector<double> PreviousNumbering = field_map_dbl["PreviousNumbering"];
        for ( int i=0; i<M_mesh_root.numNodes(); ++i )
            bamgmesh_root->PreviousNumbering[i] = PreviousNumbering[i];
    }//M_rank==0

    // mesh partitioning
    this->partitionMeshRestart();

    //set time and counters
    mesh_adapt_step = 0;
    M_nb_regrid = 0;
    if(M_rank==0)
    {
        // Set and check time
        if (!vm["restart.reset_time_counter"].as<bool>())
        {
            pcpt = misc_int[0];
            double tmp = time_init + pcpt*time_step/(24*3600.0);
            if ( time_vec[0] != tmp )
            {
                std::cout << "FiniteElement::readRestart: Time and Misc_int[0] (a.k.a pcpt) are inconsistent. \n";
                std::cout << "Time = " << time_vec[0] << " = " << to_date_time_string(time_vec[0])<<"\n";
                std::cout << "time_init + pcpt*time_step/(24*3600.0) = " << tmp << " = " << to_date_time_string(tmp)<<"\n";
                throw std::runtime_error("Inconsistent time information in restart file");
            }

            //set other counters from the restart file
            mesh_adapt_step = misc_int[2];
            M_nb_regrid     = misc_int[3];
	    }
	    else
	    {
            if ( time_vec[0] != time_init )
            {
                std::cout << "FiniteElement::readRestart: Restart Time and time_init are inconsistent. \n";
                std::cout << "Time = " << time_vec[0] << " = " << to_date_time_string(time_vec[0])<<"\n";
                std::cout << "time_init = " << time_init << " = " << to_date_time_string(time_init) <<"\n";
                throw std::runtime_error("Inconsistent time information in restart file");
            }
	    }
    }

    //transfer scalars from "Misc_int" from root to all processors
    boost::mpi::broadcast(M_comm, pcpt, 0);
    boost::mpi::broadcast(M_comm, M_flag_fix, 0);
    boost::mpi::broadcast(M_comm, mesh_adapt_step, 0);
    boost::mpi::broadcast(M_comm, M_nb_regrid, 0);
    M_current_time = time_init + pcpt*time_step/(24*3600.0);
    if(M_use_drifters)
    {
        // if current time is ahead of init-drifter time
        // (as it could be in a restart situation),
        // increase the init-drifter time
        M_drifters_time_init = std::max(M_drifters_time_init,
                std::ceil(M_current_time));
    }

    // set all variables to 0
    // - best to do this early on so M_tice has the right number of components
    //   as early as possible
    // Restart variables are reset on root to be the values from the restart file,
    // then resized back later on
    this->initVariables();

    if (M_rank == 0)
    {
        int num_elements_root = M_mesh_root.numTriangles();

        for (int i=0; i<M_tice.size(); ++i)
            LOG(DEBUG)<<"size M_tice["<<i<<"]= "<< (M_tice[i]).size() <<"\n";

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_h_thin.resize(num_elements_root);
            M_conc_thin.resize(num_elements_root);
            M_hs_thin.resize(num_elements_root);
            M_tsurf_thin.resize(num_elements_root);
        }

        M_conc          = field_map_dbl["M_conc"];
        M_thick         = field_map_dbl["M_thick"];
        M_snow_thick    = field_map_dbl["M_snow_thick"];
        M_sigma         = field_map_dbl["M_sigma"];
        M_damage        = field_map_dbl["M_damage"];
        M_ridge_ratio   = field_map_dbl["M_ridge_ratio"];
        M_random_number = field_map_dbl["M_random_number"];
        M_sst           = field_map_dbl["M_sst"];
        M_sss           = field_map_dbl["M_sss"];
        for (int i=0; i<M_tice.size(); i++)
            M_tice[i] = field_map_dbl["M_tice_"+std::to_string(i)];

        // Pre-processing
        M_VT   = field_map_dbl["M_VT"];
        M_VTM  = field_map_dbl["M_VTM"];
        M_VTMM = field_map_dbl["M_VTMM"];
        M_UM   = field_map_dbl["M_UM"];
        M_UT   = field_map_dbl["M_UT"];
        if(vm["restart.restart_at_rest"].as<bool>())
        {
            // reset M_sigma, M_VT[,M,MM] = 0
            // NB don't reset M_UT = 0 (for drifters)
            // TODO should M_UM = 0 ? - this is the mesh displacement (not part of the rheology)
            for (int i=0; i < M_sigma.size(); i++)
                M_sigma[i] = 0.;

            for (int i=0; i < M_VT.size(); i++)
            {
                M_VT[i]   = 0.;
                M_VTM[i]  = 0.;
                M_VTMM[i] = 0.;
                M_UM[i]   = 0.;
            }
        }

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_h_thin     = field_map_dbl["M_h_thin"];
            M_conc_thin  = field_map_dbl["M_conc_thin"];
            M_hs_thin    = field_map_dbl["M_hs_thin"];
            M_tsurf_thin = field_map_dbl["M_tsurf_thin"];
        }

        if (M_use_iabp_drifters)
        {
            std::vector<int>    drifter_no = field_map_int["Drifter_no"];
            std::vector<double> drifter_x  = field_map_dbl["Drifter_x"];
            std::vector<double> drifter_y  = field_map_dbl["Drifter_y"];

            if (drifter_no.size() == 0)
            {
                // Do nothing now - wait till checkDrifters() to init the IABP drifters in the normal way
                if(M_rank==0)
                {
                    LOG(WARNING) << "Warning: Couldn't read drifter positions from restart file."
                        << " Drifter positions initialised as if there was no restart.\n";
                }
            }
            else
            {
                // init the input/output files
                this->initIabpDrifterFiles();

                // insert the drifter positions from the restart file
                for ( int i=0; i<drifter_no.size(); ++i )
                    M_iabp_drifters.emplace(
                            drifter_no[i], std::array<double,2>{drifter_x[i], drifter_y[i]});

                // move the drifters
                // NB updateIabpDrifterPosition uses M_UT_root 
                // (don't need to gather since we have it already from
                // restart file)
                M_UT_root = M_UT;
                this->updateIabpDrifterPosition();

                // M_UT, M_UT_root can be set to zero now we have moved the drifters
                std::fill(M_UT.begin(), M_UT.end(), 0.);
                std::fill(M_UT_root.begin(), M_UT_root.end(), 0.);

                // still need to get M_conc at these positions
                M_conc_root = M_conc;
                M_UM_root = M_UM;
                auto movedmesh = M_mesh_root;
                movedmesh.move(M_UM_root, 1.);
                this->updateIabpDrifterConc(movedmesh);

                // Save the initial positions to the output file
                this->outputIabpDrifters();
            }
        }//M_use_iabp_drifters
    }//M_rank==0


    // get names of the variables in the restart file,
    // and set data (pointers to the corresponding vectors)
    // and num_components (number of components in those vectors)
    std::vector<std::vector<double>*> data_elements;
    std::vector<int> num_components_elements;
    std::vector<std::string> names = this->getRestartVariableNames();
    this->getVariablesIO(data_elements, num_components_elements, names);

    // get the variables (only on the root processor so far)
    // from data and put it not interp_elt_out
    // TODO do something similar for the nodes
    std::vector<double> interp_elt_out;
    std::vector<double> interp_nd_out;
    this->collectRootRestart(interp_elt_out, interp_nd_out,
            data_elements, num_components_elements);

    // Scatter elemental fields from root and put it in data_elements
    // inside a loop
    // - data_elements is a vector of pointers so the required
    //  variables are now set
    this->scatterFieldsElementIO(interp_elt_out,
            data_elements, num_components_elements);

    // Scatter nodal fields from root
    this->scatterFieldsNode(&interp_nd_out[0]);
}//readRestart
    
    
//------------------------------------------------------------------------------------------------------
//! Partitions the mesh during a restart.
//! Called by the readRestart() function.
void
FiniteElement::partitionMeshRestart()
{
    M_comm.barrier();

    if (M_rank == 0)
    {
        std::cout<<"------------------------------version       = "<< M_mesh_root.version() <<"\n";
        std::cout<<"------------------------------ordering      = "<< M_mesh_root.ordering() <<"\n";
        std::cout<<"------------------------------format        = "<< M_mesh_fileformat <<"\n";
        std::cout<<"------------------------------space         = "<< vm["mesh.partitioner-space"].as<std::string>() <<"\n";
        std::cout<<"------------------------------partitioner   = "<< vm["mesh.partitioner"].as<std::string>() <<"\n";

        // Environment::logMemoryUsage("before partitioning...");
        timer["savemesh"].first.restart();
        LOG(DEBUG) <<"Saving mesh starts\n";
        if (M_partition_space == mesh::PartitionSpace::MEMORY)
            M_mesh_root.writeToGModel();
        else if (M_partition_space == mesh::PartitionSpace::DISK)
            M_mesh_root.writeToFile(M_partitioned_mesh_filename);

        std::cout <<"Saving mesh done in "<< timer["savemesh"].first.elapsed() <<"s\n";

        // partition the mesh on root process (rank 0)
        timer["meshpartition"].first.restart();
        LOG(DEBUG) <<"Partitioning mesh starts\n";
        M_mesh_root.partition(M_partitioned_mesh_filename,
                M_partitioner, M_partition_space, M_mesh_fileformat);
        std::cout <<"Partitioning mesh done in "<< timer["meshpartition"].first.elapsed() <<"s\n";
    }

    M_prv_local_ndof = M_local_ndof;
    M_prv_num_nodes = M_num_nodes;
    M_prv_num_elements = M_local_nelements;
    M_prv_global_num_nodes = M_mesh.numGlobalNodes();
    M_prv_global_num_elements = M_mesh.numGlobalElements();

    this->distributedMeshProcessing(true);
}//partitionMeshRestart

    
//------------------------------------------------------------------------------------------------------
//! Gets the variables (only on the root processor so far) from data and store it in a structure (interp_elt_out)
//! Called by the readRestart() function.
void
FiniteElement::collectRootRestart(std::vector<double>& interp_elt_out,
        std::vector<double>& interp_nd_out,
        std::vector<std::vector<double>*> &data_elements,
        std::vector<int> &num_components_elements)
{
    // * output: interp_elt_out is vector containing all the variables
    //   on the elements to be scattered from root during readRestart
    // * output: interp_nd_out is vector containing all the variables
    //   on the nodes to be scattered from root during readRestart
    // * data_elements is a vector of pointers to the elemental variables to go
    //   into interp_elt_out
    // * num_components_elements is a vector with the number of components in
    //   each elemental variable (usually 1, but can be 3 eg for M_sigma)
    int const nb_var_element = std::accumulate(
            num_components_elements.begin(), num_components_elements.end(), 0);

    if (M_rank == 0)
    {
        int num_elements_root = M_mesh_root.numTriangles();
        int tice_size = M_tice.size();

        interp_elt_out.resize(nb_var_element*num_elements_root);
        for (int i=0; i<num_elements_root; ++i)
        {
            int tmp_nb_var=0;
            for(int j=0; j<data_elements.size(); j++)
            {
                int num_comp = num_components_elements[j];
                for (int k=0; k<num_comp; k++)
                {
                    interp_elt_out[nb_var_element*i+tmp_nb_var]
                        = (*(data_elements[j]))[num_comp*i+k];
                    tmp_nb_var++;
                }//loop over each component of variables
            }//loop over variables
            if(tmp_nb_var!=nb_var_element)
                throw std::logic_error("tmp_nb_var not equal to nb_var_element");
        }
    }//M_rank == 0: collect elemental variables

    M_nb_var_node = 10;
    if (M_rank == 0)
    {
        int num_nodes_root = M_mesh_root.numNodes();
        interp_nd_out.resize(M_nb_var_node*num_nodes_root,0.);

        int tmp_nb_var = 0;

        for (int i=0; i<num_nodes_root; ++i)
        {
            tmp_nb_var = 0;

            // VT
            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_VT[i];
            tmp_nb_var++;

            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_VT[i+num_nodes_root];
            tmp_nb_var++;

            // VTM
            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_VTM[i];
            tmp_nb_var++;

            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_VTM[i+num_nodes_root];
            tmp_nb_var++;

            // VTMM
            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_VTMM[i];
            tmp_nb_var++;

            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_VTMM[i+num_nodes_root];
            tmp_nb_var++;

            // UM
            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_UM[i];
            tmp_nb_var++;

            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_UM[i+num_nodes_root];
            tmp_nb_var++;

            // UT
            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_UT[i];
            tmp_nb_var++;

            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_UT[i+num_nodes_root];
            tmp_nb_var++;

            if(tmp_nb_var>M_nb_var_node)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }
    }
}//collectRootRestart
    

//------------------------------------------------------------------------------------------------------
//! Updates the ice velocity and total displacement (that is used for the drifters).
//! Called by the step() function.
void
FiniteElement::updateVelocity()
{
    M_VTMM = M_VTM;
    M_VTM = M_VT;
    M_VT = M_solution->container();

    // increment M_UT that is used for the drifters
    for (int nd=0; nd<M_UT.size(); ++nd)
    {
        M_UT[nd] += time_step*M_VT[nd]; // Total displacement (for drifters)
    }
}//updateVelocity


//------------------------------------------------------------------------------------------------------
//! Calculates-updates the free drift velocity (no rheology term), if option DynamicsType is set to FREE_DRIFT.
//! Called by the step() function.
void
FiniteElement::updateFreeDriftVelocity()
{
    M_VTMM = M_VTM;
    M_VTM  = M_VT;

    int index_u, index_v;
    double norm_Voce_ice, coef_Voce, norm_Vair_ice, coef_Vair;

    double norm_Voce_ice_min= 0.01; // minimum value to avoid 0 water drag term.
    double norm_Vair_ice_min= 0.01; // minimum value to avoid 0 water drag term.

    for (int nd=0; nd<M_num_nodes; ++nd)
    {
        if(M_mask_dirichlet[nd]==false)
        {
            index_u = nd;
            index_v = nd+M_num_nodes;

            // Free drift case
            norm_Voce_ice = std::hypot(M_VT[index_u]-M_ocean[index_u],M_VT[index_v]-M_ocean[index_v]);
            norm_Voce_ice = (norm_Voce_ice > norm_Voce_ice_min) ? (norm_Voce_ice):norm_Voce_ice_min;

            coef_Voce = (vm["dynamics.lin_drag_coef_water"].as<double>()+(quad_drag_coef_water*norm_Voce_ice));
            coef_Voce *= physical::rhow;

            norm_Vair_ice = std::hypot(M_VT[index_u]-M_wind [index_u],M_VT[index_v]-M_wind [index_v]);
            norm_Vair_ice = (norm_Vair_ice > norm_Vair_ice_min) ? (norm_Vair_ice):norm_Vair_ice_min;

            coef_Vair = (vm["dynamics.lin_drag_coef_air"].as<double>()+(quad_drag_coef_air*norm_Vair_ice));
            coef_Vair *= (physical::rhoa);

            M_VT[index_u] = ( coef_Vair*M_wind [index_u] + coef_Voce*M_ocean [index_u] ) / ( coef_Vair+coef_Voce );
            M_VT[index_v] = ( coef_Vair*M_wind [index_v] + coef_Voce*M_ocean [index_v] ) / ( coef_Vair+coef_Voce );

            // increment M_UT that is used for the drifters
            M_UT[index_u] += time_step*M_VT[index_u]; // Total displacement (for drifters)
            M_UT[index_v] += time_step*M_VT[index_v]; // Total displacement (for drifters)
        }
    }
}//updateFreeDriftVelocity

    
//------------------------------------------------------------------------------------------------------
//! ??
//! !Does not seem to be used!
#if 0
void
FiniteElement::speedScaling(std::vector<double>& speed_scaling)
{
    std::vector<int> sizes_elements = M_sizes_elements;
    std::vector<double> conc_local(M_local_nelements);

    for (int i=0; i<M_local_nelements; ++i)
    {
        // concentration
        conc_local[i] = M_conc[i];
    }

    std::vector<double> conc_root;

    if (M_rank == 0)
    {
        conc_root.resize(M_mesh_root.numTriangles());
        boost::mpi::gatherv(M_comm, conc_local, &conc_root[0], sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, conc_local, 0);
    }

    std::vector<double> speed_scaling_vec;

    if (M_rank == 0)
    {
        auto conc_root_nrd = conc_root;
        int gsize = M_mesh_root.numTriangles();

        for (int i=0; i<gsize; ++i)
        {
            //int ri = rmap_elements.left.find(i+1)->second-1;
            int ri = M_rmap_elements[i];
            // concentration
            conc_root[i] = conc_root_nrd[ri];
        }

        speed_scaling_vec.resize(bamgmesh_root->NodalElementConnectivitySize[0]);

        int elt_num, i, j;
        double c_max_nodal_neighbour;
        double speed_c_scaling;

        std::vector<double> cloc_elts(bamgmesh_root->NodalElementConnectivitySize[1]);

        for (i=0; i<bamgmesh_root->NodalElementConnectivitySize[0]; ++i)
        {
            for (j=0; j<bamgmesh_root->NodalElementConnectivitySize[1]; ++j)
            {
                elt_num = bamgmesh_root->NodalElementConnectivity[bamgmesh_root->NodalElementConnectivitySize[1]*i+j]-1;

                if ((0 <= elt_num) && (elt_num < gsize) && (elt_num != NAN))
                {
                    cloc_elts[j] = conc_root[elt_num];
                    //cloc_elts.push_back( conc_root[elt_num] );
                }
                else
                {
                    break;
                }
            }

            c_max_nodal_neighbour = *std::max_element(cloc_elts.begin(),cloc_elts.begin()+j-1);
            //c_max_nodal_neighbour = *std::max_element(cloc_elts.begin(),cloc_elts.end());

            c_max_nodal_neighbour = c_max_nodal_neighbour/vm["simul.drift_limit_concentration"].as<double>();
            speed_c_scaling = std::min(1.,c_max_nodal_neighbour);
            //std::cout<<"c_max_nodal_neighbour["<< i <<"]= "<< c_max_nodal_neighbour <<"\n";
            //std::cout<<"speed_c_scaling["<< i <<"]= "<< speed_c_scaling <<"\n";
            speed_scaling_vec[i] = speed_c_scaling;
        }

        auto speed_scaling_vec_nrd = speed_scaling_vec;
        speed_scaling_vec.resize(M_id_nodes.size());

        for (int i=0; i<M_id_nodes.size(); ++i)
        {
            int ri = M_id_nodes[i]-1;
            speed_scaling_vec[i] = speed_scaling_vec_nrd[ri];
        }
    }

    speed_scaling.resize(M_num_nodes);
    std::vector<int> sizes_nodes = M_sizes_nodes_with_ghost;

    if (M_rank == 0)
    {
        boost::mpi::scatterv(M_comm, speed_scaling_vec, sizes_nodes, &speed_scaling[0], 0);
    }
    else
    {
        boost::mpi::scatterv(M_comm, &speed_scaling[0], M_num_nodes, 0);
    }
}//speedScaling
#endif


//------------------------------------------------------------------------------------------------------
//! Sets the physical variables relevant to the atmosphere according to the chosen atmospheric forcing data (CONSTANT, ASR, ERAi, ...)
//! Called by the initForcings() function.
void
FiniteElement::forcingAtmosphere()
{
    double air_temperature_correction=vm["forecast.air_temperature_correction"].as<double>();

    switch (M_atmosphere_type)
    {
        case setup::AtmosphereType::CONSTANT:
            M_wind=ExternalData(
                vm["ideal_simul.constant_wind_u"].as<double>(),
                vm["ideal_simul.constant_wind_v"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());

            M_tair=ExternalData(vm["ideal_simul.constant_tair"].as<double>());
            M_mixrat=ExternalData(vm["ideal_simul.constant_mixrat"].as<double>());
            M_mslp=ExternalData(vm["ideal_simul.constant_mslp"].as<double>());
            M_Qsw_in=ExternalData(vm["ideal_simul.constant_Qsw_in"].as<double>());
            M_Qlw_in=ExternalData(vm["ideal_simul.constant_Qlw_in"].as<double>());
            M_snowfr=ExternalData(vm["ideal_simul.constant_snowfr"].as<double>());
            M_precip=ExternalData(vm["ideal_simul.constant_precip"].as<double>());
            M_dair=ExternalData(vm["ideal_simul.constant_dair"].as<double>());
        break;

        case setup::AtmosphereType::ASR:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_mixrat=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            M_snowfr=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
        break;

        case setup::AtmosphereType::ERAi:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0,true ,
                time_init, vm["simul.spinup_duration"].as<double>());

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_snowfall=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
        break;

        case setup::AtmosphereType::EC:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,air_temperature_correction,false,time_init);
            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,air_temperature_correction,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);

            // Syl: The following two lines should be removed when approxSW will be implemented in Thermo()
            M_Qsw_in=ExternalData(vm["ideal_simul.constant_Qsw_in"].as<double>());
            M_precip=ExternalData(0.);
        break;

        case setup::AtmosphereType::EC2:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,air_temperature_correction,false,time_init);
            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,air_temperature_correction,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
        break;

        case setup::AtmosphereType::EC_ERAi:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_bis_elements_dataset,M_mesh,3,false,time_init);
            M_precip=ExternalData(&M_atmosphere_bis_elements_dataset,M_mesh,5,false,time_init);
        break;

        case setup::AtmosphereType::CFSR_HI:
        case setup::AtmosphereType::CFSR:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_sphuma=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_snowfr=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
        break;

        default:
            std::cout << "invalid wind forcing"<<"\n";
            throw std::logic_error("invalid wind forcing");
    }

    // add the external data objects to M_external_data_nodes or M_external_data_elements
    // for looping
    // - these are (or should be) common to all the forcings
    // TODO raise error if not?
    M_external_data_nodes.push_back(&M_wind);
    M_external_data_elements.push_back(&M_tair);
    M_external_data_elements.push_back(&M_mslp);
    M_external_data_elements.push_back(&M_precip);
    M_external_data_elements.push_back(&M_Qsw_in);

    // either need the long wave input, or cloud cover to parameterise it
    if(M_Qlw_in.isInitialized())
        M_external_data_elements.push_back(&M_Qlw_in);
    else if(M_tcc.isInitialized())
        M_external_data_elements.push_back(&M_tcc);
    else
        throw std::runtime_error("forcingAtmosphere: One of M_Qlw_in or M_tcc should be initialised");

    // - snowfall can come from M_snowfall, M_snowfr*M_precip, or just M_precip (if M_tair<0)
    if(M_snowfr.isInitialized())
        M_external_data_elements.push_back(&M_snowfr);
    if(M_snowfall.isInitialized())
        M_external_data_elements.push_back(&M_snowfall);

    if(M_sphuma.isInitialized())
        // have specific humidity from the forcing
        M_external_data_elements.push_back(&M_sphuma);
    else if(M_mixrat.isInitialized())
        // have mixing ratio (simple relationship to specific humidity) from the forcing
        M_external_data_elements.push_back(&M_mixrat);
    else if(M_dair.isInitialized())
        // need to estimate the specific humidity from the dew point
        M_external_data_elements.push_back(&M_dair);
    else
        throw std::runtime_error("forcingAtmosphere: One of M_sphuma, M_mixrat or M_dair should be initialised");

}//forcingAtmosphere


//------------------------------------------------------------------------------------------------------
//! Nesting of forcing data.
//! !Does not seem to be used!
void
FiniteElement::forcingNesting()//(double const& u, double const& v)
{
    M_ice_thick=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 0,false,time_init);
    M_external_data_elements.push_back(&M_ice_thick);
    M_ice_conc=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 1,false,time_init);
    M_external_data_elements.push_back(&M_ice_conc);
    M_ice_snow_thick=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 2,false,time_init);
    M_external_data_elements.push_back(&M_ice_snow_thick);
    if ( Environment::vm()["thermo.newice_type"].as<int>() == 4 ) {
        M_ice_h_thin=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 3,false,time_init);
        M_external_data_elements.push_back(&M_ice_h_thin);
        M_ice_conc_thin=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 4,false,time_init);
        M_external_data_elements.push_back(&M_ice_conc_thin);
        M_ice_hs_thin=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 5,false,time_init);
        M_external_data_elements.push_back(&M_ice_hs_thin);
    }
    M_nesting_dist_elements=ExternalData(&M_nesting_distance_elements_dataset, M_mesh, 0,false,time_init);
    M_external_data_elements.push_back(&M_nesting_dist_elements);
    M_nesting_dist_nodes=ExternalData(&M_nesting_distance_nodes_dataset, M_mesh, 0,false,time_init);
    M_external_data_nodes.push_back(&M_nesting_dist_nodes);
    M_nesting_VT1=ExternalData(&M_nesting_nodes_dataset, M_mesh, 0,false,time_init);
    M_external_data_nodes.push_back(&M_nesting_VT1);
    M_nesting_VT2=ExternalData(&M_nesting_nodes_dataset, M_mesh, 1,false,time_init);
    M_external_data_nodes.push_back(&M_nesting_VT2);
    M_nesting_sigma1=ExternalData(&M_nesting_dynamics_elements_dataset, M_mesh, 0,false,time_init);
    M_external_data_elements.push_back(&M_nesting_sigma1);
    M_nesting_sigma2=ExternalData(&M_nesting_dynamics_elements_dataset, M_mesh, 1,false,time_init);
    M_external_data_elements.push_back(&M_nesting_sigma2);
    M_nesting_sigma3=ExternalData(&M_nesting_dynamics_elements_dataset, M_mesh, 2,false,time_init);
    M_external_data_elements.push_back(&M_nesting_sigma3);
    M_nesting_damage=ExternalData(&M_nesting_dynamics_elements_dataset, M_mesh, 3,false,time_init);
    M_external_data_elements.push_back(&M_nesting_damage);
    M_nesting_ridge_ratio=ExternalData(&M_nesting_dynamics_elements_dataset, M_mesh, 4,false,time_init);
    M_external_data_elements.push_back(&M_nesting_ridge_ratio);
}//forcingNesting

    
//------------------------------------------------------------------------------------------------------
//! Sets the physical variables relevant to the ocean according to the chosen ocean state and data (CONSTANT, TOPAZR, ...)
//! Called by the initForcings() function.
void
FiniteElement::forcingOcean()//(double const& u, double const& v)
{

    if(M_use_nesting)
    {
        if(M_use_ocean_nesting)
        {
            M_ocean_temp=ExternalData(&M_nesting_ocean_elements_dataset, M_mesh, 0,false,time_init);
            M_external_data_elements.push_back(&M_ocean_temp);
            M_ocean_salt=ExternalData(&M_nesting_ocean_elements_dataset, M_mesh, 1,false,time_init);
            M_external_data_elements.push_back(&M_ocean_salt);
        }
    }

    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
            M_ocean=ExternalData(
                vm["ideal_simul.constant_ocean_u"].as<double>(),
                vm["ideal_simul.constant_ocean_v"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());

            M_ssh=ExternalData(vm["ideal_simul.constant_ssh"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());

            if ( (!M_use_nesting) || ( (M_use_nesting) && (!M_use_ocean_nesting) ) )
            {
                M_ocean_temp=ExternalData(physical::ocean_freezing_temp);
                M_ocean_salt=ExternalData(physical::ocean_freezing_temp/physical::mu);
            }

            M_mld=ExternalData(vm["ideal_simul.constant_mld"].as<double>());
            break;

        case setup::OceanType::TOPAZR: case setup::OceanType::TOPAZF:
            M_ocean=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 0, true,
                time_init, vm["simul.spinup_duration"].as<double>());

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
                time_init, vm["simul.spinup_duration"].as<double>());

            if ( (!M_use_nesting) || ( (M_use_nesting) && (!M_use_ocean_nesting) ) )
            {
                M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
                M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
            }

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
    		break;

        case setup::OceanType::TOPAZR_ALTIMETER:
            M_ocean=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 0, true,
            time_init, vm["simul.spinup_duration"].as<double>());

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
            time_init, vm["simul.spinup_duration"].as<double>());

            if ( (!M_use_nesting) || ( (M_use_nesting) && (!M_use_ocean_nesting) ) )
            {
                M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
                M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
            }

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            break;

        case setup::OceanType::TOPAZR_atrest:
            M_ocean=ExternalData(
                vm["ideal_simul.constant_ocean_u"].as<double>(),
                vm["ideal_simul.constant_ocean_v"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
                time_init, vm["simul.spinup_duration"].as<double>());

            if ( (!M_use_nesting) || ( (M_use_nesting) && (!M_use_ocean_nesting) ) )
            {
                M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
                M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
            }

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            // SYL: there was a capping of the mld at minimum vm["ideal_simul.constant_mld"].as<double>()
            // but Einar said it is not necessary, so it is not implemented
    		break;

        default:
            std::cout << "invalid ocean forcing"<<"\n";
            throw std::logic_error("invalid ocean forcing");
    }

    // add the external data objects to M_external_data_nodes or M_external_data_elements
    // for looping
    M_external_data_nodes.push_back(&M_ocean);
    M_external_data_nodes.push_back(&M_ssh);
    M_external_data_elements.push_back(&M_ocean_temp);
    M_external_data_elements.push_back(&M_ocean_salt);
    M_external_data_elements.push_back(&M_mld);
}//forcingOcean


//------------------------------------------------------------------------------------------------------
//! Initializes variables relevant to a slab ocean (sss and sst).
//! Called by the initModelState() function.
void
FiniteElement::initSlabOcean()
{
    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
            //std::fill(M_sst.begin(), M_sst.end(), -1.8);
            std::fill(M_sst.begin(), M_sst.end(), 1.);
            std::fill(M_sss.begin(), M_sss.end(),  1.8/physical::mu);
            break;
        case setup::OceanType::TOPAZR:
        case setup::OceanType::TOPAZR_atrest:
        case setup::OceanType::TOPAZF:
        case setup::OceanType::TOPAZR_ALTIMETER:
            for ( int i=0; i<M_num_elements; ++i)
            {
                // Make sure the erroneous salinity and temperature don't screw up the initialisation too badly
                // This can still be done much better!
                M_sss[i] = std::max(physical::si, M_ocean_salt[i]);
                M_sst[i] = std::max(-M_sss[i]*physical::mu, M_ocean_temp[i]);
            }

            break;
        default:
            std::cout << "invalid ocean initialisation"<<"\n";
            throw std::logic_error("invalid ocean forcing");
    }
}//initSlabOcean
    
    
//------------------------------------------------------------------------------------------------------
//! Performs data assimilation (of sss and sst) in the slab ocean.
//! Called by the DataAssimilation() function.
void
FiniteElement::assimilateSlabOcean()
{
    double sigma_mod=1.;
    double sigma_obs=1.;
    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
            //std::fill(M_sst.begin(), M_sst.end(), -1.8);
            for ( int i=0; i<M_num_elements; ++i)
            {
                M_sss[i]=(sigma_obs*M_sss[i]+sigma_mod*1.8/physical::mu)/(sigma_obs+sigma_mod);
                M_sst[i]=(sigma_obs*M_sst[i]+sigma_mod*1.)/(sigma_obs+sigma_mod);
            }
            break;
        case setup::OceanType::TOPAZR:
        case setup::OceanType::TOPAZR_atrest:
        case setup::OceanType::TOPAZF:
        case setup::OceanType::TOPAZR_ALTIMETER:
            double sss_obs, sst_obs;
            for ( int i=0; i<M_num_elements; ++i)
            {
                // Make sure the erroneous salinity and temperature don't screw up the initialisation too badly
                // This can still be done much better!
                sss_obs=std::max(physical::si, M_ocean_salt[i]);
                sst_obs=std::max(-sss_obs*physical::mu, M_ocean_temp[i]);

                M_sss[i] = (sigma_obs*M_sss[i]+sigma_mod*sss_obs)/(sigma_obs+sigma_mod);
                M_sst[i] = (sigma_obs*M_sst[i]+sigma_mod*sst_obs)/(sigma_obs+sigma_mod);

                M_sst[i] = std::max(-M_sss[i]*physical::mu, M_sst[i]);
            }
            break;
        default:
            std::cout << "invalid ocean data assimilation"<<"\n";
            throw std::logic_error("invalid ocean forcing");
    }
}//assimilateSlabOcean

    
//------------------------------------------------------------------------------------------------------
//! Initializes the ice state (CONSTANT, TOPAZ4, PIOMAS, SMOS, ...).
//! Called by the initModelState() function.
void
FiniteElement::initIce()
{
    switch (M_ice_type)
    {
        case setup::IceType::CONSTANT:
        case setup::IceType::CONSTANT_PARTIAL:
            this->constantIce();
            break;
        case setup::IceType::TARGET:
            this->targetIce();
            break;
        case setup::IceType::TOPAZ4:
            this->topazIce();
            break;
        case setup::IceType::TOPAZ4OSISAFICESAT:
            this->topazIceOsisafIcesat();
            break;
        case setup::IceType::TOPAZ4F:
            this->topazForecastIce();
            break;
        case setup::IceType::TOPAZ4FAMSR2:
            this->topazForecastAmsr2Ice();
            break;
        case setup::IceType::TOPAZ4FAMSR2OSISAF:
            this->topazForecastAmsr2OsisafIce();
            break;
        case setup::IceType::TOPAZ4FAMSR2OSISAFNIC:
        case setup::IceType::TOPAZ4FAMSR2OSISAFNICWEEKLY:
            this->topazForecastAmsr2OsisafNicIce(M_ice_type==setup::IceType::TOPAZ4FAMSR2OSISAFNICWEEKLY);
            break;
        case setup::IceType::PIOMAS:
            this->piomasIce();
            break;
        case setup::IceType::AMSRE:
            this->topazAmsreIce();
            break;
        case setup::IceType::AMSR2:
            this->topazAmsr2Ice();
            break;
        case setup::IceType::CS2_SMOS:
            this->cs2SmosIce();
            break;
        case setup::IceType::CS2_SMOS_AMSR2:
            this->cs2SmosAmsr2Ice();
            break;
        case setup::IceType::SMOS:
            this->smosIce();
            break;
        default:
            std::cout << "invalid initialization of the ice"<<"\n";
            throw std::logic_error("invalid initialization of the ice");
    }

    // check consistency of fields after initialisation
    // - init ice temp everywhere
    this->checkConsistency();
}//initIce


//------------------------------------------------------------------------------------------------------
//! Checks on the consistency of fields and initial ice temperature after initialization and assimilation.
//! Called by the InitIce() and AssimilateIce() functions.
//  TODO If assimilating, we only need to initialize the ice temperature if new ice is created by the assimilation.
void
FiniteElement::checkConsistency()
{
    //check consistency of fields after init/assimilation
    //1. set things to zero if conc/abs thickness are too small
    //2. check SST is consistent
    //3. Initialise M_tice[i]
    for ( int i=0; i<M_num_elements; i++ )
    {
        //if conc or absolute thickness too small, set all fields to 0
        if (       (M_conc[i] < physical::cmin)
                || (M_thick[i] < M_conc[i]*physical::hmin))
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
            M_damage[i]=0.;
            M_ridge_ratio[i]=0.;
        }
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            if (       (M_conc_thin[i] < physical::cmin)
                    || (M_h_thin[i] < M_conc_thin[i]*physical::hmin))
            {
                M_conc_thin[i]=0.;
                M_h_thin[i]=0.;
                M_hs_thin[i]=0.;
            }
        }

        // freezing points of ice and water needed for init of ice temp
        // and to check SST
        double const Tfr_wtr = -physical::mu*M_sss[i];    //freezing point for water
        double const Tfr_ice = -physical::mu*physical::si;//freezing point for ice salinity

        // check SST is consistent
        double conc_tot = M_conc[i];
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            conc_tot += M_conc_thin[i];
        double const weight_conc = std::min(1., conc_tot*100.);//increases linearly from 0 to 1 as conc goes from 0 to .01
        if(conc_tot>0.)
            M_sst[i] = Tfr_wtr*weight_conc + M_sst[i]*(1-weight_conc);

        // init M_tice[j], M_tsurf_thin where it is needed
        // - initialization: wherever there is ice
        // - assimilation: wherever there is ice
        //   TODO make this just where new ice is added

        //init surface temp
        double tsurf = 0.;
        if ( M_snow_thick[i] > 0. )
        {
            // some snow
            // - can't be greater than melting point of snow (0degC)
            tsurf = std::min(0., M_tair[i]);
        }
        else
        {
            // no snow
            // - can't be greater than melting point of ice
            tsurf = std::min(Tfr_ice, M_tair[i]);
        }

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            //init thin ice temperature
            M_tsurf_thin[i] = tsurf;
        }

        //init thick ice temperature
        M_tice[0][i] = tsurf;

        //if using Winton, init T1 and T2
        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            if ( M_thick[i] <= 0. )
            {
                // Just set them to freezing point if no ice
                M_tice[1][i] = Tfr_ice;
                M_tice[2][i] = Tfr_ice;
            }
            else
            {
                // Calculate the temp at the top of the ice
                double Ti = M_tice[0][i];
                if(M_snow_thick[i]>0)
                {
                    // Calculate the temp of the ice-snow interface (Ti) using a zero-layer model
                    // => ki*(Tfr_bot - Ti)/hi = ks(Ti - Ts)/hs
                    // => ki*hs*(Tfr_bot - Ti) = ks*hi*(Ti - Ts)
                    // => (ki*hs+ks*hi)*Ti = ki*hs*Tfr_bot + ks*hi*Ts
                    // => a*Ti = b
                    double const a = physical::ki*M_snow_thick[i]
                        + physical::ks*M_thick[i];
                    double const b = physical::ki*M_snow_thick[i]*Tfr_wtr
                        + physical::ks*M_thick[i]*M_tice[0][i];
                    Ti = std::min(b/a, Tfr_ice);//make sure it is not higher than freezing point
                }

                // Then use linear interpolation between bottom and top of ice
                M_tice[1][i] = Tfr_wtr + .75*(Ti - Tfr_wtr);
                M_tice[2][i] = Tfr_wtr + .25*(Ti - Tfr_wtr);
            }
        }//Winton
    }
}//checkConsistency

    
//------------------------------------------------------------------------------------------------------
//! Performs the sea ice data assimilation. Different data sets available. Calls the appropriate function depending on the dataset.
//! Called by the DataAssimilation() function.
void
FiniteElement::assimilateIce()
{
    switch (M_ice_type)
    {
        case setup::IceType::TOPAZ4FAMSR2OSISAF:
            this->assimilate_topazForecastAmsr2OsisafIce();
            break;
        case setup::IceType::TOPAZ4FAMSR2OSISAFNIC:
        case setup::IceType::TOPAZ4FAMSR2OSISAFNICWEEKLY:
            this->assimilate_topazForecastAmsr2OsisafNicIce(M_ice_type==setup::IceType::TOPAZ4FAMSR2OSISAFNICWEEKLY);
            break;
        default:
            std::cout << "invalid choice for data assimilation of the ice"<<"\n";
            throw std::logic_error("invalid choice for data assimilation of the ice");
    }

    // Check consistency of fields and init ice temperature.
    // We only need to initialize the ice temperature
    // - if new ice is created by the assimilation
    // - TODO determine where to do this in the individual assimilation routines
    this->checkConsistency();
}//assimilateIce


//------------------------------------------------------------------------------------------------------
//! Sets the ice cover to a homogeneous state.
//! Called by the initIce() function.
void
FiniteElement::constantIce()
{
    LOG(DEBUG) <<"Constant Ice\n";
    double c_const = vm["ideal_simul.init_concentration"].as<double>();
    double h_const = vm["ideal_simul.init_thickness"].as<double>();
    double hs_const = vm["ideal_simul.init_snow_thickness"].as<double>();
    std::fill(M_conc.begin(), M_conc.end(), c_const);
    std::fill(M_thick.begin(), M_thick.end(), c_const*h_const);//M_thick is ice volume
    std::fill(M_snow_thick.begin(), M_snow_thick.end(), hs_const);
    std::fill(M_damage.begin(), M_damage.end(), 0.);

    if (M_ice_type==setup::IceType::CONSTANT_PARTIAL)
    {
        double xmin, xmax;
        if(M_rank==0)
        {
            auto Bx = M_mesh_root.coordX();//xmin,xmax from nodes of global mesh
            xmin = *std::min_element(Bx.begin(), Bx.end());
            xmax = *std::max_element(Bx.begin(), Bx.end());
        }
        boost::mpi::broadcast(M_comm, xmin, 0);
        boost::mpi::broadcast(M_comm, xmax, 0);
        double xedge = xmin + 0.3*(xmax-xmin);

        std::cout<<"In constantIce (partial cover)\n";
        std::cout<<"M_ice_type "<< (int)M_ice_type<<"\n";
        std::cout<<"Min conc = "<< *std::min_element(M_conc.begin(),M_conc.end()) <<"\n";
        std::cout<<"Max conc = "<< *std::max_element(M_conc.begin(),M_conc.end()) <<"\n";
        std::cout<<"Min thick = "<< *std::min_element(M_thick.begin(),M_thick.end()) <<"\n";
        std::cout<<"Max thick = "<< *std::max_element(M_thick.begin(),M_thick.end()) <<"\n";
        std::cout<<"xmin="<<xmin<<"\n";
        std::cout<<"xmax="<<xmax<<"\n";
        std::cout<<"xedge="<<xedge<<"\n";

        auto Bx = M_mesh.bCoordX();//set conc, etc on elements
        for (int i=0; i<M_conc.size(); ++i)
        {
            if (Bx[i] < xedge)
            {
                M_conc[i]       = 0.;
                M_thick[i]      = 0.;
                M_snow_thick[i] = 0.;
            }
        }
        std::cout<<"New min conc = "<< *std::min_element(M_conc.begin(),M_conc.end()) <<"\n";
        std::cout<<"New max conc = "<< *std::max_element(M_conc.begin(),M_conc.end()) <<"\n";
        std::cout<<"New min thick = "<< *std::min_element(M_thick.begin(),M_thick.end()) <<"\n";
        std::cout<<"New max thick = "<< *std::max_element(M_thick.begin(),M_thick.end()) <<"\n";
        //std::abort();
    }//partial ice cover

}//constantIce


//------------------------------------------------------------------------------------------------------
//! Sets the ice cover to a target state.
//! Called by the initIce() function.
void
FiniteElement::targetIce()
{
    /*
    double y_max=300000.;
    double y_min=150000.;
    double x_max1=240000.;
    double x_min1=150000.;
    double x_max2=300000.;
    double x_min2=260000.;
    */

    double y_max=350000.;
    double y_min=0000.;
    double x_max=400000.;
    double x_min=200000.;

    double transition=(x_max-x_min)/10.;

    double tmp_var;

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    double cmin= 0.;

    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_var = (RY[i]<=y_max)*(RY[i]>=y_min)*(RX[i]<=x_max)*(RX[i]>=x_min);

/*        tmp_var = (RY[i]<=y_max)*(RY[i]>=y_min)*(RX[i]<=x_max1)*(RX[i]>=x_min1)
                  + (RY[i]<=y_max)*(RY[i]>=y_min)*(RX[i]<=x_max2)*(RX[i]>=x_min2);
*/
            /*
            +     (RY[i]<=y_max)*(RY[i]>=y_min)*(RX[i]<x_min)*std::max(cmin,(1.-std::hypot(RX[i]-x_min,0.         )/transition))
            +     (RY[i]<=y_max)*(RY[i]>=y_min)*(RX[i]>x_max)*std::max(cmin,(1.-std::hypot(RX[i]-x_max,0.         )/transition))
            +     (RX[i]<=x_max)*(RX[i]>=x_min)*(RY[i]<y_min)*std::max(cmin,(1.-std::hypot(RY[i]-y_min,0.         )/transition))
            +     (RX[i]<=x_max)*(RX[i]>=x_min)*(RY[i]>y_max)*std::max(cmin,(1.-std::hypot(RY[i]-y_max,0.         )/transition))
            +     (RY[i]< y_min)*(RX[i]< x_min)              *std::max(cmin,(1.-std::hypot(RX[i]-x_min,RY[i]-y_min)/transition))
            +     (RY[i]< y_min)*(RX[i]> x_max)              *std::max(cmin,(1.-std::hypot(RX[i]-x_max,RY[i]-y_min)/transition))
            +     (RY[i]> y_max)*(RX[i]< x_min)              *std::max(cmin,(1.-std::hypot(RX[i]-x_min,RY[i]-y_max)/transition))
            +     (RY[i]> y_max)*(RX[i]> x_max)              *std::max(cmin,(1.-std::hypot(RX[i]-x_max,RY[i]-y_max)/transition));
            */
        std::cout<<"RX: "<< RX[i] << "RY: "<< RY[i] << "tmp_var: " << tmp_var << "\n";

        //M_conc[i]  = std::max(vm["ideal_simul.init_concentration"].as<double>()*tmp_var,cmin);
    //    M_thick[i] = vm["ideal_simul.init_thickness"].as<double>()*M_conc[i];
    //    M_snow_thick[i] = vm["ideal_simul.init_snow_thickness"].as<double>()*M_conc[i];

        M_conc[i]  = vm["ideal_simul.init_concentration"].as<double>();

    // if(i==10)
    //     M_conc[i]=0.;

        M_thick[i] = vm["ideal_simul.init_thickness"].as<double>()*M_conc[i];
        M_snow_thick[i] = vm["ideal_simul.init_snow_thickness"].as<double>()*M_conc[i];

        M_damage[i]=0.;

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_conc_thin[i]  = vm["ideal_simul.init_thin_conc"].as<double>();

            M_h_thin[i]     = (vm["thermo.h_thin_min"].as<double>()+(vm["thermo.h_thin_max"].as<double>()-vm["thermo.h_thin_min"].as<double>())/2.)*M_conc_thin[i];

            M_hs_thin[i]    = vm["ideal_simul.init_snow_thickness"].as<double>()*M_conc_thin[i];
        }

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
            M_damage[i]=0.;
        }

        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
            M_damage[i]=0.;
        }
    }
}//targetIce

    
//------------------------------------------------------------------------------------------------------
//! Initializes the ice and snow states from Topaz outputs.
//! Called by the initIce() function.
void
FiniteElement::topazIce()
{
    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_init_conc);
    external_data_tmp.push_back(&M_init_thick);
    external_data_tmp.push_back(&M_init_snow_thick);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - TOPAZ ice ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);
    external_data_tmp.resize(0);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
        // - TOPAZ puts very small values instead of 0.
        // - uses absolute thickness not effective thickness
        tmp_var=std::min(1.,M_init_conc[i]);
        M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.;
        tmp_var=M_init_thick[i];
        M_thick[i] = (tmp_var>1e-14) ? tmp_var*M_conc[i] : 0.; // TOPAZ puts very small values instead of 0.
        tmp_var=M_init_snow_thick[i];
        M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var*M_conc[i] : 0.; // TOPAZ puts very small values instead of 0.

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }

        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

        M_damage[i]=0.;
    }
}//topazIce


// -----------------------------------------------------------------------------------------------------------
//! Initializes ice state from Topaz and Osisaf data.
//! Called by the initIce() function.
void
FiniteElement::topazIceOsisafIcesat()
{
    //topaz
    external_data M_topaz_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    external_data M_topaz_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    external_data M_topaz_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);

    //obs
    external_data M_osisaf_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);
    external_data M_osisaf_conc=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,0,false,time_init);
    external_data M_icesat_thick=ExternalData(&M_ice_icesat_elements_dataset,M_mesh,0,false,time_init);
    external_data M_amsre_conc=ExternalData(&M_ice_amsre_elements_dataset,M_mesh,0,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_topaz_conc);
    external_data_tmp.push_back(&M_topaz_thick);
    external_data_tmp.push_back(&M_topaz_snow_thick);
    external_data_tmp.push_back(&M_osisaf_type);
    external_data_tmp.push_back(&M_osisaf_conc);
    external_data_tmp.push_back(&M_icesat_thick);
    external_data_tmp.push_back(&M_amsre_conc);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - TOPAZ/OSISAF/Icesat ice ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    for (int i=0; i<M_num_elements; ++i)
    {
        // - TOPAZ puts very small values instead of 0.
        // - uses absolute thickness not effective thickness
        double tmp_var=std::min(1.,M_topaz_conc[i]);
        M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.;
        tmp_var=M_topaz_thick[i];
        double hi = (tmp_var>1e-14) ? tmp_var : 0.;
        tmp_var=M_topaz_snow_thick[i];
        double hs = (tmp_var>1e-14) ? tmp_var : 0.;

        if(M_conc[i]>0.) // use osisaf only where topaz says there is ice to avoid near land issues and fake concentration over the ocean
            M_conc[i]=M_osisaf_conc[i];

        //M_type[i]==1. // No ice
        //M_type[i]==2. // First-Year ice
        //M_type[i]==3. // Multi-Year ice
        //M_type[i]==4. // Mixed
        double ratio_FYI=0.3;
        double ratio_MYI=0.9;
        double ratio_Mixed=0.5*(ratio_FYI+ratio_MYI);

        double thick_FYI=hi;
        double thick_MYI=std::max(M_icesat_thick[i],hi);//NB icesat outputs absolute thickness
        double thick_Mixed=0.5*(thick_FYI+thick_MYI);

        if((hi>0.)&&(M_conc[i])>0.2)
        {
            if(M_mesh_filename.find("kara") != std::string::npos)
            {
                LOG(DEBUG) <<"Type information is not used for the kara mesh, "
                    <<"we assume there is only FYI\n";
                M_ridge_ratio[i]=ratio_FYI;
                hi = thick_FYI;
            } 
            else
            {
                if(M_osisaf_type[i]<1.5)//1. no ice in OSISAF
                {
                    M_ridge_ratio[i]=0.;
                    hi=thick_FYI;
                }
                else if(M_osisaf_type[i]<2.5)//2. FYI
                {
                    M_ridge_ratio[i]=ratio_FYI;
                    hi=thick_FYI;
                }
                else if(M_osisaf_type[i]<3.5)//3. MYI
                {
                    M_ridge_ratio[i]=ratio_MYI;
                    hi  = thick_MYI;
                }
                else if(M_osisaf_type[i]<4.5)//4. mixed
                {
                    M_ridge_ratio[i]=ratio_Mixed;
                    hi      =thick_Mixed;
                }
                else//can't happen, can it?
                {
                    M_ridge_ratio[i]=ratio_Mixed;
                    hi=thick_Mixed;
                }
            }//not Kara
        }//ice present
        else
        {
            M_ridge_ratio[i]=0.;
        }

        if ( M_conc[i] < 0.01 || hi < physical::hmin )
        {
            //if either c or h equal zero, we set the others to zero as well
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
            M_ridge_ratio[i]=0.;
        }
        else
        {
            //convert from absolute to effective thickness
            M_ridge_ratio[i]=M_ridge_ratio[i]*M_conc[i];
            M_thick[i]=hi*M_conc[i];
            M_snow_thick[i]=hs*M_conc[i];
        }

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_conc_thin[i]=std::max(M_amsre_conc[i]-M_conc[i],0.);
            M_h_thin[i]=M_conc_thin[i]*(h_thin_min+0.5*(h_thin_max-h_thin_min));
        }

        M_damage[i]=0.;
    }
}//topazIceOsisafIcesat
    
    
// -----------------------------------------------------------------------------------------------------------
//! Initializes ice state from Topaz forecasts.
//! Called by the initIce() function.
void
FiniteElement::topazForecastIce()
{
    external_data topaz_conc = ExternalData(&M_ocean_elements_dataset, M_mesh, 3, false, time_init);
    // NB TOPAZ gives absolute ice/snow thickness
    // - cell_methods = "mean_where_ice"
    external_data topaz_thick = ExternalData(&M_ocean_elements_dataset, M_mesh, 4, false, time_init);
    external_data topaz_snow_thick = ExternalData(&M_ocean_elements_dataset, M_mesh, 5, false, time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&topaz_conc);
    external_data_tmp.push_back(&topaz_thick);
    external_data_tmp.push_back(&topaz_snow_thick);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - TOPAZ ice forecast ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
        // - TOPAZ puts very small values instead of 0.
        // - uses absolute thickness not effective thickness
        tmp_var = std::min(1., topaz_conc[i]);
        M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.;
        tmp_var = topaz_thick[i];
        M_thick[i] = (tmp_var>1e-14) ? tmp_var*M_conc[i] : 0.;
        tmp_var = topaz_snow_thick[i];
        M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var*M_conc[i] : 0.;

        //if either c or h equal zero, we set the others to zero as well
        if ( M_conc[i] < 0.01 || M_thick[i] < (M_conc[i]*physical::hmin) )
        {
            M_conc[i]        = 0.;
            M_thick[i]       = 0.;
            M_snow_thick[i]  = 0.;
        }

        //init damage and ridge ratio to 0.
        M_ridge_ratio[i] = 0.;
        M_damage[i]=0.;
    }
}//topazForecastIce
    

// -----------------------------------------------------------------------------------------------------------
//! Initializes ice state from Topaz-AMSR2.
//! Called by the initIce() function.
void
FiniteElement::topazForecastAmsr2Ice()
{
    double init_conc_tmp;

    external_data M_conc_amsr2=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init-0.5);
    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    external_data M_init_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_conc_amsr2);
    external_data_tmp.push_back(&M_init_conc);
    external_data_tmp.push_back(&M_init_thick);
    external_data_tmp.push_back(&M_init_snow_thick);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - TOPAZ ice forecast/AMSR2 ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
        double uncertainty;
        if(M_conc_amsr2[i]<0.1)
            uncertainty=0.1;
        else
            uncertainty=0.05;

        double diff_mod_obs = M_conc_amsr2[i]-M_init_conc[i];
        if(std::abs(diff_mod_obs)>=uncertainty && M_conc_amsr2[i]<=1.)
            // NB missing value for AMSR2 when not over land is 1.15
            // move towards AMSR2 value by the amount uncertainty/2
            M_conc[i] = std::min(1., M_conc_amsr2[i]-(diff_mod_obs/std::abs(diff_mod_obs))*uncertainty/2.);
        else
            M_conc[i] = std::min(1., M_init_conc[i]);

        // TOPAZ puts very small values instead of 0.
        tmp_var=M_init_conc[i];
        init_conc_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
        tmp_var=M_init_thick[i];
        M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
        tmp_var=M_init_snow_thick[i];
        M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

        // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
        if(init_conc_tmp>0.05)
        {
            //convert from absolute to effective thickness
            M_thick[i]*=M_conc[i];
            M_snow_thick[i]*=M_conc[i];
        }

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_thick[i]=0.;
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<0.1*M_conc[i])
            //10cm min ice thickness
            M_thick[i]=0.1*M_conc[i];

        M_damage[i]=1.-M_conc[i];
    }
}//topazForecastAmsr2Ice

    
// -----------------------------------------------------------------------------------------------------------
//! Get the maximum and minimum ice concentration corresponding to the original
//! NIC ice charts
//! called by <FiniteElement::topazForecastAmsr2OsisafNicIce>() and
//! <FiniteElement::assimilate_topazForecastAmsr2OsisafNicIce>()
void
FiniteElement::concBinsNic(double &thin_conc_obs_min, double &thin_conc_obs_max,
        double ci, bool use_weekly_nic)
{

    if(ci<=0.)
    {
        thin_conc_obs_min = 0.;
        thin_conc_obs_max = 0.;
    }
    else if(!use_weekly_nic)
    {
        if(ci<=0.45) // CT18
        {
            thin_conc_obs_min = 0.1;
            thin_conc_obs_max = 0.8;
        }
        else if(ci<=0.9) // CT81
        {
            thin_conc_obs_min = 0.8;
            thin_conc_obs_max = 1.;
        }
    }
    else
    {
        if(ci<=0.2) // CT13
        {
            thin_conc_obs_min = 0.1;
            thin_conc_obs_max = 0.3;
        }
        else if(ci<=0.30) // CT24
        {
            thin_conc_obs_min = 0.2;
            thin_conc_obs_max = 0.4;
        }
        else if(ci<=0.50) // CT46
        {
            thin_conc_obs_min = 0.4;
            thin_conc_obs_max = 0.6;
        }
        else if(ci<=0.70) //CT68
        {
            thin_conc_obs_min = 0.6;
            thin_conc_obs_max = 0.8;
        }
        else if(ci<=0.90) // CT81
        {
            thin_conc_obs_min = 0.8;
            thin_conc_obs_max = 1.0;
        }
        else if(ci<=1.) // CT92
        {
            thin_conc_obs_min = 0.9;
            thin_conc_obs_max = 1.0;
        }
    }
}//concBinsNic

    
// -----------------------------------------------------------------------------------------------------------
//! Assimilates Topaz forecasts, Amsr2, Osisaf and Nic ice state data
//! Called by the assimilateIce() function.
void
FiniteElement::assimilate_topazForecastAmsr2OsisafNicIce(bool use_weekly_nic)
{
    double real_thickness, init_conc_tmp;


    external_data_vec external_data_tmp;
    external_data M_nic_conc;
    if(use_weekly_nic)
        M_nic_conc = ExternalData(&M_ice_nic_weekly_elements_dataset,
                M_mesh, 0, false, time_init-0.5);
    else
        M_nic_conc = ExternalData(&M_ice_nic_elements_dataset,
                M_mesh, 0, false, time_init-0.5);
    external_data_tmp.push_back(&M_nic_conc);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"assimilate - NIC ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    double tmp_var;
    double sigma_mod=0.1;
    double sigma_amsr2=0.1;
    double sigma_osisaf=0.1;

    double topaz_conc, topaz_thick;
    double h_model, c_model;
    double hi;
    for (int i=0; i<M_num_elements; ++i)
    {
        h_model=M_thick[i];
        c_model=M_conc[i];

        if(c_model>0.01)
        {
            M_thick[i]=(h_model/c_model)*M_conc[i];            
            M_ridge_ratio[i]=(M_ridge_ratio[i]/c_model)*M_conc[i]; 
            M_damage[i]=(M_damage[i]/c_model)*M_conc[i];
        }
        else
        {
            M_thick[i]=0.;
            M_ridge_ratio[i]=0.;
        }

        //if either c or h equal zero, we set the others to zero as well
        hi=M_thick[i]/M_conc[i];
        if(M_conc[i]<0.1)
            hi=M_thick[i];

        if ( M_conc[i] < 0.01 || hi < physical::hmin )
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
            M_ridge_ratio[i]=0.;
        }

        // don't change model if NIC conc > 1
        // - then it is masked
        // - unfortunately, applying the mask in datasets led to masked values being treated as real values of 0.0
        //   so we have to do it manually here
        if(M_nic_conc[i]>1.)
            continue;


        // Use the NIC ice charts
        // - get conc bins from NIC dataset
        double thin_conc_obs = 0.;
        double thin_conc_obs_min = 0.;
        double thin_conc_obs_max = 0.;
        this->concBinsNic(thin_conc_obs_min, thin_conc_obs_max, M_nic_conc[i], use_weekly_nic);

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {

            if((M_conc[i]+M_conc_thin[i])<thin_conc_obs_min)
            {
                thin_conc_obs = thin_conc_obs_min-M_conc[i];//always >=0?
                if(thin_conc_obs>=0.)
                {
                    //increase thin ice conc so total conc = thin_conc_obs_min
                    if(thin_conc_obs>M_conc_thin[i])
                        //increase thin ice vol
                        M_h_thin[i] = M_h_thin[i]+(h_thin_min + (h_thin_max/2.-h_thin_min)*0.5)*(thin_conc_obs-M_conc_thin[i]);
                    else
                        //reduce thin ice vol
                        M_h_thin[i] = M_h_thin[i]*thin_conc_obs/M_conc_thin[i];

                    M_conc_thin[i] = thin_conc_obs;
                }
#if 0
                else
                {
                    //not possible?
                    M_conc_thin[i]=0.;
                    M_h_thin[i]=0.;

                    M_thick[i]=M_thick[i]/(M_conc[i])*(M_conc[i]+thin_conc_obs);
                    M_conc[i]=M_conc[i]+thin_conc_obs;
                }
#endif
            }
            else if((M_conc[i]+M_conc_thin[i])>thin_conc_obs_max)
            {
                thin_conc_obs = thin_conc_obs_max-M_conc[i];
                if(thin_conc_obs>=0.)
                {
                    //some thin ice
                    if(thin_conc_obs>M_conc_thin[i])
                        M_h_thin[i] = M_h_thin[i]+(h_thin_min + (h_thin_max/2.-h_thin_min)*0.5)*(thin_conc_obs-M_conc_thin[i]);
                    else
                        M_h_thin[i] = M_h_thin[i]*thin_conc_obs/M_conc_thin[i];

                    M_conc_thin[i] = thin_conc_obs;
                }
                else
                {
                    //no thin ice
                    M_conc_thin[i]=0.;
                    M_h_thin[i]=0.;

                    // reduce thick ice to max value
                    M_thick[i]=M_thick[i]*(M_conc[i]+thin_conc_obs)/(M_conc[i]);
                    M_conc[i]=M_conc[i]+thin_conc_obs;
                }
            }

            /* Two cases: Thin ice fills the cell or not */
            double min_h_thin = h_thin_min*M_conc_thin[i];
            if ( M_h_thin[i] < min_h_thin )
                M_h_thin[i] = min_h_thin;

            double max_h_thin=(h_thin_min+(h_thin_max+h_thin_min)/2.)*M_conc_thin[i];
            if ( M_h_thin[i] > max_h_thin)
                M_h_thin[i] = max_h_thin;
        }//using thin ice
        else
        {
            if(M_conc[i]<thin_conc_obs_min)
            {
                //thin_conc_obs = .25*thin_conc_obs_max + .75*thin_conc_obs_min;
                thin_conc_obs = ( thin_conc_obs_min + (thin_conc_obs_min+thin_conc_obs_max)/2.) /2.;
                M_thick[i] = M_thick[i] + std::max(hi,0.5)*(thin_conc_obs-M_conc[i]); // 50 cm minimum for the added ice
                M_conc[i] = thin_conc_obs;
            }
            else if(M_conc[i]>thin_conc_obs_max)
            {
                //thin_conc_obs = .75*thin_conc_obs_max + .25*thin_conc_obs_min;
                thin_conc_obs = ( thin_conc_obs_max + (thin_conc_obs_min+thin_conc_obs_max)/2.) /2.;
                M_thick[i] = M_thick[i]*thin_conc_obs/M_conc[i];
                M_conc[i] = thin_conc_obs;
            }
        }//not using thin ice
    }//loop over elements
}//assimilate_topazForecastAmsr2OsisafNicIce

    
// -----------------------------------------------------------------------------------------------------------
//! Assimilates Topaz forecast, Amsr2 and Osisaf ice data.
//! Called by the assimilateIce() function.
void
FiniteElement::assimilate_topazForecastAmsr2OsisafIce()
{
    double real_thickness, init_conc_tmp;

    external_data M_osisaf_conc=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,0,false,time_init-0.5);

    external_data M_osisaf_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init-0.5);

    external_data M_amsr2_conc=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init-0.5);

    external_data M_topaz_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);

    external_data M_topaz_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);

    external_data M_topaz_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);
    DataSet dist2coast_elements_dataset=DataSet("dist2coast_elements");
    external_data M_dist2coast = ExternalData(&dist2coast_elements_dataset,M_mesh,0,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_osisaf_conc);
    external_data_tmp.push_back(&M_osisaf_type);
    external_data_tmp.push_back(&M_amsr2_conc);
    external_data_tmp.push_back(&M_dist2coast);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"assimilate - OSISAF/AMSR2/dist2coast ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init-0.5, RX, RY);

    external_data_tmp.resize(0);
    external_data_tmp.push_back(&M_topaz_conc);
    external_data_tmp.push_back(&M_topaz_thick);
    external_data_tmp.push_back(&M_topaz_snow_thick);
    if(M_rank==0)
        LOG(DEBUG)<<"assimilate - TOPAZ ice forecast ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    double tmp_var;
    double sigma_mod=1.;
    double sigma_amsr2=0.5;
    double sigma_osisaf=2.;

    double topaz_conc, topaz_thick;
    double h_model, c_model;

    for (int i=0; i<M_num_elements; ++i)
    {
        //initial fields
        h_model=M_thick[i];
        c_model=M_conc[i];

        topaz_conc = (M_topaz_conc[i]>1e-14) ? M_topaz_conc[i] : 0.; // TOPAZ puts very small values instead of 0.
        topaz_thick = (M_topaz_thick[i]>1e-14) ? M_topaz_thick[i] : 0.; // TOPAZ puts very small values instead of 0.

        if(((topaz_conc>0.)||(M_conc[i]>0.))
                && (M_osisaf_conc[i]>.15)
                && (M_dist2coast[i]>25.e3))
            // use osisaf only
            // - where topaz or the model says there is ice to avoid near land issues and fake concentration over the ocean
            // - where its conc is > .15 (can be trusted)
            // - also take into account distance to coast
            M_conc[i] = (sigma_osisaf*M_conc[i]+sigma_mod*M_osisaf_conc[i])/(sigma_osisaf+sigma_mod);

        if((M_amsr2_conc[i]<M_conc[i])
                && (M_amsr2_conc[i]>.15))
            // AMSR2 is higher resolution and sees small opening that would not be see in OSISAF
            // NB AMSR2 = 1.15 if missing data over ocean, however this will not affect the example here
            M_conc[i] = M_amsr2_conc[i];

        //tmp_var=M_topaz_snow_thick[i];
        //M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
        //if(M_conc[i]<M_topaz_conc[i])
         //   M_snow_thick[i] *= M_conc[i]/M_topaz_conc[i];


        if(c_model>0.01)
        {
            M_thick[i]=(h_model/c_model)*M_conc[i];            
            M_ridge_ratio[i]=(M_ridge_ratio[i]/c_model)*M_conc[i]; 
            M_damage[i]=(M_damage[i]/c_model)*M_conc[i];
        }
        else
        {
            M_thick[i]=0.;
            M_ridge_ratio[i]=0.;
        }

        //if either c or h equal zero, we set the others to zero as well
        if ( M_conc[i] < 0.01 || M_thick[i] < (M_conc[i]*physical::hmin) )
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
            M_ridge_ratio[i]=0.;
        }

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            double thin_conc_obs  = std::max(M_amsr2_conc[i]-M_conc[i],0.);

            M_conc_thin[i] = (sigma_osisaf*M_conc_thin[i]+sigma_mod*thin_conc_obs)/(sigma_amsr2+sigma_mod);

            /* Two cases: Thin ice fills the cell or not */
            double min_h_thin = h_thin_min*M_conc_thin[i];
            if ( M_h_thin[i] < min_h_thin )
                M_h_thin[i] = min_h_thin;

            double max_h_thin=(h_thin_min+(h_thin_max+h_thin_min)/2.)*M_conc_thin[i];
            if ( M_h_thin[i] > max_h_thin)
                M_h_thin[i] = max_h_thin;
        }
    }
}//assimilate_topazForecastAmsr2OsisafIce

    
// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from Topaz forecast, AMSR2 and Osisaf data.
//! Called by the initIce() function.
void
FiniteElement::topazForecastAmsr2OsisafIce()
{
    double real_thickness, init_conc_tmp;

    // observations
    external_data M_osisaf_conc=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,0,false,time_init-0.5);
    external_data M_osisaf_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init-0.5);
    external_data M_amsr2_conc=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init-0.5);

    // topaz
    external_data M_topaz_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    external_data M_topaz_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    external_data M_topaz_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_osisaf_conc);
    external_data_tmp.push_back(&M_osisaf_type);
    external_data_tmp.push_back(&M_amsr2_conc);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - OSISAF/AMSR2 ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init-0.5, RX, RY);

    external_data_tmp.resize(0);
    external_data_tmp.push_back(&M_topaz_conc);
    external_data_tmp.push_back(&M_topaz_thick);
    external_data_tmp.push_back(&M_topaz_snow_thick);
    if(M_rank==0)
        LOG(DEBUG)<<"init - TOPAZ ice forecast ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    for (int i=0; i<M_num_elements; ++i)
    {
        // TOPAZ puts very small values instead of 0,
        // so set things to zero if below a threshold

        // get absolute ice and snow thicknesses
        double tmp_var=M_topaz_thick[i];
        double hi  = (tmp_var>1e-14) ? tmp_var : 0.;// absolute thickness
        tmp_var=M_topaz_snow_thick[i];
        double hs  = (tmp_var>1e-14) ? tmp_var : 0.;// absolute snow thickness

        tmp_var=std::min(1.,M_topaz_conc[i]);
        M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.;
        if(M_conc[i]>0.)
            // use osisaf only where topaz says there is ice
            // to avoid near land issues and fake concentration over the ocean
            M_conc[i] = M_osisaf_conc[i];

        if(M_amsr2_conc[i]<M_conc[i])
            // AMSR2 is higher resolution and sees
            // small openings that would not be see in OSISAF
            // NB AMSR2 = 1.15 if missing data over ocean, however this will not affect the example here
            M_conc[i] = M_amsr2_conc[i];


        //M_type[i]==1. // No ice
        //M_type[i]==2. // First-Year ice
        //M_type[i]==3. // Multi-Year ice
        //M_type[i]==4. // Mixed
        double ratio_FYI=0.3;
        double ratio_MYI=0.9;
        double ratio_Mixed=0.5*(ratio_FYI+ratio_MYI);

        double thickfac_FYI=1.;
        double thickfac_MYI=1.5;
        double thickfac_Mixed=0.5*(thickfac_FYI+thickfac_MYI);

        if( (hi>0.) && (M_conc[i])>0.2 )
        {

            if(M_mesh_basename.find("kara") != std::string::npos)
            {
                LOG(DEBUG) <<"Type information is not used for the kara meshes, "
                    <<"we assume there is only FYI\n";
                M_ridge_ratio[i]=ratio_FYI;
                hi*=thickfac_FYI;
            }
            else
            {
                if(M_osisaf_type[i]<1.5)//1. no ice in OSISAF
                {
                    M_ridge_ratio[i]=0.;
                    hi*=thickfac_FYI;
                }
                else if(M_osisaf_type[i]<2.5)//2. FYI
                {
                    M_ridge_ratio[i]=ratio_FYI;
                    hi*=thickfac_FYI;
                }
                else if(M_osisaf_type[i]<3.5)//3. MYI
                {
                    M_ridge_ratio[i]=ratio_MYI;
                    hi*=thickfac_MYI;
                }
                else if(M_osisaf_type[i]<4.5)//4. mixed
                {
                    M_ridge_ratio[i]=ratio_Mixed;
                    hi*=thickfac_Mixed;
                }
                else if(M_osisaf_type[i]>4.)//can't happen, can it?
                {
                    M_ridge_ratio[i]=ratio_Mixed;
                    hi*=thickfac_Mixed;
                }
            }//not Kara
        }//ice present
        else
        {
            M_ridge_ratio[i]=0.;
            hi=0.;
        }

        if ( M_conc[i] < 0.01 || hi < physical::hmin )
        {
            //if either c or h equal zero, we set the others to zero as well
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
            M_ridge_ratio[i]=0.;
            M_damage[i]=0.;
        }
        else
        {
            // convert from absolute to effective thickness
            M_ridge_ratio[i]=M_ridge_ratio[i]*M_conc[i];
            M_snow_thick[i] = M_conc[i]*hs;
            M_thick[i] = M_conc[i]*hi;
            M_damage[i]=1.-M_conc[i];
        }

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_conc_thin[i]=std::max(M_amsr2_conc[i]-M_conc[i],0.);
            M_h_thin[i]=M_conc_thin[i]*(h_thin_min+0.5*(h_thin_max-h_thin_min));
        }

    }//loop over elements
}//topazForecastAmsr2OsisafIce

    
// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from Topaz forecast, AMSR2, Osisaf and NIC ice charts data.
//! Called by the initIce() function.
void
FiniteElement::topazForecastAmsr2OsisafNicIce(bool use_weekly_nic)
{
    //observations
    external_data M_osisaf_conc=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,0,false,time_init-0.5);
    external_data M_osisaf_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init-0.5);
    external_data M_amsr2_conc=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init-0.5);

    //topaz
    external_data M_topaz_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    external_data M_topaz_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    external_data M_topaz_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_osisaf_conc);
    external_data_tmp.push_back(&M_osisaf_type);
    external_data_tmp.push_back(&M_amsr2_conc);

    external_data M_nic_conc;
    if(use_weekly_nic)
        M_nic_conc = ExternalData(&M_ice_nic_weekly_elements_dataset,
                M_mesh, 0, false, time_init-0.5);
    else
        M_nic_conc = ExternalData(&M_ice_nic_elements_dataset,
                M_mesh, 0, false, time_init-0.5);
    external_data_tmp.push_back(&M_nic_conc);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - OSISAF/AMSR2/NIC ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init-0.5, RX, RY);

    external_data_tmp.resize(0);
    external_data_tmp.push_back(&M_topaz_conc);
    external_data_tmp.push_back(&M_topaz_thick);
    external_data_tmp.push_back(&M_topaz_snow_thick);
    if(M_rank==0)
        LOG(DEBUG)<<"init - TOPAZ ice forecast ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    for (int i=0; i<M_num_elements; ++i)
    {
        // TOPAZ puts very small values instead of 0,
        // so set things to zero if below a threshold

        // get absolute ice and snow thicknesses
        double tmp_var=M_topaz_thick[i];
        double hi  = (tmp_var>1e-14) ? tmp_var : 0.;// absolute thickness
        tmp_var=M_topaz_snow_thick[i];
        double hs  = (tmp_var>1e-14) ? tmp_var : 0.;// absolute snow thickness

        tmp_var=std::min(1.,M_topaz_conc[i]);
        M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.;
        if(     (M_conc[i]>0.)
             && (M_amsr2_conc[i]>.15)
             && (M_amsr2_conc[i]<=1.))
            // use amsr2 only where
            // - topaz says there is ice to avoid near land issues and fake concentration over the ocean
            // - it is large enough to be trusted
            // - it is not masked (mask value = 1.15, but using masking from dataset fills missing values to 0)
            M_conc[i]=M_amsr2_conc[i];

        double ratio_FYI=0.3;
        double ratio_MYI=0.9;
        double ratio_Mixed=0.5*(ratio_FYI+ratio_MYI);

        double thickfac_FYI=1.;
        double thickfac_MYI=1.5;
        double thickfac_Mixed=0.5*(thickfac_FYI+thickfac_MYI);

        if( (hi>0.) && (M_conc[i])>0.2 )
        {

            if(M_mesh_basename.find("kara") != std::string::npos)
            {
                LOG(DEBUG) <<"Type information is not used for the kara meshes,"
                    << " we assume there is only FYI\n";
                M_ridge_ratio[i]=ratio_FYI;
                //M_thick[i]=thick_FYI;
                hi *= thickfac_FYI;
            }
            else
            {
                if(M_osisaf_type[i]<1.5)//1. no ice in OSISAF
                {
                    M_ridge_ratio[i]=0.;
                    hi*=thickfac_FYI;
                }
                else if(M_osisaf_type[i]<2.5)//2. FYI
                {
                    M_ridge_ratio[i]=ratio_FYI;
                    hi*=thickfac_FYI;
                }
                else if(M_osisaf_type[i]<3.5)//3. MYI
                {
                    M_ridge_ratio[i]=ratio_MYI;
                    hi*=thickfac_MYI;
                }
                else if(M_osisaf_type[i]<4.5)//4. mixed
                {
                    M_ridge_ratio[i]=ratio_Mixed;
                    hi*=thickfac_Mixed;
                }
                else//can't happen, can it?
                {
                    M_ridge_ratio[i]=ratio_Mixed;
                    hi*=thickfac_Mixed;
                }
            }//not Kara
        }//ice present
        else
        {
            M_ridge_ratio[i]=0.;
            hi = 0.;
        }

        //if either c or h equal zero, we set the others to zero as well
        if ( M_conc[i] < 0.01 || hi < physical::hmin )
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
            M_ridge_ratio[i]=0.;
            hi = 0.;
        }
        else
        {
            // convert from absolute to effective thickness
            M_ridge_ratio[i]=M_ridge_ratio[i]*M_conc[i];
            M_snow_thick[i] = M_conc[i]*hs;
            M_thick[i] = M_conc[i]*hi;
        }

        // skip application of NIC if its conc is >1
        // - then it is masked
        // - unfortunately, applying the mask in datasets led to masked values being treated as real values of 0.0
        //   so we have to do it manually here
        if(M_nic_conc[i]<=1.)
        {
            // Use the NIC ice charts
            // - get conc bins from NIC dataset
            double thin_conc_obs = 0.;
            double thin_conc_obs_min = 0.;
            double thin_conc_obs_max = 0.;
            this->concBinsNic(thin_conc_obs_min, thin_conc_obs_max, M_nic_conc[i], use_weekly_nic);

            if((M_amsr2_conc[i]>=thin_conc_obs_min) && (M_amsr2_conc[i]<=thin_conc_obs_max))
            {
                thin_conc_obs_min=M_amsr2_conc[i];
                thin_conc_obs_max=M_amsr2_conc[i];
            }
            else
            {
                thin_conc_obs_min=0.5*(thin_conc_obs_min+thin_conc_obs_max);
                thin_conc_obs_max=thin_conc_obs_min;
            }

            if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                M_conc_thin[i]=0.;

                thin_conc_obs = thin_conc_obs_min-M_conc[i];
                if(thin_conc_obs>=0.)
                {
                    //if(thin_conc_obs>M_conc_thin[i])
                    //    M_h_thin[i] = M_h_thin[i]+(h_thin_min + (h_thin_max/2.-h_thin_min)*0.5)*(thin_conc_obs-M_conc_thin[i]);
                    //else
                    //    M_h_thin[i] = M_h_thin[i]*thin_conc_obs/M_conc_thin[i];

                    M_conc_thin[i] = thin_conc_obs;
                    M_h_thin[i] = (h_thin_min + (h_thin_max/2.-h_thin_min)*0.5)*M_conc_thin[i];
                }
                else
                {
                    M_conc_thin[i]=0.;
                    M_h_thin[i]=0.;

                    M_conc[i]=M_conc[i]+thin_conc_obs;
                    M_thick[i]=hi*M_conc[i];
                }

            }//thin ice
            else
            {
                if(M_conc[i]<thin_conc_obs_min)
                {
                    M_thick[i] = M_thick[i] + std::max(hi,0.5)*(thin_conc_obs_min-M_conc[i]); // 50 cm minimum for the added ice
                    M_conc[i] = thin_conc_obs_min;
                }
                else if(M_conc[i]>thin_conc_obs_max)
                {
                    M_conc[i] = thin_conc_obs_max;
                    M_thick[i]=hi*M_conc[i];
                }
            }//no thin ice
        }//use NIC
    }//loop over elements
}//topazForecastAmsr2OsisafNicIce

    
// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from PIOMAS outputs.
//! Called by the initIce() function.
void
FiniteElement::piomasIce()
{
    external_data M_init_conc=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,0,false,time_init);
    external_data M_init_thick=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,1,false,time_init);
    external_data M_init_snow_thick=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,2,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_init_conc);
    external_data_tmp.push_back(&M_init_thick);
    external_data_tmp.push_back(&M_init_snow_thick);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - PIOMAS ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    for (int i=0; i<M_num_elements; ++i)
    {
        M_conc[i] = std::min(1.,M_init_conc[i]);
        M_thick[i] = M_init_thick[i];
        M_snow_thick[i] = M_init_snow_thick[i];

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

        M_damage[i]=0.;
    }
}//piomasIce
    
    
// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from Topaz, AMSRE data.
//! Called by the initIce() function.
void
FiniteElement::topazAmsreIce()
{
    double real_thickness, init_conc_tmp;

    //obs
    external_data M_conc_amsre=ExternalData(&M_ice_amsre_elements_dataset,M_mesh,0,false,time_init);

    //topaz
    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_conc_amsre);
    external_data_tmp.push_back(&M_init_conc);
    external_data_tmp.push_back(&M_init_thick);
    external_data_tmp.push_back(&M_init_snow_thick);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - TOPAZ/AMSR-E ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
        M_conc[i] = std::min(1.,M_conc_amsre[i]);

        // TOPAZ puts very small values instead of 0.
        tmp_var=M_init_conc[i];
        init_conc_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
        tmp_var=M_init_thick[i];
        M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
        tmp_var=M_init_snow_thick[i];
        M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

        // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
        if(init_conc_tmp>0.05)
        {
            //convert from absolute to effective thickness (with AMSRE conc)
            M_thick[i] *= M_conc[i];
            M_snow_thick[i] *= M_conc[i];
        }

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_thick[i]=0.;
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

        M_damage[i]=0.;
    }
}//topazAmsreIce TODO no thin ice; logic needs checking


// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from Topaz, AMSR2 data.
//! Called by the initIce() function.
void
FiniteElement::topazAmsr2Ice()
{
    double real_thickness, init_conc_tmp;

    external_data M_conc_amsr2=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_conc_amsr2);
    external_data_tmp.push_back(&M_init_conc);
    external_data_tmp.push_back(&M_init_thick);
    external_data_tmp.push_back(&M_init_snow_thick);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - TOPAZ/AMSR2 ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
        double uncertainty;
        if(M_conc_amsr2[i]<0.1)
            uncertainty=0.1;
        else
            uncertainty=0.05;

        double diff_mod_obs = M_conc_amsr2[i]-M_init_conc[i];
        if(std::abs(diff_mod_obs)>=uncertainty && M_conc_amsr2[i]<=1.)
            // NB missing value for AMSR2 when not over land is 1.15
            // move towards AMSR2 value by the amount uncertainty/2
            M_conc[i] = std::min(1., M_conc_amsr2[i]-(diff_mod_obs/std::abs(diff_mod_obs))*uncertainty/2.);
        else
            M_conc[i] = std::min(1., M_init_conc[i]);

        // TOPAZ puts very small values instead of 0.
        tmp_var=M_init_conc[i];
        init_conc_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
        tmp_var=M_init_thick[i];
        M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
        tmp_var=M_init_snow_thick[i];
        M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

        // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
        if(init_conc_tmp>0.05)
        {
            //convert from absolute to effective thickness
            M_thick[i]*=M_conc[i];
            M_snow_thick[i]*=M_conc[i];
        }

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_thick[i]=0.;
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

        M_damage[i]=0.;
    }
}//topazAmsr2Ice TODO no thin ice; logic needs checking

    
// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from CS2 SMOS data.
//! Called by the initIce() function.
void
FiniteElement::cs2SmosIce()
{
    external_data M_init_conc=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,0,false,time_init);
    external_data M_init_thick=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,1,false,time_init);
    external_data M_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_init_conc);
    external_data_tmp.push_back(&M_init_thick);
    external_data_tmp.push_back(&M_type);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - CS2/SMOS ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    this->warrenClimatology();

    double tmp_var, correction_factor_warren;
    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_var=std::min(1.,M_init_conc[i]);
        M_conc[i] = tmp_var;
        tmp_var=M_init_thick[i];
        M_thick[i] = tmp_var ;

        double ratio_FYI=0.3;
        double ratio_MYI=0.9;
        double ratio_Mixed=0.6;

        if((M_thick[i]>0.)&&(M_conc[i])>0.2)
        {
            if(M_type[i]<=1.)
                M_ridge_ratio[i]=0.;
            if(M_type[i]>1. && M_type[i]<=2.)
                M_ridge_ratio[i]=(M_type[i]-1.)*ratio_FYI;
            if(M_type[i]>2. && M_type[i]<=3.)
                M_ridge_ratio[i]=(1.-(M_type[i]-2.))*ratio_FYI + (M_type[i]-2.)*ratio_MYI;
            if(M_type[i]>3. && M_type[i]<=4.)
                M_ridge_ratio[i]=(1.-(M_type[i]-3.))*ratio_MYI + (M_type[i]-3.)*ratio_Mixed;
            if(M_type[i]>4.)
                M_ridge_ratio[i]=ratio_Mixed;

            M_ridge_ratio[i]=M_ridge_ratio[i]*std::exp(ridging_exponent*(1.-M_conc[i]));
        }
        else
            M_ridge_ratio[i]=0.;


        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

        // Correction of the value given by Warren as a function of the ice type
        //M_type[i]==1. // No ice
        //M_type[i]==2. // First-Year ice
        //M_type[i]==3. // Multi-Year ice
        //M_type[i]==4. // Mixed
        correction_factor_warren=std::max(0.,std::min(1.,(M_type[i]-1.)*0.5)); // == 1. for MY, and mixed, 0.5 for FY, 0. for No ice

        M_snow_thick[i]=correction_factor_warren*M_snow_thick[i]*M_conc[i];

        M_damage[i]=0.;

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            //M_conc_thin[i]=std::max(M_conc_amsre[i]-M_conc[i],0.);
            M_conc_thin[i]=std::min(1.-M_conc[i], 0.2*M_conc[i]);
            M_h_thin[i]=M_conc_thin[i]*(h_thin_min+0.5*(h_thin_max-h_thin_min));
        }

        // Check that the snow is not so thick that the ice is flooded
        double max_snow = M_thick[i]*(physical::rhow-physical::rhoi)/physical::rhos;
        M_snow_thick[i] = std::min(max_snow, M_snow_thick[i]);

    }
}//cs2SmosIce

    
// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from CS2, SMOS, AMSR2 data.
//! Called by the initIce() function.
void
FiniteElement::cs2SmosAmsr2Ice()
{
    external_data M_init_conc=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,0,false,time_init);
    external_data M_amsr2_conc=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
    external_data M_init_thick=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,1,false,time_init);
    external_data M_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_init_conc);
    external_data_tmp.push_back(&M_init_thick);
    external_data_tmp.push_back(&M_type);
    external_data_tmp.push_back(&M_amsr2_conc);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - CS2/SMOS/AMSR2 ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    this->warrenClimatology();

    double tmp_var, correction_factor_warren;
    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_var=std::min(1.,M_init_conc[i]);
        M_conc[i] = tmp_var;
        tmp_var=M_init_thick[i];
        M_thick[i] = tmp_var ;
        if(M_amsr2_conc[i]<M_conc[i])
            // AMSR2 is higher resolution and see small opening that would not be see in OSISAF
            // NB AMSR2 = 1.15 if missing data over ocean, however this will not affect the example here
            M_conc[i]=M_amsr2_conc[i];

        double ratio_FYI=0.3;
        double ratio_MYI=0.9;
        double ratio_Mixed=0.6;

        if((M_thick[i]>0.)&&(M_conc[i])>0.2)
        {
            if(M_type[i]<=1.)
                M_ridge_ratio[i]=0.;
            if(M_type[i]>1. && M_type[i]<=2.)
                M_ridge_ratio[i]=(M_type[i]-1.)*ratio_FYI;
            if(M_type[i]>2. && M_type[i]<=3.)
                M_ridge_ratio[i]=(1.-(M_type[i]-2.))*ratio_FYI + (M_type[i]-2.)*ratio_MYI;
            if(M_type[i]>3. && M_type[i]<=4.)
                M_ridge_ratio[i]=(1.-(M_type[i]-3.))*ratio_MYI + (M_type[i]-3.)*ratio_Mixed;
            if(M_type[i]>4.)
                M_ridge_ratio[i]=ratio_Mixed;

            M_ridge_ratio[i]=M_ridge_ratio[i]*std::exp(ridging_exponent*(1.-M_conc[i]));
        }
        else
            M_ridge_ratio[i]=0.;


        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

        // Correction of the value given by Warren as a function of the ice type
        //M_type[i]==1. // No ice
        //M_type[i]==2. // First-Year ice
        //M_type[i]==3. // Multi-Year ice
        //M_type[i]==4. // Mixed
        correction_factor_warren=std::max(0.,std::min(1.,(M_type[i]-1.)*0.5)); // == 1. for MY, and mixed, 0.5 for FY, 0. for No ice

        M_snow_thick[i]=correction_factor_warren*M_snow_thick[i]*M_conc[i];

        M_damage[i]=0.;

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            //M_conc_thin[i]=std::max(M_conc_amsre[i]-M_conc[i],0.);
            M_conc_thin[i]=std::min(1.-M_conc[i], 0.2*M_conc[i]);
            M_h_thin[i]=M_conc_thin[i]*(h_thin_min+0.5*(h_thin_max-h_thin_min));
        }

        // Check that the snow is not so thick that the ice is flooded
        double max_snow = M_thick[i]*(physical::rhow-physical::rhoi)/physical::rhos;
        M_snow_thick[i] = std::min(max_snow, M_snow_thick[i]);

    }
}//cs2SmosAmsr2Ice
    
    
// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from SMOS data.
//! Called by the initIce() function.
void
FiniteElement::smosIce()
{
    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    external_data M_init_thick=ExternalData(&M_ice_smos_elements_dataset,M_mesh,0,false,time_init);

    boost::gregorian::date dt = Nextsim::parse_date(time_init);
    int month_id=dt.month().as_number(); // 1 for January, 2 for February, and so on. This will be used to compute the snow from Warren climatology

    std::cout << "month_id: " << month_id <<"\n";

    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_init_conc);
    external_data_tmp.push_back(&M_init_thick);
    external_data_tmp.push_back(&M_init_snow_thick);

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    if(M_rank==0)
        LOG(DEBUG)<<"init - SMOS ExternalData objects\n";
    this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_var=std::min(1.,M_init_conc[i]);
        M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
        tmp_var=M_init_thick[i];
        M_thick[i] = tmp_var ;
        tmp_var=M_init_snow_thick[i];
        M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

        M_damage[i]=0.;
    }
}//smosIce

    
// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from CS2, SMOS, AMSR2 data.
//! Called by the cs2SmosIce() and cs2SmosAmsr2Ice() functions.
void
FiniteElement::warrenClimatology()
{
    // Coefficients for the fit from Warren et al '99
    std::vector<std::vector<double>> coeffs(12);

    coeffs[0].push_back(28.01);
    coeffs[0].push_back( 0.1270);
    coeffs[0].push_back(-1.1833);
    coeffs[0].push_back(-0.1164);
    coeffs[0].push_back(-0.0051);
    coeffs[0].push_back( 0.0243);

    coeffs[1].push_back(30.28);
    coeffs[1].push_back( 0.1056);
    coeffs[1].push_back(-0.5908);
    coeffs[1].push_back(-0.0263);
    coeffs[1].push_back(-0.0049);
    coeffs[1].push_back( 0.0044);

    coeffs[2].push_back(33.89);
    coeffs[2].push_back( 0.5486);
    coeffs[2].push_back(-0.1996);
    coeffs[2].push_back( 0.0280);
    coeffs[2].push_back( 0.0216);
    coeffs[2].push_back(-0.0176);

    coeffs[3].push_back(36.80);
    coeffs[3].push_back( 0.4046);
    coeffs[3].push_back(-0.4005);
    coeffs[3].push_back( 0.0256);
    coeffs[3].push_back( 0.0024);
    coeffs[3].push_back(-0.0641);

    coeffs[4].push_back(36.93);
    coeffs[4].push_back( 0.0214);
    coeffs[4].push_back(-1.1795);
    coeffs[4].push_back(-0.1076);
    coeffs[4].push_back(-0.0244);
    coeffs[4].push_back(-0.0142);

    coeffs[5].push_back(36.59);
    coeffs[5].push_back( 0.7021);
    coeffs[5].push_back(-1.4819);
    coeffs[5].push_back(-0.1195);
    coeffs[5].push_back(-0.0009);
    coeffs[5].push_back(-0.0603);

    coeffs[6].push_back(11.02);
    coeffs[6].push_back( 0.3008);
    coeffs[6].push_back(-1.2591);
    coeffs[6].push_back(-0.0811);
    coeffs[6].push_back(-0.0043);
    coeffs[6].push_back(-0.0959);

    coeffs[7].push_back( 4.64);
    coeffs[7].push_back( 0.3100);
    coeffs[7].push_back(-0.6350);
    coeffs[7].push_back(-0.0655);
    coeffs[7].push_back( 0.0059);
    coeffs[7].push_back(-0.0005);

    coeffs[8].push_back(15.81);
    coeffs[8].push_back( 0.2119);
    coeffs[8].push_back(-1.0292);
    coeffs[8].push_back(-0.0868);
    coeffs[8].push_back(-0.0177);
    coeffs[8].push_back(-0.0723);

    coeffs[9].push_back(22.66);
    coeffs[9].push_back( 0.3594);
    coeffs[9].push_back(-1.3483);
    coeffs[9].push_back(-0.1063);
    coeffs[9].push_back( 0.0051);
    coeffs[9].push_back(-0.0577);

    coeffs[10].push_back(25.57);
    coeffs[10].push_back( 0.1496);
    coeffs[10].push_back(-1.4643);
    coeffs[10].push_back(-0.1409);
    coeffs[10].push_back(-0.0079);
    coeffs[10].push_back(-0.0258);

    coeffs[11].push_back(26.67);
    coeffs[11].push_back(-0.1876);
    coeffs[11].push_back(-1.4229);
    coeffs[11].push_back(-0.1413);
    coeffs[11].push_back(-0.0316);
    coeffs[11].push_back(-0.0029);

    // Time interpolation
    boost::gregorian::date t = Nextsim::parse_date(time_init);
    int month  = t.month().as_number(); // 1 for January, 2 for February, and so on.
    int day    = t.day().as_number(); // Day of the month.
    int eomday = t.end_of_month().day().as_number(); // Last day of the month.

    int month2; // The other month to interpolate from (one before or after)
    double dt;
    if ( day < eomday/2. )
    { // We're in the early part of the month and interpolate to the previous month
        month2 = month-1;
        dt     = eomday/2+day;
        if(month2==0)
            month2=12;
    }
    else
    { // We're in the late part of the month and interpolate to the next month
        month2 = month+1;
        dt     = eomday/2 + eomday - day;
        if(month2==13)
            month2=1;
    }

    // Now calculate snow thickness for the current day as an inexact temporal interpolation
    // I just assume all months are as long as the current one ... but the error incurred is small
    std::vector<double> lon = M_mesh.meanLon();
    std::vector<double> lat = M_mesh.meanLat();

    for ( int i=0; i<M_num_elements; i++ )
    {
        const double pi = std::acos(-1);

        double x = (90 - lat[i]) * std::cos(lon[i]*pi/180.);
        double y = (90 - lat[i]) * std::sin(lon[i]*pi/180.);

        M_snow_thick[i] = 1e-2*dt/eomday*std::max(0.,
                coeffs[month-1][0] + coeffs[month-1][1]*x + coeffs[month-1][2]*y + coeffs[month-1][3]*x*y + coeffs[month-1][4]*x*x + coeffs[month-1][5]*y*y)
            + 1e-2*(eomday-dt)/eomday*std::max(0.,
                coeffs[month2-1][0] + coeffs[month2-1][1]*x + coeffs[month2-1][2]*y + coeffs[month2-1][3]*x*y + coeffs[month2-1][4]*x*x + coeffs[month2-1][5]*y*y);
    }

}//warrenClimatology


// -----------------------------------------------------------------------------------------------------------
//! Initializes the type of drifters used (equally spaced, IABP, RGPS, OSISAF, SIDFEX).
//  NB needs to be done before readRestart()
//! Called by the init() function.
void
FiniteElement::initDrifterOpts()
{
    // init drifters after spinup
    // NB use ceil to make sure the init time is 0:00
    M_drifters_time_init = std::ceil(time_init
        + vm["simul.spinup_duration"].as<double>());

    //to be run at start time
    std::vector<double> drifters_timesteps;
    std::vector<std::string> drifters_names;//for debugging

    // use IABP drifters (12h check)
    M_use_iabp_drifters = vm["drifters.use_iabp_drifters"].as<bool>();
    M_iabp_drifters_output_time_step = vm["drifters.iabp_drifters_output_time_step"].as<double>();
    M_iabp_drifters_input_time_step = 0.5;
    if (M_use_iabp_drifters)
    {
        // check the IABP in-/out-put time steps are consistent with each other
        // - IABP input time step should be greater than the output time step
        //   (same frequency of outputs or more frequent outputs),
        // - IABP input time step should be a multiple of the output time step
        if( M_iabp_drifters_output_time_step > M_iabp_drifters_input_time_step )
        {
            std::string msg = "IABP drifters output timestep";
            msg += " should be <= IABP input timestep";
            throw std::runtime_error(msg);
        }
        else if ( std::fmod(M_iabp_drifters_input_time_step,
                    M_iabp_drifters_output_time_step) != 0 )
        {
            throw std::runtime_error(
                    "IABP drifter input timestep should be a multiple of the IABP output timestep");
        }

        // to check all the drifter timesteps later
        drifters_timesteps.push_back(M_iabp_drifters_input_time_step);
        drifters_names.push_back("IABP (input)");
        drifters_timesteps.push_back(M_iabp_drifters_output_time_step);
        drifters_names.push_back("IABP (output)");
    }

    // use OSISAF drifters
    // - every day at 12:00 start a new set of drifters which run for 48h
    // - output time step is an option
    M_use_osisaf_drifters = vm["drifters.use_osisaf_drifters"].as<bool>();
    M_osisaf_drifters_output_time_step = vm["drifters.osisaf_drifters_output_time_step"].as<double>();
    if(M_use_osisaf_drifters)
    {
        M_osisaf_drifters.resize(2);

        if( M_osisaf_drifters_output_time_step > 2.0 )
        {
            // currently OSISAF output time step has to fit into 2 days
            throw std::runtime_error("OSISAF drifters output timestep should be <= 2.0");
        }
        else if( fmod(2.0, M_osisaf_drifters_output_time_step) != 0 )
        {
            // output timestep should fit into .5
            throw std::runtime_error("2.0 should be a multiple of the OSISAF drifters output timestep");
        }

        // to check all the drifter timesteps later
        drifters_timesteps.push_back(M_osisaf_drifters_output_time_step);
        drifters_names.push_back("OSISAF");

        // add drifters to the list of ordinary drifters
        M_ordinary_drifters.push_back(&M_osisaf_drifters[0]);
        M_ordinary_drifters.push_back(&M_osisaf_drifters[1]);
    }

    // equally spaced drifters
    M_use_equally_spaced_drifters = vm["drifters.use_equally_spaced_drifters"].as<bool>();
    M_equally_spaced_drifters_output_time_step = vm["drifters.equally_spaced_drifters_output_time_step"].as<double>();
    if (M_use_equally_spaced_drifters)
    {
        // to check all the drifter timesteps later
        drifters_timesteps.push_back(M_equally_spaced_drifters_output_time_step);
        drifters_names.push_back("Equally-spaced");

        // add drifter to the list of ordinary drifters
        M_ordinary_drifters.push_back(&M_equally_spaced_drifters);
    }

    // RGPS drifters
    M_use_rgps_drifters = vm["drifters.use_rgps_drifters"].as<bool>();
    M_rgps_drifters_output_time_step = vm["drifters.rgps_drifters_output_time_step"].as<double>();
    if (M_use_rgps_drifters)
    {
        std::string time_str = vm["drifters.RGPS_time_init"].as<std::string>();
        M_rgps_time_init = Nextsim::from_date_time_string(time_str);
        M_rgps_file = Environment::nextsimDataDir().string()
            + "/RGPS_" + time_str + ".txt";

        // make sure RGPS init time is not too early
        if(M_rgps_time_init<M_drifters_time_init)
        {
            std::string msg = "";
            msg += "M_rgps_time_init should be >= M_drifters_time_init\n";
            msg += "(ie after simulation time_init";
            msg += "(rounded up to the nearest day) and spinup period)";
            throw std::runtime_error(msg);
        }

        // to check all the drifter timesteps later
        drifters_timesteps.push_back(M_rgps_drifters_output_time_step);
        drifters_names.push_back("RGPS");

        // add drifter to the list of ordinary drifters
        M_ordinary_drifters.push_back(&M_rgps_drifters);
    }

    // SIDFEX drifters
    M_use_sidfex_drifters = vm["drifters.use_sidfex_drifters"].as<bool>();
    M_sidfex_drifters_output_time_step = vm["drifters.sidfex_drifters_output_time_step"].as<double>();
    if (M_use_sidfex_drifters)
    {
        drifters_timesteps.push_back(M_sidfex_drifters_output_time_step);
        drifters_names.push_back("SIDFEX");

        // add drifter to the list of ordinary drifters
        M_ordinary_drifters.push_back(&M_sidfex_drifters);
    }

    // can now tell if we are using any drifters and set the init time
    // (after spinup)
    M_use_drifters = ( M_use_iabp_drifters || M_ordinary_drifters.size()>0 );

    if(M_use_drifters)
    {
        // check consistency of drifter output time steps with model time step
        for(int i=0; i<drifters_timesteps.size(); i++)
        {
            if( fmod(drifters_timesteps[i]*24*3600, time_step) != 0 )
            {
                std::string msg = drifters_names[i]+" drifters' timestep not a multiple of model time step";
                throw std::runtime_error(msg);
            }
        }
    }
}//initDrifterOpts


// -----------------------------------------------------------------------------------------------------------
//! Initialise the drifters if it's the right time
//! * Input init_names  is set by initialisingDrifters() 
//!   - determines which (if any) need to be initialised this time step
//! Called by the checkDrifters() function.
void
FiniteElement::initDrifters(mesh_type_root const& movedmesh_root,
        std::vector<std::string> const& init_names)
{
    //! -1) do we initialise all the drifters that start at the usual time?
    //! * Equally-spaced, SIDFEX, IABP
    if(std::count(init_names.begin(), init_names.end(), "main")>0)
    {
        if(M_use_equally_spaced_drifters)
            this->initEquallySpacedDrifters(movedmesh_root);

        if(M_use_sidfex_drifters)
            this->initSidfexDrifters(movedmesh_root);

        if(M_use_iabp_drifters)
        {
            if(M_iabp_drifters.size()==0)
                //only init if not already initialised (eg from restart)
                this->initIabpDrifters(movedmesh_root);
        }
    }

    //! -2) initialise the RGPS drifters that start at their own time?
    if(std::count(init_names.begin(), init_names.end(), "rgps")>0)
        this->initRGPSDrifters(movedmesh_root);

    if(std::count(init_names.begin(), init_names.end(), "osisaf")>0)
        //! start a new set of OSISAF drifters?
        // NB do this after outputting netcdf,
        // to make sure the last time record is present
        // for the drifters that will be terminated
        this->initOsisafDrifters(movedmesh_root);
    
}//initDrifters


// -----------------------------------------------------------------------------------------------------------
//! Initialise the IABP drifters
//! Called by the initDrifters() function.
void
FiniteElement::initIabpDrifters(mesh_type_root const& movedmesh_root)
{
    // init the input/output files
    this->initIabpDrifterFiles();

    // Get:
    // - the 1st drifter positions (if any)
    // - M_conc at these positions
    this->updateIabpDrifters(movedmesh_root);

    // Save the initial positions to the output file
    this->outputIabpDrifters();
}//initIabpDrifters


// -----------------------------------------------------------------------------------------------------------
//! Initialise the IABP drifter input and output text files
//! Called by the initIabpDrifters() and readRestart() functions.
void
FiniteElement::initIabpDrifterFiles()
{
    // OUTPUT:
    // We should tag the file name with the init time in case of a re-start.
    std::stringstream filename_out;
    filename_out << M_export_path << "/IABP_drifters_simulated_"
        << to_date_time_string_for_filename(M_current_time) << ".txt";
    M_iabp_outfile = filename_out.str();
    std::fstream iabp_out(M_iabp_outfile, std::fstream::out );
    if ( ! iabp_out.good() )
        throw std::runtime_error("Cannot write to file: " + M_iabp_outfile);

    //write the header and close
    iabp_out << "Year Month Day Hour BuoyID Lat Lon Concentration\n";
    iabp_out.close();

    // INPUT:
    //new buoy file has a header
    std::string filename_in = Environment::nextsimDataDir().string()
        + "/IABP_drifters.txt";
    M_iabp_infile_fstream.open(filename_in, std::fstream::in);
    if ( ! M_iabp_infile_fstream.good() )
        throw std::runtime_error("File not found: " + filename_in);

    //skip header
    std::string header;
    std::getline(M_iabp_infile_fstream, header);
    std::cout<<"open IABP drifter file: "<<filename_in<<"\n";
    std::cout<<"header: "<<header<<"\n";

    int pos;    // To be able to rewind one line
    double time = from_date_string("1979-01-01");
    while ( time < M_drifters_time_init )
    {
        // Remember where we were
        pos = M_iabp_infile_fstream.tellg();

        // Read the next line
        int year, month, day, hour, number;
        double lat, lon;
        M_iabp_infile_fstream >> year >> month >> day >> hour >> number >> lat >> lon;
        std::string date = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);

        time = from_date_string(date) + hour/24.;
    }

    // We must rewind one line so that updateIabpDrifters works correctly
    M_iabp_infile_fstream.seekg(pos);
}//initIabpDrifterFiles


// -----------------------------------------------------------------------------------------------------------
//! Outputs the IABP drifter positions and conc
//! Called by the readRestart(), initDrifters() and checkDrifters() functions.
void
FiniteElement::outputIabpDrifters()
{
    if (M_rank==0)
    {
        if (M_iabp_drifters.size()==0)
            //do nothing if drifters are empty
            return;

        // Initialize the map
        mapx_class *map;
        std::string configfile = Environment::nextsimMppfile();
        std::vector<char> str(configfile.begin(), configfile.end());
        str.push_back('\0');
        map = init_mapx(&str[0]);

        //open output file for appending
        std::fstream iabp_out(M_iabp_outfile, std::fstream::out | std::fstream::app );
        if ( ! iabp_out.good() )
            throw std::runtime_error("Cannot write to file: " + M_iabp_outfile);

        // Loop over the map and output
        int j=0;
        boost::gregorian::date           date = Nextsim::parse_date( M_current_time );
        boost::posix_time::time_duration time = Nextsim::parse_time( M_current_time );
        for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); ++it )
        {
            double lat, lon;
            double conc = M_iabp_conc[j];
            inverse_mapx(map, it->second[0], it->second[1], &lat, &lon);
            j++;

            iabp_out
                << std::setw(4) << date.year()
                << " " << std::setw( 2) << date.month().as_number()
                << " " << std::setw( 2) << date.day().as_number()
                << " " << std::setw( 2) << time.hours()
                << " " << std::setw(16) << it->first
                << std::fixed << std::setprecision(5)
                << " " << std::setw( 8) << lat
                << " " << std::setw(10) << lon
                << " " << conc
                << "\n";
        }

        //close output file
        iabp_out.close();

        // close the map
        close_mapx(map);
    }
}//outputIabpDrifters


// -----------------------------------------------------------------------------------------------------------
//! Updates the IABP buoy in time, by adding the buoys that have been put into the ice and removing the dead ones.
//! * Also updates the conc at the IABP drifters
//! Called by the updateDrifterPosition() function.
void
FiniteElement::updateIabpDrifters(mesh_type_root const& movedmesh_root)
{

    if(M_rank==0)
    {
        // Initialize the map
        mapx_class *map;
        std::string configfile = Environment::nextsimMppfile();
        std::vector<char> str(configfile.begin(), configfile.end());
        str.push_back('\0');
        map = init_mapx(&str[0]);

        // Read the current buoys from file
        double time = M_current_time;
        std::vector<int> keepers;
        while ( time == M_current_time )
        {
            // Read the next line
            int year, month, day, hour, number;
            double lat, lon;
            M_iabp_infile_fstream >> year >> month >> day >> hour >> number >> lat >> lon;
            std::string date = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);
            time = from_date_string(date) + hour/24.;

            // Remember which buoys are in the ice according to IABP
            keepers.push_back(number);

            // Project and add the buoy to the map if it's missing
            if ( M_iabp_drifters.count(number) == 0 )
            {
                double x, y;
                forward_mapx(map,lat,lon,&x,&y);
                M_iabp_drifters.emplace(number, std::array<double,2>{x, y});
                LOG(DEBUG)<<"new IABP buoy: time, index, lon, lat, x, y: "
                    <<to_date_time_string(time)<<", "<<number<<", "
                    <<lon<<", "<<lat<<", "<<x<<", "<<y<<"\n";
            }
        }
        close_mapx(map);

        // update M_iabp_conc (get model conc at all drifters)
        this->updateIabpDrifterConc(movedmesh_root);

        // Rebuild the M_iabp_drifters map
        double clim = vm["drifters.concentration_limit"].as<double>();

        // Go through the M_iabp_drifters map and throw out:
        // (i) the ones which IABP doesn't report as being in the ice anymore
        // (ii) the ones which have a low conc according to the model
        int j =0;
        auto model_conc = M_iabp_conc;
        M_iabp_conc.resize(0);
        for ( auto model = M_iabp_drifters.begin(); model != M_iabp_drifters.end(); j++)// NB ++model is not allowed here, because we use 'erase'
        {
            double conc = model_conc[j];
            bool keep = (conc > clim)                                               // in ice (model)
                && ( std::count(keepers.begin(), keepers.end(), model->first) >0 ); // in ice (IABP)

            // Delete or advance the iterator
            if ( ! keep )
                model = M_iabp_drifters.erase(model);
            else
            {
                ++model;
                M_iabp_conc.push_back(conc);
            }
        }//remove drifters not in ice according to IABP or model
    }
}//updateIabpDrifters


// -----------------------------------------------------------------------------------------------------------
//! Updates the concentration variable at the IABP buoys
//! Called by the updateIabpDrifters() and readRestart() function.
void
FiniteElement::updateIabpDrifterConc(mesh_type_root const& movedmesh_root)
{
    if( M_rank ==0 )
    {
        // Assemble the coordinates from the unordered_map
        int Ndrifters = M_iabp_drifters.size();
        std::vector<double> drifter_X(Ndrifters);
        std::vector<double> drifter_Y(Ndrifters);
        int j=0;
        for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); ++it )
        {
            drifter_X[j] = it->second[0];
            drifter_Y[j] = it->second[1];
            ++j;
        }

        // Interpolate the concentration
        double* interp_drifter_c_out;
        InterpFromMeshToMesh2dx(&interp_drifter_c_out,
                                &movedmesh_root.indexTr()[0], &movedmesh_root.coordX()[0], &movedmesh_root.coordY()[0],
                                movedmesh_root.numNodes(), movedmesh_root.numTriangles(),
                                &M_conc_root[0],
                                movedmesh_root.numTriangles(), 1,
                                &drifter_X[0], &drifter_Y[0],
                                Ndrifters,
                                true, 0.);

        // Finally retrieve M_iabp_conc and delete the pointer
        M_iabp_conc.resize(M_iabp_drifters.size());
        for( int j=0; j<M_iabp_drifters.size(); j++ )
            M_iabp_conc[j] = interp_drifter_c_out[j];
        xDelete<double>(interp_drifter_c_out);
    }

}//updateIabpDrifterConc


// -----------------------------------------------------------------------------------------------------------
//! Updates the positions of the IABP buoys
//! Called by the updateIabpDrifters() function.
void
FiniteElement::updateIabpDrifterPosition()
{
    if(M_rank==0)
    {
        chrono.restart();
        LOG(DEBUG) <<"IABP Drifter moving starts\n";

        // Assemble the coordinates from the unordered_map
        std::vector<double> drifter_X(M_iabp_drifters.size());
        std::vector<double> drifter_Y(M_iabp_drifters.size());
        int j=0;
        for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); ++it )
        {
            drifter_X[j] = it->second[0];
            drifter_Y[j] = it->second[1];
            ++j;
        }

        // Interpolate the total displacement onto the drifter positions
        int nb_var=2;
        std::vector<double> interp_drifter_in(nb_var*M_mesh_root.numNodes());
        for (int i=0; i<M_mesh_root.numNodes(); ++i)
        {
            interp_drifter_in[nb_var*i]   = M_UT_root[i];
            interp_drifter_in[nb_var*i+1] = M_UT_root[i+M_mesh_root.numNodes()];
        }

        double* interp_drifter_out;
        InterpFromMeshToMesh2dx(&interp_drifter_out,
                                &M_mesh_root.indexTr()[0], &M_mesh_root.coordX()[0], &M_mesh_root.coordY()[0],
                                M_mesh_root.numNodes(), M_mesh_root.numTriangles(),
                                &interp_drifter_in[0],
                                M_mesh_root.numNodes(), nb_var,
                                &drifter_X[0], &drifter_Y[0], M_iabp_drifters.size(),
                                true, 0.);

        // Rebuild the M_iabp_drifters map
        // (move the IABP drifters)
        j = 0;
        for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); ++it)
        {
            M_iabp_drifters[it->first] = std::array<double,2> {
                it->second[0]+interp_drifter_out[nb_var*j],
                it->second[1]+interp_drifter_out[nb_var*j+1]
            };
            ++j;
        }

        xDelete<double>(interp_drifter_out);

        LOG(DEBUG) <<"IABP Drifter move done in "<< chrono.elapsed() <<"s\n";
    }
}//updateIabpDrifterPosition


// -----------------------------------------------------------------------------------------------------------
//! This function tests which (if any) drifters we are initialising
//! at the current time step
//! It also calculates if we are doing any initialising at this timestep
//! Called by checkDrifters()
void
FiniteElement::initialisingDrifters(
        std::vector<std::string> &init_names,
        bool &init_any
        )
{
    init_names.resize(0);
    init_any = false;
    if(M_rank != 0)
        return;

    // inputting/outputting main drifters?
    // - those that start at the usual time
    //   - IABP
    //   - equally-spaced
    //   - SIDFEX
    if (       M_use_iabp_drifters
            || M_use_equally_spaced_drifters
            || M_use_sidfex_drifters
            )
        if( M_current_time == M_drifters_time_init)
            init_names.push_back("main");

    // initialising RGPS drifters?
    if(M_use_rgps_drifters)
        if( M_current_time == M_rgps_time_init )
            init_names.push_back("rgps");

    // initialising OSISAF drifters?
    // - is it 12:00 and past the spinup duration?
    if ( M_use_osisaf_drifters )
        if(M_current_time > M_drifters_time_init
                && std::fmod( M_current_time + 0.5, 1.0) == 0)
            init_names.push_back("osisaf");

    // calculate if we need to initialise ANY,
    // since then we need to ensure M_UT = 0
    init_any = (init_names.size()>0);
}//initialisingDrifters


// -----------------------------------------------------------------------------------------------------------
//! This function tests which (if any) drifters we are outputting (or inputting for IABP)
//! at the current time step
//! It also calculates if we are doing any outputting or inputting
//! Called by checkDrifters()
void
FiniteElement::outputtingDrifters(
        bool &input_iabp,
        bool &output_iabp,
        bool &io_any
        )
{

    input_iabp = false;
    output_iabp = false;
    io_any = false;
    if(M_rank != 0)
        return;

    // we only output after init time (not at init time)
    // - netcdf/text file outputs are already written at init time of drifters
    if(M_current_time>M_drifters_time_init)
    {
        // inputting/outputting IABP drifters?
        if ( M_use_iabp_drifters )
        {
            input_iabp = (std::fmod(
                    M_current_time - M_drifters_time_init,
                    M_iabp_drifters_input_time_step
                    )==0);
            output_iabp = (std::fmod(
                    M_current_time - M_drifters_time_init,
                    M_iabp_drifters_output_time_step
                    )==0);
        }
    }

    // Check the non-IABP drifters by looping over M_ordinary_drifters
    io_any = input_iabp || output_iabp;
    for(auto it = M_ordinary_drifters.begin(); it!= M_ordinary_drifters.end(); it++)
        io_any = io_any || (*it)->isOutputTime(M_current_time);

}//outputtingDrifters()


// -----------------------------------------------------------------------------------------------------------
//! Checks if we have to init, update, output or move any drifters
//! Called by the <FiniteElement::checkOutputs>() function.
void
FiniteElement::checkDrifters()
{
    // although these are only needed on the root processor,
    // we need these in several scopes, so declare them here
    bool io_any = false;                    // input/output any drifters (then we need to move)?
    bool input_iabp = false;                // need to update IABP?
    bool output_iabp = false;               // need to output IABP?
    bool init_any = false;                  //init any drifters (then we need to move)?
    std::vector<std::string> init_names;    // names of the drifters to init this time step
    bool move_drifters = false;             //needed on all processors (broadcast later)

    if(M_rank==0)
    {
        //! - 1) check if need to:
        //! * move the drifters (if needed)
        //! * get new inputs (IABP, OSISAF) (if needed)
        //! * output to text file or netcdf (if needed)

        // do we need to init any drifters?
        this->initialisingDrifters(
                init_names,
                init_any);

        // do we need to output any drifters?
        this->outputtingDrifters(
                input_iabp,
                output_iabp,
                io_any);

        //! - 2) Do we need to move any drifters?
        //! * If we are outputting ANY, we move ALL the drifters, since they all use the same
        //!   total displcement (M_UT)
        //! * Also if we are initialising any, we also need to move the drifters
        //!   - this is because we need M_UT=0 at init time,
        //!     so the next move is relative to their initial position
        move_drifters = (init_any || io_any);
    }

    // let all the processors know if we need to gather the three vectors
    // or if we need to reset M_UT
    boost::mpi::broadcast(M_comm, move_drifters, 0);

    if(move_drifters)
    {
        //! - 3) Gather the fields needed by the drifters if we are moving
        this->gatherNodalField(M_UT, M_UT_root);
        this->gatherNodalField(M_UM, M_UM_root);
        this->gatherElementField(M_conc, M_conc_root);

        // can now reset M_UT to 0
        std::fill(M_UT.begin(), M_UT.end(), 0.);

        if(M_rank==0)
        {
            //! - 4) Move the drifters
            // NB M_UT is relative to the fixed mesh, not the moved mesh
            // NB to update the conc we need to move the mesh
            //! * IABP drifters 
            if ( M_use_iabp_drifters )
                if(M_current_time>M_drifters_time_init)
                    this->updateIabpDrifterPosition();

            //! * other drifters 
            for(auto it=M_ordinary_drifters.begin(); it!=M_ordinary_drifters.end(); it++)
                if ((*it)->isInitialised())
                    (*it)->move(M_mesh_root, M_UT_root);
            
            // reset M_UT_root
            // TODO do we need to do this?
            // - it isn't used and seems to be gathered each time it's needed
            // (eg in writeRestart)
            // - needs checking though
            LOG(DEBUG)<<"Reset M_UT_root: "
                <<to_date_time_string(M_current_time)<<"\n";
            std::fill(M_UT_root.begin(), M_UT_root.end(), 0.);

            auto movedmesh_root = M_mesh_root;
            movedmesh_root.move(M_UM_root, 1.);

            //! - 5) Do we need to input/output any drifters?
            if (input_iabp)
                // check if we need to add new IABP drifters
                // NB do this after moving
                // NB this updates M_iabp_conc
                this->updateIabpDrifters(movedmesh_root);

            if (output_iabp)
            {
                // output IABP drifters
                // NB do this after moving
                if(!input_iabp)
                    //still need to update M_iabp_conc
                    this->updateIabpDrifterConc(movedmesh_root);
                this->outputIabpDrifters();
            }

            for(auto it=M_ordinary_drifters.begin(); it!=M_ordinary_drifters.end(); it++)
                if ((*it)->isOutputTime(M_current_time))
                {
                    (*it)->updateConc(movedmesh_root, M_conc_root);
                    (*it)->appendNetCDF(M_current_time);
                }
             
            //! - 6) Do we need to initialise any drifters?
            //  * Do initialisation after moving,
            //    to avoid moving drifters that are only just initialised
            //  * Also do it after output - this is mainly for the OSISAF drifters
            //    (to output the last record of M_osisaf_drifters[1], before it is overwritten)
            if(init_any)
                this->initDrifters(movedmesh_root, init_names);

        }//M_rank==0
    }//if(move_drifters)

}//checkDrifters


// -----------------------------------------------------------------------------------------------------------
//! Calculates the planetary vorticity (f)
//! Called by the step() function.
void
FiniteElement::calcCoriolis()
{
    // Interpolation of the latitude
    std::vector<double> lat = M_mesh.meanLat();

    for (int i=0; i<M_fcor.size(); ++i)
        M_fcor[i] = 2*(physical::omega)*std::sin(lat[i]*PI/180.);
}//calcCoriolis()


// -----------------------------------------------------------------------------------------------------------
//! Initializes the bathymetry according to the type used (constant or ETOPO).
//! Called by the step() function.
void
FiniteElement::initBathymetry()//(double const& u, double const& v)
{
    // This always needs to be done, regardless of which bathymetry type we
    // have as the object M_bathymetry_elements_dataset must be initialised. But
    // if we use CONSTANT then we don't put any data into the object.
    M_bathymetry_elements_dataset=DataSet("etopo_elements");

    switch (M_bathymetry_type)
    {
        case setup::BathymetryType::CONSTANT:
            M_element_depth=ExternalData(vm["ideal_simul.constant_bathymetry"].as<double>());
            M_external_data_elements.push_back(&M_element_depth);
            break;
        case setup::BathymetryType::ETOPO:
            M_element_depth=ExternalData(&M_bathymetry_elements_dataset,M_mesh,0,false,time_init);
            M_external_data_elements.push_back(&M_element_depth);
            break;
        default:
            std::cout << "invalid bathymetry"<<"\n";
            throw std::logic_error("invalid bathymetry");
    }
    M_datasets_regrid.push_back(&M_bathymetry_elements_dataset);//this needs to be reloaded if we are regridding
}//initBathymetry


void
FiniteElement::nodesToElements(double const* depth, std::vector<double>& v)
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


// -----------------------------------------------------------------------------------------------------------
//! Initializes the equally-spaced drifters, and saves the first netcdf time record
//! Called by the initDrifters() function.
void
FiniteElement::initEquallySpacedDrifters(mesh_type_root const& movedmesh_root)
{
    if (M_rank == 0)
    {
        M_equally_spaced_drifters = Drifters(1e3*vm["drifters.spacing"].as<double>(),
                movedmesh_root, M_conc_root, vm["drifters.concentration_limit"].as<double>(),
                M_current_time, M_equally_spaced_drifters_output_time_step);
        M_equally_spaced_drifters.initNetCDF(M_export_path+"/Drifters_", M_current_time);
        M_equally_spaced_drifters.appendNetCDF(M_current_time);
    }
}//initEquallySpacedDrifters


// -----------------------------------------------------------------------------------------------------------
//! Initializes the RGPS drifters, and saves the first netcdf time record
//! Called by the initDrifters() function.
void
FiniteElement::initRGPSDrifters(mesh_type_root const& movedmesh_root)
{
    if (M_rank == 0)
    {
        //called once when M_current_time == M_rgps_time_init
        M_rgps_drifters = Drifters(M_rgps_file, movedmesh_root, M_conc_root,
                -1, //assume that RGPS drifters' initial positions are OK and don't need masking due to low concentrations
                M_rgps_time_init, M_rgps_drifters_output_time_step);

        //save first outputs
        M_rgps_drifters.initNetCDF(M_export_path+"/RGPS_Drifters_", M_current_time);
        M_rgps_drifters.appendNetCDF(M_current_time);
    }
}//initRGPSDrifters
    

// -----------------------------------------------------------------------------------------------------------
//! Initializes the SIDFEX drifters, and saves the first netcdf time record
//! Called by the checkDrifters() function.
void
FiniteElement::initSidfexDrifters(mesh_type_root const& movedmesh_root)
{
    if (M_rank == 0)
    {
        std::string filename = Environment::nextsimDataDir().string() +"/"
            + vm["drifters.sidfex_filename"].as<std::string>();

        M_sidfex_drifters = Drifters(filename, movedmesh_root, M_conc_root,
                -1, //assume that SIDFEX drifters' initial positions are OK and don't need masking due to low concentrations
                M_current_time, M_sidfex_drifters_output_time_step);
        
        M_sidfex_drifters.initNetCDF(M_export_path+"/SIDFEx_Drifters_", M_current_time);
        M_sidfex_drifters.appendNetCDF(M_current_time);
    }
}//initSidfexDrifters


// -----------------------------------------------------------------------------------------------------------
//! Initializes a new set of OSISAF drifters, and saves the first netcdf time record
//! Called by the checkDrifters() function.
void
FiniteElement::initOsisafDrifters(mesh_type_root const& movedmesh_root)
{
    // OSISAF drift is calculated as a drifter displacement over 48 hours
    // - start a new one each day at 12:00
    // We have two sets of drifters in the field at all times
    // - [0] is the newest
    // Here we start a new set of drifters.

    if(M_rank != 0)
        return;

    LOG(DEBUG)<<"Initialize a new set of OSISAF drifters\n";

    // Flip the vector so we move [0] to be [1]
    std::reverse(M_osisaf_drifters.begin(), M_osisaf_drifters.end());

    // Create a new M_drifters instance in [0], with a properly initialised netCDF file
    std::string osi_grid_file = "ice_drift_nh_polstere-625_multi-oi.nc";
    std::string osi_output_path = M_export_path+"/OSISAF_";
    if(vm["drifters.use_refined_osisaf_grid"].as<bool>())
    {
        // use grid refined by a factor of 9
        // - can then compare averaged drift to the observations
        // - using an odd number in the refinement means the original grid points are a sub-sample of the refined grid
        osi_grid_file = "ice_drift_nh_polstere-625_multi-grid_refined_9.nc";
        osi_output_path = M_export_path+"/OSISAF_refined9_";
    }
        
    M_osisaf_drifters[0] = Drifters(osi_grid_file,
            "xc", "yc",
            "lat", "lon", movedmesh_root,
            M_conc_root, vm["drifters.concentration_limit"].as<double>(),
            M_current_time, M_osisaf_drifters_output_time_step);

    M_osisaf_drifters[0].initNetCDF(osi_output_path, M_current_time);
    M_osisaf_drifters[0].appendNetCDF(M_current_time);
}//initOsisafDrifters()

    
// -------------------------------------------------------------------------------------
//! Imports a BAMG mesh grid.
//! Called by the readRestart() and adaptMesh functions.
void
FiniteElement::importBamg(BamgMesh const* bamg_mesh)
{
    std::vector<point_type> mesh_nodes;
    std::vector<element_type> mesh_triangles;
    std::vector<double> coords(3,0);

    mesh_nodes.resize(bamg_mesh->VerticesSize[0]);

    for (int id=0; id<bamg_mesh->VerticesSize[0]; ++id)
    {
        coords[0] = bamg_mesh->Vertices[3*id];
        coords[1] = bamg_mesh->Vertices[3*id+1];

        mesh_nodes[id].id = id+1;
        mesh_nodes[id].coords = coords;
    }

    int type = 2;
    int physical = 0;
    int elementary = 0;

    int numVertices = 3;
    std::vector<int> indices(numVertices);

    mesh_triangles.resize(bamg_mesh->TrianglesSize[0]);

    for (int tr=0; tr<bamg_mesh->TrianglesSize[0]; ++tr)
    {
        indices[0] = bamg_mesh->Triangles[4*tr];
        indices[1] = bamg_mesh->Triangles[4*tr+1];
        indices[2] = bamg_mesh->Triangles[4*tr+2];

        element_type gmshElt( tr,
                              type,
                              physical,
                              elementary,
                              numVertices,
                              indices );

        mesh_triangles[tr] = gmshElt;
    }

    LOG(DEBUG) <<"\n";
    LOG(DEBUG) <<"Previous  NumNodes     = "<< M_mesh.numNodes() <<"\n";
    LOG(DEBUG) <<"Previous  NumTriangles = "<< M_mesh.numTriangles() <<"\n";

    M_mesh_previous_root = M_mesh_root;
    M_mesh_root.update(mesh_nodes, mesh_triangles);

    LOG(DEBUG) <<"\n";
    LOG(DEBUG) <<"Current  NumNodes      = "<< M_mesh_root.numNodes() <<"\n";
    LOG(DEBUG) <<"Current  NumTriangles  = "<< M_mesh_root.numTriangles() <<"\n";
    LOG(DEBUG) <<"\n";
}//importBamg

    
// -------------------------------------------------------------------------------------
//! Creates a new graphmpi_type object.
//! Called by the distributeMeshProcessing() function.
void
FiniteElement::createGraph()
{
    auto M_local_ghost = M_mesh.localGhost();
    auto M_transfer_map = M_mesh.transferMap();

    int Nd = bamgmesh->NodalConnectivitySize[1];
    std::vector<int> dz;
    std::vector<int> ddz_j;
    std::vector<int> ddz_i;

    std::vector<int> d_nnz;
    std::vector<int> o_nnz;

    for (int i=0; i<bamgmesh->NodalConnectivitySize[0]; ++i)
    {

        int counter_dnnz = 0;
        int counter_onnz = 0;

        int Ncc = bamgmesh->NodalConnectivity[Nd*(i+1)-1];
        int gid = M_transfer_map.right.find(i+1)->second;
        //if (std::find(M_local_ghost.begin(),M_local_ghost.end(),gid) == M_local_ghost.end())
        if (!std::binary_search(M_local_ghost.begin(),M_local_ghost.end(),gid))
        {
            for (int j=0; j<Ncc; ++j)
            {
                int currentr = bamgmesh->NodalConnectivity[Nd*i+j];

                int gid2 = M_transfer_map.right.find(currentr)->second;
                //if (std::find(M_local_ghost.begin(),M_local_ghost.end(),gid2) == M_local_ghost.end())
                if (!std::binary_search(M_local_ghost.begin(),M_local_ghost.end(),gid2))
                {
                    ++counter_dnnz;
                }
                else
                {
                    ++counter_onnz;
                }

                //std::cout<<"Connect["<< j <<"]= "<< currentr << " or "<< M_transfer_map.right.find(currentr)->second <<"\n";
            }

            d_nnz.push_back(2*(counter_dnnz+1));
            o_nnz.push_back(2*(counter_onnz));
        }
    }

    auto d_nnz_count = d_nnz.size();
    d_nnz.resize(2*d_nnz_count);
    std::copy_n(d_nnz.begin(), d_nnz_count, d_nnz.begin() + d_nnz_count);

    auto o_nnz_count = o_nnz.size();
    o_nnz.resize(2*o_nnz_count);
    std::copy_n(o_nnz.begin(), o_nnz_count, o_nnz.begin() + o_nnz_count);

    int sM = M_mesh.numNodes();

    std::vector<int> global_indices_with_ghost = M_mesh.localDofWithGhost();
    int glsize = global_indices_with_ghost.size();

    for (int gl=0; gl<glsize; ++gl)
        global_indices_with_ghost[gl] = global_indices_with_ghost[gl]-1;

    std::vector<int> global_indices_without_ghost = M_mesh.localDofWithoutGhost();
    glsize = global_indices_without_ghost.size();

    for (int gl=0; gl<glsize; ++gl)
        global_indices_without_ghost[gl] = global_indices_without_ghost[gl]-1;

    M_graphmpi = graphmpi_type(d_nnz, o_nnz, global_indices_without_ghost, global_indices_with_ghost);

#if 0
    std::cout<<"\n";
    std::cout<<"["<< M_comm.rank() <<"] GRAPHCSR INFO: MIN NZ ON-DIAGONAL (per row)     = "<< *std::min_element(d_nnz.begin(),d_nnz.end()) <<"\n";
    std::cout<<"["<< M_comm.rank() <<"] GRAPHCSR INFO: MAX NZ ON-DIAGONAL (per row)     = "<< *std::max_element(d_nnz.begin(),d_nnz.end()) <<"\n";
    std::cout<<"["<< M_comm.rank() <<"] GRAPHCSR INFO: MIN NZ OFF-DIAGONAL (per row)    = "<< *std::min_element(o_nnz.begin(),o_nnz.end()) <<"\n";
    std::cout<<"["<< M_comm.rank() <<"] GRAPHCSR INFO: MAX NZ OFF-DIAGONAL (per row)    = "<< *std::max_element(o_nnz.begin(),o_nnz.end()) <<"\n";
    std::cout<<"\n";
#endif

#if 0
    M_comm.barrier();

    if (M_rank == 1)
    {
        std::cout<<"************00************\n";
        for (int const& index : global_indices_without_ghost)
            std::cout<<"WITHOUT GHOST "<< index+1 <<"\n";

        std::cout<<"************01************\n";
        for (int const& index : global_indices_with_ghost)
            std::cout<<"WITH GHOST    "<< index+1 <<"\n";

        std::cout<<"************02************\n";
        for (int const& index : M_local_ghost)
            std::cout<<"GHOST         "<< index <<"\n";
    }
#endif

}//createGraph
    

// -------------------------------------------------------------------------------------
//! Exports the model outputs.
//! Called by the checkOutputs() and step() functions.
void
FiniteElement::exportResults(bool const& export_mesh,
        bool const& export_fields, bool const& apply_displacement)
{
    // determine the filename of the output files.
    // if output.datetime_in_filename = false,
    //   an integer "output_step_num" goes into the filename
    // else,
    //   filename is [mesh,field]_%Y%m%dT%H%M%SZ.[bin,dat]

    std::string name_str;
    if (vm["output.datetime_in_filename"].as<bool>())
        name_str = to_date_time_string_for_filename(M_current_time);
    else
    {
        int const output_step_num = pcpt*time_step/output_time_step;
        name_str = (boost::format( "%1%" )
                               % output_step_num ).str();
    }

    this->exportResults(name_str, export_mesh, export_fields, apply_displacement);
}//exportResults


// -------------------------------------------------------------------------------------
//! Exports the model outputs.
//! Called by the other exportResults() interface, init(), step() functions.
void
FiniteElement::exportResults(std::string const& name_str, bool const& export_mesh,
        bool const& export_fields, bool const& apply_displacement)
{
    //define filenames from name_str
    std::string meshfile = (boost::format( "%1%/mesh_%2%" )
                               % M_export_path
                               % name_str ).str();

    std::string fieldfile = (boost::format( "%1%/field_%2%" )
                               % M_export_path
                               % name_str ).str();

    std::vector<std::string> filenames = {meshfile,fieldfile};
    this->exportResults(filenames, export_mesh, export_fields, apply_displacement);
}//exportResults


// -------------------------------------------------------------------------------------
//! Exports the model outputs.
//! Called by the other exportResults() function.
void
FiniteElement::exportResults(std::vector<std::string> const& filenames, bool const& export_mesh,
        bool const& export_fields, bool const& apply_displacement)
{
    std::vector<double> M_VT_root;
    this->gatherNodalField(M_VT,M_VT_root);

    std::vector<double> M_UM_root;
    if (apply_displacement)
    {
        this->gatherNodalField(M_UM,M_UM_root);
    }

    // fields defined on mesh elements

    M_prv_local_ndof = M_local_ndof;
    M_prv_num_nodes = M_num_nodes;
    M_prv_num_elements = M_local_nelements;
    M_prv_global_num_nodes = M_mesh.numGlobalNodes();
    M_prv_global_num_elements = M_mesh.numGlobalElements();

    M_nb_var_element = 15 + M_tice.size();//15;
    int nb_var_element = M_nb_var_element;
    if (M_ice_cat_type!=setup::IceCategoryType::THIN_ICE)
    {
        nb_var_element -= 4;
    }

    if (vm["output.save_diagnostics"].as<bool>())
    {
        nb_var_element += 7;
    }

    std::vector<double> interp_in_elements;
    this->gatherFieldsElementIO(interp_in_elements,M_ice_cat_type==setup::IceCategoryType::THIN_ICE);


    M_comm.barrier();
#if 1
    if (M_rank == 0)
    {
        int num_elements_root = M_mesh_root.numTriangles();
        int tice_size = M_tice.size();

        std::vector<double> M_conc_root(num_elements_root);
        std::vector<double> M_thick_root(num_elements_root);
        std::vector<double> M_snow_thick_root(num_elements_root);
        std::vector<double> M_sigma_root(3*num_elements_root);
        std::vector<double> M_damage_root(num_elements_root);
        std::vector<double> M_ridge_ratio_root(num_elements_root);
        std::vector<double> M_random_number_root(num_elements_root);
        std::vector<double> M_sss_root(num_elements_root);
        std::vector<double> M_sst_root(num_elements_root);
        //std::vector<double> M_tice_root(tice_size*num_elements_root);
        std::vector<double> M_tice0_root(num_elements_root);
        std::vector<double> M_tice1_root;
        std::vector<double> M_tice2_root;

        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            M_tice1_root.resize(num_elements_root);
            M_tice2_root.resize(num_elements_root);
        }

        std::vector<double> M_h_thin_root;
        std::vector<double> M_conc_thin_root;
        std::vector<double> M_hs_thin_root;
        std::vector<double> M_tsurf_thin_root;

        if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_h_thin_root.resize(num_elements_root);
            M_conc_thin_root.resize(num_elements_root);
            M_hs_thin_root.resize(num_elements_root);
            M_tsurf_thin_root.resize(num_elements_root);
        }

        std::vector<double> D_Qa_root;
        std::vector<double> D_Qsw_root;
        std::vector<double> D_Qlw_root;
        std::vector<double> D_Qsh_root;
        std::vector<double> D_Qlh_root;
        std::vector<double> D_Qo_root;
        std::vector<double> D_delS_root;

        if (vm["output.save_diagnostics"].as<bool>())
        {
            D_Qa_root.resize(num_elements_root);
            D_Qsw_root.resize(num_elements_root);
            D_Qlw_root.resize(num_elements_root);
            D_Qsh_root.resize(num_elements_root);
            D_Qlh_root.resize(num_elements_root);
            D_Qo_root.resize(num_elements_root);
            D_delS_root.resize(num_elements_root);
        }

        int tmp_nb_var = 0;

        for (int i=0; i<M_mesh_root.numTriangles(); ++i)
        {
            tmp_nb_var=0;
            int ri = M_rmap_elements[i];

            // concentration
            M_conc_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // thickness
            M_thick_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // snow thickness
            M_snow_thick_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress1
            M_sigma_root[3*i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress2
            M_sigma_root[3*i+1] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress3
            M_sigma_root[3*i+2] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // damage
            M_damage_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // ridge_ratio
            M_ridge_ratio_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // ridge_ratio
            M_random_number_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // sss
            M_sss_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // sst
            M_sst_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // tice1
            M_tice0_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            if ( M_thermo_type == setup::ThermoType::WINTON )
            {
                // tice2
                M_tice1_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;

                // tice3
                M_tice2_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;
            }

            if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                // We also add the "thin" quantities to the "thick" ones
                // h_thin
                M_h_thin_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                M_thick_root[i] += M_h_thin_root[i];
                tmp_nb_var++;

                // conc_thin
                M_conc_thin_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                M_conc_root[i] += M_conc_thin_root[i];
                tmp_nb_var++;

                // hs_thin
                M_hs_thin_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                M_snow_thick_root[i] += M_hs_thin_root[i];
                tmp_nb_var++;

                // tsurf_thin
                M_tsurf_thin_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;
            }

            if (vm["output.save_diagnostics"].as<bool>())
            {
                // Qatm
                D_Qa_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;

                // Qsw
                D_Qsw_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;

                // Qlw
                D_Qlw_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;

                // Qsh
                D_Qsh_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;

                // Qlh
                D_Qlh_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;

                // Qocean
                D_Qo_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;

                // Saltflux
                D_delS_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;
            }

            if(tmp_nb_var>nb_var_element)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }

        Exporter exporter(vm["output.exporter_precision"].as<std::string>());
        std::string fileout;

        if (export_mesh)
        {
            fileout = filenames[0]+".bin";
            LOG(INFO) <<"MESH BINARY: Exporter Filename= "<< fileout <<"\n";

            if(apply_displacement)
            {
                // move the mesh for the export
                M_mesh_root.move(M_UM_root,1.);
            }

            std::fstream meshbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
            if ( !meshbin.good() )
                throw std::runtime_error("Cannot write to file: " + fileout);

            exporter.writeMesh(meshbin, M_mesh_root);
            meshbin.close();

            if(apply_displacement)
            {
                // move it back after the export
                M_mesh_root.move(M_UM_root,-1.);
            }

            fileout = filenames[0]+".dat";

            LOG(INFO) <<"RECORD MESH: Exporter Filename= "<< fileout <<"\n";

            std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
            if ( !outrecord.good() )
                throw std::runtime_error("Cannot write to file: " + fileout);

            exporter.writeRecord(outrecord,"mesh");
            outrecord.close();
        }

        if (export_fields)
        {
            fileout = filenames[1]+".bin";
            LOG(INFO) <<"BINARY: Exporter Filename= "<< fileout <<"\n";

            std::fstream outbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
            if ( !outbin.good() )
                throw std::runtime_error("Cannot write to file: " + fileout);

            std::vector<double> timevec = {M_current_time};
            std::vector<int> regridvec = {M_nb_regrid};

            exporter.writeField(outbin, timevec, "Time");
            exporter.writeField(outbin, regridvec, "M_nb_regrid");
            exporter.writeField(outbin, M_surface_root, "Element_area");
            exporter.writeField(outbin, M_VT_root, "M_VT");
            exporter.writeField(outbin, M_dirichlet_flags_root, "M_dirichlet_flags");
            exporter.writeField(outbin, M_conc_root, "Concentration");
            exporter.writeField(outbin, M_thick_root, "Thickness");
            exporter.writeField(outbin, M_snow_thick_root, "Snow");
            exporter.writeField(outbin, M_damage_root, "Damage");
            exporter.writeField(outbin, M_ridge_ratio_root, "Ridge_ratio");

            // std::vector<double> AllMinAngle = this->AllMinAngle(M_mesh, M_UM, 0.);
            // exporter.writeField(outbin, AllMinAngle, "AllMinAngle");

            // for (int i=0; i<num_elements_root; ++i)
            // {
            //     M_tice[0][i] = M_tice_root[tice_size*i];

            //     if ( M_thermo_type == setup::ThermoType::WINTON )
            //     {
            //         M_tice[1][i] = M_tice_root[tice_size*i+1];
            //         M_tice[2][i] = M_tice_root[tice_size*i+2];
            //     }
            // }

            // int i=0;
            // for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
            // {
            //     exporter.writeField(outbin, *it, "Tice_"+std::to_string(i));
            //     i++;
            // }

            exporter.writeField(outbin, M_tice0_root, "Tice_0");
            if ( M_thermo_type == setup::ThermoType::WINTON )
            {
                exporter.writeField(outbin, M_tice1_root, "Tice_1");
                exporter.writeField(outbin, M_tice2_root, "Tice_2");
            }

            exporter.writeField(outbin, M_sst_root, "SST");
            exporter.writeField(outbin, M_sss_root, "SSS");

            if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                exporter.writeField(outbin, M_h_thin_root, "Thin_ice");
                exporter.writeField(outbin, M_hs_thin_root, "Snow_thin_ice");
                exporter.writeField(outbin, M_tsurf_thin_root, "Tsurf_thin_ice");
                exporter.writeField(outbin, M_conc_thin_root, "Concentration_thin_ice");
            }

#if 1
            // EXPORT sigma1 sigma2
            std::vector<double> sigma1(M_mesh_root.numTriangles());
            std::vector<double> sigma2(M_mesh_root.numTriangles());
            double sigma_s, sigma_n;
            std::vector<double> sigma_pred(3);

            for ( int i=0; i<M_mesh_root.numTriangles(); ++i )
            {
                sigma_pred[0]=M_sigma_root[3*i];
                sigma_pred[1]=M_sigma_root[3*i+1];
                sigma_pred[2]=M_sigma_root[3*i+2];

                sigma_s=std::hypot((sigma_pred[0]-sigma_pred[1])/2.,sigma_pred[2]);
                sigma_n= -         (sigma_pred[0]+sigma_pred[1])/2.;

                sigma1[i] = sigma_n+sigma_s;
                sigma2[i] = sigma_n-sigma_s;
            }
            exporter.writeField(outbin, sigma1, "Sigma1");
            exporter.writeField(outbin, sigma2, "Sigma2");
#endif

            if (vm["output.save_diagnostics"].as<bool>())
            {
                exporter.writeField(outbin, D_Qa_root, "Qatm");
                exporter.writeField(outbin, D_Qsw_root, "Qsw");
                exporter.writeField(outbin, D_Qlw_root, "Qlw");
                exporter.writeField(outbin, D_Qsh_root, "Qsh");
                exporter.writeField(outbin, D_Qlh_root, "Qlh");
                exporter.writeField(outbin, D_Qo_root,  "Qocean");
                exporter.writeField(outbin, D_delS_root, "Saltflux");
            }
            outbin.close();

            fileout = filenames[1]+".dat";
            LOG(INFO) <<"RECORD FIELD: Exporter Filename= "<< fileout <<"\n";

            std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
            if ( !outrecord.good() )
                throw std::runtime_error("Cannot write to file: " + fileout);

            exporter.writeRecord(outrecord);
            outrecord.close();
        }
    }

#endif
}// exportResults()

    
// -------------------------------------------------------------------------------------
//! Gets GitHub revision version of the model code.
//! Called by the writeLogFile() function.
std::string
FiniteElement::gitRevision()
{
    //std::string command = "git rev-parse HEAD";
    return this->system("git rev-parse HEAD");
}//gitRevision

    
// -------------------------------------------------------------------------------------
//! Run a system command
//! Called by the writeLogFile(), gitRevision(), createGmshMesh() functions.
std::string
FiniteElement::system(std::string const& command)
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


// -------------------------------------------------------------------------------------
//! Gets the location of the libraries.
//! Called by the writeLogFile() function.
std::string
FiniteElement::getEnv(std::string const& envname)
{
    const char* senv = ::getenv(envname.c_str());
    if ( senv == NULL )
        senv = "NULL";
    return std::string(senv);
}//getEnv
    

// -------------------------------------------------------------------------------------
//! Writes a log file with the location of libraries and system information.
//! Called by the run() function.
void
FiniteElement::writeLogFile()
{
    std::string logfilename = "";
    if ((vm["output.logfile"].as<std::string>()).empty())
    {
        logfilename = "nextsim.log";
    }
    else
    {
        logfilename = vm["output.logfile"].as<std::string>();
    }

    std::string fileout = (boost::format( "%1%/%2%" )
               % M_export_path
               % logfilename ).str();

    std::fstream logfile(fileout, std::ios::out | std::ios::trunc);
    std::cout << "Writing log file " << fileout << "...\n";

    int log_width = 55;
    if (logfile.is_open())
    {
        logfile << "#----------Info\n";
        logfile << std::setw(log_width) << std::left << "Build date "  << Nextsim::current_time_local() <<"\n";
        logfile << std::setw(log_width) << std::left << "Git revision "  << gitRevision() <<"\n";

        logfile << "#----------Compilers\n";
        logfile << std::setw(log_width) << std::left << "C "  << system("which gcc") << " (version "<< system("gcc -dumpversion") << ")" <<"\n";
        logfile << std::setw(log_width) << std::left << "C++ "  << system("which g++") << " (version "<< system("g++ -dumpversion") << ")" <<"\n";

        logfile << "#----------Environment variables\n";
        logfile << std::setw(log_width) << std::left << "NEXTSIM_DATA_DIR "  << getEnv("NEXTSIM_DATA_DIR") <<"\n";
        logfile << std::setw(log_width) << std::left << "NEXTSIM_MESH_DIR "  << getEnv("NEXTSIM_MESH_DIR") <<"\n";

        logfile << "#----------Program options\n";

        for (po::variables_map::iterator it = vm.begin(); it != vm.end(); it++)
        {
            // ignore wim options if no coupling
#if !defined (WAVES)
            if ((it->first.find("nextwim.") != std::string::npos) || (it->first.find("wim.") != std::string::npos))
            {
                continue;
            }
#endif

            logfile << std::setw(log_width) << std::left << it->first;

#if 0
            if (((boost::any)it->second.value()).empty())
            {
                std::cout << "(empty)";
            }
            if (vm[it->first].defaulted() || it->second.defaulted()) {
                std::cout << "(default)";
            }
#endif

            bool is_char;
            try
            {
                boost::any_cast<const char *>(it->second.value());
                is_char = true;
            }
            catch (const boost::bad_any_cast &)
            {
                is_char = false;
            }

            bool is_str;
            try
            {
                boost::any_cast<std::string>(it->second.value());
                is_str = true;
            }
            catch (const boost::bad_any_cast &)
            {
                is_str = false;
            }

            if (((boost::any)it->second.value()).type() == typeid(int))
            {
                logfile << vm[it->first].as<int>() <<"\n";
            }
            else if (((boost::any)it->second.value()).type() == typeid(bool))
            {
                logfile << vm[it->first].as<bool>() <<"\n";
            }
            else if (((boost::any)it->second.value()).type() == typeid(double))
            {
                logfile << vm[it->first].as<double>() <<"\n";
            }
            else if (is_char)
            {
                logfile << vm[it->first].as<const char * >() <<"\n";
            }
            else if (is_str)
            {
                std::string temp = vm[it->first].as<std::string>();

                logfile << temp <<"\n";
#if 0
                if (temp.size())
                {
                    logfile << temp <<"\n";
                }
                else
                {
                    logfile << "true" <<"\n";
                }
#endif
            }
            else
            { // Assumes that the only remainder is vector<string>
                try
                {
                    std::vector<std::string> vect = vm[it->first].as<std::vector<std::string> >();
                    if(vect.size()==0)
                        //make sure we go to the next line if vector is empty
                        logfile<<"\n";

                    uint i = 0;
                    for (auto oit=vect.begin(); oit != vect.end(); oit++, ++i)
                    {
                        //logfile << it->first << "[" << i << "]=" << (*oit) <<"\n";
                        if (i > 0)
                            logfile << std::setw(log_width) << std::right<<" ";

                        logfile << "[" << i << "]=" << (*oit) <<"\n";
                    }
                }
                catch (const boost::bad_any_cast &)
                {
                    std::cout << "UnknownType(" << ((boost::any)it->second.value()).type().name() << ")" <<"\n";
                }
            }
        }
    }

    //move git_changes.txt from current dir to output dir
    fs::path path1("git_changes.txt");
    if ( fs::exists(path1) )
    {
        // try to rename, otherwise copy and remove
        fs::path path2(M_export_path+"/git_changes.txt");
        try
        {
            fs::rename(path1,path2);
        }
        catch (const boost::filesystem::filesystem_error &)
        {
            std::ifstream in (path1.string());
            std::ofstream out(path2.string());
            out << in.rdbuf();
            out.close();
            in.close();
            fs::remove(path1);
        }
    }
}//writeLogFile

    
// -------------------------------------------------------------------------------------
//! Checks fields for NaNs and for too big ice thickness values.
//! Called by the step() function.
void
FiniteElement::checkFields()
{

    int itest = vm["debugging.test_element_number"].as<int>();
    bool printout = (
            M_rank == vm["debugging.test_proc_number"].as<int>()
            && itest>0);
    double xtest = 0.;
    double ytest = 0.;
    double lat_test = 0.;
    double lon_test = 0.;
    if(printout)
    {
        auto movedmesh = M_mesh;
        movedmesh.move(M_UM, 1.);
        xtest = movedmesh.bCoordX()[itest];
        ytest = movedmesh.bCoordY()[itest];

        // get lon, lat at test position
        mapx_class *map;
        std::string configfile = Environment::nextsimMppfile();
        std::vector<char> str(configfile.begin(), configfile.end());
        str.push_back('\0');
        map = init_mapx(&str[0]);
        inverse_mapx(map, xtest, ytest, &lat_test, &lon_test);
        close_mapx(map);
    }
    std::vector< std::vector<double>* > vecs_to_check;
    std::vector< std::string > vec_names;
    std::vector< ExternalData* > forcings_to_check;
    std::vector< std::string > forcing_names;

    vecs_to_check.push_back(&M_conc);
    vec_names.push_back("M_conc");
    vecs_to_check.push_back(&M_thick);
    vec_names.push_back("M_thick");
    vecs_to_check.push_back(&(M_tice[0]));
    vec_names.push_back("M_tice[0]");
    vecs_to_check.push_back(&(M_sss));
    vec_names.push_back("M_sss");
    vecs_to_check.push_back(&(M_sst));
    vec_names.push_back("M_sst");

    forcings_to_check.push_back(&M_ocean_temp);
    forcing_names.push_back("SST_forcing");
    forcings_to_check.push_back(&M_ocean_salt);
    forcing_names.push_back("SSS_forcing");

    for(int i=0; i<M_num_elements; i++)
    {
        int j = 0;
        if(printout && i==itest)
        {
            std::cout<<"In checkFields\n";
            std::cout<<"pcpt =  "<<pcpt<<"\n";
            std::cout<<"date =  "<<to_date_time_string(M_current_time)<<"\n";
            std::cout<<"M_nb_regrid = "<<M_nb_regrid<<"\n";
            std::cout<<"M_rank = "<<M_rank<<"\n";
            std::cout<<"element number = "<<i<<"\n";
            std::cout<<"x,y = " <<xtest <<"," <<ytest <<"\n";
            std::cout<<"lon,lat = " <<lon_test <<"," <<lat_test <<"\n";
            for (int j=0; j<vec_names.size(); j++)
            {
                double val = ( *(vecs_to_check[j]) )[i];//vecs_to_check[j] is a pointer, so dereference
                std::string name = vec_names[j];
                std::cout<<name <<" = "<< val <<"\n";
            }
            for (int j=0; j<forcing_names.size(); j++)
            {
                double val = forcings_to_check[j]->get(i);//forcings_to_check[j] is a pointer
                std::string name = forcing_names[j];
                std::cout<<name <<" = "<< val <<"\n";
            }
            std::cout<<"\n";
        }
        for (int j=0; j<vec_names.size(); j++)
        {
            double val = ( *(vecs_to_check[j]) )[i];//vecs_to_check[j] is a pointer, so dereference
            std::string name = vec_names[j];
            if(std::isnan(val))
            {
                std::cout<<"NaN in "<<name<<"\n";
                std::cout<<"pcpt =  "<<pcpt<<"\n";
                std::cout<<"M_nb_regrid = "<<M_nb_regrid<<"\n";
                std::cout<<"date =  "<<to_date_time_string(M_current_time)<<"\n";
                std::cout<<"M_rank = "<<M_rank<<"\n";
                std::cout<<"element number = "<<i<<"\n";
                throw std::runtime_error("found NaN");
            }
            if(name=="M_thick" && val>25.)
            {
                std::cout<<"M_thick too big(" <<val <<")\n";
                std::cout<<"pcpt =  "<<pcpt<<"\n";
                std::cout<<"M_nb_regrid = "<<M_nb_regrid<<"\n";
                std::cout<<"date =  "<<to_date_time_string(M_current_time)<<"\n";
                std::cout<<"M_rank = "<<M_rank<<"\n";
                std::cout<<"element number = "<<i<<"\n";
                throw std::runtime_error("M_thick too big");
            }
        }
    }
}//checkFields

    
// -------------------------------------------------------------------------------------
//! Finalizes the run: clears meshes and some matrices used by the solver.
//! Called by the step() function.
void
FiniteElement::finalise()
{
    // Don't forget to close the iabp file!
    if (M_use_iabp_drifters)
    {
        M_iabp_infile_fstream.close();
    }

    // Clear ponters etc
    M_comm.barrier();

    delete bamgmesh;
    delete bamggeom;

    if (M_mesh.comm().rank() == 0)
    {
        delete bamgopt;
        delete bamggeom_root;
        delete bamgmesh_root;

        // We need to point these to NULL because 'delete bamgopt' clears the
        // memory they were pointing to before
        bamgopt_previous->hminVertices      = NULL;
        bamgopt_previous->hmaxVertices      = NULL;

        delete bamgmesh_previous;
        delete bamggeom_previous;
        delete bamgopt_previous;

        // clear GModel from mesh data structure
        if (M_partition_space == mesh::PartitionSpace::MEMORY)
        {
            M_mesh_root.clear();
        }
    }

    M_matrix->clear();
    M_vector->clear();
    M_solution->clear();
    M_solver->clear();
}//finalise

} // Nextsim
