/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   finiteelement.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
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
FiniteElement::FiniteElement()
    :
    vm(Environment::vm()),
    M_mesh(),
    M_solver(),
    M_matrix(),
    M_vector(),
    timer()
{}

// Initialisation of the mesh and forcing
void
FiniteElement::initMesh()
{
    switch (M_domain_type)
    {
        case setup::DomainType::DEFAULT:
            M_flag_fix = 10000; // free = [10001 10002];
            break;
        case setup::DomainType::KARA:
            M_flag_fix = 17; // free = [15 16];
            break;
        case setup::DomainType::BERINGKARA:
            M_flag_fix = 1; // free = [];
            break;
        case setup::DomainType::BIGKARA:
            M_flag_fix = 158; // free = 157;
            break;
        case setup::DomainType::ARCTIC:
            M_flag_fix = 174; // free = [172 173];
            break;
        case setup::DomainType::BIGARCTIC:
            M_flag_fix = 161; // free = 158:160;
            break;
        default:
            std::cout << "invalid domain type"<<"\n";
            throw std::logic_error("invalid domain type");
    }

    this->initBamg();

    M_comm = M_mesh.comm();
    M_rank = M_comm.rank();

    this->rootMeshProcessing();

    if (!M_use_restart)
    {
        this->distributedMeshProcessing(true);
    }
}

void
FiniteElement::distributedMeshProcessing(bool start)
{
    // if (M_mesh_filename.substr(3,1) != std::to_string(Environment::comm().size()))
    // {
    //     std::cout<<"the number of processor cores does not match the number of mesh partitions: "
    //              << std::to_string(Environment::comm().size()) <<" != " << M_mesh_filename.substr(3,1) <<"\n";

    //     throw std::logic_error("invalid number of processor cores or number of mesh partitions");
    // }

    M_comm.barrier();

    if (!start)
    {
        M_mesh = mesh_type();
    }

    M_mesh.setOrdering("gmsh");

    if (M_rank == 0)
    {
        //LOG(INFO) <<"["<< M_rank <<"] " <<"filename= "<< M_mesh_filename <<"\n";
        std::cout << "[INFO]: " <<"["<< M_rank <<"] " <<"filename= "<< M_mesh_filename <<"\n";
    }

    timer["meshread"].first.restart();
    M_mesh.readFromFile(M_mesh_filename, M_mesh_fileformat);
    //M_mesh.readFromFileBinary(M_mesh_filename);
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

}

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
}

void
FiniteElement::rootMeshProcessing()
{
    if (M_rank == 0)
    {
        M_mesh_root.setOrdering("bamg");
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
        case setup::MeshType::FROM_GMSH:
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
            M_mesh_filename = (boost::format( "par%1%%2%" ) % M_comm.size() % M_mesh_filename ).str();
            LOG(DEBUG) <<"["<< M_rank <<"] " <<"filename= "<< M_mesh_filename <<"\n";

            std::cout<<"------------------------------version       = "<< M_mesh_root.version() <<"\n";
            std::cout<<"------------------------------ordering      = "<< M_mesh_root.ordering() <<"\n";
            std::cout<<"------------------------------format        = "<< M_mesh_fileformat <<"\n";
            std::cout<<"------------------------------space         = "<< vm["mesh.partition-space"].as<std::string>() <<"\n";
            std::cout<<"------------------------------partitioner   = "<< vm["mesh.partitioner"].as<std::string>() <<"\n";


            // save mesh (only root process)
            chrono.restart();
            if (M_partition_space == mesh::PartitionSpace::MEMORY)
            {
                // Environment::logMemoryUsage("before gmodel...");
                M_mesh_root.initGModel();
                M_mesh_root.writeToGModel(M_mesh_filename);
                // Environment::logMemoryUsage("before after...");
            }
            else if (M_partition_space == mesh::PartitionSpace::DISK)
            {
                M_mesh_root.writeTofile(M_mesh_filename);
            }
            //LOG(DEBUG) <<"Saving mesh done in "<< chrono.elapsed() <<"s\n";
            std::cout <<"Writing mesh done in "<< chrono.elapsed() <<"s\n";

            // partition the mesh on root process (rank 0)
            chrono.restart();
            M_mesh_root.partition(M_mesh_filename,M_partitioner,M_partition_space, M_mesh_fileformat);
            //LOG(DEBUG) <<"Partitioning mesh done in "<< chrono.elapsed() <<"s\n";
            std::cout <<"Partitioning mesh done in "<< chrono.elapsed() <<"s\n";
        }
    }
    else
    {
        // Add information on the number of partition to mesh filename
        M_mesh_filename = (boost::format( "par%1%%2%" ) % M_comm.size() % M_mesh_filename ).str();
        // LOG(DEBUG) <<"["<< M_rank <<"] " <<"filename= "<< M_mesh_filename <<"\n";
    }
}

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
}

// Initialise size of all physical variables with values set to zero
void
FiniteElement::initVariables()
{
    chrono_tot.restart();

    M_solver = solver_ptrtype(new solver_type());
    M_matrix = matrix_ptrtype(new matrix_type());
    M_vector = vector_ptrtype(new vector_type());
    M_solution = vector_ptrtype(new vector_type());

    M_reuse_prec = true;
    M_nb_regrid = 0;

    M_VT.resize(2*M_num_nodes,0.);
    M_VTM.resize(2*M_num_nodes,0.);
    M_VTMM.resize(2*M_num_nodes,0.);

    M_sst.resize(M_num_elements);
    M_sss.resize(M_num_elements);

    M_h_thin.assign(M_num_elements,0.);
    M_conc_thin.assign(M_num_elements,0.);
    M_hs_thin.assign(M_num_elements,0.);
    M_tsurf_thin.assign(M_num_elements,0.);

    // stresses
    M_sigma.resize(3*M_num_elements,0.);

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

    // if (M_rank == 0)
    // {
    //     std::cout<<"Random MIN= "<< *std::min_element(M_random_number_root.begin(), M_random_number_root.end()) <<"\n";
    //     std::cout<<"Random MAX= "<< *std::max_element(M_random_number_root.begin(), M_random_number_root.end()) <<"\n";
    // }

    boost::mpi::broadcast(M_comm, &M_random_number_root[0], M_mesh.numGlobalElements(), 0);

    M_random_number.resize(M_num_elements);
    auto id_elements = M_mesh.trianglesIdWithGhost();

    for (int i=0; i<M_random_number.size(); ++i)
    {
        M_random_number[i] = M_random_number_root[id_elements[i]-1];
    }

    M_conc.resize(M_num_elements);
    M_thick.resize(M_num_elements);
    M_damage.resize(M_num_elements);
    M_ridge_ratio.resize(M_num_elements);
    M_snow_thick.resize(M_num_elements);

    switch (M_thermo_type)
    {
    case (setup::ThermoType::ZERO_LAYER):
        M_tice.resize(1);
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
        if ((M_conc[i] <= 0.) || (M_thick[i] <= 0.) )
        {
            M_conc[i] = 0.;
            M_thick[i] = 0.;
        }
    }

    // Diagnostics
    D_Qa.resize(M_num_elements);
    D_Qsh.resize(M_num_elements);
    D_Qlh.resize(M_num_elements);
    D_Qlw.resize(M_num_elements);
    D_Qsw.resize(M_num_elements);
    D_Qo.resize(M_num_elements);
    D_delS.resize(M_num_elements);

    this->assignVariables();
}

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
    M_UT.assign(2*M_num_nodes,0.);
    M_fcor.resize(M_num_elements);

    M_atmosphere_nodes_dataset.target_size=M_num_nodes;
    M_atmosphere_elements_dataset.target_size=M_num_elements;
    M_atmosphere_bis_elements_dataset.target_size=M_num_elements;
    M_ocean_nodes_dataset.target_size=M_num_nodes;
    M_ocean_elements_dataset.target_size=M_num_elements;

    M_ice_topaz_elements_dataset.target_size=M_num_elements;
    M_ice_piomas_elements_dataset.target_size=M_num_elements;
    M_ice_amsre_elements_dataset.target_size=M_num_elements;
    M_ice_osisaf_elements_dataset.target_size=M_num_elements;
    M_ice_osisaf_type_elements_dataset.target_size=M_num_elements;
    M_ice_amsr2_elements_dataset.target_size=M_num_elements;
    M_ice_cs2_smos_elements_dataset.target_size=M_num_elements;
    M_ice_smos_elements_dataset.target_size=M_num_elements;
    M_bathymetry_elements_dataset.target_size=M_num_elements;

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

    M_Cohesion.resize(M_num_elements);
    M_Compressive_strength.resize(M_num_elements);
    M_time_relaxation_damage.resize(M_num_elements,time_relaxation_damage);

    // root
    M_UM_root.assign(2*M_mesh.numGlobalNodes(),0.);

    // number of variables to interpolate
    M_nb_var_element = /*11*/15 + M_tice.size();
}

void
FiniteElement::initModelState()
{
    // Initialise the physical state of the model

    this->initIce();

    this->initSlabOcean();

    this->initDrifter();
}

void
FiniteElement::initDatasets()
{
    // Definition of the datasets
    switch(M_atmosphere_type){
        case setup::AtmosphereType::CONSTANT:
            break;

        case setup::AtmosphereType::ASR:
            M_atmosphere_nodes_dataset=DataSet("asr_nodes",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("asr_elements",M_num_elements);
            break;

        case setup::AtmosphereType::ERAi:
            M_atmosphere_nodes_dataset=DataSet("ERAi_nodes",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("ERAi_elements",M_num_elements);
            break;

        case setup::AtmosphereType::EC:
            M_atmosphere_nodes_dataset=DataSet("ec_nodes",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("ec_elements",M_num_elements);
            break;

        case setup::AtmosphereType::EC2:
            M_atmosphere_nodes_dataset=DataSet("ec2_nodes",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("ec2_elements",M_num_elements);
            break;

        case setup::AtmosphereType::EC_ERAi:
            M_atmosphere_nodes_dataset=DataSet("ec_nodes",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("ec_elements",M_num_elements);
            M_atmosphere_bis_elements_dataset=DataSet("ERAi_elements",M_num_elements);
            break;

        case setup::AtmosphereType::CFSR:
            M_atmosphere_nodes_dataset=DataSet("cfsr_nodes",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("cfsr_elements",M_num_elements);
            break;

        case setup::AtmosphereType::CFSR_HI:
            M_atmosphere_nodes_dataset=DataSet("cfsr_nodes_hi",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("cfsr_elements",M_num_elements);
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
            M_ocean_nodes_dataset=DataSet("topaz_nodes",M_num_nodes);
            M_ocean_elements_dataset=DataSet("topaz_elements",M_num_elements);
            break;

        case setup::OceanType::TOPAZR_ALTIMETER:
            M_ocean_nodes_dataset=DataSet("ocean_currents_nodes",M_num_nodes);
            M_ocean_elements_dataset=DataSet("topaz_elements",M_num_elements);
            break;

        case setup::OceanType::TOPAZF:
            M_ocean_nodes_dataset=DataSet("topaz_forecast_nodes",M_num_nodes);
            M_ocean_elements_dataset=DataSet("topaz_forecast_elements",M_num_elements);
            break;

        default:
            std::cout << "invalid ocean forcing"<<"\n";throw std::logic_error("invalid wind forcing");
    }

    M_ice_topaz_elements_dataset=DataSet("ice_topaz_elements",M_num_elements);

    M_ice_icesat_elements_dataset=DataSet("ice_icesat_elements",M_num_elements);

    M_ice_piomas_elements_dataset=DataSet("ice_piomas_elements",M_num_elements);

    M_ice_amsre_elements_dataset=DataSet("ice_amsre_elements",M_num_elements);

    M_ice_osisaf_elements_dataset=DataSet("ice_osisaf_elements",M_num_elements);

    M_ice_amsr2_elements_dataset=DataSet("ice_amsr2_elements",M_num_elements);

    M_ice_cs2_smos_elements_dataset=DataSet("ice_cs2_smos_elements",M_num_elements);

    M_ice_smos_elements_dataset=DataSet("ice_smos_elements",M_num_elements);

    M_bathymetry_elements_dataset=DataSet("etopo_elements",M_num_elements);//M_num_nodes);

    // datasets that need to be re-interpolated after regridding
    // - not needed if only used at initialisation, or if not interpolated onto
    // mesh (eg wave datasets are interpolated onto a fixed grid)
    M_datasets_regrid.push_back(&M_atmosphere_nodes_dataset);
    M_datasets_regrid.push_back(&M_atmosphere_elements_dataset);
    M_datasets_regrid.push_back(&M_atmosphere_bis_elements_dataset);
    M_datasets_regrid.push_back(&M_ocean_nodes_dataset);
    M_datasets_regrid.push_back(&M_ocean_elements_dataset);
    M_datasets_regrid.push_back(&M_bathymetry_elements_dataset);
}

void
FiniteElement::checkReloadDatasets(external_data_vec const& ext_data_vec,
        double const& CRtime, std::string const& printout)
{
    //std::cout<<"size of external data vector = "<<ext_data_vec.size()<<"\n";
    if ( ext_data_vec.size()==0 )
    {
        LOG(DEBUG) <<"check_and_reload ("<<printout<<"):\n";
        LOG(DEBUG) <<"nothing to do\n";
        return;
    }

    //loop over M_external_data and call check and reload for each:
    chrono.restart();
    LOG(DEBUG) <<"check_and_reload ("<<printout<<") starts\n";

    //don't rotate yet since rotation angle not always defined yet
    auto RX_nod = M_mesh.coordX ();//nodes
    auto RY_nod = M_mesh.coordY ();
    auto RX_el  = M_mesh.bCoordX();//elements
    auto RY_el  = M_mesh.bCoordY();

    for ( auto it = ext_data_vec.begin(); it != ext_data_vec.end(); ++it )
    {
        if ( (*it)->M_is_constant )
        {
            std::vector<double> RX,RY;//not needed for constant forcings
            (*it)->check_and_reload(RX,RY,CRtime);
        }
        else
        {
            LOG(DEBUG)<<" dataset = "<<(*it)->M_dataset->name<<"\n";
            //dataset & interpolation etc needed
            if ( (*it)->M_dataset->grid.target_location=="mesh_nodes" )
            {
                LOG(DEBUG)<<"in nodes: dataset = "<<(*it)->M_dataset->name<<"\n";
                (*it)->check_and_reload(RX_nod,RY_nod,CRtime);
            }
            else if ( (*it)->M_dataset->grid.target_location=="mesh_elements" )
            {
                LOG(DEBUG)<<"in elements: dataset = "<<(*it)->M_dataset->name<<"\n";
                (*it)->check_and_reload(RX_el,RY_el,CRtime);
            }
#if defined (WAVES)
            else if ( (*it)->M_dataset->grid.target_location=="wim_grid" )
            {
                //LOG(DEBUG)<<"in wim_grid: dataset = "<<(*it)->M_dataset->name<<"\n";
                std::cout<<"in wim_grid: dataset = "<<(*it)->M_dataset->name<<"\n";
                //interp to WIM grid
                (*it)->check_and_reload(wim_grid.X,wim_grid.Y,CRtime);
            }
#endif
            else
            {
                std::cout<<"Bad value for dataset->grid.target_location: "<<(*it)->M_dataset->grid.target_location<<"\n";
                std::cout<<"- set to \"mesh_nodes\", or \"mesh_elements\"\n";
#if defined (WAVES)
                std::cout<<"or \"wim_grid\"\n";
#endif
                std::abort();
            }
        }
    }

    RX_el.resize(0);
    RY_el.resize(0);
    RX_nod.resize(0);
    RY_nod.resize(0);
    LOG(DEBUG) <<"check_and_reload in "<< chrono.elapsed() <<"s\n";
}//checkReloadDatasets


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
    bamgopt->verbose           = vm["simul.verbose"].as<int>();

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
}

void
FiniteElement::initConstant()
{
    nu0 = vm["simul.nu0"].as<double>();
    young = vm["simul.young"].as<double>();
    rhoi = physical::rhoi;
    rhos = physical::rhos;

    days_in_sec = 24.0*3600.0;
    time_init = from_date_time_string(vm["simul.time_init"].as<std::string>());
    output_time_step =  days_in_sec/vm["simul.output_per_day"].as<int>();
    ptime_step =  days_in_sec/vm["simul.ptime_per_day"].as<int>();

    time_step = vm["simul.timestep"].as<double>();
    duration = (vm["simul.duration"].as<double>())*days_in_sec;


    mooring_output_time_step =  vm["simul.mooring_output_timestep"].as<double>()*days_in_sec;
    mooring_time_factor = time_step/mooring_output_time_step;
    //drifter_output_time_step =  vm["simul.drifter_output_timestep"].as<double>()*days_in_sec;

    restart_time_step =  vm["setup.restart_time_step"].as<double>()*days_in_sec;
    M_use_restart   = vm["setup.use_restart"].as<bool>();
    M_write_restart = vm["setup.write_restart"].as<bool>();
    if ( fmod(restart_time_step,time_step) != 0)
    {
        std::cout << restart_time_step << " " << time_step << "\n";
        throw std::logic_error("restart_time_step not an integer multiple of time_step");
    }

    ocean_turning_angle_rad = (PI/180.)*vm["simul.oceanic_turning_angle"].as<double>();
    ridging_exponent = vm["simul.ridging_exponent"].as<double>();

    quad_drag_coef_water = vm["simul.quad_drag_coef_water"].as<double>();

    basal_k2 = vm["simul.Lemieux_basal_k2"].as<double>();
    basal_u_0 = vm["simul.Lemieux_basal_u_0"].as<double>();
    basal_Cb = vm["simul.Lemieux_basal_Cb"].as<double>();

    time_relaxation_damage = vm["simul.time_relaxation_damage"].as<double>()*days_in_sec;
    deltaT_relaxation_damage = vm["simul.deltaT_relaxation_damage"].as<double>()/**days_in_sec*/;

    h_thin_min = vm["simul.h_thin_min"].as<double>();
    h_thin_max = vm["simul.h_thin_max"].as<double>();

    compr_strength = vm["simul.compr_strength"].as<double>();
    tract_coef = vm["simul.tract_coef"].as<double>();
    scale_coef = vm["simul.scale_coef"].as<double>();
    alea_factor = vm["simul.alea_factor"].as<double>();
    cfix = vm["simul.cfix"].as<double>();

    C_fix    = cfix*scale_coef;          // C_fix;...  : cohesion (mohr-coulomb) in MPa (40000 Pa)
    C_alea   = alea_factor*C_fix;        // C_alea;... : alea sur la cohesion (Pa)
    tan_phi = vm["simul.tan_phi"].as<double>();
    //ridge_h = vm["simul.ridge_h"].as<double>();

    if ( vm["simul.newice_type"].as<int>() == 4 )
        M_ice_cat_type = setup::IceCategoryType::THIN_ICE;
    else
        M_ice_cat_type = setup::IceCategoryType::CLASSIC;

    const boost::unordered_map<const std::string, setup::ThermoType> str2thermo = boost::assign::map_list_of
        ("zero-layer", setup::ThermoType::ZERO_LAYER)
        ("winton", setup::ThermoType::WINTON);
    M_thermo_type = str2thermo.find(vm["setup.thermo-type"].as<std::string>())->second;

    const boost::unordered_map<const std::string, setup::AtmosphereType> str2atmosphere = boost::assign::map_list_of
        ("constant", setup::AtmosphereType::CONSTANT)
        ("asr", setup::AtmosphereType::ASR)
        ("erai", setup::AtmosphereType::ERAi)
        ("ec", setup::AtmosphereType::EC)
        ("ec2", setup::AtmosphereType::EC2)
        ("ec_erai", setup::AtmosphereType::EC_ERAi)
        ("cfsr", setup::AtmosphereType::CFSR)
        ("cfsr_hi", setup::AtmosphereType::CFSR_HI);
    M_atmosphere_type = str2atmosphere.find(vm["setup.atmosphere-type"].as<std::string>())->second;

    switch(M_atmosphere_type)
    {
        case setup::AtmosphereType::CONSTANT:
            quad_drag_coef_air = vm["simul.ASR_quad_drag_coef_air"].as<double>();
            break;
        case setup::AtmosphereType::ASR:
            quad_drag_coef_air = vm["simul.ASR_quad_drag_coef_air"].as<double>();
            break;
        case setup::AtmosphereType::CFSR_HI:
        case setup::AtmosphereType::CFSR:
            quad_drag_coef_air = vm["simul.CFSR_quad_drag_coef_air"].as<double>();
            break;
        case setup::AtmosphereType::ERAi:
            quad_drag_coef_air = vm["simul.ERAi_quad_drag_coef_air"].as<double>();
            break;
        case setup::AtmosphereType::EC:
        case setup::AtmosphereType::EC2:
        case setup::AtmosphereType::EC_ERAi:
            quad_drag_coef_air = vm["simul.ECMWF_quad_drag_coef_air"].as<double>();
            break;
        default:
            std::cout << "invalid wind forcing"<<"\n";throw std::logic_error("invalid wind forcing");
    }

    const boost::unordered_map<const std::string, setup::OceanType> str2ocean = boost::assign::map_list_of
        ("constant", setup::OceanType::CONSTANT)
        ("topaz", setup::OceanType::TOPAZR)
        ("topaz_atrest", setup::OceanType::TOPAZR_atrest)
        ("topaz_forecast", setup::OceanType::TOPAZF)
        ("topaz_altimeter", setup::OceanType::TOPAZR_ALTIMETER);
    M_ocean_type = str2ocean.find(vm["setup.ocean-type"].as<std::string>())->second;

    const boost::unordered_map<const std::string, setup::IceType> str2conc = boost::assign::map_list_of
        ("constant", setup::IceType::CONSTANT)
        ("constant_partial", setup::IceType::CONSTANT_PARTIAL)
        ("target", setup::IceType::TARGET)
        ("binary", setup::IceType::BINARY)
        ("topaz", setup::IceType::TOPAZ4)
        ("topaz_forecast", setup::IceType::TOPAZ4F)
        ("topaz_forecast_amsr2", setup::IceType::TOPAZ4FAMSR2)
        ("topaz_forecast_amsr2_osisaf", setup::IceType::TOPAZ4FAMSR2OSISAF)
        ("amsre", setup::IceType::AMSRE)
        ("amsr2", setup::IceType::AMSR2)
        ("osisaf", setup::IceType::OSISAF)
        ("piomas", setup::IceType::PIOMAS)
        ("cs2_smos", setup::IceType::CS2_SMOS)
        ("cs2_smos_amsr2", setup::IceType::CS2_SMOS_AMSR2)
        ("smos", setup::IceType::SMOS)
        ("topaz_osisaf_icesat", setup::IceType::TOPAZ4OSISAFICESAT);
    M_ice_type = str2conc.find(vm["setup.ice-type"].as<std::string>())->second;

    const boost::unordered_map<const std::string, setup::DynamicsType> str2dynamics = boost::assign::map_list_of
        ("default", setup::DynamicsType::DEFAULT)
        ("no_motion", setup::DynamicsType::NO_MOTION)
        ("free_drift", setup::DynamicsType::FREE_DRIFT);
    M_dynamics_type = str2dynamics.find(vm["setup.dynamics-type"].as<std::string>())->second;

    const boost::unordered_map<const std::string, setup::BathymetryType> str2bathymetry = boost::assign::map_list_of
        ("constant", setup::BathymetryType::CONSTANT)
        ("etopo", setup::BathymetryType::ETOPO);
    M_bathymetry_type = str2bathymetry.find(vm["setup.bathymetry-type"].as<std::string>())->second;

    const boost::unordered_map<const std::string, setup::BasalStressType> str2basal_stress= boost::assign::map_list_of
        ("lemieux", setup::BasalStressType::LEMIEUX)
        ("bouillon", setup::BasalStressType::BOUILLON);
    M_basal_stress_type = str2basal_stress.find(vm["simul.basal_stress-type"].as<std::string>())->second;

    M_use_iabp_drifters=vm["simul.use_iabp_drifters"].as<bool>();
    M_equallyspaced_drifters_output_time_step=vm["simul.equallyspaced_drifters_output_time_step"].as<double>();
    M_rgps_drifters_output_time_step=vm["simul.rgps_drifters_output_time_step"].as<double>();
    M_use_osisaf_drifters=vm["simul.use_osisaf_drifters"].as<bool>();

    M_use_equallyspaced_drifters=false;
    M_use_rgps_drifters=false;

    if(M_equallyspaced_drifters_output_time_step>0.)
        M_use_equallyspaced_drifters=true;
    if(M_rgps_drifters_output_time_step>0.)
        M_use_rgps_drifters=true;


    // option for enabling/disabling the moorings
    M_use_moorings =  vm["simul.use_moorings"].as<bool>();
    M_moorings_snapshot =  vm["simul.mooring_snapshot"].as<bool>();
    const boost::unordered_map<const std::string, GridOutput::fileLength> str2mooringsfl = boost::assign::map_list_of
        ("inf", GridOutput::fileLength::inf)
        ("daily", GridOutput::fileLength::daily)
        ("weekly", GridOutput::fileLength::weekly)
        ("monthly", GridOutput::fileLength::monthly)
        ("yearly", GridOutput::fileLength::yearly);
    M_moorings_file_length = str2mooringsfl.find(vm["simul.mooring_file_length"].as<std::string>())->second;

    const boost::unordered_map<const std::string, mesh::Partitioner> str2partitioner = boost::assign::map_list_of
        ("chaco", mesh::Partitioner::CHACO)
        ("metis", mesh::Partitioner::METIS);
    M_partitioner = str2partitioner.find(vm["mesh.partitioner"].as<std::string>())->second;

    const boost::unordered_map<const std::string, mesh::PartitionSpace> str2partitionspace = boost::assign::map_list_of
        ("memory", mesh::PartitionSpace::MEMORY)
        ("disk", mesh::PartitionSpace::DISK);
    M_partition_space = str2partitionspace.find(vm["mesh.partition-space"].as<std::string>())->second;

    const boost::unordered_map<const std::string, LogLevel> str2log = boost::assign::map_list_of
        ("info", INFO)
        ("warning", WARNING)
        ("debug", DEBUG)
        ("error", ERROR);
    M_log_level = str2log.find(vm["simul.log-level"].as<std::string>())->second;

    M_mesh.setOrdering("bamg");

    M_mesh_filename = vm["mesh.filename"].as<std::string>();

    if (M_mesh_filename.find("split") != std::string::npos)
    {
        M_domain_type = setup::DomainType::DEFAULT;
        M_mesh_type = setup::MeshType::FROM_SPLIT;
    }
    else
    {
        M_domain_type = setup::DomainType::BIGARCTIC;
        M_mesh_type = setup::MeshType::FROM_GMSH;
    }

    M_mesh_fileformat = vm["mesh.fileformat"].as<std::string>();
}

void
FiniteElement::createGMSHMesh(std::string const& geofilename)
{
    std::string gmshgeofile = Environment::nextsimDir().string() + "/mesh/" + geofilename;

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
}

double
FiniteElement::jacobian(element_type const& element, mesh_type const& mesh) const
{
    std::vector<double> vertex_0 = mesh.nodes().find(element.indices[0])->second.coords;
    std::vector<double> vertex_1 = mesh.nodes().find(element.indices[1])->second.coords;
    std::vector<double> vertex_2 = mesh.nodes().find(element.indices[2])->second.coords;

    double jac = (vertex_1[0]-vertex_0[0])*(vertex_2[1]-vertex_0[1]);
    jac -= (vertex_2[0]-vertex_0[0])*(vertex_1[1]-vertex_0[1]);

    return  jac;
}

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
}

double
FiniteElement::jacobian(element_type const& element, mesh_type_root const& mesh) const
{
    std::vector<double> vertex_0 = mesh.nodes()[element.indices[0]-1].coords;
    std::vector<double> vertex_1 = mesh.nodes()[element.indices[1]-1].coords;
    std::vector<double> vertex_2 = mesh.nodes()[element.indices[2]-1].coords;

    double jac = (vertex_1[0]-vertex_0[0])*(vertex_2[1]-vertex_0[1]);
    jac -= (vertex_2[0]-vertex_0[0])*(vertex_1[1]-vertex_0[1]);

    return  jac;
}

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
}

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
}

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
}

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
}

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
}

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
}

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
}

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
}


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
}

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
}

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
}

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
}


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
        hmin[i] = std::sqrt(2.)*std::sqrt(*std::min_element(measure.begin(),measure.end()))*0.9;
    }

    return hmin;
}

std::vector<double>
FiniteElement::hmaxVertices(mesh_type_root const& mesh, BamgMesh const* bamg_mesh) const
{
    std::vector<double> hmax = this->hminVertices(mesh,bamg_mesh);

    std::for_each(hmax.begin(), hmax.end(), [&](double& f){ f = 1.2*f; });

    return hmax;
}

template<typename FEMeshType>
double
FiniteElement::measure(element_type const& element, FEMeshType const& mesh) const
{
    return (1./2)*std::abs(jacobian(element,mesh));
}

template<typename FEMeshType>
double
FiniteElement::measure(element_type const& element, FEMeshType const& mesh,
                       std::vector<double> const& um, double factor) const
{
    return (1./2)*std::abs(jacobian(element,mesh,um,factor));
}

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
}

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
}

void
FiniteElement::interpVertices()
{
    chrono.restart();
    LOG(DEBUG) <<"Interpolate hminVertices starts\n";
    // Interpolate hminVertices and hmaxVertices onto the current mesh

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
}

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
        int out_size = std::accumulate(sizes_nodes.begin(),sizes_nodes.end(),0);
        M_id_nodes.resize(out_size);

        boost::mpi::gatherv(M_comm, M_mesh.localDofWithGhostInit(), &M_id_nodes[0], sizes_nodes, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, M_mesh.localDofWithGhostInit(), 0);
    }

    // -------------------------------------------------------------
    if (M_rank == 0)
    {
        M_rmap_nodes = M_mesh.mapNodes();
        M_rmap_elements = M_mesh.mapElements();
    }
    // -------------------------------------------------------------
}

void
FiniteElement::collectVariables(std::vector<double>& interp_elt_in_local, bool slab, bool ghosts)
{
    int num_elements = M_local_nelements;
    if (ghosts)
    {
        num_elements = M_num_elements;
    }

    interp_elt_in_local.resize(M_nb_var_element*num_elements);
    M_interp_method.resize(M_nb_var_element);
    M_diffusivity_parameters.resize(M_nb_var_element);

    int tmp_nb_var=0;
    for (int i=0; i<num_elements; ++i)
    {
        tmp_nb_var=0;

        // concentration
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_conc[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // thickness
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_thick[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // snow thickness
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_snow_thick[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // integrated_stress1
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sigma[3*i]*M_thick[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // integrated_stress2
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sigma[3*i+1]*M_thick[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // integrated_stress3
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sigma[3*i+2]*M_thick[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // damage
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_damage[i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // ridge_ratio
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_ridge_ratio[i]*M_thick[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // random_number
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_random_number[i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // diffusivity_sss
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sss[i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=vm["simul.diffusivity_sss"].as<double>();
        tmp_nb_var++;

        // diffusivity_sst
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sst[i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=vm["simul.diffusivity_sst"].as<double>();
        tmp_nb_var++;

        // Ice temperature
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_tice[0][i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = ( M_tice[1][i] - physical::mu*physical::si*physical::Lf/(physical::C*M_tice[1][i]) ) * M_thick[i]; // (39) times volume with f1=1
            M_interp_method[tmp_nb_var] = 1;
            M_diffusivity_parameters[tmp_nb_var]=0.;
            tmp_nb_var++;

            interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = ( M_tice[2][i] ) * M_thick[i]; // (39) times volume with f1=0
            M_interp_method[tmp_nb_var] = 1;
            M_diffusivity_parameters[tmp_nb_var]=0.;
            tmp_nb_var++;
        }

        // thin ice thickness
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_h_thin[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // thin ice thickness
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_conc_thin[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // snow on thin ice
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_hs_thin[i];
        M_interp_method[tmp_nb_var] = 1;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        // Ice surface temperature for thin ice
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_tsurf_thin[i];
        M_interp_method[tmp_nb_var] = 0;
        M_diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        if(tmp_nb_var>M_nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }
    }
}

void
FiniteElement::redistributeVariables(std::vector<double> const& out_elt_values, bool slab)
{
    int tmp_nb_var=0;

    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_nb_var=0;

        // concentration
        //M_conc[i] = std::max(0., std::min(1.,out_elt_values[M_nb_var_element*i+tmp_nb_var]));
        M_conc[i] = std::max(0., out_elt_values[M_nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        // thickness
        M_thick[i] = std::max(0., out_elt_values[M_nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        // snow thickness
        M_snow_thick[i] = std::max(0., out_elt_values[M_nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        if (M_thick[i] != 0.)
        {
            // integrated_stress1
            M_sigma[3*i] = out_elt_values[M_nb_var_element*i+tmp_nb_var]/M_thick[i];
            tmp_nb_var++;

            // integrated_stress2
            M_sigma[3*i+1] = out_elt_values[M_nb_var_element*i+tmp_nb_var]/M_thick[i];
            tmp_nb_var++;

            // integrated_stress3
            M_sigma[3*i+2] = out_elt_values[M_nb_var_element*i+tmp_nb_var]/M_thick[i];
            tmp_nb_var++;

            // damage
            M_damage[i] = std::max(0., std::min(1.,out_elt_values[M_nb_var_element*i+tmp_nb_var]));
            tmp_nb_var++;

            // ridge_ratio
            M_ridge_ratio[i] = std::max(0., std::min(1.,out_elt_values[M_nb_var_element*i+tmp_nb_var]/M_thick[i]));
            tmp_nb_var++;
        }
        else
        {
            tmp_nb_var += 3;

            // damage
            M_damage[i] = 0.;
            tmp_nb_var++;

            // ridge_ratio
            M_ridge_ratio[i] = 0;
            tmp_nb_var++;
        }

        // random_number
        M_random_number[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        //M_random_number[i] = std::max(0., std::min(1.,out_elt_values[M_nb_var_element*i+tmp_nb_var]));
        tmp_nb_var++;

        // SSS
        M_sss[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // SST
        M_sst[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // Ice temperature
        M_tice[0][i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            if(M_thick[i]>0.)
            {
                double tmp = out_elt_values[M_nb_var_element*i+tmp_nb_var]/M_thick[i];
                M_tice[1][i] = 0.5*( tmp - std::sqrt(tmp*tmp + 4*physical::mu*physical::si*physical::Lf/physical::C) ); // (38) divided with volume with f1=1
                tmp_nb_var++;

                M_tice[2][i] = out_elt_values[M_nb_var_element*i+tmp_nb_var]/M_thick[i]; // (40) divided with volume with f1=0
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
        M_h_thin[i] = std::max(0., out_elt_values[M_nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        // thin ice concentration
        M_conc_thin[i] = std::max(0., out_elt_values[M_nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        // snow on thin ice
        M_hs_thin[i] = std::max(0., out_elt_values[M_nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        // Ice surface temperature for thin ice
        M_tsurf_thin[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

#if 0
        // slab ocean (Msst and M_sss)
        if (slab)
        {
            // Sea surface temperature
            M_sst[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
            tmp_nb_var++;

            // Sea surface salinity
            M_sss[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
            tmp_nb_var++;
        }
#endif

        if(tmp_nb_var!=M_nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }
    }
}

void
FiniteElement::advect(std::vector<double> const& interp_elt_in, std::vector<double>& interp_elt_out)
{
    M_comm.barrier();

    int ALE_smoothing_step_nb = vm["simul.ALE_smoothing_step_nb"].as<int>();
    // ALE_smoothing_step_nb<0 is the diffusive eulerian case where M_UM is not changed and then =0.
    // ALE_smoothing_step_nb=0 is the purely Lagrangian case where M_UM is updated with M_VT
    // ALE_smoothing_step_nb>0 is the ALE case where M_UM is updated with a smoothed version of M_VT

    std::vector<double> vt_root;
    std::vector<double> M_VT_smoothed;
    std::vector<double> M_VT_smoothed_root;
    this->gatherNodalField(M_VT,vt_root);

    // get the global number of nodes
    int num_nodes = M_mesh_root.numNodes();

    interp_elt_out.resize(M_nb_var_element*M_num_elements);

    if((ALE_smoothing_step_nb>=0) && (M_rank == 0))
    {
        M_VT_smoothed_root = vt_root;
        std::vector<double> M_VT_tmp = M_VT_smoothed_root;
        int Nd = bamgmesh_root->NodalConnectivitySize[1];

        for (int k=0; k<ALE_smoothing_step_nb; ++k)
        {
            M_VT_tmp = M_VT_smoothed_root;

            for (int i=0; i<M_ndof; ++i)
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
            x_ind = (M_elements[cpt]).indices[i]-1;
            y_ind = (M_elements[cpt]).indices[i]-1+M_num_nodes;

            x[i] = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[0];
            y[i] = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[1];

            /* old and new positions of the mesh */
            x_new[i] = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[0]+M_UM[x_ind];
            y_new[i] = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[1]+M_UM[y_ind];
            x[i]     = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[0]+UM_P[x_ind];
            y[i]     = M_nodes.find((M_elements[cpt]).indices[i])->second.coords[1]+UM_P[y_ind];

            if(ALE_smoothing_step_nb==-2)
            {
                VC_x[i] =0.;
                VC_y[i] =0.;
            }
            else
            {
                VC_x[i] = M_VT[x_ind]-(M_UM[x_ind]-UM_P[x_ind])/time_step;
                VC_y[i] = M_VT[y_ind]-(M_UM[y_ind]-UM_P[y_ind])/time_step;
            }
        }

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
                surface = this->measure(M_elements[cpt],M_mesh, UM_P);
                outer_fluxes_area[i] = std::min(surface/time_step/3.,outer_fluxes_area[i]);
                fluxes_source_id[i] = cpt;
            }
            else
            {
                neighbour_double = bamgmesh->ElementConnectivity[cpt*3+i];
                neighbour_int = (int)bamgmesh->ElementConnectivity[cpt*3+i];

                if (!std::isnan(neighbour_double) && neighbour_int>0)
                {
                    surface = this->measure(M_elements[neighbour_int-1],M_mesh, UM_P);
                    outer_fluxes_area[i] = -std::min(surface/time_step/3.,-outer_fluxes_area[i]);
                    fluxes_source_id[i] = neighbour_int-1;
                }
                else // open boundary with incoming fluxes
                {
                    fluxes_source_id[i]=cpt;
                }
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
}

void
FiniteElement::gatherFieldsElement(std::vector<double>& interp_in_elements)
{
    //M_comm.barrier();

    timer["gather"].first.restart();

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------GATHER ELEMENT starts\n";

    // ELEMENT INTERPOLATION With Cavities
    std::vector<int> sizes_elements = M_sizes_elements;
    //M_nb_var_element = 13 + M_tice.size();
    //std::cout<<"------------------------------------------------------------------------------------M_nb_var_element= "<< M_nb_var_element <<"\n";
    std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = M_nb_var_element*f; });

    std::vector<double> interp_elt_in_local;
    this->collectVariables(interp_elt_in_local, true);

#if 0
    int prv_num_nodes = M_local_ndof;
    int prv_num_elements = M_local_nelements;

    std::vector<double> interp_elt_in_local(M_nb_var_element*prv_num_elements);

    //LOG(DEBUG) <<"ELEMENT: Interp starts\n";

    int tmp_nb_var=0;
    for (int i=0; i<prv_num_elements; ++i)
    {
        tmp_nb_var=0;

        // concentration
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_conc[i];
        tmp_nb_var++;

        // thickness
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_thick[i];
        tmp_nb_var++;

        // snow thickness
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_snow_thick[i];
        tmp_nb_var++;

        // integrated_stress1
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sigma[3*i]*M_thick[i];
        tmp_nb_var++;

        // integrated_stress2
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sigma[3*i+1]*M_thick[i];
        tmp_nb_var++;

        // integrated_stress3
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sigma[3*i+2]*M_thick[i];
        tmp_nb_var++;

        // compliance
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = 1./(1.-M_damage[i]);
        tmp_nb_var++;

#if 0
        // divergence_rate
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_divergence_rate[i];
        tmp_nb_var++;

        // h_ridged_thin_ice
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_h_ridged_thin_ice[i];
        tmp_nb_var++;

        // h_ridged_thick_ice
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_h_ridged_thick_ice[i];
        tmp_nb_var++;
#endif

        // random_number
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_random_number[i];
        tmp_nb_var++;

        // // Ice surface temperature
        // interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_tsurf[i];
        // tmp_nb_var++;

        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_tice[0][i];
        tmp_nb_var++;
        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = ( M_tice[1][i] - physical::mu*physical::si*physical::Lf/(physical::C*M_tice[1][i]) ) * M_thick[i]; // (39) times volume with f1=1
            tmp_nb_var++;
            interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = ( M_tice[2][i] ) * M_thick[i]; // (39) times volume with f1=0
            tmp_nb_var++;
        }

        // thin ice thickness
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_h_thin[i];
        tmp_nb_var++;

        // snow on thin ice
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_hs_thin[i];
        tmp_nb_var++;

        // Ice surface temperature for thin ice
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_tsurf_thin[i];
        tmp_nb_var++;

        // Sea surface temperature
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sst[i];
        tmp_nb_var++;

        // Sea surface salinity
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sss[i];
        tmp_nb_var++;

        if(tmp_nb_var>M_nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }
    }
#endif

    if (M_rank == 0)
    {
        interp_in_elements.resize(M_nb_var_element*M_mesh_previous_root.numTriangles());
        boost::mpi::gatherv(M_comm, interp_elt_in_local, &interp_in_elements[0], sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, interp_elt_in_local, 0);
    }

    if (M_rank == 0)
    {
        //auto rmap_elements = M_mesh.mapElements();
        auto interp_in_elements_nrd = interp_in_elements;

        for (int i=0; i<M_mesh_previous_root.numTriangles(); ++i)
        {
            //int ri = rmap_elements.left.find(i+1)->second-1;
            int ri = M_rmap_elements[i];

            for (int j=0; j<M_nb_var_element; ++j)
            {
                interp_in_elements[M_nb_var_element*i+j] = interp_in_elements_nrd[M_nb_var_element*ri+j];
            }
        }
    }

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------GATHER ELEMENT done in "<< timer["gather"].first.elapsed() <<"s\n";
}

void
FiniteElement::scatterFieldsElement(double* interp_elt_out)
{
    timer["scatter"].first.restart();

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------SCATTER ELEMENT starts\n";

    std::vector<int> sizes_elements = M_sizes_elements_with_ghost;
    std::vector<int> id_elements;
    int out_size;

    chrono.restart();
    if (M_rank == 0)
    {
        out_size = std::accumulate(sizes_elements.begin(),sizes_elements.end(),0);
        id_elements.resize(out_size);

        boost::mpi::gatherv(M_comm, M_mesh.trianglesIdWithGhost(), &id_elements[0], sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, M_mesh.trianglesIdWithGhost(), 0);
    }

    std::vector<double> in_elt_values;

    if (M_rank == 0)
    {
        in_elt_values.resize(M_nb_var_element*id_elements.size());

        for (int i=0; i<id_elements.size(); ++i)
        {
            //int ri = rmap_elements.right.find(id_elements[i])->second-1;
            int ri = id_elements[i]-1;

            for (int j=0; j<M_nb_var_element; ++j)
            {
                in_elt_values[M_nb_var_element*i+j] = interp_elt_out[M_nb_var_element*ri+j];
            }
        }
    }

    std::vector<double> out_elt_values(M_nb_var_element*M_num_elements);

    if (M_rank == 0)
    {
        std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = M_nb_var_element*f; });
        boost::mpi::scatterv(M_comm, in_elt_values, sizes_elements, &out_elt_values[0], 0);
    }
    else
    {
        boost::mpi::scatterv(M_comm, &out_elt_values[0], M_nb_var_element*M_num_elements, 0);
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

    this->redistributeVariables(out_elt_values, true);

#if 0
    int tmp_nb_var=0;

    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_nb_var=0;

        // concentration
        M_conc[i] = std::max(0., std::min(1.,out_elt_values[M_nb_var_element*i+tmp_nb_var]));
        tmp_nb_var++;

        // thickness
        M_thick[i] = std::max(0., out_elt_values[M_nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        // snow thickness
        M_snow_thick[i] = std::max(0., out_elt_values[M_nb_var_element*i+tmp_nb_var]);
        tmp_nb_var++;

        if (M_thick[i] != 0.)
        {
            // integrated_stress1
            M_sigma[3*i] = out_elt_values[M_nb_var_element*i+tmp_nb_var]/M_thick[i];
            tmp_nb_var++;

            // integrated_stress2
            M_sigma[3*i+1] = out_elt_values[M_nb_var_element*i+tmp_nb_var]/M_thick[i];
            tmp_nb_var++;

            // integrated_stress3
            M_sigma[3*i+2] = out_elt_values[M_nb_var_element*i+tmp_nb_var]/M_thick[i];
            tmp_nb_var++;
        }
        else
        {
            tmp_nb_var+=3;
        }

        // compliance
        if (out_elt_values[M_nb_var_element*i+tmp_nb_var] != 0.)
        {
            M_damage[i] = std::max(0., std::min(1.,1.-(1./out_elt_values[M_nb_var_element*i+tmp_nb_var])));
            tmp_nb_var++;
        }
        else
        {
            M_damage[i] = 0.;
            tmp_nb_var++;
        }

#if 0
        // divergence_rate
        M_divergence_rate[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // h_ridged_thin_ice
        M_h_ridged_thin_ice[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // h_ridged_thick_ice
        M_h_ridged_thick_ice[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;
#endif

        // random_number
        M_random_number[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        //M_random_number[i] = std::max(0., std::min(1.,out_elt_values[M_nb_var_element*i+tmp_nb_var]));
        tmp_nb_var++;

        // // Ice surface temperature
        // M_tsurf[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        // tmp_nb_var++;

        // Ice temperature
        M_tice[0][i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;
        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            double tmp = out_elt_values[M_nb_var_element*i+tmp_nb_var]/M_thick[i];
            M_tice[1][i] = 0.5*( tmp - std::sqrt(tmp*tmp + 4*physical::mu*physical::si*physical::Lf/physical::C) ); // (38) divided with volume with f1=1
            tmp_nb_var++;
            M_tice[2][i] = out_elt_values[M_nb_var_element*i+tmp_nb_var]/M_thick[i]; // (40) divided with volume with f1=0
            tmp_nb_var++;
        }

        // thin ice thickness
        M_h_thin[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // snow on thin ice
        M_hs_thin[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // Ice surface temperature for thin ice
        M_tsurf_thin[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // Sea surface temperature
        M_sst[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // Sea surface salinity
        M_sss[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        if(tmp_nb_var!=M_nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }
    }
#endif

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------SCATTER ELEMENT done in "<< timer["scatter"].first.elapsed() <<"s\n";
}

void
FiniteElement::interpFieldsElement()
{
    std::vector<double> interp_in_elements;

    timer["gather"].first.restart();
    this->gatherFieldsElement(interp_in_elements);
    if (M_rank == 0)
        std::cout<<"-------------------GATHER done in "<< timer["gather"].first.elapsed() <<"s\n";

    double* interp_elt_out;

    if (M_rank == 0)
    {
        std::vector<double> surface_previous(M_mesh_previous_root.numTriangles());
        std::vector<double> surface(M_mesh_root.numTriangles());

        int cpt = 0;
        for (auto it=M_mesh_previous_root.triangles().begin(), end=M_mesh_previous_root.triangles().end(); it!=end; ++it)
        {
            surface_previous[cpt] = this->measure(*it,M_mesh_previous_root);
            ++cpt;
        }

        cpt = 0;
        for (auto it=M_mesh_root.triangles().begin(), end=M_mesh_root.triangles().end(); it!=end; ++it)
        {
            surface[cpt] = this->measure(*it,M_mesh_root);
            ++cpt;
        }

        // The interpolation with the cavities still needs to be tested on a long run.
        // By default, we then use the non-conservative MeshToMesh interpolation

        timer["cavities"].first.restart();
        InterpFromMeshToMesh2dCavities(&interp_elt_out,&interp_in_elements[0],
                                       &M_interp_method[0], M_nb_var_element,
                                       &surface_previous[0], &surface[0], bamgmesh_previous, bamgmesh_root);

        if (M_rank == 0)
            std::cout<<"-------------------CAVITIES done in "<< timer["cavities"].first.elapsed() <<"s\n";

#if 0

        // std::cout<<"M_mesh_previous_root.indexTr().size()= "<< M_mesh_previous_root.indexTr().size() <<"\n";
        // std::cout<<"M_mesh_previous_root.numTriangles()  = "<< M_mesh_previous_root.numTriangles() <<"\n";

        // auto indextr_ = M_mesh_previous_root.indexTr();
        // std::cout<<"["<< M_rank <<"]: " <<"Min index= "<< *std::min_element(indextr_.begin(), indextr_.end()) <<"\n";
        // std::cout<<"["<< M_rank <<"]: " <<"Max index= "<< *std::max_element(indextr_.begin(), indextr_.end()) <<"\n";


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
                                &M_mesh_root.bcoordX()[0],&M_mesh_root.bcoordY()[0],M_mesh_root.numTriangles(),
                                false);

        // std::cout<<"InterpFromMeshToMesh2dx done in "<< chrono.elapsed() <<"\n";
#endif

    } // rank 0

    timer["distributed"].first.restart();
    this->distributedMeshProcessing();
    if (M_rank == 0)
        std::cout<<"-------------------DISTRIBUTED done in "<< timer["distributed"].first.elapsed() <<"s\n";

    timer["scatter"].first.restart();
    this->scatterFieldsElement(interp_elt_out);
    if (M_rank == 0)
        std::cout<<"-------------------SCATTER done in "<< timer["scatter"].first.elapsed() <<"s\n";

    if (M_rank == 0)
    {
        xDelete<double>(interp_elt_out);
    }
}

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
}

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
}

void
FiniteElement::interpFieldsNode(std::vector<int> const& rmap_nodes, std::vector<int> sizes_nodes)
{
    std::vector<double> interp_in_nodes;

    this->gatherFieldsNode(interp_in_nodes, rmap_nodes, sizes_nodes);

    double* interp_nd_out;

    if (M_rank == 0)
    {
        InterpFromMeshToMesh2dx(&interp_nd_out,
                                &M_mesh_previous_root.indexTr()[0],&M_mesh_previous_root.coordX()[0],&M_mesh_previous_root.coordY()[0],
                                M_mesh_previous_root.numNodes(),M_mesh_previous_root.numTriangles(),
                                &interp_in_nodes[0],
                                M_mesh_previous_root.numNodes(),M_nb_var_node,
                                &M_mesh_root.coordX()[0],&M_mesh_root.coordY()[0],M_mesh_root.numNodes(),
                                false);
    }


    this->scatterFieldsNode(interp_nd_out);

    if (M_rank == 0)
    {
        xDelete<double>(interp_nd_out);
    }

    //carefull here (merging)
    //M_comm.barrier();
}

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
        //auto rmap_nodes = M_mesh.mapNodes();
        int global_num_nodes = M_mesh.numGlobalNodes();

        auto field_root_nrd = field_root;

        for (int i=0; i<global_num_nodes; ++i)
        {
            //int ri =  rmap_nodes.left.find(i+1)->second-1;
            int ri =  M_rmap_nodes[i];

            field_root[i] = field_root_nrd[2*ri];
            field_root[i+global_num_nodes] = field_root_nrd[2*ri+1];
        }
    }
}

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
}

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
}
#endif

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
    double minang = 0.;

    std::vector<double> um_root;

    //this->gatherUM(um_root);
    this->gatherNodalField(M_UM,um_root);

    if (M_rank == 0)
    {
        chrono.restart();
        LOG(DEBUG) <<"Flip starts\n";

        while (flip || (minang<(vm["simul.regrid_angle"].as<double>())/10.))
        {
            ++substep;
            displacement_factor /= 2.;
            step_order++;

            flip = this->flip(M_mesh_root,um_root,displacement_factor);
            minang = this->minAngle(M_mesh,M_UM,displacement_factor,true);

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
            if(step && (vm["simul.regrid_output_flag"].as<bool>()))
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
            if(step && (vm["simul.regrid_output_flag"].as<bool>()))
            {
                std::string tmp_string2    = (boost::format( "after_adaptMesh_%1%_mesh_adapt_step_%2%_substep_%3%" )
                                              % step
                                              % mesh_adapt_step
                                              % substep ).str();

                this->exportResults(tmp_string2,true,false,false);
            }
#endif

            // save mesh (only root process)

#if 0
            std::string src_fname = Environment::nextsimDir().string() + "/mesh/" + M_mesh_filename;
            std::string desc_fname = Environment::nextsimDir().string() + "/mesh/" + "prev_" + M_mesh_filename;
            fs::copy_file(fs::path(src_fname), fs::path(desc_fname), fs::copy_option::overwrite_if_exists);

            src_fname = Environment::nextsimDir().string() + "/mesh/" + "seq_" + M_mesh_filename;
            desc_fname = Environment::nextsimDir().string() + "/mesh/" + "seq_prev_" + M_mesh_filename;
            fs::copy_file(fs::path(src_fname), fs::path(desc_fname), fs::copy_option::overwrite_if_exists);
#endif
            std::cout<<"------------------------------version       = "<< M_mesh_root.version() <<"\n";
            std::cout<<"------------------------------ordering      = "<< M_mesh_root.ordering() <<"\n";
            std::cout<<"------------------------------format        = "<< M_mesh_fileformat <<"\n";
            std::cout<<"------------------------------space         = "<< vm["mesh.partition-space"].as<std::string>() <<"\n";
            std::cout<<"------------------------------partitioner   = "<< vm["mesh.partitioner"].as<std::string>() <<"\n";

            // Environment::logMemoryUsage("before partitioning...");
            timer["savemesh"].first.restart();
            LOG(DEBUG) <<"Saving mesh starts\n";
            if (M_partition_space == mesh::PartitionSpace::MEMORY)
            {
                M_mesh_root.writeToGModel(M_mesh_filename);
            }
            else if (M_partition_space == mesh::PartitionSpace::DISK)
            {
                M_mesh_root.writeTofile(M_mesh_filename);
            }
            std::cout <<"Saving mesh done in "<< timer["savemesh"].first.elapsed() <<"s\n";

#if 0
            src_fname = Environment::nextsimDir().string() + "/mesh/" + M_mesh_filename;
            desc_fname = Environment::nextsimDir().string() + "/mesh/" + "seq_" + M_mesh_filename;
            fs::copy_file(fs::path(src_fname), fs::path(desc_fname), fs::copy_option::overwrite_if_exists);
#endif

            // partition the mesh on root process (rank 0)
            timer["meshpartition"].first.restart();
            LOG(DEBUG) <<"Partitioning mesh starts\n";
            M_mesh_root.partition(M_mesh_filename,M_partitioner,M_partition_space, M_mesh_fileformat);
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

    timer["felt"].first.restart();
    this->interpFieldsElement();
    if (M_rank == 0)
        std::cout <<"interpFieldsElement done in "<< timer["felt"].first.elapsed() <<"s\n";

    timer["fnd"].first.restart();
    this->interpFieldsNode(prv_rmap_nodes, sizes_nodes);
    if (M_rank == 0)
        std::cout <<"interpFieldsNode done in "<< timer["fnd"].first.elapsed() <<"s\n";

    // --------------------------------END-------------------------------

    if (M_rank == 0)
        std::cout <<"TIMER REGRIDDING= "<< timer["regrid"].first.elapsed() <<"s\n";

    this->assignVariables();

#if 0
    M_comm.barrier();

    if (step)
    {
        M_solver->clear();
        M_solver.reset();
        M_solver = solver_ptrtype(new solver_type());
    }
#endif
}

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
    }

    timer["bamgmesh"].first.restart();
    Bamgx(bamgmesh_root,bamggeom_root,bamgmesh_previous,bamggeom_previous,bamgopt_previous);
    std::cout <<"---BAMGMESH done in "<< timer["bamgmesh"].first.elapsed() <<"s\n";

    this->importBamg(bamgmesh_root);

    this->updateBoundaryFlags();
}

void
FiniteElement::updateBoundaryFlags()
{
    LOG(DEBUG) <<"CLOSED: FLAGS SIZE BEFORE= "<< M_dirichlet_flags_root.size() <<"\n";
    LOG(DEBUG) <<"OPEN  : FLAGS SIZE BEFORE= "<< M_neumann_flags_root.size() <<"\n";

    // update dirichlet nodes
    M_dirichlet_flags_root.resize(0);
    M_neumann_flags_root.resize(0);

    // get the global number of nodes
    int num_nodes = M_mesh_root.numNodes();

    // We mask out the boundary nodes
    M_mask_root.assign(num_nodes,false);
    M_mask_dirichlet_root.assign(num_nodes,false);

    for (int edg=0; edg<bamgmesh_root->EdgesSize[0]; ++edg)
    {
        if (bamgmesh_root->Edges[3*edg+2] == M_flag_fix)
        {
            M_dirichlet_flags_root.push_back(bamgmesh_root->Edges[3*edg]/*-1*/);

            // new addition for masking the dirichlet nodes
            M_mask_dirichlet_root[edg] = true;
        }
        else
        {
            M_neumann_flags_root.push_back(bamgmesh_root->Edges[3*edg]/*-1*/);
        }

        M_mask_root[edg] = true;
    }

    M_dirichlet_nodes_root.resize(2*(M_dirichlet_flags_root.size()));
    for (int i=0; i<M_dirichlet_flags_root.size(); ++i)
    {
        M_dirichlet_nodes_root[2*i] = M_dirichlet_flags_root[i];
        M_dirichlet_nodes_root[2*i+1] = M_dirichlet_flags_root[i]+num_nodes;
    }

    M_neumann_nodes_root.resize(2*(M_neumann_flags_root.size()));
    for (int i=0; i<M_neumann_flags_root.size(); ++i)
    {
        M_neumann_nodes_root[2*i] = M_neumann_flags_root[i];
        M_neumann_nodes_root[2*i+1] = M_neumann_flags_root[i]+num_nodes;
    }


    LOG(DEBUG) <<"CLOSED: FLAGS SIZE AFTER= "<< M_dirichlet_flags_root.size() <<"\n";
    LOG(DEBUG) <<"OPEN  : FLAGS SIZE AFTER= "<< M_neumann_flags_root.size() <<"\n";
}

void
FiniteElement::assemble(int pcpt)
{
    M_comm.barrier();
    LOG(DEBUG) << "Reinitialize matrix and vector to zero starts\n";
    M_matrix->zero();
    M_vector->zero();
    LOG(DEBUG) << "Reinitialize matrix and vector to zero done\n";

    double coef_V, coef_Voce, coef_Vair, coef_basal, coef_X, coef_Y, coef_C;
    double coef = 0;
    double coef_P = 0.;
    double mass_e = 0.;
    double surface_e = 0.;
    double g_ssh_e_x = 0.;
    double g_ssh_e = 0.;
    double g_ssh_e_y = 0.;
    double tmp_thick, tmp_conc;
    double coef_sigma;

    double welt_oce_ice = 0.;
    double welt_air_ice = 0.;
    double welt_ice = 0.;
    double welt_ssh = 0.;
    double element_ssh = 0.;
    double critical_h = 0.;

    double max_keel_height;
    double ice_to_keel_factor;
    double keel_height_estimate;
    double critical_h_mod;

    double norm_Voce_ice = 0.;
    double norm_Voce_ice_min = 0.;
    double norm_Vair_ice = 0.;
    double norm_Vair_ice_min = 0.;
    double norm_Vice = 0.;

    double Vcor_index_v, Vcor_index_u;

    double b0tj_sigma_hu = 0.;
    double b0tj_sigma_hv = 0.;
    double mloc = 0.;
    double dloc = 0.;

    double undamaged_time_relaxation_sigma= 0.;
    double exponent_relaxation_sigma = 0.;;

    double time_viscous = 0.;
    double multiplicator = 0.;


    double duu, dvu, duv, dvv;
    int index_u, index_v;

    int nind;

    //double epsilon_veloc_i;
    //double divergence_rate;

    std::vector<int> extended_dirichlet_nodes = M_dirichlet_nodes;

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

    // if (M_rank == 0)
    //     std::cout<<"Assembling starts\n";

    timer["assembly"].first.restart();
    LOG(DEBUG) <<"Loop starts...\n";
    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
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
        total_thickness =       (vm["simul.min_h"].as<double>()>total_thickness)        ? vm["simul.min_h"].as<double>() : total_thickness;
        total_concentration =   (vm["simul.min_c"].as<double>()>total_concentration)    ? vm["simul.min_c"].as<double>() : total_concentration;

        double coef_min = 10.;
        double coef_drag = 0.;


        coef = coef_min;
        mass_e = 0.;
        coef_C = 0.;
        coef_V = 0.;
        coef_X = 0.;
        coef_Y = 0.;
        coef_sigma = 0.;

        coef_Vair = 0.;
        coef_Voce = 0.;
        coef_basal = 0.;
        coef_P = 0.;

        mloc = 0.;
        dloc = 0.;

        b0tj_sigma_hu = 0.;
        b0tj_sigma_hv = 0.;

        undamaged_time_relaxation_sigma = vm["simul.undamaged_time_relaxation_sigma"].as<double>();
        exponent_relaxation_sigma = vm["simul.exponent_relaxation_sigma"].as<double>();

        time_viscous = undamaged_time_relaxation_sigma*std::pow(1.-M_damage[cpt],exponent_relaxation_sigma-1.);
        multiplicator = time_viscous/(time_viscous+time_step);

        norm_Voce_ice = 0.;
        norm_Voce_ice_min = 0.01; // minimum value to avoid 0 water drag term.

        norm_Vair_ice = 0.;
        norm_Vair_ice_min = 0.01; // minimum value to avoid 0 water drag term.

        norm_Vice = 0.;

        element_ssh = 0.;
        critical_h = 0.;
        max_keel_height=28; // [m] from "A comprehensive analysis of the morphology of first-year sea ice ridges"
        ice_to_keel_factor=19.28; // from "A comprehensive analysis of the morphology of first-year sea ice ridges"
        keel_height_estimate;
        critical_h_mod=0.;

        //if(tmp_conc > vm["simul.min_c"].as<double>())
        if( (M_conc[cpt] > vm["simul.min_c"].as<double>()) && (M_thick[cpt] > vm["simul.min_h"].as<double>()) )
        {
            /* Compute the value that only depends on the element */
            welt_ice = 0.;
            welt_ssh = 0.;
            nind = 0;

            for (int i=0; i<3; ++i)
            {
                nind = it->indices[i]-1;
                welt_ssh += M_ssh[nind];
            }

            element_ssh = welt_ssh/3.;

            if(M_conc[cpt]>vm["simul.min_c"].as<double>())
            {
                switch ( M_basal_stress_type )
                {
                case setup::BasalStressType::BOUILLON:
                    // Sylvain's grounding scheme
                    keel_height_estimate = ice_to_keel_factor*std::pow(M_thick[cpt]/M_conc[cpt],0.5);
                    keel_height_estimate = ( keel_height_estimate > max_keel_height ) ? max_keel_height : keel_height_estimate;
                    critical_h      = M_conc[cpt]*std::pow((M_element_depth[cpt]+element_ssh)/ice_to_keel_factor,2.);
                    critical_h_mod  = M_conc[cpt]*std::pow(keel_height_estimate/ice_to_keel_factor,2.);
                    break;
                case setup::BasalStressType::LEMIEUX:
                    // JF Lemieux's grounding
                    keel_height_estimate = vm["simul.Lemieux_basal_k1"].as<double>()*M_thick[cpt]/M_conc[cpt];
                    keel_height_estimate = ( keel_height_estimate > max_keel_height ) ? max_keel_height : keel_height_estimate;

                    critical_h      = M_conc[cpt]*(M_element_depth[cpt]+element_ssh)/(vm["simul.Lemieux_basal_k1"].as<double>());
                    critical_h_mod  = M_conc[cpt]*keel_height_estimate/(vm["simul.Lemieux_basal_k1"].as<double>());
                    break;
                }
            }

            if(young>0.) // EB rheology
            {
                coef = multiplicator*young*(1.-M_damage[cpt])*M_thick[cpt]*std::exp(ridging_exponent*(1.-M_conc[cpt]));
            }
            else // Linear viscous rheology where nominal viscosity is defined as -young*time_step
            {
                double norm_factor=vm["simul.cohesion_thickness_normalisation"].as<double>();
                double exponent=vm["simul.cohesion_thickness_exponent"].as<double>();
                double mult_factor = std::pow(M_thick[cpt]/norm_factor,exponent);
                coef = -young*M_thick[cpt]*mult_factor*std::exp(ridging_exponent*(1.-M_conc[cpt]));
            }

            coef = (coef<coef_min) ? coef_min : coef ;

            if (vm["simul.use_coriolis"].as<bool>())
            {
                mass_e = (rhoi*total_thickness + rhos*total_snow)/total_concentration;
            }
            else
            {
                mass_e = 0.;
            }

            /* compute the x and y derivative of g*ssh */
            g_ssh_e_x = 0.;
            g_ssh_e_y = 0.;
            g_ssh_e;
            for(int i=0; i<3; i++)
            {
                g_ssh_e = (physical::gravity)*M_ssh[it->indices[i]-1] /*g_ssh*/;   /* g*ssh at the node k of the element e */
                g_ssh_e_x += M_shape_coeff[cpt][i]*g_ssh_e; /* x derivative of g*ssh */
                g_ssh_e_y += M_shape_coeff[cpt][i+3]*g_ssh_e; /* y derivative of g*ssh */
            }

            coef_drag  = 1.;
            coef_C = mass_e*M_fcor[cpt];              /* for the Coriolis term */
            coef_V = mass_e/time_step;             /* for the inertial term */
            coef_X = - mass_e*g_ssh_e_x;              /* for the ocean slope */
            coef_Y = - mass_e*g_ssh_e_y;              /* for the ocean slope */
            coef_sigma = M_thick[cpt]*multiplicator;

        } // test assemble

        surface_e = M_surface[cpt];


        //double coef_sigma = tmp_thick*multiplicator;
        //coef_sigma = (tmp_conc > vm["simul.min_c"].as<double>()) ? (coef_sigma):0.;

        // coef_Voce = (vm["simul.lin_drag_coef_water"].as<double>()+(quad_drag_coef_water*norm_Voce_ice));
        // coef_Voce *= physical::rhow;

        // coef_Vair = (vm["simul.lin_drag_coef_air"].as<double>()+(quad_drag_coef_air*norm_Vair_ice));
        // coef_Vair *= (physical::rhoa);


        if (0)
        {
            std::cout<<"************************\n";
            std::cout<<"Coef_C    = "<< coef_C <<"\n";
            std::cout<<"Coef_V    = "<< coef_V <<"\n";
            std::cout<<"Coef_X    = "<< coef_X <<"\n";
            std::cout<<"Coef_Y    = "<< coef_Y <<"\n";
            std::cout<<"coef_Vair = "<< coef_Vair <<"\n";
            std::cout<<"coef_Voce = "<< coef_Voce <<"\n";
            std::cout<<"coef_basal= "<< coef_basal <<"\n";
        }

        std::vector<int> rindices(6); //new
        std::vector<int> cindices(6);
        int vs = 0;

        for (int s=0; s<3; ++s)
        {
            int index_u = it->indices[s]-1;

            if (!it->ghostNodes[s])
            {
                rindices[2*vs] = index_u;
                rindices[2*vs+1] = index_u+M_local_ndof;
                ++vs;
            }

            cindices[2*s] = index_u;
            cindices[2*s+1] = index_u+M_local_ndof_ghost;
        }

        rindices.resize(2*vs);

        int nlrow = rindices.size();
        int nlcol = cindices.size();

        std::vector<double> data(nlrow*nlcol,0.);
        std::vector<double> fvdata(nlcol,0.);


        int l_j = -1;

        /* Loop over the 6 by 6 components of the finite element intergrale
         * this is done smartly by looping over j=0:2 and i=0:2
         * col = (mwIndex)it[2*j]-1  , row = (mwIndex)it[2*i]-1;
         * col  , row   -> UU component
         * col  , row+1 -> VU component
         * col+1, row   -> VV component
         * col+1, row+1 -> UV component */

        for(int j=0; j<3; j++)
        {
            /* Column corresponding to indice j (we also assemble terms in col+1) */
            //col = (mwIndex)it[2*j]-1; /* -1 to use the indice convention of C */

            index_u = it->indices[j]-1;
            index_v = it->indices[j]-1+M_num_nodes;

#if 0
            Vcor_index_v=beta0*M_VT[index_v] + beta1*M_VTM[index_v] + beta2*M_VTMM[index_v];
            Vcor_index_u=beta0*M_VT[index_u] + beta1*M_VTM[index_u] + beta2*M_VTMM[index_u];

            norm_Voce_ice = std::hypot(M_VT[index_u]-M_ocean[index_u],M_VT[index_v]-M_ocean[index_v]);
            norm_Voce_ice = (norm_Voce_ice > norm_Voce_ice_min) ? (norm_Voce_ice):norm_Voce_ice_min;

            norm_Vair_ice = std::hypot(M_VT[index_u]-M_wind [index_u],M_VT[index_v]-M_wind [index_v]);
            norm_Vair_ice = (norm_Vair_ice > norm_Vair_ice_min) ? (norm_Vair_ice):norm_Vair_ice_min;

            norm_Vice = std::hypot(M_VT[index_u],M_VT[index_v]);
            norm_Vice = (norm_Vice > basal_u_0) ? (norm_Vice):basal_u_0;
#else
            double vt_u = M_VT[index_u];
            double vt_v = M_VT[index_v];

            double ocean_u = M_ocean[index_u];
            double ocean_v = M_ocean[index_v];

            double wind_u = M_wind[index_u];
            double wind_v = M_wind[index_v];


            Vcor_index_v = beta0*vt_v + beta1*M_VTM[index_v] + beta2*M_VTMM[index_v];
            Vcor_index_u = beta0*vt_u + beta1*M_VTM[index_u] + beta2*M_VTMM[index_u];

            // norm_Voce_ice = std::hypot(vt_u-ocean_u,vt_v-ocean_v);
            // norm_Voce_ice = (norm_Voce_ice > norm_Voce_ice_min) ? (norm_Voce_ice):norm_Voce_ice_min;

            // norm_Vair_ice = std::hypot(vt_u-wind_u,vt_v-wind_v);
            // norm_Vair_ice = (norm_Vair_ice > norm_Vair_ice_min) ? (norm_Vair_ice):norm_Vair_ice_min;

            // norm_Vice = std::hypot(vt_u,vt_v);
            // norm_Vice = (norm_Vice > basal_u_0) ? (norm_Vice):basal_u_0;
#endif

            // coef_basal = basal_k2/norm_Vice;
            // coef_basal *= std::max(0., M_thick[cpt]-critical_h)*std::exp(-basal_Cb*(1.-M_conc[cpt]));

            // if(tmp_conc == vm["simul.min_c"].as<double>())
            // {
            //     coef_Voce = 0.;
            //     coef_Vair = 0.;
            //     coef_basal = 0.;
            // }


            if (!it->ghostNodes[j])
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

                    // --------------------------------------------------------------------------------------------------------
                    norm_Voce_ice = std::hypot(M_VT[index_u]-M_ocean[index_u],M_VT[index_v]-M_ocean[index_v]);
                    norm_Voce_ice = (norm_Voce_ice > norm_Voce_ice_min) ? (norm_Voce_ice):norm_Voce_ice_min;

                    coef_Voce = (vm["simul.lin_drag_coef_water"].as<double>()+(quad_drag_coef_water*norm_Voce_ice));
                    coef_Voce *= coef_drag*physical::rhow; //(vm["simul.rho_water"].as<double>());

                    norm_Vair_ice = std::hypot(M_VT[index_u]-M_wind [index_u],M_VT[index_v]-M_wind [index_v]);
                    norm_Vair_ice = (norm_Vair_ice > norm_Vair_ice_min) ? (norm_Vair_ice):norm_Vair_ice_min;

                    coef_Vair = (vm["simul.lin_drag_coef_air"].as<double>()+(quad_drag_coef_air*norm_Vair_ice));
                    coef_Vair *= coef_drag*(physical::rhoa);

                    norm_Vice = std::hypot(M_VT[index_u],M_VT[index_v]);
                    norm_Vice = (norm_Vice > basal_u_0) ? (norm_Vice):basal_u_0;

                    coef_basal = basal_k2/norm_Vice;
                    coef_basal *= coef_drag*std::max(0., critical_h_mod-critical_h)*std::exp(-basal_Cb*(1.-M_conc[cpt]));
                    // --------------------------------------------------------------------------------------------------------

                    /* ---------- UU component */
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

#if 1
                    data[(2*l_j  )*6+2*i  ] = duu;
                    data[(2*l_j+1)*6+2*i  ] = dvu;
                    data[(2*l_j  )*6+2*i+1] = duv;
                    data[(2*l_j+1)*6+2*i+1] = dvv;

                    fvdata[2*i] += surface_e*( mloc*(
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
#endif
                }
            }
        }

        // update matrix
        M_matrix->addMatrix(&rindices[0], rindices.size(), &cindices[0], cindices.size(), &data[0]);

        // update vector
        M_vector->addVector(&cindices[0], cindices.size(), &fvdata[0]);

        if ((cpt < 0))
        {
            std::cout<<"-------------------------------------------\n";

            for (int i=0; i<nlrow; ++i)
            {
                for (int j=0; j<nlcol; ++j)
                {
                    std::cout<< std::left << std::setw(12) << data[nlcol*i+j] <<"  ";
                }
                std::cout<<"\n";
            }
        }
        ++cpt;
    }

    // close petsc matrix
    LOG(DEBUG) <<"Closing matrix starts\n";
    M_matrix->close();
    LOG(DEBUG) <<"Closing matrix done\n";

#if 0
    if (pcpt == 0)
    {
        M_matrix->sameNonZeroPattern();
    }
#endif

    // close petsc vector
    LOG(DEBUG) <<"Closing vector starts\n";
    M_vector->close();
    LOG(DEBUG) <<"Closing vector done\n";

    if (M_rank == 0)
    {
        //std::cout<<"TIMER ASSEMBLY= " << timer["assembly"].first.elapsed() <<"s\n";
    }
    //std::cout<<"[" << M_rank <<"] " <<" TIMER ASSEMBLY= " << chrono.elapsed() <<"s\n";

    // extended dirichlet nodes (add nodes where M_conc <= 0)
    timer["dirichlet"].first.restart();

#if 0
    for (int i=0; i<M_local_ndof; ++i)
    {
        if(!M_valid_conc[i])
        {
            extended_dirichlet_nodes.push_back(i);
            extended_dirichlet_nodes.push_back(i+M_num_nodes);
        }
    }
    std::sort(extended_dirichlet_nodes.begin(), extended_dirichlet_nodes.end());
    extended_dirichlet_nodes.erase(std::unique( extended_dirichlet_nodes.begin(), extended_dirichlet_nodes.end() ), extended_dirichlet_nodes.end());

    std::cout<<"["<< M_rank <<"] EXTENDED = "<< extended_dirichlet_nodes.size() <<"\n";
    std::cout<<"["<< M_rank <<"] DIRICHLET= "<< M_dirichlet_nodes.size() <<"\n";
#endif

    M_matrix->on(extended_dirichlet_nodes,*M_vector);

    // if (M_rank==0)
    //     std::cout <<"TIMER DBCA= " << timer["dirichlet"].first.elapsed() <<"s\n";

#if 0
    int S1 = M_matrix->size1();
    int S2 = M_matrix->size2();
    double _norm = M_matrix->linftyNorm();
    double _l2norm = M_vector->l2Norm();

    if (M_rank == 0)
    {
        std::cout<<"[PETSC MATRIX] CLOSED          = "<< M_matrix->closed() <<"\n";
        //std::cout<<"[PETSC MATRIX] SIZE        = "<< M_matrix->size1() << " " << M_matrix->size2() <<"\n";
        std::cout<<"[PETSC MATRIX] SIZE            = "<< S1 << " " << S2 <<"\n";
        //std::cout<<"[PETSC MATRIX] SYMMETRIC   = "<< M_matrix->isSymmetric() <<"\n";
        //std::cout<<"[PETSC MATRIX] NORM        = "<< M_matrix->linftyNorm() <<"\n";
        std::cout<<"[PETSC MATRIX] NORM            = "<< _norm <<"\n";
        std::cout<<"[PETSC VECTOR] NORM            = "<< _l2norm <<"\n";
    }
#endif

    //M_matrix->printMatlab("stiffness.m");
    //M_vector->printMatlab("rhs.m");
}

void
FiniteElement::tensors()
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
        {
            for(int j=0; j<i; j++)
            {
                B0T_Dunit_B0T[i*6+j] = B0T_Dunit_B0T[j*6+i];
            }
        }

        M_shape_coeff[cpt]        = shapecoeff;
        M_B0T[cpt]                = B0T;
        M_B0T_Dunit_B0T[cpt]      = B0T_Dunit_B0T;

        ++cpt;
    }
}

void
FiniteElement::cohesion()
{
    for (int i=0; i<M_Cohesion.size(); ++i)
        M_Cohesion[i] = C_fix+C_alea*(M_random_number[i]-0.5);

    for (int i=0; i<M_Compressive_strength.size(); ++i)
        M_Compressive_strength[i] = compr_strength*scale_coef;
}

void
FiniteElement::update()
{
    // collect the variables into a single structure
    std::vector<double> interp_elt_in_local;
    this->collectVariables(interp_elt_in_local, false, true);

    // advect
    std::vector<double> interp_elt_out;
    this->advect(interp_elt_in_local, interp_elt_out);

    // redistribute the interpolated values
    this->redistributeVariables(interp_elt_out);

    double old_damage;
    /* deformation, deformation rate and internal stress tensor and temporary variables */
    double epsilon_veloc_i;
    std::vector<double> epsilon_veloc(3);

    std::vector<double> sigma_pred(3);
    double sigma_dot_i;

    /* divergence */
    double divergence_rate;

    /* shear*/
    double shear_rate;

    // ridging scheme
    double delta_ridging;
    double G_star=0.15;
    double e_factor=2.;

    /* invariant of the internal stress tensor and some variables used for the damaging process*/
    double sigma_s, sigma_n, sigma_1, sigma_2;
    double tract_max, sigma_t, sigma_c, q;
    double tmp, sigma_target;


    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        // Temporary memory
        old_damage = M_damage[cpt];

        /*======================================================================
         * Diagnostic:
         * Elastic deformation and instantaneous deformation rate
         *======================================================================
         */

        /* Compute the elastic deformation and the instantaneous deformation rate */
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
         * Ridging scheme
         * After the advection the concentration can be higher than 1, meaning that ridging should have occured.
         *======================================================================
         */
        double open_water_concentration=1.-M_conc[cpt];

        /* Thin ice category */
        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            open_water_concentration-=M_conc_thin[cpt];
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
                if( (M_conc[cpt] > vm["simul.min_c"].as<double>()) && (M_thick[cpt] > vm["simul.min_h"].as<double>()) && (new_conc_thin < M_conc_thin[cpt] ))
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

                M_conc_thin[cpt]= new_conc_thin;
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




        /*======================================================================
         * Update the internal stress
         *======================================================================
         */

        if( (M_conc[cpt] > vm["simul.min_c"].as<double>()) && (M_thick[cpt] > vm["simul.min_h"].as<double>()) && (young>0.))
        {

            double damaging_exponent = ridging_exponent;
            double undamaged_time_relaxation_sigma=vm["simul.undamaged_time_relaxation_sigma"].as<double>();
            double exponent_relaxation_sigma=vm["simul.exponent_relaxation_sigma"].as<double>();

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
            sigma_pred[i] = (M_conc[cpt] > vm["simul.min_c"].as<double>()) ? (sigma_pred[i]):0.;

            M_sigma[3*cpt+i] = (M_sigma[3*cpt+i]+time_step*sigma_dot_i)*multiplicator;
            M_sigma[3*cpt+i] = (M_conc[cpt] > vm["simul.min_c"].as<double>()) ? (M_sigma[3*cpt+i]):0.;
        }

        /*======================================================================
         * Correct the internal stress and the damage
         *======================================================================
         */

        /* Compute the shear and normal stress, which are two invariants of the internal stress tensor */

        sigma_s = std::hypot((sigma_pred[0]-sigma_pred[1])/2.,sigma_pred[2]);
        sigma_n = -          (sigma_pred[0]+sigma_pred[1])/2.;

        sigma_1 = sigma_n+sigma_s; // max principal component following convention (positive sigma_n=pressure)
        sigma_2 = sigma_n-sigma_s; // max principal component following convention (positive sigma_n=pressure)

        double ridge_to_normal_cohesion_ratio=vm["simul.ridge_to_normal_cohesion_ratio"].as<double>();
        double norm_factor=vm["simul.cohesion_thickness_normalisation"].as<double>();
        double exponent=vm["simul.cohesion_thickness_exponent"].as<double>();

        double hi=0.;
        if(M_conc[cpt]>0.1)
        {
            hi = M_thick[cpt]/M_conc[cpt];
        }
        else
        {
            hi = M_thick[cpt]/0.1;
        }

        double mult_factor = std::pow(hi/norm_factor,exponent)*(1. + M_ridge_ratio[cpt]*(ridge_to_normal_cohesion_ratio-1.) );

        double effective_cohesion = mult_factor * M_Cohesion[cpt];
        double effective_compressive_strength = mult_factor * M_Compressive_strength[cpt];

        q = std::pow(std::pow(std::pow(tan_phi,2.)+1,.5)+tan_phi,2.);
        sigma_c=2.*effective_cohesion/(std::pow(std::pow(tan_phi,2.)+1,.5)-tan_phi);
        sigma_t=-sigma_c/q;

        /* minimum and maximum normal stress */
        tract_max=-tract_coef*effective_cohesion/tan_phi;

        /* Correction of the damage */
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

        /*
         * Diagnostic:
         * Recompute the internal stress
         */

        for(int i=0;i<3;i++)
        {
            if(old_damage<1.0)
            {
                M_sigma[3*cpt+i] = (1.-M_damage[cpt])/(1.-old_damage)*M_sigma[3*cpt+i] ;
            }
            else
            {
                M_sigma[3*cpt+i] = 0. ;
            }
        }

        }
        else // if M_conc or M_thick too low, set sigma to 0.
        {
            for(int i=0;i<3;i++)
            {
                M_sigma[3*cpt+i] = 0.;
                M_sigma[3*cpt+i] = 0.;
            }
        }

        /*======================================================================
         * Update:
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
        {
            tmp=0.;
        }

        M_damage[cpt]=((tmp>0.)?(tmp):(0.));
    }
}

#if 0
void
FiniteElement::update()
{
    M_comm.barrier();

    std::vector<double> UM_P = M_UM;

    for (int nd=0; nd<M_UM.size(); ++nd)
    {
        M_UM[nd] += time_step*M_VT[nd];
    }

    for (const int& nd : M_neumann_nodes)
    {
        M_UM[nd] = UM_P[nd];
    }

    //std::cout<<"[" << M_rank <<"] " <<" Divergence_Rate MIN= "<< *std::min_element(M_divergence_rate.begin(),M_divergence_rate.end()) <<"\n";
    //std::cout<<"[" << M_rank <<"] " <<" Divergence_Rate MAX= "<< *std::max_element(M_divergence_rate.begin(),M_divergence_rate.end()) <<"\n";

    double old_thick;
    double old_snow_thick;
    double old_conc;
    double old_damage;
    double old_h_ridged_thick_ice;

    double old_h_thin;
    double old_hs_thin;
    double old_h_ridged_thin_ice;

    /* deformation, deformation rate and internal stress tensor and temporary variables */
    double epsilon_veloc_i;
    std::vector<double> epsilon_veloc(3);
    std::vector<double> sigma_pred(3);
    double sigma_dot_i;

    /* divergence */
    double divergence_rate;

    /* some variables used for the advection*/
    double surface, surface_new, ar, ar_new;
    double ice_surface, ice_volume, snow_volume;
    double thin_ice_surface, thin_ice_volume, thin_snow_volume;
    double ridging_thin_ice, ridging_thick_ice, ridging_snow_thin_ice;
    double ridged_thin_ice_volume, ridged_thick_ice_volume;

    /* invariant of the internal stress tensor and some variables used for the damaging process*/
    double sigma_s, sigma_n, sigma_1, sigma_2;
    double tract_max, sigma_t, sigma_c, q;
    double tmp, sigma_target;

    /* some variables used for the ice redistribution*/
    double tanalpha, rtanalpha, del_v, del_c, del_vs, new_v_thin;

    /* Indices loop*/
    int i,j;

    bool to_be_updated;

    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        // std::cout<<"-------------------"<< cpt << "-------------------"<<"\n";
        if (0)//(M_rank == 5)
        {
            std::cout<<"--------------------------------------"<<"\n";
            std::cout<<"(M_elements[cpt]).rank                = "<< M_rank <<"\n";
            std::cout<<"(M_elements[cpt]).number              = "<< (M_elements[cpt]).number <<"\n";
            std::cout<<"(M_elements[cpt]).type                = "<< (M_elements[cpt]).type <<"\n";
            std::cout<<"(M_elements[cpt]).physical            = "<< (M_elements[cpt]).physical <<"\n";
            std::cout<<"(M_elements[cpt]).elementary          = "<< (M_elements[cpt]).elementary <<"\n";
            std::cout<<"(M_elements[cpt]).numPartitions       = "<< (M_elements[cpt]).numPartitions <<"\n";
            std::cout<<"(M_elements[cpt]).partition           = "<< (M_elements[cpt]).partition <<"\n";
            std::cout<<"(M_elements[cpt]).is_ghost            = "<< (M_elements[cpt]).is_ghost <<"\n";
            std::cout<<"(M_elements[cpt]).ghosts              = "<<"\n";

            for (int k=0; k<(M_elements[cpt]).ghosts.size(); ++k)
            {
                std::cout<<"                    ghosts["<< k <<"]= "<< (M_elements[cpt]).ghosts[k] <<"\n";
            }

            std::cout<<"(M_elements[cpt]).is_on_processor     = "<< (M_elements[cpt]).is_on_processor <<"\n";
            std::cout<<"(M_elements[cpt]).ghost_partition_id  = "<< (M_elements[cpt]).ghost_partition_id <<"\n";
        }



        /* set constants for the ice redistribution */
        tanalpha  = h_thin_max/c_thin_max;
        rtanalpha = 1./tanalpha;

        // beginning of the original code (without openMP)
        // Temporary memory
        old_thick = M_thick[cpt];
        old_snow_thick = M_snow_thick[cpt];
        old_conc = M_conc[cpt];
        old_damage = M_damage[cpt];
        //old_h_ridged_thick_ice = M_h_ridged_thick_ice[cpt];

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            old_h_thin = M_h_thin[cpt];
            old_hs_thin = M_hs_thin[cpt];
            //old_h_ridged_thin_ice = M_h_ridged_thin_ice[cpt];
        }

        /*======================================================================
         * Diagnostic:
         * Elastic deformation and instantaneous deformation rate
         *======================================================================
         */

        /* Compute the elastic deformation and the instantaneous deformation rate */
        for(i=0;i<3;i++)
        {
            epsilon_veloc_i = 0.0;
            for(j=0;j<3;j++)
            {
                /* deformation */
                epsilon_veloc_i += M_B0T[cpt][i*6 + 2*j]*M_VT[(M_elements[cpt]).indices[j]-1];
                epsilon_veloc_i += M_B0T[cpt][i*6 + 2*j + 1]*M_VT[(M_elements[cpt]).indices[j]-1+M_num_nodes];
            }

            epsilon_veloc[i] = epsilon_veloc_i;
        }

        divergence_rate[cpt]= (epsilon_veloc[0]+epsilon_veloc[1]);

        /*======================================================================
         * Update the internal stress
         *======================================================================
         */

        for(i=0;i<3;i++)
        {
            sigma_dot_i = 0.0;
            for(j=0;j<3;j++)
            {
                sigma_dot_i += std::exp(ridging_exponent*(1.-old_conc))*young*(1.-old_damage)*M_Dunit[3*i + j]*epsilon_veloc[j];
            }

            M_sigma[3*cpt+i] += time_step*sigma_dot_i;
            sigma_pred[i]    = M_sigma[3*cpt+i] + time_step*sigma_dot_i;
        }

        /*======================================================================
         * Correct the internal stress and the damage
         *======================================================================
         */

        /* Compute the shear and normal stress, which are two invariants of the internal stress tensor */

        sigma_s=std::hypot((sigma_pred[0]-sigma_pred[1])/2.,sigma_pred[2]);
        sigma_n=-          (sigma_pred[0]+sigma_pred[1])/2.;

        sigma_1 = sigma_n+sigma_s; // max principal component following convention (positive sigma_n=pressure)
        sigma_2 = sigma_n-sigma_s; // max principal component following convention (positive sigma_n=pressure)

        q=std::pow(std::pow(std::pow(tan_phi,2.)+1,.5)+tan_phi,2.);
        sigma_c=2.*M_Cohesion[cpt]/(std::pow(std::pow(tan_phi,2.)+1,.5)-tan_phi);

        sigma_t=-sigma_c/q;

        /* minimum and maximum normal stress */
        tract_max=-tract_coef*M_Cohesion[cpt]/tan_phi;

        /* Correction of the damage */
        if(sigma_n>M_Compressive_strength[cpt])
        {
            sigma_target=M_Compressive_strength[cpt];

            tmp=1.0-sigma_target/sigma_n*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }


        if(sigma_1-q*sigma_2>sigma_c)
        {
            sigma_target=sigma_c;

            tmp=1.0-sigma_target/(sigma_1-q*sigma_2)*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }

        if(sigma_n<tract_max)
        {
            sigma_target=tract_max;

            tmp=1.0-sigma_target/sigma_n*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }


        /*
         * Diagnostic:
         * Recompute the internal stress
         */
        for(i=0;i<3;i++)
        {
            if(old_damage<1.0)
            {
                M_sigma[3*cpt+i] = (1.-M_damage[cpt])/(1.-old_damage)*M_sigma[3*cpt+i] ;
            }
            else
            {
                M_sigma[3*cpt+i] = 0. ;
            }
        }


        /*======================================================================
         * Update:
         * Ice damage
         * We use now a constant healing rate defined as 1/time_recovery_damage
         * so that we are now able to reset the damage to 0.
         * otherwise, it will never heal completely.
         * time_recovery_damage still depends on the temperature when themodynamics is activated.
         *======================================================================
         */
        tmp=1./(1.-M_damage[cpt]);
        tmp-= 1000*time_step/M_time_relaxation_damage[cpt];
        tmp=((tmp>1.)?(tmp):(1.));
        M_damage[cpt]=-1./tmp + 1.;

        /*======================================================================
         * Update:
         * Ice and snow thickness, and concentration using a Lagrangian or an Eulerian scheme
         *======================================================================
         */

        to_be_updated=false;
        if( divergence_rate[cpt]!=0.)
            to_be_updated=true;


        if(std::binary_search(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[0]-1) ||
           std::binary_search(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[1]-1) ||
           std::binary_search(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[2]-1))
            to_be_updated=false;

        if((old_conc>0.)  && (to_be_updated))
        {
            surface = this->measure(M_elements[cpt],M_mesh, UM_P);
            surface_new = this->measure(M_elements[cpt],M_mesh,M_UM);

            ice_surface = old_conc*surface;
            ice_volume = old_thick*surface;
            snow_volume = old_snow_thick*surface;
            ridged_thick_ice_volume = old_h_ridged_thick_ice*surface;

            M_conc[cpt]    = ice_surface/surface_new;
            M_thick[cpt]   = ice_volume/surface_new; // Hold on! Isn't M_thick the effective thickness?
            M_snow_thick[cpt]   = snow_volume/surface_new;
            //M_h_ridged_thick_ice[cpt]   =   ridged_thick_ice_volume/surface_new;


            if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                thin_ice_volume = old_h_thin*surface;
                thin_snow_volume = old_hs_thin*surface;
                ridged_thin_ice_volume = old_h_ridged_thin_ice*surface;

                M_h_thin[cpt]        = thin_ice_volume/surface_new;
                M_hs_thin[cpt]   = thin_snow_volume/surface_new;
                //M_h_ridged_thin_ice[cpt]    =   ridged_thin_ice_volume/surface_new;
            }

            /* Ridging scheme */
            if(surface_new<surface)
            {
                if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
                {
                    ridging_thin_ice=(surface-surface_new)/surface_new*old_h_thin;
                    ridging_snow_thin_ice=(surface-surface_new)/surface_new*old_hs_thin;

                    ridging_thin_ice=((ridging_thin_ice<old_h_thin)?(ridging_thin_ice):(old_h_thin));
                    ridging_snow_thin_ice=((ridging_snow_thin_ice<old_hs_thin)?(ridging_snow_thin_ice):(old_hs_thin)) ;

                    M_thick[cpt] += ridging_thin_ice ;
                    M_h_thin[cpt] -= ridging_thin_ice ;

                    M_snow_thick[cpt] += ridging_snow_thin_ice ;
                    M_hs_thin[cpt] -= ridging_snow_thin_ice ;

                    //M_h_ridged_thin_ice[cpt] += ridging_thin_ice;
                    M_conc[cpt] += ridging_thin_ice/ridge_h;

                    /* upper bounds (only for the concentration) */
                    ridging_thin_ice = ((M_conc[cpt]<1.)?(0.):(M_h_thin[cpt])) ;
                    ridging_snow_thin_ice = ((M_conc[cpt]<1.)?(0.):(M_hs_thin[cpt])) ;

                    M_snow_thick[cpt] += ridging_snow_thin_ice ;
                    M_hs_thin[cpt] -= ridging_snow_thin_ice ;

                    M_thick[cpt] += ridging_thin_ice;
                    M_h_thin[cpt] -= ridging_thin_ice;
                }
                /* upper bounds (only for the concentration) */
                ridging_thick_ice=((M_conc[cpt]<1.)?(0.):(M_thick[cpt]*(M_conc[cpt]-1.)));
                M_conc[cpt] = ((M_conc[cpt]<1.)?(M_conc[cpt]):(1.)) ;
            }

            /* lower bounds */
            M_conc[cpt] = ((M_conc[cpt]>0.)?(M_conc[cpt] ):(0.)) ;
            M_thick[cpt]        = ((M_thick[cpt]>0.)?(M_thick[cpt]     ):(0.)) ;
            M_snow_thick[cpt]   = ((M_snow_thick[cpt]>0.)?(M_snow_thick[cpt]):(0.)) ;

            if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                M_h_thin[cpt]    = ((M_h_thin[cpt]>0.)?(M_h_thin[cpt] ):(0.)) ;
                M_hs_thin[cpt]   = ((M_hs_thin[cpt]>0.)?(M_hs_thin[cpt]):(0.)) ;
            }
        }

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            /* Compute the redistribution of thin ice. */
            /* Returns the change in volume and concentration of thick ice as well as the
             * change in volume of thin ice. It is called after the
             * dynamics are done. */

            if(M_h_thin[cpt]>0.)
            {
                thin_ice_redistribute(M_h_thin[cpt], M_hs_thin[cpt], 0., M_conc[cpt],
                                      tanalpha, rtanalpha, h_thin_max, &new_v_thin, &del_v, &del_c, &del_vs);

                M_conc[cpt]       += del_c;

                M_thick[cpt]           += del_v;
                M_h_thin[cpt]      -= del_v;

                M_snow_thick[cpt]      += del_vs;
                M_hs_thin[cpt] -= del_vs;
            }
            else
            {
                M_snow_thick[cpt] += M_hs_thin[cpt];
                M_hs_thin[cpt] = 0. ;
            }
        }
    }
}
#endif

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
                    _reuse_prec=vm["solver.ksp-reuse-prec"].as<bool>()/*false*/,
                    _rebuild=M_regrid
                    );


    M_solution->close();

    // if (M_rank==0)
    //     std::cout<<"TIMER SOLUTION= " << timer["solution"].first.elapsed() <<"s\n";

    //M_solution->printMatlab("solution.m");
}

// Routine for the 1D thermodynamical model
// No thin ice for now
// No stability dependent drag for now
void
FiniteElement::thermo()
{
    M_comm.barrier();

    // There is now only one big loop for the thermodynamics so that we can use multithreading.

    // constant variables
    // Set local variable to values defined by options
    double const timeT = vm["simul.ocean_nudge_timeT"].as<double>();
    double const timeS = vm["simul.ocean_nudge_timeS"].as<double>();
    double const Qdw_const = vm["simul.constant_Qdw"].as<double>();
    double const Fdw_const = vm["simul.constant_Fdw"].as<double>();

    double const ocean_albedo = vm["simul.albedoW"].as<double>();
    double const drag_ocean_t = vm["simul.drag_ocean_t"].as<double>();
    double const drag_ocean_q = vm["simul.drag_ocean_q"].as<double>();

    double const rh0   = 1./vm["simul.hnull"].as<double>();
    double const rPhiF = 1./vm["simul.PhiF"].as<double>();

    double const qi = physical::Lf * physical::rhoi;// J m^{-3}
    double const qs = physical::Lf * physical::rhos;// J m^{-3}

    int const newice_type = vm["simul.newice_type"].as<int>();
    int const melt_type = vm["simul.melt_type"].as<int>();
    double const PhiM = vm["simul.PhiM"].as<double>();
    double const PhiF = vm["simul.PhiF"].as<double>();

    const double aw=6.1121e2, bw=18.729, cw=257.87, dw=227.3;
    const double Aw=7.2e-4, Bw=3.20e-6, Cw=5.9e-10;

    const double alpha=0.62197, beta=0.37803;

    // initialisation of the multithreading
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    // #pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int i=0; i < M_num_elements; ++i)
    {
        // -------------------------------------------------
        // 1) Initialization of the temporary variables

        double  hi=0.;     // Ice thickness (slab)
        double  hi_old=0.; // Ice thickness at the start of the time step (slab)
        double  hs=0.;     // Snow thickness (slab)

        double  hi_thin=0.;     // Thin ice thickness (slab)
        double  hi_thin_old=0.; // Thin ice thickness at the start of the time step (slab)
        double  hs_thin=0.;     // Snow thickness on thin ice (slab)

        double  del_hi=0.; // Change in ice thickness (slab only)
        double  del_hi_thin=0.; // Change in thin ice thickness (slab only)

        double  evap=0.;   // Evaporation

        double  Qdw=0.;    // Heat flux from ocean nudging
        double  Fdw=0.;    // Fresh water flux from ocean nudging

        double  Qio=0.;    // Ice-ocean heat flux
        double  Qio_thin=0.;    // Ice-ocean heat flux through thin ice
        double  Qai=0.;    // Atmosphere-ice heat flux
        double  Qai_thin=0.;    // Atmosphere-ice heat flux over thin ice
        double  Qswi=0.;    // Atmosphere-ice heat flux
        double  Qsw_thin=0.;    // Atmosphere-ice heat flux over thin ice
        double  Qlwi=0.;    // Atmosphere-ice heat flux
        double  Qlw_thin=0.;    // Atmosphere-ice heat flux over thin ice
        double  Qshi=0.;    // Atmosphere-ice heat flux
        double  Qsh_thin=0.;    // Atmosphere-ice heat flux over thin ice
        double  Qlhi=0.;    // Atmosphere-ice heat flux
        double  Qlh_thin=0.;    // Atmosphere-ice heat flux over thin ice
        double  Qow=0.;    // Open water heat flux

        // Save old _volumes_ and concentration and calculate wind speed
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
        if(M_snowfr.M_initialized)
            tmp_snowfall=M_precip[i]*M_snowfr[i];
        else if (M_snowfall.M_initialized)
            tmp_snowfall=M_snowfall[i];
        else
        {
            if(M_tair[i]<0)
                tmp_snowfall=M_precip[i];
            else
                tmp_snowfall=0.;
        }

        double Qsw_in;
        if(M_Qsw_in.M_initialized)
            Qsw_in=M_Qsw_in[i];
        else
        {
            throw std::logic_error("The function approxSW not yet implemented, you need to initialized M_Qsw_in");
            //Qsw_in=approxSW();
        }

        double mld=( M_mld[i] > vm["simul.constant_mld"].as<double>() ) ? M_mld[i] : vm["simul.constant_mld"].as<double>();

        // -------------------------------------------------
        // 2) We calculate or set the flux due to nudging
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
                Fdw = delS * mld * physical::rhow /(timeS*M_sss[i] - time_step*delS);
            } else {
                Qdw = Qdw_const;
                Fdw = Fdw_const;
            }
        }

        // SYL: Calculate the drag coefficients missing??

        // -------------------------------------------------
        // 3) Calculate fluxes in the open water portion (openWaterFlux in old c code, Qow_mex in matlab)
        double sphuma, sphumw;

        // Calculate atmospheric fluxes

        /* Out-going long-wave flux */
        double Qlw_out = physical::eps*physical::sigma_sb*std::pow(M_sst[i]+physical::tfrwK,4.);

        // -------------------------------------------------
        // 3.1) Specific humidity - atmosphere (calcSphumA in matlab)
        /* There are two ways to calculate this. We decide which one by
         * checking mixrat - the calling routine must set this to a negative
         * value if the dewpoint should be used. */
        if ( M_sphuma.M_initialized )
            sphuma = M_sphuma[i];
        else if ( M_dair.M_initialized )
        {
            double fa     = 1. + Aw + M_mslp[i]*1e-2*( Bw + Cw*M_dair[i]*M_dair[i] );
            double esta   = fa*aw*std::exp( (bw-M_dair[i]/dw)*M_dair[i]/(M_dair[i]+cw) );
            sphuma = alpha*fa*esta/(M_mslp[i]-beta*fa*esta) ;
        }
        else if ( M_mixrat.M_initialized )
            sphuma = M_mixrat[i]/(1.+M_mixrat[i]) ;
        else
            throw std::logic_error("Neither M_dair nor M_mixrat have been initialized. I cannot calculate sphuma.");

        // -------------------------------------------------
        // 3.2) Specific humidity - ocean surface (calcSphumW in matlab)
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
        if(M_Qlw_in.M_initialized)
            Qlw_in=M_Qlw_in[i];
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
        // 4) Thickness change of the ice slab (thermoIce0 in matlab)

        switch ( M_thermo_type )
        {
            case setup::ThermoType::ZERO_LAYER:
                this->thermoIce0(i, wspeed, sphuma, M_conc[i], M_thick[i], M_snow_thick[i],
                        Qlw_in, Qsw_in, mld, tmp_snowfall, hi, hs, hi_old, Qio, del_hi, M_tice[0][i],
                        Qai, Qswi, Qlwi, Qshi, Qlhi);
                break;
            case setup::ThermoType::WINTON:
                this->thermoWinton(i, time_step, wspeed, sphuma, M_conc[i], M_thick[i], M_snow_thick[i],
                        Qlw_in, Qsw_in, mld, tmp_snowfall, hi, hs, hi_old, Qio, del_hi,
                        M_tice[0][i], M_tice[1][i], M_tice[2][i],
                        Qai, Qswi, Qlwi, Qshi, Qlhi);
                break;
        }

        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            this->thermoIce0(i, wspeed, sphuma, old_conc_thin, M_h_thin[i], M_hs_thin[i],
                    Qlw_in, Qsw_in, mld, tmp_snowfall, hi_thin, hs_thin, hi_thin_old, Qio_thin, del_hi_thin, M_tsurf_thin[i],
                        Qai_thin, Qsw_thin, Qlw_thin, Qsh_thin, Qlh_thin);
            M_h_thin[i]  = hi_thin * old_conc_thin;
            M_hs_thin[i] = hs_thin * old_conc_thin;
        }

        // -------------------------------------------------
        // 5) Ice growth over open water and lateral melt (thermoOW in matlab)

        /* Local variables */
        double tw_new, tfrw, newice, del_c, newsnow, h0;

        /* dT/dt due to heatflux ocean->atmosphere */
        tw_new = M_sst[i] - Qow*time_step/(mld*physical::rhow*physical::cpw);
        tfrw   = -physical::mu*M_sss[i];

        /* Form new ice in case of super cooling, and reset Qow and evap */
        if ( tw_new < tfrw )
        {
            newice  = (1.-M_conc[i]-M_conc_thin[i])*(tfrw-tw_new)*mld*physical::rhow*physical::cpw/qi;// m
            Qow  = -(tfrw-M_sst[i])*mld*physical::rhow*physical::cpw/time_step;
            // evap = 0.;
        } else {
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
                    } else {
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
                                M_conc_thin[i] = M_h_thin[i]/h_thin_min;
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
                        del_c += PhiM*(1.-M_conc[i])*std::min(0.,Qow)*time_step/( hi*qi+hs*qs );
                        /* Deliver the fraction (1-PhiM) of Qow to the ocean */
                        Qow = (1.-PhiM)*Qow;
                    } else
                        del_c = -M_conc[i];
                    // This is handled below
                    // /* + Deliver excess energy to the ocean when there's no ice left */
                    //         + std::min(0., std::max(0.,M_conc[i]+del_c)*( hi*qi+hs*qs )/time_step);
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
                Qow = Qow + del_c*hs*qs/time_step;
            } else {
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
        if ( M_conc[i] < physical::cmin || hi < physical::hmin )
        {
            // Extract heat from the ocean corresponding to the heat in the
            // remaining ice and snow
            Qow    = Qow + M_conc[i]*hi*qi/time_step + M_conc[i]*hs*qs/time_step;
            M_conc[i]  = 0.;
            for (int j=0; j<M_tice.size(); j++)
                M_tice[j][i] = tfrw;
            M_tsurf_thin[i] = tfrw;
            hi     = 0.;
            hs     = 0.;
        }
#endif
        // -------------------------------------------------
        // 6) Calculate effective ice and snow thickness
        M_thick[i] = hi*M_conc[i];
        M_snow_thick[i] = hs*M_conc[i];

        // -------------------------------------------------
        // 7) Slab Ocean (slabOcean in matlab)

        // local variables
        double del_vi;      // Change in ice volume
        double del_vs;      // Change in snow olume
        double rain;        // Liquid precipitation
        double emp;         // Evaporation minus liquid precipitation
        double Qio_mean;    // Element mean ice-ocean heat flux
        double Qow_mean;    // Element mean open water heat flux

        // Calculate change in volume to calculate salt rejection
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
        M_sst[i] = M_sst[i] - time_step*( Qio_mean + Qow_mean - Qdw )/(physical::rhow*physical::cpw*mld);

        /* Change in salinity */
        double denominator= ( mld*physical::rhow - del_vi*physical::rhoi - ( del_vs*physical::rhos + (emp-Fdw)*time_step) );
        denominator = ( denominator > 1.*physical::rhow ) ? denominator : 1.*physical::rhow;

        double sss_old = M_sss[i];
        M_sss[i] = M_sss[i] + ( (M_sss[i]-physical::si)*physical::rhoi*del_vi + M_sss[i]*(del_vs*physical::rhos + (emp-Fdw)*time_step) )
            / denominator;

        // -------------------------------------------------
        // 8) Damage manipulation (thermoDamage in matlab)

        // local variables
        double deltaT;      // Temperature difference between ice bottom and the snow-ice interface

        // Newly formed ice is undamaged - and unridged
        // calculate damage and ridge ratio as a weighted
        // average of the old damage - ridge ratio and 0, weighted with volume.
        if ( M_thick[i] > old_vol )
        {
            M_damage[i] = M_damage[i]*old_vol/M_thick[i];
            M_ridge_ratio[i] = M_ridge_ratio[i]*old_vol/M_thick[i];
        }

        if ( vm["simul.use_temperature_dependent_healing"].as<bool>() )
        {
            // Set time_relaxation_damage to be inversely proportional to
            // temperature difference between bottom and snow-ice interface
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
                M_time_relaxation_damage[i] = std::max(time_relaxation_damage*deltaT_relaxation_damage/deltaT, time_step);
            } else {
                M_time_relaxation_damage[i] = 1e36;
            }
        }
        // -------------------------------------------------

        // Diagnostics
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
        D_delS[i] = (M_sss[i] - sss_old)*physical::rhow*mld/time_step;

    }// end for loop
}// end thermo function

// Atmospheric fluxes through bulk formula
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
}

// Ice-ocean heat flux
double
FiniteElement::iceOceanHeatflux(int cpt, double sst, double sss, double mld, double dt)
{
    /* Use all excess heat to melt or grow ice. This is not
     * accurate, but will have to do for now! */
    double const Tbot = -physical::mu*sss; // Temperature at ice base (bottom), also freezing point of sea-water
    if ( vm["simul.Qio-type"].as<std::string>() == "basic" )
    {
        return (sst-Tbot)*physical::rhow*physical::cpw*mld/dt;
    } else if ( vm["simul.Qio-type"].as<std::string>() == "exchange" ) {
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
        std::cout << "Qio-type = " << vm["simul.Qio-type"].as<std::string>() << "\n";
        throw std::logic_error("Wrong Qio-type");
    }
}

// Albedo
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
}

// Winton thermo dynamics (ice temperature, growth, and melt)
void
FiniteElement::thermoWinton(int i, double dt, double wspeed, double sphuma, double conc, double voli, double vols,
        double Qlw_in, double Qsw_in, double mld, double snowfall,
        double &hi, double &hs, double &hi_old, double &Qio, double &del_hi, double &Tsurf, double &T1, double &T2,
        double &Qai, double &Qsw, double &Qlw, double &Qsh, double &Qlh)
{
    // Constants
    double const alb_ice = vm["simul.alb_ice"].as<double>();
    double const alb_sn  = vm["simul.alb_sn"].as<double>();
    double const I_0     = vm["simul.I_0"].as<double>();
    int const alb_scheme = vm["simul.alb_scheme"].as<int>();

    double const drag_ice_t = vm["simul.drag_ice_t"].as<double>();

    bool const flooding = vm["simul.flooding"].as<bool>();

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
        double Qsw = -Qsw_in*(1.-FiniteElement::albedo(alb_scheme, Tsurf, hs, alb_sn, alb_ice, I_0))*(1.-I_0);
        // The rest is calculated by bulk formula
        FiniteElement::atmFluxBulk(i, Tsurf, sphuma, drag_ice_t, Qsw, Qlw_in, wspeed, Qai, dQaidT,subl,
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
        Qio    = FiniteElement::iceOceanHeatflux(i, M_sst[i], M_sss[i], mld, dt);
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
}

// This is Semtner zero layer
void
FiniteElement::thermoIce0(int i, double wspeed, double sphuma, double conc, double voli, double vols, double Qlw_in, double Qsw_in, double mld, double snowfall,
        double &hi, double &hs, double &hi_old, double &Qio, double &del_hi, double &Tsurf,
        double &Qai, double &Qsw, double &Qlw, double &Qsh, double &Qlh)
{

    // Constants
    double const alb_ice = vm["simul.alb_ice"].as<double>();
    double const alb_sn  = vm["simul.alb_sn"].as<double>();
    double const I_0     = vm["simul.I_0"].as<double>();
    int const alb_scheme = vm["simul.alb_scheme"].as<int>();

    double const drag_ice_t = vm["simul.drag_ice_t"].as<double>();

    bool const flooding = vm["simul.flooding"].as<bool>();

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
        double Qsw, Qout, dQaidT, subl;
        double Qic, del_hs, del_ht, del_hb, draft;

        double Qlw_out, dQlwdT;
        double tairK, sphumi;
        double rhoair, Qsh, dQshdT;
        double Qlh, dsphumidT, dQlhdT;

        double fi, esti;
        double dsphumdesti, destidT, dfidT;

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
            FiniteElement::atmFluxBulk(i, Tsurf, sphuma, drag_ice_t, Qsw, Qlw_in, wspeed, Qai, dQaidT,subl,
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
        del_hs = std::min(Qai-Qic,0.)*time_step/qs - subl*time_step/physical::rhos;
        /* Use the energy left over after snow melts to melt the ice */
        del_ht = std::min(hs+del_hs,0.)*qs/qi;
        /* Can't have negative hs! */
        del_hs = std::max(del_hs,-hs);
        // snowfall in kg/m^2/s
        hs  = hs + del_hs + snowfall/physical::rhos*time_step;

        /* Heatflux from ocean */
        Qio = FiniteElement::iceOceanHeatflux(i, M_sst[i], M_sss[i], mld, time_step);
        /* Bottom melt/growth */
        del_hb = (Qic-Qio)*time_step/qi;

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
            Qio     = Qio + hi*qi/time_step + hs*qs/time_step;

            hi      = 0.;
            hs      = 0.;
            Tsurf   = Tbot;
        }
    }
}

void
FiniteElement::init()
{
    // define export path
    M_export_path = Environment::nextsimDir().string() + "/matlab";
    // change directory for outputs if the option "output_directory" is not empty
    if ( ! (vm["simul.output_directory"].as<std::string>()).empty() )
    {
        M_export_path = vm["simul.output_directory"].as<std::string>();

        fs::path path(M_export_path);
        // add a subdirecory if needed
        // path /= "subdir";

        // create the output directory if it does not exist
        if ( (!fs::exists(path)) && (M_comm.rank()==0) )
            fs::create_directories(path);
    }
    // export path

    M_comm.barrier();


    // Initialise everything that doesn't depend on the mesh (constants, data set description, and time)
    this->initConstant();
    current_time = time_init /*+ pcpt*time_step/(24*3600.0)*/;
    this->initDatasets();

    if (M_rank==0)
    {
        LOG(DEBUG) <<"TIMESTEP= "<< time_step <<"\n";
        LOG(DEBUG) <<"DURATION= "<< duration <<"\n";
    }

    // Initialise the mesh
    this->initMesh();

    // We need to set the scale_coeff et al after initialising the mesh - this was previously done in initConstants
    // The mean resolution of the small_arctic_10km mesh is 7446.71 m. Using 74.5 gives scale_coef = 0.100022, for that mesh
    boost::mpi::broadcast(M_comm, M_res_root_mesh, 0);
    scale_coef = std::sqrt(74.5/M_res_root_mesh/*this->resolution(M_mesh)*/);
    C_fix    = cfix*scale_coef;          // C_fix;...  : cohesion (mohr-coulomb) in MPa (40000 Pa)
    C_alea   = alea_factor*C_fix;        // C_alea;... : alea sur la cohesion (Pa)
    LOG(DEBUG) << "SCALE_COEF = " << scale_coef << "\n";


    if (!M_use_restart)
    {
        // Check the minimum angle of the grid
        minang = this->minAngle(M_mesh);

        if (minang < vm["simul.regrid_angle"].as<double>())
        {
            // LOG(INFO) <<"invalid regridding angle: should be smaller than the minimal angle in the intial grid\n";
            std::cout << "[INFO]: " <<"invalid regridding angle: should be smaller than the minimal angle in the intial grid\n";
            throw std::logic_error("invalid regridding angle: should be smaller than the minimal angle in the intial grid");
        }
    }

    std::cout<<"-----------------------------------------------------------USE RESTART= "<< M_use_restart <<"\n";

    if ( M_use_restart )
    {
        LOG(DEBUG) <<"Reading restart file\n";
        pcpt = this->readRestart(vm["setup.step_nb"].as<int>());

        // Check the minimum angle of the grid
        minang = this->minAngle(M_mesh);

        if (minang < vm["simul.regrid_angle"].as<double>())
        {
            // LOG(INFO) <<"invalid regridding angle: should be smaller than the minimal angle in the intial grid\n";
            std::cout << "[INFO]: " <<"invalid regridding angle: should be smaller than the minimal angle in the intial grid\n";
            throw std::logic_error("invalid regridding angle: should be smaller than the minimal angle in the intial grid");
        }

        // current_time = time_init + pcpt*time_step/(24*3600.0);
        if(M_use_osisaf_drifters)
            this->initOSISAFDrifters();

        if(fmod(pcpt*time_step,output_time_step) == 0)
        {
            LOG(DEBUG) <<"export starts\n";
            this->exportResults((int) pcpt*time_step/output_time_step);
            LOG(DEBUG) <<"export done in " << chrono.elapsed() <<"s\n";
        }

        std::cout<< "-----------------------------------------------------------pcpt= "<< pcpt <<"\n";
    }
    else
    {
        // Do one regrid to get the mesh right
        //this->regrid(pcpt);

        // Initialise variables
        chrono.restart();
        timer["initvar"].first.restart();
        // if (M_rank == 0)
        //     std::cout <<"Initialize variables\n";
        this->initVariables();
        if (M_rank == 0)
            std::cout <<"Initialize variables done "<< timer["initvar"].first.elapsed() <<"s\n";
    }

#if 1

    timer["atmost"].first.restart();
    // if (M_rank == 0)
    //     std::cout <<"Initialize forcingAtmosphere\n";
    this->forcingAtmosphere();
    if (M_rank == 0)
        std::cout <<"Initialize forcingAtmosphere done in "<< timer["atmost"].first.elapsed() <<"s\n";

    timer["ocean"].first.restart();
    // if (M_rank == 0)
    //     std::cout <<"Initialize forcingOcean\n";
    this->forcingOcean();
    if (M_rank == 0)
        std::cout <<"Initialize forcingOcean done in "<< timer["ocean"].first.elapsed() <<"s\n";

    timer["bathy"].first.restart();
    // if (M_rank == 0)
    //     std::cout <<"Initialize bathymetry\n";
    this->bathymetry();
    //if (M_rank == 0)
    //std::cout <<"Initialize bathymetry done in "<< timer["bathy"].first.elapsed() <<"s\n";

    timer["checkload"].first.restart();
    // if (M_rank == 0)
    //     std::cout <<"check_and_reload starts\n";

    this->checkReloadDatasets(M_external_data,
                              time_init,
                              "init - time-dependant");
    if (M_rank == 0)
        std::cout <<"check_and_reload in "<< timer["checkload"].first.elapsed() <<"s\n";

    if (!M_use_restart)
    {
        timer["state"].first.restart();
        // if (M_rank == 0)
        //     std::cout <<"initModelState starts\n";
        this->initModelState();
        if (M_rank == 0)
            std::cout <<"initModelState done in "<< timer["state"].first.elapsed() <<"s\n";
    }

    // Open the output file for drifters
    // TODO: Is this the right place to open the file?
    if ( M_use_iabp_drifters )
    {
        // We should tag the file name with the init time in case of a re-start.
        std::stringstream filename;
        filename << M_export_path << "/drifters_out_" << current_time << ".txt";
        M_iabp_out.open(filename.str(), std::fstream::out);
        if ( ! M_iabp_out.good() )
            throw std::runtime_error("Cannot write to file: " + filename.str());
    }

    // Initialise the moorings - if requested
    if ( M_use_moorings )
        this->initMoorings();

#endif

}

 void
 FiniteElement::step()
 {
     this->updateDrifterPosition();

     // The first time step we behave as if we just did a regrid
     M_regrid = (pcpt==0);

     double displacement_factor = 1.;
     bool force_regrid = (!(pcpt % 30)) && (pcpt != 0);

     if (vm["simul.regrid"].as<std::string>() == "bamg")
     {
         minang = this->minAngle(M_mesh,M_UM,displacement_factor);

         if (M_rank == 0)
         {
             std::cout <<"REGRID ANGLE= "<< minang <<"\n";
         }

         if ( minang < vm["simul.regrid_angle"].as<double>() )
         {
             if ( M_use_moorings && ! M_moorings_snapshot )
                 M_moorings.updateGridMean(M_mesh);

             LOG(DEBUG) <<"Regriding starts\n";
             chrono.restart();
             this->regrid(pcpt);
             LOG(DEBUG) <<"Regriding done in "<< chrono.elapsed() <<"s\n";

             if ( M_use_moorings && ! M_moorings_snapshot )
                 M_moorings.updateGridMean(M_mesh);

             M_regrid = true;
             ++M_nb_regrid;
         }
     }

     M_comm.barrier();

     if ( M_regrid || M_use_restart )
     {
         timer["tensors"].first.restart();
         this->tensors();
         if (M_rank == 0)
             std::cout <<"---timer tensors:              "<< timer["tensors"].first.elapsed() <<"\n";

         timer["cohesion"].first.restart();
         this->cohesion();
         if (M_rank == 0)
             std::cout <<"---timer cohesion:             "<< timer["cohesion"].first.elapsed() <<"\n";

         timer["coriolis"].first.restart();
         this->coriolis();
         if (M_rank == 0)
             std::cout <<"---timer coriolis:             "<< timer["coriolis"].first.elapsed() <<"\n";
     }

     timer["reload"].first.restart();
     this->checkReloadDatasets(M_external_data,
                               current_time+time_step/(24*3600.0),
                               "init - time-dependant");
     if (M_rank == 0)
         std::cout <<"---timer check_and_reload:     "<< timer["reload"].first.elapsed() <<"s\n";

     if (pcpt == 0)
     {
         LOG(DEBUG) <<"first export starts\n";
         this->exportResults(0);
         this->writeRestart(pcpt, 0); // Write a restart before regrid - useful for debugging
         LOG(DEBUG) <<"first export done\n";
     }

    //======================================================================
    // Do the thermodynamics
    //======================================================================
    if(vm["simul.use_thermo_forcing"].as<bool>())
    {
        timer["thermo"].first.restart();
        this->thermo();
        if (M_rank == 0)
            std::cout <<"---timer thermo:               "<< timer["thermo"].first.elapsed() <<"s\n";
    }

    //======================================================================
    // Do the dynamics
    //======================================================================

    if ( M_dynamics_type == setup::DynamicsType::DEFAULT )
    {
        //======================================================================
        // Assemble the matrix
        //======================================================================
        timer["assemble"].first.restart();
        this->assemble(pcpt);
        if (M_rank == 0)
            std::cout <<"---timer assemble:             "<< timer["assemble"].first.elapsed() <<"s\n";

#if 0
        if(had_remeshed && (vm["simul.regrid_output_flag"].as<bool>()))
        {
            std::string tmp_string3    = (boost::format( "after_assemble_%1%_mesh_adapt_step_%2%" )
                                                                             % pcpt
                                          % mesh_adapt_step ).str();

            this->exportResults(tmp_string3);

            had_remeshed=false;
        }
#endif

        //======================================================================
        // Solve the linear problem
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

    ++pcpt;
    current_time = time_init + pcpt*time_step/(24*3600.0);

#if 1
#if 1
    if(fmod(pcpt*time_step,output_time_step) == 0)
    {
        chrono.restart();
        LOG(DEBUG) <<"export starts\n";
        this->exportResults((int) pcpt*time_step/output_time_step);
        LOG(DEBUG) <<"export done in " << chrono.elapsed() <<"s\n";
    }

    if ( M_use_moorings )
    {
        this->updateMoorings();
    }

    if ( fmod(pcpt*time_step,restart_time_step) == 0)
    {
        std::cout << "Writing restart file after time step " <<  pcpt-1 << "\n";
        this->writeRestart(pcpt, (int) pcpt*time_step/restart_time_step);
    }
#endif
#endif

 }

 // This is the main working function, called from main.cpp (same as perform_simul in the old code)
void
FiniteElement::run()
{
    if (M_comm.rank() == 0)
    {
        // LOG(INFO) << "-----------------------Simulation started on "<< current_time_local() <<"\n";
        std::cout << "[INFO]: " << "-----------------------Simulation started on "<< current_time_local() <<"\n";
    }

    std::string current_time_system = current_time_local();
    double displacement_factor = 1.;

    pcpt = 0;
    niter = vm["simul.maxiteration"].as<int>();
    mesh_adapt_step = 0;
    had_remeshed = false;
    minang = 0.;

    this->init();

    // write the logfile: assigned to the process master (rank 0)
    if (M_comm.rank() == 0)
    {
        this->writeLogFile();
    }

    // Debug file that records the time step
    std::fstream pcpt_file;

    if (M_rank == 0)
    {
        pcpt_file.open(M_export_path + "/Timestamp.txt", std::ios::out | std::ios::trunc);
    }

#if 1
    current_time = time_init + pcpt*time_step/(24*3600.0);
    bool is_running = true;

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
                std::cout <<" ---------- progression: ("<< 100.0*(pcpt*time_step/duration) <<"%)"
                          <<" ---------- time spent: "<< time_spent(current_time_system);
            }

            std::cout <<"\n";
        }

        std::cout<<"pcpt = "<< pcpt <<"\n";

        //is_running = ((pcpt+1)*time_step) < duration;
        is_running = (pcpt*time_step) < duration;

        if (pcpt == niter)
            is_running = false;

        //current_time = time_init + pcpt*time_step/days_in_sec;

        // **********************************************************************
        // Take one time-step
        // **********************************************************************
        this->step();

        if (M_rank == 0)
        {
            pcpt_file << pcpt << "\n";
            pcpt_file << to_date_string(current_time) << "\n";
            pcpt_file.seekp(0);
        }
    }

    if (M_rank == 0)
    {
        pcpt_file.close();
    }

    this->exportResults(1000);

    //this->clear();
    this->finalise();

#endif

    M_comm.barrier();

    if (M_rank==0)
    {
        std::cout <<"nb regrid total = " << M_nb_regrid <<"\n";

        std::cout << "[INFO]: " << "-----------------------Simulation done on "<< current_time_local() <<"\n";
        std::cout << "[INFO]: " << "-----------------------Total time spent:  "<< time_spent(current_time_system) <<"\n";
    }
}

void
FiniteElement::updateDrifterPosition()
{
    // Update the drifters position twice a day, important to keep the same frequency as the IABP data, for the moment
    if( pcpt==0 || std::fmod(current_time,0.5)==0 )
    {
        // Read in the new buoys and output
        if ( M_use_iabp_drifters )
        {
            this->updateIABPDrifter();

            chrono.restart();
            LOG(DEBUG) <<"Drifter starts\n";
            LOG(DEBUG) <<"DRIFTER: Interp starts\n";

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

            // Interpolate the total displacement and concentration onto the drifter positions
            int nb_var=2;
            std::vector<double> interp_drifter_in(nb_var*M_mesh.numNodes());

            // Interpolate the velocity
            for (int i=0; i<M_mesh.numNodes(); ++i)
            {
                interp_drifter_in[nb_var*i]   = M_UT[i];
                interp_drifter_in[nb_var*i+1] = M_UT[i+M_mesh.numNodes()];
            }

            double* interp_drifter_out;
            InterpFromMeshToMesh2dx(&interp_drifter_out,
                                    &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                                    M_mesh.numNodes(),M_mesh.numTriangles(),
                                    &interp_drifter_in[0],
                                    M_mesh.numNodes(),nb_var,
                                    &drifter_X[0],&drifter_Y[0],M_iabp_drifters.size(),
                                    true, 0.);


            // Interpolate the concentration - re-use interp_drifter_in
            interp_drifter_in.resize(M_mesh.numTriangles());
            for (int i=0; i<M_mesh.numTriangles(); ++i)
            {
                interp_drifter_in[i]   = M_conc[i];
            }

            double* interp_drifter_c_out;
            InterpFromMeshToMesh2dx(&interp_drifter_c_out,
                                    &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                                    M_mesh.numNodes(),M_mesh.numTriangles(),
                                    &interp_drifter_in[0],
                                    M_mesh.numTriangles(),1,
                                    &drifter_X[0],&drifter_Y[0],M_iabp_drifters.size(),
                                    true, 0.);

            // Rebuild the M_iabp_drifters map
            double clim = vm["simul.drifter_climit"].as<double>();
            j=0;
            for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); /* ++it is not allowed here, because we use 'erase' */ )
            {
                if ( interp_drifter_c_out[j] > clim )
                {
                    M_iabp_drifters[it->first] = std::array<double,2> {it->second[0]+interp_drifter_out[nb_var*j], it->second[1]+interp_drifter_out[nb_var*j+1]};
                    ++it;
                } else {
                    // Throw out drifters that drift out of the ice
                    it = M_iabp_drifters.erase(it);
                }
                ++j;
            }

            xDelete<double>(interp_drifter_out);
            xDelete<double>(interp_drifter_c_out);

            LOG(DEBUG) <<"DRIFTER: Interp done\n";
            LOG(DEBUG) <<"Drifter interp done in "<< chrono.elapsed() <<"s\n";

            // TODO: Do we want to output drifters at a different time interval?
            this->outputDrifter(M_iabp_out);
        }

        if ( M_use_equallyspaced_drifters )
            M_equallyspaced_drifters.move(M_mesh, M_UT);
        if ( M_use_rgps_drifters )
            M_rgps_drifters.move(M_mesh, M_UT);
        if ( M_use_osisaf_drifters )
            for (auto it=M_osisaf_drifters.begin(); it!=M_osisaf_drifters.end(); it++)
                it->move(M_mesh, M_UT);

	    for (int i=0; i<M_num_nodes; ++i)
	    {
	        // UM
		    M_UT[i] = 0.;
		    M_UT[i+M_num_nodes] = 0.;
	    }
    }

    if(pcpt>0)
    {
        if ( M_use_equallyspaced_drifters && fmod(current_time,M_equallyspaced_drifters_output_time_step) == 0 )
            M_equallyspaced_drifters.appendNetCDF(current_time, M_mesh, M_UT);

        if ( M_use_rgps_drifters )
        {
            std::string time_str = vm["simul.RGPS_time_init"].as<std::string>();
            double RGPS_time_init = from_date_time_string(time_str);

            if( !M_rgps_drifters.isInitialised() && current_time == RGPS_time_init)
                this->updateRGPSDrifters();

            if( current_time != RGPS_time_init && fmod(current_time,M_rgps_drifters_output_time_step) == 0 )
                if ( M_rgps_drifters.isInitialised() )
                    M_rgps_drifters.appendNetCDF(current_time, M_mesh, M_UT);
        }
    }

    if ( M_use_osisaf_drifters && fmod(current_time+0.5,1.) == 0 )
    {
        // OSISAF drift is calculated as a dirfter displacement over 48 hours
        // and they have two sets of drifters in the field at all times.

        // Write out the contents of [1] if it's meaningfull
        if ( M_osisaf_drifters[1].isInitialised() )
            M_osisaf_drifters[1].appendNetCDF(current_time, M_mesh, M_UT);

        // Flip the vector so we move [0] to be [1]
        std::reverse(M_osisaf_drifters.begin(), M_osisaf_drifters.end());

        // Create a new M_drifters instance in [0], with a properly initialised netCDF file
        M_osisaf_drifters[0] = Drifters("data", "ice_drift_nh_polstere-625_multi-oi.nc", "yc", "yx", "lat", "lon", M_mesh, M_conc, vm["simul.drifter_climit"].as<double>());
        M_osisaf_drifters[0].initNetCDF(M_export_path+"/OSISAF_", current_time);
        M_osisaf_drifters[0].appendNetCDF(current_time, M_mesh, M_UT);
    }
}

// Add to the mean on the mesh
void
FiniteElement::updateMeans(GridOutput &means, double time_factor)
{
    // Update elements and multiply with time_factor
    for ( auto it=means.M_elemental_variables.begin(); it!=means.M_elemental_variables.end(); ++it )
    {
        switch (it->varID)
        {
            // Prognostic variables
            case (GridOutput::variableID::conc):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_conc[i]*time_factor;
                break;

            case (GridOutput::variableID::thick):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_thick[i]*time_factor;
                break;

            case (GridOutput::variableID::damage):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_damage[i]*time_factor;
                break;

            case (GridOutput::variableID::snow):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_snow_thick[i]*time_factor;
                break;

            case (GridOutput::variableID::tsurf):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += ( M_conc[i]*M_tice[0][i] + M_conc_thin[i]*M_tsurf_thin[i] + (1-M_conc[i]-M_conc_thin[i])*M_sst[i] )*time_factor;
                break;

            case (GridOutput::variableID::sst):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_sst[i]*time_factor;
                break;

            case (GridOutput::variableID::sss):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_sss[i]*time_factor;
                break;

            case (GridOutput::variableID::tsurf_ice):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_tice[0][i]*time_factor;
                break;

            case (GridOutput::variableID::t1):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_tice[1][i]*time_factor;
                break;

            case (GridOutput::variableID::t2):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_tice[2][i]*time_factor;
                break;

            case (GridOutput::variableID::conc_thin):
                for (int i=0; i<M_num_elements; i++)
                {
                    it->data_mesh[i] += M_conc_thin[i]*time_factor;
                }
                break;

            case (GridOutput::variableID::h_thin):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_h_thin[i]*time_factor;
                break;

            case (GridOutput::variableID::hs_thin):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_hs_thin[i]*time_factor;
                break;

            // Diagnostic variables
            case (GridOutput::variableID::Qa):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += D_Qa[i]*time_factor;
                break;
            case (GridOutput::variableID::Qsw):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += D_Qsw[i]*time_factor;
                break;
            case (GridOutput::variableID::Qlw):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += D_Qlw[i]*time_factor;
                break;
            case (GridOutput::variableID::Qsh):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += D_Qsh[i]*time_factor;
                break;
            case (GridOutput::variableID::Qlh):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += D_Qlh[i]*time_factor;
                break;
            case (GridOutput::variableID::Qo):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += D_Qo[i]*time_factor;
                break;
            case (GridOutput::variableID::delS):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += D_delS[i]*time_factor;
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
                    it->data_mesh[i] += M_VT[i];
                break;

            case (GridOutput::variableID::VT_y):
                for (int i=0; i<M_num_nodes; i++)
                    it->data_mesh[i] += M_VT[i+M_num_nodes];
                break;

            default: std::logic_error("Updating of given variableID not implimented (nodes)");
        }
    }
}

// Initialise everything w.r.t. the moorings
void
FiniteElement::initMoorings()
{
    // Output and averaging grids
    std::vector<double> data_nodes(M_num_nodes);
    std::vector<double> data_elements(M_num_elements);
    std::vector<double> data_grid;

    // Output variables - elements
    GridOutput::Variable conc(GridOutput::variableID::conc, data_elements, data_grid);
    GridOutput::Variable thick(GridOutput::variableID::thick, data_elements, data_grid);
    GridOutput::Variable snow(GridOutput::variableID::snow, data_elements, data_grid);
    GridOutput::Variable tsurf(GridOutput::variableID::tsurf, data_elements, data_grid);
    GridOutput::Variable Qa(GridOutput::variableID::Qa, data_elements, data_grid);
    GridOutput::Variable Qsw(GridOutput::variableID::Qsw, data_elements, data_grid);
    GridOutput::Variable Qlw(GridOutput::variableID::Qlw, data_elements, data_grid);
    GridOutput::Variable Qsh(GridOutput::variableID::Qsh, data_elements, data_grid);
    GridOutput::Variable Qlh(GridOutput::variableID::Qlh, data_elements, data_grid);
    GridOutput::Variable Qo(GridOutput::variableID::Qo, data_elements, data_grid);
    GridOutput::Variable delS(GridOutput::variableID::delS, data_elements, data_grid);

    std::vector<GridOutput::Variable> elemental_variables(11);
    elemental_variables[0] = conc;
    elemental_variables[1] = thick;
    elemental_variables[2] = snow;
    elemental_variables[3] = tsurf;
    elemental_variables[4] = Qa;
    elemental_variables[5] = Qsw;
    elemental_variables[6] = Qlw;
    elemental_variables[7] = Qsh;
    elemental_variables[8] = Qlh;
    elemental_variables[9] = Qo;
    elemental_variables[10] = delS;
    if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
    {
        GridOutput::Variable conc_thin(GridOutput::variableID::conc_thin, data_elements, data_grid);
        GridOutput::Variable h_thin(GridOutput::variableID::h_thin, data_elements, data_grid);
        GridOutput::Variable hs_thin(GridOutput::variableID::hs_thin, data_elements, data_grid);

        elemental_variables.push_back(conc_thin);
        elemental_variables.push_back(h_thin);
        elemental_variables.push_back(hs_thin);
    }

    // Output variables - nodes
    GridOutput::Variable siu(GridOutput::variableID::VT_x, data_nodes, data_grid);

    GridOutput::Variable siv(GridOutput::variableID::VT_y, data_nodes, data_grid);

    std::vector<GridOutput::Variable> nodal_variables(2);
    nodal_variables[0] = siu;
    nodal_variables[1] = siv;

    // The vectorial variables are (always on the nodes) ...
    std::vector<int> siuv_id(2);
    siuv_id[0] = 0;
    siuv_id[1] = 1;

    GridOutput::Vectorial_Variable siuv{
        components_Id: siuv_id,
        // east_west_oriented: true
        east_west_oriented: false
    };

    std::vector<GridOutput::Vectorial_Variable> vectorial_variables(1);
    vectorial_variables[0] = siuv;

#if 1
    // Calculate the grid spacing (assuming a regular grid for now)
    auto RX = M_mesh.coordX();
    auto RY = M_mesh.coordY();
    auto xcoords = std::minmax_element( RX.begin(), RX.end() );
    auto ycoords = std::minmax_element( RY.begin(), RY.end() );

    double mooring_spacing = 1e3 * vm["simul.mooring_spacing"].as<double>();
    int ncols = (int) ( 0.5 + ( *xcoords.second - *xcoords.first )/mooring_spacing );
    int nrows = (int) ( 0.5 + ( *ycoords.second - *ycoords.first )/mooring_spacing );

    // Define the mooring dataset
    M_moorings = GridOutput(ncols, nrows, mooring_spacing, *xcoords.first, *ycoords.first, nodal_variables, elemental_variables, vectorial_variables);
#else
    // Read the grid in from file
    // Define a grid
    GridOutput::Grid grid{
        gridFile: "ice_drift_nh_polstere-625_multi-oi.nc",
        dirname: "data",
        mpp_file: Environment::vm()["mesh.mppfile"].as<std::string>(),
        dimNameX: "yc",
        dimNameY: "xc",
        latName: "lat",
        lonName: "lon"
    };

    // Define the mooring dataset
    M_moorings = GridOutput(grid, nodal_variables, elemental_variables, vectorial_variables);

    /* Just for debuging
    // Save the grid info - this is still just an ascii dump!
    std::ofstream myfile;
    myfile.open("lon_grid.dat");
    std::copy(M_moorings.M_grid.gridLON.begin(), M_moorings.M_grid.gridLON.end(), ostream_iterator<float>(myfile," "));
    myfile.close();
    myfile.open("lat_grid.dat");
    std::copy(M_moorings.M_grid.gridLAT.begin(), M_moorings.M_grid.gridLAT.end(), ostream_iterator<float>(myfile," "));
    myfile.close();

    myfile.open("x_grid.dat");
    std::copy(M_moorings.M_grid.gridX.begin(), M_moorings.M_grid.gridX.end(), ostream_iterator<float>(myfile," "));
    myfile.close();
    myfile.open("y_grid.dat");
    std::copy(M_moorings.M_grid.gridY.begin(), M_moorings.M_grid.gridY.end(), ostream_iterator<float>(myfile," "));
    myfile.close();
    */
#endif

    double output_time;
    if ( M_moorings_snapshot )
        // shift the timestamp in the file to the centre of the output interval
        output_time = current_time;
    else
        output_time = current_time - mooring_output_time_step/86400/2;

    M_moorings_file = M_moorings.initNetCDF(M_export_path + "/Moorings", M_moorings_file_length, output_time);
}

void
FiniteElement::updateMoorings()
{
    // If we're taking snapshots the we only call updateMeans before writing to file
    if ( ! M_moorings_snapshot )
        this->updateMeans(M_moorings, mooring_time_factor);

    if ( fmod(pcpt*time_step,mooring_output_time_step) == 0 )
    {
        double output_time;
        if ( M_moorings_snapshot )
        {
            // Update the snapshot
            this->updateMeans(M_moorings, 1.);
            // shift the timestamp in the file to the centre of the output interval
            output_time = current_time;
        }
        else
        {
            output_time = current_time - mooring_output_time_step/86400/2;
        }

        // If it's a new day we check if we need a new file
        double not_used;
        if ( M_moorings_file_length != GridOutput::fileLength::inf && modf(output_time, &not_used) < time_step*86400 )
        {
            boost::gregorian::date now = Nextsim::parse_date(output_time);
            switch (M_moorings_file_length)
            {
            case GridOutput::fileLength::daily:
                M_moorings_file = M_moorings.initNetCDF(M_export_path + "/Moorings", M_moorings_file_length, output_time);
                break;
            case GridOutput::fileLength::weekly:
                if ( now.day_of_week().as_number() == 1 )
                    M_moorings_file = M_moorings.initNetCDF(M_export_path + "/Moorings", M_moorings_file_length, output_time);
                break;
            case GridOutput::fileLength::monthly:
                if ( now.day().as_number() == 1 )
                    M_moorings_file = M_moorings.initNetCDF(M_export_path + "/Moorings", M_moorings_file_length, output_time);
                break;
            case GridOutput::fileLength::yearly:
                if ( now.day_of_year() == 1 )
                    M_moorings_file = M_moorings.initNetCDF(M_export_path + "/Moorings", M_moorings_file_length, output_time);
            }
        }

        M_moorings.updateGridMean(M_mesh);

        M_moorings.appendNetCDF(M_moorings_file, output_time);

        M_moorings.resetMeshMean(M_mesh);
        M_moorings.resetGridMean();
    }
}

void
FiniteElement::writeRestart(int pcpt, int step)
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

    std::cout<< "Rank "<< M_rank <<": Arrived here\n";

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

    // fields defined on mesh elements
    std::vector<int> sizes_elements = M_sizes_elements;
    int nb_var_element = 15 + M_tice.size();//15;
    std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });
    std::vector<double> interp_elt_in_local(nb_var_element*M_local_nelements);


    for (int i=0; i<M_local_nelements; ++i)
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
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_ridge_ratio[i];
        tmp_nb_var++;

        // random_number
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_random_number[i];
        tmp_nb_var++;

        // sss
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sss[i];
        tmp_nb_var++;

        // sst
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_sst[i];
        tmp_nb_var++;

        // Ice temperature
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_tice[0][i];
        tmp_nb_var++;

        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_tice[1][i];
            tmp_nb_var++;

            interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_tice[2][i];
            tmp_nb_var++;
        }

        // thin ice thickness
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_h_thin[i];
        tmp_nb_var++;

        // thin ice thickness
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_conc_thin[i];
        tmp_nb_var++;

        // snow on thin ice
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_hs_thin[i];
        tmp_nb_var++;

        // Ice surface temperature for thin ice
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_tsurf_thin[i];
        tmp_nb_var++;

        if(tmp_nb_var>nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }
    }

    std::vector<double> interp_in_elements;

    if (M_rank == 0)
    {
        interp_in_elements.resize(nb_var_element*M_mesh_root.numTriangles());
        boost::mpi::gatherv(M_comm, interp_elt_in_local, &interp_in_elements[0], sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, interp_elt_in_local, 0);
    }

    M_comm.barrier();

    if (M_rank == 0)
    {
        int num_elements_root = M_mesh_root.numTriangles();
        int tice_size = M_tice.size();

        std::cout<<"tice_size        = "<< tice_size <<"\n";

        std::cout<< "["<< M_rank << "] "<<"num_elements_root                   = "<< num_elements_root <<"\n";

        std::vector<double> M_conc_root(num_elements_root);
        std::vector<double> M_thick_root(num_elements_root);
        std::vector<double> M_snow_thick_root(num_elements_root);
        std::vector<double> M_sigma_root(3*num_elements_root);
        std::vector<double> M_damage_root(num_elements_root);
        std::vector<double> M_ridge_ratio_root(num_elements_root);
        std::vector<double> M_random_number_root(num_elements_root);
        std::vector<double> M_sss_root(num_elements_root);
        std::vector<double> M_sst_root(num_elements_root);
        std::vector<double> M_tice_root(tice_size*num_elements_root);
        std::vector<double> M_h_thin_root(num_elements_root);
        std::vector<double> M_conc_thin_root(num_elements_root);
        std::vector<double> M_hs_thin_root(num_elements_root);
        std::vector<double> M_tsurf_thin_root(num_elements_root);

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
            M_tice_root[tice_size*i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            if ( M_thermo_type == setup::ThermoType::WINTON )
            {
                // tice2
                M_tice_root[tice_size*i+1] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;

                // tice3
                M_tice_root[tice_size*i+2] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
                tmp_nb_var++;
            }

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

            if(tmp_nb_var>nb_var_element)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }


        Exporter exporter;
        std::string filename;

        // === Start with the mesh ===
        // First the data
        std::string directory = Environment::nextsimDir().string() + "/restart";
        // change directory for outputs if the option "output_directory" is not empty
        if ( ! (vm["simul.output_directory"].as<std::string>()).empty() )
            directory = vm["simul.output_directory"].as<std::string>() + "/restart";

        // create the output directory if it does not exist
        fs::path path(directory);
        if ( !fs::exists(path) )
            fs::create_directories(path);

        filename = (boost::format( "%1%/mesh_%2%.bin" )
                    % directory
                    % step ).str();

        std::fstream meshbin(filename, std::ios::binary | std::ios::out | std::ios::trunc);
        if ( ! meshbin.good() )
            throw std::runtime_error("Cannot write to file: " + filename);
        exporter.writeMesh(meshbin, M_mesh_root);
        meshbin.close();

        // Then the record
        filename = (boost::format( "%1%/mesh_%2%.dat" )
                    % directory
                    % step ).str();

        std::fstream meshrecord(filename, std::ios::out | std::ios::trunc);
        if ( ! meshrecord.good() )
            throw std::runtime_error("Cannot write to file: " + filename);
        exporter.writeRecord(meshrecord,"mesh");
        meshrecord.close();

        // === Write the prognostic variables ===
        // First the data
        filename = (boost::format( "%1%/field_%2%.bin" )
                    % directory
                    % step ).str();
        std::fstream outbin(filename, std::ios::binary | std::ios::out | std::ios::trunc );
        if ( ! outbin.good() )
            throw std::runtime_error("Cannot write to file: " + filename);

        std::vector<int> misc_int(3);

        std::cout<<"*********pcpt           = "<< pcpt <<"\n";
        std::cout<<"*********M_flag_fix     = "<< M_flag_fix <<"\n";
        std::cout<<"*********current_time   = "<< current_time <<"\n";
        std::cout<<"*********mesh_adapt_step= "<< mesh_adapt_step <<"\n";

        misc_int[0] = pcpt;
        misc_int[1] = M_flag_fix;
        //misc_int[2] = current_time;
        misc_int[2] = mesh_adapt_step;

        exporter.writeField(outbin, misc_int, "Misc_int");
        exporter.writeField(outbin, M_dirichlet_flags_root, "M_dirichlet_flags");

        exporter.writeField(outbin, M_conc_root, "M_conc");
        exporter.writeField(outbin, M_thick_root, "M_thick");
        exporter.writeField(outbin, M_snow_thick_root, "M_snow_thick");
        exporter.writeField(outbin, M_sigma_root, "M_sigma");
        exporter.writeField(outbin, M_damage_root, "M_damage");
        exporter.writeField(outbin, M_ridge_ratio_root, "M_ridge_ratio");
        exporter.writeField(outbin, M_random_number_root, "M_random_number");

        exporter.writeField(outbin, M_tice_root, "M_Tice");

        // int i=0;
        // for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
        // {
        //     std::cout<<"TICE TO STORE SIE= "<< it->size() <<"\n";
        //     exporter.writeField(outbin, *it, "M_Tice_" + std::to_string(i));
        //     i++;
        // }

        exporter.writeField(outbin, M_sst_root, "M_sst");
        exporter.writeField(outbin, M_sss_root, "M_sss");
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

        if (M_use_iabp_drifters)
        {
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

        outbin.close();

        // Then the record
        filename = (boost::format( "%1%/field_%2%.dat" )
                    % directory
                    % step ).str();

        std::fstream outrecord(filename, std::ios::out | std::ios::trunc);
        if ( ! outrecord.good() )
        {
            throw std::runtime_error("Cannot write to file: " + filename);
        }

        exporter.writeRecord(outrecord);
        outrecord.close();
    }
}

int
FiniteElement::readRestart(int step)
{
    Exporter exp_field, exp_mesh;
    std::string filename;
    boost::unordered_map<std::string, std::vector<int>>    field_map_int;
    boost::unordered_map<std::string, std::vector<double>> field_map_dbl;

    //std::fstream inbin;

    if (M_rank == 0)
    {
        M_mesh_filename = (boost::format( "par%1%%2%" ) % M_comm.size() % M_mesh_filename ).str();

        // === Read in the mesh restart files ===
        std::string restart_path;
        if ( (vm["setup.restart_path"].as<std::string>()).empty() )
        {
            //default restart path is $NEXTSIMDIR/restart
            restart_path = Environment::nextsimDir().string()+"/restart";
        }
        else
        {
            restart_path = vm["setup.restart_path"].as<std::string>();
        }

        // Start with the record
        filename = (boost::format( "%1%/mesh_%2%.dat" )
                    % restart_path
                    % step ).str();

        std::ifstream meshrecord(filename);
        if ( ! meshrecord.good() )
        {
            throw std::runtime_error("File not found: " + filename);
        }

        exp_mesh.readRecord(meshrecord);
        meshrecord.close();

        // Then onto the data itself
        filename = (boost::format( "%1%/mesh_%2%.bin" )
                    % restart_path
                    % step ).str();

        std::fstream meshbin(filename, std::ios::binary | std::ios::in );
        if ( ! meshbin.good() )
        {
            throw std::runtime_error("File not found: " + filename);
        }

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
                    % step ).str();

        std::ifstream inrecord(filename);
        if ( ! inrecord.good() )
        {
            throw std::runtime_error("File not found: " + filename);
        }

        std::cout<<"filename= "<< filename <<"\n";

        exp_field.readRecord(inrecord);
        inrecord.close();

        // Then onto the data itself
        filename = (boost::format( "%1%/field_%2%.bin" )
                    % restart_path
                    % step ).str();

        //inbin.open(filename, std::ios::binary | std::ios::in);
        std::fstream inbin(filename, std::ios::binary | std::ios::in);

        if ( ! inbin.good() )
        {
            throw std::runtime_error("File not found: " + filename);
        }

        std::cout<<"filename= "<< filename <<"\n";

        field_map_int.clear();
        field_map_dbl.clear();
        exp_field.loadFile(inbin, field_map_int, field_map_dbl);
        inbin.close();

        std::cout<<".......................END PROCESS\n";

        // === Recreate the mesh ===
        // Create bamgmesh and bamggeom
        BamgConvertMeshx(
                         bamgmesh_root,bamggeom_root,
                         &indexTr[0],&coordX[0],&coordY[0],
                         coordX.size(), indexTr.size()/3.
                         );

        std::cout<<".....................SIZE= "<< field_map_int["Misc_int"].size() <<"\n";


        // Fix boundaries
        auto map_field  = field_map_int["Misc_int"];
        int pcpt        = map_field.at(0);
        M_flag_fix      = map_field.at(1);
        mesh_adapt_step = map_field.at(2);

        std::cout<<".....................PCPT      = "<< pcpt <<"\n";
        std::cout<<".....................M_flag_fix= "<< M_flag_fix <<"\n";
        //std::cout<<".....................M_current_time= "<< current_time <<"\n";
        std::cout<<".....................M_mesh_adapt_step= "<< mesh_adapt_step <<"\n";

        current_time = time_init + (pcpt*time_step)/(24*3600);

        std::cout<<"current_time= "<< current_time <<"\n";


        //int pcpt = (current_time-time_init)*(24*3600)/time_step;

        std::vector<int> dirichlet_flags = field_map_int["M_dirichlet_flags"];

        for (int edg=0; edg<bamgmesh_root->EdgesSize[0]; ++edg)
        {
            int fnd = bamgmesh_root->Edges[3*edg];

            if ((std::binary_search(dirichlet_flags .begin(),dirichlet_flags.end(),fnd)))
            {
                bamggeom_root->Edges[3*edg+2] = M_flag_fix;
                bamgmesh_root->Edges[3*edg+2] = M_flag_fix;
            }
            // else
            // {
            //     bamggeom_root->Edges[3*edg+2] = M_flag_fix+1; // we just want it to be different than M_flag_fix
            //     bamgmesh_root->Edges[3*edg+2] = M_flag_fix+1; // we just want it to be different than M_flag_fix
            // }
        }

        // Import the bamg structs
        this->importBamg(bamgmesh_root);

        // update the boundary flags
        this->updateBoundaryFlags();
    }

    // mesh partitioning
    this->partitionMeshRestart();

    // Initialise all the variables to zero
    // this->initVariables();
#if 1

    if (M_rank == 0)
    {
        int num_elements_root = M_mesh_root.numTriangles();
        int tice_size = M_tice.size();

        M_VT.resize(2*M_ndof);
        M_VTM.resize(2*M_ndof);
        M_VTMM.resize(2*M_ndof);
        M_UM.resize(2*M_ndof);
        M_UT.resize(2*M_ndof);

        M_conc.resize(num_elements_root);
        M_thick.resize(num_elements_root);
        M_snow_thick.resize(num_elements_root);
        M_sigma.resize(3*num_elements_root);
        M_damage.resize(num_elements_root);
        M_ridge_ratio.resize(num_elements_root);
        M_random_number.resize(num_elements_root);
        M_sss.resize(num_elements_root);
        M_sst.resize(num_elements_root);


        std::vector<double> M_tice_all(tice_size*num_elements_root);

        switch (M_thermo_type)
        {
            case (setup::ThermoType::ZERO_LAYER):
                M_tice.resize(1);
                break;
            case (setup::ThermoType::WINTON):
                M_tice.resize(3);
                break;
            default:
                std::cout << "thermo_type= " << (int)M_thermo_type << "\n";
                throw std::logic_error("Wrong thermo_type");
        }

        for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
        {
            it->resize(num_elements_root);
        }

        for (int i=0; i<M_tice.size(); ++i)
        {
            std::cout<<"M_tice[i]= "<< (M_tice[i]).size() <<"\n";
        }

        M_h_thin.resize(num_elements_root);
        M_conc_thin.resize(num_elements_root);
        M_hs_thin.resize(num_elements_root);
        M_tsurf_thin.resize(num_elements_root);

        M_conc            = field_map_dbl["M_conc"];
        M_thick           = field_map_dbl["M_thick"];
        M_snow_thick      = field_map_dbl["M_snow_thick"];
        M_sigma           = field_map_dbl["M_sigma"];
        M_damage          = field_map_dbl["M_damage"];
        M_ridge_ratio     = field_map_dbl["M_ridge_ratio"];
        M_random_number   = field_map_dbl["M_random_number"];
        M_tice_all   = field_map_dbl["M_Tice"];

        for (int i=0; i<num_elements_root; ++i)
        {
            M_tice[0][i] = M_tice_all[tice_size*i];

            if ( M_thermo_type == setup::ThermoType::WINTON )
            {
                M_tice[1][i] = M_tice_all[tice_size*i+1];
                M_tice[2][i] = M_tice_all[tice_size*i+2];
            }
        }

        // int i=0;
        // for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
        // {
        //     *it = field_map_dbl["M_Tice_"+std::to_string(i)];
        //     i++;
        // }

        M_sst        = field_map_dbl["M_sst"];
        M_sss        = field_map_dbl["M_sss"];

        // Pre-processing
        if(vm["setup.restart_at_rest"].as<bool>())
        {
            for (int i=0; i < M_sigma.size(); i++)
            {
                M_sigma[i] = 0.;
            }

            for (int i=0; i < M_VT.size(); i++)
            {
                M_VT[i]   = 0.;
                M_VTM[i]  = 0.;
                M_VTMM[i] = 0.;
                M_UM[i]   = 0.;
                M_UT[i]   = 0.;
            }
        }
        else
        {
            M_VT         = field_map_dbl["M_VT"];
            M_VTM        = field_map_dbl["M_VTM"];
            M_VTMM       = field_map_dbl["M_VTMM"];
            M_UM         = field_map_dbl["M_UM"];
            M_UT         = field_map_dbl["M_UT"];
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

            this->initIABPDrifter();
            if (drifter_no.size() == 0)
            {
                LOG(WARNING) << "Warning: Couldn't read drifter positions from restart file. Drifter positions initialised as if there was no restart.\n";
                this->updateIABPDrifter();
            }
            else
            {
                for ( int i=0; i<drifter_no.size(); ++i )
                {
                    M_iabp_drifters.emplace(drifter_no[i], std::array<double,2>{drifter_x[i], drifter_y[i]});
                }
            }
        }
#if 0
#endif
    }

    std::vector<double> interp_elt_out;
    std::vector<double> interp_nd_out;
    this->collectRootRestart(interp_elt_out, interp_nd_out);

    // Initialise all the variables to zero
    this->initVariables();

    this->scatterFieldsElement(&interp_elt_out[0]);

    this->scatterFieldsNode(&interp_nd_out[0]);

#endif
    return pcpt;
}

void
FiniteElement::partitionMeshRestart()
{
    M_comm.barrier();

    if (M_rank == 0)
    {
        std::cout<<"------------------------------version       = "<< M_mesh_root.version() <<"\n";
        std::cout<<"------------------------------ordering      = "<< M_mesh_root.ordering() <<"\n";
        std::cout<<"------------------------------format        = "<< M_mesh_fileformat <<"\n";
        std::cout<<"------------------------------space         = "<< vm["mesh.partition-space"].as<std::string>() <<"\n";
        std::cout<<"------------------------------partitioner   = "<< vm["mesh.partitioner"].as<std::string>() <<"\n";

        // Environment::logMemoryUsage("before partitioning...");
        timer["savemesh"].first.restart();
        LOG(DEBUG) <<"Saving mesh starts\n";
        if (M_partition_space == mesh::PartitionSpace::MEMORY)
        {
            M_mesh_root.writeToGModel(M_mesh_filename);
        }
        else if (M_partition_space == mesh::PartitionSpace::DISK)
        {
            M_mesh_root.writeTofile(M_mesh_filename);
        }

        std::cout <<"Saving mesh done in "<< timer["savemesh"].first.elapsed() <<"s\n";

        // partition the mesh on root process (rank 0)
        timer["meshpartition"].first.restart();
        LOG(DEBUG) <<"Partitioning mesh starts\n";
        M_mesh_root.partition(M_mesh_filename,M_partitioner,M_partition_space, M_mesh_fileformat);
        std::cout <<"Partitioning mesh done in "<< timer["meshpartition"].first.elapsed() <<"s\n";
    }

    M_prv_local_ndof = M_local_ndof;
    M_prv_num_nodes = M_num_nodes;
    M_prv_num_elements = M_local_nelements;
    M_prv_global_num_nodes = M_mesh.numGlobalNodes();
    M_prv_global_num_elements = M_mesh.numGlobalElements();

    this->distributedMeshProcessing(true);
}

void
FiniteElement::collectRootRestart(std::vector<double>& interp_elt_out, std::vector<double>& interp_nd_out)
{
    M_nb_var_element = 15 + M_tice.size();//15;

    if (M_rank == 0)
    {
        int num_elements_root = M_mesh_root.numTriangles();
        int tice_size = M_tice.size();

        interp_elt_out.resize(M_nb_var_element*num_elements_root);

        // std::cout<< "["<< M_rank << "] "<<"num_elements_root                   = "<< num_elements_root <<"\n";

        int tmp_nb_var = 0;

        for (int i=0; i<num_elements_root; ++i)
        {
            tmp_nb_var=0;
            //int ri = M_rmap_elements[i];

            // concentration
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_conc[i];
            tmp_nb_var++;

            // thickness
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_thick[i];
            tmp_nb_var++;

            // snow thickness
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_snow_thick[i];
            tmp_nb_var++;

            // integrated_stress1
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_sigma[3*i];
            tmp_nb_var++;

            // integrated_stress2
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_sigma[3*i+1];
            tmp_nb_var++;

            // integrated_stress3
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_sigma[3*i+2];
            tmp_nb_var++;

            // damage
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_damage[i];
            tmp_nb_var++;

            // ridge_ratio
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_ridge_ratio[i];
            tmp_nb_var++;

            // ridge_ratio
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_random_number[i];
            tmp_nb_var++;

            // sss
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_sss[i];
            tmp_nb_var++;

            // sst
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_sst[i];
            tmp_nb_var++;

            // tice1
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_tice[0][i];
            tmp_nb_var++;

            if ( interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_thermo_type == setup::ThermoType::WINTON )
            {
                // tice2
                interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_tice[1][i];
                tmp_nb_var++;

                // tice3
                interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_tice[2][i];
                tmp_nb_var++;
            }

            // h_thin
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_h_thin[i];
            tmp_nb_var++;

            // conc_thin
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_conc_thin[i];
            tmp_nb_var++;

            // hs_thin
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_hs_thin[i];
            tmp_nb_var++;

            // tsurf_thin
            interp_elt_out[M_nb_var_element*i+tmp_nb_var] = M_tsurf_thin[i];
            tmp_nb_var++;

            if(tmp_nb_var>M_nb_var_element)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }
    }

#if 1
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
#endif
}

void
FiniteElement::updateVelocity()
{
    M_VTMM = M_VTM;
    M_VTM = M_VT;
    M_VT = M_solution->container();

#if 0
    std::vector<double> speed_scaling;
    this->speedScaling(speed_scaling);

    // linear scaling of ice velocity
    for (int i=0; i<M_num_nodes; ++i)
    {
        M_VT[i] *= speed_scaling[i];
        M_VT[i+M_num_nodes] *= speed_scaling[i];
    }
#endif

#if 0
    M_speed_scaling = speed_scaling;

    double min_elt = *std::min_element(M_VT.begin(),M_VT.end());
    double max_elt = *std::max_element(M_VT.begin(),M_VT.end());

    // std::cout<<"----------------------------[" << M_rank <<"] " <<" VT MIN= "<< min_elt <<"\n";
    // std::cout<<"----------------------------[" << M_rank <<"] " <<" VT MAX= "<< max_elt <<"\n";

    M_comm.barrier();

    double gmin = boost::mpi::all_reduce(M_comm, min_elt, boost::mpi::minimum<double>());
    double gmax = boost::mpi::all_reduce(M_comm, max_elt, boost::mpi::maximum<double>());

    if (M_comm.rank()==0)
    {
        std::cout<<"----------------------------VT MIN= "<< gmin <<"\n";
        std::cout<<"----------------------------VT MAX= "<< gmax <<"\n";
    }
#endif

    // increment M_UT that is used for the drifters
    for (int nd=0; nd<M_UM.size(); ++nd)
    {
        M_UT[nd] += time_step*M_VT[nd]; // Total displacement (for drifters)
    }
}

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

            coef_Voce = (vm["simul.lin_drag_coef_water"].as<double>()+(quad_drag_coef_water*norm_Voce_ice));
            coef_Voce *= physical::rhow;

            norm_Vair_ice = std::hypot(M_VT[index_u]-M_wind [index_u],M_VT[index_v]-M_wind [index_v]);
            norm_Vair_ice = (norm_Vair_ice > norm_Vair_ice_min) ? (norm_Vair_ice):norm_Vair_ice_min;

            coef_Vair = (vm["simul.lin_drag_coef_air"].as<double>()+(quad_drag_coef_air*norm_Vair_ice));
            coef_Vair *= (physical::rhoa);

            M_VT[index_u] = ( coef_Vair*M_wind [index_u] + coef_Voce*M_ocean [index_u] ) / ( coef_Vair+coef_Voce );
            M_VT[index_v] = ( coef_Vair*M_wind [index_v] + coef_Voce*M_ocean [index_v] ) / ( coef_Vair+coef_Voce );

            // increment M_UT that is used for the drifters
            M_UT[index_u] += time_step*M_VT[index_u]; // Total displacement (for drifters)
            M_UT[index_v] += time_step*M_VT[index_v]; // Total displacement (for drifters)
        }
    }
}

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
}

void
FiniteElement::forcingAtmosphere()
{
    switch (M_atmosphere_type)
    {
        case setup::AtmosphereType::CONSTANT:
            M_wind=ExternalData(
                vm["simul.constant_wind_u"].as<double>(),
                vm["simul.constant_wind_v"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(vm["simul.constant_tair"].as<double>());
            M_external_data.push_back(&M_tair);

            M_mixrat=ExternalData(vm["simul.constant_mixrat"].as<double>());
            M_external_data.push_back(&M_mixrat);

            M_mslp=ExternalData(vm["simul.constant_mslp"].as<double>());
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(vm["simul.constant_Qsw_in"].as<double>());
            M_external_data.push_back(&M_Qsw_in);

            M_Qlw_in=ExternalData(vm["simul.constant_Qlw_in"].as<double>());
            M_external_data.push_back(&M_Qlw_in);

            M_snowfr=ExternalData(vm["simul.constant_snowfr"].as<double>());
            M_external_data.push_back(&M_snowfr);

            M_precip=ExternalData(vm["simul.constant_precip"].as<double>());
            M_external_data.push_back(&M_precip);

            M_dair=ExternalData(vm["simul.constant_dair"].as<double>());
            M_external_data.push_back(&M_dair);
        break;

        case setup::AtmosphereType::ASR:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_tair);

            M_mixrat=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_external_data.push_back(&M_mixrat);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_Qsw_in);

            M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            M_external_data.push_back(&M_Qlw_in);

            M_snowfr=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_external_data.push_back(&M_snowfr);

            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
            M_external_data.push_back(&M_precip);
        break;

        case setup::AtmosphereType::ERAi:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            /*
            variables[0] = tair;
            variables[1] = dair;
            variables[2] = mslp;
            variables[3] = Qsw_in;
            variables[4] = tcc;
            variables[5] = precip;
            */

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_tair);

            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_external_data.push_back(&M_dair);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_Qsw_in);

            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            M_external_data.push_back(&M_tcc);

            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_external_data.push_back(&M_precip);

            M_snowfall=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
            M_external_data.push_back(&M_snowfall);
        break;

        case setup::AtmosphereType::EC:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_tair);

            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_external_data.push_back(&M_dair);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_tcc);

            // Syl: The following two lines should be removed when approxSW will be implemented in Thermo()
            M_Qsw_in=ExternalData(vm["simul.constant_Qsw_in"].as<double>());
            M_external_data.push_back(&M_Qsw_in);

            M_precip=ExternalData(0.);
            M_external_data.push_back(&M_precip);
        break;

        case setup::AtmosphereType::EC2:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_tair);

            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_external_data.push_back(&M_dair);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_tcc);

            // Syl: The following two lines should be removed when approxSW will be implemented in Thermo()
            M_Qsw_in=ExternalData(vm["simul.constant_Qsw_in"].as<double>());
            M_external_data.push_back(&M_Qsw_in);

            M_precip=ExternalData(0.);
            M_external_data.push_back(&M_precip);
        break;

        case setup::AtmosphereType::EC_ERAi:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_tair);

            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_external_data.push_back(&M_dair);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_tcc);

            M_Qsw_in=ExternalData(&M_atmosphere_bis_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_Qsw_in);

            M_precip=ExternalData(&M_atmosphere_bis_elements_dataset,M_mesh,5,false,time_init);
            M_external_data.push_back(&M_precip);
        break;

        case setup::AtmosphereType::CFSR_HI:
        case setup::AtmosphereType::CFSR:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_tair);

            M_sphuma=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_external_data.push_back(&M_sphuma);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_Qsw_in);

            M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            M_external_data.push_back(&M_Qlw_in);

            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_external_data.push_back(&M_precip);

            M_snowfr=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
            M_external_data.push_back(&M_snowfr);
        break;

        default:
            std::cout << "invalid wind forcing"<<"\n";
            throw std::logic_error("invalid wind forcing");
    }
}

void
FiniteElement::forcingOcean()
{
    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
            M_ocean=ExternalData(
                vm["simul.constant_ocean_v"].as<double>(),
                vm["simul.constant_ocean_v"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ocean);

            M_ssh=ExternalData(vm["simul.constant_ssh"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ssh);

            M_ocean_temp=ExternalData(physical::ocean_freezing_temp);
            M_external_data.push_back(&M_ocean_temp);

            M_ocean_salt=ExternalData(physical::ocean_freezing_temp/physical::mu);
            M_external_data.push_back(&M_ocean_salt);

            M_mld=ExternalData(vm["simul.constant_mld"].as<double>());
            M_external_data.push_back(&M_mld);
            break;
        case setup::OceanType::TOPAZR: case setup::OceanType::TOPAZF:
            M_ocean=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 0, true,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ocean);

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ssh);

            M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
            M_external_data.push_back(&M_ocean_temp);

            M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
            M_external_data.push_back(&M_ocean_salt);

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            M_external_data.push_back(&M_mld);
            // SYL: there was a capping of the mld at minimum vm["simul.constant_mld"].as<double>()
            // but Einar said it is not necessary, so it is not implemented
    		break;

        case setup::OceanType::TOPAZR_ALTIMETER:
            M_ocean=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 0, true,
            time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ocean);

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
            time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ssh);

            M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
            M_external_data.push_back(&M_ocean_temp);

            M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
            M_external_data.push_back(&M_ocean_salt);

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            M_external_data.push_back(&M_mld);
            // SYL: there was a capping of the mld at minimum vm["simul.constant_mld"].as<double>()
            // but Einar said it is not necessary, so it is not implemented
            break;

        case setup::OceanType::TOPAZR_atrest:
            M_ocean=ExternalData(
                vm["simul.constant_ocean_v"].as<double>(),
                vm["simul.constant_ocean_v"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ocean);

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ssh);

            M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
            M_external_data.push_back(&M_ocean_temp);

            M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
            M_external_data.push_back(&M_ocean_salt);

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            M_external_data.push_back(&M_mld);
            // SYL: there was a capping of the mld at minimum vm["simul.constant_mld"].as<double>()
            // but Einar said it is not necessary, so it is not implemented
    		break;

        default:
            std::cout << "invalid ocean forcing"<<"\n";
            throw std::logic_error("invalid ocean forcing");
    }
}

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
}

void
FiniteElement::initIce()
{
    switch (M_ice_type)
    {
        case setup::IceType::CONSTANT:
            this->constantIce();
            break;
        case setup::IceType::CONSTANT_PARTIAL:
            this->constantIce();
            break;
        case setup::IceType::TARGET:
            this->targetIce();
            break;
        case setup::IceType::BINARY:
            this->binaryIce();
            break;
        case setup::IceType::TOPAZ4:
            this->topazIce();
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
        case setup::IceType::SMOS:
            this->smosIce();
            break;
        default:
            std::cout << "invalid initialization of the ice"<<"\n";
            throw std::logic_error("invalid initialization of the ice");
    }

    double hi, conc_tot, weight_conc;
    // Consistency check for the slab ocean + initialization of the ice temperature (especially for Winton)
    for ( int i=0; i<M_num_elements; i++ )
    {
        hi=0.;
        if(M_conc[i]>0.)
            hi = M_thick[i]/M_conc[i];

        if ( M_conc[i] < physical::cmin || hi < physical::hmin)
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }

        hi = 0.;
        if(M_conc_thin[i]>0.)
            hi = M_h_thin[i]/M_conc_thin[i];

        if ( M_conc_thin[i] < physical::cmin || hi < physical::hmin)
        {
            M_conc_thin[i]=0.;
            M_h_thin[i]=0.;
            M_hs_thin[i]=0.;
        }

        conc_tot=M_conc[i]+M_conc_thin[i];
        weight_conc=std::min(1.,conc_tot*2.);
        if(conc_tot>0.)
            M_sst[i] = -M_sss[i]*physical::mu*weight_conc+M_sst[i]*(1.-weight_conc);
        //M_sst[i] = -M_sss[i]*physical::mu;//*M_conc[i]+M_ocean_temp[i]*(1.-M_conc[i]);

        if ( M_snow_thick[i] > 0. )
            M_tice[0][i] = std::min(0., M_tair[i]);
        else
            M_tice[0][i] = std::min(-physical::mu*physical::si, M_tair[i]);
    }

    if ( M_thermo_type == setup::ThermoType::WINTON )
    {
        for ( int i=0; i<M_num_elements; i++ )
        {
            double Tbot  = -physical::mu*M_sss[i];
            if ( M_thick[i] > 0. )
            {
                // Just a linear interpolation between bottom and snow-ice interface (i.e. a zero layer model)
                double deltaT = (Tbot - M_tice[0][i] ) / ( 1. + physical::ki*M_snow_thick[i]/(physical::ks*M_thick[i]) );
                double slope = -deltaT/M_thick[i];
                M_tice[1][i] = Tbot + 3*slope*M_thick[i]/4;
                M_tice[2][i] = Tbot +   slope*M_thick[i]/4;
            } else {
                M_tice[1][i] = Tbot;
                M_tice[2][i] = Tbot;
            }
        }
    }

#if 0
    // It's nice to initialise the ice temperature (especially for Winton)
    for ( int i=0; i<M_num_elements; i++ )
        if ( M_snow_thick[i] > 0. )
            M_tice[0][i] = std::min(0., M_tair[i]);
        else
            M_tice[0][i] = std::min(-physical::mu*physical::si, M_tair[i]);

    if ( M_thermo_type == setup::ThermoType::WINTON )
    {
        for ( int i=0; i<M_num_elements; i++ )
        {
            double Tbot  = -physical::mu*M_sss[i];
            if ( M_thick[i] > 0. )
            {
                // Just a linear interpolation between bottom and snow-ice interface (i.e. a zero layer model)
                double deltaT = (Tbot - M_tice[0][i] ) / ( 1. + physical::ki*M_snow_thick[i]/(physical::ks*M_thick[i]) );
                double slope = -deltaT/M_thick[i];
                M_tice[1][i] = Tbot + 3*slope*M_thick[i]/4;
                M_tice[2][i] = Tbot +   slope*M_thick[i]/4;
            } else {
                M_tice[1][i] = Tbot;
                M_tice[2][i] = Tbot;
            }
        }
    }
#endif
}

void
FiniteElement::constantIce()
{
	LOG(DEBUG) <<"Constant Ice\n";
    std::fill(M_conc.begin(), M_conc.end(), vm["simul.init_concentration"].as<double>());
    std::fill(M_thick.begin(), M_thick.end(), vm["simul.init_thickness"].as<double>());
    std::fill(M_snow_thick.begin(), M_snow_thick.end(), vm["simul.init_snow_thickness"].as<double>());
    std::fill(M_damage.begin(), M_damage.end(), 0.);

    if (M_ice_type==setup::IceType::CONSTANT_PARTIAL)
    {
        auto Bx = M_mesh.coordX();//xmin,xmax from nodes
        double xmin = *std::min_element(Bx.begin(),Bx.end());
        double xmax = *std::max_element(Bx.begin(),Bx.end());
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

        Bx = M_mesh.bCoordX();//set conc, etc on elements
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

}

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

        //M_conc[i]  = std::max(vm["simul.init_concentration"].as<double>()*tmp_var,cmin);
	//	M_thick[i] = vm["simul.init_thickness"].as<double>()*M_conc[i];
	//	M_snow_thick[i] = vm["simul.init_snow_thickness"].as<double>()*M_conc[i];

        M_conc[i]  = vm["simul.init_concentration"].as<double>();

	if(i==10)
		M_conc[i]=0.;

	M_thick[i] = vm["simul.init_thickness"].as<double>()*M_conc[i];
	M_snow_thick[i] = vm["simul.init_snow_thickness"].as<double>()*M_conc[i];

	M_damage[i]=0.;

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_conc_thin[i]  = vm["simul.init_thin_conc"].as<double>();

            M_h_thin[i]     = (vm["simul.h_thin_min"].as<double>()+(vm["simul.h_thin_max"].as<double>()-vm["simul.h_thin_min"].as<double>())/2.)*M_conc_thin[i];

            M_hs_thin[i]    = vm["simul.init_snow_thickness"].as<double>()*M_conc_thin[i];
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
}

void
FiniteElement::binaryIce()
{
    std::fstream input;
    std::string filename;
    int length;

    // First concentration
    filename = Environment::nextsimDir().string() + "/data/initConc.dat";
    input.open(filename, std::fstream::in);

    // get length of file:
    input.seekg (0, input.end);
    length = input.tellg();
    input.seekg (0, input.beg);
    if ( length != M_num_elements*sizeof(double) )
        throw std::runtime_error("Couldn't read the correct number of elements from " + filename +
                ". File length is " + std::to_string(length) + " while M_num_elements is " + std::to_string(M_num_elements) + ".\n");

    input.read((char*) &M_conc[0], M_num_elements*sizeof(double));
    input.close();

    // Then thickness
    filename = Environment::nextsimDir().string() + "/data/initThick.dat";
    input.open(filename, std::fstream::in);

    // get length of file:
    input.seekg (0, input.end);
    length = input.tellg();
    input.seekg (0, input.beg);
    if ( length != M_num_elements*sizeof(double) )
        throw std::runtime_error("Couldn't read the correct number of elements from " + filename +
                ". File length is " + std::to_string(length) + " while M_num_elements is " + std::to_string(M_num_elements) + ".\n");

    input.read((char*) &M_thick[0], M_num_elements*sizeof(double));
    input.close();

    // Finally snow thickness
    filename = Environment::nextsimDir().string() + "/data/initSnow.dat";
    input.open(filename, std::fstream::in);

    // get length of file:
    input.seekg (0, input.end);
    length = input.tellg();
    input.seekg (0, input.beg);
    if ( length != M_num_elements*sizeof(double) )
        throw std::runtime_error("Couldn't read the correct number of elements from " + filename +
                ". File length is " + std::to_string(length) + " while M_num_elements is " + std::to_string(M_num_elements) + ".\n");

    input.read((char*) &M_snow_thick[0], M_num_elements*sizeof(double));
    input.close();

    // Make sure damage is zero and do something for thin ice
    for (int i=0; i<M_num_elements; ++i)
    {
		M_damage[i]=0.;
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_conc_thin[i]=std::min(1.-M_conc[i], 0.2*M_conc[i]);
            M_h_thin[i]=M_conc_thin[i]*(h_thin_min+0.5*(h_thin_max-h_thin_min));
        }
    }
}

void
FiniteElement::topazIce()
{
    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,"init - TOPAZ ice");
    M_external_data_tmp.resize(0);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		tmp_var=std::min(1.,M_init_conc[i]);
		M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
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
}

void
FiniteElement::topazIceOsisafIcesat()
{
    external_data M_topaz_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);

    external_data M_topaz_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);

    external_data M_topaz_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);

    external_data M_osisaf_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);

    external_data M_osisaf_conc=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,0,false,time_init);

    external_data M_icesat_thick=ExternalData(&M_ice_icesat_elements_dataset,M_mesh,0,false,time_init);

    external_data M_amsre_conc=ExternalData(&M_ice_amsre_elements_dataset,M_mesh,0,false,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_topaz_conc);
    M_external_data_tmp.push_back(&M_topaz_thick);
    M_external_data_tmp.push_back(&M_topaz_snow_thick);
    M_external_data_tmp.push_back(&M_osisaf_type);
    M_external_data_tmp.push_back(&M_osisaf_conc);
    M_external_data_tmp.push_back(&M_icesat_thick);
    M_external_data_tmp.push_back(&M_amsre_conc);
    this->checkReloadDatasets(M_external_data_tmp,time_init,"init - TOPAZ ice");
    M_external_data_tmp.resize(0);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		tmp_var=std::min(1.,M_topaz_conc[i]);
		M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
		tmp_var=M_topaz_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
        double hi_topaz=M_thick[i]/M_conc[i];
        hi_topaz=std::min(hi_topaz,3.); // TOPAZ thickness is only used for FYI and then capped to 3 m, to avoid unrealistic thickness where concentration is low

        if(M_conc[i]>0.) // use osisaf only where topaz says there is ice to avoid near land issues and fake concentration over the ocean
            M_conc[i]=M_osisaf_conc[i];

		tmp_var=M_topaz_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
        if(M_conc[i]<M_topaz_conc[i])
            M_snow_thick[i] *= M_conc[i]/M_topaz_conc[i];


        //M_type[i]==1. // No ice
        //M_type[i]==2. // First-Year ice
        //M_type[i]==3. // Multi-Year ice
        //M_type[i]==4. // Mixed
        double ratio_FYI=0.3;
        double ratio_MYI=0.9;
        double ratio_Mixed=0.5*(ratio_FYI+ratio_MYI);

        double thick_FYI=hi_topaz;
        double thick_MYI=std::max(M_icesat_thick[i],hi_topaz);
        double thick_Mixed=0.5*(thick_FYI+thick_MYI);

        if((M_thick[i]>0.)&&(M_conc[i])>0.2)
        {
            if(M_osisaf_type[i]<=1.)
            {
                M_ridge_ratio[i]=0.;
                M_thick[i]=thick_FYI;
            }
            if(M_osisaf_type[i]>1. && M_osisaf_type[i]<=2.)
            {
                M_ridge_ratio[i]=(M_osisaf_type[i]-1.)*ratio_FYI;
                M_thick[i]      =thick_FYI;
            }
            if(M_osisaf_type[i]>2. && M_osisaf_type[i]<=3.)
            {
                M_ridge_ratio[i]=(1.-(M_osisaf_type[i]-2.))*ratio_FYI + (M_osisaf_type[i]-2.)*ratio_MYI;
                M_thick[i]      =(1.-(M_osisaf_type[i]-2.))*thick_FYI + (M_osisaf_type[i]-2.)*thick_MYI;
            }
            if(M_osisaf_type[i]>3. && M_osisaf_type[i]<=4.)
            {
                M_ridge_ratio[i]=(1.-(M_osisaf_type[i]-3.))*ratio_MYI + (M_osisaf_type[i]-3.)*ratio_Mixed;
                M_thick[i]      =(1.-(M_osisaf_type[i]-3.))*thick_MYI + (M_osisaf_type[i]-3.)*thick_Mixed;
            }
            if(M_osisaf_type[i]>4.)
            {
                M_ridge_ratio[i]=ratio_Mixed;
                M_thick[i]=thick_Mixed;
            }
        }
        else
        {
            M_ridge_ratio[i]=0.;
            M_thick[i]=hi_topaz;
        }

        M_ridge_ratio[i]=M_ridge_ratio[i]*M_conc[i]; // Icesat gives the actual thickness (see "Uncertainties in Arctic sea ice thickness and volume: new estimates and implications for trends")
        M_thick[i]=M_thick[i]*M_conc[i];

        //if either c or h equal zero, we set the others to zero as well
        double hi=M_thick[i]/M_conc[i];
        if ( M_conc[i] < 0.01 || hi < physical::hmin )
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
            M_ridge_ratio[i]=0.;
        }

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_conc_thin[i]=std::max(M_amsre_conc[i]-M_conc[i],0.);
            M_h_thin[i]=M_conc_thin[i]*(h_thin_min+0.5*(h_thin_max-h_thin_min));
        }

		M_damage[i]=0.;
	}
}

void
FiniteElement::topazForecastIce()
{
    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,"init - TOPAZ ice forecast");
    M_external_data_tmp.resize(0);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		tmp_var=std::min(1.,M_init_conc[i]);
		M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
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
}

void
FiniteElement::topazForecastAmsr2Ice()
{
    double real_thickness, init_conc_tmp;

    external_data M_conc_amsr2=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
    //M_conc_amsr2.check_and_reload(M_mesh,time_init);

    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_conc_amsr2);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init - TOPAZ ice forecast + AMSR2");
    M_external_data_tmp.resize(0);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
        double uncertainty;
        if(M_conc_amsr2[i]<0.1)
            uncertainty=0.1;
		else
            uncertainty=0.05;

        double diff_mod_obs=M_conc_amsr2[i]-M_init_conc[i];
        if(std::abs(diff_mod_obs)>=uncertainty)
            M_conc[i] = std::min(1.,M_conc_amsr2[i]-(diff_mod_obs)/std::abs(diff_mod_obs)*uncertainty/2.);
        else
            M_conc[i] = std::min(1.,M_init_conc[i]);

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
            real_thickness=M_thick[i]/init_conc_tmp;
            M_thick[i]=real_thickness*M_conc[i];
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
            M_thick[i]=0.1*M_conc[i];

		M_damage[i]=0.;
	}
}

void
FiniteElement::topazForecastAmsr2OsisafIce()
{
    double real_thickness, init_conc_tmp;

    external_data M_osisaf_conc=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,0,false,time_init);

    external_data M_osisaf_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);

    external_data M_amsr2_conc=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);

    external_data M_topaz_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);

    external_data M_topaz_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);

    external_data M_topaz_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_osisaf_conc);
    M_external_data_tmp.push_back(&M_osisaf_type);
    M_external_data_tmp.push_back(&M_amsr2_conc);
    M_external_data_tmp.push_back(&M_topaz_conc);
    M_external_data_tmp.push_back(&M_topaz_thick);
    M_external_data_tmp.push_back(&M_topaz_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init - TOPAZ ice forecast + AMSR2 + OSISAF");
    M_external_data_tmp.resize(0);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		tmp_var=std::min(1.,M_topaz_conc[i]);
		M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
		tmp_var=M_topaz_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
        double hi_topaz=M_thick[i]/M_conc[i];
        hi_topaz=std::min(hi_topaz,3.); // TOPAZ thickness is only used for FYI and then capped to 3 m, to avoid unrealistic thickness where concentration is low

        if(M_conc[i]>0.) // use osisaf only where topaz says there is ice to avoid near land issues and fake concentration over the ocean
            M_conc[i]=M_osisaf_conc[i];

        if(M_amsr2_conc[i]<M_conc[i]) // AMSR2 is higher resolution and see small opening that would not be see in OSISAF
            M_conc[i]=M_amsr2_conc[i];

		tmp_var=M_topaz_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
        if(M_conc[i]<M_topaz_conc[i])
            M_snow_thick[i] *= M_conc[i]/M_topaz_conc[i];


        //M_type[i]==1. // No ice
        //M_type[i]==2. // First-Year ice
        //M_type[i]==3. // Multi-Year ice
        //M_type[i]==4. // Mixed
        double ratio_FYI=0.3;
        double ratio_MYI=0.9;
        double ratio_Mixed=0.5*(ratio_FYI+ratio_MYI);

        double thick_FYI=hi_topaz;
        double thick_MYI=1.5*hi_topaz;
        double thick_Mixed=0.5*(thick_FYI+thick_MYI);

        if( (M_thick[i]>0.) && (M_conc[i])>0.2 )
        {

            if(M_mesh_filename.find("kara") != std::string::npos)
            {
                LOG(DEBUG) <<"Type information is not used for the kara mesh, we assume there is only FYI\n";
                M_ridge_ratio[i]=ratio_FYI;
                M_thick[i]=thick_FYI;
            }
            else
            {
            if(M_osisaf_type[i]<=1.)
            {
                M_ridge_ratio[i]=0.;
                M_thick[i]=thick_FYI;
            }
            if(M_osisaf_type[i]>1. && M_osisaf_type[i]<=2.)
            {
                M_ridge_ratio[i]=(M_osisaf_type[i]-1.)*ratio_FYI;
                M_thick[i]      =thick_FYI;
            }
            if(M_osisaf_type[i]>2. && M_osisaf_type[i]<=3.)
            {
                M_ridge_ratio[i]=(1.-(M_osisaf_type[i]-2.))*ratio_FYI + (M_osisaf_type[i]-2.)*ratio_MYI;
                M_thick[i]      =(1.-(M_osisaf_type[i]-2.))*thick_FYI + (M_osisaf_type[i]-2.)*thick_MYI;
            }
            if(M_osisaf_type[i]>3. && M_osisaf_type[i]<=4.)
            {
                M_ridge_ratio[i]=(1.-(M_osisaf_type[i]-3.))*ratio_MYI + (M_osisaf_type[i]-3.)*ratio_Mixed;
                M_thick[i]      =(1.-(M_osisaf_type[i]-3.))*thick_MYI + (M_osisaf_type[i]-3.)*thick_Mixed;
            }
            if(M_osisaf_type[i]>4.)
            {
                M_ridge_ratio[i]=ratio_Mixed;
                M_thick[i]=thick_Mixed;
            }
            }
        }
        else
        {
            M_ridge_ratio[i]=0.;
            M_thick[i]=hi_topaz;
        }

        M_ridge_ratio[i]=M_ridge_ratio[i]*M_conc[i];
        // Icesat gives the actual thickness (see "Uncertainties in Arctic sea ice thickness and volume: new estimates and implications for trends")
        M_thick[i]=M_thick[i]*M_conc[i];

        //if either c or h equal zero, we set the others to zero as well
        double hi=M_thick[i]/M_conc[i];
        if ( M_conc[i] < 0.01 || hi < physical::hmin )
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
            M_ridge_ratio[i]=0.;
        }

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_conc_thin[i]=std::max(M_amsr2_conc[i]-M_conc[i],0.);
            M_h_thin[i]=M_conc_thin[i]*(h_thin_min+0.5*(h_thin_max-h_thin_min));
        }

		M_damage[i]=0.;
	}
}

void
FiniteElement::piomasIce()
{
    external_data M_init_conc=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,0,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,1,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,2,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init ice - PIOMAS");
    M_external_data_tmp.resize(0);

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
}

void
FiniteElement::topazAmsreIce()
{
    double real_thickness, init_conc_tmp;

    external_data M_conc_amsre=ExternalData(&M_ice_amsre_elements_dataset,M_mesh,0,false,time_init);
    //M_conc_amsre.check_and_reload(M_mesh,time_init);

    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_conc_amsre);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init ice - TOPAZ + AMSR-E");
    M_external_data_tmp.resize(0);

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
            real_thickness=M_thick[i]/init_conc_tmp;
            M_thick[i]=real_thickness*M_conc[i];
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
}

void
FiniteElement::topazAmsr2Ice()
{
    double real_thickness, init_conc_tmp;

    external_data M_conc_amsr2=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
    //M_conc_amsr2.check_and_reload(M_mesh,time_init);

    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_conc_amsr2);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init ice - TOPAZ + AMSR2");
    M_external_data_tmp.resize(0);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		M_conc[i] = std::min(1.,M_conc_amsr2[i]);

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
            real_thickness=M_thick[i]/init_conc_tmp;
            M_thick[i]=real_thickness*M_conc[i];
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
}

void
FiniteElement::cs2SmosIce()
{
    external_data M_init_conc=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,0,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,1,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);
    //M_type.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_type);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init ice - CS2 + SMOS");
    M_external_data_tmp.resize(0);

    warrenClimatology();

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
}

void
FiniteElement::cs2SmosAmsr2Ice()
{
    external_data M_init_conc=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,0,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_amsr2_conc=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
    //M_conc_amsr2.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,1,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);
    //M_type.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_type);
    M_external_data_tmp.push_back(&M_amsr2_conc);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init ice - CS2 + SMOS + AMSR2");
    M_external_data_tmp.resize(0);

    warrenClimatology();

    double tmp_var, correction_factor_warren;
    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_var=std::min(1.,M_init_conc[i]);
        M_conc[i] = tmp_var;
        tmp_var=M_init_thick[i];
        M_thick[i] = tmp_var ;
        if(M_amsr2_conc[i]<M_conc[i]) // AMSR2 is higher resolution and see small opening that would not be see in OSISAF
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
}

void
FiniteElement::smosIce()
{
    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_smos_elements_dataset,M_mesh,0,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    boost::gregorian::date dt = Nextsim::parse_date(time_init);
    int month_id=dt.month().as_number(); // 1 for January, 2 for February, and so on. This will be used to compute the snow from Warren climatology

    std::cout << "month_id: " << month_id <<"\n";

    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init ice - SMOS");
    M_external_data_tmp.resize(0);

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
}

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

}

void
FiniteElement::initThermodynamics()
{
    // std::fill(damage.begin(), damage.end(), 0.);

#if 0
    for (int i=0; i<M_num_elements; ++i)
    {
        damage[i] = 1.0 - conc[i];
    }
#endif
}

void
FiniteElement::initDrifter()
{
    if(M_use_equallyspaced_drifters)
        this->equallySpacedDrifter();

    if(M_use_iabp_drifters)
        this->initIABPDrifter();

    if(M_use_rgps_drifters)
        this->initRGPSDrifters();

    if(M_use_osisaf_drifters)
        this->initOSISAFDrifters();
}

void
FiniteElement::coriolis()
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
}

void
FiniteElement::bathymetry()
{
    switch (M_bathymetry_type)
    {
        case setup::BathymetryType::CONSTANT:
            M_element_depth=ExternalData(vm["simul.constant_bathymetry"].as<double>());
            M_external_data.push_back(&M_element_depth);
            break;
        case setup::BathymetryType::ETOPO:
            M_element_depth=ExternalData(&M_bathymetry_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_element_depth);
            break;
        default:
            std::cout << "invalid bathymetry"<<"\n";
            throw std::logic_error("invalid bathymetry");
    }
}

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
}

// A simple function to output the drifters in the model, IABP or otherwise
// The output could well be prettier!
void
FiniteElement::outputDrifter(std::fstream &drifters_out)
{
    // Initialize the map
    mapx_class *map;
    std::string configfile = (boost::format( "%1%/%2%/%3%" )
                              % Environment::nextsimDir().string()
                              % "data"
                              % vm["mesh.mppfile"].as<std::string>()
                              ).str();

    std::vector<char> str(configfile.begin(), configfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

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

    // Interpolate the total displacement and concentration onto the drifter positions
    int nb_var=2;
    std::vector<double> interp_drifter_in(nb_var*M_mesh.numNodes());

    // Interpolate the velocity
    for (int i=0; i<M_mesh.numNodes(); ++i)
    {
        interp_drifter_in[nb_var*i]   = M_UT[i];
        interp_drifter_in[nb_var*i+1] = M_UT[i+M_mesh.numNodes()];
    }

    double* interp_drifter_out;
    InterpFromMeshToMesh2dx(&interp_drifter_out,
        &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
        M_mesh.numNodes(),M_mesh.numTriangles(),
        &interp_drifter_in[0],
        M_mesh.numNodes(),nb_var,
        &drifter_X[0],&drifter_Y[0],M_iabp_drifters.size(),
        true, 0.);

    // Loop over the map and output
    j=0;
    boost::gregorian::date           date = Nextsim::parse_date( current_time );
    boost::posix_time::time_duration time = Nextsim::parse_time( current_time );
    for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); ++it )
    {
        double lat, lon;
        inverse_mapx(map,it->second[0]+interp_drifter_out[nb_var*j],it->second[1]+interp_drifter_out[nb_var*j+1],&lat,&lon);
        j++;

        drifters_out << setw(4) << date.year()
            << " " << setw( 2) << date.month().as_number()
            << " " << setw( 2) << date.day().as_number()
            << " " << setw( 2) << time.hours()
            << " " << setw(16) << it->first
            << fixed << setprecision(5)
            << " " << setw( 8) << lat
            << " " << setw(10) << lon << "\n";
    }

    xDelete<double>(interp_drifter_out);
    close_mapx(map);
}

// Add the buoys that have been put into the ice and remove dead ones
void
FiniteElement::updateIABPDrifter()
{
    // Initialize the map
    mapx_class *map;
    std::string configfile = (boost::format( "%1%/%2%/%3%" )
                              % Environment::nextsimDir().string()
                              % "data"
                              % vm["mesh.mppfile"].as<std::string>()
                              ).str();

    std::vector<char> str(configfile.begin(), configfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    // Read the current buoys from file
    int pos;    // To be able to rewind one line
    double time = current_time;
    std::vector<int> keepers;
    while ( time == current_time )
    {
        // Remember where we were
        pos = M_iabp_file.tellg();

        // Read the next line
        int year, month, day, hour, number;
        double lat, lon;
        M_iabp_file >> year >> month >> day >> hour >> number >> lat >> lon;
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

        }
    }
    close_mapx(map);

    // Go through the M_iabp_drifters map and throw out the ones which IABP doesn't
    // report as being in the ice anymore
    for ( auto model = M_iabp_drifters.begin(); model != M_iabp_drifters.end(); /* ++model is not allowed here, because we use 'erase' */ )
    {
        bool keep = false;
        // Check against all the buoys we want to keep
        for ( auto obs = keepers.begin(); obs != keepers.end(); ++obs )
        {
            if ( model->first == *obs )
            {
                keep = true;
                break;
            }
        }

        // Delete or advance the iterator
        if ( ! keep )
            model = M_iabp_drifters.erase(model);
        else
            ++model;
    }
}

// Initialise by reading all the data from '79 up to time_init
// This is too slow, but only happens once so I won't try to optimise that for now
void
FiniteElement::initIABPDrifter()
{
    std::string filename = Environment::nextsimDir().string() + "/data/IABP_buoys.txt";
    M_iabp_file.open(filename, std::fstream::in);
    if ( ! M_iabp_file.good() )
        throw std::runtime_error("File not found: " + filename);

    int pos;    // To be able to rewind one line
    double time = from_date_string("1979-01-01");
    while ( time < time_init )
    {
        // Remember where we were
        pos = M_iabp_file.tellg();

        // Read the next line
        int year, month, day, hour, number;
        double lat, lon;
        M_iabp_file >> year >> month >> day >> hour >> number >> lat >> lon;
        std::string date = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);

        time = from_date_string(date) + hour/24.;
    }

    // We must rewind one line so that updateIABPDrifter works correctly
    M_iabp_file.seekg(pos);
}

void
FiniteElement::equallySpacedDrifter()
{
    // M_drifters = Drifters(1e3*vm["simul.drifter_spacing"].as<double>(), M_mesh, M_conc, vm["simul.drifter_climit"].as<double>());
    // M_drifters.initNetCDF(M_export_path+"/Drifters_", current_time);
    // M_drifters.appendNetCDF(current_time, M_mesh, M_UT);
    M_equallyspaced_drifters = Drifters(1e3*vm["simul.drifter_spacing"].as<double>(), M_mesh, M_conc, vm["simul.drifter_climit"].as<double>());
    M_equallyspaced_drifters.initNetCDF(M_export_path+"/Drifters_", current_time);
    M_equallyspaced_drifters.appendNetCDF(current_time, M_mesh, M_UT);
}

void
FiniteElement::initRGPSDrifters()
{
    M_rgps_drifters = Drifters();
}

void
FiniteElement::updateRGPSDrifters()
{
    std::string time_str = vm["simul.RGPS_time_init"].as<std::string>();
    double RGPS_time_init = from_date_time_string(time_str);

    std::string filename = Environment::nextsimDir().string() + "/data/RGPS_" + time_str + ".txt";
    M_rgps_drifters = Drifters(filename, M_mesh, M_conc, vm["simul.drifter_climit"].as<double>(),RGPS_time_init);

    M_rgps_drifters.initNetCDF(M_export_path+"/RGPS_Drifters_", current_time);
    M_rgps_drifters.appendNetCDF(current_time, M_mesh, M_UT);
}

void
FiniteElement::initOSISAFDrifters()
{
    M_osisaf_drifters.resize(2);
}

void
FiniteElement::importBamg(BamgMesh const* bamg_mesh)
{
    std::vector<point_type> mesh_nodes;
    //std::vector<element_type> mesh_edges;
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

#if 0
    int numVertices = 2;
    std::vector<int> edges(numVertices);

    for (int edg=0; edg<bamg_mesh->EdgesSize[0]; ++edg)
    {
        edges[0] = bamg_mesh->Edges[3*edg];
        edges[1] = bamg_mesh->Edges[3*edg+1];
        //std::cout<< "Edges= "<< bamg_mesh->Edges[3*edg+2] <<"\n";

        element_type gmshElt( edg,
                              type,
                              physical,
                              elementary,
                              numVertices,
                              edges );

        //mesh_edges.insert(std::make_pair(edg,gmshElt));
        mesh_edges.push_back(gmshElt);
    }
#endif

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

        //mesh_triangles.insert(std::make_pair(tr,gmshElt));
        //mesh_triangles.push_back(gmshElt);
        mesh_triangles[tr] = gmshElt;
    }

    std::cout<<"\n";
    std::cout<<"INFO: Previous  NumNodes     = "<< M_mesh_root.numNodes() <<"\n";
    std::cout<<"INFO: Previous  NumTriangles = "<< M_mesh_root.numTriangles() <<"\n";
    //std::cout<<"INFO: Previous  NumEdges     = "<< M_mesh_root.numEdges() <<"\n";

    M_mesh_previous_root = M_mesh_root;
    //M_mesh_root = mesh_type_root(mesh_nodes,mesh_edges,mesh_triangles);
    //M_mesh_root = mesh_type_root(mesh_nodes,mesh_triangles);
    M_mesh_root.update(mesh_nodes,mesh_triangles);
    //M_mesh.writeTofile("out.msh");

    std::cout<<"\n";
    std::cout<<"INFO: Current  NumNodes      = "<< M_mesh_root.numNodes() <<"\n";
    std::cout<<"INFO: Current  NumTriangles  = "<< M_mesh_root.numTriangles() <<"\n";
    //std::cout<<"INFO: Current  NumEdges      = "<< M_mesh_root.numEdges() <<"\n";
    std::cout<<"\n";
}

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

}

void
FiniteElement::exportResults(int step, bool export_mesh)
{
    // velocity
    std::vector<double> vt_local(2*M_local_ndof,0.);
    for (int i=0; i<M_local_ndof; ++i)
    {
        vt_local[2*i] = M_VT[i];
        vt_local[2*i+1] = M_VT[i+M_num_nodes];
    }

    std::vector<int> sizes_nodes = M_sizes_nodes;
    std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = 2*f; });

    // send displacement vector to the root process (rank 0)
    std::vector<double> vt_root;

    chrono.restart();
    if (M_rank == 0)
    {
        vt_root.resize(2*M_ndof);
        boost::mpi::gatherv(M_comm, vt_local, &vt_root[0], sizes_nodes, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, vt_local, 0);
    }

    // fields defined on mesh elements
#if 1
    std::vector<int> sizes_elements = M_sizes_elements;

    // ELEMENT INTERPOLATION With Cavities
    int nb_var_element=7;//15;
    std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });
    std::vector<double> interp_elt_in_local(nb_var_element*M_local_nelements);

    //std::cout<<"ELEMENT: Interp starts\n";

    int tmp_nb_var=0;
    for (int i=0; i<M_local_nelements; ++i)
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
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i]/**M_thick[i]*/;
        tmp_nb_var++;

        // integrated_stress2
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i+1]/**M_thick[i]*/;
        tmp_nb_var++;

        // integrated_stress3
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i+2]/**M_thick[i]*/;
        tmp_nb_var++;

        // compliance
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_damage[i];
        tmp_nb_var++;

#if 0
        // h_ridged_thin_ice
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_h_ridged_thin_ice[i];
        tmp_nb_var++;

        // h_ridged_thick_ice
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_h_ridged_thick_ice[i];
        tmp_nb_var++;

        // random_number
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_random_number[i];
        tmp_nb_var++;

        // Ice surface temperature
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_tsurf[i];
        tmp_nb_var++;

        // thin ice thickness
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_h_thin[i];
        tmp_nb_var++;

        // snow on thin ice
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_hs_thin[i];
        tmp_nb_var++;

        // Ice surface temperature for thin ice
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_tsurf_thin[i];
        tmp_nb_var++;
#endif

        if(tmp_nb_var>nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }
    }

    std::vector<double> interp_in_elements;
    //M_comm.barrier();
    if (M_rank == 0)
    {
        interp_in_elements.resize(nb_var_element*M_mesh_root.numTriangles());
        boost::mpi::gatherv(M_comm, interp_elt_in_local, &interp_in_elements[0], sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, interp_elt_in_local, 0);
    }

    std::vector<double> conc_root(M_mesh_root.numTriangles());
    std::vector<double> thick_root(M_mesh_root.numTriangles());
    std::vector<double> thick_snow(M_mesh_root.numTriangles());
    std::vector<double> stress1(M_mesh_root.numTriangles());
    std::vector<double> stress2(M_mesh_root.numTriangles());
    std::vector<double> stress3(M_mesh_root.numTriangles());
    std::vector<double> damage(M_mesh_root.numTriangles());
    // std::vector<double> ridged_thin(M_mesh_root.numTriangles());
    // std::vector<double> ridged_thick(M_mesh_root.numTriangles());
    // std::vector<double> random(M_mesh_root.numTriangles());
    // std::vector<double> tsurf(M_mesh_root.numTriangles());
    // std::vector<double> hthin(M_mesh_root.numTriangles());
    // std::vector<double> hsthin(M_mesh_root.numTriangles());
    // std::vector<double> tsurfthin(M_mesh_root.numTriangles());

    if (M_rank == 0)
    {
        tmp_nb_var=0;
        //auto rmap_elements = M_mesh.mapElements();
        //auto interp_in_elements_nrd = interp_in_elements;

        for (int i=0; i<M_mesh_root.numTriangles(); ++i)
        {
            tmp_nb_var=0;
            //int ri = rmap_elements.left.find(i+1)->second-1;
            int ri = M_rmap_elements[i];

            // for (int j=0; j<nb_var_element; ++j)
            // {
            //     interp_in_elements[nb_var_element*i+j] = interp_in_elements_nrd[nb_var_element*ri+j];
            // }

            // concentration
            conc_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // thickness
            thick_root[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // snow thickness
            thick_snow[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress1
            stress1[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress2
            stress2[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress3
            stress3[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // damage
            damage[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

#if 0
            // h_ridged_thin_ice
            ridged_thin[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // h_ridged_thick_ice
            ridged_thick[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // random_number
            random[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // Ice surface temperature
            tsurf[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // thin ice thickness
            hthin[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // snow on thin ice
            hsthin[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // Ice surface temperature for thin ice
            tsurfthin[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;
#endif

            if(tmp_nb_var>nb_var_element)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }
    }

#endif

    if (M_rank == 0)
    {
        Exporter exporter("float");
        std::string fileout;

        if (export_mesh)
        {
            fileout = (boost::format( "%1%/matlab/mesh_%2%_%3%.bin" )
                       % Environment::nextsimDir().string()
                       % M_comm.size() /*M_rank*/
                       % step ).str();

            std::cout<<"MESH BINARY: Exporter Filename= "<< fileout <<"\n";

            // move the mesh for the export
            // M_mesh.move(M_UM,1.);

            std::fstream meshbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
            if ( ! meshbin.good() )
                throw std::runtime_error("Cannot write to file: " + fileout);

#if 0
            auto rmap_nodes = M_mesh.mapNodes();
            auto rmap_elements = M_mesh.mapElements();

            auto meshr = M_mesh_root;
            meshr.reorder(rmap_nodes,rmap_elements);
            //exporter.writeMesh(meshbin, M_mesh);
            //exporter.writeMesh(meshbin, M_mesh_root);
            exporter.writeMesh(meshbin, meshr);
            meshbin.close();
#endif
            exporter.writeMesh(meshbin, M_mesh_root);
            meshbin.close();

            // move it back after the export
            //M_mesh.move(M_UM,-1.);

            fileout = (boost::format( "%1%/matlab/mesh_%2%_%3%.dat" )
                       % Environment::nextsimDir().string()
                       % M_comm.size()
                       % step ).str();

            std::cout<<"RECORD MESH: Exporter Filename= "<< fileout <<"\n";

            std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
            if ( ! outrecord.good() )
                throw std::runtime_error("Cannot write to file: " + fileout);

            exporter.writeRecord(outrecord,"mesh");
            outrecord.close();
        }


        fileout = (boost::format( "%1%/matlab/field_%2%_%3%.bin" )
                   % Environment::nextsimDir().string()
                   % M_comm.size() /*M_rank*/
                   % step ).str();

        std::cout<<"BINARY: Exporter Filename= "<< fileout <<"\n";

        std::fstream outbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
        if ( ! outbin.good() )
            throw std::runtime_error("Cannot write to file: " + fileout);


        //auto rmap_nodes = M_mesh.mapNodes();
        int prv_global_num_nodes = M_mesh.numGlobalNodes();

        auto interp_in_nodes_nrd = vt_root;

        for (int i=0; i<prv_global_num_nodes; ++i)
        {
            //int ri =  rmap_nodes.left.find(i+1)->second-1;
            int ri =  M_rmap_nodes[i];

            vt_root[i] = interp_in_nodes_nrd[2*ri];
            vt_root[i+prv_global_num_nodes] = interp_in_nodes_nrd[2*ri+1];

            //if reorder mesh
            //vt_root[i] = interp_in_nodes_nrd[2*i];
            //vt_root[i+prv_global_num_nodes] = interp_in_nodes_nrd[2*i+1];

            // if ((vt_root[i] == 0) && (vt_root[i+prv_global_num_nodes] == 0))
            // {
            //     std::cout<<"vt_root["<< 2*i+j <<"]---["<< 2*ri+j <<"]" << ": global "<< ri <<"\n";
            //     std::cout<<"vt_root["<< i <<"]---["<< 2*ri <<"]" << ": global "<< ri <<"\n";
            //     std::cout<<"vt_root["<< i+prv_global_num_nodes <<"]---["<< 2*ri+1 <<"]" << ": global "<< 2*ri+1 <<"\n";
            // }
        }

        // Add time info and export it for plotting
        std::vector<double> timevec(1,current_time);

        exporter.writeField(outbin, timevec, "Time");
        exporter.writeField(outbin, vt_root, "M_VT");
        exporter.writeField(outbin, conc_root, "Concentration");
        exporter.writeField(outbin, thick_root, "Thickness");
        exporter.writeField(outbin, thick_snow, "Snow");
        exporter.writeField(outbin, stress1, "Stress1");
        exporter.writeField(outbin, stress2, "Stress2");
        exporter.writeField(outbin, stress3, "Stress3");
        exporter.writeField(outbin, damage, "Damage");
        // exporter.writeField(outbin, ridged_thin, "Ridgedthin");
        // exporter.writeField(outbin, ridged_thick, "Ridgedthick");
        // exporter.writeField(outbin, random, "Random");
        // exporter.writeField(outbin, tsurf, "Tsurf");
        // exporter.writeField(outbin, hthin, "Hthin");
        // exporter.writeField(outbin, hsthin, "Hsthin");
        // exporter.writeField(outbin, tsurfthin, "Tsurfthin");



        // exporter.writeField(outbin, M_tsurf, "Tsurf");
        // exporter.writeField(outbin, M_sst, "SST");
        // exporter.writeField(outbin, M_sss, "SSS");

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            exporter.writeField(outbin, M_h_thin, "Thin_ice");
            exporter.writeField(outbin, M_hs_thin, "Snow_thin_ice");
            exporter.writeField(outbin, M_tsurf_thin, "Tsurf_thin_ice");

            // Re-create the diagnostic variable 'concentration of thin ice'
            std::vector<double> conc_thin(M_mesh.numTriangles());
            double const rtanalpha = 1.;//c_thin_max/h_thin_max;
            for ( int i=0; i<M_mesh.numTriangles(); ++i )
            {
                conc_thin[i] = std::min(std::min(M_h_thin[i]/physical::hmin, std::sqrt(2.*M_h_thin[i]*rtanalpha)), 1.-M_conc[i]);
            }
            exporter.writeField(outbin, conc_thin, "Concentration_thin_ice");
        }

        outbin.close();


        fileout = (boost::format( "%1%/matlab/field_%2%_%3%.dat" )
                   % Environment::nextsimDir().string()
                   % M_comm.size()
                   % step ).str();

        std::cout<<"RECORD FIELD: Exporter Filename= "<< fileout <<"\n";

        std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
        if ( ! outrecord.good() )
            throw std::runtime_error("Cannot write to file: " + fileout);

        exporter.writeRecord(outrecord);
        outrecord.close();
    }
}

std::string
FiniteElement::gitRevision()
{
    //std::string command = "git rev-parse HEAD";
    return this->system("git rev-parse HEAD");
}

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
}

std::string
FiniteElement::getEnv(std::string const& envname)
{
    char* senv = ::getenv(envname.c_str());
    return std::string(senv);
}

void
FiniteElement::writeLogFile()
{
    std::string logfilename = "";
    if ((vm["simul.logfile"].as<std::string>()).empty())
    {
        logfilename = "nextsim.log";
    }
    else
    {
        logfilename = vm["simul.logfile"].as<std::string>();
    }

    std::string fileout = (boost::format( "%1%/%2%" )
               % M_export_path
               % logfilename ).str();

    std::fstream logfile(fileout, std::ios::out | std::ios::trunc);
    std::cout << "Writing log file " << fileout << "...\n";

    if (logfile.is_open())
    {
        logfile << "#----------Info\n";
        logfile << std::setw(40) << std::left << "Build date "  << current_time_local() <<"\n";
        logfile << std::setw(40) << std::left << "Git revision "  << gitRevision() <<"\n";

        logfile << "#----------Compilers\n";
        logfile << std::setw(40) << std::left << "C "  << system("which gcc") << " (version "<< system("gcc -dumpversion") << ")" <<"\n";
        logfile << std::setw(40) << std::left << "C++ "  << system("which g++") << " (version "<< system("g++ -dumpversion") << ")" <<"\n";

        logfile << "#----------Environment variables\n";
        logfile << std::setw(40) << std::left << "NEXTSIMDIR "  << getEnv("NEXTSIMDIR") <<"\n";
        logfile << std::setw(40) << std::left << "SIMDATADIR "  << getEnv("SIMDATADIR") <<"\n";
        logfile << std::setw(40) << std::left << "SIMFORECASTDIR "  << getEnv("SIMFORECASTDIR") <<"\n";
        logfile << std::setw(40) << std::left << "PETSC_DIR "  << getEnv("PETSC_DIR") <<"\n";
        logfile << std::setw(40) << std::left << "BOOST_DIR "  << getEnv("BOOST_DIR") <<"\n";
        logfile << std::setw(40) << std::left << "GMSH_DIR "  << getEnv("GMSH_DIR") <<"\n";
        logfile << std::setw(40) << std::left << "NETCDF_DIR "  << getEnv("NETCDF_DIR") <<"\n";
        logfile << std::setw(40) << std::left << "OPENMPI_LIB_DIR "  << getEnv("OPENMPI_LIB_DIR") <<"\n";
        logfile << std::setw(40) << std::left << "OPENMPI_INCLUDE_DIR "  << getEnv("OPENMPI_INCLUDE_DIR") <<"\n";

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

            logfile << std::setw(40) << std::left << it->first;

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
                    uint i = 0;
                    for (std::vector<std::string>::iterator oit=vect.begin(); oit != vect.end(); oit++, ++i)
                    {
                        //logfile << it->first << "[" << i << "]=" << (*oit) <<"\n";
                        if (i > 0)
                            logfile << std::setw(41) << std::right;

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
}

void
FiniteElement::clear()
{
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

        delete bamgopt_previous;
        delete bamggeom_previous;
        delete bamgmesh_previous;

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
}

// Finalise everything
void
FiniteElement::finalise()
{
    // Don't forget to close the iabp file!
    if (M_use_iabp_drifters)
    {
        M_iabp_file.close();
        M_iabp_out.close();
    }

    this->clear();
}

} // Nextsim
