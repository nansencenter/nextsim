/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

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
    M_vector()
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

    this->distributedMeshProcessing(true);
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
    //M_mesh.setOrdering("bamg");

    //chrono.restart();
    LOG(INFO) <<"["<< M_rank <<"] " <<"filename= "<< M_mesh_filename <<"\n";
    M_mesh.readFromFile(M_mesh_filename);
    //M_mesh.readFromFile("par4topazreducedsplit2.msh");
    //std::cout<<"Reading mesh done in "<< chrono.elapsed() <<"s\n";

    if (!start)
    {
        delete bamggeom;
        delete bamgmesh;

        bamggeom = new BamgGeom();
        bamgmesh = new BamgMesh();
    }

    BamgConvertMeshx(
                     bamgmesh,bamggeom,
                     &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                     M_mesh.numNodes(), M_mesh.numTriangles()
                     );

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_ndof = M_mesh.numGlobalNodes();

    M_local_ndof = M_mesh.numLocalNodesWithoutGhost();
    M_local_ndof_ghost = M_mesh.numLocalNodesWithGhost();

    M_local_nelements = M_mesh.numTrianglesWithoutGhost();
    M_num_nodes = M_local_ndof_ghost;

    this->bcMarkedNodes();

    this->createGraph();

    this->gatherSizes();

    LOG(INFO) <<"["<< M_rank << "] NODES   = "<< M_mesh.numGlobalNodes() << " --- "<< M_local_ndof <<"\n";
    LOG(INFO) <<"["<< M_rank << "] ELEMENTS= "<< M_mesh.numGlobalElements() << " --- "<< M_local_nelements <<"\n";

#if 0
    int num_nodes = boost::mpi::all_reduce(M_comm, M_local_ndof, std::plus<int>());
    int num_elements = boost::mpi::all_reduce(M_comm, M_local_nelements, std::plus<int>());


    if(M_mesh.numGlobalNodesFromSarialMesh() != num_nodes)
    {
        throw std::logic_error("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Inconsistant NODAL PARTITIONS");
    }

    std::cout<<"COMPARE "<< M_mesh.numGlobalElementsFromSarialMesh() << " and "<< num_elements <<"\n";

    if(M_mesh.numGlobalElementsFromSarialMesh() != num_elements)
    {
        throw std::logic_error("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@INCONSISTANT ELEMENT PARTITIONS");
    }
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
                }
            }
            else
            {
                // including ghost nodes
                M_neumann_flags.push_back(lindex);
                //std::cout<<"["<< M_comm.rank() <<"] " << "-----------------here  = "<< flags_root[i]  <<"\n";
            }
        }
    }

    std::sort(M_dirichlet_flags.begin(), M_dirichlet_flags.end());
    M_dirichlet_flags.erase(std::unique(M_dirichlet_flags.begin(), M_dirichlet_flags.end() ), M_dirichlet_flags.end());

    std::sort(M_neumann_flags.begin(), M_neumann_flags.end());
    M_neumann_flags.erase(std::unique(M_neumann_flags.begin(), M_neumann_flags.end() ), M_neumann_flags.end());

    LOG(DEBUG) <<"["<< M_comm.rank() <<"] " << "Dirichlet flags= "<< M_dirichlet_flags.size() <<"\n";
    LOG(DEBUG) <<"["<< M_comm.rank() <<"] " << "Neumann flags  = "<< M_neumann_flags.size() <<"\n";


    // if (M_rank == 1)
    // {
    //     for (int i=0; i<flags_root.size(); ++i)
    //     {
    //         std::cout<<"flags_root["<< i <<"]= "<< flags_root[i] <<"\n";
    //     }
    // }
    // std::cout<<"M_local_ndof= "<< M_local_ndof <<"\n";

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

    // if (M_rank == 0)
    // {
    //     for (int i=0; i<M_dirichlet_nodes.size(); ++i)
    //     {
    //         std::cout<<"M_dirichlet_nodes["<< i <<"]= "<< M_dirichlet_nodes[i] <<"\n";
    //     }
    // }
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
            //if ((it->physical==158) || (it->physical==159) || (it->physical==160) || (it->physical==161))
            //if (it->elementary==M_flag_fix)
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

        LOG(DEBUG) <<"MESH: HMIN= "<< h[0] <<"\n";
        LOG(DEBUG) <<"MESH: HMAX= "<< h[1] <<"\n";
        LOG(DEBUG) <<"MESH: RES = "<< this->resolution(M_mesh_root) <<"\n";

        //M_mesh_type = setup::MeshType::FROM_SPLIT;

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
            this->adaptMesh();
            bamgopt->KeepVertices=1;
            LOG(DEBUG) <<"First adaptation done in "<< chrono.elapsed() <<"s\n";

            chrono.restart();
            // Interpolate hminVertices and hmaxVertices onto the current mesh
            this->interpVertices();
        }

        chrono.restart();
        LOG(DEBUG) <<"AdaptMesh starts\n";
        this->adaptMesh();
        LOG(DEBUG) <<"AdaptMesh done in "<< chrono.elapsed() <<"s\n";

        // Add information on the number of partition to mesh filename
        M_mesh_filename = (boost::format( "par%1%%2%" ) % M_comm.size() % M_mesh_filename ).str();
        LOG(DEBUG) <<"["<< M_rank <<"] " <<"filename= "<< M_mesh_filename <<"\n";

        std::cout<<"------------------------------version       = "<< M_mesh_root.version() <<"\n";
        std::cout<<"------------------------------ordering      = "<< M_mesh_root.ordering() <<"\n";
        std::cout<<"------------------------------space         = "<< (int)M_partition_space <<"\n";
        std::cout<<"------------------------------partitioner   = "<< (int)M_partitioner <<"\n";

        // save mesh (only root process)
        chrono.restart();
        if (M_partition_space == mesh::PartitionSpace::MEMORY)
        {
            M_mesh_root.initGModel();
            M_mesh_root.writeToGModel(M_mesh_filename);
        }
        else if (M_partition_space == mesh::PartitionSpace::DISK)
        {
            M_mesh_root.writeTofile(M_mesh_filename);
        }
        //LOG(DEBUG) <<"Saving mesh done in "<< chrono.elapsed() <<"s\n";
        std::cout <<"Writing mesh done in "<< chrono.elapsed() <<"s\n";

        // partition the mesh on root process (rank 0)
        chrono.restart();
        M_mesh_root.partition(M_mesh_filename,M_partitioner,M_partition_space);
        //LOG(DEBUG) <<"Partitioning mesh done in "<< chrono.elapsed() <<"s\n";
        std::cout <<"Partitioning mesh done in "<< chrono.elapsed() <<"s\n";
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

    M_VT.resize(2*M_num_nodes,0.);
    M_VTM.resize(2*M_num_nodes,0.);
    M_VTMM.resize(2*M_num_nodes,0.);

    M_sst.resize(M_num_elements);
    M_sss.resize(M_num_elements);

    M_h_thin.assign(M_num_elements,0.);
    M_hs_thin.assign(M_num_elements,0.);
    M_tsurf_thin.assign(M_num_elements,0.);

    M_h_ridged_thin_ice.resize(M_num_elements,0.);
    M_h_ridged_thick_ice.resize(M_num_elements,0.);

    M_divergence_rate.resize(M_num_elements,0.);
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

    M_comm.barrier();

    M_random_number.resize(M_num_elements);
    auto id_elements = M_mesh.trianglesIdWithGhost();

    for (int i=0; i<M_random_number.size(); ++i)
    {
        M_random_number[i] = M_random_number_root[id_elements[i]-1];
    }

    M_conc.resize(M_num_elements);
    M_thick.resize(M_num_elements);
    M_damage.resize(M_num_elements);
    M_snow_thick.resize(M_num_elements);
    M_tsurf.resize(M_num_elements);

    for (int i=0; i<M_num_elements; ++i)
    {
        if ((M_conc[i] <= 0.) || (M_thick[i] <= 0.) )
        {
            M_conc[i] = 0.;
            M_thick[i] = 0.;
        }
    }

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
    M_vector_reduction.resize(2*M_num_nodes,0.);
    M_valid_conc.resize(M_local_ndof,false);

    M_fcor.resize(M_num_elements);

    M_asr_nodes_dataset.target_size=M_num_nodes;
    M_asr_elements_dataset.target_size=M_num_elements;
    M_topaz_nodes_dataset.target_size=M_num_nodes;
    M_topaz_elements_dataset.target_size=M_num_elements;
    M_ice_topaz_elements_dataset.target_size=M_num_elements;
    M_ice_amsre_elements_dataset.target_size=M_num_elements;
    M_ice_osisaf_elements_dataset.target_size=M_num_elements;
    M_ice_amsr2_elements_dataset.target_size=M_num_elements;
    M_etopo_elements_dataset.target_size=M_num_elements;
    M_ERAi_nodes_dataset.target_size=M_num_nodes;
    M_ERAi_elements_dataset.target_size=M_num_elements;

    M_asr_nodes_dataset.reloaded=false;
    M_asr_elements_dataset.reloaded=false;
    M_topaz_nodes_dataset.reloaded=false;
    M_topaz_elements_dataset.reloaded=false;
    M_ice_topaz_elements_dataset.reloaded=false;
    M_ice_amsre_elements_dataset.reloaded=false;
    M_ice_osisaf_elements_dataset.reloaded=false;
    M_ice_amsr2_elements_dataset.reloaded=false;
    M_etopo_elements_dataset.reloaded=false;
    M_ERAi_nodes_dataset.reloaded=false;
    M_ERAi_elements_dataset.reloaded=false;

    M_asr_nodes_dataset.grid.loaded=false;
    M_asr_elements_dataset.grid.loaded=false;
    M_topaz_nodes_dataset.grid.loaded=false;
    M_topaz_elements_dataset.grid.loaded=false;
    M_ice_topaz_elements_dataset.grid.loaded=false;
    M_ice_amsre_elements_dataset.grid.loaded=false;
    M_ice_osisaf_elements_dataset.grid.loaded=false;
    M_ice_amsr2_elements_dataset.grid.loaded=false;
    M_etopo_elements_dataset.grid.loaded=false;
    M_ERAi_nodes_dataset.grid.loaded=false;
    M_ERAi_elements_dataset.grid.loaded=false;

    M_Cohesion.resize(M_num_elements);
    M_Compressive_strength.resize(M_num_elements);
    M_time_relaxation_damage.resize(M_num_elements,time_relaxation_damage);

    // root
    M_UM_root.assign(2*M_mesh.numGlobalNodes(),0.);
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
    // Definition of the datasets available in the code
    M_asr_nodes_dataset=DataSet("asr_nodes",M_num_nodes);

    M_asr_elements_dataset=DataSet("asr_elements",M_num_elements);

    M_topaz_nodes_dataset=DataSet("topaz_nodes",M_num_nodes);

    M_topaz_elements_dataset=DataSet("topaz_elements",M_num_elements);

    M_ice_topaz_elements_dataset=DataSet("ice_topaz_elements",M_num_elements);

    M_ice_amsre_elements_dataset=DataSet("ice_amsre_elements",M_num_elements);

    M_ice_osisaf_elements_dataset=DataSet("ice_osisaf_elements",M_num_elements);

    M_ice_amsr2_elements_dataset=DataSet("ice_amsr2_elements",M_num_elements);

    M_etopo_elements_dataset=DataSet("etopo_elements",M_num_elements);//M_num_nodes);

    M_ERAi_nodes_dataset=DataSet("ERAi_nodes",M_num_nodes);

    M_ERAi_elements_dataset=DataSet("ERAi_elements",M_num_elements);
}

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
    time_init = dateStr2Num(vm["simul.time_init"].as<std::string>());
    output_time_step =  days_in_sec/vm["simul.output_per_day"].as<int>();
    ptime_step =  days_in_sec/vm["simul.ptime_per_day"].as<int>();
    mooring_output_time_step =  vm["simul.mooring_output_timestep"].as<double>()*days_in_sec;

#if 1
    time_step = vm["simul.timestep"].as<double>();
    duration = (vm["simul.duration"].as<double>())*days_in_sec;
    restart_time_step =  vm["setup.restart_time_step"].as<double>()*days_in_sec;
    M_use_restart   = vm["setup.use_restart"].as<bool>();
    M_write_restart = vm["setup.write_restart"].as<bool>();
    if ( fmod(restart_time_step,time_step) != 0)
    {
        std::cout << restart_time_step << " " << time_step << "\n";
        throw std::logic_error("restart_time_step not an integer multiple of time_step");
    }

    divergence_min = /*(1./days_in_sec)**/vm["simul.divergence_min"].as<double>();
    compression_factor = vm["simul.compression_factor"].as<double>();
    exponent_compression_factor = vm["simul.exponent_compression_factor"].as<double>();
    ocean_turning_angle_rad = (PI/180.)*vm["simul.oceanic_turning_angle"].as<double>();
    ridging_exponent = vm["simul.ridging_exponent"].as<double>();

    quad_drag_coef_air = vm["simul.ASR_quad_drag_coef_air"].as<double>();
    quad_drag_coef_water = vm["simul.quad_drag_coef_water"].as<double>();

    basal_k2 = vm["simul.Lemieux_basal_k2"].as<double>();
    basal_drag_coef_air = vm["simul.Lemieux_drag_coef_air"].as<double>();
    basal_u_0 = vm["simul.Lemieux_basal_u_0"].as<double>();
    basal_Cb = vm["simul.Lemieux_basal_Cb"].as<double>();

    time_relaxation_damage = vm["simul.time_relaxation_damage"].as<double>()*days_in_sec;
    deltaT_relaxation_damage = vm["simul.deltaT_relaxation_damage"].as<double>()*days_in_sec;

    h_thin_max = vm["simul.h_thin_max"].as<double>();
    c_thin_max = vm["simul.c_thin_max"].as<double>();

    compr_strength = vm["simul.compr_strength"].as<double>();
    tract_coef = vm["simul.tract_coef"].as<double>();
    scale_coef = vm["simul.scale_coef"].as<double>();
    alea_factor = vm["simul.alea_factor"].as<double>();
    cfix = vm["simul.cfix"].as<double>();

    C_fix    = cfix*scale_coef;          // C_fix;...  : cohesion (mohr-coulomb) in MPa (40000 Pa)
    C_alea   = alea_factor*C_fix;        // C_alea;... : alea sur la cohesion (Pa)
    tan_phi = vm["simul.tan_phi"].as<double>();
    ridge_h = vm["simul.ridge_h"].as<double>();

    if ( vm["simul.newice_type"].as<int>() == 4 )
        M_ice_cat_type = setup::IceCategoryType::THIN_ICE;
    else
        M_ice_cat_type = setup::IceCategoryType::CLASSIC;

    const boost::unordered_map<const std::string, setup::AtmosphereType> str2atmosphere = boost::assign::map_list_of
        ("constant", setup::AtmosphereType::CONSTANT)
        ("asr", setup::AtmosphereType::ASR)
        ("erai", setup::AtmosphereType::ERAi);
    M_atmosphere_type = str2atmosphere.find(vm["setup.atmosphere-type"].as<std::string>())->second;

    //std::cout<<"AtmosphereType= "<< (int)M_atmosphere_type <<"\n";

    const boost::unordered_map<const std::string, setup::OceanType> str2ocean = boost::assign::map_list_of
        ("constant", setup::OceanType::CONSTANT)
        ("topaz", setup::OceanType::TOPAZR);
    M_ocean_type = str2ocean.find(vm["setup.ocean-type"].as<std::string>())->second;

    //std::cout<<"OCEANTYPE= "<< (int)M_ocean_type <<"\n";

    const boost::unordered_map<const std::string, setup::IceType> str2conc = boost::assign::map_list_of
        ("constant", setup::IceType::CONSTANT)
        ("target", setup::IceType::TARGET)
        ("topaz", setup::IceType::TOPAZ4)
        ("amsre", setup::IceType::AMSRE)
        ("amsr2", setup::IceType::AMSR2)
        ("osisaf", setup::IceType::OSISAF);
    M_ice_type = str2conc.find(vm["setup.ice-type"].as<std::string>())->second;

    const boost::unordered_map<const std::string, setup::BathymetryType> str2bathymetry = boost::assign::map_list_of
        ("constant", setup::BathymetryType::CONSTANT)
        ("etopo", setup::BathymetryType::ETOPO);
    M_bathymetry_type = str2bathymetry.find(vm["setup.bathymetry-type"].as<std::string>())->second;

    const boost::unordered_map<const std::string, setup::DrifterType> str2drifter = boost::assign::map_list_of
        ("none", setup::DrifterType::NONE)
        ("equallyspaced", setup::DrifterType::EQUALLYSPACED)
        ("iabp", setup::DrifterType::IABP);
    M_drifter_type = str2drifter.find(vm["setup.drifter-type"].as<std::string>())->second;

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
#endif
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
FiniteElement::minAngle(FEMeshType const& mesh, std::vector<double> const& um, double factor) const
{
    auto movedmesh = mesh;
    movedmesh.move(um,factor);

    std::vector<double> all_min_angle(movedmesh.numTriangles());

    int cpt = 0;
    for (auto it=movedmesh.triangles().begin(), end=movedmesh.triangles().end(); it!=end; ++it)
    {
        all_min_angle[cpt] = this->minAngles(*it,movedmesh);
        ++cpt;
    }

    double min_angle = *std::min_element(all_min_angle.begin(),all_min_angle.end());
    //return min_angle;
    return boost::mpi::all_reduce(M_comm, min_angle, boost::mpi::minimum<double>());
}

template<typename FEMeshType>
bool
FiniteElement::flip(FEMeshType const& mesh, std::vector<double> const& um, double factor) const
{
    auto movedmesh = mesh;
    movedmesh.move(um,factor);

    std::vector<double> area(movedmesh.numTriangles());

    int cpt = 0;
    for (auto it=movedmesh.triangles().begin(), end=movedmesh.triangles().end(); it!=end; ++it)
    {
        area[cpt] = this->jacobian(*it,movedmesh);
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

    M_comm.barrier();
}

void
FiniteElement::gatherFieldsElement(std::vector<double>& interp_in_elements)
{
    //M_comm.barrier();

    timer["gather"].first.restart();

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------GATHER ELEMENT starts\n";

    int prv_num_nodes = M_local_ndof;
    int prv_num_elements = M_local_nelements;

    // ELEMENT INTERPOLATION With Cavities
    std::vector<int> sizes_elements = M_sizes_elements;
    M_nb_var_element=17;
    std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = M_nb_var_element*f; });

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

        // divergence_rate
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_divergence_rate[i];
        tmp_nb_var++;

        // h_ridged_thin_ice
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_h_ridged_thin_ice[i];
        tmp_nb_var++;

        // h_ridged_thick_ice
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_h_ridged_thick_ice[i];
        tmp_nb_var++;

        // random_number
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_random_number[i];
        tmp_nb_var++;

        // Ice surface temperature
        interp_elt_in_local[M_nb_var_element*i+tmp_nb_var] = M_tsurf[i];
        tmp_nb_var++;

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

    //LOG(DEBUG) <<"ELEMENT: Interp done in "<< chrono.elapsed() <<"\n";

    //timer["gather2"].first.elapsed();
    //LOG(DEBUG) <<"["<< M_rank <<"]: " <<"GATHERV2 starts\n";

    if (M_rank == 0)
    {
        interp_in_elements.resize(M_nb_var_element*M_mesh_previous_root.numTriangles());
        boost::mpi::gatherv(M_comm, interp_elt_in_local, &interp_in_elements[0], sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, interp_elt_in_local, 0);
    }

    //LOG(DEBUG) <<"["<< M_rank <<"]: " <<"GATHERV2 done in "<< timer["gather2"].first.elapsed() <<"s\n";


    if (M_rank == 0)
    {
        auto rmap_elements = M_mesh.mapElements();
        auto interp_in_elements_nrd = interp_in_elements;

        for (int i=0; i<M_mesh_previous_root.numTriangles(); ++i)
        {
            int ri = rmap_elements.left.find(i+1)->second-1;

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
        auto rmap_elements = M_mesh.mapElements();
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

    M_divergence_rate.assign(M_num_elements,0.);
    M_h_ridged_thin_ice.assign(M_num_elements,0.);
    M_h_ridged_thick_ice.assign(M_num_elements,0.);
    M_random_number.resize(M_num_elements);
    M_tsurf.assign(M_num_elements,0.);
    M_h_thin.assign(M_num_elements,0.);
    M_hs_thin.assign(M_num_elements,0.);
    M_tsurf_thin.assign(M_num_elements,0.);

    M_sst.resize(M_num_elements);
    M_sss.resize(M_num_elements);


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

        // divergence_rate
        M_divergence_rate[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // // h_ridged_thin_ice
        M_h_ridged_thin_ice[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // h_ridged_thick_ice
        M_h_ridged_thick_ice[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // random_number
        M_random_number[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        //M_random_number[i] = std::max(0., std::min(1.,out_elt_values[M_nb_var_element*i+tmp_nb_var]));
        tmp_nb_var++;

        // Ice surface temperature
        M_tsurf[i] = out_elt_values[M_nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

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

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------SCATTER ELEMENT done in "<< timer["scatter"].first.elapsed() <<"s\n";
}

void
FiniteElement::interpFieldsElement()
{
    M_comm.barrier();

    std::vector<double> interp_in_elements;

    this->gatherFieldsElement(interp_in_elements);

    M_comm.barrier();

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

#if 1
        // The interpolation with the cavities still needs to be tested on a long run.
        // By default, we then use the non-conservative MeshToMesh interpolation

        //chrono.restart();
        //std::cout<<"InterpFromMeshToMesh2dCavities starts\n";

        InterpFromMeshToMesh2dCavities(&interp_elt_out,&interp_in_elements[0],M_nb_var_element,
                                       &surface_previous[0], &surface[0], bamgmesh_previous, bamgmesh_root);

        //std::cout<<"InterpFromMeshToMesh2dCavities done in "<< chrono.elapsed() <<"\n";
#endif

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

    M_comm.barrier();

    this->distributedMeshProcessing();

    M_comm.barrier();

    this->scatterFieldsElement(interp_elt_out);

    if (M_rank == 0)
    {
        xDelete<double>(interp_elt_out);
    }

    // std::cout<<"["<< M_rank <<"]: " <<"AFTER Min val= "<< *std::min_element(M_conc.begin(), M_conc.end()) <<"\n";
    // std::cout<<"["<< M_rank <<"]: " <<"AFTER Max val= "<< *std::max_element(M_conc.begin(), M_conc.end()) <<"\n";
}

void
FiniteElement::gatherFieldsNode(std::vector<double>& interp_in_nodes, bimap_type const& rmap_nodes, std::vector<int> sizes_nodes)
{
    timer["gather.node"].first.restart();

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------GATHER NODE starts\n";

    M_nb_var_node = 8;
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
            int ri =  rmap_nodes.left.find(i+1)->second-1;

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
        //auto rmap_nodes = M_mesh.mapNodes();
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
    }


    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------SCATTER NODE done in "<< timer["scatter.node"].first.elapsed() <<"s\n";
}

void
FiniteElement::interpFieldsNode(bimap_type const& rmap_nodes, std::vector<int> sizes_nodes)
{

    M_comm.barrier();

    std::vector<double> interp_in_nodes;

    this->gatherFieldsNode(interp_in_nodes, rmap_nodes, sizes_nodes);

    M_comm.barrier();

    double* interp_nd_out;

    if (M_rank == 0)
    {
        InterpFromMeshToMesh2dx(&interp_nd_out,
                                &M_mesh_previous_root.indexTr()[0],&M_mesh_previous_root.coordX()[0],&M_mesh_previous_root.coordY()[0],
                                //M_prv_global_num_nodes,M_prv_global_num_elements,
                                M_mesh_previous_root.numNodes(),M_mesh_previous_root.numTriangles(),
                                &interp_in_nodes[0],
                                //M_prv_global_num_nodes,M_nb_var_node,
                                M_mesh_previous_root.numNodes(),M_nb_var_node,
                                &M_mesh_root.coordX()[0],&M_mesh_root.coordY()[0],M_mesh_root.numNodes(),
                                false);
    }


    this->scatterFieldsNode(interp_nd_out);

    if (M_rank == 0)
    {
        xDelete<double>(interp_nd_out);
    }

    M_comm.barrier();



#if 0
    M_comm.barrier();

    std::vector<int> indextr(3*M_prv_global_num_elements);
    std::vector<double> coordX(M_prv_global_num_nodes);
    std::vector<double> coordY(M_prv_global_num_nodes);

    if (M_rank == 0)
    {
        indextr = M_mesh_previous_root.indexTr();
        coordX = M_mesh_previous_root.coordX();
        coordY = M_mesh_previous_root.coordY();
    }

    // mesh elements
    boost::mpi::broadcast(M_comm, &indextr[0], 3*M_prv_global_num_elements, 0);

    // mesh x-coordinates
    boost::mpi::broadcast(M_comm, &coordX[0], M_prv_global_num_nodes, 0);

    // mesh y-coordinates
    boost::mpi::broadcast(M_comm, &coordY[0], M_prv_global_num_nodes, 0);

    // if (M_rank == 1)
    // {
    //     for (int i=0; i<10; ++i)
    //     {
    //         std::cout<<"["<< M_rank <<"]: indextr["<< i <<"]= "<< indextr[i] <<"\n";
    //     }
    // }


    M_nb_var_node = 8;
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
    }


    std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = M_nb_var_node*f; });

    std::vector<double> interp_in_nodes(M_nb_var_node*M_prv_global_num_nodes);

    chrono.restart();

    if (M_rank == 0)
    {
        boost::mpi::gatherv(M_comm, interp_node_in_local, &interp_in_nodes[0], sizes_nodes, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, interp_node_in_local, 0);
    }

    boost::mpi::broadcast(M_comm, &interp_in_nodes[0], M_nb_var_node*M_prv_global_num_nodes, 0);

    //std::cout<<"ALL_GATHER done in "<< chrono.elapsed() <<"s\n";

    // if (M_rank == 1)
    // {
    //     for (int i=0; i<10; ++i)
    //     {
    //         std::cout<<"BEFORE out_values["<< i <<"]= "<< interp_in_nodes[i] <<"\n";
    //     }
    // }

    auto interp_in_nodes_nrd = interp_in_nodes;

    for (int i=0; i<M_prv_global_num_nodes; ++i)
    {
        int ri =  rmap_nodes.left.find(i+1)->second-1;

        for (int j=0; j<M_nb_var_node; ++j)
        {
            interp_in_nodes[M_nb_var_node*i+j] = interp_in_nodes_nrd[M_nb_var_node*ri+j];
        }
    }

    // if (M_rank == 1)
    // {
    //     for (int i=0; i<10; ++i)
    //     {
    //         std::cout<<"out_values["<< i <<"]= "<< interp_in_nodes[i] <<"\n";
    //     }
    // }

    double* interp_out;

    chrono.restart();

    InterpFromMeshToMesh2dx(&interp_out,
                            &indextr[0],&coordX[0],&coordY[0],
                            M_prv_global_num_nodes,M_prv_global_num_elements,
                            &interp_in_nodes[0],
                            M_prv_global_num_nodes,M_nb_var_node,
                            &M_mesh.coordX()[0],&M_mesh.coordY()[0],M_mesh.numNodes(),
                            false);

    //std::cout<<"NODAL INTERPOLATION done in "<< chrono.elapsed() <<"s\n";

    M_VT.assign(2*M_num_nodes,0.);
    M_VTM.assign(2*M_num_nodes,0.);
    M_VTMM.assign(2*M_num_nodes,0.);
    M_UM.assign(2*M_num_nodes,0.);


    for (int i=0; i<M_num_nodes; ++i)
    {
        // VT
        M_VT[i] = interp_out[M_nb_var_node*i];
        M_VT[i+M_num_nodes] = interp_out[M_nb_var_node*i+1];

        // VTM
        M_VTM[i] = interp_out[M_nb_var_node*i+2];
        M_VTM[i+M_num_nodes] = interp_out[M_nb_var_node*i+3];

        // VTMM
        M_VTMM[i] = interp_out[M_nb_var_node*i+4];
        M_VTMM[i+M_num_nodes] = interp_out[M_nb_var_node*i+5];

        // UM
        M_UM[i] = interp_out[M_nb_var_node*i+6];
        M_UM[i+M_num_nodes] = interp_out[M_nb_var_node*i+7];
    }

    xDelete<double>(interp_out);

    M_comm.barrier();
#endif
}

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
        auto rmap_nodes = M_mesh.mapNodes();
        int prv_global_num_nodes = M_mesh.numGlobalNodes();

        auto interp_in_nodes_nrd = um;

        for (int i=0; i<prv_global_num_nodes; ++i)
        {
            int ri =  rmap_nodes.left.find(i+1)->second-1;

            um[i] = interp_in_nodes_nrd[2*ri];
            um[i+prv_global_num_nodes] = interp_in_nodes_nrd[2*ri+1];
        }
    }
}

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

    this->gatherUM(um_root);

    if (M_rank == 0)
    {
        chrono.restart();
        LOG(DEBUG) <<"Flip starts\n";

        while (flip)
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
                LOG(DEBUG) <<"Interp vertices done in "<< timer["interpvertices"].first.elapsed() <<"\n";
            }

            timer["adaptmesh"].first.restart();
            LOG(DEBUG) <<"---TRUE AdaptMesh starts\n";
            this->adaptMesh();
            LOG(DEBUG) <<"---TRUE AdaptMesh done in "<< timer["adaptmesh"].first.elapsed() <<"s\n";

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
            std::cout<<"------------------------------space         = "<< (int)M_partition_space <<"\n";
            std::cout<<"------------------------------partitioner   = "<< (int)M_partitioner <<"\n";


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
            LOG(DEBUG) <<"Saving mesh done in "<< timer["savemesh"].first.elapsed() <<"s\n";

#if 0
            src_fname = Environment::nextsimDir().string() + "/mesh/" + M_mesh_filename;
            desc_fname = Environment::nextsimDir().string() + "/mesh/" + "seq_" + M_mesh_filename;
            fs::copy_file(fs::path(src_fname), fs::path(desc_fname), fs::copy_option::overwrite_if_exists);
#endif

            // partition the mesh on root process (rank 0)
            timer["meshpartition"].first.restart();
            LOG(DEBUG) <<"Partitioning mesh starts\n";
            M_mesh_root.partition(M_mesh_filename,M_partitioner,M_partition_space);
            LOG(DEBUG) <<"Partitioning mesh done in "<< timer["meshpartition"].first.elapsed() <<"s\n";
        }
    } // rank 0

    // --------------------------------BEGINNING-------------------------

    M_prv_local_ndof = M_local_ndof;
    M_prv_num_nodes = M_num_nodes;
    M_prv_num_elements = M_local_nelements;
    bimap_type prv_rmap_nodes = M_mesh.mapNodes();
    M_prv_global_num_nodes = M_mesh.numGlobalNodes();
    M_prv_global_num_elements = M_mesh.numGlobalElements();
    std::vector<int> sizes_nodes = M_sizes_nodes;

    M_comm.barrier();

    timer["felt"].first.restart();
    this->interpFieldsElement();
    LOG(DEBUG) <<"interpFieldsElement done in "<< timer["felt"].first.elapsed() <<"s\n";

    timer["fnd"].first.restart();
    this->interpFieldsNode(prv_rmap_nodes, sizes_nodes);
    LOG(DEBUG) <<"interpFieldsNode done in "<< timer["fnd"].first.elapsed() <<"s\n";

    // --------------------------------END-------------------------------

    LOG(DEBUG) <<"TIMER REGRIDDING: this is done in "<< timer["regrid"].first.elapsed() <<"s\n";

    M_comm.barrier();

    this->assignVariables();
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

    Bamgx(bamgmesh_root,bamggeom_root,bamgmesh_previous,bamggeom_previous,bamgopt_previous);

    this->importBamg(bamgmesh_root);

    LOG(DEBUG) <<"CLOSED: FLAGS SIZE BEFORE= "<< M_dirichlet_flags_root.size() <<"\n";
    LOG(DEBUG) <<"OPEN  : FLAGS SIZE BEFORE= "<< M_neumann_flags_root.size() <<"\n";

    // update dirichlet nodes
    M_dirichlet_flags_root.resize(0);
    M_neumann_flags_root.resize(0);

    for (int edg=0; edg<bamgmesh_root->EdgesSize[0]; ++edg)
    {
        if (bamgmesh_root->Edges[3*edg+2] == M_flag_fix)
        {
            M_dirichlet_flags_root.push_back(bamgmesh_root->Edges[3*edg]/*-1*/);
        }
        else
        {
            M_neumann_flags_root.push_back(bamgmesh_root->Edges[3*edg]/*-1*/);
        }
    }

    int num_nodes = M_mesh_root.numNodes();

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

    M_matrix->zero();
    M_vector->zero();

    double coef_V, coef_Voce, coef_Vair, coef_basal, coef_X, coef_Y, coef_C;
    double coef = 0;
    double coef_P = 0.;
    double mass_e = 0.;
    double surface_e = 0.;
    double g_ssh_e_x = 0.;
    double g_ssh_e = 0.;
    double g_ssh_e_y = 0.;
    double tmp_thick, tmp_conc;

    double welt_oce_ice = 0.;
    double welt_air_ice = 0.;
    double welt_ice = 0.;
    double welt_ssh = 0.;
    double element_ssh = 0.;
    double critical_h = 0.;

    double norm_Voce_ice = 0.;
    double norm_Vair_ice = 0.;
    double norm_Vice = 0.;

    double Vcor_index_v, Vcor_index_u;

    double b0tj_sigma_hu = 0.;
    double b0tj_sigma_hv = 0.;
    double mloc = 0.;

    double duu, dvu, duv, dvv;
    int index_u, index_v;

    int nind;

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

    double cos_ocean_turning_angle=std::cos(ocean_turning_angle_rad);
    double sin_ocean_turning_angle=std::sin(ocean_turning_angle_rad);

    // ---------- Assembling starts -----------

    // if (M_rank == 0)
    //     std::cout<<"Assembling starts\n";

    //chrono.restart();
    timer["assembly"].first.restart();

    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        /* Compute the value that only depends on the element */
        welt_oce_ice = 0.;
        welt_air_ice = 0.;
        welt_ice = 0.;
        welt_ssh = 0.;
        nind;

        for (int i=0; i<3; ++i)
        {
            nind = it->indices[i]-1;
            welt_oce_ice += std::hypot(M_VT[nind]-M_ocean[nind],M_VT[nind+M_num_nodes]-M_ocean[nind+M_num_nodes]);
            welt_air_ice += std::hypot(M_VT[nind]-M_wind [nind],M_VT[nind+M_num_nodes]-M_wind [nind+M_num_nodes]);
            welt_ice += std::hypot(M_VT[nind],M_VT[nind+M_num_nodes]);

            welt_ssh += M_ssh[nind];
        }

        norm_Voce_ice = welt_oce_ice/3.;
        norm_Vair_ice = welt_air_ice/3.;
        norm_Vice = welt_ice/3.;

        element_ssh = welt_ssh/3.;

        coef_Vair = (vm["simul.lin_drag_coef_air"].as<double>()+(quad_drag_coef_air*norm_Vair_ice));
        coef_Vair *= (vm["simul.rho_air"].as<double>());

        coef_Voce = (vm["simul.lin_drag_coef_water"].as<double>()+(quad_drag_coef_water*norm_Voce_ice));
        coef_Voce *= (vm["simul.rho_water"].as<double>());


        critical_h = M_conc[cpt]*(M_element_depth[cpt]+element_ssh)/(vm["simul.Lemieux_basal_k1"].as<double>());
        double _coef = std::max(0., M_thick[cpt]-critical_h);
        coef_basal = quad_drag_coef_air*basal_k2/(basal_drag_coef_air*(norm_Vice+basal_u_0));
        coef_basal *= _coef*std::exp(-basal_Cb*(1.-M_conc[cpt]));

        //std::vector<double> sigma_P(3,0.); /* temporary variable for the resistance to the compression */
        //std::vector<double> B0Tj_sigma_h(2,0);

        tmp_thick=(0.05>M_thick[cpt]) ? 0.05 : M_thick[cpt];
        tmp_conc=(0.01>M_conc[cpt]) ? 0.01 : M_conc[cpt];

        coef = young*(1.-M_damage[cpt])*tmp_thick*std::exp(ridging_exponent*(1.-tmp_conc));

        coef_P = 0.;
        if(M_divergence_rate[cpt] < 0.)
        {
            coef_P = compression_factor*std::pow(tmp_thick,exponent_compression_factor)*std::exp(ridging_exponent*(1.-tmp_conc));
            coef_P = coef_P/(std::abs(M_divergence_rate[cpt])+divergence_min);
        }

        mass_e = rhoi*tmp_thick + rhos*M_snow_thick[cpt];
        mass_e = (tmp_conc > 0.) ? (mass_e/tmp_conc):0.;
        surface_e = M_surface[cpt];

        // /* compute the x and y derivative of g*ssh */
        g_ssh_e_x = 0.;
        g_ssh_e_y = 0.;
        g_ssh_e;
        for(int i=0; i<3; i++)
        {
            g_ssh_e = (vm["simul.gravity"].as<double>())*M_ssh[it->indices[i]-1] /*g_ssh*/;   /* g*ssh at the node k of the element e */
            g_ssh_e_x += M_shape_coeff[cpt][i]*g_ssh_e; /* x derivative of g*ssh */
            g_ssh_e_y += M_shape_coeff[cpt][i+3]*g_ssh_e; /* y derivative of g*ssh */
        }

        coef_C     = mass_e*M_fcor[cpt];              /* for the Coriolis term */
        coef_V     = mass_e/time_step;             /* for the inertial term */
        coef_X     = - mass_e*g_ssh_e_x;              /* for the ocean slope */
        coef_Y     = - mass_e*g_ssh_e_y;              /* for the ocean slope */

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

        //std::vector<int> rindices_u;
        std::vector<int> rindices(6); //new
        std::vector<int> cindices(6);
        int vs = 0;

        for (int s=0; s<3; ++s)
        {
            int index_u = it->indices[s]-1;

            if (!it->ghostNodes[s])
            {
                // rindices_u.push_back(index_u);
                // rindices.push_back(index_u);
                // rindices.push_back(index_u+M_local_ndof);

                rindices[2*vs] = index_u;
                rindices[2*vs+1] = index_u+M_local_ndof;
                ++vs;

                if((M_conc[cpt]>0.))
                {
                    M_valid_conc[index_u] = true;
                }
            }

            //cindices.push_back(index_u);
            //cindices.push_back(index_u+M_local_ndof_ghost);
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

            Vcor_index_v=beta0*M_VT[index_v] + beta1*M_VTM[index_v] + beta2*M_VTMM[index_v];
            Vcor_index_u=beta0*M_VT[index_u] + beta1*M_VTM[index_u] + beta2*M_VTMM[index_u];

            if (!it->ghostNodes[j])
            {
                l_j = l_j + 1;

                for(int i=0; i<3; i++)
                {
                    /* Row corresponding to indice i (we also assemble terms in row+1) */
                    //row = (mwIndex)it[2*i]-1; /* -1 to use the indice convention of C */

                    /* Select the nodal weight values from M_loc */
                    mloc = M_Mass[3*j+i];

                    b0tj_sigma_hu = 0.;
                    b0tj_sigma_hv = 0.;

                    for(int k=0; k<3; k++)
                    {
                        b0tj_sigma_hu += M_B0T[cpt][k*6+2*i]*(M_sigma[3*cpt+k]*tmp_thick/*+sigma_P[k]*/);
                        b0tj_sigma_hv += M_B0T[cpt][k*6+2*i+1]*(M_sigma[3*cpt+k]*tmp_thick/*+sigma_P[k]*/);
                    }

                    /* ---------- UU component */
                    duu = surface_e*( mloc*(coef_Vair+coef_Voce*cos_ocean_turning_angle+coef_V+coef_basal)
                                      +M_B0T_Dunit_B0T[cpt][(2*i)*6+2*j]*coef*time_step+M_B0T_Dunit_comp_B0T[cpt][(2*i)*6+2*j]*coef_P);

                    /* ---------- VU component */
                    dvu = surface_e*(+M_B0T_Dunit_B0T[cpt][(2*i+1)*6+2*j]*coef*time_step+M_B0T_Dunit_comp_B0T[cpt][(2*i+1)*6+2*j]*coef_P);

                    /* ---------- UV component */
                    duv = surface_e*(+M_B0T_Dunit_B0T[cpt][(2*i)*6+2*j+1]*coef*time_step+M_B0T_Dunit_comp_B0T[cpt][(2*i)*6+2*j+1]*coef_P);

                    /* ---------- VV component */
                    dvv = surface_e*( mloc*(coef_Vair+coef_Voce*cos_ocean_turning_angle+coef_V+coef_basal)
                                      +M_B0T_Dunit_B0T[cpt][(2*i+1)*6+2*j+1]*coef*time_step+M_B0T_Dunit_comp_B0T[cpt][(2*i+1)*6+2*j+1]*coef_P);


                    data[(2*l_j  )*6+2*i  ] = duu;
                    data[(2*l_j+1)*6+2*i  ] = dvu;
                    data[(2*l_j  )*6+2*i+1] = duv;
                    data[(2*l_j+1)*6+2*i+1] = dvv;


                    fvdata[2*i] += surface_e*( mloc*( coef_Vair*M_wind[index_u]
                                                      +coef_Voce*cos_ocean_turning_angle*M_ocean[index_u]
                                                      +coef_X
                                                      +coef_V*M_VT[index_u]
                                                      -coef_Voce*sin_ocean_turning_angle*(M_ocean[index_v]-M_VT[index_v])
                                                      +coef_C*Vcor_index_v)
                                               - b0tj_sigma_hu/3);


                    fvdata[2*i+1] += surface_e*( mloc*( coef_Vair*M_wind[index_v]
                                                        +coef_Voce*cos_ocean_turning_angle*M_ocean[index_v]
                                                        +coef_Y
                                                        +coef_V*M_VT[index_v]
                                                        +coef_Voce*sin_ocean_turning_angle*(M_ocean[index_u]-M_VT[index_u])
                                                        -coef_C*Vcor_index_u)
                                                 - b0tj_sigma_hv/3);

                }
            }
        }

        // update matrix
        M_matrix->addMatrix(&rindices[0], rindices.size(),
                            &cindices[0], cindices.size(), &data[0]);

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

#if 0
        if((M_conc[cpt]>0.))
        {
            for (int const& idn : rindices_u)
                M_valid_conc[idn] = true;
        }
#endif

        ++cpt;
    }

    // close petsc matrix
    M_matrix->close();

    // close petsc vector
    M_vector->close();

    if (M_rank == 0)
    {
        //std::cout<<"Assembling done\n";
        //std::cout<<"TIMER ASSEMBLY= " << chrono.elapsed() <<"s\n";
        std::cout<<"TIMER ASSEMBLY= " << timer["assembly"].first.elapsed() <<"s\n";
    }
    //std::cout<<"[" << M_rank <<"] " <<" TIMER ASSEMBLY= " << chrono.elapsed() <<"s\n";

    // extended dirichlet nodes (add nodes where M_conc <= 0)
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

    //std::cout<<"["<< M_rank <<"] EXTENDED = "<< extended_dirichlet_nodes.size() <<"\n";
    //std::cout<<"["<< M_rank <<"] DIRICHLET= "<< M_dirichlet_nodes.size() <<"\n";

    //chrono.restart();
    //M_matrix->on(M_dirichlet_nodes,*M_vector);
    M_matrix->on(extended_dirichlet_nodes,*M_vector);

    //std::cout<<"[" << M_rank <<"] " <<"-------------------DIFF SIZE EXTENDED_DIRICHLET= " << (int)extended_dirichlet_nodes.size()-(int)M_dirichlet_nodes.size() <<"\n";

    // if (M_rank==0)
    //     std::cout <<"TIMER DBCA= " << chrono.elapsed() <<"s\n";


    //std::cout<<"[" << M_rank <<"] " <<"TIMER DBCA= " << chrono.elapsed() <<"s\n";

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
    M_Dunit_comp.assign(9,0);
    M_Mass.assign(9,0);

    for (int k=0; k<6; k+=3)
    {
        for (int kk=0; kk<2; ++kk )
        {
            M_Dunit[k+kk] = (1-((k+kk)%2)*(1-nu0))/(1-std::pow(nu0,2.));
            M_Dunit_comp[k+kk] = 1.;
        }
    }
    M_Dunit[8] = (1-nu0)/(2.*(1-std::pow(nu0,2.)));

    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            M_Mass[3*i+j] = ((i == j) ? 2.0 : 1.0)/12.0;
            //std::cout<< std::left << std::setw(12) << Mass[3*i+j] <<"  ";
        }
    }

    M_B0T.resize(M_num_elements);
    M_B0T_Dunit_B0T.resize(M_num_elements);
    M_B0T_Dunit_comp_B0T.resize(M_num_elements);
    M_shape_coeff.resize(M_num_elements);

    std::vector<double> B0T(18,0);
    std::vector<double> B0Tj_Dunit(6,0);
    std::vector<double> B0T_Dunit_B0T(36,0);
    std::vector<double> B0Tj_Dunit_comp(6,0);
    std::vector<double> B0T_Dunit_comp_B0T(36,0);

    double B0Tj_Dunit_tmp0, B0Tj_Dunit_tmp1;
    double B0Tj_Dunit_B0Ti_tmp0, B0Tj_Dunit_B0Ti_tmp1, B0Tj_Dunit_B0Ti_tmp2, B0Tj_Dunit_B0Ti_tmp3;
    double B0Tj_Dunit_comp_tmp0, B0Tj_Dunit_comp_tmp1;
    double B0Tj_Dunit_comp_B0Ti_tmp0, B0Tj_Dunit_comp_B0Ti_tmp1, B0Tj_Dunit_comp_B0Ti_tmp2, B0Tj_Dunit_comp_B0Ti_tmp3;

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

                B0Tj_Dunit_comp_tmp0 = 0.;
                B0Tj_Dunit_comp_tmp1 = 0.;

                for(int kk=0; kk<3; kk++)
                {
                    B0Tj_Dunit_tmp0 += B0T[kk*6+2*j]*M_Dunit[3*i+kk];
                    B0Tj_Dunit_tmp1 += B0T[kk*6+2*j+1]*M_Dunit[3*i+kk];

                    B0Tj_Dunit_comp_tmp0 += B0T[kk*6+2*j]*M_Dunit_comp[3*i+kk];
                    B0Tj_Dunit_comp_tmp1 += B0T[kk*6+2*j+1]*M_Dunit_comp[3*i+kk];
                }

                B0Tj_Dunit[2*i] = B0Tj_Dunit_tmp0;
                B0Tj_Dunit[2*i+1] = B0Tj_Dunit_tmp1;

                B0Tj_Dunit_comp[2*i] = B0Tj_Dunit_comp_tmp0;
                B0Tj_Dunit_comp[2*i+1] = B0Tj_Dunit_comp_tmp1;
            }

            for(int i=0; i<3; i++)
            {
                /* The rigidity matrix */
                /* scalar product of B0Ti_Dunit and the first column of B0T */
                B0Tj_Dunit_B0Ti_tmp0 = 0.;
                B0Tj_Dunit_B0Ti_tmp1 = 0.;
                B0Tj_Dunit_B0Ti_tmp2 = 0.;
                B0Tj_Dunit_B0Ti_tmp3 = 0.;

                /* scalar product of B0Ti_Dunit_comp and the first column of B0T */
                B0Tj_Dunit_comp_B0Ti_tmp0 = 0.;
                B0Tj_Dunit_comp_B0Ti_tmp1 = 0.;
                B0Tj_Dunit_comp_B0Ti_tmp2 = 0.;
                B0Tj_Dunit_comp_B0Ti_tmp3 = 0.;

                for(int kk=0; kk<3; kk++)
                {
                    B0Tj_Dunit_B0Ti_tmp0 += B0Tj_Dunit[2*kk]*B0T[kk*6+2*i];
                    B0Tj_Dunit_B0Ti_tmp1 += B0Tj_Dunit[2*kk]*B0T[kk*6+2*i+1];
                    B0Tj_Dunit_B0Ti_tmp2 += B0Tj_Dunit[2*kk+1]*B0T[kk*6+2*i];
                    B0Tj_Dunit_B0Ti_tmp3 += B0Tj_Dunit[2*kk+1]*B0T[kk*6+2*i+1];

                    B0Tj_Dunit_comp_B0Ti_tmp0 += B0Tj_Dunit_comp[2*kk]*B0T[kk*6+2*i];
                    B0Tj_Dunit_comp_B0Ti_tmp1 += B0Tj_Dunit_comp[2*kk]*B0T[kk*6+2*i+1];
                    B0Tj_Dunit_comp_B0Ti_tmp2 += B0Tj_Dunit_comp[2*kk+1]*B0T[kk*6+2*i];
                    B0Tj_Dunit_comp_B0Ti_tmp3 += B0Tj_Dunit_comp[2*kk+1]*B0T[kk*6+2*i+1];
                }

                B0T_Dunit_B0T[(2*i)*6+2*j] = B0Tj_Dunit_B0Ti_tmp0;
                B0T_Dunit_B0T[(2*i+1)*6+2*j] = B0Tj_Dunit_B0Ti_tmp1;
                B0T_Dunit_B0T[(2*i)*6+2*j+1] = B0Tj_Dunit_B0Ti_tmp2;
                B0T_Dunit_B0T[(2*i+1)*6+2*j+1] = B0Tj_Dunit_B0Ti_tmp3;

                B0T_Dunit_comp_B0T[(2*i)*6+2*j] = B0Tj_Dunit_comp_B0Ti_tmp0;
                B0T_Dunit_comp_B0T[(2*i+1)*6+2*j] = B0Tj_Dunit_comp_B0Ti_tmp1;
                B0T_Dunit_comp_B0T[(2*i)*6+2*j+1] = B0Tj_Dunit_comp_B0Ti_tmp2;
                B0T_Dunit_comp_B0T[(2*i+1)*6+2*j+1] = B0Tj_Dunit_comp_B0Ti_tmp3;
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
                B0T_Dunit_comp_B0T[i*6+j] = B0T_Dunit_comp_B0T[j*6+i];
            }
        }

        M_shape_coeff[cpt]        = shapecoeff;
        M_B0T[cpt]                = B0T;
        M_B0T_Dunit_B0T[cpt]      = B0T_Dunit_B0T;
        M_B0T_Dunit_comp_B0T[cpt] = B0T_Dunit_comp_B0T;

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
        old_h_ridged_thick_ice = M_h_ridged_thick_ice[cpt];

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            old_h_thin = M_h_thin[cpt];
            old_hs_thin = M_hs_thin[cpt];
            old_h_ridged_thin_ice = M_h_ridged_thin_ice[cpt];
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

        M_divergence_rate[cpt]= (epsilon_veloc[0]+epsilon_veloc[1]);

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
#if 0
        sigma_s=std::hypot((sigma_pred[0]-sigma_pred[1])/2.,sigma_pred[2]);
        sigma_n=           (sigma_pred[0]+sigma_pred[1])/2.;

        /* minimum and maximum normal stress */
        tract_max=tract_coef*M_Cohesion[cpt]/tan_phi;

        /* Correction of the damage */

        if((sigma_n>tract_max) || (sigma_n<(-M_Compressive_strength[cpt])))
        {
            if(sigma_n>tract_max)
            {
                sigma_target=tract_max;
            }
            else
            {
                sigma_target=-M_Compressive_strength[cpt];
            }

            tmp=1.0-sigma_target/sigma_n*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }

        if(sigma_s>M_Cohesion[cpt]-sigma_n*tan_phi)
        {
            tmp=1.0-M_Cohesion[cpt]/(sigma_s+sigma_n*tan_phi)*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }
#endif

#if 1
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
#endif

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
        if( M_divergence_rate[cpt]!=0.)
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
            M_h_ridged_thick_ice[cpt]   =   ridged_thick_ice_volume/surface_new;


            if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                thin_ice_volume = old_h_thin*surface;
                thin_snow_volume = old_hs_thin*surface;
                ridged_thin_ice_volume = old_h_ridged_thin_ice*surface;

                M_h_thin[cpt]        = thin_ice_volume/surface_new;
                M_hs_thin[cpt]   = thin_snow_volume/surface_new;
                M_h_ridged_thin_ice[cpt]    =   ridged_thin_ice_volume/surface_new;
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

                    M_h_ridged_thin_ice[cpt] += ridging_thin_ice;
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

void
FiniteElement::solve()
{
    M_comm.barrier();
    //SolverPetsc ksp;
    //ksp.solve(M_matrix, M_solution, M_vector);

    timer["solution"].first.restart();

    M_solver->solve(_matrix=M_matrix,
                    _solution=M_solution,
                    _rhs=M_vector,
                    _ksp=vm["solver.ksp-type"].as<std::string>()/*"preonly"*/,
                    _pc=vm["solver.pc-type"].as<std::string>()/*"cholesky"*/,
                    _pcfactormatsolverpackage=vm["solver.mat-package-type"].as<std::string>()/*"cholmod"*/,
                    _reuse_prec=false,
                    _rebuild=M_regrid
                    );

    if (M_rank==0)
        std::cout<<"TIMER SOLUTION= " << timer["solution"].first.elapsed() <<"s\n";

    //std::cout<<"[" << M_rank <<"] " <<"TIMER SOLUTION= " << timer["solution"].first.elapsed() <<"s\n";

    M_comm.barrier();

    M_solution->close();
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

    double const ocean_albedo=vm["simul.albedoW"].as<double>();
    double const drag_ocean_t=vm["simul.drag_ocean_t"].as<double>();
    double const drag_ocean_q=vm["simul.drag_ocean_q"].as<double>();

    double const rh0   = 1./vm["simul.hnull"].as<double>();
    double const rPhiF = 1./vm["simul.PhiF"].as<double>();

    double const tanalpha  = h_thin_max/c_thin_max;
    double const rtanalpha = 1./tanalpha;

    double const qi = physical::Lf * physical::rhoi;
    double const qs = physical::Lf * physical::rhos;

    int const newice_type = vm["simul.newice_type"].as<int>();
    int const melt_type = vm["simul.melt_type"].as<int>();
    double const PhiM = vm["simul.PhiM"].as<double>();
    double const PhiF = vm["simul.PhiF"].as<double>();

    const double aw=6.1121e2, bw=18.729, cw=257.87, dw=227.3;
    const double Aw=7.2e-4, Bw=3.20e-6, Cw=5.9e-10;

    const double alpha=0.62197, beta=0.37803;

    // initialisation of the multithreading
    //int thread_id;
    //int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    //#pragma omp parallel for num_threads(max_threads) private(thread_id)
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
        double  Qow=0.;    // Open water heat flux

        // Save old _volumes_ and concentration and calculate wind speed
        double  old_vol=M_thick[i];
        double  old_snow_vol=M_snow_thick[i];
        double  old_conc=M_conc[i];

        double  old_h_thin = 0;
        double  old_hs_thin = 0;
        double  old_conc_thin=0;
        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            old_h_thin  = M_h_thin[i];
            old_hs_thin = M_hs_thin[i];
            old_conc_thin = std::min(std::min(old_h_thin/physical::hmin, std::sqrt(2.*old_h_thin*rtanalpha)), 1-old_conc);
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

		if(M_mld[i]<=0.)
		{
            LOG(DEBUG) << "M_mld[i] = " << M_mld[i] << "\n";
            throw std::logic_error("negative or 0 mld, Oups!!");
		}

        // definition of the fraction of snow
        double tmp_snowfr;
        if(M_snowfr.M_initialized)
            tmp_snowfr=M_snowfr[i];
        else
        {
            if(M_tair[i]<0)
                tmp_snowfr=1.;
            else
                tmp_snowfr=0.;
        }

        // -------------------------------------------------
        // 2) We calculate or set the flux due to nudging
        if ( M_atmosphere_type == setup::AtmosphereType::CONSTANT || M_ocean_type == setup::OceanType::CONSTANT )
        {
            Qdw=Qdw_const;
            Fdw=Fdw_const;
        }
        else
        {
            // nudgeFlux
            if ( M_ocean_salt[i] > physical::si )
            {
                Qdw = -(M_sst[i]-M_ocean_temp[i]) * M_mld[i] * physical::rhow * physical::cpw/timeT;

                double delS = M_sss[i] - M_ocean_salt[i];
                Fdw = delS * M_mld[i] * physical::rhow /(timeS*M_sss[i] - time_step*delS);
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
        if ( M_mixrat[i] < 0. )
        {
            double fa     = 1. + Aw + M_mslp[i]*1e-2*( Bw + Cw*M_dair[i]*M_dair[i] );
            double esta   = fa*aw*std::exp( (bw-M_dair[i]/dw)*M_dair[i]/(M_dair[i]+cw) );
            sphuma = alpha*fa*esta/(M_mslp[i]-beta*fa*esta) ;
        }
        else
            sphuma = M_mixrat[i]/(1.+M_mixrat[i]) ;

        // -------------------------------------------------
        // 3.2) Specific humidity - ocean surface (calcSphumW in matlab)
        double fw     = 1. + Aw + M_mslp[i]*1e-2*( Bw + Cw*M_sst[i]*M_sst[i] );
        double estw   = aw*std::exp( (bw-M_sst[i]/dw)*M_sst[i]/(M_sst[i]+cw) )*(1-5.37e-4*M_sss[i]);
        sphumw = alpha*fw*estw/(M_mslp[i]-beta*fw*estw) ;

        // -------------------------------------------------
        /* Density of air */
        double rhoair = M_mslp[i]/(physical::Ra*(M_tair[i]+tfrwK)) * (1.+sphuma)/(1.+1.609*sphuma);

        /* Sensible heat flux */
        double Qsh = drag_ocean_t*rhoair*physical::cpa*wspeed*( M_sst[i] - M_tair[i] );

        /* Latent heat flux */
        double Lv  = physical::Lv0 - 2.36418e3*M_tair[i] + 1.58927*M_tair[i]*M_tair[i] - 6.14342e-2*std::pow(M_tair[i],3.);
        double Qlh = drag_ocean_q*rhoair*Lv*wspeed*( sphumw - sphuma );

        /* Evaporation */
        evap = Qlh/(physical::rhofw*Lv);

        // Sum them up
        double tmp_Qlw_in;
        if(M_Qlw_in.M_initialized)
            tmp_Qlw_in=M_Qlw_in[i];
        else
        {
        	double tsa = M_tsurf[i] + tfrwK;
        	double taa = M_tair[i]  + tfrwK;
        	tmp_Qlw_in = sigma_sb*pow(taa,4) \
        			*( 1. - 0.261*exp(-7.77e-4*std::pow(taa-tfrwK,2)) ) \
        			*( 1. + 0.275*M_tcc[i] );
        }

        Qow = -M_Qsw_in[i]*(1.-ocean_albedo) - tmp_Qlw_in + Qlw_out + Qsh + Qlh;

        // -------------------------------------------------
        // 4) Thickness change of the ice slab (thermoIce0 in matlab)

        this->thermoIce0(i, wspeed, sphuma, M_conc[i], M_thick[i], M_snow_thick[i], hi, hs, hi_old, Qio, del_hi, M_tsurf[i],tmp_Qlw_in,tmp_snowfr);
        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            this->thermoIce0(i, wspeed, sphuma, old_conc_thin, M_h_thin[i], M_hs_thin[i], hi_thin, hs_thin, hi_thin_old, Qio_thin, del_hi_thin, M_tsurf_thin[i],tmp_Qlw_in,tmp_snowfr);
            M_h_thin[i]  = hi_thin * old_conc_thin;
            M_hs_thin[i] = hs_thin * old_conc_thin;
        }

        // -------------------------------------------------
        // 5) Ice growth over open water and lateral melt (thermoOW in matlab)

        /* Local variables */
        double tw_new, tfrw, newice, del_c, newsnow, h0;

        /* dT/dt due to heatflux atmos.-ocean */
        tw_new = M_sst[i] - Qow*time_step/(M_mld[i]*physical::rhow*physical::cpw);
        tfrw   = physical::mu*M_sss[i];

        /* Form new ice in case of super cooling, and reset Qow and evap */
        if ( tw_new < tfrw )
        {
            newice  = (1.-M_conc[i])*(tfrw-tw_new)*M_mld[i]*physical::rhow*physical::cpw/qi;
            Qow  = -(tfrw-M_sst[i])*M_mld[i]*physical::rhow*physical::cpw/time_step;
            // evap = 0.;
        } else {
            newice  = 0.;
        }

        /* Decide the change in ice fraction (del_c) */
        /* Initialise to be safe */
        del_c = 0.;
        newsnow = 0.;
        if ( newice > 0. )
        {
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
                            del_c = 1.;
                    }
                    break;
                case 3:
                    /* Olason and Harms (09) */
                    h0    = (1.+0.1*wspeed)/15.;
                    del_c = newice/std::max(rPhiF*hi_old,h0);
                    break;
                case 4:
                    /* Thin ice category */
                    thin_ice_redistribute(M_h_thin[i], M_hs_thin[i], newice/(1.-M_conc[i]), M_conc[i],
                                      tanalpha, rtanalpha, h_thin_max, &M_h_thin[i], &newice, &del_c, &newsnow);

                    // Change the snow _thickness_ for thick ice and _volume_ for thin ice
                    M_hs_thin[i]    -= newsnow;
                    // M_snow_thick[i] += newsnow; <- this is done properly below
                    break;
                default:
                    std::cout << "newice_type = " << newice_type << "\n";
                    throw std::logic_error("Wrong newice_type");
            }
            /* Check bounds on del_c */
            del_c = std::min( 1.-M_conc[i], del_c );
        }
        else if ( del_hi < 0. )
        {
            /* Melting conditions */
            switch ( melt_type )
            {
                case 1:
                    /* Hibler (79) using PhiM to tune. PhiM = 0.5 is
                     * equivilent Hibler's (79) approach */
                    if ( M_conc[i] < 1. )
                        del_c = del_hi*M_conc[i]*PhiM/hi_old;
                    else
                        del_c = 0.;
                    break;
                case 2:
                    /* Mellor and Kantha (89) */
                    /* Use the fraction PhiM of (1-c)*Qow to melt laterally */
                    del_c = PhiM*(1.-M_conc[i])*std::min(0.,Qow)*time_step/( hi*qi+hs*qs );
                    /* Deliver the fraction (1-PhiM) of Qow to the ocean */
                    Qow = (1.-PhiM)*Qow;
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
        else
        {
            /* There is no ice growing or melting */
            del_c = 0.;
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
                hs  = ( hs*(old_conc) + newsnow )/M_conc[i];
            }
        }

        /* Check limits */
        if ( M_conc[i] < physical::cmin || hi < physical::hmin )
        {
            //Qow    = Qow + (M_conc[i]-del_c)*hi*qi/time_step + (M_conc[i]-del_c)*hs*qs/time_step;
            //SYL: is this commented line right??
            // I think irt is wrong as we already modify Qow, I would do:
            // Qow    = Qow + M_conc[i]*hi*qi/time_step + M_conc[i]*hs*qs/time_step;
            // Extract heat from the ocean corresponding to the heat in all the
            // ice and snow present at the start of the time step.
            Qow    = Qow + old_conc*hi*qi/time_step + old_conc*hs*qs/time_step;
            M_conc[i]  = 0.;
            M_tsurf[i] = 0.;
            hi     = 0.;
            hs     = 0.;
        }

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
        rain = (1.-old_conc-old_conc_thin)*M_precip[i] + (old_conc+old_conc_thin)*(1.-tmp_snowfr)*M_precip[i];
        emp  = (evap*(1.-old_conc-old_conc_thin)-rain);

        Qio_mean = Qio*old_conc + Qio_thin*old_conc_thin;
        Qow_mean = Qow*(1.-old_conc-old_conc_thin);

        /* Heat-flux */
        M_sst[i] = M_sst[i] - time_step*( Qio_mean + Qow_mean - Qdw )/(physical::rhow*physical::cpw*M_mld[i]);

        /* Change in salinity */
        M_sss[i] = M_sss[i] + ( (M_sss[i]-physical::si)*physical::rhoi*del_vi + M_sss[i]*(del_vs*physical::rhos + (emp-Fdw)*time_step) )
            / ( M_mld[i]*physical::rhow - del_vi*physical::rhoi - ( del_vs*physical::rhos + (emp-Fdw)*time_step) );

        // -------------------------------------------------
        // 8) Damage manipulation (thermoDamage in matlab)

        // local variables
        double deltaT;      // Temperature difference between ice bottom and the snow-ice interface

        // Newly formed ice is undamaged - calculate damage as a weighted
        // average of the old damage and 0, weighted with volume.
        if ( M_thick[i] > old_vol )
            M_damage[i] = M_damage[i]*old_vol/M_thick[i];

        // SYL: manipulation of the ridged ice missing??

        // Set time_relaxation_damage to be inversely proportional to
        // temperature difference between bottom and snow-ice interface
        if ( M_thick[i] > 0. )
        {
            deltaT = std::max(0., physical::mu*M_sss[i] - M_tsurf[i] )
                / ( 1. + physical::ki*M_snow_thick[i]/(physical::ks*M_thick[i]) );
            M_time_relaxation_damage[i] = std::max(time_relaxation_damage*deltaT_relaxation_damage/deltaT, time_step);
        } else {
            M_time_relaxation_damage[i] = time_relaxation_damage;
        }
        // -------------------------------------------------

    }// end for loop
}// end thermo function

// The slab ice thermodynamics model
// This is Semtner zero layer
// We should add a new one (like Winton) at some point
void
FiniteElement::thermoIce0(int i, double wspeed, double sphuma, double conc, double voli, double vols,
        double &hi, double &hs, double &hi_old, double &Qio, double &del_hi, double &Tsurf, double tmp_Qlw_in, double tmp_snowfr)
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

    const double ai=6.1115e2, bi=23.036, ci=279.82, di=333.7;
    const double Ai=2.2e-4, Bi=3.83e-6, Ci=6.4e-10;

    const double alpha=0.62197, beta=0.37803;

    // local variables
    double Qai = 0;

    /* Don't do anything if there's no ice */
    if ( conc <=0. )
    {
        hi      = 0.;
        hi_old  = 0.;
        hs      = 0.;
        Tsurf   = 0.;
        Qio     = 0.;
        del_hi  = 0.;
        Qai     = 0.;
    } else {
        /* Calculate the slab thickness */
        hi     = voli/conc;
        hi_old = hi;
        hs     = vols/conc;

        /* Local variables */
        double albedo;

        double albs, albi, frac_sn;

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
        double dtsurf   = 1.;
        double tbot     = physical::mu*M_sss[i];
        int nb_iter_while=0;
        while ( dtsurf > 1e-4 )
        {
            nb_iter_while++;
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
                    /* Albedo scheme from ccsm 3 */
                    /* The scheme is simplified using the assumption that visible solar
                     * radiation accounts for 52% of the spectrum and the near-infrared for
                     * 48% (the same split as used in cice when running in stand-alone
                     * mode). */

                    /* This is the ccsm3 scheme when alb_ice = 0.538 and alb_sn = 0.8256 */
                    if ( Tsurf > -1. )
                    {
                        albi = alb_ice - 0.075*(Tsurf+1.);
                        albs = alb_sn  - 0.124*(Tsurf+1.);
                    } else {
                        albi = alb_ice;
                        albs = alb_sn;
                    }

                    /* Snow cover fraction */
                    frac_sn = hs/(hs+0.02);

                    /* Final albedo */
                    albedo = frac_sn*albs + (1.-frac_sn)*albi;

                    break;
                default:
                    std::cout << "alb_scheme = " << alb_scheme << "\n";
                    throw std::logic_error("Wrong albedo_scheme");
            }

            /* Calculate atmospheric fluxes */
            Qsw = -M_Qsw_in[i]*(1.-albedo);

            // -------------------------------------------------
            // 4.1) BULK ICE
            /* Out-going long-wave flux and derivative */
            Qlw_out =   physical::eps * physical::sigma_sb * std::pow(Tsurf+physical::tfrwK,4);
            dQlwdT  = 4.*physical::eps * physical::sigma_sb * std::pow(Tsurf+physical::tfrwK,3);

            // -------------------------------------------------
            // calcSphumI
            /* Specific humidity - ice surface */
            fi      = 1. + Ai + M_mslp[i]*1e-2*( Bi + Ci*Tsurf*Tsurf );
            esti    = ai*std::exp( (bi-Tsurf/di)*Tsurf/(Tsurf+ci) );
            sphumi = alpha*fi*esti/(M_mslp[i]-beta*fi*esti);

            /* We need the derivative of sphumi wrt. tsurf */
            dsphumdesti = alpha/(M_mslp[i]-beta*fi*esti)*( 1. + beta*fi*esti/(M_mslp[i]-beta*fi*esti) );
            destidT     = ( bi*ci*di-Tsurf*( 2.*ci+Tsurf) )/( di*std::pow(ci+Tsurf,2) )*esti;
            dfidT       = 2.*Ci*Bi*Tsurf;
            dsphumidT   = dsphumdesti*(fi*destidT+esti*dfidT);

            // -------------------------------------------------

            /* Density of air */
            tairK  = M_tair[i] + physical::tfrwK;
            rhoair = M_mslp[i]/(physical::Ra*tairK) * (1.+sphuma)/(1.+1.609*sphuma);

            /* Sensible heat flux and derivative */
            Qsh    = drag_ice_t * rhoair * physical::cpa * wspeed*( Tsurf - M_tair[i] );
            dQshdT = drag_ice_t * rhoair * physical::cpa * wspeed;

            /* Latent heat flux and derivative */
            Qlh    = drag_ice_t*rhoair*(physical::Lf+physical::Lv0)*wspeed*( sphumi - sphuma );
            dQlhdT = drag_ice_t*(physical::Lf+physical::Lv0)*rhoair*wspeed*dsphumidT;

            /* Sum them up */
            Qout    = Qlw_out + Qsh + Qlh;
            dQaidT = dQlwdT + dQshdT + dQlhdT;

            /* Sublimation */
            subl    = Qlh/(physical::Lf+physical::Lv0);

            /* Sum them up */
            Qai = Qsw - tmp_Qlw_in + Qout;

            // -------------------------------------------------
            /* Recalculate Tsurf */
            dtsurf = Tsurf;
            Qic    = physical::ks*( tbot-Tsurf )/( hs + physical::ks*hi/physical::ki );
            Tsurf = Tsurf + ( Qic - Qai )/
                ( physical::ks/(hs+physical::ks*hi/physical::ki) + dQaidT );

            /* Set Tsurf to the freezing point of snow or ice */
            if ( hs > 0. )
                Tsurf = std::min(0., Tsurf);
            else
                Tsurf = std::min(physical::mu*physical::si, Tsurf);

            /* Re-evaluate the exit condition */
            dtsurf = std::abs(dtsurf-Tsurf);
        }

        if(nb_iter_while>10)
        {
            LOG(DEBUG) << "nb_iter_while = " << nb_iter_while << "\n";
            throw std::logic_error("nb_iter_while larger than 10");
        }

        /* Conductive flux through the ice */
        Qic = physical::ks*( tbot-Tsurf )/( hs + physical::ks*hi/physical::ki );

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
        hs  = hs + del_hs + M_precip[i]*tmp_snowfr/physical::rhos*time_step;

        /* Heatflux from ocean */
        /* Use all excess heat to melt or grow ice. This is not
         * accurate, but will have to do for now! */
        Qio = (M_sst[i]-tbot)*physical::rhow*physical::cpw*M_mld[i]/time_step;
        /* Bottom melt/growth */
        del_hb = (Qic-Qio)*time_step/qi;

        /* Combine top and bottom */
        del_hi = del_ht+del_hb;
        hi     = hi + del_hi;

        /* Make sure we don't get too small hi_new */
        if ( hi < physical::hmin )
        {
            del_hi  = del_hi-hi;
            Qio     = Qio + hi*qi/time_step + hs*qs/time_step;

            hi      = 0.;
            hs      = 0.;
            Tsurf   = 0.;
        }

        /* Snow-to-ice conversion */
        draft = ( hi*physical::rhoi + hs*physical::rhos ) / physical::rhow;
        if ( flooding && draft > hi )
        {
            /* Subtract the mass of snow converted to ice from hs_new */
            hs = hs - ( draft - hi )*physical::rhoi/physical::rhos;
            hi = draft;
        }
    }
}

// This is the main working function, called from main.cpp (same as perform_simul in the old code)
void
FiniteElement::run()
{
    if (M_comm.rank() == 0)
        LOG(INFO) << "-----------------------Simulation started on "<< current_time_local() <<"\n";

    std::string current_time_system = current_time_local();

    // Initialise everything that doesn't depend on the mesh (constants, data set description, and time)
    this->initConstant();
    current_time = time_init /*+ pcpt*time_step/(24*3600.0)*/;
    this->initDatasets();

    int pcpt = 0;
    int niter = vm["simul.maxiteration"].as<int>();
    //std::cout<<"MAXITER= "<< maxiter <<"\n";

    if (M_rank==0)
    {
        LOG(DEBUG) <<"TIMESTEP= "<< time_step <<"\n";
        LOG(DEBUG) <<"DURATION= "<< duration <<"\n";
    }

    double displacement_factor = 1.;
    double minang = 0.;
    bool is_running = true;

    // Initialise the mesh
    this->initMesh();

#if 1
    // Check the minimum angle of the grid
    minang = this->minAngle(M_mesh);

    if (minang < vm["simul.regrid_angle"].as<double>())
    {
        LOG(INFO) <<"invalid regridding angle: should be smaller than the minimal angle in the intial grid\n";
        throw std::logic_error("invalid regridding angle: should be smaller than the minimal angle in the intial grid");
    }

    bool use_restart   = false;//vm["setup.use_restart"].as<bool>();
    bool write_restart = false;//vm["setup.write_restart"].as<bool>();
    if (0)//( use_restart )
    {
        //pcpt = this->readRestart(vm["setup.step_nb"].as<int>());
        //current_time = time_init + pcpt*time_step/(24*3600.0);

        LOG(DEBUG) <<"Initialize forcingAtmosphere\n";
        this->forcingAtmosphere();

        LOG(DEBUG) <<"Initialize forcingOcean\n";
        this->forcingOcean();

        LOG(DEBUG) <<"Initialize bathymetry\n";
        this->bathymetry();

        chrono.restart();
        LOG(DEBUG) <<"check_and_reload starts\n";
        for ( auto it = M_external_data.begin(); it != M_external_data.end(); ++it )
            (*it)->check_and_reload(M_mesh,time_init);
        LOG(DEBUG) <<"check_and_reload in "<< chrono.elapsed() <<"s\n";
    }
    else
    {
        // Do one regrid to get the mesh right
        //this->regrid(pcpt);

        // this->initVariables();
        // this->initModelState();

        // Initialise variables
        chrono.restart();
        LOG(DEBUG) <<"Initialize variables\n";
        this->initVariables();

        LOG(DEBUG) <<"Initialize forcingAtmosphere\n";
        this->forcingAtmosphere();

        LOG(DEBUG) <<"Initialize forcingOcean\n";
        this->forcingOcean();

        LOG(DEBUG) <<"Initialize bathymetry\n";
        this->bathymetry();

        chrono.restart();
        LOG(DEBUG) <<"check_and_reload starts\n";
        for ( auto it = M_external_data.begin(); it != M_external_data.end(); ++it )
            (*it)->check_and_reload(M_mesh,time_init);
        LOG(DEBUG) <<"check_and_reload in "<< chrono.elapsed() <<"s\n";

        this->initModelState();
        LOG(DEBUG) <<"initSimulation done in "<< chrono.elapsed() <<"s\n";
    }

    // for (int i=0; i<M_num_elements; ++i)
    // {
    //     if (M_conc[i] > 0. )
    //     {
    //         M_conc[i] = 1.;
    //         //M_thick[i] = 0.;
    //     }
    // }

    // Open the output file for drifters
    // TODO: Is this the right place to open the file?
    if (M_drifter_type == setup::DrifterType::IABP )
    {
        // We should tag the file name with the init time in case of a re-start.
        std::stringstream filename;
        filename << Environment::nextsimDir().string() << "/matlab/drifters_out_" << current_time << ".txt";
        M_drifters_out.open(filename.str(), std::fstream::out);
        if ( ! M_drifters_out.good() )
            throw std::runtime_error("Cannot write to file: " + filename.str());
    }

    // // Initialise the moorings - if requested
    // if ( M_use_moorings )
    //     M_grid_size = this->initMoorings(M_ncols, M_nrows);

#endif

    int rg_cpt = 0;

#if 1
    // main loop for nextsim program
    while (is_running)
    {
        M_comm.barrier();

        if (M_rank == 0)
        {
            std::cout <<"---------------------- TIME STEP "<< pcpt << " : "
                      << time_init << " + "<< pcpt*time_step/days_in_sec;

            if(fmod(pcpt*time_step, ptime_step) == 0)
            {
                std::cout <<" ---------- progression: ("<< 100.0*(pcpt*time_step/duration) <<"%)"
                          <<" ---------- time spent: "<< time_spent(current_time_system);
            }

            std::cout <<"\n";
        }


        is_running = ((pcpt+1)*time_step) < duration;

        if (pcpt == niter)
            is_running = false;


        current_time = time_init + pcpt*time_step/days_in_sec;

        // step 0: preparation
        // remeshing and remapping of the prognostic variables

#if 0
        minang = 0.;

        if (M_rank == 0)
            minang = this->minAngle(M_mesh_root,M_UM_root,1.);

        double min_angle = minang;

        minang = boost::mpi::all_reduce(M_comm, min_angle, boost::mpi::maximum<double>());
#endif
        // The first time step we behave as if we just did a regrid
        M_regrid = (pcpt==0);

        bool force_regrid = (!(pcpt % 30)) && (pcpt != 0);

        if (vm["simul.regrid"].as<std::string>() == "bamg")
        {
            minang = this->minAngle(M_mesh,M_UM,displacement_factor);
            //std::cout<<"[" << M_rank <<"] " <<" REGRID ANGLE= "<< minang <<"\n";

            if (M_rank == 0)
            {
                //std::cout <<"---------------------- REGRID ANGLE= "<< minang <<"\n";
                std::cout <<"REGRID ANGLE= "<< minang <<"\n";
            }

            if ( minang < vm["simul.regrid_angle"].as<double>() )
            {
                //this->exportResults(2000+rg_cpt);
                // std::cout<<"UUMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\n";
                //return;

                // if ( M_use_moorings )
                //     this->updateMoorings(M_grid_size, M_ncols, M_nrows);

                LOG(DEBUG) <<"Regriding starts\n";
				chrono.restart();
                this->regrid(pcpt);
                LOG(DEBUG) <<"Regriding done in "<< chrono.elapsed() <<"s\n";

                M_regrid = true;

                //this->exportResults(3000+rg_cpt);
                //return;
                ++rg_cpt;
            }
        }

        M_comm.barrier();

        // Read in the new buoys and output
        // if (M_drifter_type == setup::DrifterType::IABP && std::fmod(current_time,0.5) == 0)
        // {
        //     this->updateIABPDrifter();
        //     // TODO: Do we want to output drifters at a different time interval?
        //     this->outputDrifter(M_drifters_out);
        // }

        if ( M_regrid || use_restart )
        {
            this->tensors();
            //this->tensorsOnRoot();
            this->cohesion();
            this->coriolis();
        }

        for ( auto it = M_external_data.begin(); it != M_external_data.end(); ++it )
            (*it)->check_and_reload(M_mesh,current_time+time_step/(24*3600.0));


        if (pcpt == 0)
        {
            this->exportResults(0);
        }

        //======================================================================
        // Do the thermodynamics
        //======================================================================
        if(vm["simul.use_thermo_forcing"].as<bool>())
        {
            chrono.restart();
            LOG(DEBUG) <<"thermo starts\n";
            this->thermo();
            LOG(DEBUG) <<"thermo done in "<< chrono.elapsed() <<"s\n";
        }

        //======================================================================
        // Assemble the matrix
        //======================================================================
        this->assemble(pcpt);

        //======================================================================
        // Solve the linear problem
        //======================================================================
        this->solve();


        this->updateVelocity();


        this->update();
        //this->updateOnRoot();

        if(fmod((pcpt+1)*time_step,output_time_step) == 0)
        {
            this->exportResults((int) (pcpt+1)*time_step/output_time_step);
        }

        ++pcpt;

        //this->exportResults(pcpt);

#if 0
        if(fmod((pcpt+1)*time_step,output_time_step) == 0)
        {
            chrono.restart();
            std::cout<<"export starts\n";
            this->exportResults((int) (pcpt+1)*time_step/output_time_step);
            std::cout<<"export done in " << chrono.elapsed() <<"s\n";
        }

        if ( fmod(pcpt*time_step,restart_time_step) == 0)
        {
            std::cout << "Writing restart file after time step " <<  pcpt-1 << endl;
            this->writeRestart(pcpt, (int) pcpt*time_step/restart_time_step);
        }
#endif
    }

    this->exportResults(1000);

    // if (M_rank==0)
    //     std::cout<<"TIMER total = " << chrono_tot.elapsed() <<"s\n";

    // // Don't forget to close the iabp file!
    // if (M_drifter_type == setup::DrifterType::IABP)
    // {
    //     M_iabp_file.close();
    //     drifters_out.close();
    // }

    this->clear();

#endif

    if (M_rank==0)
        LOG(INFO) << "-----------------------Simulation done on "<< current_time_local() <<"\n";
}

void
FiniteElement::writeRestart(int pcpt, int step)
{
}

void
FiniteElement::readRestart(int &pcpt, int step)
{
}

void
FiniteElement::updateVelocity()
{
    M_comm.barrier();

    // timer["updatevelocity"].first.restart();

    // if (M_rank == 0)
    //     std::cout<<"UPDATEVELOCITY STARTS\n";

    M_VTMM = M_VTM;
    M_VTM = M_VT;
    M_VT = M_solution->container();

    std::vector<double> speed_scaling;
    this->speedScaling(speed_scaling);

    // linear scaling of ice velocity
    for (int i=0; i<M_num_nodes; ++i)
    {
        M_VT[i] *= speed_scaling[i];
        M_VT[i+M_num_nodes] *= speed_scaling[i];
    }

    // if (M_rank == 0)
    //     std::cout<<"TIMER UPDATEVELOCITY= " << timer["updatevelocity"].first.elapsed() <<"s\n";

    //M_speed_scaling = speed_scaling;

    // double min_elt = *std::min_element(M_VT.begin(),M_VT.end());
    // double max_elt = *std::max_element(M_VT.begin(),M_VT.end());

    // // std::cout<<"----------------------------[" << M_rank <<"] " <<" VT MIN= "<< min_elt <<"\n";
    // // std::cout<<"----------------------------[" << M_rank <<"] " <<" VT MAX= "<< max_elt <<"\n";

    // M_comm.barrier();

    // double gmin = boost::mpi::all_reduce(M_comm, min_elt, boost::mpi::minimum<double>());
    // double gmax = boost::mpi::all_reduce(M_comm, max_elt, boost::mpi::maximum<double>());

    // if (M_comm.rank()==0)
    // {
    //     std::cout<<"----------------------------VT MIN= "<< gmin <<"\n";
    //     std::cout<<"----------------------------VT MAX= "<< gmax <<"\n";
    // }
}

void
FiniteElement::speedScaling(std::vector<double>& speed_scaling)
{
    M_comm.barrier();

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
        //auto rmap_nodes = M_mesh.mapNodes();
        auto rmap_elements = M_mesh.mapElements();

        auto conc_root_nrd = conc_root;
        int gsize = M_mesh_root.numTriangles();

        for (int i=0; i<gsize; ++i)
        {
            int ri = rmap_elements.left.find(i+1)->second-1;
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
                }
                else
                {
                    break;
                }
            }

            c_max_nodal_neighbour = *std::max_element(cloc_elts.begin(),cloc_elts.begin()+j-1);
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

    M_comm.barrier();
}

void
FiniteElement::error()
{
    double l2_error = 0;
    double sh1_error = 0;

    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        double area = measure(*it,M_mesh);
        std::vector<double> x(3);
        std::vector<double> y(3);

        double l2_contrib = 0;
        double sh1_contrib = 0;
        double entry_contrib = 0;

        for (int i=0; i<3; ++i)
        {
            // x[i] = M_nodes.find(it->second.indices[i])->second.coords[0];
            // y[i] = M_nodes.find(it->second.indices[i])->second.coords[1];

            x[i] = M_nodes[it->indices[i]-1].coords[0];
            y[i] = M_nodes[it->indices[i]-1].coords[1];
        }

        int lc = 0;
        for (int j=0; j<3; ++j)
        {
            l2_contrib = 0;
            sh1_contrib = 0;
            // x-axis
            int jp1 = (j+1)%3;
            int jp2 = (j+2)%3;

            for (int k=0; k<3; ++k)
            {
                // y-axis
                int kp1 = (k+1)%3;
                int kp2 = (k+2)%3;

                // semi h1 error
                entry_contrib = (y[jp1]-y[jp2])*(y[kp1]-y[kp2])+(x[jp1]-x[jp2])*(x[kp1]-x[kp2]);
                sh1_contrib += entry_contrib*M_exact->operator()(it->indices[k]-1)/(4.0*area);

                // l2 error
                l2_contrib += M_exact->operator()(it->indices[k]-1)*((j == k) ? 2.0 : 1.0)*area/12.0;
            }

            l2_error += l2_contrib*M_exact->operator()(it->indices[j]-1);
            sh1_error += sh1_contrib*M_exact->operator()(it->indices[j]-1);
        }
    }

    std::cout<<"||u-uh||_L2  = "<< std::sqrt(l2_error) <<"\n";
    std::cout<<"||u-uh||_H1  = "<< std::sqrt(l2_error+sh1_error) <<"\n";
}

void
FiniteElement::forcingAtmosphere()//(double const& u, double const& v)
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
                &M_asr_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(&M_asr_elements_dataset,M_mesh,0,false);
            M_external_data.push_back(&M_tair);

            M_mixrat=ExternalData(&M_asr_elements_dataset,M_mesh,1,false);
            M_external_data.push_back(&M_mixrat);

            M_mslp=ExternalData(&M_asr_elements_dataset,M_mesh,2,false);
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(&M_asr_elements_dataset,M_mesh,3,false);
            M_external_data.push_back(&M_Qsw_in);

            M_Qlw_in=ExternalData(&M_asr_elements_dataset,M_mesh,4,false);
            M_external_data.push_back(&M_Qlw_in);

            M_snowfr=ExternalData(&M_asr_elements_dataset,M_mesh,5,false);
            M_external_data.push_back(&M_snowfr);

            M_precip=ExternalData(&M_asr_elements_dataset,M_mesh,6,false);
            M_external_data.push_back(&M_precip);

            M_dair=ExternalData(-1.);
            M_external_data.push_back(&M_dair);
        break;

        case setup::AtmosphereType::ERAi:
            M_wind=ExternalData(
                &M_ERAi_nodes_dataset,M_mesh,0,true ,
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

            M_tair=ExternalData(&M_ERAi_elements_dataset,M_mesh,0,false);
            M_external_data.push_back(&M_tair);

            M_dair=ExternalData(&M_ERAi_elements_dataset,M_mesh,1,false);
            M_external_data.push_back(&M_dair);

            M_mslp=ExternalData(&M_ERAi_elements_dataset,M_mesh,2,false);
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(&M_ERAi_elements_dataset,M_mesh,3,false);
            M_external_data.push_back(&M_Qsw_in);

            M_tcc=ExternalData(&M_ERAi_elements_dataset,M_mesh,4,false);
            M_external_data.push_back(&M_tcc);

            M_precip=ExternalData(&M_ERAi_elements_dataset,M_mesh,5,false);
            M_external_data.push_back(&M_precip);

            M_mixrat=ExternalData(-1.);
            M_external_data.push_back(&M_mixrat);

        break;

        default:
            std::cout << "invalid wind forcing"<<"\n";
            throw std::logic_error("invalid wind forcing");
    }
}

void
FiniteElement::forcingOcean()//(double const& u, double const& v)
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

            M_ocean_temp=ExternalData(vm["simul.constant_ocean_temp"].as<double>());
            M_external_data.push_back(&M_ocean_temp);

            M_ocean_salt=ExternalData(vm["simul.constant_ocean_salt"].as<double>());
            M_external_data.push_back(&M_ocean_salt);

            M_mld=ExternalData(vm["simul.constant_mld"].as<double>());
            M_external_data.push_back(&M_mld);
            break;
        case setup::OceanType::TOPAZR:
            M_ocean=ExternalData(
                &M_topaz_nodes_dataset, M_mesh, 0, true,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ocean);

            M_ssh=ExternalData(
                &M_topaz_nodes_dataset, M_mesh, 2, false,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ssh);

            M_ocean_temp=ExternalData(&M_topaz_elements_dataset, M_mesh, 0,false);
            M_external_data.push_back(&M_ocean_temp);

            M_ocean_salt=ExternalData(&M_topaz_elements_dataset, M_mesh, 1,false);
            M_external_data.push_back(&M_ocean_salt);

            M_mld=ExternalData(&M_topaz_elements_dataset, M_mesh, 2,false);
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
            std::fill(M_sst.begin(), M_sst.end(), -1.8);
            std::fill(M_sss.begin(), M_sss.end(), -1.8/physical::mu);
            break;
        case setup::OceanType::TOPAZR:
            for ( int i=0; i<M_num_elements; ++i)
            {
                // Make sure the erroneous salinity and temperature don't screw up the initialisation too badly
                // This can still be done much better!
                M_sss[i] = std::max(physical::si, M_ocean_salt[i]);
                M_sst[i] = std::max(M_sss[i]*physical::mu, M_ocean_temp[i]);
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
        case setup::IceType::TARGET:
            this->targetIce();
            break;
        case setup::IceType::TOPAZ4:
            this->topazIce();
            break;
        case setup::IceType::AMSRE:
            this->amsreIce();
            break;
        case setup::IceType::OSISAF:
            this->osisaf2Ice();
            break;
        case setup::IceType::AMSR2:
            this->amsr2Ice();
            break;

        default:
            std::cout << "invalid initialization of the ice"<<"\n";
            throw std::logic_error("invalid initialization of the ice");
    }
}

void
FiniteElement::constantIce()
{
	LOG(DEBUG) <<"Constant Ice\n";
    std::fill(M_conc.begin(), M_conc.end(), vm["simul.init_concentration"].as<double>());
    std::fill(M_thick.begin(), M_thick.end(), vm["simul.init_thickness"].as<double>());

    std::fill(M_snow_thick.begin(), M_snow_thick.end(), vm["simul.init_snow_thickness"].as<double>());
    std::fill(M_damage.begin(), M_damage.end(), 0.);
}

void
FiniteElement::targetIce()
{
    double y_max=300000.;
    double x_max=350000.;
    double x_min=200000.;

	double tmp_var;

    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();

    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_var = (RY[i]<=y_max)*(RX[i]<=x_max)*(RX[i]>=x_min);

        std::cout<<"RX: "<< RX[i] << "RY: "<< RY[i] << "tmp_var: " << tmp_var << "\n";

        M_conc[i]  = vm["simul.init_concentration"].as<double>()*tmp_var;
		M_thick[i] = vm["simul.init_thickness"].as<double>()*tmp_var;
		M_snow_thick[i] = vm["simul.init_snow_thickness"].as<double>()*tmp_var;
        M_damage[i]=0.;

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
            M_damage[i]=1.;
        }
        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
            M_damage[i]=1.;
        }
    }
}

void
FiniteElement::topazIce()
{
    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false);
    M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		tmp_var=M_init_conc[i];
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
FiniteElement::amsreIce()
{
    double real_thickness, init_conc_topaz_tmp;

    external_data M_init_conc=ExternalData(&M_ice_amsre_elements_dataset,M_mesh,0,false);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_conc_topaz=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false);
    M_init_conc_topaz.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false);
    M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		M_conc[i] = M_init_conc[i];

        // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_conc_topaz[i];
		init_conc_topaz_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

        // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
        if(init_conc_topaz_tmp>0.05)
        {
            real_thickness=M_thick[i]/init_conc_topaz_tmp;
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
FiniteElement::osisaf2Ice()
{
    double real_thickness, init_conc_topaz_tmp;

    external_data M_init_conc=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,0,false);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_conc_topaz=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false);
    M_init_conc_topaz.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false);
    M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		M_conc[i] = M_init_conc[i];

        // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_conc_topaz[i];
		init_conc_topaz_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

        // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
        if(init_conc_topaz_tmp>0.05)
        {
            real_thickness=M_thick[i]/init_conc_topaz_tmp;
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
FiniteElement::amsr2Ice()
{
    double real_thickness, init_conc_topaz_tmp;

    external_data M_init_conc=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_conc_topaz=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false);
    M_init_conc_topaz.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false);
    M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		M_conc[i] = M_init_conc[i];

        // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_conc_topaz[i];
		init_conc_topaz_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

        // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
        if(init_conc_topaz_tmp>0.05)
        {
            real_thickness=M_thick[i]/init_conc_topaz_tmp;
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
    switch (M_drifter_type)
    {
        case setup::DrifterType::NONE:
            break;

        case setup::DrifterType::EQUALLYSPACED:
            this->equallySpacedDrifter();
            break;

        case setup::DrifterType::IABP:
            this->initIABPDrifter();
            break;

        default:
            std::cout << "invalid initialization of drifter"<<"\n";
            throw std::logic_error("invalid initialization of drifter");
    }
}

void
FiniteElement::coriolis()
{
    // Interpolation of the latitude
    std::vector<double> lat = M_mesh.meanLat();

    for (int i=0; i<M_fcor.size(); ++i)
    {
        M_fcor[i] = 2*(vm["simul.omega"].as<double>())*std::sin(lat[i]*PI/180.);
        //std::cout <<"Coeff= "<< M_fcor[i] <<"\n";
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
            M_element_depth=ExternalData(&M_etopo_elements_dataset,M_mesh,0,false);
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
                              % "NpsNextsim.mpp"
                              ).str();

    std::vector<char> str(configfile.begin(), configfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    // Assemble the coordinates from the unordered_map
    std::vector<double> drifter_X(M_drifter.size());
    std::vector<double> drifter_Y(M_drifter.size());
    int j=0;
    for ( auto it = M_drifter.begin(); it != M_drifter.end(); ++it )
    {
        drifter_X[j] = it->second[0];
        drifter_Y[j] = it->second[1];
        ++j;
    }

    // Interpolate the velocity onto the drifter positions
    int nb_var=2;
    std::vector<double> interp_drifter_in(nb_var*M_mesh.numNodes());
    double* interp_drifter_out;

    for (int i=0; i<M_mesh.numNodes(); ++i)
    {
        interp_drifter_in[nb_var*i]   = M_UM[i];
        interp_drifter_in[nb_var*i+1] = M_UM[i+M_mesh.numNodes()];
    }

    // Interpolate the velocity
    InterpFromMeshToMesh2dx(&interp_drifter_out,
        &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
        M_mesh.numNodes(),M_mesh.numTriangles(),
        &interp_drifter_in[0],
        M_mesh.numNodes(),nb_var,
        &drifter_X[0],&drifter_Y[0],M_drifter.size(),
        false);

    // Loop over the map and output
    j=0;
    boost::gregorian::date           date = Nextsim::parse_date( current_time );
    boost::posix_time::time_duration time = Nextsim::parse_time( current_time );
    for ( auto it = M_drifter.begin(); it != M_drifter.end(); ++it )
    {
        double lat, lon;
        inverse_mapx(map,it->second[0]+interp_drifter_out[j],it->second[1]+interp_drifter_out[j+M_drifter.size()],&lat,&lon);
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
                              % "NpsNextsim.mpp"
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
        time = dateStr2Num(date) + hour/24.;

        // Remember which buoys are in the ice according to IABP
        keepers.push_back(number);

        // Project and add the buoy to the map if it's missing
        if ( M_drifter.count(number) == 0 )
        {
            double x, y;
            forward_mapx(map,lat,lon,&x,&y);
            M_drifter.emplace(number, std::array<double,2>{x, y});

        }
    }
    close_mapx(map);

    // Go through the M_drifter map and throw out the ones which IABP doesn't
    // report as being in the ice anymore
    for ( auto model = M_drifter.begin(); model != M_drifter.end(); /* ++model is not allowed here, because we use 'erase' */ )
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
            model = M_drifter.erase(model);
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
    double time = dateStr2Num("1979-01-01");
    while ( time < time_init )
    {
        // Remember where we were
        pos = M_iabp_file.tellg();

        // Read the next line
        int year, month, day, hour, number;
        double lat, lon;
        M_iabp_file >> year >> month >> day >> hour >> number >> lat >> lon;
        std::string date = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);

        time = dateStr2Num(date) + hour/24.;
    }

    // We must rewind one line so that updateIABPDrifter works correctly
    M_iabp_file.seekg(pos);
}

void
FiniteElement::equallySpacedDrifter()
{
    std::cout << "equally spaced drifters not yet implemented" << endl;
    throw std::logic_error("equally spaced drifters not yet implemented");
}

void
FiniteElement::importBamg(BamgMesh const* bamg_mesh)
{
    //mesh_type_root mesh;
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
        mesh_triangles.push_back(gmshElt);
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
FiniteElement::createGraph()//(BamgMesh const* bamg_mesh)
{
    M_comm.barrier();

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
        if (std::find(M_local_ghost.begin(),M_local_ghost.end(),gid) == M_local_ghost.end())
        {
            for (int j=0; j<Ncc; ++j)
            {
                int currentr = bamgmesh->NodalConnectivity[Nd*i+j];

                int gid2 = M_transfer_map.right.find(currentr)->second;
                if (std::find(M_local_ghost.begin(),M_local_ghost.end(),gid2) == M_local_ghost.end())
                    ++counter_dnnz;
                else
                    ++counter_onnz;

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

    M_comm.barrier();

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

    //std::vector<int> sizes_nodes(M_comm.size());
    //boost::mpi::gather(M_comm, M_local_ndof, sizes_nodes, 0);

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
    std::vector<int> sizes_elements = M_sizes_elements;//(M_comm.size());
    //boost::mpi::gather(M_comm, M_local_nelements, sizes_elements, 0);

    // ELEMENT INTERPOLATION With Cavities
    int nb_var_element=8;//15;
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

        // divergence_rate
        interp_elt_in_local[nb_var_element*i+tmp_nb_var] = M_divergence_rate[i];
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
    std::vector<double> divergence(M_mesh_root.numTriangles());
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
        auto rmap_elements = M_mesh.mapElements();
        //auto interp_in_elements_nrd = interp_in_elements;

        for (int i=0; i<M_mesh_root.numTriangles(); ++i)
        {
            tmp_nb_var=0;
            int ri = rmap_elements.left.find(i+1)->second-1;

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

            // divergence
            divergence[i] = interp_in_elements[nb_var_element*ri+tmp_nb_var];
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
        Exporter exporter;
        std::string fileout;

        if (export_mesh)
        {
            fileout = (boost::format( "%1%/matlab/mesh_%2%_%3%.bin" )
                       % Environment::nextsimDir().string()
                       % M_rank
                       % step ).str();

            std::cout<<"MESH BINARY: Exporter Filename= "<< fileout <<"\n";

            // move the mesh for the export
            //M_mesh.move(M_UM,1.);

            std::fstream meshbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
            if ( ! meshbin.good() )
                throw std::runtime_error("Cannot write to file: " + fileout);

            if (0)
            {
                auto rmap_nodes = M_mesh.mapNodes();
                auto rmap_elements = M_mesh.mapElements();

                auto meshr = M_mesh_root;
                meshr.reorder(rmap_nodes,rmap_elements);
                //exporter.writeMesh(meshbin, M_mesh);
                //exporter.writeMesh(meshbin, M_mesh_root);
                exporter.writeMesh(meshbin, meshr);
                meshbin.close();
            }

            exporter.writeMesh(meshbin, M_mesh_root);
            meshbin.close();

            // move it back after the export
            //M_mesh.move(M_UM,-1.);

            fileout = (boost::format( "%1%/matlab/mesh_%2%_%3%.dat" )
                       % Environment::nextsimDir().string()
                       % M_rank
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
                   % M_rank
                   % step ).str();

        std::cout<<"BINARY: Exporter Filename= "<< fileout <<"\n";

        std::fstream outbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
        if ( ! outbin.good() )
            throw std::runtime_error("Cannot write to file: " + fileout);


        auto rmap_nodes = M_mesh.mapNodes();
        int prv_global_num_nodes = M_mesh.numGlobalNodes();

        auto interp_in_nodes_nrd = vt_root;

        for (int i=0; i<prv_global_num_nodes; ++i)
        {
            int ri =  rmap_nodes.left.find(i+1)->second-1;

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

        exporter.writeField(outbin, vt_root, "M_VT");
        exporter.writeField(outbin, conc_root, "Concentration");
        exporter.writeField(outbin, thick_root, "Thickness");
        exporter.writeField(outbin, thick_snow, "Snow");
        exporter.writeField(outbin, stress1, "Stress1");
        exporter.writeField(outbin, stress2, "Stress2");
        exporter.writeField(outbin, stress3, "Stress3");
        exporter.writeField(outbin, damage, "Damage");
        exporter.writeField(outbin, divergence, "Divergence");
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
            double const rtanalpha = c_thin_max/h_thin_max;
            for ( int i=0; i<M_mesh.numTriangles(); ++i )
            {
                conc_thin[i] = std::min(std::min(M_h_thin[i]/physical::hmin, std::sqrt(2.*M_h_thin[i]*rtanalpha)), 1.-M_conc[i]);
            }
            exporter.writeField(outbin, conc_thin, "Concentration_thin_ice");
        }

        outbin.close();


        fileout = (boost::format( "%1%/matlab/field_%2%_%3%.dat" )
                   % Environment::nextsimDir().string()
                   % M_rank
                   % step ).str();

        std::cout<<"RECORD FIELD: Exporter Filename= "<< fileout <<"\n";

        std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
        if ( ! outrecord.good() )
            throw std::runtime_error("Cannot write to file: " + fileout);

        exporter.writeRecord(outrecord);
        outrecord.close();
    }
}


void
FiniteElement::tensorsOnRoot()
{
    if (M_rank == 0)
    {
        M_B0T_root.resize(M_mesh_root.numTriangles());
        std::vector<double> B0T(18,0.);

        //int num_elements = M_mesh_root.numTriangles();
        auto elements_root = M_mesh_root.triangles();

        int cpt = 0;
        for (auto it=elements_root.begin(), end=elements_root.end(); it!=end; ++it)
        {
            std::vector<double> shapecoeff = this->shapeCoeff(*it,M_mesh_root);

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

            M_B0T_root[cpt]                = B0T;
            ++cpt;
        }
    }
}

void
FiniteElement::updateOnRoot()
{
    M_comm.barrier();

    int nb_var_node = 2;
    std::vector<double> VT_local(nb_var_node*M_local_ndof,0.);

    for (int i=0; i<M_local_ndof; ++i)
    {
        // VT
        VT_local[nb_var_node*i] = M_VT[i];
        VT_local[nb_var_node*i+1] = M_VT[i+M_num_nodes];
    }

    std::vector<int> sizes_nodes = M_sizes_nodes;
    //boost::mpi::gather(M_comm, M_local_ndof, sizes_nodes, 0);
    std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = nb_var_node*f; });

    auto rmap_nodes = M_mesh.mapNodes();
    int global_num_nodes = M_mesh.numGlobalNodes();

    std::vector<double> VT_root;

    if (M_rank == 0)
    {
        VT_root.resize(nb_var_node*global_num_nodes);
        boost::mpi::gatherv(M_comm, VT_local, &VT_root[0], sizes_nodes, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, VT_local, 0);
    }

    std::vector<double> UM_P;

    if (M_rank == 0)
    {
        auto VT_root_nrd = VT_root;

        for (int i=0; i<global_num_nodes; ++i)
        {
            int ri =  rmap_nodes.left.find(i+1)->second-1;

            VT_root[i] = VT_root_nrd[nb_var_node*ri+0];
            VT_root[i+global_num_nodes] = VT_root_nrd[nb_var_node*ri+1];
        }

        double vnorm = std::accumulate(VT_root.begin(), VT_root.end(), 0.);

        // move all node except to neumann
        UM_P.resize(M_UM_root.size());
        UM_P = M_UM_root;

        for (int nd=0; nd<M_UM_root.size(); ++nd)
        {
            M_UM_root[nd] += time_step*VT_root[nd];
        }

        for (const int& nd : M_neumann_nodes_root)
        {
            M_UM_root[nd] = UM_P[nd];
        }

        double minang = this->minAngle(M_mesh_root,M_UM_root,1.);

        std::cout<<"----------------------------[" << M_rank <<"] " <<" VT MIN   = "<< *std::min_element(VT_root.begin(),VT_root.end()) <<"\n";
        std::cout<<"----------------------------[" << M_rank <<"] " <<" VT MAX   = "<< *std::max_element(VT_root.begin(),VT_root.end()) <<"\n";
        std::cout<<"----------------------------[" << M_rank <<"] " <<" VT NORM  = "<< vnorm <<"\n";
        std::cout<<"----------------------------[" << M_rank <<"] " <<" MIN ANGLE= "<< minang <<"\n";
    }

#if 1

    std::vector<int> sizes_elements = M_sizes_elements;
    //boost::mpi::gather(M_comm, M_local_nelements, sizes_elements, 0);
    int nb_var_element=12;
    std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });

    std::vector<double> fields_local(nb_var_element*M_local_nelements);

    int tmp_nb_var=0;
    for (int i=0; i<M_local_nelements; ++i)
    {
        tmp_nb_var=0;

        // concentration
        fields_local[nb_var_element*i+tmp_nb_var] = M_conc[i];
        tmp_nb_var++;

        // thickness
        fields_local[nb_var_element*i+tmp_nb_var] = M_thick[i];
        tmp_nb_var++;

        // snow thickness
        fields_local[nb_var_element*i+tmp_nb_var] = M_snow_thick[i];
        tmp_nb_var++;

        // integrated_stress1
        fields_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i]/**M_thick[i]*/;
        tmp_nb_var++;

        // integrated_stress2
        fields_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i+1]/**M_thick[i]*/;
        tmp_nb_var++;

        // integrated_stress3
        fields_local[nb_var_element*i+tmp_nb_var] = M_sigma[3*i+2]/**M_thick[i]*/;
        tmp_nb_var++;

        // damage
        fields_local[nb_var_element*i+tmp_nb_var] = M_damage[i];
        tmp_nb_var++;

        // divergence_rate
        fields_local[nb_var_element*i+tmp_nb_var] = M_divergence_rate[i];
        tmp_nb_var++;

        // h_ridged_thin_ice
        fields_local[nb_var_element*i+tmp_nb_var] = M_h_ridged_thin_ice[i];
        tmp_nb_var++;

        // h_ridged_thick_ice
        fields_local[nb_var_element*i+tmp_nb_var] = M_h_ridged_thick_ice[i];
        tmp_nb_var++;

        // thin ice thickness
        fields_local[nb_var_element*i+tmp_nb_var] = M_h_thin[i];
        tmp_nb_var++;

        // snow on thin ice
        fields_local[nb_var_element*i+tmp_nb_var] = M_hs_thin[i];
        tmp_nb_var++;

        if(tmp_nb_var>nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }
    }

    std::vector<double> fields_root;

    if (M_rank == 0)
    {
        fields_root.resize(nb_var_element*M_mesh_root.numTriangles());
        boost::mpi::gatherv(M_comm, fields_local, &fields_root[0], sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, fields_local, 0);
    }


    auto rmap_elements = M_mesh.mapElements();

    std::vector<double> conc_root;
    std::vector<double> thick_root;
    std::vector<double> snow_thick_root;
    std::vector<double> sigma_root;
    std::vector<double> damage_root;
    std::vector<double> divergence_rate_root;
    std::vector<double> h_ridged_thin_ice_root;
    std::vector<double> h_ridged_thick_ice_root;
    std::vector<double> h_thin_root;
    std::vector<double> hs_thin_root;

    if (M_rank == 0)
    {
        int gsize = M_mesh_root.numTriangles();
        conc_root.resize(gsize);
        thick_root.resize(gsize);
        snow_thick_root.resize(gsize);
        sigma_root.resize(3*gsize);
        damage_root.resize(gsize);
        divergence_rate_root.resize(gsize);
        h_ridged_thin_ice_root.resize(gsize);
        h_ridged_thick_ice_root.resize(gsize);
        h_thin_root.resize(gsize);
        hs_thin_root.resize(gsize);


        tmp_nb_var=0;

        for (int i=0; i<gsize; ++i)
        {
            tmp_nb_var=0;
            int ri = rmap_elements.left.find(i+1)->second-1;

            // concentration
            conc_root[i] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // thickness
            thick_root[i] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // snow thickness
            snow_thick_root[i] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress1
            sigma_root[3*i] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress2
            sigma_root[3*i+1] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // integrated_stress3
            sigma_root[3*i+2] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // damage
            damage_root[i] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // divergence
            divergence_rate_root[i] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // h_ridged_thin_ice
            h_ridged_thin_ice_root[i] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // h_ridged_thick_ice
            h_ridged_thick_ice_root[i] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // thin ice thickness
            h_thin_root[i] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            // snow on thin ice
            hs_thin_root[i] = fields_root[nb_var_element*ri+tmp_nb_var];
            tmp_nb_var++;

            if(tmp_nb_var>nb_var_element)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }
    }

    //std::cout<<"[" << M_rank <<"] " <<" Divergence_Rate MIN= "<< *std::min_element(M_divergence_rate.begin(),M_divergence_rate.end()) <<"\n";
    //std::cout<<"[" << M_rank <<"] " <<" Divergence_Rate MAX= "<< *std::max_element(M_divergence_rate.begin(),M_divergence_rate.end()) <<"\n";



    if (M_rank == 0)
    {
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

        /* some variables used for the advection*/
        double surface, surface_new, ar, ar_new;
        double ice_surface, ice_volume, snow_volume;
        double thin_ice_surface, thin_ice_volume, thin_snow_volume;
        double ridging_thin_ice, ridging_thick_ice, ridging_snow_thin_ice;
        double ridged_thin_ice_volume, ridged_thick_ice_volume;

        /* invariant of the internal stress tensor and some variables used for the damaging process*/
        double sigma_s, sigma_n;
        double tract_max;
        double tmp, sigma_target;

        /* some variables used for the ice redistribution*/
        double tanalpha, rtanalpha, del_v, del_c, del_vs, new_v_thin;

        /* Indices loop*/
        int i,j;

        bool to_be_updated;

        double Compressive_strength = compr_strength*scale_coef;

        int num_elements = M_mesh_root.numTriangles();
        auto elements_root = M_mesh_root.triangles();

        for (int cpt=0; cpt < num_elements; ++cpt)
        {
            /* set constants for the ice redistribution */
            tanalpha  = h_thin_max/c_thin_max;
            rtanalpha = 1./tanalpha;

            old_thick = thick_root[cpt];
            old_snow_thick = snow_thick_root[cpt];
            old_conc = conc_root[cpt];
            old_damage = damage_root[cpt];
            old_h_ridged_thick_ice=h_ridged_thick_ice_root[cpt];

            if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                old_h_thin = h_thin_root[cpt];
                old_hs_thin=hs_thin_root[cpt];
                old_h_ridged_thin_ice=h_ridged_thin_ice_root[cpt];
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
                    //if (!(M_elements[cpt]).ghostNodes[j])
                    {
                        epsilon_veloc_i += M_B0T_root[cpt][i*6 + 2*j]*VT_root[(elements_root[cpt]).indices[j]-1];
                        epsilon_veloc_i += M_B0T_root[cpt][i*6 + 2*j + 1]*VT_root[(elements_root[cpt]).indices[j]-1+global_num_nodes];
                    }
                }

                epsilon_veloc[i] = epsilon_veloc_i;
            }

            divergence_rate_root[cpt]= (epsilon_veloc[0]+epsilon_veloc[1]);

            /*======================================================================
             * Update the internal stress
             *======================================================================
             */

            for(i=0;i<3;i++)
            {
                sigma_dot_i = 0.0;
                for(j=0;j<3;j++)
                {
                    sigma_dot_i += std::exp(ridging_exponent*(1.-old_conc))*young*(1.-old_damage)*M_Dunit[i*3 + j]*epsilon_veloc[j];
                }

                sigma_root[3*cpt+i] += time_step*sigma_dot_i;
                sigma_pred[i]    = sigma_root[3*cpt+i] + time_step*sigma_dot_i;
            }

            /*======================================================================
             * Correct the internal stress and the damage
             *======================================================================
             */

            /* Compute the shear and normal stress, which are two invariants of the internal stress tensor */

            sigma_s=std::hypot((sigma_pred[0]-sigma_pred[1])/2.,sigma_pred[2]);
            sigma_n=           (sigma_pred[0]+sigma_pred[1])/2.;

            /* minimum and maximum normal stress */
            tract_max=tract_coef*C_fix/tan_phi;

            /* Correction of the damage */

            if((sigma_n>tract_max) || (sigma_n<(-Compressive_strength)))
            {
                if(sigma_n>tract_max)
                {
                    sigma_target=tract_max;
                }
                else
                {
                    sigma_target=-Compressive_strength;
                }

                tmp=1.0-sigma_target/sigma_n*(1.0-old_damage);

                if(tmp>damage_root[cpt])
                {
                    damage_root[cpt]=tmp;
                }
            }

            if(sigma_s>C_fix-sigma_n*tan_phi)
            {
                tmp=1.0-C_fix/(sigma_s+sigma_n*tan_phi)*(1.0-old_damage);

                if(tmp>damage_root[cpt])
                {
                    damage_root[cpt]=tmp;
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
                    sigma_root[3*cpt+i] = (1.-damage_root[cpt])/(1.-old_damage)*sigma_root[3*cpt+i] ;
                }
                else
                {
                    sigma_root[3*cpt+i] = 0. ;
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
            tmp=1./(1.-damage_root[cpt]);
            tmp-= 1000*time_step/time_relaxation_damage;
            tmp=((tmp>1.)?(tmp):(1.));
            damage_root[cpt]=-1./tmp + 1.;

            /*======================================================================
             * Update:
             * Ice and snow thickness, and concentration using a Lagrangian or an Eulerian scheme
             *======================================================================
             */

            to_be_updated=false;
            if( divergence_rate_root[cpt]!=0.)
                to_be_updated=true;

            std::sort(M_neumann_flags_root.begin(), M_neumann_flags_root.end());

            if(std::binary_search(M_neumann_flags_root.begin(),M_neumann_flags_root.end(),(elements_root[cpt]).indices[0]-1) ||
               std::binary_search(M_neumann_flags_root.begin(),M_neumann_flags_root.end(),(elements_root[cpt]).indices[1]-1) ||
               std::binary_search(M_neumann_flags_root.begin(),M_neumann_flags_root.end(),(elements_root[cpt]).indices[2]-1))
                to_be_updated=false;

            if((old_conc>0.)  && (to_be_updated))
            {
                surface = this->measure(elements_root[cpt],M_mesh_root, UM_P);
                surface_new = this->measure(elements_root[cpt],M_mesh_root,M_UM_root);

                // for (int i=0; i<M_UM_root.size(); ++i)
                // {
                //     std::cout<<"@@@@@@@@@@@@@@@@@@@@ "<< UM_P[i] <<" --- "<< M_UM_root[i] <<"\n";
                // }

                // std::cout<<"----------------------------[" << M_rank <<"] " <<" UM_P MIN   = "<< *std::min_element(UM_P.begin(),UM_P.end()) <<"\n";
                // std::cout<<"----------------------------[" << M_rank <<"] " <<" UM_P MAX   = "<< *std::max_element(UM_P.begin(),UM_P.end()) <<"\n";

                // std::cout<<"----------------------------[" << M_rank <<"] " <<" UM ROOT MIN= "<< *std::min_element(M_UM_root.begin(),M_UM_root.end()) <<"\n";
                // std::cout<<"----------------------------[" << M_rank <<"] " <<" UM ROOT MAX= "<< *std::max_element(M_UM_root.begin(),M_UM_root.end()) <<"\n";


                //std::cout<<"---------surface    ["<< cpt <<"]= "<< surface <<"\n";
                //std::cout<<"---------surface_new["<< cpt <<"]= "<< surface_new <<"\n";

                ice_surface = old_conc*surface;
                ice_volume = old_thick*surface;
                snow_volume = old_snow_thick*surface;
                ridged_thick_ice_volume = old_h_ridged_thick_ice*surface;

                conc_root[cpt]    = ice_surface/surface_new;
                thick_root[cpt]   = ice_volume/surface_new; // Hold on! Isn't thick the effective thickness?
                snow_thick_root[cpt]   = snow_volume/surface_new;
                h_ridged_thick_ice_root[cpt]   =   ridged_thick_ice_volume/surface_new;

                //if (((M_elements[cpt]).ghosts.size() != 0))
                //std::cout<<"[" << M_rank <<"] " <<"M_thick["<< (M_elements[cpt]).number <<"]= "<< M_thick[cpt] <<"\n";


                if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
                {
                    thin_ice_volume = old_h_thin*surface;
                    thin_snow_volume = old_hs_thin*surface;
                    ridged_thin_ice_volume = old_h_ridged_thin_ice*surface;

                    h_thin_root[cpt]        = thin_ice_volume/surface_new;
                    hs_thin_root[cpt]   = thin_snow_volume/surface_new;
                    h_ridged_thin_ice_root[cpt]    =   ridged_thin_ice_volume/surface_new;
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

                        thick_root[cpt] += ridging_thin_ice ;
                        h_thin_root[cpt] -= ridging_thin_ice ;

                        snow_thick_root[cpt] += ridging_snow_thin_ice ;
                        hs_thin_root[cpt] -= ridging_snow_thin_ice ;

                        h_ridged_thin_ice_root[cpt] += ridging_thin_ice;
                        conc_root[cpt] += ridging_thin_ice/ridge_h;

                        /* upper bounds (only for the concentration) */
                        ridging_thin_ice = ((conc_root[cpt]<1.)?(0.):(h_thin_root[cpt])) ;
                        ridging_snow_thin_ice = ((conc_root[cpt]<1.)?(0.):(hs_thin_root[cpt])) ;

                        snow_thick_root[cpt] += ridging_snow_thin_ice ;
                        hs_thin_root[cpt] -= ridging_snow_thin_ice ;

                        thick_root[cpt] += ridging_thin_ice;
                        h_thin_root[cpt] -= ridging_thin_ice;
                    }
                    /* upper bounds (only for the concentration) */
                    ridging_thick_ice=((conc_root[cpt]<1.)?(0.):(thick_root[cpt]*(conc_root[cpt]-1.)));
                    conc_root[cpt] = ((conc_root[cpt]<1.)?(conc_root[cpt]):(1.)) ;
                }

                /* lower bounds */
                conc_root[cpt] = ((conc_root[cpt]>0.)?(conc_root[cpt] ):(0.)) ;
                thick_root[cpt]        = ((thick_root[cpt]>0.)?(thick_root[cpt]     ):(0.)) ;
                snow_thick_root[cpt]   = ((snow_thick_root[cpt]>0.)?(snow_thick_root[cpt]):(0.)) ;

                if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
                {
                    h_thin_root[cpt]    = ((h_thin_root[cpt]>0.)?(h_thin_root[cpt] ):(0.)) ;
                    hs_thin_root[cpt]   = ((hs_thin_root[cpt]>0.)?(hs_thin_root[cpt]):(0.)) ;
                }
            }

            if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                /* Compute the redistribution of thin ice. */
                /* Returns the change in volume and concentration of thick ice as well as the
                 * change in volume of thin ice. It is called after the
                 * dynamics are done. */

                if(h_thin_root[cpt]>0.)
                {
                    thin_ice_redistribute(h_thin_root[cpt], hs_thin_root[cpt], 0., conc_root[cpt],
                                          tanalpha, rtanalpha, h_thin_max, &new_v_thin, &del_v, &del_c, &del_vs);

                    conc_root[cpt]       += del_c;

                    thick_root[cpt]           += del_v;
                    h_thin_root[cpt]      -= del_v;

                    snow_thick_root[cpt]      += del_vs;
                    hs_thin_root[cpt] -= del_vs;
                }
                else
                {
                    snow_thick_root[cpt] += hs_thin_root[cpt];
                    hs_thin_root[cpt] = 0. ;
                }
            }
        } // loop on elements

        int gsize = M_mesh_root.numTriangles();
        fields_root.resize(nb_var_element*gsize);

        tmp_nb_var=0;
        for (int i=0; i<gsize; ++i)
        {
            tmp_nb_var=0;

            // concentration
            fields_root[nb_var_element*i+tmp_nb_var] = conc_root[i];
            tmp_nb_var++;

            // thickness
            fields_root[nb_var_element*i+tmp_nb_var] = thick_root[i];
            tmp_nb_var++;

            // snow thickness
            fields_root[nb_var_element*i+tmp_nb_var] = snow_thick_root[i];
            tmp_nb_var++;

            // integrated_stress1
            fields_root[nb_var_element*i+tmp_nb_var] = sigma_root[3*i]/**thick[i]*/;
            tmp_nb_var++;

            // integrated_stress2
            fields_root[nb_var_element*i+tmp_nb_var] = sigma_root[3*i+1]/**thick[i]*/;
            tmp_nb_var++;

            // integrated_stress3
            fields_root[nb_var_element*i+tmp_nb_var] = sigma_root[3*i+2]/**thick[i]*/;
            tmp_nb_var++;

            // damage
            fields_root[nb_var_element*i+tmp_nb_var] = damage_root[i];
            tmp_nb_var++;

            // divergence_rate
            fields_root[nb_var_element*i+tmp_nb_var] = divergence_rate_root[i];
            tmp_nb_var++;

            // h_ridged_thin_ice
            fields_root[nb_var_element*i+tmp_nb_var] = h_ridged_thin_ice_root[i];
            tmp_nb_var++;

            // h_ridged_thick_ice
            fields_root[nb_var_element*i+tmp_nb_var] = h_ridged_thick_ice_root[i];
            tmp_nb_var++;

            // thin ice thickness
            fields_root[nb_var_element*i+tmp_nb_var] = h_thin_root[i];
            tmp_nb_var++;

            // snow on thin ice
            fields_root[nb_var_element*i+tmp_nb_var] = hs_thin_root[i];
            tmp_nb_var++;

            if(tmp_nb_var>nb_var_element)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }
    } // rank 0


    std::vector<int> g_sizes_elements = M_sizes_elements_with_ghost;
    //boost::mpi::gather(M_comm, M_num_elements, g_sizes_elements, 0);
    std::vector<int> id_elements;
    int out_size;

    if (M_rank == 0)
    {
        // for (int i=0; i<g_sizes_elements.size(); ++i)
        // {
        //     std::cout<<"GSIZE["<< i <<"]= "<< g_sizes_elements[i] <<"\n";
        // }
        // std::cout<<"Coincidence= "<< M_mesh_root.numTriangles() <<" and "<< M_mesh.numGlobalElements() <<"\n";

        out_size = std::accumulate(g_sizes_elements.begin(),g_sizes_elements.end(),0);
        id_elements.resize(out_size);

        boost::mpi::gatherv(M_comm, M_mesh.trianglesIdWithGhost(), &id_elements[0], g_sizes_elements, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, M_mesh.trianglesIdWithGhost(), 0);
    }

    if (M_rank == 0)
    {
        auto fields_root_nrd = fields_root;
        fields_root.resize(nb_var_element*id_elements.size());

        for (int i=0; i<id_elements.size(); ++i)
        {
            int ri = id_elements[i]-1;

            for (int j=0; j<nb_var_element; ++j)
            {
                fields_root[nb_var_element*i+j] = fields_root_nrd[nb_var_element*ri+j];
            }
        }
    }

    fields_local.resize(nb_var_element*M_num_elements);

    if (M_rank == 0)
    {
        std::for_each(g_sizes_elements.begin(), g_sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });
        boost::mpi::scatterv(M_comm, fields_root, g_sizes_elements, &fields_local[0], 0);
    }
    else
    {
        boost::mpi::scatterv(M_comm, &fields_local[0], nb_var_element*M_num_elements, 0);
    }

    M_comm.barrier();

    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_nb_var=0;

        // concentration
        M_conc[i] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // thickness
        M_thick[i] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // snow thickness
        M_snow_thick[i] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // integrated_stress1
        M_sigma[3*i] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // integrated_stress2
        M_sigma[3*i+1] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // integrated_stress3
        M_sigma[3*i+2] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // damage
        M_damage[i] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // divergence
        M_divergence_rate[i] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // h_ridged_thin_ice
        M_h_ridged_thin_ice[i] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // h_ridged_thick_ice
        M_h_ridged_thick_ice[i] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // thin ice thickness
        M_h_thin[i] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // snow on thin ice
        M_hs_thin[i] = fields_local[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        if(tmp_nb_var>nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }
    }
#endif
}

void
FiniteElement::clear()
{
    M_comm.barrier();

    // delete[] M_topaz_grid.pfindex;

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

} // Nextsim
