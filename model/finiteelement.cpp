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
            throw std::runtime_error("No \"coast\" marker in mesh file. Check your input file (" + M_mesh_filename + ")");
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
            M_mesh_filename = (boost::format( "par%1%%2%" ) % M_comm.size() % M_mesh_filename ).str();
            LOG(DEBUG) <<"["<< M_rank <<"] " <<"filename= "<< M_mesh_filename <<"\n";

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

    boost::mpi::broadcast(M_comm, M_flag_fix, 0);

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

    M_nb_regrid = 0;

    M_solver = solver_ptrtype(new solver_type());
    M_matrix = matrix_ptrtype(new matrix_type());
    M_vector = vector_ptrtype(new vector_type());
    M_solution = vector_ptrtype(new vector_type());

    M_reuse_prec = true;

    M_VT.resize(2*M_num_nodes,0.);
    M_VTM.resize(2*M_num_nodes,0.);
    M_VTMM.resize(2*M_num_nodes,0.);

    M_h_thin.assign(M_num_elements,0.);
    M_conc_thin.assign(M_num_elements,0.);
    M_hs_thin.assign(M_num_elements,0.);
    M_tsurf_thin.assign(M_num_elements,0.);

    // stresses
    M_sigma.assign(3*M_num_elements,0.);

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

    M_conc.resize(M_num_elements);
    M_thick.resize(M_num_elements);
    M_damage.resize(M_num_elements);
    M_ridge_ratio.assign(M_num_elements,0.);
    M_snow_thick.resize(M_num_elements);

    M_sst.resize(M_num_elements);
    M_sss.resize(M_num_elements);

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
        M_sigma[3*i]=0.;
        M_sigma[3*i+1]=0.;
        M_sigma[3*i+2]=0.;

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

    M_UT.assign(2*M_num_nodes,0.);

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
}//end initVariables

void
FiniteElement::DataAssimilation()
{
    // Initialise the physical state of the model
    LOG(DEBUG) << "assimilateSlabOcean\n";
    this->assimilateSlabOcean();

    LOG(DEBUG) << "assimilateIce\n";
    this->assimilateIce();
}//DataAssimilation

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

    M_fcor.resize(M_num_elements);

#if 1
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

    M_Cohesion.resize(M_num_elements);
    M_Compressive_strength.resize(M_num_elements);
    M_time_relaxation_damage.resize(M_num_elements,time_relaxation_damage);

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
}

void
FiniteElement::initModelState()
{
    // Initialise the physical state of the model
    LOG(DEBUG) << "initSlabOcean\n";
    this->initSlabOcean();

    LOG(DEBUG) << "initIce\n";
    this->initIce();
}//initModelState

void
FiniteElement::initForcings()
{
    // First we initialise the dataset definitions
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
            std::cout << "invalid ocean forcing"<<"\n";throw std::logic_error("invalid ocean forcing");
    }
    if (M_use_nesting)
    {
        M_nesting_nodes_dataset=DataSet("nesting_nodes",M_num_nodes);
        M_nesting_ocean_elements_dataset=DataSet("nesting_ocean_elements",M_num_elements);
        M_nesting_ice_elements_dataset=DataSet("nesting_ice_elements",M_num_elements);
        M_nesting_dynamics_elements_dataset=DataSet("nesting_dynamics_elements",M_num_elements);
        M_nesting_distance_elements_dataset=DataSet("nesting_distance_elements",M_num_elements);
        M_nesting_distance_nodes_dataset=DataSet("nesting_distance_nodes",M_num_nodes);
    }

    M_ice_topaz_elements_dataset=DataSet("ice_topaz_elements",M_num_elements);

    M_ice_icesat_elements_dataset=DataSet("ice_icesat_elements",M_num_elements);

    M_ice_piomas_elements_dataset=DataSet("ice_piomas_elements",M_num_elements);

    M_ice_amsre_elements_dataset=DataSet("ice_amsre_elements",M_num_elements);

    M_ice_osisaf_elements_dataset=DataSet("ice_osisaf_elements",M_num_elements);

    M_ice_osisaf_type_elements_dataset=DataSet("ice_osisaf_type_elements",M_num_elements);

    M_ice_amsr2_elements_dataset=DataSet("ice_amsr2_elements",M_num_elements);

    M_ice_nic_elements_dataset=DataSet("ice_nic_elements",M_num_elements);

    M_ice_nic_weekly_elements_dataset=DataSet("ice_nic_weekly_elements",M_num_elements);

    M_ice_cs2_smos_elements_dataset=DataSet("ice_cs2_smos_elements",M_num_elements);

    M_ice_smos_elements_dataset=DataSet("ice_smos_elements",M_num_elements);

    // datasets that need to be re-interpolated after regridding
    // - not needed if only used at initialisation, or if not interpolated onto
    // mesh (eg wave datasets are interpolated onto a fixed grid)
    M_datasets_regrid.push_back(&M_atmosphere_nodes_dataset);
    M_datasets_regrid.push_back(&M_atmosphere_elements_dataset);
    M_datasets_regrid.push_back(&M_atmosphere_bis_elements_dataset);
    M_datasets_regrid.push_back(&M_ocean_nodes_dataset);
    M_datasets_regrid.push_back(&M_ocean_elements_dataset);

    // Then we populate the forcing variables
    LOG(DEBUG) <<"Initialize forcingAtmosphere\n";
    this->forcingAtmosphere();

    LOG(DEBUG) <<"Initialize forcingOcean\n";
    this->forcingOcean();
}//initForcings

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

    //loop over ext_data_vec and call check and reload for each:
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
            else if ( (*it)->M_dataset->grid.target_location=="wim_elements" )
            {
                //LOG(DEBUG)<<"in wim_elements: dataset = "<<(*it)->M_dataset->name<<"\n";
                std::cout<<"in wim_elements: dataset = "<<(*it)->M_dataset->name<<"\n";
                //interp to WIM elements
                (*it)->check_and_reload(M_wim.getX(),M_wim.getY(),CRtime);
            }
#endif
            else
            {
                std::cout<<"in dataset = "<<(*it)->M_dataset->name<<"\n";
                std::cout<<"Bad value for dataset->grid.target_location: "<<(*it)->M_dataset->grid.target_location<<"\n";
                std::cout<<"- set to \"mesh_nodes\", or \"mesh_elements\"\n";
#if defined (WAVES)
                std::cout<<"or \"wim_elements\"\n";
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

void
FiniteElement::initOptAndParam()
{
    // log
    const boost::unordered_map<const std::string, LogLevel> str2log = boost::assign::map_list_of
        ("info", INFO)
        ("warning", WARNING)
        ("debug", DEBUG)
        ("error", ERROR);

    M_log_level = str2log.find(vm["debugging.log-level"].as<std::string>())->second;

    // define export path
    M_export_path = Environment::nextsimDir().string() + "/matlab";
    // change directory for outputs if the option "output.exporter_path" is not empty
    if ( ! (vm["output.exporter_path"].as<std::string>()).empty() )
    {
        M_export_path = vm["output.exporter_path"].as<std::string>();
        fs::path path(M_export_path);

        // create the output directory if it does not exist
        if ( (!fs::exists(path)) && (M_comm.rank()==0) )
            fs::create_directories(path);
    }

    nu0 = vm["dynamics.nu0"].as<double>();
    young = vm["dynamics.young"].as<double>();
    rhoi = physical::rhoi;
    rhos = physical::rhos;

    days_in_sec = 24.0*3600.0;
    if (vm["simul.time_init"].as<std::string>() == "")
        throw std::runtime_error("Please provide simul.time_init option (start time)\n");
    else
        time_init = Nextsim::from_date_time_string(vm["simul.time_init"].as<std::string>());
    ptime_step =  days_in_sec/vm["debugging.ptime_per_day"].as<int>();

    time_step = vm["simul.timestep"].as<double>();
    thermo_timestep = vm["simul.thermo_timestep"].as<double>();
    if ( fmod(thermo_timestep,time_step) != 0)
    {
        std::cout << thermo_timestep << " " << time_step << "\n";
        throw std::runtime_error("thermo_timestep is not an integer multiple of time_step");
    }
    // Temporarly disabling super-stepping of the thermodynamcis. The model hangs randomly when it's enabled
    thermo_timestep = time_step;

    output_time_step =  (vm["output.output_per_day"].as<int>()<0) ? time_step : time_step * floor(days_in_sec/vm["output.output_per_day"].as<int>()/time_step);
    mooring_output_time_step =  vm["moorings.output_timestep"].as<double>()*days_in_sec;
    mooring_time_factor = time_step/mooring_output_time_step;
    if ( fmod(mooring_output_time_step,time_step) != 0)
    {
        std::cout << mooring_output_time_step << " " << time_step << "\n";
        throw std::runtime_error("mooring_output_time_step is not an integer multiple of time_step");
    }

    duration = (vm["simul.duration"].as<double>())*days_in_sec;
    restart_time_step =  vm["restart.output_time_step"].as<double>()*days_in_sec;
    M_use_assimilation   = vm["setup.use_assimilation"].as<bool>();
    M_use_restart   = vm["restart.start_from_restart"].as<bool>();
    M_write_restart = vm["restart.write_restart"].as<bool>();
    if ( fmod(restart_time_step,time_step) != 0)
    {
        std::cout << restart_time_step << " " << time_step << "\n";
        throw std::runtime_error("restart_time_step not an integer multiple of time_step");
    }

    ocean_turning_angle_rad = 0.;
    if (vm["dynamics.use_coriolis"].as<bool>())
        ocean_turning_angle_rad = (PI/180.)*vm["dynamics.oceanic_turning_angle"].as<double>();
    ridging_exponent = vm["dynamics.ridging_exponent"].as<double>();

    quad_drag_coef_water = vm["dynamics.quad_drag_coef_water"].as<double>();

    basal_k2 = vm["dynamics.Lemieux_basal_k2"].as<double>();
    basal_u_0 = vm["dynamics.Lemieux_basal_u_0"].as<double>();
    basal_Cb = vm["dynamics.Lemieux_basal_Cb"].as<double>();

    time_relaxation_damage = vm["dynamics.time_relaxation_damage"].as<double>()*days_in_sec;
    deltaT_relaxation_damage = vm["dynamics.deltaT_relaxation_damage"].as<double>();

    h_thin_max = vm["thermo.h_thin_max"].as<double>();
    h_thin_min = vm["thermo.h_thin_min"].as<double>();

    compr_strength = vm["dynamics.compr_strength"].as<double>();
    tract_coef = vm["dynamics.tract_coef"].as<double>();
    // scale_coef is now set after initialising the mesh
    // scale_coef = vm["dynamics.scale_coef"].as<double>();
    alea_factor = vm["dynamics.alea_factor"].as<double>();
    cfix = vm["dynamics.cfix"].as<double>();
    // C_fix    = cfix*scale_coef;          // C_fix;...  : cohesion (mohr-coulomb) in MPa (40000 Pa)
    // C_alea   = alea_factor*C_fix;        // C_alea;... : alea sur la cohesion (Pa)
    tan_phi = vm["dynamics.tan_phi"].as<double>();

    if ( vm["thermo.newice_type"].as<int>() == 4 )
        M_ice_cat_type = setup::IceCategoryType::THIN_ICE;
    else
        M_ice_cat_type = setup::IceCategoryType::CLASSIC;

    const boost::unordered_map<const std::string, setup::ThermoType> str2thermo = boost::assign::map_list_of
        ("zero-layer", setup::ThermoType::ZERO_LAYER)
        ("winton", setup::ThermoType::WINTON);
    M_thermo_type = str2thermo.find(vm["setup.thermo-type"].as<std::string>())->second;
    LOG(DEBUG)<<"ThermoType= "<< (int)M_thermo_type <<"\n";

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

    M_use_nesting= vm["nesting.use_nesting"].as<bool>();

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
    M_ocean_type = str2ocean.find(vm["setup.ocean-type"].as<std::string>())->second;
    LOG(DEBUG) <<"OCEANTYPE= "<< (int)M_ocean_type <<"\n";

    const boost::unordered_map<const std::string, setup::IceType> str2conc = boost::assign::map_list_of
        ("constant", setup::IceType::CONSTANT)
        ("constant_partial", setup::IceType::CONSTANT_PARTIAL)
        ("target", setup::IceType::TARGET)
        ("binary", setup::IceType::BINARY)
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
    M_dynamics_type = str2dynamics.find(vm["setup.dynamics-type"].as<std::string>())->second;
    LOG(DEBUG) <<"DYNAMICSTYPE= "<< (int)M_dynamics_type <<"\n";

    const boost::unordered_map<const std::string, setup::BathymetryType> str2bathymetry = boost::assign::map_list_of
        ("constant", setup::BathymetryType::CONSTANT)
        ("etopo", setup::BathymetryType::ETOPO);
    M_bathymetry_type = str2bathymetry.find(vm["setup.bathymetry-type"].as<std::string>())->second;
    LOG(DEBUG) <<"BATHYMETRYTYPE= "<< (int) M_bathymetry_type <<"\n";

    const boost::unordered_map<const std::string, setup::BasalStressType> str2basal_stress= boost::assign::map_list_of
        ("none", setup::BasalStressType::NONE)
        ("lemieux", setup::BasalStressType::LEMIEUX)
        ("bouillon", setup::BasalStressType::BOUILLON);
    M_basal_stress_type = str2basal_stress.find(vm["setup.basal_stress-type"].as<std::string>())->second;
    LOG(DEBUG) <<"BASALSTRESTYPE= "<< (int) M_basal_stress_type <<"\n";

    // Are we using drifters?
    M_use_iabp_drifters=vm["drifters.use_iabp_drifters"].as<bool>();

    M_use_equallyspaced_drifters=vm["drifters.use_equally_spaced_drifters"].as<bool>();
    M_equallyspaced_drifters_output_time_step=vm["drifters.equallyspaced_drifters_output_time_step"].as<double>();

    M_use_rgps_drifters=vm["drifters.use_rgps_drifters"].as<bool>();
    M_rgps_drifters_output_time_step=vm["drifters.rgps_drifters_output_time_step"].as<double>();

    M_use_osisaf_drifters=vm["drifters.use_osisaf_drifters"].as<bool>();

    M_use_drifters = (M_use_iabp_drifters) || (M_use_osisaf_drifters) || (M_use_equallyspaced_drifters) || (M_use_rgps_drifters);

    M_use_drifters = (M_use_iabp_drifters) || (M_use_osisaf_drifters) || (M_use_equallyspaced_drifters) || (M_use_rgps_drifters);
    LOG(DEBUG) <<"USEDRIFTERS= " << M_use_drifters << " "
        << (int) M_use_iabp_drifters << " "
        << (int) M_use_osisaf_drifters << " "
        << (int) M_use_equallyspaced_drifters << " "
        << (int) M_use_rgps_drifters << "\n";

    // mesh type
    const boost::unordered_map<const std::string, setup::MeshType> str2mesh = boost::assign::map_list_of
        ("from_unref", setup::MeshType::FROM_UNREF)
        ("from_split", setup::MeshType::FROM_SPLIT);
    M_mesh_type = str2mesh.find(vm["mesh.type"].as<std::string>())->second;
    LOG(DEBUG) <<"MESHTYPE= "<< (int) M_mesh_type <<"\n";

    M_mesh_filename = vm["mesh.filename"].as<std::string>();
    M_mesh_fileformat = vm["mesh.partitioner-fileformat"].as<std::string>();
    M_mesh.setOrdering("bamg");

    // Moorings
    M_use_moorings =  vm["moorings.use_moorings"].as<bool>();
    M_moorings_snapshot =  vm["moorings.snapshot"].as<bool>();
    M_moorings_parallel_output =  vm["moorings.parallel_output"].as<bool>();
    const boost::unordered_map<const std::string, GridOutput::fileLength> str2mooringsfl = boost::assign::map_list_of
        ("inf", GridOutput::fileLength::inf)
        ("daily", GridOutput::fileLength::daily)
        ("weekly", GridOutput::fileLength::weekly)
        ("monthly", GridOutput::fileLength::monthly)
        ("yearly", GridOutput::fileLength::yearly);
    M_moorings_file_length = str2mooringsfl.find(vm["moorings.file_length"].as<std::string>())->second;

    const boost::unordered_map<const std::string, mesh::Partitioner> str2partitioner = boost::assign::map_list_of
        ("chaco", mesh::Partitioner::CHACO)
        ("metis", mesh::Partitioner::METIS);
    M_partitioner = str2partitioner.find(vm["mesh.partitioner"].as<std::string>())->second;

    const boost::unordered_map<const std::string, mesh::PartitionSpace> str2partitionspace = boost::assign::map_list_of
        ("memory", mesh::PartitionSpace::MEMORY)
        ("disk", mesh::PartitionSpace::DISK);
    M_partition_space = str2partitionspace.find(vm["mesh.partitioner-space"].as<std::string>())->second;
}//initOptAndParam

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
        hmin[i] = std::sqrt(2.)*std::sqrt(*std::min_element(measure.begin(),measure.end()))*0.8;
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
}

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
}

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
}

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
}

void
FiniteElement::redistributeVariablesIO(std::vector<double> const& out_elt_values, bool thin_ice)
{
    int nb_var_element = M_nb_var_element;
    if (!thin_ice)
    {
        nb_var_element -= 4;
    }

    int tmp_nb_var=0;

    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_nb_var=0;

        // concentration
        M_conc[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // thickness
        M_thick[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // snow thickness
        M_snow_thick[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // integrated_stress1
        M_sigma[3*i] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // integrated_stress2
        M_sigma[3*i+1] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // integrated_stress3
        M_sigma[3*i+2] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // damage
        M_damage[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

        // ridge_ratio
        M_ridge_ratio[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
        tmp_nb_var++;

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
            M_tice[1][i] = out_elt_values[nb_var_element*i+tmp_nb_var];
            tmp_nb_var++;

            M_tice[2][i] = out_elt_values[nb_var_element*i+tmp_nb_var];
            tmp_nb_var++;
        }

        if (thin_ice)
        {
            // thin ice thickness
            M_h_thin[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
            tmp_nb_var++;

            // thin ice concentration
            M_conc_thin[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
            tmp_nb_var++;

            // snow on thin ice
            M_hs_thin[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
            tmp_nb_var++;

            // Ice surface temperature for thin ice
            M_tsurf_thin[i] = out_elt_values[nb_var_element*i+tmp_nb_var];
            tmp_nb_var++;
        }

        if(tmp_nb_var!=nb_var_element)
        {
            throw std::logic_error("tmp_nb_var not equal to nb_var");
        }
    }
}

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

            if(outer_fluxes_area[i]>0)
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
}//advect()

#if 1//advectRoot not used - looks weird too
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
    }

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
}
#endif

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
}

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
}

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
}

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
}

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
}

void
FiniteElement::scatterFieldsElementIO(std::vector<double> const& interp_elt_out, bool thin_ice)
{
    timer["scatter"].first.restart();

    LOG(DEBUG) <<"["<< M_rank <<"]: " <<"----------SCATTER ELEMENT starts\n";

    int nb_var_element = M_nb_var_element;
    if (M_ice_cat_type != setup::IceCategoryType::THIN_ICE)
    {
        nb_var_element -= 4;
    }

    if (vm["output.save_diagnostics"].as<bool>())
    {
        nb_var_element += 7;
    }

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

    this->redistributeVariablesIO(out_elt_values, thin_ice);

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

        // The interpolation with the cavities still needs to be tested on a long run.
        // By default, we then use the non-conservative MeshToMesh interpolation

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
        int global_num_nodes = M_mesh.numGlobalNodes();

        auto field_root_nrd = field_root;

        for (int i=0; i<global_num_nodes; ++i)
        {
            int ri =  M_rmap_nodes[i];

            field_root[i] = field_root_nrd[2*ri];
            field_root[i+global_num_nodes] = field_root_nrd[2*ri+1];
        }
    }
}

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

}

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
            std::cout<<"------------------------------space         = "<< vm["mesh.partitioner-space"].as<std::string>() <<"\n";
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
}//regrid

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

    // Save the old id_node before redefining it
    std::vector<int> old_node_id = M_mesh_root.id();

    // Import the mesh from bamg
    this->importBamg(bamgmesh_root);

    this->updateBoundaryFlags(old_node_id);
}//adaptMesh()

void
FiniteElement::updateBoundaryFlags(std::vector<int> const& old_node_id)
{
    LOG(DEBUG) <<"CLOSED: FLAGS SIZE BEFORE= "<< M_dirichlet_flags_root.size() <<"\n";
    LOG(DEBUG) <<"OPEN  : FLAGS SIZE BEFORE= "<< M_neumann_flags_root.size() <<"\n";

    // update dirichlet nodes
    M_dirichlet_flags_root.resize(0);
    M_neumann_flags_root.resize(0);

    // get the global number of nodes
    int num_nodes = M_mesh_root.numNodes();

#if 0
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
#endif

#if 1
    // We mask out the boundary nodes
    M_mask_root.assign(bamgmesh_root->VerticesSize[0],false) ;
    M_mask_dirichlet_root.assign(bamgmesh_root->VerticesSize[0],false) ;

    for (int vert=0; vert<bamgmesh_root->VerticesOnGeomVertexSize[0]; ++vert)
        M_mask_root[bamgmesh_root->VerticesOnGeomVertex[2*vert]-1] = true; // The factor 2 is because VerticesOnGeomVertex has 2 dimensions in bamg

    // Recompute the node ids
    if(bamgopt->KeepVertices)
    {
        std::vector<int> new_nodes_id=M_mesh_root.id();

        int Boundary_id  = 0;
        int nb_new_nodes = 0;

        // The new id will have values higher than the previous one
        int first_new_node = *std::max_element(old_node_id.begin(),old_node_id.end())+1;

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
                    new_nodes_id[vert] = old_node_id[bamgmesh_root->PreviousNumbering[vert]-1];
                }
            }
        }
        M_mesh_root.setId(new_nodes_id);
    }

    std::vector<int> boundary_flags_root;

    for (int edg=0; edg<bamgmesh_root->EdgesSize[0]; ++edg)
    {
        boundary_flags_root.push_back(bamgmesh_root->Edges[3*edg]);

        if (bamgmesh_root->Edges[3*edg+2] == M_flag_fix)
        {
            M_dirichlet_flags_root.push_back(bamgmesh_root->Edges[3*edg]);
        }
    }

    std::sort(M_dirichlet_flags_root.begin(), M_dirichlet_flags_root.end());

    std::sort(boundary_flags_root.begin(), boundary_flags_root.end());

    std::set_difference(boundary_flags_root.begin(), boundary_flags_root.end(),
                        M_dirichlet_flags_root.begin(), M_dirichlet_flags_root.end(),
                        std::back_inserter(M_neumann_flags_root));

#endif


    M_dirichlet_nodes_root.resize(2*(M_dirichlet_flags_root.size()));
    for (int i=0; i<M_dirichlet_flags_root.size(); ++i)
    {
        M_dirichlet_nodes_root[2*i] = M_dirichlet_flags_root[i];
        M_dirichlet_nodes_root[2*i+1] = M_dirichlet_flags_root[i]+num_nodes;
        M_mask_dirichlet_root[M_dirichlet_flags_root[i]] = true;
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

            /* compute the x and y derivative of g*ssh */
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
            coef_V     = mass_e/time_step;             /* for the inertial term */
            coef_X     = - mass_e*g_ssh_e_x;              /* for the ocean slope */
            coef_Y     = - mass_e*g_ssh_e_y;              /* for the ocean slope */
            coef_sigma = M_thick[cpt]*multiplicator;
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

            /* Skip gohst nodes */
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

void
FiniteElement::calcCohesion()
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
    this->collectVariables(interp_elt_in_local, true);

    // advect
    std::vector<double> interp_elt_out;
    this->advectRoot(interp_elt_in_local, interp_elt_out);

    // redistribute the interpolated values
    this->redistributeVariables(interp_elt_out);

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
         * Ridging scheme and mechanical redistribution
         * After the advection the concentration can be higher than 1, meaning that ridging should have occured.
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
         * Update the internal stress
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
         * Correct the internal stress and the damage
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
}

// Routine for nesting the ice variables () from outer domain
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
}

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
}

// Routine for the 1D thermodynamical model
// No stability dependent drag for now
void
FiniteElement::thermo(double dt)
{
    M_comm.barrier();

    // There is now only one big loop for the thermodynamics so that we can use multithreading.

    // constant variables
    // Set local variable to values defined by options
    double const timeT = vm["thermo.ocean_nudge_timeT"].as<double>();
    double const timeS = vm["thermo.ocean_nudge_timeS"].as<double>();
    double const Qdw_const = vm["ideal_simul.constant_Qdw"].as<double>();
    double const Fdw_const = vm["ideal_simul.constant_Fdw"].as<double>();

    double const ocean_albedo = vm["thermo.albedoW"].as<double>();
    double const drag_ocean_t = vm["thermo.drag_ocean_t"].as<double>();
    double const drag_ocean_q = vm["thermo.drag_ocean_q"].as<double>();

    double const rh0   = 1./vm["thermo.hnull"].as<double>();
    double const rPhiF = 1./vm["thermo.PhiF"].as<double>();

    double const qi = physical::Lf * physical::rhoi;// J m^{-3}
    double const qs = physical::Lf * physical::rhos;// J m^{-3}

    int const newice_type = vm["thermo.newice_type"].as<int>();
    int const melt_type = vm["thermo.melt_type"].as<int>();
    double const PhiM = vm["thermo.PhiM"].as<double>();
    double const PhiF = vm["thermo.PhiF"].as<double>();

    const double aw=6.1121e2, bw=18.729, cw=257.87, dw=227.3;
    const double Aw=7.2e-4, Bw=3.20e-6, Cw=5.9e-10;

    const double alpha=0.62197, beta=0.37803;

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
        {
            Qsw_in=M_Qsw_in[i];
        }
        else
        {
            throw std::logic_error("The function approxSW not yet implemented, you need to initialized M_Qsw_in");
            //Qsw_in=approxSW();
        }

        double mld=( M_mld[i] > vm["ideal_simul.constant_mld"].as<double>() ) ? M_mld[i] : vm["ideal_simul.constant_mld"].as<double>();

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
        {
            sphuma = M_sphuma[i];
        }
        else if ( M_dair.M_initialized )
        {
            double fa     = 1. + Aw + M_mslp[i]*1e-2*( Bw + Cw*M_dair[i]*M_dair[i] );
            double esta   = fa*aw*std::exp( (bw-M_dair[i]/dw)*M_dair[i]/(M_dair[i]+cw) );
            sphuma = alpha*fa*esta/(M_mslp[i]-beta*fa*esta) ;
        }
        else if ( M_mixrat.M_initialized )
        {
            sphuma = M_mixrat[i]/(1.+M_mixrat[i]) ;
        }
        else
        {
            throw std::logic_error("Neither M_dair nor M_mixrat have been initialized. I cannot calculate sphuma.");
        }

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
        // 4) Thickness change of the ice slab (thermoIce0 in matlab)

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
        // 5) Ice growth over open water and lateral melt (thermoOW in matlab)

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
        M_sst[i] = M_sst[i] - dt*( Qio_mean + Qow_mean - Qdw )/(physical::rhow*physical::cpw*mld);

        /* Change in salinity */
        double denominator= ( mld*physical::rhow - del_vi*physical::rhoi - ( del_vs*physical::rhos + (emp-Fdw)*dt) );
        denominator = ( denominator > 1.*physical::rhow ) ? denominator : 1.*physical::rhow;

        double sss_old = M_sss[i];
        M_sss[i] = M_sss[i] + ( (M_sss[i]-physical::si)*physical::rhoi*del_vi + M_sss[i]*(del_vs*physical::rhos + (emp-Fdw)*dt) )
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

        if ( vm["dynamics.use_temperature_dependent_healing"].as<bool>() )
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
                M_time_relaxation_damage[i] = std::max(time_relaxation_damage*deltaT_relaxation_damage/deltaT, dt);
            }
            else
            {
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
        D_delS[i] = (M_sss[i] - sss_old)*physical::rhow*mld/dt;
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
}

// This is Semtner zero layer
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
}

// Initialise everything
void
FiniteElement::init()
{
    // Initialise everything that doesn't depend on the mesh (constants, data set description, and time)

    M_comm.barrier();

    pcpt = 0;
    mesh_adapt_step=0;
    had_remeshed=false;

    this->initOptAndParam();
    M_current_time = time_init /*+ pcpt*time_step/(24*3600.0)*/;

    // Initialise the mesh
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

    // Check the minimum angle of the grid
    double minang = this->minAngle(M_mesh);

    if (minang < vm["numerics.regrid_angle"].as<double>())
    {
        LOG(INFO) <<"invalid regridding angle: should be smaller than the minimal angle in the initial grid\n";
        throw std::logic_error("invalid regridding angle: should be smaller than the minimal angle in the intial grid");
    }

    if ( M_use_restart )
    {
        LOG(DEBUG) <<"Reading restart file\n";

        std::string res_str = vm["restart.restart_string"].as<std::string>();
        if ( !res_str.empty() )
            pcpt = this->readRestart(res_str);
        else
            pcpt = this->readRestart(vm["restart.step_nb"].as<int>());
        M_current_time = time_init + pcpt*time_step/(24*3600.0);

        if(fmod(pcpt*time_step,output_time_step) == 0)
        {
            LOG(DEBUG) <<"export starts\n";
            //this->exportResults((int) pcpt*time_step/output_time_step);
            this->exportResults("restart");
            LOG(DEBUG) <<"export done in " << chrono.elapsed() <<"s\n";
        }
    }
    else
    {
        // Do one regrid to get the mesh right
        //this->regrid(pcpt);

        // Initialise variables
        chrono.restart();
        LOG(DEBUG) <<"Initialize variables\n";
        this->initVariables();
    }

    // Initialise atmospheric and oceanic forcing
    this->initForcings();

    // Initialise bathymetry
    LOG(DEBUG) <<"Initialize bathymetry\n";
    this->initBathymetry();

    // Load data from the datasets we just initialised
    timer["reload"].first.restart();
    this->checkReloadDatasets(M_external_data,
                              M_current_time,
                              "init - time-dependant");
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

    // Open the output file for drifters
    // TODO: Is this the right place to open the file?
    if ( (M_rank == 0) && (M_use_iabp_drifters) )
    {
        // We should tag the file name with the init time in case of a re-start.
        std::stringstream filename;
        filename << M_export_path << "/drifters_out_" << M_current_time << ".txt";

        M_iabp_out.open(filename.str(), std::fstream::out);
        if ( ! M_iabp_out.good() )
            throw std::runtime_error("Cannot write to file: " + filename.str());
    }

    // Initialise the moorings - if requested
    LOG(DEBUG) << "initMoorings\n";
    if ( M_use_moorings )
        this->initMoorings();

    // Initialise drifters - if requested
    LOG(DEBUG) << "initDrifters\n";
    if (M_use_drifters)
        this->initDrifters();
}

// Take one time step
void
FiniteElement::step()
{

    if (vm["debugging.check_fields"].as<bool>())
        // check fields for nans and if thickness is too big
        this->checkFields();

    this->updateDrifterPosition();

#if 1
    if (pcpt == 0)
    {
        // Write results/restart before regrid - useful for debugging
        // NB this only helps if starting from a restart,
        // otherwise regridding has already happened in init(),
        // so won't happen this time step
        // TODO just write restart in init() as well?
        LOG(DEBUG) <<"first export starts\n";
        if (vm["output.datetime_in_filename"].as<bool>())
            this->exportResults(M_current_time);
        else
        {
            int ostep = 0;//need to declare as an int, to make sure it's not interpreted as a double
            this->exportResults(ostep);
        }
        // this->writeRestart(pcpt, 0); // Write a restart before regrid - useful for debugging
        LOG(DEBUG) <<"first export done\n";
    }
#endif


    // remeshing and remapping of the prognostic variables

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

        timer["calcCoriolis"].first.restart();
        this->calcCoriolis();
        if (M_rank == 0)
            std::cout <<"---timer calcCoriolis:             "<< timer["calcCoriolis"].first.elapsed() <<"\n";
    }

    timer["reload"].first.restart();
    this->checkReloadDatasets(M_external_data,
                              M_current_time+time_step/(24*3600.0),
                              "step - time-dependant");
    if (M_rank == 0)
        std::cout <<"---timer check_and_reload:     "<< timer["reload"].first.elapsed() <<"s\n";

    //======================================================================
    // Do the thermodynamics
    //======================================================================
    if ( vm["thermo.use_thermo_forcing"].as<bool>() && (fmod(pcpt*time_step,thermo_timestep) == 0) )
    {
        timer["thermo"].first.restart();
        this->thermo(thermo_timestep);
        if (M_rank == 0)
            std::cout <<"---timer thermo:               "<< timer["thermo"].first.elapsed() <<"s\n";
    }

    //======================================================================
    // Do the nesting of the ice
    //======================================================================
    if( M_use_nesting )
    {
        chrono.restart();
        LOG(DEBUG) <<"nestingIce starts\n";
        this->nestingIce();
        LOG(DEBUG) <<"nestingIce done in "<< chrono.elapsed() <<"s\n";

        //======================================================================
        // Do the nesting of the dynamics
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
        if(had_remeshed && (vm["numerics.regrid_output_flag"].as<bool>()))
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


    //======================================================================
    // Update time, do post-processing and check output
    //======================================================================

    ++pcpt;
    M_current_time = time_init + pcpt*time_step/(24*3600.0);

#if 1
    //======================================================================
    // Output (export and moorings)
    //======================================================================

    if(fmod(pcpt*time_step,output_time_step) == 0)
    {
        chrono.restart();
        LOG(DEBUG) <<"export starts\n";
        if (vm["output.datetime_in_filename"].as<bool>())
            this->exportResults(M_current_time);
        else
        {
            int ostep = pcpt*time_step/output_time_step;//need to declare as an int, to make sure it's not interpreted as a double
            this->exportResults(ostep);
        }
        LOG(DEBUG) <<"export done in " << chrono.elapsed() <<"s\n";
    }

    if ( M_use_moorings )
    {
        this->updateMoorings();
    }

    // check if writing restart each timestep
    bool write_restart = vm["restart.debugging"].as<bool>();
    if(!write_restart)
        // else check if we've reached the restart timestep
        write_restart = (fmod(pcpt*time_step,restart_time_step) == 0);
    if (write_restart)
    {
        std::cout << "Writing restart file after time step " <<  pcpt-1 << "\n";
        if (vm["output.datetime_in_filename"].as<bool>())
            this->writeRestart(pcpt, M_current_time);
        else
        {
            int rstep = pcpt*time_step/restart_time_step;//need to declare as an int, to make sure it's not interpreted as a double
            this->writeRestart(pcpt, rstep );
        }
    }
#endif
#if 0
#endif
 }

// This is the main working function, called from main.cpp (same as perform_simul in the old code)
void
FiniteElement::run()
{
    std::string current_time_system = current_time_local();
    this->init();
    int niter = vm["debugging.maxiteration"].as<int>();

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

        if (pcpt == niter-1)
            is_running = false;

        // **********************************************************************
        // Take one time-step
        // **********************************************************************
        this->step();
    }

    this->exportResults("final");

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
    //bool enable_drifters = (M_use_iabp_drifters) || (M_use_osisaf_drifters) || (M_use_equallyspaced_drifters) || (M_use_rgps_drifters);
    // std::vector<double> M_VT_root;

    if (M_use_drifters)
    {
        this->gatherNodalField(M_UT,M_UT_root);

        this->gatherElementField(M_conc, M_conc_root);
    }

    //std::cout<<"fine"

    // Update the drifters position twice a day, important to keep the same frequency as the IABP data, for the moment
    if ((M_rank == 0) && (M_use_drifters))
    {
        if ( (pcpt==0) || (std::fmod(M_current_time,0.5)==0) )
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
                int ndof_root = M_mesh_root.numNodes();
                std::vector<double> interp_drifter_in(nb_var*ndof_root);

                // Interpolate the velocity
                for (int i=0; i<ndof_root; ++i)
                {
                    interp_drifter_in[nb_var*i]   = M_UT_root[i];
                    interp_drifter_in[nb_var*i+1] = M_UT_root[i+ndof_root];
                }

                double* interp_drifter_out;
                InterpFromMeshToMesh2dx(&interp_drifter_out,
                                        &M_mesh_root.indexTr()[0],&M_mesh_root.coordX()[0],&M_mesh_root.coordY()[0],
                                        M_mesh_root.numNodes(),M_mesh_root.numTriangles(),
                                        &interp_drifter_in[0],
                                        M_mesh_root.numNodes(),nb_var,
                                        &drifter_X[0],&drifter_Y[0],M_iabp_drifters.size(),
                                        true, 0.);


                // Interpolate the concentration - re-use interp_drifter_in
                interp_drifter_in.resize(M_mesh_root.numTriangles());
                for (int i=0; i<M_mesh_root.numTriangles(); ++i)
                {
                    interp_drifter_in[i]   = M_conc_root[i];
                }

                double* interp_drifter_c_out;
                InterpFromMeshToMesh2dx(&interp_drifter_c_out,
                                        &M_mesh_root.indexTr()[0],&M_mesh_root.coordX()[0],&M_mesh_root.coordY()[0],
                                        M_mesh_root.numNodes(),M_mesh_root.numTriangles(),
                                        &interp_drifter_in[0],
                                        M_mesh_root.numTriangles(),1,
                                        &drifter_X[0],&drifter_Y[0],M_iabp_drifters.size(),
                                        true, 0.);

                // Rebuild the M_iabp_drifters map
                double clim = vm["drifters.concentration_limit"].as<double>();
                j=0;
                for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); /* ++it is not allowed here, because we use 'erase' */ )
                {
                    if ( interp_drifter_c_out[j] > clim )
                    {
                        M_iabp_drifters[it->first] = std::array<double,2> {it->second[0]+interp_drifter_out[nb_var*j],
                                                                           it->second[1]+interp_drifter_out[nb_var*j+1]};
                        ++it;
                    }
                    else
                    {
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
            {
                M_equallyspaced_drifters.move(M_mesh_root, M_UT_root);
            }

            if ( M_use_rgps_drifters )
            {
                M_rgps_drifters.move(M_mesh_root, M_UT_root);
            }

            if ( M_use_osisaf_drifters )
            {
                for (auto it=M_osisaf_drifters.begin(); it!=M_osisaf_drifters.end(); it++)
                {
                    it->move(M_mesh_root, M_UT_root);
                }
            }

            for (int i=0; i<M_UT_root.size(); ++i)
            {
                // UT_root
                M_UT_root[i] = 0.;
            }
        }// (pcpt == 0)

        if(pcpt>0)
        {
            if ( M_use_equallyspaced_drifters && fmod(M_current_time,M_equallyspaced_drifters_output_time_step) == 0 )
            {
                M_equallyspaced_drifters.appendNetCDF(M_current_time, M_mesh_root, M_UT_root);
            }

            if ( M_use_rgps_drifters )
            {
                std::string time_str = vm["drifters.RGPS_time_init"].as<std::string>();
                double RGPS_time_init = from_date_time_string(time_str);

                if( (!M_rgps_drifters.isInitialised()) && (M_current_time == RGPS_time_init))
                {
                    this->updateRGPSDrifters();
                }

                if( (M_current_time != RGPS_time_init) && (fmod(M_current_time,M_rgps_drifters_output_time_step) == 0) )
                {
                    if ( M_rgps_drifters.isInitialised() )
                    {
                        M_rgps_drifters.appendNetCDF(M_current_time, M_mesh_root, M_UT_root);
                    }
                }
            }
        }

        if ( (M_use_osisaf_drifters && fmod(M_current_time+0.5,1.) == 0) )
        {
            // OSISAF drift is calculated as a dirfter displacement over 48 hours
            // and they have two sets of drifters in the field at all times.

            // Write out the contents of [1] if it's meaningfull
            if ( M_osisaf_drifters[1].isInitialised() )
            {
                M_osisaf_drifters[1].appendNetCDF(M_current_time, M_mesh_root, M_UT_root);
            }

            // Flip the vector so we move [0] to be [1]
            std::reverse(M_osisaf_drifters.begin(), M_osisaf_drifters.end());

            // Create a new M_drifters instance in [0], with a properly initialised netCDF file
            M_osisaf_drifters[0] = Drifters("data", "ice_drift_nh_polstere-625_multi-oi.nc",
                                            "xc", "yc",
                                            "lat", "lon",
                                            M_mesh_root, M_conc_root, vm["drifters.concentration_limit"].as<double>());

            M_osisaf_drifters[0].initNetCDF(M_export_path+"/OSISAF_", M_current_time);
            M_osisaf_drifters[0].appendNetCDF(M_current_time, M_mesh_root, M_UT_root);
        }
    }

    if ( (M_use_drifters) && ((pcpt == 0) || (std::fmod(M_current_time,0.5)==0)) )
    {
        for (int i=0; i<M_num_nodes; ++i)
	    {
	        // UT
		    M_UT[i] = 0.;
		    M_UT[i+M_num_nodes] = 0.;
	    }
    }
}

// Add to the mean on the mesh
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
}

// Initialise everything w.r.t. the moorings
void
FiniteElement::initMoorings()
{
    bool use_ice_mask = false;

    // Output variables - elements
    std::vector<GridOutput::Variable> elemental_variables(0);

    // Output variables - nodes
    std::vector<GridOutput::Variable> nodal_variables(0);

    // The vectorial variables are (always on the nodes) ...
    std::vector<GridOutput::Vectorial_Variable> vectorial_variables(0);

    std::vector<std::string> names = vm["moorings.variables"].as<std::vector<std::string>>();

    for ( auto it=names.begin(); it!=names.end(); ++it )
    {
        // Element variables
        if ( *it == "conc" )
        {
            GridOutput::Variable conc(GridOutput::variableID::conc);
            elemental_variables.push_back(conc);
        }
        else if ( *it == "thick" )
        {
            GridOutput::Variable thick(GridOutput::variableID::thick);
            elemental_variables.push_back(thick);
        }
        else if ( *it == "snow" )
        {
            GridOutput::Variable snow(GridOutput::variableID::snow);
            elemental_variables.push_back(snow);
        }
        else if ( *it == "tsurf" )
        {
            GridOutput::Variable tsurf(GridOutput::variableID::tsurf);
            elemental_variables.push_back(tsurf);
        }
        else if ( *it == "Qa" )
        {
            GridOutput::Variable Qa(GridOutput::variableID::Qa);
            elemental_variables.push_back(Qa);
        }
        else if ( *it == "Qsw" )
        {
            GridOutput::Variable Qsw(GridOutput::variableID::Qsw);
            elemental_variables.push_back(Qsw);
        }
        else if ( *it == "Qlw" )
        {
            GridOutput::Variable Qlw(GridOutput::variableID::Qlw);
            elemental_variables.push_back(Qlw);
        }
        else if ( *it == "Qsh" )
        {
            GridOutput::Variable Qsh(GridOutput::variableID::Qsh);
            elemental_variables.push_back(Qsh);
        }
        else if ( *it == "Qlh" )
        {
            GridOutput::Variable Qlh(GridOutput::variableID::Qlh);
            elemental_variables.push_back(Qlh);
        }
        else if ( *it == "Qo" )
        {
            GridOutput::Variable Qo(GridOutput::variableID::Qo);
            elemental_variables.push_back(Qo);
        }
        else if ( *it == "delS" )
        {
            GridOutput::Variable delS(GridOutput::variableID::delS);
            elemental_variables.push_back(delS);
        }
        else if ( *it == "conc_thin" & M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            GridOutput::Variable conc_thin(GridOutput::variableID::conc_thin);
            elemental_variables.push_back(conc_thin);
        }
        else if ( *it == "h_thin" & M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            GridOutput::Variable h_thin(GridOutput::variableID::h_thin);
            elemental_variables.push_back(h_thin);
        }
        else if ( *it == "hs_thin" & M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            GridOutput::Variable hs_thin(GridOutput::variableID::hs_thin);
            elemental_variables.push_back(hs_thin);
        }
        // Nodal variables and vectors
        else if ( *it == "velocity_xy" | *it == "velocity_uv" )
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
            if ( *it == "velocity_xy" )
                siuv.east_west_oriented = false;
            else
                siuv.east_west_oriented = true;

            vectorial_variables.push_back(siuv);
        }
        // Error
        else
        {
            std::cout << "Invalid mooring name: " << *it << std::endl;
            std::cout << "Available names are ";
            std::cout << "conc, ";
            std::cout << "thick, ";
            std::cout << "snow, ";
            std::cout << "tsurf, ";
            std::cout << "Qa, ";
            std::cout << "Qsw, ";
            std::cout << "Qlw, ";
            std::cout << "Qsh, ";
            std::cout << "Qlh, ";
            std::cout << "Qo, ";
            std::cout << "delS, ";
            if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
            {
                std::cout << "conc_thin, ";
                std::cout << "h_thin, ";
                std::cout << "hs_thin, ";
            }
            std::cout << "velocity_xy, ";
            std::cout << "velocity_uv";

            throw std::runtime_error("Invalid mooring name");
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
        M_moorings = GridOutput(M_mesh, ncols, nrows, mooring_spacing, xmin, ymin, nodal_variables, elemental_variables, vectorial_variables);
    }
    else if(vm["moorings.grid_type"].as<std::string>()=="from_file")
    {
        // Read the grid in from file
        GridOutput::Grid grid{
            gridFile: Environment::vm()["moorings.grid_file"].as<std::string>(),
            dirname: "data",
            mpp_file: Environment::vm()["mesh.mppfile"].as<std::string>(),
            dimNameX: "y",
            dimNameY: "x",
            latName: "latitude",
            lonName: "longitude"
        };

        // Define the mooring dataset
        M_moorings = GridOutput(M_mesh, grid, nodal_variables, elemental_variables, vectorial_variables);
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
            // shift the timestamp in the file to the centre of the output interval
            output_time = M_current_time;
        else
            output_time = M_current_time - mooring_output_time_step/86400/2;

        std::string filename_root;
        if ( M_moorings_parallel_output )
            filename_root = M_export_path + "/Moorings_" + std::to_string(M_rank);
        else
            filename_root = M_export_path + "/Moorings";

        M_moorings_file = M_moorings.initNetCDF(filename_root, M_moorings_file_length, output_time);
    }

} //initMoorings

void
FiniteElement::updateMoorings()
{
    // If we're taking snapshots the we only call updateMeans before writing to file
    if ( !M_moorings_snapshot )
        this->updateMeans(M_moorings, mooring_time_factor);

    if ( fmod(pcpt*time_step,mooring_output_time_step) == 0 )
    {
        double output_time;
        if ( M_moorings_snapshot )
        {
            // Update the snapshot
            this->updateMeans(M_moorings, 1.);
            // shift the timestamp in the file to the centre of the output interval
            output_time = M_current_time;
        }
        else
        {
            output_time = M_current_time - mooring_output_time_step/86400/2;
        }

        // If it's a new day we check if we need a new file
        double not_used;
        if ( (M_rank==0||M_moorings_parallel_output) && (M_moorings_file_length != GridOutput::fileLength::inf) && (modf(output_time, &not_used) < time_step*86400) )
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

        M_moorings.updateGridMean(M_mesh);

        if ( ! M_moorings_parallel_output )
        {
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

        if ( (M_rank==0) || M_moorings_parallel_output )
            M_moorings.appendNetCDF(M_moorings_file, output_time);

        M_moorings.resetMeshMean(M_mesh);
        M_moorings.resetGridMean();
    }
}//updateMoorings

void
FiniteElement::writeRestart(int pcpt, int step)
{
    std::string tmp = (boost::format( "%1%" ) % step).str();
    this->writeRestart(pcpt,tmp);
}

void
FiniteElement::writeRestart(int pcpt, double date_time)
{
    std::string tmp = to_date_time_string_for_filename(date_time);
    this->writeRestart(pcpt,tmp);
}

void
FiniteElement::writeRestart(int pcpt, std::string step)
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
        std::vector<double> M_tice_root(tice_size*num_elements_root);

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
        std::string directory = Environment::nextsimDir().string() + "/restart";
        // change directory for outputs if the option "output.exporter_path" is not empty
        if ( ! (vm["output.exporter_path"].as<std::string>()).empty() )
            directory = vm["output.exporter_path"].as<std::string>() + "/restart";

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

        exporter.writeField(outbin, M_tice_root, "M_Tice");

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
            throw std::runtime_error("Cannot write to file: " + filename);

        exporter.writeRecord(outrecord);
        outrecord.close();
    }
}

int
FiniteElement::readRestart(int step)
{
    std::string tmp = (boost::format( "%1%" ) % step).str();
    return this->readRestart(tmp);
}

int
FiniteElement::readRestart(std::string step)
{
    Exporter exp_field("double"), exp_mesh("double");
    std::string filename;
    boost::unordered_map<std::string, std::vector<int>>    field_map_int;
    boost::unordered_map<std::string, std::vector<double>> field_map_dbl;

    if (M_rank == 0)
    {
        M_mesh_filename = (boost::format( "par%1%%2%" ) % M_comm.size() % M_mesh_filename ).str();

        // === Read in the mesh restart files ===
        std::string restart_path;
        if ( (vm["restart.input_path"].as<std::string>()).empty() )
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
        LOG(DEBUG)<<"restart file = "<<filename<<"\n";

        std::ifstream meshrecord(filename);
        if ( ! meshrecord.good() )
            throw std::runtime_error("File not found: " + filename);

        exp_mesh.readRecord(meshrecord);
        meshrecord.close();

        // Then onto the data itself
        filename = (boost::format( "%1%/mesh_%2%.bin" )
                    % restart_path
                    % step ).str();

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
                    % step ).str();

        std::ifstream inrecord(filename);
        if ( ! inrecord.good() )
            throw std::runtime_error("File not found: " + filename);

        exp_field.readRecord(inrecord);
        inrecord.close();

        // Then onto the data itself
        filename = (boost::format( "%1%/field_%2%.bin" )
                    % restart_path
                    % step ).str();

        std::fstream inbin(filename, std::ios::binary | std::ios::in);

        if ( ! inbin.good() )
            throw std::runtime_error("File not found: " + filename);

        field_map_int.clear();
        field_map_dbl.clear();
        exp_field.loadFile(inbin, field_map_int, field_map_dbl);
        inbin.close();

        // === Recreate the mesh ===
        // Create bamgmesh and bamggeom
        BamgConvertMeshx(
                         bamgmesh_root,bamggeom_root,
                         &indexTr[0],&coordX[0],&coordY[0],
                         coordX.size(), indexTr.size()/3.
                         );

        // Set and check time
        int pcpt = field_map_int["Misc_int"].at(0);
        std::vector<double> time = field_map_dbl["Time"];
        if (!vm["restart.reset_time_counter"].as<bool>())
        {
            double tmp = time_init + pcpt*time_step/(24*3600.0);
            if ( time[0] != tmp )
            {
                std::cout << "FiniteElement::readRestart: Time and Misc_int[0] (a.k.a pcpt) are inconsistent. \n";
                std::cout << "Time = " << time[0] << " = " << to_date_time_string(time[0])<<"\n";
                std::cout << "time_init + pcpt*time_step/(24*3600.0) = " << tmp << " = " << to_date_time_string(tmp)<<"\n";
                throw std::runtime_error("Inconsistent time information in restart file");
            }
	    }
	    else
	    {
            pcpt = 0;
            if ( time[0] != time_init )
            {
                std::cout << "FiniteElement::readRestart: Restart Time and time_init are inconsistent. \n";
                std::cout << "Time = " << time[0] << " = " << to_date_time_string(time[0])<<"\n";
                std::cout << "time_init = " << time_init << " = " << to_date_time_string(time_init) <<"\n";
                throw std::runtime_error("Inconsistent time information in restart file");
            }
	    }

        // Fix boundaries
        M_flag_fix   = field_map_int["Misc_int"].at(1);

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

        // Import the bamg structs
        this->importBamg(bamgmesh_root);

        // update the boundary flags (to uncomment)
        // this->updateBoundaryFlags();
    }

    // mesh partitioning
    this->partitionMeshRestart();

    // Initialise all the variables to zero
    // done below
    // this->initVariables();

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

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_h_thin.resize(num_elements_root);
            M_conc_thin.resize(num_elements_root);
            M_hs_thin.resize(num_elements_root);
            M_tsurf_thin.resize(num_elements_root);
        }

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
    }

    std::vector<double> interp_elt_out;
    std::vector<double> interp_nd_out;
    this->collectRootRestart(interp_elt_out, interp_nd_out);

    // Initialise all the variables to zero
    this->initVariables();

    this->scatterFieldsElementIO(interp_elt_out, M_ice_cat_type==setup::IceCategoryType::THIN_ICE);

    this->scatterFieldsNode(&interp_nd_out[0]);

    return pcpt;
}//readRestart

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
    int nb_var_element = M_nb_var_element;
    if (M_ice_cat_type != setup::IceCategoryType::THIN_ICE)
    {
        nb_var_element -= 4;
    }


    if (M_rank == 0)
    {
        int num_elements_root = M_mesh_root.numTriangles();
        int tice_size = M_tice.size();

        interp_elt_out.resize(nb_var_element*num_elements_root);

        // std::cout<< "["<< M_rank << "] "<<"num_elements_root                   = "<< num_elements_root <<"\n";

        int tmp_nb_var = 0;

        for (int i=0; i<num_elements_root; ++i)
        {
            tmp_nb_var=0;
            //int ri = M_rmap_elements[i];

            // concentration
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_conc[i];
            tmp_nb_var++;

            // thickness
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_thick[i];
            tmp_nb_var++;

            // snow thickness
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_snow_thick[i];
            tmp_nb_var++;

            // integrated_stress1
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_sigma[3*i];
            tmp_nb_var++;

            // integrated_stress2
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_sigma[3*i+1];
            tmp_nb_var++;

            // integrated_stress3
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_sigma[3*i+2];
            tmp_nb_var++;

            // damage
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_damage[i];
            tmp_nb_var++;

            // ridge_ratio
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_ridge_ratio[i];
            tmp_nb_var++;

            // ridge_ratio
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_random_number[i];
            tmp_nb_var++;

            // sss
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_sss[i];
            tmp_nb_var++;

            // sst
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_sst[i];
            tmp_nb_var++;

            // tice1
            interp_elt_out[nb_var_element*i+tmp_nb_var] = M_tice[0][i];
            tmp_nb_var++;

            if ( interp_elt_out[nb_var_element*i+tmp_nb_var] = M_thermo_type == setup::ThermoType::WINTON )
            {
                // tice2
                interp_elt_out[nb_var_element*i+tmp_nb_var] = M_tice[1][i];
                tmp_nb_var++;

                // tice3
                interp_elt_out[nb_var_element*i+tmp_nb_var] = M_tice[2][i];
                tmp_nb_var++;
            }

            if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                // h_thin
                interp_elt_out[nb_var_element*i+tmp_nb_var] = M_h_thin[i];
                tmp_nb_var++;

                // conc_thin
                interp_elt_out[nb_var_element*i+tmp_nb_var] = M_conc_thin[i];
                tmp_nb_var++;

                // hs_thin
                interp_elt_out[nb_var_element*i+tmp_nb_var] = M_hs_thin[i];
                tmp_nb_var++;

                // tsurf_thin
                interp_elt_out[nb_var_element*i+tmp_nb_var] = M_tsurf_thin[i];
                tmp_nb_var++;
            }

            if(tmp_nb_var>nb_var_element)
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

    // increment M_UT that is used for the drifters
    for (int nd=0; nd<M_UT.size(); ++nd)
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
}

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
}
#endif

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
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(vm["ideal_simul.constant_tair"].as<double>());
            M_external_data.push_back(&M_tair);

            M_mixrat=ExternalData(vm["ideal_simul.constant_mixrat"].as<double>());
            M_external_data.push_back(&M_mixrat);

            M_mslp=ExternalData(vm["ideal_simul.constant_mslp"].as<double>());
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(vm["ideal_simul.constant_Qsw_in"].as<double>());
            M_external_data.push_back(&M_Qsw_in);

            M_Qlw_in=ExternalData(vm["ideal_simul.constant_Qlw_in"].as<double>());
            M_external_data.push_back(&M_Qlw_in);

            M_snowfr=ExternalData(vm["ideal_simul.constant_snowfr"].as<double>());
            M_external_data.push_back(&M_snowfr);

            M_precip=ExternalData(vm["ideal_simul.constant_precip"].as<double>());
            M_external_data.push_back(&M_precip);

            M_dair=ExternalData(vm["ideal_simul.constant_dair"].as<double>());
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

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,air_temperature_correction,false,time_init);
            M_external_data.push_back(&M_tair);

            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,air_temperature_correction,false,time_init);
            M_external_data.push_back(&M_dair);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_tcc);

            // Syl: The following two lines should be removed when approxSW will be implemented in Thermo()
            M_Qsw_in=ExternalData(vm["ideal_simul.constant_Qsw_in"].as<double>());
            M_external_data.push_back(&M_Qsw_in);

            M_precip=ExternalData(0.);
            M_external_data.push_back(&M_precip);
        break;

        case setup::AtmosphereType::EC2:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,air_temperature_correction,false,time_init);
            M_external_data.push_back(&M_tair);

            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,air_temperature_correction,false,time_init);
            M_external_data.push_back(&M_dair);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_Qsw_in);

            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_external_data.push_back(&M_tcc);

            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
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
FiniteElement::forcingNesting()//(double const& u, double const& v)
{
    M_ice_thick=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 0,false,time_init);
    M_external_data.push_back(&M_ice_thick);
    M_ice_conc=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 1,false,time_init);
    M_external_data.push_back(&M_ice_conc);
    M_ice_snow_thick=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 2,false,time_init);
    M_external_data.push_back(&M_ice_snow_thick);
    if ( Environment::vm()["thermo.newice_type"].as<int>() == 4 ) {
        M_ice_h_thin=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 3,false,time_init);
        M_external_data.push_back(&M_ice_h_thin);
        M_ice_conc_thin=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 4,false,time_init);
        M_external_data.push_back(&M_ice_conc_thin);
        M_ice_hs_thin=ExternalData(&M_nesting_ice_elements_dataset, M_mesh, 5,false,time_init);
        M_external_data.push_back(&M_ice_hs_thin);
    }
    M_nesting_dist_elements=ExternalData(&M_nesting_distance_elements_dataset, M_mesh, 0,false,time_init);
    M_external_data.push_back(&M_nesting_dist_elements);
    M_nesting_dist_nodes=ExternalData(&M_nesting_distance_nodes_dataset, M_mesh, 0,false,time_init);
    M_external_data.push_back(&M_nesting_dist_nodes);
    M_nesting_VT1=ExternalData(&M_nesting_nodes_dataset, M_mesh, 0,false,time_init);
    M_external_data.push_back(&M_nesting_VT1);
    M_nesting_VT2=ExternalData(&M_nesting_nodes_dataset, M_mesh, 1,false,time_init);
    M_external_data.push_back(&M_nesting_VT2);
    M_nesting_sigma1=ExternalData(&M_nesting_dynamics_elements_dataset, M_mesh, 0,false,time_init);
    M_external_data.push_back(&M_nesting_sigma1);
    M_nesting_sigma2=ExternalData(&M_nesting_dynamics_elements_dataset, M_mesh, 1,false,time_init);
    M_external_data.push_back(&M_nesting_sigma2);
    M_nesting_sigma3=ExternalData(&M_nesting_dynamics_elements_dataset, M_mesh, 2,false,time_init);
    M_external_data.push_back(&M_nesting_sigma3);
    M_nesting_damage=ExternalData(&M_nesting_dynamics_elements_dataset, M_mesh, 3,false,time_init);
    M_external_data.push_back(&M_nesting_damage);
    M_nesting_ridge_ratio=ExternalData(&M_nesting_dynamics_elements_dataset, M_mesh, 4,false,time_init);
    M_external_data.push_back(&M_nesting_ridge_ratio);
}

void
FiniteElement::forcingOcean()//(double const& u, double const& v)
{

    if(M_use_nesting)
    {
        if(M_use_ocean_nesting)
        {
            M_ocean_temp=ExternalData(&M_nesting_ocean_elements_dataset, M_mesh, 0,false,time_init);
            M_external_data.push_back(&M_ocean_temp);
            M_ocean_salt=ExternalData(&M_nesting_ocean_elements_dataset, M_mesh, 1,false,time_init);
            M_external_data.push_back(&M_ocean_salt);
        }
    }

    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
            M_ocean=ExternalData(
                vm["ideal_simul.constant_ocean_u"].as<double>(),
                vm["ideal_simul.constant_ocean_v"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ocean);

            M_ssh=ExternalData(vm["ideal_simul.constant_ssh"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ssh);

            if ( (!M_use_nesting) || ( (M_use_nesting) && (!M_use_ocean_nesting) ) )
            {
                M_ocean_temp=ExternalData(physical::ocean_freezing_temp);
                M_external_data.push_back(&M_ocean_temp);

                M_ocean_salt=ExternalData(physical::ocean_freezing_temp/physical::mu);
                M_external_data.push_back(&M_ocean_salt);
            }

            M_mld=ExternalData(vm["ideal_simul.constant_mld"].as<double>());
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

            if ( (!M_use_nesting) || ( (M_use_nesting) && (!M_use_ocean_nesting) ) )
            {
                M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
                M_external_data.push_back(&M_ocean_temp);

                M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
                M_external_data.push_back(&M_ocean_salt);
            }

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            M_external_data.push_back(&M_mld);
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

            if ( (!M_use_nesting) || ( (M_use_nesting) && (!M_use_ocean_nesting) ) )
            {
                M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
                M_external_data.push_back(&M_ocean_temp);

                M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
                M_external_data.push_back(&M_ocean_salt);
            }

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            M_external_data.push_back(&M_mld);
            break;

        case setup::OceanType::TOPAZR_atrest:
            M_ocean=ExternalData(
                vm["ideal_simul.constant_ocean_u"].as<double>(),
                vm["ideal_simul.constant_ocean_v"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ocean);

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ssh);

            if ( (!M_use_nesting) || ( (M_use_nesting) && (!M_use_ocean_nesting) ) )
            {
                M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
                M_external_data.push_back(&M_ocean_temp);

                M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
                M_external_data.push_back(&M_ocean_salt);
            }

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            M_external_data.push_back(&M_mld);
            // SYL: there was a capping of the mld at minimum vm["ideal_simul.constant_mld"].as<double>()
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

    double hi, conc_tot, weight_conc;
    // Consistency check for the slab ocean + initialization of the ice temperature (especially for Winton)
    for ( int i=0; i<M_num_elements; i++ )
    {
        hi=0.;
        if(M_conc[i]>0.)
            hi = M_thick[i]/M_conc[i];

        if ( (M_conc[i] < physical::cmin) || (hi < physical::hmin) )
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }

        hi = 0.;
        if(M_conc_thin[i]>0.)
            hi = M_h_thin[i]/M_conc_thin[i];

        if ( (M_conc_thin[i] < physical::cmin) || (hi < physical::hmin) )
        {
            M_conc_thin[i]=0.;
            M_h_thin[i]=0.;
            M_hs_thin[i]=0.;
        }

        conc_tot=M_conc[i]+M_conc_thin[i];
        weight_conc=std::min(1.,conc_tot*100.);
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
}

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
        weight_conc=std::min(1.,conc_tot*100.);
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
}

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

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_init_conc);
    external_data_tmp.push_back(&M_init_thick);
    external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(external_data_tmp,time_init,"init - TOPAZ ice");
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
    this->checkReloadDatasets(external_data_tmp,time_init,"init - TOPAZ ice");
    external_data_tmp.resize(0);

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

void
FiniteElement::topazForecastIce()
{
    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    external_data M_init_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_init_conc);
    external_data_tmp.push_back(&M_init_thick);
    external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(external_data_tmp,time_init,"init - TOPAZ ice forecast");
    external_data_tmp.resize(0);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
        // - TOPAZ puts very small values instead of 0.
        // - uses absolute thickness not effective thickness
        tmp_var=std::min(1.,M_init_conc[i]);
        M_conc[i] = (tmp_var>1e-14) ? tmp_var*M_conc[i] : 0.;
        tmp_var=M_init_thick[i];
        M_thick[i] = (tmp_var>1e-14) ? tmp_var*M_conc[i] : 0.;
        tmp_var=M_init_snow_thick[i];
        M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

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
}//topazForecastIce

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
    this->checkReloadDatasets(external_data_tmp,time_init,
            "init - TOPAZ ice forecast + AMSR2");
    external_data_tmp.resize(0);

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

void
FiniteElement::concBinsNic(double &thin_conc_obs_min,double &thin_conc_obs_max,double ci,bool use_weekly_nic)
{

    double alpha_up = 0.;
    if(!use_weekly_nic)
    {
        if(ci<=1.e-3)
        {
            thin_conc_obs_min=0.;
            thin_conc_obs_max=0.1;
        }
        if(ci<=0.45)
        {
            // CT18
            // .45: .1 - .8
            thin_conc_obs_min=0.1;
            thin_conc_obs_max=0.8;
        }
        else if(ci<=0.9)
        {
            // CT81
            // .9: .8 - 1.
            thin_conc_obs_min=0.8;
            thin_conc_obs_max=1.;
        }
        else if(ci<=1.)
        {
            thin_conc_obs_min=0.8;
            thin_conc_obs_max=1.;
        }
        else // should not happen
        {
            thin_conc_obs_min=1.;
            thin_conc_obs_max=1.;
        }
    }
    else
    {
        if(ci<=0.05) // CT01
        {
            // 0 - .1
            thin_conc_obs_min=0.;
            thin_conc_obs_max=0.1;
        }
        else if(ci<=0.10) // CT02
        {
            // 0 - .2
            thin_conc_obs_min=0.;
            thin_conc_obs_max=0.2;
        }
        else if(ci<=0.25) // CT14, CT13, CT24
        {
            // .1 - .4
            thin_conc_obs_min=0.1;
            thin_conc_obs_max=0.4;
        }
        else if(ci<=0.50) // CT46
        {
            // .4 - .6
            thin_conc_obs_min=0.4;
            thin_conc_obs_max=0.6;
        }
        else if(ci<=0.70) //CT68
        {
            // .6 - .8
            thin_conc_obs_min=0.6;
            thin_conc_obs_max=0.8;
        }
        else if(ci<=0.90) // CT81
        {
            // .8 - 1.
            thin_conc_obs_min=0.8;
            thin_conc_obs_max=1.0;
        }
        else if(ci<=1.) // CT92
        {
            // .9 - 1.
            thin_conc_obs_min=0.9;
            thin_conc_obs_max=1.0;
        }
        else // should not happen
        {
            thin_conc_obs_min=1.;
            thin_conc_obs_max=1.;
        }
    }
}//concBinsNic

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

    this->checkReloadDatasets(external_data_tmp, time_init-0.5,
            "assimilate - NIC");

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
    DataSet dist2coast_elements_dataset=DataSet("dist2coast_elements",M_num_elements);
    external_data M_dist2coast = ExternalData(&dist2coast_elements_dataset,M_mesh,0,false,time_init);

    external_data_vec external_data_tmp;
    external_data_tmp.push_back(&M_osisaf_conc);
    external_data_tmp.push_back(&M_osisaf_type);
    external_data_tmp.push_back(&M_amsr2_conc);
    external_data_tmp.push_back(&M_dist2coast);
    this->checkReloadDatasets(external_data_tmp,time_init-0.5,
            "assim - OSISAF - AMSR2 -dist2coast");

    external_data_tmp.resize(0);
    external_data_tmp.push_back(&M_topaz_conc);
    external_data_tmp.push_back(&M_topaz_thick);
    external_data_tmp.push_back(&M_topaz_snow_thick);
    this->checkReloadDatasets(external_data_tmp,time_init,
            "assim - TOPAZ ice forecast");
    external_data_tmp.resize(0);

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

        if((M_amsr2_conc[i]<M_conc[i]) // AMSR2 is higher resolution and sees small opening that would not be see in OSISAF
                && (M_amsr2_conc[i]>.15))
            M_conc[i]=M_amsr2_conc[i];

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
}
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
    this->checkReloadDatasets(external_data_tmp,time_init-0.5,
            "init - OSISAF - AMSR2");

    external_data_tmp.resize(0);
    external_data_tmp.push_back(&M_topaz_conc);
    external_data_tmp.push_back(&M_topaz_thick);
    external_data_tmp.push_back(&M_topaz_snow_thick);
    this->checkReloadDatasets(external_data_tmp,time_init,
            "init - TOPAZ ice forecast");
    external_data_tmp.resize(0);

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

            if(M_mesh_filename.find("kara") != std::string::npos)
            {
                LOG(DEBUG) <<"Type information is not used for the kara mesh, "
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

    this->checkReloadDatasets(external_data_tmp, time_init-0.5,
            "init - OSISAF - AMSR2 - NIC");

    external_data_tmp.resize(0);
    external_data_tmp.push_back(&M_topaz_conc);
    external_data_tmp.push_back(&M_topaz_thick);
    external_data_tmp.push_back(&M_topaz_snow_thick);
    this->checkReloadDatasets(external_data_tmp, time_init,
            "init - TOPAZ ice forecast");
    external_data_tmp.resize(0);

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

            if(M_mesh_filename.find("kara") != std::string::npos)
            {
                LOG(DEBUG) <<"Type information is not used for the kara mesh, we assume there is only FYI\n";
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
    this->checkReloadDatasets(external_data_tmp,time_init,
            "init ice - PIOMAS");
    external_data_tmp.resize(0);

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
    this->checkReloadDatasets(external_data_tmp,time_init,
            "init ice - TOPAZ + AMSR-E");
    external_data_tmp.resize(0);

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
    this->checkReloadDatasets(external_data_tmp,time_init,
            "init ice - TOPAZ + AMSR2");
    external_data_tmp.resize(0);

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
    this->checkReloadDatasets(external_data_tmp,time_init,
            "init ice - CS2 + SMOS");
    external_data_tmp.resize(0);

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
}//cs2SmosIce

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
    this->checkReloadDatasets(external_data_tmp,time_init,
            "init ice - CS2 + SMOS + AMSR2");
    external_data_tmp.resize(0);

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
}//cs2SmosAmsr2Ice

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
    this->checkReloadDatasets(external_data_tmp,time_init,
            "init ice - SMOS");
    external_data_tmp.resize(0);

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
FiniteElement::initDrifters()
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
FiniteElement::calcCoriolis()
{
    // Interpolation of the latitude
    std::vector<double> lat = M_mesh.meanLat();

    for (int i=0; i<M_fcor.size(); ++i)
    {
        if (vm["dynamics.use_coriolis"].as<bool>())
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
FiniteElement::initBathymetry()//(double const& u, double const& v)
{
    // This always needs to be done, regardless of which bathymetry type we
    // have as the object M_bathymetry_elements_dataset must be initialised. But
    // if we use CONSTANT then we don't put any data into the object.
    M_bathymetry_elements_dataset=DataSet("etopo_elements",M_num_elements);//M_num_nodes);

    switch (M_bathymetry_type)
    {
        case setup::BathymetryType::CONSTANT:
            M_element_depth=ExternalData(vm["ideal_simul.constant_bathymetry"].as<double>());
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
    M_datasets_regrid.push_back(&M_bathymetry_elements_dataset);//this needs to be reloaded if we are regridding
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
FiniteElement::outputDrifter(std::fstream& drifters_out)
{
    if (M_rank == 0)
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
        int ndof_root = M_mesh_root.numNodes();
        std::vector<double> interp_drifter_in(nb_var*ndof_root);

        // Interpolate the velocity
        for (int i=0; i<ndof_root; ++i)
        {
            interp_drifter_in[nb_var*i]   = M_UT_root[i];
            interp_drifter_in[nb_var*i+1] = M_UT_root[i+ndof_root];
        }

        double* interp_drifter_out;
        InterpFromMeshToMesh2dx(&interp_drifter_out,
                                &M_mesh_root.indexTr()[0],&M_mesh_root.coordX()[0],&M_mesh_root.coordY()[0],
                                M_mesh_root.numNodes(),M_mesh_root.numTriangles(),
                                &interp_drifter_in[0],
                                M_mesh_root.numNodes(),nb_var,
                                &drifter_X[0],&drifter_Y[0],M_iabp_drifters.size(),
                                true, 0.);

        // Loop over the map and output
        j=0;
        boost::gregorian::date           date = Nextsim::parse_date( M_current_time );
        boost::posix_time::time_duration time = Nextsim::parse_time( M_current_time );
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
}

// Add the buoys that have been put into the ice and remove dead ones
void
FiniteElement::updateIABPDrifter()
{
    if (M_rank == 0)
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
        double time = M_current_time;
        std::vector<int> keepers;
        while ( time == M_current_time )
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
}

// Initialise by reading all the data from '79 up to time_init
// This is too slow, but only happens once so I won't try to optimise that for now
void
FiniteElement::initIABPDrifter()
{
    if (M_rank == 0)
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
}

void
FiniteElement::equallySpacedDrifter()
{
    if (M_rank == 0)
    {
        M_equallyspaced_drifters = Drifters(1e3*vm["drifters.spacing"].as<double>(), M_mesh_root, M_conc_root, vm["drifters.concentration_limit"].as<double>());
        M_equallyspaced_drifters.initNetCDF(M_export_path+"/Drifters_", M_current_time);
        M_equallyspaced_drifters.appendNetCDF(M_current_time, M_mesh_root, M_UT_root);
    }
}

void
FiniteElement::initRGPSDrifters()
{
    if (M_rank == 0)
    {
        M_rgps_drifters = Drifters();
    }
}

void
FiniteElement::updateRGPSDrifters()
{
    if (M_rank == 0)
    {
        std::string time_str = vm["drifters.RGPS_time_init"].as<std::string>();
        double RGPS_time_init = Nextsim::from_date_time_string(time_str);

        std::string filename = Environment::nextsimDir().string() + "/data/RGPS_" + time_str + ".txt";
        M_rgps_drifters = Drifters(filename, M_mesh_root, M_conc_root, vm["drifters.concentration_limit"].as<double>(),RGPS_time_init);

        M_rgps_drifters.initNetCDF(M_export_path+"/RGPS_Drifters_", M_current_time);
        M_rgps_drifters.appendNetCDF(M_current_time, M_mesh_root, M_UT_root);
    }
}

void
FiniteElement::initOSISAFDrifters()
{
    if (M_rank == 0)
    {
        M_osisaf_drifters.resize(2);
    }
}

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
    //M_mesh_root = mesh_type_root(mesh_nodes,mesh_edges,mesh_triangles);
    M_mesh_root.update(mesh_nodes,mesh_triangles);
    //M_mesh.writeTofile("out.msh");

    LOG(DEBUG) <<"\n";
    LOG(DEBUG) <<"Current  NumNodes      = "<< M_mesh_root.numNodes() <<"\n";
    LOG(DEBUG) <<"Current  NumTriangles  = "<< M_mesh_root.numTriangles() <<"\n";
    LOG(DEBUG) <<"\n";
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
FiniteElement::exportResults(int step, bool export_mesh, bool export_fields, bool apply_displacement)
{
    //define name_str from step
    std::string name_str    = (boost::format( "%1%" )
                               % step ).str();

    this->exportResults(name_str, export_mesh, export_fields, apply_displacement);
}

void
FiniteElement::exportResults(double date_time, bool export_mesh, bool export_fields, bool apply_displacement)
{
    //define name_str from date_time
    std::string name_str = to_date_time_string_for_filename(date_time);
    this->exportResults(name_str, export_mesh, export_fields, apply_displacement);
}

void
FiniteElement::exportResults(std::string const& name_str, bool export_mesh, bool export_fields, bool apply_displacement)
{
    //define filenames from iname_str
    std::string meshfile    = (boost::format( "%1%/mesh_%2%" )
                               % M_export_path
                               % name_str ).str();

    std::string fieldfile   = (boost::format( "%1%/field_%2%" )
                               % M_export_path
                               % name_str ).str();

    std::vector<std::string> filenames = {meshfile,fieldfile};
    this->exportResults(filenames, export_mesh, export_fields, apply_displacement);
}

void
FiniteElement::exportResults(std::vector<std::string> const& filenames, bool export_mesh, bool export_fields, bool apply_displacement)
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
    const char* senv = ::getenv(envname.c_str());
    if ( senv == NULL )
        senv = "NULL";
    return std::string(senv);
}

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
        logfile << std::setw(log_width) << std::left << "NEXTSIMDIR "  << getEnv("NEXTSIMDIR") <<"\n";
        logfile << std::setw(log_width) << std::left << "SIMDATADIR "  << getEnv("SIMDATADIR") <<"\n";
        logfile << std::setw(log_width) << std::left << "SIMFORECASTDIR "  << getEnv("SIMFORECASTDIR") <<"\n";
        logfile << std::setw(log_width) << std::left << "PETSC_DIR "  << getEnv("PETSC_DIR") <<"\n";
        logfile << std::setw(log_width) << std::left << "BOOST_INCDIR "  << getEnv("BOOST_INCDIR") <<"\n";
        logfile << std::setw(log_width) << std::left << "BOOST_LIBDIR "  << getEnv("BOOST_LIBDIR") <<"\n";
        logfile << std::setw(log_width) << std::left << "GMSH_DIR "  << getEnv("GMSH_DIR") <<"\n";
        logfile << std::setw(log_width) << std::left << "NETCDF_DIR "  << getEnv("NETCDF_DIR") <<"\n";
        logfile << std::setw(log_width) << std::left << "OPENMPI_LIB_DIR "  << getEnv("OPENMPI_LIB_DIR") <<"\n";
        logfile << std::setw(log_width) << std::left << "OPENMPI_INCLUDE_DIR "  << getEnv("OPENMPI_INCLUDE_DIR") <<"\n";

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
        std::string configfile = (boost::format( "%1%/%2%/%3%" )
                                  % Environment::nextsimDir().string()
                                  % "data"
                                  % vm["mesh.mppfile"].as<std::string>()
                                  ).str();

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
}

} // Nextsim
