/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   meshhandler.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 17 Aug 2022 07:59:06 CEST
 */

#include <meshhandler.hpp>
#include <environment.hpp>
#include <debug.hpp>

namespace Nextsim
{

void
MeshHandler::initOptAndParam()
{
    bamg_verbose = vm["debugging.bamg_verbose"].as<int>();

    //! Defines the export (output) path.
    std::string export_path = vm["output.exporter_path"].as<std::string>(); //! \param export_path (string) Path of the export files
    // Changes directory for outputs if the option "output.exporter_path" is not empty
    fs::path output_path(export_path);

    // Creates the output directory if it does not exist
    if ( (!fs::exists(output_path)) && (M_comm.rank()==0) )
    fs::create_directories(output_path);

    // mesh ordering
    std::vector<std::string> order_opts = {"gmsh", "bamg"};
    M_mesh_ordering = OptionHandler::getAllowedOption(vm, "mesh.ordering", order_opts);
        //! \param M_mesh_ordering (std::string) Mesh ordering ("gmsh" or "bamg")
    LOG(DEBUG) <<"MESH_ORDERING = "<< M_mesh_ordering <<"\n";

    //! Sets the type and format of the mesh and the mesh filename
    const boost::unordered_map<const std::string, setup::MeshType> str2mesh = boost::assign::map_list_of
        ("from_unref", setup::MeshType::FROM_UNREF)
        ("from_split", setup::MeshType::FROM_SPLIT);
    M_mesh_type = OptionHandler::getOptionFromMap(vm, "mesh.type", str2mesh);
        //! \param M_mesh_type (enum) Mesh type (unref or split)
    LOG(DEBUG) <<"MESHTYPE= "<< (int) M_mesh_type <<"\n";

    M_mesh_basename = vm["mesh.filename"].as<std::string>(); //! \param M_mesh_basename (string) Mesh filename
    M_mesh_filename = (boost::format( "%1%/%2%" ) // \param M_mesh_filename (string) Mesh filename (with path)
                       % Environment::nextsimMeshDir().string()
            % M_mesh_basename
            ).str();
    M_partitioned_mesh_filename = (boost::format( "%1%/par%2%%3%" )
            % export_path
            % M_comm.size()
            % M_mesh_basename
            ).str();
    M_mesh_fileformat = vm["mesh.partitioner-fileformat"].as<std::string>(); //! \param M_mesh_fileformat (string) Format of the partitioned mesh file (used if mesh.partitioner-space=="disk")

    //! Sets the type of partitioner and partition space
    const boost::unordered_map<const std::string, mesh::Partitioner> str2partitioner = boost::assign::map_list_of
        ("chaco", mesh::Partitioner::CHACO)
        ("metis", mesh::Partitioner::METIS);
    M_partitioner = OptionHandler::getOptionFromMap(vm, "mesh.partitioner", str2partitioner);
        //! \param M_partitioner (string) Sets the type of partioner (CHACO or METIS)
    LOG(DEBUG) << "MeshPartitioner: "<< (int)M_partitioner<<"\n";

    const boost::unordered_map<const std::string, mesh::PartitionSpace> str2partitionspace = boost::assign::map_list_of
        ("memory", mesh::PartitionSpace::MEMORY)
        ("disk", mesh::PartitionSpace::DISK);

    M_partition_space = OptionHandler::getOptionFromMap(vm, "mesh.partitioner-space", str2partitionspace);
        //! \param M_partition_space (string) Sets the space for partitions (memory or disk)
    LOG(DEBUG) << "MeshPartitionerSpace:" << (int)M_partition_space<<"\n";

    M_use_restart = vm["restart.start_from_restart"].as<bool>(); //! \param M_use_restart (boolean) Option on using starting simulation from a restart file
}

//------------------------------------------------------------------------------------------------------
//! Initialisation of the mesh.
//! Called by the init() function.
void
MeshHandler::initMesh()
{
    this->initBamg();
    this->rootMeshProcessing();
    if (!M_use_restart)
        this->distributedMeshProcessing(true);
}//initMesh

//------------------------------------------------------------------------------------------------------
//! Initializes a Bamg mesh grid.
//! Called by the initMesh() function.
void
MeshHandler::initBamg()
{
    bamgopt = new BamgOpts();
    bamggeom = new BamgGeom();
    bamgmesh = new BamgMesh();

    if ( M_comm.rank() == 0 )
    {
        bamggeom_root = new BamgGeom();
        bamgmesh_root = new BamgMesh();

        bamggeom_previous = new BamgGeom();
        bamgmesh_previous = new BamgMesh();
    }

    bamgopt->Crack             = 0;
    bamgopt->anisomax          = 1e30;
    bamgopt->coeff             = 1;
    bamgopt->cutoff            = 1e-5;
    // bamgopt->err               = 0.01;
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
    bamgopt->splitcorners      = 1; //the Devil!  Changed to 0, original 1 Phil
    bamgopt->geometricalmetric = 0;
    bamgopt->random            = true;
    bamgopt->verbose           = bamg_verbose;

    bamgopt->Check();

}//initBamg

//------------------------------------------------------------------------------------------------------
//! Reads, converts and applies boundary conditions to the mesh.
//! Called by the initMesh() function.
void
MeshHandler::rootMeshProcessing()
{
    if (M_comm.rank() == 0)
    {

        // read the original input mesh
        M_mesh_root.setOrdering(M_mesh_ordering);
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

        LOG(DEBUG) <<"MESH: HMIN= "<< h[0] <<"\n";
        LOG(DEBUG) <<"MESH: HMAX= "<< h[1] <<"\n";
        LOG(DEBUG) <<"MESH: RESOLUTION= "<< M_res_root_mesh <<"\n";

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

            LOG(DEBUG) <<"HMIN MIN= "<< *std::min_element(M_hminVertices.begin(), M_hminVertices.end()) <<"\n";
            LOG(DEBUG) <<"HMIN MAX= "<< *std::max_element(M_hminVertices.begin(), M_hminVertices.end()) <<"\n";
            LOG(DEBUG) <<"HMAX MIN= "<< *std::min_element(M_hmaxVertices.begin(), M_hmaxVertices.end()) <<"\n";
            LOG(DEBUG) <<"HMAX MAX= "<< *std::max_element(M_hmaxVertices.begin(), M_hmaxVertices.end()) <<"\n";

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
            LOG(DEBUG) <<"filename= "<< M_partitioned_mesh_filename <<"\n";

            LOG(DEBUG)<<"------------------------------version       = "<< M_mesh_root.version() <<"\n";
            LOG(DEBUG)<<"------------------------------ordering      = "<< M_mesh_root.ordering() <<"\n";
            LOG(DEBUG)<<"------------------------------format        = "<< M_mesh_fileformat <<"\n";
            LOG(DEBUG)<<"------------------------------space         = "<< vm["mesh.partitioner-space"].as<std::string>() <<"\n";
            LOG(DEBUG)<<"------------------------------partitioner   = "<< vm["mesh.partitioner"].as<std::string>() <<"\n";


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
            LOG(DEBUG) <<"Writing mesh done in "<< chrono.elapsed() <<"s\n";

            // partition the mesh on root process (rank 0)
            chrono.restart();
            M_mesh_root.partition(M_partitioned_mesh_filename,
                    M_partitioner, M_partition_space, M_mesh_fileformat);
            //LOG(DEBUG) <<"Partitioning mesh done in "<< chrono.elapsed() <<"s\n";
            LOG(DEBUG) <<"Partitioning mesh done in "<< chrono.elapsed() <<"s\n";
        }
    }

    boost::mpi::broadcast(M_comm, M_flag_fix, 0);

}//rootMeshProcessing

//------------------------------------------------------------------------------------------------------
//! Distribution of mesh processing for parallel computing.
//! Called by the interpFields(), initMesh() and distributedMeshProcessing() functions.
void
MeshHandler::distributedMeshProcessing(bool start)
{
    M_comm.barrier();

    if (!start)
    {
        M_mesh = mesh_type();
    }

    M_mesh.setOrdering("gmsh");

    LOG(VERBOSE) <<"filename= "<< M_partitioned_mesh_filename <<"\n";

    chrono.restart();
    M_mesh.readFromFile(M_partitioned_mesh_filename, M_mesh_fileformat);
    LOG(DEBUG)<<"-------------------MESHREAD done in "<< chrono.elapsed() <<"s\n";

    chrono.restart();
    delete bamgmesh;
    delete bamggeom;
    bamgmesh = new BamgMesh();
    bamggeom = new BamgGeom();
    BamgConvertMeshx(
                     bamgmesh,bamggeom,
                     &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                     M_mesh.numNodes(), M_mesh.numTriangles());

    LOG(DEBUG)<<"-------------------CREATEBAMG done in "<< chrono.elapsed() <<"s\n";

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_ndof = M_mesh.numGlobalNodes();

    M_local_ndof = M_mesh.numLocalNodesWithoutGhost();
    M_local_ndof_ghost = M_mesh.numLocalNodesWithGhost();

    M_local_nelements = M_mesh.numTrianglesWithoutGhost();
    M_num_nodes = M_local_ndof_ghost;

    chrono.restart();
    this->bcMarkedNodes();
    LOG(DEBUG)<<"-------------------BCMARKER done in "<< chrono.elapsed() <<"s\n";

    chrono.restart();
    this->createGraph();
    LOG(DEBUG)<<"-------------------CREATEGRAPH done in "<< chrono.elapsed() <<"s\n";

    chrono.restart();
    this->gatherSizes();
    LOG(DEBUG)<<"-------------------GATHERSIZE done in "<< chrono.elapsed() <<"s\n";

    chrono.restart();
    this->scatterElementConnectivity();
    LOG(DEBUG)<<"-------------------CONNECTIVITY done in "<< chrono.elapsed() <<"s\n";

    chrono.restart();
    this->initUpdateGhosts();
    LOG(DEBUG)<<"-------------------INITUPDATEGHOSTS done in "<< chrono.elapsed() <<"s\n";

#if 0
    // LOG(DEBUG) << NODES   = "<< M_mesh.numGlobalNodes() << " --- "<< M_local_ndof <<"\n";
    // LOG(DEBUG) << ELEMENTS= "<< M_mesh.numGlobalElements() << " --- "<< M_local_nelements <<"\n";

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

// -------------------------------------------------------------------------------------
//! Creates a new graphmpi_type object.
//! Called by the distributeMeshProcessing() function.
void
MeshHandler::createGraph()
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

//------------------------------------------------------------------------------------------------------
//! Interpolates hminVertices and hmaxVertices onto the current mesh.
//! Called by the distributedMeshProcessing() function.
void
MeshHandler::gatherSizes()
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

    if (M_comm.rank() == 0)
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

    if (M_comm.rank() == 0)
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
    if (M_comm.rank() == 0)
    {
        M_rmap_nodes = M_mesh.mapNodes();
        M_rmap_elements = M_mesh.mapElements();
    }
    // -------------------------------------------------------------
}//gatherSizes

//------------------------------------------------------------------------------------------------------
//! ?? Has to do with the parallel computing.
//! Called by distributedMeshProcessing(), initMesh and  functions.
void
MeshHandler::scatterElementConnectivity()
{
    auto transfer_map_local = M_mesh.transferMapElt();
    std::vector<int> sizes_elements = M_sizes_elements_with_ghost;

    int nb_var_element = 3;
    std::vector<double> connectivity_root;

    if (M_comm.rank() == 0)
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

    if (M_comm.rank() == 0)
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

// -------------------------------------------------------------------------------------
//! Initialise maps to update ghost nodes
//! Called after remeshing and at init by distributedMeshProcessing
void
MeshHandler::initUpdateGhosts()
{
    auto M_transfer_map = M_mesh.transferMap();
    auto M_local_ghost = M_mesh.localGhost();

    std::vector<std::vector<int>> local_ghosts_global_index(M_comm.size());

    for (int i=0; i<M_local_ghost.size(); i++)
    {
        int currentid = M_local_ghost[i];

        // std::cout<<"["<< M_comm.rank() <<"]: Global id= "
        //          << currentid
        //          << " <-----> " << M_transfer_map.left.find(currentid)->second
        //          <<" true ghost= "<< globalNumToprocId(currentid) <<"\n";

        local_ghosts_global_index[this->globalNumToprocId(currentid)].push_back(currentid);
    }

    M_local_ghosts_proc_id.resize(0);
    M_local_ghosts_local_index.resize(M_comm.size());

    for (int i=0; i<M_comm.size(); i++)
    {
        if (local_ghosts_global_index[i].size() != 0)
        {
            M_local_ghosts_proc_id.push_back(i);
        }

        M_local_ghosts_local_index[i].resize(local_ghosts_global_index[i].size());

        // std::cout<<"            ["<< M_rank <<"] <---> "<< i <<"             \n";

        for (int j=0; j<local_ghosts_global_index[i].size(); j++)
        {
            int currentindex = local_ghosts_global_index[i][j];
            M_local_ghosts_local_index[i][j] = M_transfer_map.left.find(currentindex)->second-1;
            // std::cout<<" "<< M_local_ghosts_local_index[i][j] <<"  ";
        }
        // std::cout<<"\n";
    }

    std::vector<std::vector<int>> recipients_proc_id_extended;
    boost::mpi::all_gather(M_comm, M_local_ghosts_proc_id, recipients_proc_id_extended);

    M_recipients_proc_id.resize(0);

    for (int i=0; i<M_comm.size(); i++)
    {
        if (M_comm.rank() == 0)
        {
            LOG(DEBUG) << "["<< i <<"]  ";
            for (int j=0; j<recipients_proc_id_extended[i].size(); j++)
                LOG(DEBUG) << recipients_proc_id_extended[i][j] <<"  ";

            LOG(DEBUG) << "\n";
        }

        for (int j=0; j<recipients_proc_id_extended[i].size(); j++)
        {
            if (recipients_proc_id_extended[i][j] == M_comm.rank())
                M_recipients_proc_id.push_back(i);
        }
    }

    std::vector<std::vector<int>> extract_global_index(M_comm.size());

    for (int const& proc : M_local_ghosts_proc_id)
        M_comm.send(proc, M_comm.rank(), local_ghosts_global_index[proc]);

    for (int const& proc : M_recipients_proc_id)
        M_comm.recv(proc, proc, extract_global_index[proc]);

    M_extract_local_index.resize(M_comm.size());

    for (int i=0; i<extract_global_index.size(); i++)
    {
        int srl = extract_global_index[i].size();
        M_extract_local_index[i].resize(srl);

        for (int j=0; j<extract_global_index[i].size(); j++)
        {
            M_extract_local_index[i][j] = M_transfer_map.left.find(extract_global_index[i][j])->second-1;
        }
    }
} //initUpdateGhosts

// -------------------------------------------------------------------------------------
//! Called by initUpdateGhosts
int
MeshHandler::globalNumToprocId(int global_num)
{
    int cpt = 0;
    for (int i=0; i<M_comm.size(); i++)
    {
        if ((cpt < global_num) && (global_num <= cpt+M_sizes_nodes[i]))
            return i;

        cpt += 2*M_sizes_nodes[i];
    }
} //globalNumToprocId

// -------------------------------------------------------------------------------------
//! Updates the ghost nodes
//! Called by the explicit solver
void
MeshHandler::updateGhosts(std::vector<double>& mesh_nodal_vec)
{
    std::vector<std::vector<double>> extract_local_values(M_comm.size());

    for (int i=0; i<M_extract_local_index.size(); i++)
    {
        int const srl = M_extract_local_index[i].size();
        extract_local_values[i].resize(2*srl);

        for (int j=0; j<M_extract_local_index[i].size(); j++)
        {
            extract_local_values[i][j] = mesh_nodal_vec[M_extract_local_index[i][j]];
            extract_local_values[i][j+srl] = mesh_nodal_vec[M_extract_local_index[i][j]+M_num_nodes];
        }
    }

    std::vector<std::vector<double>> ghost_update_values(M_comm.size());

    for (int const& proc : M_recipients_proc_id)
        M_comm.send(proc, M_comm.rank(), extract_local_values[proc]);

    for (int const& proc : M_local_ghosts_proc_id)
        M_comm.recv(proc, proc, ghost_update_values[proc]);

    for (int i=0; i<M_local_ghosts_local_index.size(); i++)
    {
        for (int j=0; j<M_local_ghosts_local_index[i].size(); j++)
        {
            int const srl = M_local_ghosts_local_index[i].size();
            mesh_nodal_vec[M_local_ghosts_local_index[i][j]] = ghost_update_values[i][j];
            mesh_nodal_vec[M_local_ghosts_local_index[i][j]+M_num_nodes] = ghost_update_values[i][j+srl];
        }
    }
} //updateGhosts

//------------------------------------------------------------------------------------------------------
//! Marks the mesh nodes where boundary conditions should apply.
//! Called by the distributedMeshProcessing() function.
void
MeshHandler::bcMarkedNodes()
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
    if (M_comm.rank() == 0)
    {
        flags_size_root = {(int)M_dirichlet_flags_root.size(), (int)M_neumann_flags_root.size()};
    }

    boost::mpi::broadcast(M_comm, &flags_size_root[0], 2, 0);

    int dir_size = flags_size_root[0];
    int nmn_size = flags_size_root[1];

    std::vector<int> flags_root;
    if (M_comm.rank() == 0)
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

    if (M_comm.rank() != 0)
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

    LOG(DEBUG) << "Dirichlet flags= "<< M_dirichlet_flags.size() <<"\n";
    LOG(DEBUG) << "Neumann flags  = "<< M_neumann_flags.size() <<"\n";


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
//! Adapts the mesh grid.
//! Called by the regrid() function.
void
MeshHandler::adaptMesh()
{
    delete bamgmesh_previous;
    delete bamggeom_previous;
    bamgmesh_previous = bamgmesh_root;
    bamggeom_previous = bamggeom_root;

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

    /* We need bamgopt_modified because Bamgx modifies (at least) bamgopt->hmin
     * and hmax and we don't want anyone to use the modified values! */
    BamgOpts *bamgopt_modified;
    bamgopt_modified = new BamgOpts();
    *bamgopt_modified = *bamgopt;

    /* We also reset bamgmesh_root & bamggeom_root for output
     * (bamgmesh_previous and bamggeom_prevous are pointing to the data, so
     * nothing's lost here). */
    bamgmesh_root = NULL; // probably only needed to convince valgrind that nothing's lost
    bamggeom_root = NULL;
    bamgmesh_root = new BamgMesh();
    bamggeom_root = new BamgGeom();

    chrono.restart();
    Bamgx(bamgmesh_root,bamggeom_root,bamgmesh_previous,bamggeom_previous,bamgopt_modified);
    delete bamgopt_modified;

    LOG(DEBUG) <<"---BAMGMESH done in "<< chrono.elapsed() <<"s\n";

    //! Imports the mesh from bamg, updates the boundary flags and node ID's
    this->importBamg(bamgmesh_root);
    this->updateBoundaryFlags();
    if(bamgopt->KeepVertices)
        this->updateNodeIds();
}//adaptMesh

//------------------------------------------------------------------------------------------------------
//! Interpolates hminVertices and hmaxVertices onto the current mesh.
//! Called by the rootMeshProcessing() and regrid() functions.
void
MeshHandler::interpVertices()
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

// -------------------------------------------------------------------------------------
//! Imports a BAMG mesh grid.
//! Called by the readRestart() and adaptMesh functions.
void
MeshHandler::importBamg(BamgMesh const* bamg_mesh)
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

//------------------------------------------------------------------------------------------------------
//! Updates the boundary flags (Neumann vs Dirichlet) after regriding and mesh adaptation.
//! Called by the adaptMesh() function.
void
MeshHandler::updateBoundaryFlags()
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
//!  Updates the node ID's after regriding and mesh adaptation. Called by the adaptMesh() function.
void
MeshHandler::updateNodeIds()
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
//! Calculates the minimum height (i.e., dimension) of each triangular mesh element.
//! Called by the rootMeshProcessing() function.
std::vector<double>
MeshHandler::hminVertices(mesh_type_root const& mesh, BamgMesh const* bamg_mesh) const
{
    std::vector<double> hmin(bamg_mesh->NodalElementConnectivitySize[0]);

    for (int i=0; i<bamg_mesh->NodalElementConnectivitySize[0]; ++i)
    {
        std::vector<double> measure(bamg_mesh->NodalElementConnectivitySize[1]);
        int k = 0;
        for (int j=0; j<bamg_mesh->NodalElementConnectivitySize[1]; ++j)
        {
            int elt_num = bamg_mesh->NodalElementConnectivity[bamg_mesh->NodalElementConnectivitySize[1]*i+j]-1;

            if ((0 <= elt_num) && (elt_num < mesh.numTriangles()) && (elt_num != NAN))
            {
                measure[k] = mesh.measure(mesh.triangles()[elt_num]);
                k++;
            }
            else
            {
                continue;
            }
        }

        measure.resize(k);
        hmin[i] = std::sqrt(2.)*std::sqrt(*std::min_element(measure.begin(),measure.end()))*0.8;
    }

    return hmin;
}//hminVertices


//------------------------------------------------------------------------------------------------------
//! Returns the maximum height (i.e., dimension) of all mesh elements
//! Called by the rootMeshProcessing() function.
std::vector<double>
MeshHandler::hmaxVertices(mesh_type_root const& mesh, BamgMesh const* bamg_mesh) const
{
    std::vector<double> hmax = this->hminVertices(mesh,bamg_mesh);

    std::for_each(hmax.begin(), hmax.end(), [&](double& f){ f = 1.2*f; });

    return hmax;
}//hmaxVertices


//------------------------------------------------------------------------------------------------------
//! Calculates the maximum and minimum sizes of the vertices of the triangular mesh elements
//! Called by the rootMeshProcessing function.
std::vector<double>
MeshHandler::minMaxSide(mesh_type_root const& mesh) const
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
//! Calculates the mean resolution of the non-structured triangular elements mesh grid.
//! Called by the rootMeshProcessing() function.
double
MeshHandler::resolution(mesh_type_root const& mesh) const
{
    std::vector<double> all_min_measure(mesh.numTriangles());

    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        all_min_measure[cpt] = mesh.measure(*it);
        ++cpt;
    }

    double resol = std::accumulate(all_min_measure.begin(),all_min_measure.end(),0.)/(all_min_measure.size());
    resol = std::pow(resol,0.5);

    return resol;
}//resolution

//------------------------------------------------------------------------------------------------------
//! Calculates the length of the vertices of the triangular mesh elements.
//! Called by the minMaxSides() and minAngles() functions.
std::vector<double>
MeshHandler::sides(element_type const& element, mesh_type const& mesh) const
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
MeshHandler::sides(element_type const& element, mesh_type const& mesh,
                     std::vector<double> const& um, double factor) const
{
    std::vector<double> vertex_0 = mesh.nodes().find(element.indices[0])->second.coords;
    std::vector<double> vertex_1 = mesh.nodes().find(element.indices[1])->second.coords;
    std::vector<double> vertex_2 = mesh.nodes().find(element.indices[2])->second.coords;

    int num_nodes = mesh.numLocalNodesWithGhost();
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


std::vector<double>
MeshHandler::sides(element_type const& element, mesh_type_root const& mesh) const
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
MeshHandler::sides(element_type const& element, mesh_type_root const& mesh,
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
//! Calculates the minimum angle of each triangular mesh element.
//! Called by the minAngle() function.
double
MeshHandler::minAngles(element_type const& element, mesh_type const& mesh) const
{
    std::vector<double> side = this->sides(element,mesh);
    //std::for_each(side.begin(), side.end(), [&](double& f){ f = 1000.*f; });
    std::sort(side.begin(),side.end());
    double minang = std::acos( (std::pow(side[1],2.) + std::pow(side[2],2.) - std::pow(side[0],2.) )/(2*side[1]*side[2]) );
    minang = minang*45.0/std::atan(1.0);

    return minang;
}//minAngles


double
MeshHandler::minAngles(element_type const& element, mesh_type const& mesh,
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
double
MeshHandler::minAngle(mesh_type const& mesh) const
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


double
MeshHandler::minAngle(mesh_type const& mesh, std::vector<double> const& um, double factor, bool root) const
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


#if 0
//------------------------------------------------------------------------------------------------------
//! Performs a re-numbering of the mesh nodes and elements
//! !Does not seem to be used!
void
MeshHandler::rootMeshRenumbering()
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
}//rootMeshRenumbering

//------------------------------------------------------------------------------------------------------
//! Creates a GMSH mesh grid.
//! !Does not seem to be used!
void
MeshHandler::createGMSHMesh(std::string const& geofilename)
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

        LOG(DEBUG) << "[Gmsh::generate] execute '" <<  gmshstr.str() << "'\n";
        auto err = ::system( gmshstr.str().c_str() );
    }
    else
    {
        throw std::runtime_error("Cannot find " + gmshgeofile + "\n");
    }
}//createGMSHMesh
#endif

}



