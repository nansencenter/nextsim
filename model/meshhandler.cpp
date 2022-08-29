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

}



