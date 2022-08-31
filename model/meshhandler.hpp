/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   meshhandler.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 17 Aug 2022 07:59:06 CEST
 */

#ifndef __MeshHandler_HPP
#define __MeshHandler_HPP 1

#include <boost/program_options.hpp>
#include <gmshmesh.hpp>
#include <gmshmeshseq.hpp>
#include <BamgConvertMeshx.h>
#include <BamgTriangulatex.h>
#include <Bamgx.h>
#include <enums.hpp>
#include <optionhandler.hpp>
#include <InterpFromMeshToMesh2dx.h>
#include <graphcsrmpi.hpp>
#include <boost/serialization/vector.hpp>

namespace Nextsim
{

    class MeshHandler
    {
        // Types
        protected:
            typedef typename GmshMesh::point_type point_type;
            typedef typename GmshMesh::element_type element_type;
            typedef GmshMeshSeq mesh_type_root;
            typedef GmshMesh mesh_type;

            typedef GraphCSRMPI graphmpi_type;

        // Methods
        public:
            MeshHandler(Communicator const& comm)
                : vm(Environment::vm()),
                    M_mesh(mesh_type(comm)),
                    M_log_level(Environment::logLevel()),
                    M_log_all(Environment::logAll()),
                    M_comm(comm)
            {
                this->initOptAndParam();
            };

            void initMesh();

        protected:
            void initBamg();

            void rootMeshProcessing();

            void updateBoundaryFlags();

            void importBamg(BamgMesh const* bamg_mesh);

            void interpVertices();

            void adaptMesh();

            std::vector<double> sides(element_type const& element, mesh_type const& mesh) const;
            std::vector<double> sides(element_type const& element, mesh_type const& mesh,
                                      std::vector<double> const& um, double factor = 1.) const;

            std::vector<double> sides(element_type const& element, mesh_type_root const& mesh) const;
            std::vector<double> sides(element_type const& element, mesh_type_root const& mesh,
                                      std::vector<double> const& um, double factor = 1.) const;

            double minAngles(element_type const& element, mesh_type const& mesh) const;
            double minAngles(element_type const& element, mesh_type const& mesh,
                             std::vector<double> const& um, double factor) const;

            double minAngle(mesh_type const& mesh) const;
            double minAngle(mesh_type const& mesh, std::vector<double> const& um, double factor, bool root = false) const;

            void distributedMeshProcessing(bool start = false);

            void bcMarkedNodes();
            void createGraph();//(BamgMesh const* bamg_mesh);
            void gatherSizes();
            void scatterElementConnectivity();
            void initUpdateGhosts();
            int globalNumToprocId(int global_num);
            void updateGhosts(std::vector<double>& mesh_nodal_vec);

        private:
            void initOptAndParam();

            std::vector<double> minMaxSide(mesh_type_root const& mesh) const;
            double resolution(mesh_type_root const& mesh) const;

            std::vector<double> hminVertices(mesh_type_root const& mesh, BamgMesh const* bamg_mesh) const;
            std::vector<double> hmaxVertices(mesh_type_root const& mesh, BamgMesh const* bamg_mesh) const;

            void updateNodeIds();

        // Variables
        public:

        protected:
            mesh_type M_mesh;
            mesh_type M_mesh_previous;
            mesh_type_root M_mesh_previous_root;

            std::vector<int> M_neumann_flags_root;
            std::vector<int> M_neumann_nodes_root;
            std::vector<int> M_dirichlet_nodes_root;

            // local masks
            std::vector<bool> M_mask;
            std::vector<bool> M_mask_dirichlet;

            // global masks on root
            std::vector<bool> M_mask_root;
            std::vector<bool> M_mask_dirichlet_root;

            /* These need to be pointers because the bamg routines expect
             * pointers. We get memory leaks if we pass them references
             * instead. */
            // Global bamg options
            BamgOpts *bamgopt;

            // Per-processor mesh
            BamgMesh *bamgmesh;
            BamgGeom *bamggeom;

            // Root mesh
            BamgMesh *bamgmesh_root;
            BamgGeom *bamggeom_root;

            // Old root mesh
            BamgMesh *bamgmesh_previous;
            BamgGeom *bamggeom_previous;

            mesh_type_root M_mesh_root;
            mesh_type_root M_mesh_init_root;

            std::string M_mesh_ordering;
            setup::MeshType M_mesh_type;

            std::string M_mesh_basename;
            std::string M_mesh_filename;
            std::string M_partitioned_mesh_filename;
            std::string M_mesh_fileformat;

            mesh::Partitioner M_partitioner;
            mesh::PartitionSpace M_partition_space;

            int M_flag_fix;
            std::vector<int> M_dirichlet_flags_root;
            double M_res_root_mesh;

            std::vector<double> M_hminVertices;
            std::vector<double> M_hmaxVertices;

            int M_num_elements;
            int M_ndof;
            int M_local_ndof;
            int M_local_ndof_ghost;
            int M_local_nelements;

            int M_num_nodes;

            std::vector<element_type> M_elements;
            std::map<int, point_type > M_nodes;

            std::vector<int> M_dirichlet_flags;
            std::vector<int> M_dirichlet_nodes;
            std::vector<int> M_neumann_flags;
            std::vector<int> M_neumann_nodes;

            std::vector<int> M_sizes_nodes;
            std::vector<int> M_sizes_nodes_with_ghost;
            std::vector<int> M_sizes_elements;
            std::vector<int> M_sizes_elements_with_ghost;

            std::vector<int> M_id_elements;
            std::vector<int> M_id_nodes;
            std::vector<int> M_rmap_nodes;

            std::vector<int> M_rmap_elements;
            std::vector<double> M_element_connectivity;

            graphmpi_type M_graphmpi;

            // update solution from explicit solver
            std::vector<std::vector<int>> M_extract_local_index;
            std::vector<int> M_recipients_proc_id;
            std::vector<int> M_local_ghosts_proc_id;
            std::vector<std::vector<int>> M_local_ghosts_local_index;

        private:
            int bamg_verbose;
            po::variables_map vm;

            LogLevel M_log_level;
            bool M_log_all;
            Communicator M_comm;

            boost::mpi::timer chrono, chrono_tot;

            bool M_use_restart;
    };
}
#endif
