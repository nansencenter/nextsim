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

            void initBamg();
            void rootMeshProcessing();

        protected:
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

            template<typename FEMeshType>
            double measure(element_type const& element, FEMeshType const& mesh) const;

            template<typename FEMeshType>
            double measure(element_type const& element, FEMeshType const& mesh,
                           std::vector<double> const& um, double factor = 1.) const;

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

            BamgOpts bamgopt;
            BamgMesh bamgmesh;
            BamgGeom bamggeom;

            BamgMesh bamgmesh_root;
            BamgGeom bamggeom_root;

            BamgOpts bamgopt_previous;
            BamgMesh bamgmesh_previous;
            BamgGeom bamggeom_previous;

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
