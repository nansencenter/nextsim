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

namespace Nextsim
{

    class MeshHandler
    {
        // Methods
        public:

            MeshHandler()
                : vm(Environment::vm())
            {
                bamg_verbose = vm["debugging.bamg_verbose"].as<int>();
            };

            void initBamg(BamgOpts *bamgopt, BamgGeom *bamggeom, BamgMesh *bamgmesh,
                    BamgGeom *bamggeom_root = NULL, BamgMesh *bamgmes_root = NULL,
                    BamgGeom *bamgeom_previous = NULL, BamgMesh *bamgmesh_previous = NULL, int rank = 0);

        protected:

        // Variables
        public:

        private:
            int bamg_verbose;
            po::variables_map vm;

    };
}
#endif
