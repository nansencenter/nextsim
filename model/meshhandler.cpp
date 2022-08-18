/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   meshhandler.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 17 Aug 2022 07:59:06 CEST
 */

#include <meshhandler.hpp>

namespace Nextsim
{

//------------------------------------------------------------------------------------------------------
//! Initializes a Bamg mesh grid.
//! Called by the initMesh() function.
void
MeshHandler::initBamg(BamgOpts &bamgopt)
{

    bamgopt.Crack             = 0;
    bamgopt.anisomax          = 1e30;
    bamgopt.coeff             = 1;
    bamgopt.cutoff            = 1e-5;
    // bamgopt.err               = 0.01;
    bamgopt.errg              = 0.1;
    bamgopt.field             = NULL;
    bamgopt.gradation         = 1.5;
    bamgopt.Hessiantype       = 0;
    bamgopt.hmin              = 1e-100;
    bamgopt.hmax              = 1e100;
    bamgopt.hminVertices      = NULL;
    bamgopt.hmaxVertices      = NULL;
    bamgopt.hVertices         = NULL;
    bamgopt.KeepVertices      = 1;
    bamgopt.MaxCornerAngle    = 10;
    bamgopt.maxnbv            = 1e7;
    bamgopt.maxsubdiv         = 10;
    bamgopt.metric            = NULL;
    bamgopt.Metrictype        = 0;
    bamgopt.nbjacobi          = 1;
    bamgopt.nbsmooth          = 3;
    bamgopt.omega             = 1.8;
    bamgopt.power             = 1.;
    bamgopt.splitcorners      = 1; //the Devil!  Changed to 0, original 1 Phil
    bamgopt.geometricalmetric = 0;
    bamgopt.random            = true;
    bamgopt.verbose           = bamg_verbose;

    bamgopt.Check();

}//initBamg

}



