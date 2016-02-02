/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   enums.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Mon Oct 12 12:04:59 2015
 */

namespace Nextsim
{
namespace setup
{
    enum class InterpolationType
    {
        InterpFromGridToMesh = 0,
        InterpFromMeshToMesh2dx = 1,
        InterpFromMeshToMesh2dCavities = 2,
    };

    enum class AtmosphereType
    {
        CONSTANT = 0,
        ASR      = 1,
        ERA      = 2,
        EC       = 3,
        ASRFINAL = 4
    };

	enum class OceanType
	{
		CONSTANT = 0,
		TOPAZR   = 1,
		TOPAZF   = 2,
		MITGCM   = 3
    };

    enum class ConcentrationType
	{
		CONSTANT  = 0,
		AMSRE     = 1,
		TOPAZ4    = 2,
        ARBITRARY = 4,
        AMSRE2    = 5,
        TOPAZ4F   = 6,
        MITGCM    = 7
	};

    enum class ThicknessType
    {
        CONSTANT  = 0,
        SMOS      = 1,
        ARBITRARY = 2,
        ICESAT    = 3,
        TOPAZ4    = 4,
        MITGCM    = 5
    };

    enum class DamageType
    {
        CONSTANT = 0,
        RGPS     = 1,
        RCONC    = 2
    };

    enum class IceCategoryType
    {
        CLASSIC     = 0,
        THIN_ICE    = 1,
        MULTI       = 2
    };

    enum class SnowThicknessType
    {
        CONSTANT  = 0,
        ARBITRARY = 1,
        TOPAZ4    = 2,
        MITGCM    = 3
    };

    enum class DrifterType
    {
        EQUALLYSPACED = 0,
        IABP          = 1
    };

    enum class DomainType
    {
        DEFAULT        = 0,
        KARA           = 1,
        BERINGKARA     = 2,
        BIGKARA        = 3,
        ARCTIC         = 4,
        BIGARCTIC      = 5
    };

    enum class MeshType
    {
        FROM_GMSH      = 0,
        FROM_SPLIT     = 1
    };

} // setup
} // Nextsim
