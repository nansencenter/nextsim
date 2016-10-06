/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   enums.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Mon Oct 12 12:04:59 2015
 */

namespace Nextsim
{
namespace setup
{

    enum class AtmosphereType
    {
        CONSTANT = 0,
        ASR      = 1,
        ERAi     = 2,
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

    enum class IceType
	{
		CONSTANT  = 0,
		AMSRE     = 1,
		TOPAZ4    = 2,
        ARBITRARY = 4,
        AMSR2     = 5,
        TOPAZ4F   = 6,
        MITGCM    = 7,
        TARGET    = 8,
        OSISAF    = 9
	};

    enum class WaveType
    {
        CONSTANT = 0,
        WW3A     = 1,
        ERAI_WAVES_1DEG = 2
    };

    enum class BathymetryType
    {
        CONSTANT = 0,
        ETOPO    = 1
    };

    enum class IceCategoryType
    {
        CLASSIC     = 0,
        THIN_ICE    = 1,
        MULTI       = 2
    };

    enum class DrifterType
    {
        NONE          = 0,
        EQUALLYSPACED = 1,
        IABP          = 2
    };

    enum class DomainType
    {
        DEFAULT        = 0,
        KARA           = 1,
        BERINGKARA     = 2,
        BIGKARA        = 3,
        ARCTIC         = 4,
        BIGARCTIC      = 5,
        WIM            = 6
    };

    enum class MeshType
    {
        FROM_GMSH      = 0,
        FROM_SPLIT     = 1
    };

    enum class ThermoType
    {
        ZERO_LAYER = 0,
        WINTON     = 1
    };

} // setup
} // Nextsim
