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
        EC2      = 4,
        EC_ERAi  = 5,
        CFSR     = 6,
        CFSR_HI  = 7
    };

	enum class OceanType
	{
		CONSTANT = 0,
		TOPAZR   = 1,
		TOPAZF   = 2,
        MITGCM   = 3,
        TOPAZR_atrest   = 4,
        TOPAZR_ALTIMETER   = 5
    };

    enum class IceType
	{
		CONSTANT            = 0,
		CONSTANT_PARTIAL    = 1,
		AMSRE               = 2,
		TOPAZ4              = 3,
        ARBITRARY           = 5,
        AMSR2               = 6,
        TOPAZ4F             = 7,
        MITGCM              = 8,
        TARGET              = 9,
        OSISAF              = 10,
        PIOMAS              = 11,
        TOPAZ4FAMSR2        = 12,
        TOPAZ4FAMSR2OSISAF  = 13,
        CS2_SMOS            = 14,
        CS2_SMOS_AMSR2      = 15,
        SMOS                = 16,
        BINARY              = 17,
        TOPAZ4OSISAFICESAT  = 18,
        TOPAZ4FAMSR2OSISAFNIC= 19,
        TOPAZ4FAMSR2OSISAFNICWEEKLY= 20
	};

    enum class WaveType
    {
        SET_IN_WIM          = 0,
        CONSTANT            = 1,
        CONSTANT_PARTIAL    = 2,
        WW3A                = 3,
        ERAI_WAVES_1DEG     = 4
    };

    enum class BathymetryType
    {
        CONSTANT = 0,
        ETOPO    = 1
    };

    enum class BasalStressType
    {
        NONE     = 0,
        LEMIEUX  = 1,
        BOUILLON = 2
    };

    enum class IceCategoryType
    {
        CLASSIC     = 0,
        THIN_ICE    = 1,
        MULTI       = 2
    };

    enum class DomainType
    {
        DEFAULT        = 0,
        KARA           = 1,
        BERINGKARA     = 2,
        BIGKARA        = 3,
        ARCTIC         = 4,
        BIGARCTIC      = 5,
        WIM            = 6,
        FROM_SPLIT     = 7,
        UNREF          = 8
    };

    enum class MeshType
    {
        FROM_UNREF     = 0,
        FROM_SPLIT     = 1
    };

    enum class ThermoType
    {
        ZERO_LAYER = 0,
        WINTON     = 1
    };

    enum class DynamicsType
    {
        DEFAULT         = 0,
        NO_MOTION       = 1,
        FREE_DRIFT      = 2
    };

} // setup
} // Nextsim
