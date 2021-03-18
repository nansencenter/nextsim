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
        CONSTANT           = 0,
        GENERIC_PS         = 1,
        ASR                = 2,
        ERAi               = 3,
        EC2                = 4,
        EC_ERAi            = 5,
        CFSR               = 6,
        CFSR_HI            = 7,
        EC2_AROME          = 8,
        EC2_AROME_ENSEMBLE = 9,
        ERA5               = 10
    };

    enum class OceanType
    {
        CONSTANT           = 0,
        TOPAZR             = 1,
        TOPAZF             = 2,
        MITGCM             = 3,
        TOPAZR_atrest      = 4,
        TOPAZR_ALTIMETER   = 5,
        COUPLED            = 6
    };

    enum class IceType
    {
        CONSTANT                    = 0,
        CONSTANT_PARTIAL            = 1,
        AMSRE                       = 2,
        TOPAZ4                      = 3,
        ARBITRARY                   = 5,
        AMSR2                       = 6,
        TOPAZ4F                     = 7,
        MITGCM                      = 8,
        TARGET                      = 9,
        OSISAF                      = 10,
        PIOMAS                      = 11,
        TOPAZ4FAMSR2                = 12,
        TOPAZ4FAMSR2OSISAF          = 13,
        CS2_SMOS                    = 14,
        CS2_SMOS_AMSR2              = 15,
        SMOS                        = 16,
        BINARY                      = 17,
        TOPAZ4OSISAFICESAT          = 18,
        TOPAZ4FAMSR2OSISAFNIC       = 19,
        TOPAZ4FAMSR2OSISAFNICWEEKLY = 20,
        CREG                        = 21,
        AMSR2CSTTHICK               = 22
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

    enum class WeldingType
    {
        NONE     = 0,
        ROACH    = 1
    };

    enum class FSDType
    {
        CONSTANT_SIZE   = 0,
        CONSTANT_AREA   = 1
    };
    enum class BreakupType
    {
        NONE  = 0,
        UNIFORM_SIZE = 1,
        ZHANG   =2,
        DUMONT  =3
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

    enum class FreezingPointType
    {
        LINEAR     = 0,
        NON_LINEAR = 1,
        UNESCO     = 2
    };


    enum class OceanHeatfluxScheme
    {
        BASIC      = 0,
        EXCHANGE   = 1
    };

    enum class DynamicsType
    {
        BBM        = 0,
        NO_MOTION  = 1,
        FREE_DRIFT = 2,
        EVP        = 3,
        mEVP       = 4
    };

} // setup

namespace schemes
{
    enum class damageDiscretisation
    {
        EXPLICIT    = 0,
        IMPLICIT    = 1,
        RECURSIVE   = 2
    };

    enum class tdType
    {
        FIXED            = 0,
        DAMAGE_DEPENDENT = 1,
    };

    enum class specificHumidity
    {
        ATMOSPHERE = 0,
        WATER      = 1,
        ICE        = 2
    };
} // schemes

} // Nextsim
