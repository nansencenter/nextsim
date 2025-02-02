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
        CONSTANT                 = 0,
        GENERIC_PS               = 1,
        ASR                      = 2,
        ECMWF_NRT                = 3,
        CFSR                     = 4,
        CFSR_HI                  = 5,
        ECMWF_NRT_AROME          = 6,
        ECMWF_NRT_AROME_ENSEMBLE = 7,
        ERA5                     = 8
    };

    enum class OceanType
    {
        CONSTANT          = 0,
        TOPAZ4R           = 1,
        TOPAZ4NRT         = 2,
        TOPAZ5NRT         = 3,
        MITGCM            = 4,
        TOPAZ4R_ATREST    = 5,
        TOPAZ4R_ALTIMETER = 6,
        COUPLED           = 7,
        GLORYS12R         = 8
    };

    enum class IceType
    {
        CONSTANT                         = 0,
        CONSTANT_PARTIAL                 = 1,
        AMSRE                            = 2,
        TOPAZ4R                          = 3,
        ARBITRARY                        = 5,
        AMSR2                            = 6,
        TOPAZ4NRT                        = 7,
        TOPAZ5NRT                        = 8,
        MITGCM                           = 9,
        TARGET                           = 10,
        OSISAF                           = 11,
        PIOMAS                           = 12,
        TOPAZ4NRT_AMSR2                  = 13,
        TOPAZ4NRT_AMSR2_OSISAF           = 14,
        CS2_SMOS                         = 15,
        CS2_SMOS_AMSR2                   = 16,
        SMOS                             = 17,
        BINARY                           = 18,
        TOPAZ4R_OSISAF_ICESAT            = 19,
        TOPAZ4NRT_AMSR2_OSISAF_NIC       = 20,
        TOPAZ4NRT_AMSR2_OSISAF_NICWEEKLY = 21,
        NEMO                             = 22,
        CICE                             = 23,
        AMSR2_CSTTHICK                   = 24,
        GLORYS12R                        = 25
    };

    enum class WaveType
    {
        SET_IN_WIM       = 0,
        CONSTANT         = 1,
        CONSTANT_PARTIAL = 2,
        WW3A             = 3,
        ERAI_WAVES_1DEG  = 4
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
        CLASSIC   = 0,
        YOUNG_ICE = 1,
        MULTI     = 2
    };

    enum class WeldingType
    {
        NONE  = 0,
        ROACH = 1
    };

    enum class FSDType
    {
        CONSTANT_SIZE = 0,
        CONSTANT_AREA = 1
    };
    enum class BreakupType
    {
        NONE         = 0,
        UNIFORM_SIZE = 1,
        ZHANG        = 2,
        DUMONT       = 3
    };
    enum class MeshType
    {
        FROM_UNREF = 0,
        FROM_SPLIT = 1
    };

    enum class ThermoType
    {
        ZERO_LAYER = 0,
        WINTON     = 1
    };

    enum class FreezingPointType
    {
        LINEAR = 0,
        UNESCO = 1
    };


    enum class OceanHeatfluxScheme
    {
        BASIC    = 0,
        EXCHANGE = 1
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
    enum class specificHumidity
    {
        ATMOSPHERE = 0,
        WATER      = 1,
        ICE        = 2
    };
} // schemes

} // Nextsim
