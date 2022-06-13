/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gridoutput.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Thu Aug  4 09:47:27 CEST 2016
 */

#ifndef __ModelVariable_H
#define __ModelVariable_H 1

#include <environment.hpp>
#include <constants.hpp>

/**
 * @class GridOutput
 * @brief Manage output of gridded data
 *
 * @see
 *
 */

namespace Nextsim
{

class ModelVariable: public std::vector<double> // inherit from std::vector<double>
{

public:

    //////////////////////////////////////////////////////////////////////
    // Misc. enumerations and type definitions
    //////////////////////////////////////////////////////////////////////


    enum variableKind
    {
        nodal      =  0,
        elemental  =  1
    };

    enum interpTransformation
    {
        // transformation to be done before (and inverted after) interpolation at regrid time (and advection if using ALE)
        none     = 0,//do nothing (conservative variable already)
        conc     = 1,//multiply by M_conc
        thick    = 2,//multiply by M_thick
        enthalpy = 3,//convert to enthalpy (for M_tice[1])
    };

    enum variableID
    {
        // Prognostic variables
        M_conc            =  1,
        M_thick           =  2,
        M_damage          =  3,
        M_snow_thick      =  4,
        M_ridge_ratio     =  5,
        M_tice            =  6,
        M_sigma           =  7,
        M_sst             =  8,
        M_sss             =  9,
        M_tsurf_young     = 10,
        M_h_young         = 11,
        M_hs_young        = 12,
        M_conc_young      = 13,
        M_random_number   = 14,
        M_conc_fsd        = 15,
        M_fyi_fraction    = 16,
        M_age_det         = 17,
        M_age             = 18,
        M_conc_upd        = 19,
        M_conc_mech_fsd   = 20,
        M_cum_damage      = 21,
        M_cum_wave_damage = 22,
        M_divergence      = 23,
        M_conc_myi        = 24,
        M_melt_seconds    = 25,
        M_conc_summer     = 26,
        M_freeze_seconds  = 27,
        M_melt_onset      = 28,
        M_freeze_onset    = 29,
        M_del_hi_tend     = 30,
        M_thick_myi       = 31,
        M_thick_summer    = 32,

        // Diagnostic variables
        D_conc         = 100,
        D_thick        = 101,
        D_snow_thick   = 102,
        D_tsurf        = 103,
        D_sigma        = 104,
        D_Qa           = 105,
        D_Qsw          = 106,
        D_Qlw          = 107,
        D_Qsh          = 108,
        D_Qlh          = 109,
        D_Qo           = 110,
        D_Qnosun       = 111,
        D_Qsw_ocean    = 112,
        D_tau_ow       = 113,
        D_delS         = 114,
        D_fwflux       = 115,
        D_brine        = 116,
        D_evap         = 117,
        D_rain         = 118,
        D_dmax         = 119,
        D_dmean        = 120,
        D_Qassim       = 121,
        D_fwflux_ice   = 122,
        D_dcrit        = 123,
        D_sigma_p      = 124,
        D_vice_melt    = 125,
        D_newice       = 126,
        D_del_hi       = 127,
        D_del_hi_young = 128,
        D_mlt_bot      = 129,
        D_mlt_top      = 130,
        D_snow2ice     = 131,
        D_del_vi_young = 132,
        D_divergence   = 133,
        D_albedo       = 134,
        D_del_ci_mlt_myi   = 135,
        D_del_vi_mlt_myi   = 136,
        D_del_ci_rplnt_myi = 137,
        D_del_vi_rplnt_myi = 138,
        D_del_ci_ridge_myi = 139,
    };


private:
    variableID M_varID;
    variableKind M_var_kind;
    int M_component_number;
    std::string M_name;//variable name and name in restart file
    std::string M_export_name;//name in binary file
    bool M_prognostic;//is it a prognostic variable (in the restart file, and is it regridded)?
    bool M_exporting;//export to binary output? (Default which could be overridden by config file in future)
    interpTransformation M_interp_transformation;//transformation to use during interpolation/advection (ALE)
    double M_diffusivity;//diffusivity parameter
        // 0. for non added diffusion;
        // positive value for active diffusion in [m^2/s] (only non conservative implementation available)
    bool M_has_min = false; //does the variable have a strict min (eg >0)?
    bool M_has_max = false; //does the variable have a strict max (eg <1)?
    bool M_has_value_no_thick_ice = false; //does the variable have a designated value if no thick ice (eg 0, freezing temp)?
    double M_min_val;
    double M_max_val;
    double M_value_no_thick_ice;
    double M_tfr_ice;


public:
    ModelVariable(){}  // empty constructor
    ~ModelVariable(){} // destructor

    // main constructor
    ModelVariable(variableID id, int comp_num=-1)
        : M_varID(id), M_component_number(comp_num)
    {

        M_tfr_ice = -physical::mu*physical::si;
        // freezing point of ice - value for M_tice[i] in open water

        // see if it is an elemental variable
        bool elemental = this->initElemental();
        bool nodal = false;

        if(elemental)
            M_var_kind = variableKind::elemental;
        else
            // if it's not an elemental variable, see if is nodal
            nodal = this->initNodal();

        // TODO add the nodes
        if(nodal)
            M_var_kind = variableKind::nodal;
        else if(!elemental)
            //if it's neither elemental nor nodal raise an error
            throw std::logic_error("ModelVariable: variable ID not defined: "+std::to_string(M_varID));

        if(M_component_number>-1)
        {
            //we export vector components separately with modified name eg M_tice_[0,1,2], M_sigma_[0,1,2]
            M_name += "_" + std::to_string(M_component_number);
            M_export_name += "_" + std::to_string(M_component_number);
        }

    }

private:
    bool initElemental();
    bool initNodal() { return false; } // TODO implement nodal variables

public:

    // get attributes
    variableID varID() { return M_varID; }
    variableKind varKind() { return M_var_kind; }
    int componentNumber() { return M_component_number; }
    std::string name() { return M_name; }
    std::string exportName() { return M_export_name; }
    bool isPrognostic() {return M_prognostic; }
    bool exporting() { return M_exporting; }
    interpTransformation getInterpTransformation() { return M_interp_transformation; }
    double diffusivity() { return M_diffusivity; }
    bool hasMinVal() { return M_has_min; }
    bool hasMaxVal() { return M_has_max; }
    bool hasValueNoThickIce() { return M_has_value_no_thick_ice; }
    double minVal() { return M_min_val; }
    double maxVal() { return M_max_val; }
    double valueNoThickIce() { return M_value_no_thick_ice; }

    // set attributes
    void setExporting(bool const& do_export) { M_exporting = do_export; }

    void testCall();
};

} // Nextsim
#endif // __ModelVariable_H
