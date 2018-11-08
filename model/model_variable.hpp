/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gridoutput.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Thu Aug  4 09:47:27 CEST 2016
 */

#ifndef __ModelVariable_H
#define __ModelVariable_H 1

#include <environment.hpp>

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

    enum InterpMethod
    {
        meshToMesh   = 0,
        conservative = 1
    };

    enum InterpTransformation
    {
        none     = 0,//do nothing (conservative variable already)
        conc     = 1,//multiply by M_conc
        thick    = 2,//multiply by M_thick
        enthalpy = 3,//convert to enthalpy (for M_tice[1])
    };

    enum variableID
    {
        // Prognostic variables
        M_conc          =  1,
        M_thick         =  2,
        M_damage        =  3,
        M_snow_thick    =  4,
        M_ridge_ratio   =  5,
        M_tice          =  6,
        M_sigma         =  7,
        M_sst           =  8,
        M_sss           =  9,
        M_tsurf_thin    = 10,
        M_h_thin        = 11,
        M_hs_thin       = 12,
        M_conc_thin     = 13,
        M_random_number = 14,
        //M_fyi_fraction  = 15,
        //M_age_obs       = 16,
        //M_age           = 17,

        // Diagnostic variables
        D_conc       = 100,
        D_thick      = 101,
        D_snow_thick = 102,
        D_tsurf      = 103,
        D_sigma      = 104,
        D_Qa         = 105,
        D_Qsw        = 106,
        D_Qlw        = 107,
        D_Qsh        = 108,
        D_Qlh        = 109,
        D_Qo         = 110,
        D_Qnosun     = 111,
        D_Qsw_ocean  = 112,
        D_delS       = 113,
        D_emp        = 114,
        D_brine      = 115
    };


private:
    variableID M_varID;
    variableKind M_var_kind;
    int M_component_number;
    std::string M_name;//variable name and name in restart file
    std::string M_export_name;//name in binary file
    bool M_prognostic;//is it a prognostic variable (in the restart file, and is it regridded)?
    bool M_exporting;//export to binary output? (Default which could be overridden by config file in future)
    InterpMethod M_interp_method;//interpolation method during regridding
    InterpTransformation M_interp_transformation;//transformation to use during interpolation/advection (ALE)
    double M_diffusivity;//diffusivity parameter
        // 0. for non added diffusion;
        // positive value for active diffusion in [m^2/s] (only non conservative implementation available)
    bool M_has_min = false;
    bool M_has_max = false;
    double M_min_val;
    double M_max_val;


public:
    ModelVariable(){}  // empty constructor
    ~ModelVariable(){} // destructor

    // main constructor
    ModelVariable(variableID id, int comp_num=-1)
        : M_varID(id), M_component_number(comp_num)
    {

        // see if it is an elemental variable
        bool elemental = this->initElemental();
        bool nodal = false;

        if(elemental)
            M_var_kind = variableKind::elemental;
        else
            // if it's not an elemental variable, see if is nodal
            nodal = this->initNodal();

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
    InterpMethod interpMethod() {return M_interp_method; }
    InterpTransformation interpTransformation() { return M_interp_transformation; }
    double diffusivity() { return M_diffusivity; }
    bool hasMinVal() { return M_has_min; }
    bool hasMaxVal() { return M_has_max; }
    double minVal() { return M_min_val; }
    double maxVal() { return M_max_val; }

    // set attributes
    void setExporting(bool const& do_export) { M_exporting = do_export; }

    void testCall();
};

} // Nextsim
#endif // __ModelVariable_H
