/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gridoutput.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Thu Aug  4 09:46:28 CEST 2016
 */

#include <model_variable.hpp>

/**
 * @class GridOutput
 * @brief Manage output of gridded data
 *
 * @see
 *
*/

namespace Nextsim
{

bool
ModelVariable::initElemental()
{
    bool elemental = true;
    switch (varID)
    {

        case (variableID::M_conc):
            // concentration of thick ice
            name = "M_conc";
            export_name = "Concentration_thick_ice";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::conservative;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;

        case (variableID::M_thick):
            // thickness of thick ice [m]
            name = "M_thick";
            export_name = "Thick_ice";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::conservative;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;

        case (variableID::M_snow_thick):
            // snow thickness on thick ice [m]
            name = "M_snow_thick";
            export_name = "Snow_thick_ice";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::conservative;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;

        case (variableID::M_sigma):
            // stress tensor [Pa]
            name = "M_sigma";
            export_name = "M_sigma";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::conservative;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            if(component_number<0 || component_number>2)
                throw std::runtime_error(
                        "Unauthorised component number for M_sigma: "
                        +std::to_string(component_number));
            break;

        case (variableID::M_damage):
            // damage
            name = "M_damage";
            export_name = "Damage";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::meshToMesh;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;

        case (variableID::M_ridge_ratio):
            // volume ratio of ridged ice
            name = "M_ridge_ratio";
            export_name = "Ridge_ratio";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::conservative;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;

        case (variableID::M_tice):
            // ice temperature [degC]
            name = "M_tice";
            export_name = "Tice";
            in_restart = true;
            regrid = true;
            diffusivity = 0;
            switch(component_number)
            {
                case(0):
                    // surface temperature
                    interp_method = interpMethod::meshToMesh;
                    interp_transformation = interpTransformation::none;
                    break;

                case(1):
                    // mean temperature of top layer
                    interp_method = interpMethod::conservative;
                    interp_transformation = interpTransformation::enthalpy; // (Winton, 2000: eq 39) times volume with f1=1
                    break;

                case(2):
                    // mean temperature of lower layer
                    interp_method = interpMethod::conservative;
                    interp_transformation = interpTransformation::thick; // (Winton, 2000: eq 39) times volume with f1=0
                    break;

                default:
                    throw std::runtime_error(
                            "Unauthorised component number for M_tice: "
                            +std::to_string(component_number));
            }
            break;

        case (variableID::M_sst):
            // sea surface temperature of slab ocean [deg C]
            name = "M_sst";
            export_name = "SST";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::meshToMesh;
            interp_transformation = interpTransformation::none;
            diffusivity = Environment::vm()["thermo.diffusivity_sst"].as<double>();
            break;

        case (variableID::M_sss):
            // sea surface salinity of slab ocean [ppt]
            name = "M_sss";
            export_name = "SSS";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::meshToMesh;
            interp_transformation = interpTransformation::none;
            diffusivity = Environment::vm()["thermo.diffusivity_sss"].as<double>();
            break;

        case (variableID::M_tsurf_thin_ice):
            // surface temperature above thin ice [deg C]
            name = "M_tsurf_thin_ice";
            export_name = "Tsurf_thin_ice";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::meshToMesh;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;

        case (variableID::M_h_thin):
            // thickness of thin ice [m]
            name = "M_h_thin";
            export_name = "Thin_ice";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::conservative;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;

        case (variableID::M_hs_thin):
            // snow thickness on thin ice [m]
            name = "M_hs_thin";
            export_name = "Snow_thin_ice";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::conservative;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;

        case (variableID::M_conc_thin):
            // concentration of thin ice
            name = "M_conc_thin";
            export_name = "Concentration_thin_ice";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::conservative;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;

        case (variableID::M_random_number):
            // random number (uniform in [0,1]) used in the definition of the cohesion as follows:
            // M_Cohesion[i] = C_fix+C_alea*(M_random_number[i]-0.5);
            name = "M_random_number";
            export_name = "M_random_number";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::meshToMesh;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;

#if 0
        case (variableID::M_fyi_fraction):
            name = "M_fyi_fraction";
            export_name = "Fyi_fraction";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::conservative;
            interp_transformation = interpTransformation::conc;
            diffusivity = 0;
            break;

        case (variableID::M_age_obs):
            name = "M_age_obs";
            export_name = "Age_o";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::conservative;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;

        case (variableID::M_age):
            name = "M_age";
            export_name = "Age";
            in_restart = true;
            regrid = true;
            interp_method = interpMethod::conservative;
            interp_transformation = interpTransformation::none;
            diffusivity = 0;
            break;
#endif

        // Diagnostic variables
        case (variableID::D_conc):
            // total concentration (thin + thick ice)
            name = "D_conc";
            export_name = "Concentration";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_thick):
            // total thickness (thin + thick ice) [m]
            name = "D_thick";
            export_name = "Thickness";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_snow_thick):
            // total snow thickness (thin + thick ice) [m]
            name = "D_snow_thick";
            export_name = "Snow";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_tsurf):
            // mean surface temperature (ocean + thin ice + thick ice) [deg C]
            name = "D_tsurf";
            export_name = "Tsurf";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_sigma):
            // principal stresses [Pa]
            name = "D_sigma";
            export_name = "Sigma";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            if(component_number<0 || component_number>1)
                throw std::runtime_error(
                        "Unauthorised component number for D_sigma: "
                        +std::to_string(component_number));
            break;

        case (variableID::D_Qa):
            // Total heat flux to the atmosphere [W/m^2]
            name = "D_Qa";
            export_name = "Qa";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_Qsw):
            // Short wave flux to the atmosphere [W/m^2]
            name = "D_Qsw";
            export_name = "Qsw";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_Qlw):
            // Long wave flux to the atmosphere [W/m^2]
            name = "D_Qlw";
            export_name = "Qlw";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_Qsh):
            // Sensible heat flux to the atmosphere [W/m^2]
            name = "D_Qsh";
            export_name = "Qsh";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_Qlh):
            // Latent heat flux to the atmosphere [W/m^2]
            name = "D_Qlh";
            export_name = "Qlh";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_Qo):
            // Total heat lost by ocean [W/m^2]
            name = "D_Qo";
            export_name = "Qocean";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_delS):
            // Salt balance of the ocean (all sources) [kg/day]
            name = "D_delS";
            export_name = "Saltflux";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_emp):
            // Freshwater balance at the surface [kg/m^2/s]
            name = "D_emp";
            export_name = "D_emp";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        case (variableID::D_brine):
            // Brine release - kg/m^2/s
            name = "D_brine";
            export_name = "D_brine";
            in_restart = false;
            regrid = false;
            diffusivity = 0;
            break;

        default:
            elemental = false;
    }

    return elemental;
}//initElemental


void ModelVariable::testCall()
{
    // test basic attributes
    std::cout<<"name = "<<this->name<<"\n";
    std::cout<<"export_name = "<<this->export_name<<"\n";

#if 0
    // test inherited methods
    // - assign
    // - begin
    // - end
    int N = 5;
    int i = 0;
    this->assign(N, 4.);
    for(auto it=this->begin(); it!=this->end(); it++, i++)
        std::cout<<"i<"<<N<<"? "<<i <<" value = 4? "<<*it<<"\n";

    // - resize
    // - std::fill
    // - begin
    // - end
    N = 4;
    i = 0;
    this->resize(N);
    std::fill(this->begin(), this->end(), 2.5);
    for(auto it=this->begin(); it!=this->end(); it++, i++)
        std::cout<<"i<"<<N<<"? "<<i <<" value = 2.5? "<<*it<<"\n";

    // - size
    std::cout<<"size = "<<N<<"? "<<this->size()<<"\n";
#endif

    std::abort();
}//testCall


} // Nextsim
