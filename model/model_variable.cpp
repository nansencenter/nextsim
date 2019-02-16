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
    switch (M_varID)
    {

        case (variableID::M_conc):
            // concentration of thick ice
            M_name = "M_conc";
            M_export_name = "Concentration_thick_ice";
            M_prognostic = true;
            M_exporting = false;
            M_interp_method = interpMethod::conservative;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_thick):
            // thickness of thick ice [m]
            M_name = "M_thick";
            M_export_name = "Thick_ice";
            M_prognostic = true;
            M_exporting = false;
            M_interp_method = interpMethod::conservative;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_snow_thick):
            // snow thickness on thick ice [m]
            M_name = "M_snow_thick";
            M_export_name = "Snow_thick_ice";
            M_prognostic = true;
            M_exporting = false;
            M_interp_method = interpMethod::conservative;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_sigma):
            // stress tensor [Pa]
            M_name = "M_sigma";
            M_export_name = "M_sigma";
            M_prognostic = true;
            M_exporting = false;
            M_interp_method = interpMethod::conservative;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            if(M_component_number<0 || M_component_number>2)
                throw std::runtime_error(
                        "Unauthorised component number for M_sigma: "
                        +std::to_string(M_component_number));
            break;

        case (variableID::M_damage):
            // damage
            M_name = "M_damage";
            M_export_name = "Damage";
            M_prognostic = true;
            M_exporting = true;
            M_interp_method = interpMethod::nearest_neighbour;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_ridge_ratio):
            // volume ratio of ridged ice
            M_name = "M_ridge_ratio";
            M_export_name = "Ridge_ratio";
            M_prognostic = true;
            M_exporting = true;
            M_interp_method = interpMethod::conservative;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_tice):
            // ice temperature [degC]
            M_name = "M_tice";
            M_export_name = "Tice";
            M_prognostic = true;
            M_exporting = true;
            M_diffusivity = 0;
            switch(M_component_number)
            {
                case(0):
                    // surface temperature
                    M_interp_method = interpMethod::nearest_neighbour;
                    M_interp_transformation = interpTransformation::none;
                    break;

                case(1):
                    // mean temperature of top layer
                    M_interp_method = interpMethod::conservative;
                    M_interp_transformation = interpTransformation::enthalpy; // (Winton, 2000: eq 39) times volume with f1=1
                    break;

                case(2):
                    // mean temperature of lower layer
                    M_interp_method = interpMethod::conservative;
                    M_interp_transformation = interpTransformation::thick; // (Winton, 2000: eq 39) times volume with f1=0
                    break;

                default:
                    throw std::runtime_error(
                            "Unauthorised component number for M_tice: "
                            +std::to_string(M_component_number));
            }
            break;

        case (variableID::M_sst):
            // sea surface temperature of slab ocean [deg C]
            M_name = "M_sst";
            M_export_name = "SST";
            M_prognostic = true;
            M_exporting = true;
            M_interp_method = interpMethod::nearest_neighbour;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = Environment::vm()["thermo.diffusivity_sst"].as<double>();
            break;

        case (variableID::M_sss):
            // sea surface salinity of slab ocean [ppt]
            M_name = "M_sss";
            M_export_name = "SSS";
            M_prognostic = true;
            M_exporting = true;
            M_interp_method = interpMethod::nearest_neighbour;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = Environment::vm()["thermo.diffusivity_sss"].as<double>();
            break;

        case (variableID::M_tsurf_thin):
            // surface temperature above thin ice [deg C]
            M_name = "M_tsurf_thin";
            M_export_name = "Tsurf_thin_ice";
            M_prognostic = true;
            M_exporting = true;
            M_interp_method = interpMethod::nearest_neighbour;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_h_thin):
            // thickness of thin ice [m]
            M_name = "M_h_thin";
            M_export_name = "Thin_ice";
            M_prognostic = true;
            M_exporting = true;
            M_interp_method = interpMethod::conservative;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_hs_thin):
            // snow thickness on thin ice [m]
            M_name = "M_hs_thin";
            M_export_name = "Snow_thin_ice";
            M_prognostic = true;
            M_exporting = true;
            M_interp_method = interpMethod::conservative;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_conc_thin):
            // concentration of thin ice
            M_name = "M_conc_thin";
            M_export_name = "Concentration_thin_ice";
            M_prognostic = true;
            M_exporting = true;
            M_interp_method = interpMethod::conservative;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_random_number):
            // random number (uniform in [0,1]) used in the definition of the cohesion as follows:
            // M_Cohesion[i] = C_fix+C_alea*(M_random_number[i]-0.5);
            M_name = "M_random_number";
            M_export_name = "M_random_number";
            M_prognostic = true;
            M_exporting = false;
            M_interp_method = interpMethod::nearest_neighbour;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_fyi_fraction):
            M_name = "M_fyi_fraction";
            M_export_name = "Fyi_fraction";
            M_prognostic = true;
            M_exporting = true;
            M_interp_method = interpMethod::conservative;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_age_det):
            M_name = "M_age_det";
            M_export_name = "Age_d";
            M_prognostic = true;
            M_exporting = true;
            M_interp_method = interpMethod::conservative;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_age):
            M_name = "M_age";
            M_export_name = "Age";
            M_prognostic = true;
            M_exporting = true;
            M_interp_method = interpMethod::conservative;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        // Diagnostic variables
        case (variableID::D_conc):
            // total concentration (thin + thick ice)
            M_name = "D_conc";
            M_export_name = "Concentration";
            M_prognostic = false;
            M_exporting = true;
            M_diffusivity = 0;
            break;

        case (variableID::D_thick):
            // total thickness (thin + thick ice) [m]
            M_name = "D_thick";
            M_export_name = "Thickness";
            M_prognostic = false;
            M_exporting = true;
            M_diffusivity = 0;
            break;

        case (variableID::D_snow_thick):
            // total snow thickness (thin + thick ice) [m]
            M_name = "D_snow_thick";
            M_export_name = "Snow";
            M_prognostic = false;
            M_exporting = true;
            M_diffusivity = 0;
            break;

        case (variableID::D_tsurf):
            // mean surface temperature (ocean + thin ice + thick ice) [deg C]
            M_name = "D_tsurf";
            M_export_name = "Tsurf";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        case (variableID::D_sigma):
            // principal stresses [Pa]
            M_name = "D_sigma";
            M_export_name = "Sigma";
            M_prognostic = false;
            M_exporting = true;
            M_diffusivity = 0;
            if(M_component_number<0 || M_component_number>1)
                throw std::runtime_error(
                        "Unauthorised component number for D_sigma: "
                        +std::to_string(M_component_number));
            break;

        case (variableID::D_Qa):
            // Total heat flux to the atmosphere [W/m^2]
            M_name = "D_Qa";
            M_export_name = "Qatm";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        case (variableID::D_Qsw):
            // Short wave flux to the atmosphere [W/m^2]
            M_name = "D_Qsw";
            M_export_name = "Qsw";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        case (variableID::D_Qlw):
            // Long wave flux to the atmosphere [W/m^2]
            M_name = "D_Qlw";
            M_export_name = "Qlw";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        case (variableID::D_Qsh):
            // Sensible heat flux to the atmosphere [W/m^2]
            M_name = "D_Qsh";
            M_export_name = "Qsh";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        case (variableID::D_Qlh):
            // Latent heat flux to the atmosphere [W/m^2]
            M_name = "D_Qlh";
            M_export_name = "Qlh";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        case (variableID::D_Qo):
            // Total heat lost by ocean [W/m^2]
            M_name = "D_Qo";
            M_export_name = "Qocean";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        case (variableID::D_delS):
            // Salt balance of the ocean (all sources) [kg/day]
            M_name = "D_delS";
            M_export_name = "Saltflux";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        case (variableID::D_emp):
            // Freshwater balance at the surface [kg/m^2/s]
            M_name = "D_emp";
            M_export_name = "D_emp";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        case (variableID::D_brine):
            // Brine release - kg/m^2/s
            M_name = "D_brine";
            M_export_name = "D_brine";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        default:
            elemental = false;
    }

    return elemental;
}//initElemental


void ModelVariable::testCall()
{
    // test basic attributes
    std::cout<<"name = "<<this->M_name<<"\n";
    std::cout<<"export_name = "<<this->M_export_name<<"\n";

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
