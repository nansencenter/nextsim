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
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            M_has_max = true;
            M_max_val = 1.;
            break;

        case (variableID::M_thick):
            // thickness of thick ice [m]
            M_name = "M_thick";
            M_export_name = "Thick_ice";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            break;

        case (variableID::M_snow_thick):
            // snow thickness on thick ice [m]
            M_name = "M_snow_thick";
            M_export_name = "Snow_thick_ice";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            break;

        case (variableID::M_sigma):
            // stress tensor [Pa]
            M_name = "M_sigma";
            M_export_name = "M_sigma";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
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
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            M_has_max = true;
            M_max_val = 1.0 - 1e-10;
            break;

        case (variableID::M_cum_damage):
            // damage
            M_name = "M_cum_damage";
            M_export_name = "Cumulated_damage";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            break;

        case (variableID::M_cum_wave_damage):
            // damage
            M_name = "M_cum_wave_damage";
            M_export_name = "Cumulated_wave_damage";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            break;

        case (variableID::M_ridge_ratio):
            // volume ratio of ridged ice
            M_name = "M_ridge_ratio";
            M_export_name = "Ridge_ratio";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            M_has_max = true;
            M_max_val = 1.;
            break;

        case (variableID::M_tice):
            // ice temperature [degC]
            M_name = "M_tice";
            M_export_name = "Tice";
            M_prognostic = true;
            M_exporting = false;
            M_diffusivity = 0.;
            switch(M_component_number)
            {
                case(0):
                    // surface temperature
                    M_interp_transformation = interpTransformation::none;
                    break;

                case(1):
                    // mean temperature of top layer
                    M_interp_transformation = interpTransformation::enthalpy; // (Winton, 2000: eq 39) times volume with f1=1
                    break;

                case(2):
                    // mean temperature of lower layer
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
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = Environment::vm()["thermo.diffusivity_sst"].as<double>();
            break;

        case (variableID::M_sss):
            // sea surface salinity of slab ocean [ppt]
            M_name = "M_sss";
            M_export_name = "SSS";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = Environment::vm()["thermo.diffusivity_sss"].as<double>();
            break;

        case (variableID::M_tsurf_young):
            // surface temperature above young ice [deg C]
            M_name = "M_tsurf_young";
            M_export_name = "Tsurf_young_ice";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_h_young):
            // thickness of young ice [m]
            M_name = "M_h_young";
            M_export_name = "Thickness_young_ice";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            break;

        case (variableID::M_hs_young):
            // snow thickness on young ice [m]
            M_name = "M_hs_young";
            M_export_name = "Snow_young_ice";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            break;

        case (variableID::M_conc_young):
            // concentration of young ice
            M_name = "M_conc_young";
            M_export_name = "Concentration_young_ice";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            M_has_max = true;
            M_max_val = 1.;
            break;

        case (variableID::M_random_number):
            // random number (uniform in [0,1]) used in the definition of the cohesion as follows:
            // M_Cohesion[i] = C_fix+C_alea*M_random_number[i];
            M_name = "M_random_number";
            M_export_name = "M_random_number";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0;
            break;

        case (variableID::M_conc_fsd):
            // concentration of young ice
            M_name = "M_conc_fsd";
            M_export_name = "Concentration_fsd";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            M_has_max = true;
            M_max_val = 1.;
            break;

        case (variableID::M_conc_mech_fsd):
            // concentration of young ice
            M_name = "M_conc_mech_fsd";
            M_export_name = "Concentration_mech_fsd";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            M_has_max = true;
            M_max_val = 1.;
            break;

        case (variableID::M_fyi_fraction):
            M_name = "M_fyi_fraction";
            M_export_name = "Fyi_fraction";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            M_has_max = true;
            M_max_val = 1.;
            break;

        case (variableID::M_age_det):
            M_name = "M_age_det";
            M_export_name = "Age_d";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            break;

        case (variableID::M_age):
            M_name = "M_age";
            M_export_name = "Age";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            break;

        case (variableID::M_divergence):
            // update of concentration by assimilation
            M_name = "M_divergence";
            M_export_name = "Divergence";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            break;

        case (variableID::M_conc_upd):
            // update of concentration by assimilation
            M_name = "M_conc_upd";
            M_export_name = "Concentration_update";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = -1.;
            M_has_max = true;
            M_max_val = 1.;
            break;
        // MYI Variables    
        case (variableID::M_conc_myi):
            M_name = "M_conc_myi";
            M_export_name = "Concentration_myi";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            M_has_max = true;
            M_max_val = 1.;
            break;

        case (variableID::M_thick_myi):
            M_name = "M_thick_myi";
            M_export_name = "Thickness_myi";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            break;

        case (variableID::M_conc_summer):
            M_name = "M_conc_summer";
            M_export_name = "Concentration_summer";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            M_has_max = true;
            M_max_val = 1.;
            break;

        case (variableID::M_thick_summer):
            M_name = "M_thick_summer";
            M_export_name = "Thickness_summer";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            break;

        case (variableID::M_freeze_days):
            M_name = "M_freeze_days";
            M_export_name = "Ice_freezing_days";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            break;
        
        case (variableID::M_freeze_onset):
            M_name = "M_freeze_onset";
            M_export_name = "Ice_freezing_onset";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            M_has_min = true;
            M_min_val = 0.;
            M_has_max = true;
            M_max_val = 1.;
            break;

        case (variableID::M_del_vi_tend):
            M_name = "M_del_vi_tend";
            M_export_name = "Ice_volume_tendency_daily";
            M_prognostic = true;
            M_exporting = false;
            M_interp_transformation = interpTransformation::none;
            M_diffusivity = 0.;
            break;
        
        // Diagnostic variables
        case (variableID::D_conc):
            // total concentration (young + thick ice)
            M_name = "D_conc";
            M_export_name = "Concentration";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_thick):
            // total thickness (young + thick ice) [m]
            M_name = "D_thick";
            M_export_name = "Thickness";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_snow_thick):
            // total snow thickness (young + thick ice) [m]
            M_name = "D_snow_thick";
            M_export_name = "Snow";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_tsurf):
            // mean surface temperature (ocean + young ice + thick ice) [deg C]
            M_name = "D_tsurf";
            M_export_name = "Tsurf";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_sigma):
            // principal stresses [Pa]
            M_name = "D_sigma";
            M_export_name = "Sigma";
            M_prognostic = false;
            M_exporting = false;
            if(M_component_number<0 || M_component_number>1)
                throw std::runtime_error(
                        "Unauthorised component number for D_sigma: "
                        +std::to_string(M_component_number));
            break;

        case (variableID::D_divergence):
            // Divergence of ice flow [1/s]
            M_name = "D_divergence";
            M_export_name = "Divergence";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_Qa):
            // Total heat flux to the atmosphere [W/m^2]
            M_name = "D_Qa";
            M_export_name = "Qatm";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_Qsw):
            // Short wave flux to the atmosphere [W/m^2]
            M_name = "D_Qsw";
            M_export_name = "Qsw";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_Qlw):
            // Long wave flux to the atmosphere [W/m^2]
            M_name = "D_Qlw";
            M_export_name = "Qlw";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_Qsh):
            // Sensible heat flux to the atmosphere [W/m^2]
            M_name = "D_Qsh";
            M_export_name = "Qsh";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_Qlh):
            // Latent heat flux to the atmosphere [W/m^2]
            M_name = "D_Qlh";
            M_export_name = "Qlh";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_Qo):
            // Total heat lost by ocean [W/m^2]
            M_name = "D_Qo";
            M_export_name = "Qocean";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_Qnosun):
            // Non-solar fluxes to ocean [W/m^2]
            M_name = "D_Qnosun";
            M_export_name = "Qnosun";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_Qsw_ocean):
            // Short-wave fluxes to ocean [W/m^2]
            M_name = "D_Qsw_ocean";
            M_export_name = "Qsw_ocean";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_Qassim):
            // flux from assim [W/m^2]
            M_name = "D_Qassim";
            M_export_name = "Qassim";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_delS):
            // Salt balance of the ocean (all sources) [kg/day]
            M_name = "D_delS";
            M_export_name = "Saltflux";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_fwflux):
            // Freshwater balance at the surface [kg/m^2/s]
            M_name = "D_fwflux";
            M_export_name = "D_fwflux";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_fwflux_ice):
            // Freshwater flux at the surface due to ice processes [kg/m^2/s]
            M_name = "D_fwflux_ice";
            M_export_name = "D_fwflux_ice";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_vice_melt):
            // ice volume gain/loss by freezing/melt [m/day]
            M_name = "D_vice_melt";
            M_export_name = "D_vice_melt";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_del_vi_young):
            // Young ice volume gain/loss by freezing/melt [m/day]
            M_name = "D_del_vi_young";
            M_export_name = "D_del_vi_young";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_del_hi):
            // ice growth/melt rate [m/day]
            M_name = "D_del_hi";
            M_export_name = "D_del_hi";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_del_hi_young):
            //young ice growth/melt rate [m/day]
            M_name = "D_del_hi_young";
            M_export_name = "D_del_hi_young";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_newice):
            // ice growth in open water [m/day]
            M_name = "D_newice";
            M_export_name = "D_newice";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_mlt_bot):
            // ice melt at bottom [m/day]
            M_name = "D_mlt_bot";
            M_export_name = "D_mlt_bot";
            M_prognostic = false;
            M_exporting = false;
            break;
        
        case (variableID::D_mlt_top):
            // ice melt at top [m/day]
            M_name = "D_mlt_top";
            M_export_name = "D_mlt_top";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_snow2ice):
            // ice formed from snow by flooding [m/day]
            M_name = "D_snow2ice";
            M_export_name = "D_snow2ice";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_brine):
            // Brine release - kg/m^2/s
            M_name = "D_brine";
            M_export_name = "D_brine";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_dmax):
            // mean floe diameter
            M_name = "D_dmax";
            M_export_name = "Dmax";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_dmean):
            // mean floe diameter
            M_name = "D_dmean";
            M_export_name = "Dmean";
            M_prognostic = false;
            M_exporting = false;

        case (variableID::M_drag_ui):
            // mean floe diameter
            M_name = "M_drag_ui";
            M_export_name = "ice-atmosphere_drag";
            M_prognostic = true;
            M_exporting = false;

        case (variableID::M_drag_ti):
            // mean floe diameter
            M_name = "M_drag_ti";
            M_export_name = "ice-atmosphere_thermo_drag";
            M_prognostic = true;
            M_exporting = false;

        case (variableID::M_drag_ui_young):
            // mean floe diameter
            M_name = "M_drag_ui_young";
            M_export_name = "ice-atmosphere_drag_young";
            M_prognostic = true;
            M_exporting = false;

        case (variableID::M_drag_ti_young):
            // mean floe diameter
            M_name = "M_drag_ti_young";
            M_export_name = "ice-atmosphere_thermo_drag_young";
            M_prognostic = true;
            M_exporting = false;

        case (variableID::D_tau_ow):
            // Ocean atmosphere drag coefficient - still needs to be multiplied with the wind [Pa/s/m]
            M_name = "D_tau_ow";
            M_export_name = "D_tau_ow";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_evap):
            // Evaporation - kg/m^2/s
            M_name = "D_evap";
            M_export_name = "D_evap";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        case (variableID::D_rain):
            // Rain - kg/m^2/s
            M_name = "D_rain";
            M_export_name = "D_rain";
            M_prognostic = false;
            M_exporting = false;
            M_diffusivity = 0;
            break;

        case (variableID::D_del_ci_ridge_myi):
            M_name = "D_del_ci_ridge_myi";
            M_export_name = "MYI_area_ridging_change_daily_rate";
            M_prognostic = false;
            M_exporting = false;
            break;
        
        case (variableID::D_del_vi_mlt_myi):
            M_name = "D_del_vi_mlt_myi";
            M_export_name = "MYI_volume_melt_change_daily_rate";
            M_prognostic = false;
            M_exporting = false;
            break;
        
        case (variableID::D_del_ci_mlt_myi):
            M_name = "D_del_ci_mlt_myi";
            M_export_name = "MYI_area_melt_change_daily_rate";
            M_prognostic = false;
            M_exporting = false;
            break;
        
        case (variableID::D_del_vi_rplnt_myi):
            M_name = "D_del_vi_rplnt_myi";
            M_export_name = "MYI_volume_replenishment_change_daily_rate";
            M_prognostic = false;
            M_exporting = false;
            break;
        
        case (variableID::D_del_ci_rplnt_myi):
            M_name = "D_del_ci_rplnt_myi";
            M_export_name = "MYI_area_replenishment_change_daily_rate";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_albedo):
            // Surface albedo
            M_name = "D_albedo";
            M_export_name = "albedo";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_sialb):
            // Sea ice albedo
            M_name = "D_sialb";
            M_export_name = "sialb";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_dcrit):
            // diagnostic variable to tell use how far outside the
            // Mohr-Coulomb envelope we are
            M_name = "D_dcrit";
            M_export_name = "d_crit";
            M_prognostic = false;
            M_exporting = false;
            break;

        case (variableID::D_sigma_p):
            // stress tensor originally based on the "pressure" term
            // of Rampal et al (2016)
            // - now it is more of a visco-plastic term
            //   - viscous for small deformations
            //   - plastic (independent of strain rate) for large deformations
            M_name = "D_sigma_p";
            M_export_name = "D_sigma_p";
            M_prognostic = false;
            M_exporting = false;
            if(M_component_number<0 || M_component_number>2)
                throw std::runtime_error(
                        "Unauthorised component number for D_sigma_p: "
                        +std::to_string(M_component_number));
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

    std::abort();
}//testCall


} // Nextsim
