
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   meshhandler.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 17 Aug 2022 07:59:06 CEST
 */

#include <variablemanager.hpp>
#include <environment.hpp>

namespace Nextsim
{

//------------------------------------------------------------------------------------------------------
//! Initializes the size of all physical variables with values set to zero.
//! Called by the init() and readRestart() functions.
//! \note
//! - The prefix D_ is used for diagnostic variables (outputs),
//! - The prefix M_ is used for global variables of the finite element class, accessible for all functions defined in finiteelement.cpp.
//! - The suffix M is used for a quantity at the previous (nth) time step (e.g., VTM)
//! - The suffix M is used for a quantity at the second-previous (n-1 th) time step (e.g., VTMM)

void
VariableManager::initVariables(int num_elements, int num_nodes)
{
    //! Global variables are assigned the prefix M_
    this->initModelVariables(num_elements);
    this->assignVariables(num_elements, num_nodes);

}//initVariables

// ==============================================================================
//! Initialise the ModelVariable objects
//! - 1) Instantiate them
//! - 2) Loop over them and sort them into prognostic variables
//!      and ones we will export to binary files
//!
//! \note Steps to create a new elemental variable (nodal var's not implemented yet):
//!     - i) define a variable ID in model_variable.hpp
//!     - ii) define ModelVariable attributes for the variable in model_variable.cpp
//!     - iii) instantiate and add it to M_variables_elt
//!     - iv) initialise it in initModelState() (initIce() or initSlabOcean())
void
VariableManager::initModelVariables(int num_elements)
{

    //! - 1) init all ModelVariable's and put them in M_variables_elt
    // Prognostic variables
    M_conc = ModelVariable(ModelVariable::variableID::M_conc);//! \param M_conc (double) Concentration of thick ice
    M_variables_elt.push_back(&M_conc);
    M_thick = ModelVariable(ModelVariable::variableID::M_thick);//! \param M_thick (double) Thickness of thick ice [m]
    M_variables_elt.push_back(&M_thick);
    M_damage = ModelVariable(ModelVariable::variableID::M_damage);//! \param M_damage (double) Level of damage
    M_variables_elt.push_back(&M_damage);
    M_snow_thick = ModelVariable(ModelVariable::variableID::M_snow_thick);//! \param M_snow_thick (double) Snow thickness (on top of thick ice) [m]
    M_variables_elt.push_back(&M_snow_thick);
    M_ridge_ratio = ModelVariable(ModelVariable::variableID::M_ridge_ratio);//! \param M_ridge_ratio (double) Ratio of ridged vs unridged ice
    M_variables_elt.push_back(&M_ridge_ratio);
    M_conc_upd = ModelVariable(ModelVariable::variableID::M_conc_upd);//! \param M_conc_upd (double) Concentration update by assimilation
    M_variables_elt.push_back(&M_conc_upd);

    const boost::unordered_map<const std::string, setup::ThermoType> str2thermo = boost::assign::map_list_of
        ("zero-layer", setup::ThermoType::ZERO_LAYER)
        ("winton", setup::ThermoType::WINTON);
    setup::ThermoType thermo_type = OptionHandler::getOptionFromMap(vm, "setup.thermo-type", str2thermo);

    switch (thermo_type)
    {
        case (setup::ThermoType::ZERO_LAYER):
            M_tice.resize(1);   //! \param M_tice (double) Ice surface temperature [C]
            break;
        case (setup::ThermoType::WINTON):
            M_tice.resize(3);//! \param M_tice (double) 0: Ice surface temperature [C]; 1,2: Temperatures of upper and lower ice layers [C]
            break;
        default:
            std::cout << "thermo_type= " << (int) thermo_type << "\n";
            throw std::logic_error("Wrong thermo_type");
    }
    for(int k=0; k<M_tice.size(); k++)
    {
        M_tice[k] = ModelVariable(ModelVariable::variableID::M_tice, k);
        M_variables_elt.push_back(&(M_tice[k]));
    }

    M_sigma.resize(3);//! \param M_sigma (double) Tensor components of stress [Pa]
    for(int k=0; k<M_sigma.size(); k++)
    {
        M_sigma[k] = ModelVariable(ModelVariable::variableID::M_sigma, k);
        M_variables_elt.push_back(&(M_sigma[k]));
    }

    M_sst = ModelVariable(ModelVariable::variableID::M_sst);//! \param M_sst (double) Sea surface temperature (slab ocean) [C]
    M_variables_elt.push_back(&M_sst);
    M_sss = ModelVariable(ModelVariable::variableID::M_sss);//! \param M_sss (double) Sea surface salinity (slab ocean) [C]
    M_variables_elt.push_back(&M_sss);
    if ( vm["thermo.newice_type"].as<int>() == 4 )
    {
        M_tsurf_young = ModelVariable(ModelVariable::variableID::M_tsurf_young);//! \param M_tsurf_young (double) Young ice surface temperature [C]
        M_variables_elt.push_back(&M_tsurf_young);
        M_h_young = ModelVariable(ModelVariable::variableID::M_h_young);//! \param M_h_young (double) Thickness of young ice [m]
        M_variables_elt.push_back(&M_h_young);
        M_hs_young = ModelVariable(ModelVariable::variableID::M_hs_young);//! \param M_hs_young (double) Snow thickness (on top of young ice) [m]
        M_variables_elt.push_back(&M_hs_young);
        M_conc_young = ModelVariable(ModelVariable::variableID::M_conc_young);//! \param M_conc (double) Concentration of young ice
        M_variables_elt.push_back(&M_conc_young);
    }
    M_random_number = ModelVariable(ModelVariable::variableID::M_random_number);//! \param M_random_number (double) Random component of cohesion
    M_variables_elt.push_back(&M_random_number);

#ifdef OASIS
    // FSD
    M_conc_fsd.resize(M_num_fsd_bins);
    for(int k=0; k<M_num_fsd_bins; k++)
    {
        M_conc_fsd[k] = ModelVariable(ModelVariable::variableID::M_conc_fsd, k);
        M_variables_elt.push_back(&(M_conc_fsd[k]));
    }
    M_conc_mech_fsd.resize(M_num_fsd_bins);
    for(int k=0; k<M_num_fsd_bins; k++)
    {
        M_conc_mech_fsd[k] = ModelVariable(ModelVariable::variableID::M_conc_mech_fsd, k);
        M_variables_elt.push_back(&(M_conc_mech_fsd[k]));
    }

    M_cum_damage = ModelVariable(ModelVariable::variableID::M_cum_damage);//! \param M_cum_damage (double) Level of accumulated damage (no healing accounted)
    M_variables_elt.push_back(&M_cum_damage);
    M_cum_wave_damage = ModelVariable(ModelVariable::variableID::M_cum_wave_damage);//! \param M_cum_wave_damage (double) Level of accumulated damage (no healing accounted)
    M_variables_elt.push_back(&M_cum_wave_damage);
#endif

    M_fyi_fraction = ModelVariable(ModelVariable::variableID::M_fyi_fraction);//! \param M_fyi_fraction (double) Fraction of FYI
    M_variables_elt.push_back(&M_fyi_fraction);
    M_age_det = ModelVariable(ModelVariable::variableID::M_age_det);//! \param M_age_det (double) Sea ice age observable/detectable from space [s]
    M_variables_elt.push_back(&M_age_det);
    M_age = ModelVariable(ModelVariable::variableID::M_age);//! \param M_age (double) Sea ice age (volumetric) [s]
    M_variables_elt.push_back(&M_age);
    M_conc_myi = ModelVariable(ModelVariable::variableID::M_conc_myi);//! \param M_conc_myi (double) Concentration of MYI
    M_variables_elt.push_back(&M_conc_myi);
    M_thick_myi = ModelVariable(ModelVariable::variableID::M_thick_myi);//! \param M_thick_myi (double) Thickness of MYI
    M_variables_elt.push_back(&M_thick_myi);
    M_freeze_days = ModelVariable(ModelVariable::variableID::M_freeze_days);//! \param M_freeze_days (double) Counter of time (days) of ice freezeing for myi reset
    M_variables_elt.push_back(&M_freeze_days);
    M_conc_summer = ModelVariable(ModelVariable::variableID::M_conc_summer);//! \param M_conc_summer (double) Concentration at end of summer
    M_variables_elt.push_back(&M_conc_summer);
    M_thick_summer = ModelVariable(ModelVariable::variableID::M_thick_summer);//! \param M_thick_summer (double) Thickness at end of summer
    M_variables_elt.push_back(&M_thick_summer);
    M_freeze_onset = ModelVariable(ModelVariable::variableID::M_freeze_onset);//! \param M_freeze_onset (double) Binary for if freeze onset has occurred for this year or not
    M_variables_elt.push_back(&M_freeze_onset);
    M_del_vi_tend = ModelVariable(ModelVariable::variableID::M_del_vi_tend);//! \param M_del_vi_tend (double) Counter of daily total del_hi to deduce melt/freeze day
    M_variables_elt.push_back(&M_del_vi_tend);

    // Diagnostic variables are assigned the prefix D_
    D_conc = ModelVariable(ModelVariable::variableID::D_conc);//! \param D_conc (double) Total concentration of ice
    M_variables_elt.push_back(&D_conc);
    D_thick = ModelVariable(ModelVariable::variableID::D_thick);//! \param D_thick (double) Total thickness of ice [m]
    M_variables_elt.push_back(&D_thick);
    D_snow_thick = ModelVariable(ModelVariable::variableID::D_snow_thick);//! \param D_snow_thick (double) Total snow thickness [m]
    M_variables_elt.push_back(&D_snow_thick);
    D_tsurf = ModelVariable(ModelVariable::variableID::D_tsurf);//! \param D_tsurf (double) Mean surface temperature [C]
    M_variables_elt.push_back(&D_tsurf);
    D_sigma.resize(2);//! \param D_sigma (double) Principal components of stress [Pa]
    for(int k=0; k<D_sigma.size(); k++)
    {
        D_sigma[k] = ModelVariable(ModelVariable::variableID::D_sigma, k);
        M_variables_elt.push_back(&(D_sigma[k]));
    }
    D_divergence = ModelVariable(ModelVariable::variableID::D_divergence);//! \param D_divergence (double) Divergence of ice flow
    M_variables_elt.push_back(&D_divergence);
    D_Qa = ModelVariable(ModelVariable::variableID::D_Qa);//! \param D_Qa (double) Total heat flux to the atmosphere
    M_variables_elt.push_back(&D_Qa);
    D_Qsw = ModelVariable(ModelVariable::variableID::D_Qsw);//! \param D_Qsw (double) Short wave heat flux to the atmosphere
    M_variables_elt.push_back(&D_Qsw);
    D_Qlw = ModelVariable(ModelVariable::variableID::D_Qlw);//! \param D_Qlw (double) Long wave heat flux to the atmosphere
    M_variables_elt.push_back(&D_Qlw);
    D_Qsh = ModelVariable(ModelVariable::variableID::D_Qsh);//! \param D_Qsh (double) Sensible heat flux to the atmosphere
    M_variables_elt.push_back(&D_Qsh);
    D_Qlh = ModelVariable(ModelVariable::variableID::D_Qlh);//! \param D_Qlh (double) Latent heat flux to the atmosphere
    M_variables_elt.push_back(&D_Qlh);
    D_Qo = ModelVariable(ModelVariable::variableID::D_Qo);//! \param D_Qo (double) Total heat lost by the ocean
    M_variables_elt.push_back(&D_Qo);
    D_delS = ModelVariable(ModelVariable::variableID::D_delS);//! \param D_delS (double) Virtual salt flux to the ocean [kg/m^2/s]
    M_variables_elt.push_back(&D_delS);
    D_Qnosun = ModelVariable(ModelVariable::variableID::D_Qnosun);//! \param D_Qnosun (double) Non-solar heat loss from ocean [W/m2]
    M_variables_elt.push_back(&D_Qnosun);
    D_Qsw_ocean = ModelVariable(ModelVariable::variableID::D_Qsw_ocean);//! \param D_Qsw_ocean (double) SW flux out of the ocean [W/m2]
    M_variables_elt.push_back(&D_Qsw_ocean);
    D_vice_melt = ModelVariable(ModelVariable::variableID::D_vice_melt);//! \param D_vice_melt (double) Ice volume formed/melted per element area [m/day]
    M_variables_elt.push_back(&D_vice_melt);
    D_del_vi_young = ModelVariable(ModelVariable::variableID::D_del_vi_young);//! \param D_del_vi_young (double) Young Ice volume formed/melted per element area [m/day]
    M_variables_elt.push_back(&D_del_vi_young);
    D_newice = ModelVariable(ModelVariable::variableID::D_newice);//! \param D_newice (double) Ice volume formed in open water  per element area [m/day]
    M_variables_elt.push_back(&D_newice);
    D_mlt_top = ModelVariable(ModelVariable::variableID::D_mlt_top);//! \param D_mlt_top (double) Ice volume melted at top  per element area [m/day]
    M_variables_elt.push_back(&D_mlt_top);
    D_mlt_bot = ModelVariable(ModelVariable::variableID::D_mlt_bot);//! \param D_mlt_bot (double) Ice volume melted at bottom  per element area [m/day]
    M_variables_elt.push_back(&D_mlt_bot);
    D_snow2ice = ModelVariable(ModelVariable::variableID::D_snow2ice);//! \param D_snow2ice (double) Ice volume formed in from snow flooding per element area [m/day]
    M_variables_elt.push_back(&D_snow2ice);
    D_del_hi_young = ModelVariable(ModelVariable::variableID::D_del_hi_young);//! \param D_del_hi_young (double) Young growth/melt rate [m/day]
    M_variables_elt.push_back(&D_del_hi_young);
    D_del_hi = ModelVariable(ModelVariable::variableID::D_del_hi);//! \param D_del_hi (double) Ice growth/melt rate  [m/day]
    M_variables_elt.push_back(&D_del_hi);
    D_del_vi_mlt_myi = ModelVariable(ModelVariable::variableID::D_del_vi_mlt_myi);//! \param D_del_vi_mlt_myi (double) Ice growth/melt rate  [m/day]
    M_variables_elt.push_back(&D_del_vi_mlt_myi);
    D_del_ci_mlt_myi = ModelVariable(ModelVariable::variableID::D_del_ci_mlt_myi);//! \param D_del_ci_mlt_myi (double) Ice growth/melt rate  [m/day]
    M_variables_elt.push_back(&D_del_ci_mlt_myi);
    D_del_ci_rplnt_myi = ModelVariable(ModelVariable::variableID::D_del_ci_rplnt_myi);//! \param D_del_ci_rplnt_myi (double) Ice growth/melt rate  [m/day]
    M_variables_elt.push_back(&D_del_ci_rplnt_myi);
    D_del_vi_rplnt_myi = ModelVariable(ModelVariable::variableID::D_del_vi_rplnt_myi);//! \param D_del_vi_rplnt_myi (double) Ice growth/melt rate  [m/day]
    M_variables_elt.push_back(&D_del_vi_rplnt_myi);
    D_del_ci_ridge_myi = ModelVariable(ModelVariable::variableID::D_del_ci_ridge_myi);//! \param D_del_ci_ridge_myi (double) Ice growth/melt rate  [m/day]
    M_variables_elt.push_back(&D_del_ci_ridge_myi);
    D_fwflux = ModelVariable(ModelVariable::variableID::D_fwflux);//! \param D_fwflux (double) Fresh-water flux at ocean surface [kg/m2/s]
    M_variables_elt.push_back(&D_fwflux);
    D_fwflux_ice = ModelVariable(ModelVariable::variableID::D_fwflux_ice);//! \param D_fwflux_ice (double) Fresh-water flux at ocean surface due to ice processes [kg/m2/s]
    M_variables_elt.push_back(&D_fwflux_ice);
    D_Qassim = ModelVariable(ModelVariable::variableID::D_Qassim);//! \param D_Qassim (double) flux from assimilation [W/m2]
    M_variables_elt.push_back(&D_Qassim);
    D_brine = ModelVariable(ModelVariable::variableID::D_brine);//! \param D_brine (double) Brine release into the ocean [kg/m2/s]
    M_variables_elt.push_back(&D_brine);
    D_tau_ow = ModelVariable(ModelVariable::variableID::D_tau_ow);//! \param D_tau_ow (double) Ocean atmosphere drag coefficient - still needs to be multiplied with the wind [Pa/s/m] (for the coupled ice-ocean system)
    M_variables_elt.push_back(&D_tau_ow);
    D_evap = ModelVariable(ModelVariable::variableID::D_evap);//! \param D_evap (double) Evaporation from the ocean
    M_variables_elt.push_back(&D_evap);
    D_rain = ModelVariable(ModelVariable::variableID::D_rain);//! \param D_rain (double) Rain into the ocean
    M_variables_elt.push_back(&D_rain);
    D_albedo = ModelVariable(ModelVariable::variableID::D_albedo);//! \param D_albedo (double) Surface albedo - area-weighted average of ocean, young and old ice albedo
    M_variables_elt.push_back(&D_albedo);
    D_sialb = ModelVariable(ModelVariable::variableID::D_albedo);//! \param D_sialb (double) Sea ice albedo - mean albedo where ice
    M_variables_elt.push_back(&D_sialb);

    D_dmax = ModelVariable(ModelVariable::variableID::D_dmax);
    M_variables_elt.push_back(&D_dmax);
    D_dmean = ModelVariable(ModelVariable::variableID::D_dmean);
    M_variables_elt.push_back(&D_dmean);

    //! - 2) loop over M_variables_elt in order to sort them
    //!     for restart/regrid/export
    M_prognostic_variables_elt.resize(0);
    M_export_variables_elt.resize(0);
    M_export_names_elt.resize(0);
    for(auto ptr: M_variables_elt)
    {
        if(ptr->isPrognostic())
        {
            // restart, regrid variables
            M_prognostic_variables_elt.push_back(ptr);
        }
        else if (vm["output.save_diagnostics"].as<bool>())
        {
            // export all diagnostic variables to binary
            ptr->setExporting(true);
        }

        if(ptr->exporting())
        {
            // export variables
            M_export_variables_elt.push_back(ptr);
            M_export_names_elt.push_back(ptr->exportName());
        }

        // resize - note need to make sure variables are initialised properly in initModelState()
        // otherwise model will crash
        ptr->assign(num_elements,0.);
    }// loop over M_variables_elt

    //! - 3) finally sort the prognostic variables into M_prognostic_variables_elt_indices
    //! using ModelVariable::interpTransformation
    //! \note need to do the ModelVariable::interpTransformation::none type variables
    //! first, since the others need M_conc or M_thick
    this->sortPrognosticVars();
}//initModelVariables

//------------------------------------------------------------------------------------------------------
//! Assigns variables in the context of remeshing : the size of variables needs to be update when remeshing because the nb of elements/nodes has changed.
//! Called by the regrid() and initVariables() functions.
void
VariableManager::assignVariables(int num_elements, int num_nodes)
{
    //nodal var's
    M_VT.resize(2*num_nodes); //! \param M_VT (double) Instantaneous velocity vector at the (n+1)th (current) t-step [m/s]
    M_UM.resize(2*num_nodes); //! \param M_UM (double) Total mesh displacement [m]
    D_tau_w.resize(2*num_nodes); //! \param D_tau_w (double) Ice-ocean drag [Pa]
    D_tau_a.resize(2*num_nodes); //! \param D_tau_a (double) Ice-atmosphere drag [Pa]

    M_delta_x.resize(num_nodes);

    M_UM.assign(2*num_nodes,0.);

    // For drifters:
    M_UT.assign(2*num_nodes,0.); //! \param M_UT (double) Total ice displacement (M_UT[] += time_step*M_VT[]) [m]

    // --------------------------------------------------------------

    // //loop over vector of pointers to datasets defined in initExternalData()
    // for (auto it=M_datasets_regrid.begin(), end=M_datasets_regrid.end(); it!=end; ++it)
    // {
    //     if (M_rank == 0)
    //         std::cout<<"REGRIDDING: need to re-interpolate dataset "<<(*it)->name<<"\n";

    //     (*it)->interpolated=false;

    //     // for the parallel code, it will be necessary to add these lines
    //     // as the domain covered by the partitions changes at each remeshing/partitioning
    //     // (*it)->grid.interpolated=false;
    //     (*it)->grid.loaded=false;
    // }

    M_surface.resize(num_elements);
    M_Cohesion.resize(num_elements); // \param M_Cohesion (double) Ice cohesive strength [N/m2]
    const double time_relaxation_damage = vm["dynamics.time_relaxation_damage"].as<double>()*days_in_sec; //! \param time_relaxation_damage (double) Characteristic healing time [s]
    M_time_relaxation_damage.resize(num_elements,time_relaxation_damage); // \param M_time_relaxation_damage (double) Characteristic time for healing [s]

    // root
    // M_UM_root.assign(2*M_mesh.numGlobalNodes(),0.);

}//assignVariables

//------------------------------------------------------------------------------------------------------
//! Sort the prognostic variables by the transformation needed to be done at interpolation (regrid) time.
//! The main reason this is done is to make sure the "none" type is redistributed first,
//! since the "conc" type needs M_conc to be defined, and the "thick" and "enthalpy" types
//! need M_thick. This ordering also ensures collectVariables and redistributeVariables() insert and
//! read variables from the interpolation function in the same order.
//! Called by collectVariables() and redistributeVariables().
void
VariableManager::sortPrognosticVars()
{
    //! - 1) sort variables into indices corresponding to ModelVariable::interpTransformation
    std::vector<std::vector<int>> prognostic_variables_elt_indices(4);
    prognostic_variables_elt_indices.resize(4);
    int nb_var_element = M_prognostic_variables_elt.size();
    for(int j=0; j<M_prognostic_variables_elt.size(); j++)
    {
        auto vptr = M_prognostic_variables_elt[j];
        int k = vptr->getInterpTransformation();
        prognostic_variables_elt_indices[k].push_back(j);
    }

    //! - 2) reorder M_prognostic_variables_elt and:
    //!    * set M_restart_names_elt
    auto tmp_prog_vars = M_prognostic_variables_elt;
    M_prognostic_variables_elt.resize(0);
    M_restart_names_elt.resize(0);
    for (auto inds : prognostic_variables_elt_indices)
        for (int j : inds)
        {
            auto vptr = tmp_prog_vars[j];//this is a pointer to a ModelVariable object
            M_prognostic_variables_elt.push_back(vptr);
            M_restart_names_elt.push_back(vptr->name());
        }
}//sortPrognosticVars

}
