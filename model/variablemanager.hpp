/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   variablemanager.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 17 Aug 2022 07:59:06 CEST
 */

#ifndef __VariableManager_HPP
#define __VariableManager_HPP 1

#include <model_variable.hpp>
#include <boost/program_options.hpp>
#include <enums.hpp>
#include <optionhandler.hpp>

namespace Nextsim
{

    class VariableManager
    {
        // Methods
        public:
            VariableManager() : vm(Environment::vm()) {}

            void initVariables(int num_elements, int num_nodes);
            void assignVariables(int num_elements, int num_nodes);

        private:

            void initModelVariables(int num_elements);
            void sortPrognosticVars();

        // Variables
        protected:
            //Nodal variables
            std::vector<double> M_delta_x;

            std::vector<double> M_UM;
            std::vector<double> M_UT;
            std::vector<double> M_VT;

            std::vector<double> D_tau_w; // Ice-ocean drag [Pa]
            std::vector<double> D_tau_a; // Ice-atmosphere drag [Pa]

            // Elemental variables
            ModelVariable M_conc;               // Ice concentration
            ModelVariable M_thick;              // Effective ice thickness [m]
            ModelVariable M_damage;             // Ice damage
#ifdef OASIS
            ModelVariable M_cum_damage;         // Ice cumulated damage
            ModelVariable M_cum_wave_damage;    // Ice cumulated damage due to wave
#endif
            ModelVariable M_snow_thick;         // Effective snow thickness [m]
            ModelVariable M_ridge_ratio;
            std::vector<ModelVariable> M_tice;  // Ice temperature - 0 for surface and higher ordinals for layers in the ice
            std::vector<ModelVariable> M_sigma; // Internal stress tensor
            ModelVariable M_sst;                // Sea-surface (slab ocean) temperature [C]
            ModelVariable M_sss;                // Sea-surface (slab ocean) salinity [psu]
            ModelVariable M_tsurf_young;        // Ice surface temperature of young ice [C]
            ModelVariable M_h_young;
            ModelVariable M_hs_young;
            ModelVariable M_conc_young;
            ModelVariable M_random_number;
            ModelVariable M_fyi_fraction;
            ModelVariable M_age_det;
            ModelVariable M_age;
            ModelVariable M_conc_upd;           // Ice concentration update by assimilation
            ModelVariable M_conc_myi;           // Concentration of multiyear ice
            ModelVariable M_thick_myi;          // Thickness of multiyear ice
            ModelVariable M_freeze_days;        // Days of consecutive freezing
            ModelVariable M_freeze_onset;       // Onset of freezing
            ModelVariable M_conc_summer;        // Concentration at end of summer
            ModelVariable M_thick_summer;       // Thickness at end of summer
            ModelVariable M_del_vi_tend;        // Daily sum of ice volume tendency

#ifdef OASIS
            // Following variables are related to floe size distribution
            std::vector<ModelVariable> M_conc_fsd;
            //std::vector<ModelVariable> M_conc_fsd_thick;
            //std::vector<ModelVariable> M_conc_fsd_young ;
            std::vector<ModelVariable> M_conc_mech_fsd;
            int M_num_fsd_bins;
            std::vector<double> M_fsd_bin_widths;
            double M_fsd_bin_cst_width;
            double M_fsd_min_floe_size;
            std::vector<double> M_fsd_bin_centres;
            std::vector<double> M_fsd_bin_low_limits;
            std::vector<double> M_fsd_bin_up_limits;
            double M_fsd_unbroken_floe_size     ;
            // Non-circularity of floes
            double M_floe_shape                              ;
            // Lettie's variables
            std::vector<double> M_floe_area_up             ;
            std::vector<double> M_floe_area_low            ;
            std::vector<double> M_floe_area_centered       ;
            std::vector<double> M_floe_area_binwidth       ;

            std::vector<double> M_fsd_area_scaled_up       ;
            std::vector<double> M_fsd_area_scaled_low      ;
            std::vector<double> M_fsd_area_scaled_centered ;
            std::vector<double> M_fsd_area_scaled_binwidth ;
            std::vector<double> M_fsd_area_lims            ;
            std::vector<double> M_fsd_area_lims_scaled     ;

            std::vector<std::vector<int> > M_alpha_fsd_merge ;
            // In namelist
            bool   M_distinguish_mech_fsd             ;
            bool   M_debug_fsd                        ;
            int    M_fsd_damage_type                  ;
            double M_floes_flex_strength              ;
            double M_floes_flex_young                 ;
            double M_welding_kappa                    ;
            bool   M_fsd_welding_use_scaled_area      ;
            double M_dmax_c_threshold                 ;
            double M_breakup_thick_min                ;
            bool   M_breakup_in_dt                    ;
            bool   M_breakup_cell_average_thickness   ;
            // Horvat et Tziperman (2015) lead fraction, lat. surf and lead width
            // double M_lead_width    ;
#endif


            // Diagnostic variables
            ModelVariable D_conc; //total concentration
            ModelVariable D_thick; // total thickness [m]
            ModelVariable D_snow_thick;// total snow thickness [m]
            ModelVariable D_tsurf; //mean surface temperature (thick + young ice + slab ocean) [deg C]
            std::vector<ModelVariable> D_sigma; //principal stresses [Pa]
            ModelVariable D_divergence; // ice velocity divergence [1/s]
            ModelVariable D_Qa; // Heat loss to atmosphere [W/m2]
            ModelVariable D_Qsw; // Total short wave at surface [W/m2]
            ModelVariable D_Qlw; // Total long wave at surface [W/m2]
            ModelVariable D_Qsh; // Total sensible heat flux at surface [W/m2]
            ModelVariable D_Qlh; // Total latent heat flux at surface [W/m2]
            ModelVariable D_Qo; // Heat loss from ocean [W/m2]
            ModelVariable D_Qnosun; // Non-solar heat loss from ocean [W/m2]
            ModelVariable D_Qsw_ocean; // SW flux out of the ocean [W/m2]
            ModelVariable D_Qassim; // flux from assim [W/m2]
            ModelVariable D_vice_melt; // ice volume (/element_area) melted/formed [m/day]
            ModelVariable D_del_vi_young; // young ice volume (/element_area) melted/formed [m/day]
            ModelVariable D_del_hi; // ice growth/melt rate [m/sdat
            ModelVariable D_del_hi_young; // young ice growth/melt rate [m/day]
            ModelVariable D_newice; // ice volume (/element_area) formed in open water [m/day]
            ModelVariable D_mlt_top; // ice volume (/element_area) melted at top [m/day]
            ModelVariable D_mlt_bot; // ice volume (/element_area) melted at bottom [m/day]
            ModelVariable D_del_vi_mlt_myi;   //myi ice volume (/element_area) melted [m/day]
            ModelVariable D_del_vi_rplnt_myi;  //myi ice vol change (/element_area) due to replenishment  [./day]
            ModelVariable D_del_ci_rplnt_myi;  //myi ice area change (/element_area) due to replenishment  [./day]
            ModelVariable D_del_ci_mlt_myi;    //myi ice area (/element_area) melted   [m/day]
            ModelVariable D_del_ci_ridge_myi; //myi ice area change (/element_area) due to ridging [./day]
            ModelVariable D_snow2ice; // ice volume (/element_area) melted at bottom [m/day]
            ModelVariable D_delS; // Salt flux to ocean
            ModelVariable D_fwflux; // Fresh-water flux at ocean surface [kg/m2/s]
            ModelVariable D_fwflux_ice; // Fresh-water flux at ocean surface due to ice processes [kg/m2/s]
            ModelVariable D_brine; // Brine release into the ocean [kg/m2/s]
            ModelVariable D_dmax; //max floe size [m]
            ModelVariable D_dmean; //mean floe size [m]
            ModelVariable D_tau_ow; // Ocean atmosphere drag coefficient - still needs to be multiplied with the wind [Pa/s/m] (for the coupled ice-ocean system)
            ModelVariable D_evap; // Evaporation out of the ocean [kg/m2/s]
            ModelVariable D_rain; // Rain into the ocean [kg/m2/s]
            ModelVariable D_albedo; // surface albedo
            ModelVariable D_sialb; // sea ice albedo

            std::vector<double> M_surface;
            std::vector<double> M_Cohesion;
            std::vector<double> M_time_relaxation_damage;

            // vectors of pointers to variables (for looping)
            std::vector<ModelVariable*> M_variables_elt;
            std::vector<ModelVariable*> M_prognostic_variables_elt;//for restart, regrid
            std::vector<ModelVariable*> M_export_variables_elt;

            // other vectors related to export/restart
            std::vector<std::string> M_restart_names_elt;
            std::vector<std::string> M_export_names_elt;

        private:
            po::variables_map vm;
            double const days_in_sec  = 86400.;
    };
}

#endif
