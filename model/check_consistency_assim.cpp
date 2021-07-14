void
FiniteElement::checkConsistency_assim()
{
    //check consistency of fields after init/assimilation
    //1. set things to zero if conc/abs thickness are too small
    //2. check SST is consistent
    //3. Initialise M_tice[i]
    for ( int i=0; i<M_num_elements; i++ )
    {
        // freezing points of ice and water needed for init of ice temp
        // and to check SST
        double const Tfr_wtr = this->freezingPoint(M_sss[i]);   //freezing point for water
        double const Tfr_ice = -physical::mu*physical::si;      //freezing point for ice salinity

        // get model variables after assimilation
        double sic_mod(M_analysis_conc[i]);
        double sit_mod(M_analysis_thick[i]);
        double rir_mod(M_analysis_ridge_ratio[i]);
        double snt_mod(M_analysis_snow_thick[i]);
        double sst_mod(M_sst[i]);
        double it0_mod(M_tice_0[i]);
        
        double sic_mod_tot(sic_mod);
        double sit_mod_tot(sit_mod);
        double snt_mod_tot(snt_mod);

        // original variables (before assimilation)
        double sic_orig(M_conc[i]);
        double sic_orig_tot(M_conc[i]);
        double sit_orig_tot(M_thick[i]);
        double snt_orig_tot(M_snow_thick[i]);

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            double sic_mod_thin(M_analysis_conc_thin[i]);
            double sit_new_thin(M_analysis_h_thin[i]);
            double snt_new_thin(M_analysis_hs_thin[i]);

            sic_mod_tot = sic_mod + sic_mod_thin;
            // remove excess of young ice
            double sic_new_thin(sic_mod_thin);
            if (sic_mod_tot > 1)
            {
                sic_new_thin -= (sic_mod_tot - 1);

                // reduce thin ice thickness proportionally to concentration reduction
                double factor(sic_new_thin / sic_mod_thin);
                sit_new_thin *= factor;
                snt_new_thin *= factor;

            }
            // recalculate total conc and thickness
            sic_mod_tot = sic_mod + sic_new_thin;
            sit_mod_tot = sit_mod + sit_new_thin;
            snt_mod_tot = snt_mod + snt_new_thin;

            sic_orig_tot += M_conc_thin[i];
            sit_orig_tot += M_h_thin[i];
            snt_orig_tot += M_hs_thin[i];

            double ist_new_thin(M_analysis_tsurf_thin[i]);
            if ( (sic_orig_tot < physical::cmin) && 
                 (sic_mod_tot >= physical::cmin) )
            {
                ist_new_thin = std::min(ist_new_thin, Tfr_ice);
            }
        }

        // Update ice surface temperature
        // - new ice only
        // - new surface temp is old ice temp, but cannot exceed the melting point
        // set ice temperature for new ice
        double it0_new(it0_mod);
        // if assimilation add new ice
        if ( (sic_orig_tot < physical::cmin) && 
             (sic_mod_tot >= physical::cmin) )
        {
            it0_new = std::min(it0_new, Tfr_ice);
            sst_mod = Tfr_wtr;
        }

        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            double it1_new(M_tice_1[i]);
            double it2_new(M_tice_2[i]);
            // Calculate the temp at the top of the ice
            double Ti = M_tice[0][i];
            if ( (sic_orig < physical::cmin) &&
                 (sic_mod >= physical::cmin) )
            {
                // Calculate the temp of the ice-snow interface (Ti) using a zero-layer model
                // => ki*(Tfr_bot - Ti)/hi = ks(Ti - Ts)/hs
                // => ki*hs*(Tfr_bot - Ti) = ks*hi*(Ti - Ts)
                // => (ki*hs+ks*hi)*Ti = ki*hs*Tfr_bot + ks*hi*Ts
                // => a*Ti = b
                double const a = physical::ki*snt_mod
                    + physical::ks*sit_mod;
                double const b = physical::ki*snt_mod*Tfr_wtr
                    + physical::ks*sit_mod*it0_new;
                Ti = std::min(b/a, Tfr_ice);//make sure it is not higher than freezing point

                // Then use linear interpolation between bottom and top of ice
                it1_new = Tfr_wtr + .75*(Ti - Tfr_wtr);
                it2_new = Tfr_wtr + .25*(Ti - Tfr_wtr);
            }

            // no thick ice
            if ((sic_mod < physical::cmin) ||
                (sit_mod < sic_mod*physical::hmin))
            {
                it1_new = Tfr_ice;
                it2_new = Tfr_ice;
                
            }
        }//Winton

        // no thick ice
        if ((sic_mod < physical::cmin) ||
            (sit_mod < sic_mod*physical::hmin))
        {
            sic_mod = 0;
            sit_mod = 0;
            snt_mod = 0;
            it0_new = Tfr_ice;

        }

        // open water
        if ((sic_mod_tot < physical::cmin) ||
            (sit_mod_tot < sic_mod_tot * physical::hmin) )
        {
            rir_mod = 0;
            if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                // young ice
                sic_new_thin = 0;
                sit_new_thin = 0;
                snt_new_thin = 0;
            }

        }

        M_conc[i] = sic_mod;
        M_thick[i] = sit_mod;
        M_ridge_ratio[i] = rir_mod;
        M_snow_thick[i] = snt_mod;
        M_tice_0[i] = it0_new;
        M_sst[i] = sst_mod;

        if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_tsurf_thin[i] = ist_new_thin;
            M_conc_thin[i] =  sic_new_thin;
            M_h_thin[i] = sit_new_thin;
            M_hs_thin[i] = snt_new_thin;
        }

        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            M_tice_1[i] = it1_new;
            M_tice_2[i] = it2_new;
        }//Winton

        M_conc[i] = std::max(M_conc[i], 0);
        M_thick[i] = std::max(M_thick[i], 0);
        M_ridge_ratio[i] = std::max(M_ridge_ratio[i], 0);
        M_snow_thick[i] = std::max(M_snow_thick[i], 0);
        M_conc_thin[i] = std::max(M_conc_thin[i], 0);
        M_h_thin[i] = std::max(M_h_thin[i], 0);
        M_hs_thin[i] = std::max(M_hs_thin[i], 0);
        M_damage[i] = std::max(M_damage[i], 0.5);
        M_conc_upd[i] = std::max(M_damage[i], -1);
       
        M_conc[i] = std::min(M_conc[i], 1);
        M_conc_thin[i] = std::min(M_conc_thin[i], 1);
        M_ridge_ratio[i] = std::min(M_ridge_ratio[i], 1);
        M_damage[i] = std::min(M_damage[i], 1);
        M_conc_upd[i] = std::min(M_conc_upd[i], 1);
    }

}//checkConsistency_assim


void
FiniteElement::AssimConc()
{
    // Young (thin) ice fraction
    double const _YIF = 0.2;
    // Thickness of new ice
    double const _HNULL = 0.25;

    for ( int i=0; i<M_num_elements; i++ )
    {
        // get model variables
        double sic_mod(M_conc[i]);
        double sit_mod(M_thick[i]);
        double rir_mod(M_ridge_ratio[i]);
        double snt_mod(M_snow_thick[i]);
        double sic_upd_mod(M_conc_upd[i]);

        // set total concentration
        double sic_mod_tot(sic_mod);
        double sit_mod_tot(sic_mod);
        double sic_tot_est(M_analysis_conc[i]);

        // in case of thin (young) ice, total concentration is sum of old and young ice conc
        if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            double sic_mod_thin(M_conc_thin[i]);
            double sit_mod_thin(M_h_thin[i]);
            double snt_mod_thin(M_hs_thin[i]);
            sic_mod_tot += sic_mod_thin;
            sit_mod_tot += sit_mod_thin;
        }

        // add OBSERVED concentration
        double sic_new_tot(sic_mod_tot);
        bool add_ice = ( (sic_tot_est >= 0.15) && (sic_mod_tot < 0.15) &&
                         (sic_mod_tot < sic_tot_est) );

        double sic_added = add_ice ? sic_tot_est - sic_mod_tot : 0.;

        sic_new_tot += sic_added;
        sic_new_tot = sic_tot_est < 0.15 ? 0. : sic_new_tot; // tot conc drops to zero if OBSERVATION says so

        double update_factor = sic_mod_tot < physical::cmin ? 0. : 1.;

        if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            // Split updated total concentration into old and young
            double sic_new(sic_mod);
            double sic_new_thin(sic_mod_thin);
            sic_new_thin += sic_added;
            if (sic_new_tot < physical::cmin)
            {
                sic_new = 0;
                sic_new_thin = 0;
            }

        }
        else
        {
            double sic_new = sic_new_tot < physical::cmin ? 0. : sic_new_tot;
        }   

        // Use code from AssimThick to separate updated SIT
        double sit_new_tot(D_analysis_thick[i]);
        update_factor = (sic_mod_tot >= physical::cmin) ? sit_new_tot/sit_mod_tot : 0;

        // update ice thickness
        double h_thin_new(0);
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            // for two ice categories
            double sit_new(sit_mod);
            double sit_new_thin(sit_mod_thin);
            // where ice was present
            if (sic_mod_tot >= physical::cmin)
            {
                sit_new = sit_mod * update_factor;
                sit_new_thin = sit_mod_thin * update_factor;
            }
            // where new ice was added
            // similar to initialiation:
            // SIT YOUNG is 20% of total SIT until it reaches _HNULL
            // SIT OLDER - remaining part
            if ((sic_mod_tot < physical::cmin) &&
                (sic_new_tot >= physical::cmin))
            {
                sit_new_thin = sit_new_tot * _YIF;
                sit_new_thin = std::min(sit_new_thin, _HNULL); 
                sit_new = sit_new_tot - sit_new_thin;
            }
            
            sit_mod_thin = sit_new_thin;
            sit_mod = sit_new;
            sit_new_tot = sit_mod_thin + sit_mod;
        }
        else
        {   
            // for one ice category
            sit_mod = sit_new_tot;
        }

        // give these value to sit_mod and sit_mod_thin
        if (sic_added > 0)
        {
            h_thin_new = (sit_new_tot - sit_mod_tot)/sic_added;
        }
        // continue AssimConc
        // Update ice thickness, ridge ratio and snow thickness proportionaly to SIC
        // where ice was present:
        double rir_new(rir_mod * update_factor);
        sit_new = sit_mod * update_factor;
        double snt_new(snt_mod * update_factor);
        if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            sit_new_thin = sit_mod_thin * update_factor;
            double snt_new_thin(snt_mod_thin * update_factor);
            sit_new_thin += h_thin_new * sic_added; // all new ice is thin
        }
        else
        {
            sit_new += h_thin_new * sic_added;
        }

        // in brand new ice:
        if ((sic_mod_tot < physical::cmin) &&
            (sic_new_tot >= physical::cmin))
        {
            rir_new = 0;
            snt_new = 0;
            if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                snt_new_thin = 0;
            }
        }

        // compute total ice and snow thickness before and after assimilation
        sit_new_tot = sit_new;
        double snt_new_tot(snt_new);
        // How much concentration was added/removed (positive - concentration added by assimilation)
        double sic_upd_new(sic_new - sic_mod);
        if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            // also add young ice
            sit_new_tot += sit_new_thin;
            snt_new_tot += snt_new_thin;
            sic_upd_new += sic_new_thin - sic_mod_thin;
        }
            
        // weighted average with previous sic_upd
        sic_upd_new = sic_upd_mod * 0.25 + sic_upd_new * 0.75

        M_analysis_conc[i]=sic_new;
        M_analysis_thick[i]=sit_new;
        M_analysis_snow_thick[i]=snt_new;
        M_analysis_ridge_ratio[i]=rir_new;

        if (M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_conc_thin[i]=sic_new_thin;
            M_h_thin[i]=sit_new_thin;
            M_hs_thin[i]=snt_new_thin;
        }

    }
}