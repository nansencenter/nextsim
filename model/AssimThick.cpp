void FiniteElement::AssimThick()
{   
    // post processing of the update ice thickness from enkf-c
    // refer to pynextsimf/assimilation.py->AssimThick()
    // Assimilate sea ice thickness from CS2SMOS  
    // Where ice in nextsim is present:
    //     *  new total thickness is calculated using enkf-c
    //     *  thick and thin/young ice thickness is changed proportionally
    // Young (thin) ice fraction
    double const _YIF = 0.2;
    // Thickness of new ice
    double const _HNULL = 0.25;
    double update_factor;
    double sit_mod, sic_mod, snt_mod, rir_mod, sit_mod_thin, sic_mod_thin, snt_mod_thin;
    double sit_new, sic_new, snt_new, rir_new, sit_new_thin, sic_new_thin, snt_new_thin;
    bool ice00,ice01,ice10,ice11;

    for ( int i=0; i<M_num_elements; i++ )
    {
        //
        sit_mod     = M_thick[i];
        sit_mod_tot = sit_mod;
        sit_new_tot = D_analysis_thick[i];
        //  
        sic_mod     = M_conc[i];
        sic_mod_tot = sic_mod;
        sic_new_tot = D_analysis_conc[i];
        //
        snt_mod = M_snow_thick[i];
        // Will we use analysis for rir?
        rir_mod = M_analysis_ridge_ratio[i];

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            sit_mod_thin= M_h_thin[i];
            sic_mod_thin= M_conc_thin[i];
            snt_mod_thin = M_hs_thin[i];
            sit_mod_tot += sit_mod_thin;
            sic_mod_tot += sic_mod_thin;
        }

        // codes identifying where ice was present before and after assim
        ice00 = (sit_mod_tot < physical::hmin) && (sit_new_tot < physical::hmin);  // is the grammar correct ?
        ice01 = (sit_mod_tot < physical::hmin) && (sit_new_tot >= physical::hmin);
        ice10 = (sit_mod_tot >= physical::hmin) && (sit_new_tot < physical::hmin);
        ice11 = (sit_mod_tot >= physical::hmin) && (sit_new_tot >= physical::hmin);

        // calculate update factor 
        update_factor = (ice10 || ice11) ? sit_new_tot/sit_mod_tot : 0;

        // update ice thickness 
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            // for two ice categories
            sit_new = sit_mod;
            sit_new_thin = sit_mod_thin;
            // where ice was present
            if (ice10 || ice11)
            {
                sit_new = sit_mod * update_factor;
                sit_new_thin = sit_mod_thin * update_factor;           
            }
            // where new ice was added
            // similar to initialiation:
            // SIT YOUNG is 20% of total SIT until it reaches _HNULL
            // SIT OLDER - remaining part
            if (ice01)
            {
                sit_new_thin = sit_new_tot * _YIF;
                sit_new_thin = std::min(sit_new_thin, _HNULL); 
                sit_new = sit_new_tot - sit_new_thin;
            }            
        }
        else
        {   
            // for one ice category
            sit_new = sit_new_tot;
        }

        // tune other ice properties based on new ice thickness (sit_new, sit_new_thin) obtained above
        // REMOVE ice where it has disappeared (total thickness below threshold)
        sic_new = sic_mod;
        snt_new = snt_mod;
        rir_new = rir_mod;
        if(sit_new < _HMIN){
            sit_new = 0;
            sic_new = 0;
            snt_new = 0;
            rir_new = 0;           
        }
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            sic_new_thin = sic_mod_thin;
            snt_new_thin = snt_mod_thin;
            if (sit_new_thin < _HMIN)
            {
                sit_new_thin = 0;
                sic_new_thin = 0;
                snt_new_thin = 0;
            }
        }

        // // increase SIC where ice appeared: linear interpolation from ice to no ice
        // new_ice_mask = (sic_mod < _CMIN)*(sic_new >0)
        // sic_new = self.interpolate_gap(sic_new, new_ice_mask)
        // if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        //     new_thin_ice_mask = (sic_mod_thin < _CMIN)*(sic_new_thin>0)
        //     sic_new_thin = self.interpolate_gap(sic_new_thin, new_thin_ice_mask)

        // update variables
        M_analysis_conc[i]=sic_new;
        M_analysis_thick[i]=sit_new;
        M_analysis_snow_thick[i]=snt_new;
        M_analysis_ridge_ratio[i]=rir_new;
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE){
            M_analysis_conc_thin[i]=sic_new_thin;
            M_analysis_h_thin[i]=sit_new_thin;
            M_analysis_hs_thin[i]=snt_new_thin;
        }
    }
}