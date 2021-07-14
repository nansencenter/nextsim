void FiniteElement::AssimThick()
{   
    // post processing of the update ice thickness from enkf-c
    // refer to pynextsimf/assimilation.py->AssimThick()
    // Assimilate sea ice thickness from CS2SMOS  
    // Where ice in nextsim is present:
    //     *  new total thickness is calculated using enkf-c
    //     *  thick and thin/young ice thickness is changed proportionally

    // Minimum ice thickness allowed [m]
    double _HMIN = 0.01;
    // Young (thin) ice fraction
    _YIF = 0.2;
    // Thickness of new ice
    _HNULL = 0.25;
    double update_factor;
    double sit_mod, sic_mod, snt_mod, rir_mod, sit_mod_thin, sic_mod_thin, snt_mod_thin;
    double sit_new, sic_new, snt_new, rir_new, sit_new_thin, sic_new_thin, snt_new_thin;
    int ice00,ice01,ice10,ice11;

    for ( int i=0; i<M_num_elements; i++ )
    {
        //
        sit_mod     = M_thick[i];
        sit_mod_thin= M_h_thin[i];
        sit_mod_tot = D_thick[i];        
        sit_new_tot = D_analysis_thick[i];
        //  
        sic_mod     = M_conc[i];          
        sic_mod_thin= M_conc_thin[i];
        sic_mod_tot = D_conc[i];      
        sic_new_tot = D_analysis_conc[i];
        //
        snt_mod = M_snow_thick[i];
        snt_mod_thin = M_hs_thin[i];
        rir_mod = M_analysis_ridge_ratio[i];

        // codes identifying where ice was present before and after assim
        ice00 = (sit_mod_tot < _HMIN) * (sit_new_tot < _HMIN);  // is the grammar correct ?
        ice01 = (sit_mod_tot < _HMIN) * (sit_new_tot >= _HMIN);
        ice10 = (sit_mod_tot >= _HMIN) * (sit_new_tot < _HMIN);
        ice11 = (sit_mod_tot >= _HMIN) * (sit_new_tot >= _HMIN);

        // calculate update factor 
        if(ice10 + ice11>0) //(sit_mod_tot>=_HMIN)
            update_factor = sit_new_tot/sit_mod_tot;

        // update ice thickness 
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {// for two ice categories
            sit_new = sit_mod;
            sit_new_thin = sit_mod_thin;
            // where ice was present
            if (ice10 + ice11>0){
                sit_new = sit_mod * update_factor;
                sit_new_thin = sit_mod_thin * update_factor;           
            }
            // where new ice was added
            // similar to initialiation:
            // SIT YOUNG is 20% of total SIT until it reaches _HNULL
            // SIT OLDER - remaining part
            if(ice01>0){
                sit_new_thin = sit_new_tot * _YIF;
                sit_new_thin = std::min(sit_new_thin, _HNULL); 
                sit_new = sit_new_tot - sit_new_thin;
            }            
        }
        else{// for one ice category
            sit_mod = sit_new_tot;
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
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE){
            sic_new_thin = sic_mod_thin
            snt_new_thin = snt_mod_thin
            if (sit_new_thin < _HMIN){
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
        M_conc[i]=sic_new;
        M_thick[i]=sit_new;
        M_snow_thick[i]=snt_new;
        M_ridge_ratio[i]=rir_new;
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE){
            M_conc_thin[i]=sic_new_thin;
            M_h_thin[i]=sit_new_thin;
            M_hs_thin[i]=snt_new_thin;
        }
    }
}