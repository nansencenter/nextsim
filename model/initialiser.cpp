//
// Created by Einar Ólason on 15/09/2022.
//

#include "initialiser.hpp"

namespace Nextsim
{

    void
    Initialiser::initOptAndParams()
    {
        const boost::unordered_map<const std::string, setup::IceType> str2conc = boost::assign::map_list_of
                ("constant", setup::IceType::CONSTANT)
                ("constant_partial", setup::IceType::CONSTANT_PARTIAL)
                ("topaz", setup::IceType::TOPAZ4)
                ("topaz_forecast", setup::IceType::TOPAZ4F)
                ("topaz_forecast_amsr2", setup::IceType::TOPAZ4FAMSR2)
                ("topaz_forecast_amsr2_osisaf", setup::IceType::TOPAZ4FAMSR2OSISAF)
                ("topaz_forecast_amsr2_osisaf_nic", setup::IceType::TOPAZ4FAMSR2OSISAFNIC)
                ("topaz_forecast_amsr2_osisaf_nic_weekly", setup::IceType::TOPAZ4FAMSR2OSISAFNICWEEKLY)
                ("amsre", setup::IceType::AMSRE)
                ("amsr2", setup::IceType::AMSR2)
                ("amsr2_cst_thick", setup::IceType::AMSR2CSTTHICK)
                ("piomas", setup::IceType::PIOMAS)
                ("nemo", setup::IceType::NEMO)
                ("cice", setup::IceType::CICE)
                ("cs2_smos", setup::IceType::CS2_SMOS)
                ("cs2_smos_amsr2", setup::IceType::CS2_SMOS_AMSR2)
                ("smos", setup::IceType::SMOS)
                ("topaz_osisaf_icesat", setup::IceType::TOPAZ4OSISAFICESAT)
                ("glorys12", setup::IceType::GLORYS12R);
        M_ice_type = OptionHandler::getOptionFromMap(vm, "setup.ice-type", str2conc);
        //! \param M_ice_type (enum) Option on the type of ice initialisation
        LOG(DEBUG) <<"IceType= "<< (int)M_ice_type <<"\n";

    }

//------------------------------------------------------------------------------------------------------
//! Initializes the physical state of the model - ocean and ice - by calling the initSlabOcean() and initIce() functions.
//! Called by the init() function.
    void
    Initialiser::initModelState()
    {

        LOG(DEBUG) << "initSlabOcean\n";
        this->initSlabOcean();

        LOG(DEBUG) << "initIce\n";
        this->initIce();
    }//initModelState

//------------------------------------------------------------------------------------------------------
//! Initializes variables relevant to a slab ocean (sss and sst).
//! Called by the initModelState() function.
    void
    Initialiser::initSlabOcean()
    {
        switch (M_ocean_type)
        {
            case setup::OceanType::CONSTANT:
                //std::fill(M_sst.begin(), M_sst.end(), -1.8);
                std::fill(M_sst.begin(), M_sst.end(), 1.);
                std::fill(M_sss.begin(), M_sss.end(),  1.8/physical::mu);
                break;
            case setup::OceanType::TOPAZR:
            case setup::OceanType::TOPAZR_atrest:
            case setup::OceanType::TOPAZF:
            case setup::OceanType::TOPAZR_ALTIMETER:
            case setup::OceanType::GLORYS12R:
                for ( int i=0; i<M_num_elements; ++i)
                {
                    // Make sure the erroneous salinity and temperature don't screw up the initialisation too badly
                    // This can still be done much better!
                    M_sss[i] = std::max(physical::si, M_ocean_salt[i]);
                    //M_sst[i] = std::max(this->freezingPoint(M_sss[i]), M_ocean_temp[i]);
                    M_sst[i] = std::max(-physical::mu*M_sss[i], M_ocean_temp[i]);
                }

                break;
#if OASIS
                case setup::OceanType::COUPLED:
            for ( int i=0; i<M_num_elements; ++i)
            {
                M_sss[i] = M_ocean_salt[i];
                M_sst[i] = M_ocean_temp[i];
            }
            break;
#endif
            default:
                std::cout << "invalid ocean initialisation"<<"\n";
                throw std::logic_error("invalid ocean forcing");
        }
    }//initSlabOcean

//------------------------------------------------------------------------------------------------------
//! Initializes the ice state (CONSTANT, TOPAZ4, PIOMAS, SMOS, ...).
//! Called by the initModelState() function.
    void
    Initialiser::initIce()
    {
        switch (M_ice_type)
        {
            case setup::IceType::CONSTANT:
            case setup::IceType::CONSTANT_PARTIAL:
                this->constantIce();
                break;
            case setup::IceType::TOPAZ4:
                this->topazIce();
                break;
            case setup::IceType::TOPAZ4OSISAFICESAT:
                this->topazIceOsisafIcesat();
                break;
            case setup::IceType::TOPAZ4F:
                this->topazForecastIce();
                break;
            case setup::IceType::TOPAZ4FAMSR2:
                this->topazForecastAmsr2Ice();
                break;
            case setup::IceType::TOPAZ4FAMSR2OSISAF:
                this->topazForecastAmsr2OsisafIce();
                break;
            case setup::IceType::TOPAZ4FAMSR2OSISAFNIC:
            case setup::IceType::TOPAZ4FAMSR2OSISAFNICWEEKLY:
                this->topazForecastAmsr2OsisafNicIce(M_ice_type==setup::IceType::TOPAZ4FAMSR2OSISAFNICWEEKLY);
                break;
            case setup::IceType::PIOMAS:
                this->piomasIce();
                break;
            case setup::IceType::NEMO:
                this->nemoIce();
                break;
            case setup::IceType::CICE:
                this->ciceIce();
                break;
            case setup::IceType::AMSRE:
                this->topazAmsreIce();
                break;
            case setup::IceType::AMSR2:
                this->topazAmsr2Ice();
                break;
            case setup::IceType::AMSR2CSTTHICK:
                this->amsr2ConstThickIce();
                break;
            case setup::IceType::CS2_SMOS:
                this->cs2SmosIce();
                break;
            case setup::IceType::CS2_SMOS_AMSR2:
                this->cs2SmosAmsr2Ice();
                break;
            case setup::IceType::SMOS:
                this->smosIce();
                break;
            case setup::IceType::GLORYS12R:
                this->glorys12Ice();
                break;
            default:
                std::cout << "invalid initialization of the ice"<<"\n";
                throw std::logic_error("invalid initialization of the ice");
        }

        // random numbers
        // - 1st set on root
        std::vector<double> random_number_root(M_mesh.numGlobalElements());
        if (M_comm.rank() == 0)
        {
            boost::random::minstd_rand intgen;
            boost::random::uniform_01<boost::random::minstd_rand> gen(intgen);

            for (int i=0; i<random_number_root.size(); ++i)
                random_number_root[i] = gen();
        }
        boost::mpi::broadcast(M_comm, &random_number_root[0], M_mesh.numGlobalElements(), 0);

        // - now set on each processor
        auto id_elements = M_mesh.trianglesIdWithGhost();
        for (int i=0; i<M_random_number.size(); ++i)
            M_random_number[i] = random_number_root[id_elements[i]-1];

        // init ice age variables and stress to 0
        std::vector<ModelVariable*> vars_to_zero;
        vars_to_zero.push_back(&M_fyi_fraction);
        vars_to_zero.push_back(&M_age_det);
        vars_to_zero.push_back(&M_age);
        vars_to_zero.push_back(&M_conc_upd);
        vars_to_zero.push_back(&M_thick_myi);
        vars_to_zero.push_back(&M_conc_myi);
        vars_to_zero.push_back(&M_freeze_days);
        vars_to_zero.push_back(&M_thick_summer);
        vars_to_zero.push_back(&M_conc_summer);
        vars_to_zero.push_back(&M_freeze_onset);
        vars_to_zero.push_back(&M_del_vi_tend);
        for (int k=0; k<3; k++)
            vars_to_zero.push_back(&(M_sigma[k]));
        for (auto ptr: vars_to_zero)
            std::fill(ptr->begin(), ptr->end(), 0.);

        // init velocities, displacements to 0.
        std::vector<std::vector<double>*> vecs_to_zero;
        vecs_to_zero.push_back(&M_VT);
        vecs_to_zero.push_back(&M_UT);
        vecs_to_zero.push_back(&M_UM);
        for (auto ptr: vecs_to_zero)
            std::fill(ptr->begin(), ptr->end(), 0.);

        // check consistency of fields after initialisation
        // - init ice temp everywhere
        // - init FSD
        // this->checkConsistency(); // Now done in FiniteElement
    }//initIce

//------------------------------------------------------------------------------------------------------
//! Sets the ice cover to a homogeneous state.
//! Called by the initIce() function.
    void
    Initialiser::constantIce()
    {
        LOG(DEBUG) <<"Constant Ice\n";
        double c_const = vm["ideal_simul.init_concentration"].as<double>();
        double h_const = vm["ideal_simul.init_thickness"].as<double>();
        double hs_const = vm["ideal_simul.init_snow_thickness"].as<double>();
        std::fill(M_conc.begin(), M_conc.end(), c_const);
        std::fill(M_thick.begin(), M_thick.end(), c_const*h_const);//M_thick is ice volume
        std::fill(M_snow_thick.begin(), M_snow_thick.end(), hs_const);

        std::vector<ModelVariable*> vars_to_zero;
        vars_to_zero.push_back(&M_damage);
        vars_to_zero.push_back(&M_ridge_ratio);
        for (auto ptr: vars_to_zero)
            std::fill(ptr->begin(), ptr->end(), 0.);

        // No ice where SST is over init_SST_limit
        double SST_limit = vm["ideal_simul.init_SST_limit"].as<double>();
        double init_young_conc = vm["ideal_simul.init_young_conc"].as<double>();
        double h_young_min = vm["thermo.h_young_min"].as<double>();
        double h_young_max = vm["thermo.h_young_max"].as<double>();
        int cnt=0;
        bool const use_thermo = vm["thermo.use_thermo_forcing"].as<bool>();
        for (int i=0; i<M_sst.size(); ++i)
        {
            bool set_young = vm["thermo.newice_type"].as<int>()==4;
            if ( use_thermo
                 && M_sst[i] > -physical::mu*M_sss[i] + SST_limit )
                 // && M_sst[i] > this->freezingPoint(M_sss[i]) + SST_limit )
            {
                M_conc[i]       = 0;
                M_thick[i]      = 0;
                M_snow_thick[i] = 0;
                cnt++;
                set_young = false;
            }
            if(set_young)
            {
                M_conc_young[i] = init_young_conc;
                M_h_young[i]    = (h_young_min+(h_young_max-h_young_min)/2.)*M_conc_young[i];
                M_hs_young[i]   = hs_const*M_conc_young[i];
            }
        }
        LOG(DEBUG) << (double)cnt/(double)M_sst.size() * 100 << "% ice covered cells cleared because of SST limit\n";

        if (M_ice_type==setup::IceType::CONSTANT_PARTIAL)
        {
            double xmin, xmax;
            if(M_comm.rank()==0)
            {
                auto Bx = M_mesh_root.coordX();//xmin,xmax from nodes of global mesh
                xmin = *std::min_element(Bx.begin(), Bx.end());
                xmax = *std::max_element(Bx.begin(), Bx.end());
            }
            boost::mpi::broadcast(M_comm, xmin, 0);
            boost::mpi::broadcast(M_comm, xmax, 0);
            double xedge = xmin + 0.3*(xmax-xmin);

            LOG(DEBUG)<<"In constantIce (partial cover)\n";
            LOG(DEBUG)<<"M_ice_type "<< (int)M_ice_type<<"\n";
            LOG(DEBUG)<<"Min conc = "<< *std::min_element(M_conc.begin(),M_conc.end()) <<"\n";
            LOG(DEBUG)<<"Max conc = "<< *std::max_element(M_conc.begin(),M_conc.end()) <<"\n";
            LOG(DEBUG)<<"Min thick = "<< *std::min_element(M_thick.begin(),M_thick.end()) <<"\n";
            LOG(DEBUG)<<"Max thick = "<< *std::max_element(M_thick.begin(),M_thick.end()) <<"\n";
            LOG(DEBUG)<<"xmin="<<xmin<<"\n";
            LOG(DEBUG)<<"xmax="<<xmax<<"\n";
            LOG(DEBUG)<<"xedge="<<xedge<<"\n";

            auto Bx = M_mesh.bCoordX();//set conc, etc on elements
            for (int i=0; i<M_conc.size(); ++i)
            {
                if (Bx[i] < xedge)
                {
                    M_conc[i]       = 0.;
                    M_thick[i]      = 0.;
                    M_snow_thick[i] = 0.;
                }
            }
            LOG(DEBUG)<<"New min conc = "<< *std::min_element(M_conc.begin(),M_conc.end()) <<"\n";
            LOG(DEBUG)<<"New max conc = "<< *std::max_element(M_conc.begin(),M_conc.end()) <<"\n";
            LOG(DEBUG)<<"New min thick = "<< *std::min_element(M_thick.begin(),M_thick.end()) <<"\n";
            LOG(DEBUG)<<"New max thick = "<< *std::max_element(M_thick.begin(),M_thick.end()) <<"\n";
            //std::abort();
        }//partial ice cover

    }//constantIce

//------------------------------------------------------------------------------------------------------
//! Initializes the ice and snow states from Topaz outputs.
//! Called by the initIce() function.
    void
    Initialiser::topazIce()
    {
        external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
        external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
        external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_init_conc);
        external_data_tmp.push_back(&M_init_thick);
        external_data_tmp.push_back(&M_init_snow_thick);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - TOPAZ ice ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);
        external_data_tmp.resize(0);

        double tmp_var;
        for (int i=0; i<M_num_elements; ++i)
        {
            // - TOPAZ puts very small values instead of 0.
            // - uses absolute thickness not effective thickness
            tmp_var=std::min(1.,M_init_conc[i]);
            M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.;
            tmp_var=M_init_thick[i];
            M_thick[i] = (tmp_var>1e-14) ? tmp_var*M_conc[i] : 0.; // TOPAZ puts very small values instead of 0.
            tmp_var=M_init_snow_thick[i];
            M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var*M_conc[i] : 0.; // TOPAZ puts very small values instead of 0.

            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }

            if(M_thick[i]<=0.)
            {
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }

            M_damage[i]=0.;
        }
    }//topazIce

// -----------------------------------------------------------------------------------------------------------
//! Initializes ice state from Topaz and Osisaf data.
//! Called by the initIce() function.
    void
    Initialiser::topazIceOsisafIcesat()
    {
        //topaz
        external_data M_topaz_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
        external_data M_topaz_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
        external_data M_topaz_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);

        //obs
        external_data M_osisaf_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);
        external_data M_osisaf_conc=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,0,false,time_init);
        external_data M_icesat_thick=ExternalData(&M_ice_icesat_elements_dataset,M_mesh,0,false,time_init);
        external_data M_amsre_conc=ExternalData(&M_ice_amsre_elements_dataset,M_mesh,0,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_topaz_conc);
        external_data_tmp.push_back(&M_topaz_thick);
        external_data_tmp.push_back(&M_topaz_snow_thick);
        external_data_tmp.push_back(&M_osisaf_type);
        external_data_tmp.push_back(&M_osisaf_conc);
        external_data_tmp.push_back(&M_icesat_thick);
        external_data_tmp.push_back(&M_amsre_conc);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - TOPAZ/OSISAF/Icesat ice ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        const double SST_limit = vm["ideal_simul.init_SST_limit"].as<double>();
        const double h_young_min = vm["thermo.h_young_min"].as<double>();
        const double h_young_max = vm["thermo.h_young_max"].as<double>();
        for (int i=0; i<M_num_elements; ++i)
        {
            // - TOPAZ puts very small values instead of 0.
            // - uses absolute thickness not effective thickness
            double tmp_var=std::min(1.,M_topaz_conc[i]);
            M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.;
            tmp_var=M_topaz_thick[i];
            double hi = (tmp_var>1e-14) ? tmp_var : 0.;
            tmp_var=M_topaz_snow_thick[i];
            double hs = (tmp_var>1e-14) ? tmp_var : 0.;

            if(M_conc[i]>0.) // use osisaf only where topaz says there is ice to avoid near land issues and fake concentration over the ocean
                M_conc[i]=M_osisaf_conc[i];

            //M_type[i]==1. // No ice
            //M_type[i]==2. // First-Year ice
            //M_type[i]==3. // Multi-Year ice
            //M_type[i]==4. // Mixed
            double ratio_FYI=0.3;
            double ratio_MYI=0.9;
            double ratio_Mixed=0.5*(ratio_FYI+ratio_MYI);

            double thick_FYI=hi;
            double thick_MYI=std::max(M_icesat_thick[i],hi);//NB icesat outputs absolute thickness
            double thick_Mixed=0.5*(thick_FYI+thick_MYI);

            if((hi>0.)&&(M_conc[i])>0.2)
            {
                if(M_mesh_filename.find("kara") != std::string::npos)
                {
                    LOG(DEBUG) <<"Type information is not used for the kara mesh, "
                               <<"we assume there is only FYI\n";
                    M_ridge_ratio[i]=ratio_FYI;
                    hi = thick_FYI;
                }
                else
                {
                    if(M_osisaf_type[i]<1.5)//1. no ice in OSISAF
                    {
                        M_ridge_ratio[i]=0.;
                        hi=thick_FYI;
                    }
                    else if(M_osisaf_type[i]<2.5)//2. FYI
                    {
                        M_ridge_ratio[i]=ratio_FYI;
                        hi=thick_FYI;
                    }
                    else if(M_osisaf_type[i]<3.5)//3. MYI
                    {
                        M_ridge_ratio[i]=ratio_MYI;
                        hi  = thick_MYI;
                    }
                    else if(M_osisaf_type[i]<4.5)//4. mixed
                    {
                        M_ridge_ratio[i]=ratio_Mixed;
                        hi      =thick_Mixed;
                    }
                    else//can't happen, can it?
                    {
                        M_ridge_ratio[i]=ratio_Mixed;
                        hi=thick_Mixed;
                    }
                }//not Kara
            }//ice present
            else
            {
                M_ridge_ratio[i]=0.;
            }

            if ( M_conc[i] < 0.01 || hi < physical::hmin )
            {
                //if either c or h equal zero, we set the others to zero as well
                M_conc[i]=0.;
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
                M_ridge_ratio[i]=0.;
            }
            else
            {
                //convert from absolute to effective thickness
                M_ridge_ratio[i]=M_ridge_ratio[i]*M_conc[i];
                M_thick[i]=hi*M_conc[i];
                M_snow_thick[i]=hs*M_conc[i];
            }

            if(bool set_young = vm["thermo.newice_type"].as<int>()==4)
            {
                M_conc_young[i]=std::min(1., std::max(M_amsre_conc[i]-M_conc[i],0.));
                M_h_young[i]=M_conc_young[i]*(h_young_min+0.5*(h_young_max-h_young_min));
            }

            M_damage[i]=0.;
        }
    }//topazIceOsisafIcesat

// -----------------------------------------------------------------------------------------------------------
//! Initializes ice state from Topaz forecasts.
//! Called by the initIce() function.
    void
    Initialiser::topazForecastIce()
    {
        external_data topaz_conc = ExternalData(&M_ocean_elements_dataset, M_mesh, 3, false, time_init);
        // NB TOPAZ gives absolute ice/snow thickness
        // - cell_methods = "mean_where_ice"
        external_data topaz_thick = ExternalData(&M_ocean_elements_dataset, M_mesh, 4, false, time_init);
        external_data topaz_snow_thick = ExternalData(&M_ocean_elements_dataset, M_mesh, 5, false, time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&topaz_conc);
        external_data_tmp.push_back(&topaz_thick);
        external_data_tmp.push_back(&topaz_snow_thick);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - TOPAZ ice forecast ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        double tmp_var;
        for (int i=0; i<M_num_elements; ++i)
        {
            // - TOPAZ puts very small values instead of 0.
            // - uses absolute thickness not effective thickness
            tmp_var = std::min(1., topaz_conc[i]);
            M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.;
            tmp_var = topaz_thick[i];
            M_thick[i] = (tmp_var>1e-14) ? tmp_var*M_conc[i] : 0.;
            tmp_var = topaz_snow_thick[i];
            M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var*M_conc[i] : 0.;

            //if either c or h equal zero, we set the others to zero as well
            if ( M_conc[i] < 0.01 || M_thick[i] < (M_conc[i]*physical::hmin) )
            {
                M_conc[i]        = 0.;
                M_thick[i]       = 0.;
                M_snow_thick[i]  = 0.;
            }

            //init damage and ridge ratio to 0.
            M_ridge_ratio[i] = 0.;
            M_damage[i]=0.;
        }
    }//topazForecastIce

// -----------------------------------------------------------------------------------------------------------
//! Initializes ice state from Topaz-AMSR2.
//! Called by the initIce() function.
    void
    Initialiser::topazForecastAmsr2Ice()
    {
        double init_conc_tmp;

        external_data M_conc_amsr2=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init-0.5);
        external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
        external_data M_init_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
        external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_conc_amsr2);
        external_data_tmp.push_back(&M_init_conc);
        external_data_tmp.push_back(&M_init_thick);
        external_data_tmp.push_back(&M_init_snow_thick);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - TOPAZ ice forecast/AMSR2 ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        double tmp_var;
        for (int i=0; i<M_num_elements; ++i)
        {
            double uncertainty;
            if(M_conc_amsr2[i]<0.1)
                uncertainty=0.1;
            else
                uncertainty=0.05;

            double diff_mod_obs = M_conc_amsr2[i]-M_init_conc[i];
            if(std::abs(diff_mod_obs)>=uncertainty && M_conc_amsr2[i]<=1.)
                // NB missing value for AMSR2 when not over land is 1.15
                // move towards AMSR2 value by the amount uncertainty/2
                M_conc[i] = std::min(1., M_conc_amsr2[i]-(diff_mod_obs/std::abs(diff_mod_obs))*uncertainty/2.);
            else
                M_conc[i] = std::min(1., M_init_conc[i]);

            // TOPAZ puts very small values instead of 0.
            tmp_var=M_init_conc[i];
            init_conc_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
            tmp_var=M_init_thick[i];
            M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
            tmp_var=M_init_snow_thick[i];
            M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

            // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
            if(init_conc_tmp>0.05)
            {
                //convert from absolute to effective thickness
                M_thick[i]*=M_conc[i];
                M_snow_thick[i]*=M_conc[i];
            }

            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_conc[i]=0.;
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }
            if(M_thick[i]<=0.)
            {
                M_thick[i]=0.;
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }
            if(M_thick[i]<0.1*M_conc[i])
                //10cm min ice thickness
                M_thick[i]=0.1*M_conc[i];

            M_damage[i]=0.;
            M_ridge_ratio[i]=0.;
        }
    }//topazForecastAmsr2Ice

// -----------------------------------------------------------------------------------------------------------
//! Get the maximum and minimum ice concentration corresponding to the original
//! NIC ice charts
//! called by <Initialiser::topazForecastAmsr2OsisafNicIce>() and
//! <Initialiser::assimilate_topazForecastAmsr2OsisafNicIce>()
    void
    Initialiser::concBinsNic(double &young_conc_obs_min, double &young_conc_obs_max,
                               double ci, bool use_weekly_nic)
    {

        if(ci<=0.)
        {
            young_conc_obs_min = 0.;
            young_conc_obs_max = 0.;
        }
        else if(!use_weekly_nic)
        {
            if(ci<=0.45) // CT18
            {
                young_conc_obs_min = 0.1;
                young_conc_obs_max = 0.8;
            }
            else if(ci<=0.9) // CT81
            {
                young_conc_obs_min = 0.8;
                young_conc_obs_max = 1.;
            }
        }
        else
        {
            if(ci<=0.2) // CT13
            {
                young_conc_obs_min = 0.1;
                young_conc_obs_max = 0.3;
            }
            else if(ci<=0.30) // CT24
            {
                young_conc_obs_min = 0.2;
                young_conc_obs_max = 0.4;
            }
            else if(ci<=0.50) // CT46
            {
                young_conc_obs_min = 0.4;
                young_conc_obs_max = 0.6;
            }
            else if(ci<=0.70) //CT68
            {
                young_conc_obs_min = 0.6;
                young_conc_obs_max = 0.8;
            }
            else if(ci<=0.90) // CT81
            {
                young_conc_obs_min = 0.8;
                young_conc_obs_max = 1.0;
            }
            else if(ci<=1.) // CT92
            {
                young_conc_obs_min = 0.9;
                young_conc_obs_max = 1.0;
            }
        }
    }//concBinsNic

// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from Topaz forecast, AMSR2 and Osisaf data.
//! Called by the initIce() function.
    void
    Initialiser::topazForecastAmsr2OsisafIce()
    {
        double real_thickness, init_conc_tmp;

        // observations
        external_data M_osisaf_conc=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,0,false,time_init-0.5);
        external_data M_osisaf_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init-0.5);
        external_data M_amsr2_conc=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init-0.5);

        // topaz
        external_data M_topaz_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
        external_data M_topaz_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
        external_data M_topaz_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_osisaf_conc);
        external_data_tmp.push_back(&M_osisaf_type);
        external_data_tmp.push_back(&M_amsr2_conc);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - OSISAF/AMSR2 ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init-0.5, RX, RY);

        external_data_tmp.resize(0);
        external_data_tmp.push_back(&M_topaz_conc);
        external_data_tmp.push_back(&M_topaz_thick);
        external_data_tmp.push_back(&M_topaz_snow_thick);
        LOG(DEBUG)<<"init - TOPAZ ice forecast ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        const double SST_limit = vm["ideal_simul.init_SST_limit"].as<double>();
        const double h_young_min = vm["thermo.h_young_min"].as<double>();
        const double h_young_max = vm["thermo.h_young_max"].as<double>();
        for (int i=0; i<M_num_elements; ++i)
        {
            // TOPAZ puts very small values instead of 0,
            // so set things to zero if below a threshold

            // get absolute ice and snow thicknesses
            double tmp_var=M_topaz_thick[i];
            double hi  = (tmp_var>1e-14) ? tmp_var : 0.;// absolute thickness
            tmp_var=M_topaz_snow_thick[i];
            double hs  = (tmp_var>1e-14) ? tmp_var : 0.;// absolute snow thickness

            tmp_var=std::min(1.,M_topaz_conc[i]);
            M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.;
            if(M_conc[i]>0.)
                // use osisaf only where topaz says there is ice
                // to avoid near land issues and fake concentration over the ocean
                M_conc[i] = M_osisaf_conc[i];

            if(M_amsr2_conc[i]<M_conc[i])
                // AMSR2 is higher resolution and sees
                // small openings that would not be see in OSISAF
                // NB AMSR2 = 1.15 if missing data over ocean, however this will not affect the example here
                M_conc[i] = M_amsr2_conc[i];


            //M_type[i]==1. // No ice
            //M_type[i]==2. // First-Year ice
            //M_type[i]==3. // Multi-Year ice
            //M_type[i]==4. // Mixed
            double ratio_FYI=0.3;
            double ratio_MYI=0.9;
            double ratio_Mixed=0.5*(ratio_FYI+ratio_MYI);

            double thickfac_FYI=1.;
            double thickfac_MYI=1.5;
            double thickfac_Mixed=0.5*(thickfac_FYI+thickfac_MYI);

            if( (hi>0.) && (M_conc[i])>0.2 )
            {

                if(M_mesh_basename.find("kara") != std::string::npos)
                {
                    LOG(DEBUG) <<"Type information is not used for the kara meshes, "
                               <<"we assume there is only FYI\n";
                    M_ridge_ratio[i]=ratio_FYI;
                    hi*=thickfac_FYI;
                }
                else
                {
                    if(M_osisaf_type[i]<1.5)//1. no ice in OSISAF
                    {
                        M_ridge_ratio[i]=0.;
                        hi*=thickfac_FYI;
                    }
                    else if(M_osisaf_type[i]<2.5)//2. FYI
                    {
                        M_ridge_ratio[i]=ratio_FYI;
                        hi*=thickfac_FYI;
                    }
                    else if(M_osisaf_type[i]<3.5)//3. MYI
                    {
                        M_ridge_ratio[i]=ratio_MYI;
                        hi*=thickfac_MYI;
                    }
                    else if(M_osisaf_type[i]<4.5)//4. mixed
                    {
                        M_ridge_ratio[i]=ratio_Mixed;
                        hi*=thickfac_Mixed;
                    }
                    else if(M_osisaf_type[i]>4.)//can't happen, can it?
                    {
                        M_ridge_ratio[i]=ratio_Mixed;
                        hi*=thickfac_Mixed;
                    }
                }//not Kara
            }//ice present
            else
            {
                M_ridge_ratio[i]=0.;
                hi=0.;
            }

            if ( M_conc[i] < 0.01 || hi < physical::hmin )
            {
                //if either c or h equal zero, we set the others to zero as well
                M_conc[i]=0.;
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
                M_ridge_ratio[i]=0.;
                M_damage[i]=0.;
            }
            else
            {
                // convert from absolute to effective thickness
                M_ridge_ratio[i]=M_ridge_ratio[i]*M_conc[i];
                M_snow_thick[i] = M_conc[i]*hs;
                M_thick[i] = M_conc[i]*hi;
                M_damage[i]=0.;
            }

            if(bool set_young = vm["thermo.newice_type"].as<int>()==4)
            {
                M_conc_young[i]=std::max(M_amsr2_conc[i]-M_conc[i],0.);
                M_h_young[i]=M_conc_young[i]*(h_young_min+0.5*(h_young_max-h_young_min));
            }

        }//loop over elements
    }//topazForecastAmsr2OsisafIce

// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from Topaz forecast, AMSR2, Osisaf and NIC ice charts data.
//! Called by the initIce() function.
    void
    Initialiser::topazForecastAmsr2OsisafNicIce(bool use_weekly_nic)
    {
        //observations
        external_data M_osisaf_conc=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,0,false,time_init-0.5);
        external_data M_osisaf_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init-0.5);
        external_data M_amsr2_conc=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init-0.5);

        //topaz
        external_data M_topaz_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
        external_data M_topaz_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
        external_data M_topaz_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_osisaf_conc);
        external_data_tmp.push_back(&M_osisaf_type);
        external_data_tmp.push_back(&M_amsr2_conc);

        external_data M_nic_conc;
        if(use_weekly_nic)
            M_nic_conc = ExternalData(&M_ice_nic_weekly_elements_dataset,
                                      M_mesh, 0, false, time_init-0.5);
        else
            M_nic_conc = ExternalData(&M_ice_nic_elements_dataset,
                                      M_mesh, 0, false, time_init-0.5);
        external_data_tmp.push_back(&M_nic_conc);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - OSISAF/AMSR2/NIC ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init-0.5, RX, RY);

        external_data_tmp.resize(0);
        external_data_tmp.push_back(&M_topaz_conc);
        external_data_tmp.push_back(&M_topaz_thick);
        external_data_tmp.push_back(&M_topaz_snow_thick);
        LOG(DEBUG)<<"init - TOPAZ ice forecast ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        const double SST_limit = vm["ideal_simul.init_SST_limit"].as<double>();
        const double h_young_min = vm["thermo.h_young_min"].as<double>();
        const double h_young_max = vm["thermo.h_young_max"].as<double>();
        for (int i=0; i<M_num_elements; ++i)
        {
            // TOPAZ puts very small values instead of 0,
            // so set things to zero if below a threshold

            // get absolute ice and snow thicknesses
            double tmp_var=M_topaz_thick[i];
            double hi  = (tmp_var>1e-14) ? tmp_var : 0.;// absolute thickness
            tmp_var=M_topaz_snow_thick[i];
            double hs  = (tmp_var>1e-14) ? tmp_var : 0.;// absolute snow thickness

            tmp_var=std::min(1.,M_topaz_conc[i]);
            M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.;
            if(     (M_conc[i]>0.)
                    && (M_amsr2_conc[i]>.15)
                    && (M_amsr2_conc[i]<=1.))
                // use amsr2 only where
                // - topaz says there is ice to avoid near land issues and fake concentration over the ocean
                // - it is large enough to be trusted
                // - it is not masked (mask value = 1.15, but using masking from dataset fills missing values to 0)
                M_conc[i]=M_amsr2_conc[i];

            double ratio_FYI=0.3;
            double ratio_MYI=0.9;
            double ratio_Mixed=0.5*(ratio_FYI+ratio_MYI);

            double thickfac_FYI=1.;
            double thickfac_MYI=1.5;
            double thickfac_Mixed=0.5*(thickfac_FYI+thickfac_MYI);

            if( (hi>0.) && (M_conc[i])>0.2 )
            {

                if(M_mesh_basename.find("kara") != std::string::npos)
                {
                    LOG(DEBUG) <<"Type information is not used for the kara meshes,"
                               << " we assume there is only FYI\n";
                    M_ridge_ratio[i]=ratio_FYI;
                    //M_thick[i]=thick_FYI;
                    hi *= thickfac_FYI;
                }
                else
                {
                    if(M_osisaf_type[i]<1.5)//1. no ice in OSISAF
                    {
                        M_ridge_ratio[i]=0.;
                        hi*=thickfac_FYI;
                    }
                    else if(M_osisaf_type[i]<2.5)//2. FYI
                    {
                        M_ridge_ratio[i]=ratio_FYI;
                        hi*=thickfac_FYI;
                    }
                    else if(M_osisaf_type[i]<3.5)//3. MYI
                    {
                        M_ridge_ratio[i]=ratio_MYI;
                        hi*=thickfac_MYI;
                    }
                    else if(M_osisaf_type[i]<4.5)//4. mixed
                    {
                        M_ridge_ratio[i]=ratio_Mixed;
                        hi*=thickfac_Mixed;
                    }
                    else//can't happen, can it?
                    {
                        M_ridge_ratio[i]=ratio_Mixed;
                        hi*=thickfac_Mixed;
                    }
                }//not Kara
            }//ice present
            else
            {
                M_ridge_ratio[i]=0.;
                hi = 0.;
            }

            //if either c or h equal zero, we set the others to zero as well
            if ( M_conc[i] < 0.01 || hi < physical::hmin )
            {
                M_conc[i]=0.;
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
                M_ridge_ratio[i]=0.;
                hi = 0.;
            }
            else
            {
                // convert from absolute to effective thickness
                M_ridge_ratio[i]=M_ridge_ratio[i]*M_conc[i];
                M_snow_thick[i] = M_conc[i]*hs;
                M_thick[i] = M_conc[i]*hi;
            }

            // skip application of NIC if its conc is >1
            // - then it is masked
            // - unfortunately, applying the mask in datasets led to masked values being treated as real values of 0.0
            //   so we have to do it manually here
            if(M_nic_conc[i]<=1.)
            {
                // Use the NIC ice charts
                // - get conc bins from NIC dataset
                double young_conc_obs = 0.;
                double young_conc_obs_min = 0.;
                double young_conc_obs_max = 0.;
                this->concBinsNic(young_conc_obs_min, young_conc_obs_max, M_nic_conc[i], use_weekly_nic);

                if((M_amsr2_conc[i]>=young_conc_obs_min) && (M_amsr2_conc[i]<=young_conc_obs_max))
                {
                    young_conc_obs_min=M_amsr2_conc[i];
                    young_conc_obs_max=M_amsr2_conc[i];
                }
                else
                {
                    young_conc_obs_min=0.5*(young_conc_obs_min+young_conc_obs_max);
                    young_conc_obs_max=young_conc_obs_min;
                }

                if(bool set_young = vm["thermo.newice_type"].as<int>()==4)
                {
                    M_conc_young[i]=0.;

                    young_conc_obs = young_conc_obs_min-M_conc[i];
                    if(young_conc_obs>=0.)
                    {
                        //if(young_conc_obs>M_conc_young[i])
                        //    M_h_young[i] = M_h_young[i]+(h_young_min + (h_young_max/2.-h_young_min)*0.5)*(young_conc_obs-M_conc_young[i]);
                        //else
                        //    M_h_young[i] = M_h_young[i]*young_conc_obs/M_conc_young[i];

                        M_conc_young[i] = young_conc_obs;
                        M_h_young[i] = (h_young_min + (h_young_max/2.-h_young_min)*0.5)*M_conc_young[i];
                    }
                    else
                    {
                        M_conc_young[i]=0.;
                        M_h_young[i]=0.;

                        M_conc[i]=M_conc[i]+young_conc_obs;
                        M_thick[i]=hi*M_conc[i];
                    }

                }//young ice
                else
                {
                    if(M_conc[i]<young_conc_obs_min)
                    {
                        M_thick[i] = M_thick[i] + std::max(hi,0.5)*(young_conc_obs_min-M_conc[i]); // 50 cm minimum for the added ice
                        M_conc[i] = young_conc_obs_min;
                    }
                    else if(M_conc[i]>young_conc_obs_max)
                    {
                        M_conc[i] = young_conc_obs_max;
                        M_thick[i]=hi*M_conc[i];
                    }
                }//no young ice
            }//use NIC

            M_damage[i]=0.;
        }//loop over elements
    }//topazForecastAmsr2OsisafNicIce

// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from PIOMAS outputs.
//! Called by the initIce() function.
    void
    Initialiser::piomasIce()
    {
        external_data M_init_conc=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,0,false,time_init);
        external_data M_init_thick=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,1,false,time_init);
        external_data M_init_snow_thick=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,2,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_init_conc);
        external_data_tmp.push_back(&M_init_thick);
        external_data_tmp.push_back(&M_init_snow_thick);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - PIOMAS ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        for (int i=0; i<M_num_elements; ++i)
        {
            M_conc[i] = std::min(1.,M_init_conc[i]);
            M_thick[i] = M_init_thick[i];
            M_snow_thick[i] = M_init_snow_thick[i];

            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }
            if(M_thick[i]<=0.)
            {
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }

            M_damage[i]=0.;
        }
    }//piomasIce

// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from NEMO outputs.
//! Called by the initIce() function.
    void
    Initialiser::nemoIce()
    {
        Dataset nemo = DataSet("ice_nemo_elements");
        external_data init_conc=ExternalData(&nemo,M_mesh,0,false,time_init);
        external_data init_thick=ExternalData(&nemo,M_mesh,1,false,time_init);
        external_data init_snow_thick=ExternalData(&nemo,M_mesh,2,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&init_conc);
        external_data_tmp.push_back(&init_thick);
        external_data_tmp.push_back(&init_snow_thick);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - NEMO ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);
        // Surface temperature over which we consider there is no ice when init.
        // There is only ice if sst <= t_freez + sst_limit (tunable)
        double const SST_limit = vm["ideal_simul.init_SST_limit"].as<double>();
        // In nemo code default value is 2, just like in neXtSIM.
        for (int i=0; i<M_num_elements; ++i)
        {
            M_conc[i] = std::min(1.,init_conc[i]);
            M_thick[i] = init_thick[i];
            M_snow_thick[i] = init_snow_thick[i];

            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }
            if(M_thick[i]<=0.)
            {
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }
            // if (M_sst[i] > this->freezingPoint(M_sss[i]) + SST_limit )
            if (M_sst[i] > -physical::mu*M_sss[i] + SST_limit )
            {
                M_thick[i]=0.;
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }
            M_damage[i]=0.;
        }
    }//nemoIce

// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from HYCOM-CICE outputs.
//! Called by the initIce() function.
    void
    Initialiser::ciceIce()
    {
        Dataset cice = DataSet("ice_cice_elements");
        external_data init_conc=ExternalData(&cice,M_mesh,0,false,time_init);
        external_data init_thick=ExternalData(&cice,M_mesh,1,false,time_init);
        external_data init_snow_thick=ExternalData(&cice,M_mesh,2,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&init_conc);
        external_data_tmp.push_back(&init_thick);
        external_data_tmp.push_back(&init_snow_thick);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - CICE ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);
        // Surface temperature over which we consider there is no ice when init.
        // There is only ice if sst <= t_freez + sst_limit (tunable)
        double const SST_limit = vm["ideal_simul.init_SST_limit"].as<double>();
        // In nemo code default value is 2, just like in neXtSIM.
        for (int i=0; i<M_num_elements; ++i)
        {
            M_conc[i] = std::min(1.,init_conc[i]);
            M_thick[i] = init_thick[i];
            M_snow_thick[i] = init_snow_thick[i];

            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }
            if(M_thick[i]<=0.)
            {
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }
            //if (M_sst[i] > this->freezingPoint(M_sss[i]) + SST_limit )
            if (M_sst[i] > -physical::mu*M_sss[i] + SST_limit )
            {
                M_thick[i]=0.;
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }
            M_damage[i]=0.;
        }
    }//ciceIce

// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from Topaz, AMSRE data.
//! Called by the initIce() function.
    void
    Initialiser::topazAmsreIce()
    {
        double real_thickness, init_conc_tmp;

        //obs
        external_data M_conc_amsre=ExternalData(&M_ice_amsre_elements_dataset,M_mesh,0,false,time_init);

        //topaz
        external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
        external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
        external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_conc_amsre);
        external_data_tmp.push_back(&M_init_conc);
        external_data_tmp.push_back(&M_init_thick);
        external_data_tmp.push_back(&M_init_snow_thick);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - TOPAZ/AMSR-E ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        double tmp_var;
        for (int i=0; i<M_num_elements; ++i)
        {
            M_conc[i] = std::min(1.,M_conc_amsre[i]);

            // TOPAZ puts very small values instead of 0.
            tmp_var=M_init_conc[i];
            init_conc_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
            tmp_var=M_init_thick[i];
            M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
            tmp_var=M_init_snow_thick[i];
            M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

            // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
            if(init_conc_tmp>0.05)
            {
                //convert from absolute to effective thickness (with AMSRE conc)
                M_thick[i] *= M_conc[i];
                M_snow_thick[i] *= M_conc[i];
            }

            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_conc[i]=0.;
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }
            if(M_thick[i]<=0.)
            {
                M_thick[i]=0.;
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }

            M_damage[i]=0.;
            M_ridge_ratio[i]=0.;
        }
    }//topazAmsreIce TODO no young ice; logic needs checking; no ice-type option for this

// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from Topaz, AMSR2 data.
//! Called by the initIce() function.
    void
    Initialiser::topazAmsr2Ice()
    {
        double real_thickness, init_conc_tmp;

        external_data M_conc_amsr2=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
        external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
        external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
        external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_conc_amsr2);
        external_data_tmp.push_back(&M_init_conc);
        external_data_tmp.push_back(&M_init_thick);
        external_data_tmp.push_back(&M_init_snow_thick);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - TOPAZ/AMSR2 ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        double tmp_var;
        for (int i=0; i<M_num_elements; ++i)
        {
            double uncertainty;
            if(M_conc_amsr2[i]<0.1)
                uncertainty=0.1;
            else
                uncertainty=0.05;

            double diff_mod_obs = M_conc_amsr2[i]-M_init_conc[i];
            if(std::abs(diff_mod_obs)>=uncertainty && M_conc_amsr2[i]<=1.)
                // NB missing value for AMSR2 when not over land is 1.15
                // move towards AMSR2 value by the amount uncertainty/2
                M_conc[i] = std::min(1., M_conc_amsr2[i]-(diff_mod_obs/std::abs(diff_mod_obs))*uncertainty/2.);
            else
                M_conc[i] = std::min(1., M_init_conc[i]);

            // TOPAZ puts very small values instead of 0.
            tmp_var=M_init_conc[i];
            init_conc_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
            tmp_var=M_init_thick[i];
            M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
            tmp_var=M_init_snow_thick[i];
            M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

            // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
            if(init_conc_tmp>0.05)
            {
                //convert from absolute to effective thickness
                M_thick[i]*=M_conc[i];
                M_snow_thick[i]*=M_conc[i];
            }

            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_conc[i]=0.;
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }
            if(M_thick[i]<=0.)
            {
                M_thick[i]=0.;
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }

            M_damage[i]=0.;
            M_ridge_ratio[i]=0.;
        }
    }//topazAmsr2Ice TODO no young ice; logic needs checking; no ice-type option for this

//! Initializes the ice state from  AMSR2 data with a constant sea ice thickness
//! Called by the initIce() function.
    void
    Initialiser::amsr2ConstThickIce()
    {
        double real_thickness, init_conc_tmp;

        double h_const = vm["ideal_simul.init_thickness"].as<double>();
        double hs_const = vm["ideal_simul.init_snow_thickness"].as<double>();

        external_data M_conc_amsr2=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
        external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_conc_amsr2);
        external_data_tmp.push_back(&M_init_conc);
        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - TOPAZ/AMSR2 ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        double tmp_var;
        for (int i=0; i<M_num_elements; ++i)
        {
            double uncertainty;
            if(M_conc_amsr2[i]<0.1)
                uncertainty=0.1;
            else
                uncertainty=0.05;
            // TEST : take AMSR2 Raw value
            M_conc[i] = std::min(1., M_conc_amsr2[i] ) ;
//        double diff_mod_obs = M_conc_amsr2[i]-M_init_conc[i];
//        if(std::abs(diff_mod_obs)>=uncertainty && M_conc_amsr2[i]<=1.)
//            // NB missing value for AMSR2 when not over land is 1.15
//            // move towards AMSR2 value by the amount uncertainty/2
//            M_conc[i] = std::min(1., M_conc_amsr2[i]-(diff_mod_obs/std::abs(diff_mod_obs))*uncertainty/2.);
//        else
//            M_conc[i] = std::min(1., M_init_conc[i]);

            // TOPAZ puts very small values instead of 0.
            tmp_var=M_init_conc[i];
            init_conc_tmp = (tmp_var>1e-14) ? tmp_var : 0.;

            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_conc[i]=0.;
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }
            else
            {
                M_thick[i]=h_const*M_conc[i];
                M_snow_thick[i]=hs_const*M_conc[i];
            }
            M_damage[i]=0.;
            M_ridge_ratio[i]=0.;
        }
    }//topaz] = std::min(1., M_conc_amsr2[i]Amsr2Ice TODO no young ice; logic needs checking; no ice-type option for this

// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from CS2 SMOS data.
//! Called by the initIce() function.
    void
    Initialiser::cs2SmosIce()
    {
        external_data M_init_conc=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,0,false,time_init);
        external_data M_init_thick=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,1,false,time_init);
        external_data M_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_init_conc);
        external_data_tmp.push_back(&M_init_thick);
        external_data_tmp.push_back(&M_type);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - CS2/SMOS ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        this->warrenClimatology();

        double tmp_var, correction_factor_warren;
        // Observational operator for correcting young ice concentrations
        // Coefficients are fitted using data (figure) from Thomas Lavergne
        // Atot = Aosisaf / OO(H)
        // OO(H) = 2y(exp(H/s + 0) / (1 + exp(H/s + 0)) - y
        // OO is applied only if H < 0.5

        // These coefs are fitting exactly data from Thomas. 5% gap for thickness above 0.5 m.
        // double oo_y = 0.95691881, oo_s = 0.06787237, oo_o = 0.42548199;

        // These coeffs are fitted to slightly modified data from Thomas - last 3 points are set to 1.
        // No gap  for thickness above 0.5 m.
        double oo_y = 0.99853900, oo_s = 0.09128844, oo_o = 0.65327181;

        const double SST_limit = vm["ideal_simul.init_SST_limit"].as<double>();
        const double h_young_min = vm["thermo.h_young_min"].as<double>();
        const double h_young_max = vm["thermo.h_young_max"].as<double>();
        const double compaction_param = vm["dynamics.compaction_param"].as<double>();

        // fraction of young ice
        double fy0 = 0.10;
        // mean allowed thickness of young ice
        double ty0 = h_young_min + 0.5 * (h_young_max - h_young_min);
        //double ty0 = 0.05;

        for (int i=0; i<M_num_elements; ++i)
        {
            M_thick[i] = M_init_thick[i];
            M_conc[i] = std::min(1., M_init_conc[i]);

            // correction of concentration by thickness using OO
            tmp_var = std::exp(M_thick[i] / oo_s + oo_o);
            tmp_var = 2 * oo_y * tmp_var / (1 + tmp_var) - oo_y;
            M_conc[i] /= tmp_var;
            M_conc[i] = std::min(1., M_conc[i]);

            double ratio_FYI=0.3;
            double ratio_MYI=0.9;
            double ratio_Mixed=0.6;

            if((M_thick[i]>0.)&&(M_conc[i])>0.2)
            {
                if(M_type[i]<=1.)
                    M_ridge_ratio[i]=0.;
                if(M_type[i]>1. && M_type[i]<=2.)
                    M_ridge_ratio[i]=(M_type[i]-1.)*ratio_FYI;
                if(M_type[i]>2. && M_type[i]<=3.)
                    M_ridge_ratio[i]=(1.-(M_type[i]-2.))*ratio_FYI + (M_type[i]-2.)*ratio_MYI;
                if(M_type[i]>3. && M_type[i]<=4.)
                    M_ridge_ratio[i]=(1.-(M_type[i]-3.))*ratio_MYI + (M_type[i]-3.)*ratio_Mixed;
                if(M_type[i]>4.)
                    M_ridge_ratio[i]=ratio_Mixed;

                M_ridge_ratio[i]=M_ridge_ratio[i]*std::exp(compaction_param*(1.-M_conc[i]));
            }
            else
                M_ridge_ratio[i]=0.;


            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }
            if(M_thick[i]<=0.)
            {
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }

            // Correction of the value given by Warren as a function of the ice type
            //M_type[i]==1. // No ice
            //M_type[i]==2. // First-Year ice
            //M_type[i]==3. // Multi-Year ice
            //M_type[i]==4. // Mixed
            correction_factor_warren=std::max(0.,std::min(1.,(M_type[i]-1.)*0.5)); // == 1. for MY, and mixed, 0.5 for FY, 0. for No ice

            M_snow_thick[i]=correction_factor_warren*M_snow_thick[i]*M_conc[i];

            M_damage[i]=0.;

            if(bool set_young = vm["thermo.newice_type"].as<int>()==4)
            {
                // In case of young ice
                // Concentration is split into fy0 and 1 - fy0 for young and old ice
                M_conc_young[i] = fy0 * M_conc[i];
                M_conc[i]  *= (1. - fy0);
                if (M_thick[i] / M_conc[i] <= ty0) {
                    // if absolute thickness (To = Ho / Ao) is below ty0 (see ty0 above):
                    // thickness is also split into fy0 and 1 - fy0
                    M_h_young[i] = fy0 * M_thick[i];
                    M_thick[i] *= (1. - fy0);
                } else {
                    // if absolute thickness (To = Ho / Ao) is above ty0:
                    // young ice thickness is maximum: h = ty0 * a
                    // old ice thickness is the remaining part: H = H0 - h
                    M_h_young[i] = ty0 * M_conc_young[i];
                    M_thick[i] -= M_h_young[i];
                }
            }

            // Check that the snow is not so thick that the ice is flooded
            double max_snow = M_thick[i]*(physical::rhow-physical::rhoi)/physical::rhos;
            M_snow_thick[i] = std::min(max_snow, M_snow_thick[i]);

        }
    }//cs2SmosIce

// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from CS2, SMOS, AMSR2 data.
//! Called by the initIce() function.
    void
    Initialiser::cs2SmosAmsr2Ice()
    {
        external_data M_init_conc=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,0,false,time_init);
        external_data M_amsr2_conc=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
        external_data M_init_thick=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,1,false,time_init);
        external_data M_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_init_conc);
        external_data_tmp.push_back(&M_init_thick);
        external_data_tmp.push_back(&M_type);
        external_data_tmp.push_back(&M_amsr2_conc);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - CS2/SMOS/AMSR2 ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        this->warrenClimatology();

        const double SST_limit = vm["ideal_simul.init_SST_limit"].as<double>();
        const double h_young_min = vm["thermo.h_young_min"].as<double>();
        const double h_young_max = vm["thermo.h_young_max"].as<double>();
        const double compaction_param = vm["dynamics.compaction_param"].as<double>();
        double tmp_var, correction_factor_warren;
        for (int i=0; i<M_num_elements; ++i)
        {
            tmp_var=std::min(1.,M_init_conc[i]);
            M_conc[i] = tmp_var;
            tmp_var=M_init_thick[i];
            M_thick[i] = tmp_var ;
            if(M_amsr2_conc[i]<M_conc[i])
                // AMSR2 is higher resolution and see small opening that would not be see in OSISAF
                // NB AMSR2 = 1.15 if missing data over ocean, however this will not affect the example here
                M_conc[i]=M_amsr2_conc[i];

            double ratio_FYI=0.3;
            double ratio_MYI=0.9;
            double ratio_Mixed=0.6;

            if((M_thick[i]>0.)&&(M_conc[i])>0.2)
            {
                if(M_type[i]<=1.)
                    M_ridge_ratio[i]=0.;
                if(M_type[i]>1. && M_type[i]<=2.)
                    M_ridge_ratio[i]=(M_type[i]-1.)*ratio_FYI;
                if(M_type[i]>2. && M_type[i]<=3.)
                    M_ridge_ratio[i]=(1.-(M_type[i]-2.))*ratio_FYI + (M_type[i]-2.)*ratio_MYI;
                if(M_type[i]>3. && M_type[i]<=4.)
                    M_ridge_ratio[i]=(1.-(M_type[i]-3.))*ratio_MYI + (M_type[i]-3.)*ratio_Mixed;
                if(M_type[i]>4.)
                    M_ridge_ratio[i]=ratio_Mixed;

                M_ridge_ratio[i]=M_ridge_ratio[i]*std::exp(compaction_param*(1.-M_conc[i]));
            }
            else
                M_ridge_ratio[i]=0.;


            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }
            if(M_thick[i]<=0.)
            {
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }

            // Correction of the value given by Warren as a function of the ice type
            //M_type[i]==1. // No ice
            //M_type[i]==2. // First-Year ice
            //M_type[i]==3. // Multi-Year ice
            //M_type[i]==4. // Mixed
            correction_factor_warren=std::max(0.,std::min(1.,(M_type[i]-1.)*0.5)); // == 1. for MY, and mixed, 0.5 for FY, 0. for No ice

            M_snow_thick[i]=correction_factor_warren*M_snow_thick[i]*M_conc[i];

            M_damage[i]=0.;

            if(bool set_young = vm["thermo.newice_type"].as<int>()==4)
            {
                // we assume the young ice is not seen by passive microwave
                // so arbitrarily increase the total conc by 20%, putting this into young ice
                // at a constant absolute thickness (default=27.5cm)
                M_conc_young[i]=std::min(1.-M_conc[i], 0.2*M_conc[i]);
                M_h_young[i]=M_conc_young[i]*(h_young_min+0.5*(h_young_max-h_young_min));
            }

            // Check that the snow is not so thick that the ice is flooded
            double max_snow = M_thick[i]*(physical::rhow-physical::rhoi)/physical::rhos;
            M_snow_thick[i] = std::min(max_snow, M_snow_thick[i]);

        }
    }//cs2SmosAmsr2Ice

// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from SMOS data.
//! Called by the initIce() function.
    void
    Initialiser::smosIce()
    {
        external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
        external_data M_init_thick=ExternalData(&M_ice_smos_elements_dataset,M_mesh,0,false,time_init);

        boost::gregorian::date dt = Nextsim::parse_date(time_init);
        int month_id=dt.month().as_number(); // 1 for January, 2 for February, and so on. This will be used to compute the snow from Warren climatology

        LOG(DEBUG) << "month_id: " << month_id <<"\n";

        external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);

        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_init_conc);
        external_data_tmp.push_back(&M_init_thick);
        external_data_tmp.push_back(&M_init_snow_thick);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - SMOS ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        double tmp_var;
        for (int i=0; i<M_num_elements; ++i)
        {
            tmp_var=std::min(1.,M_init_conc[i]);
            M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
            tmp_var=M_init_thick[i];
            M_thick[i] = tmp_var ;
            tmp_var=M_init_snow_thick[i];
            M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.

            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }
            if(M_thick[i]<=0.)
            {
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }

            M_damage[i]=0.;
            M_ridge_ratio[i]=0.;
        }
    }//smosIce

    void
    Initialiser::glorys12Ice()
    {
        // Initialize the sea ice conc. and thick. from GLORYS12 reanalysis or forecast.
        // Snow thickness is 0 by default (used in southern ocean, with no equivalent of Warren Climatology)
        Dataset ice_glorys12_elements_dataset = DataSet("glorys12_elements");
        external_data M_init_conc   = ExternalData(&ice_glorys12_elements_dataset,M_mesh,3,false,time_init);
        external_data M_init_thick  = ExternalData(&ice_glorys12_elements_dataset,M_mesh,4,false,time_init);

        boost::gregorian::date dt = Nextsim::parse_date(time_init);


        external_data_vec external_data_tmp;
        external_data_tmp.push_back(&M_init_conc);
        external_data_tmp.push_back(&M_init_thick);

        auto RX = M_mesh.bCoordX();
        auto RY = M_mesh.bCoordY();
        LOG(DEBUG)<<"init - Glorys Ice ExternalData objects\n";
        this->checkReloadDatasets(external_data_tmp, time_init, RX, RY);

        double tmp_var;
        for (int i=0; i<M_num_elements; ++i)
        {
            tmp_var=std::min(1.,M_init_conc[i]);
            M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.;
            tmp_var=M_init_thick[i];
            M_thick[i] = tmp_var ;
            tmp_var=0.;//This is a test for the southern ocean
            M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

            //if either c or h equal zero, we set the others to zero as well
            if(M_conc[i]<=0.)
            {
                M_thick[i]=0.;
                M_snow_thick[i]=0.;
            }
            if(M_thick[i]<=0.)
            {
                M_conc[i]=0.;
                M_snow_thick[i]=0.;
            }

            M_damage[i]=0.;
            M_ridge_ratio[i]=0.;
        }
    }//glorysIce

// -----------------------------------------------------------------------------------------------------------
//! Initializes the ice state from CS2, SMOS, AMSR2 data.
//! Called by the cs2SmosIce() and cs2SmosAmsr2Ice() functions.
    void
    Initialiser::warrenClimatology()
    {
        // Coefficients for the fit from Warren et al '99
        std::vector<std::vector<double>> coeffs(12);

        coeffs[0].push_back(28.01);
        coeffs[0].push_back( 0.1270);
        coeffs[0].push_back(-1.1833);
        coeffs[0].push_back(-0.1164);
        coeffs[0].push_back(-0.0051);
        coeffs[0].push_back( 0.0243);

        coeffs[1].push_back(30.28);
        coeffs[1].push_back( 0.1056);
        coeffs[1].push_back(-0.5908);
        coeffs[1].push_back(-0.0263);
        coeffs[1].push_back(-0.0049);
        coeffs[1].push_back( 0.0044);

        coeffs[2].push_back(33.89);
        coeffs[2].push_back( 0.5486);
        coeffs[2].push_back(-0.1996);
        coeffs[2].push_back( 0.0280);
        coeffs[2].push_back( 0.0216);
        coeffs[2].push_back(-0.0176);

        coeffs[3].push_back(36.80);
        coeffs[3].push_back( 0.4046);
        coeffs[3].push_back(-0.4005);
        coeffs[3].push_back( 0.0256);
        coeffs[3].push_back( 0.0024);
        coeffs[3].push_back(-0.0641);

        coeffs[4].push_back(36.93);
        coeffs[4].push_back( 0.0214);
        coeffs[4].push_back(-1.1795);
        coeffs[4].push_back(-0.1076);
        coeffs[4].push_back(-0.0244);
        coeffs[4].push_back(-0.0142);

        coeffs[5].push_back(36.59);
        coeffs[5].push_back( 0.7021);
        coeffs[5].push_back(-1.4819);
        coeffs[5].push_back(-0.1195);
        coeffs[5].push_back(-0.0009);
        coeffs[5].push_back(-0.0603);

        coeffs[6].push_back(11.02);
        coeffs[6].push_back( 0.3008);
        coeffs[6].push_back(-1.2591);
        coeffs[6].push_back(-0.0811);
        coeffs[6].push_back(-0.0043);
        coeffs[6].push_back(-0.0959);

        coeffs[7].push_back( 4.64);
        coeffs[7].push_back( 0.3100);
        coeffs[7].push_back(-0.6350);
        coeffs[7].push_back(-0.0655);
        coeffs[7].push_back( 0.0059);
        coeffs[7].push_back(-0.0005);

        coeffs[8].push_back(15.81);
        coeffs[8].push_back( 0.2119);
        coeffs[8].push_back(-1.0292);
        coeffs[8].push_back(-0.0868);
        coeffs[8].push_back(-0.0177);
        coeffs[8].push_back(-0.0723);

        coeffs[9].push_back(22.66);
        coeffs[9].push_back( 0.3594);
        coeffs[9].push_back(-1.3483);
        coeffs[9].push_back(-0.1063);
        coeffs[9].push_back( 0.0051);
        coeffs[9].push_back(-0.0577);

        coeffs[10].push_back(25.57);
        coeffs[10].push_back( 0.1496);
        coeffs[10].push_back(-1.4643);
        coeffs[10].push_back(-0.1409);
        coeffs[10].push_back(-0.0079);
        coeffs[10].push_back(-0.0258);

        coeffs[11].push_back(26.67);
        coeffs[11].push_back(-0.1876);
        coeffs[11].push_back(-1.4229);
        coeffs[11].push_back(-0.1413);
        coeffs[11].push_back(-0.0316);
        coeffs[11].push_back(-0.0029);

        // Time interpolation
        boost::gregorian::date t = Nextsim::parse_date(time_init);
        int month  = t.month().as_number(); // 1 for January, 2 for February, and so on.
        int day    = t.day().as_number(); // Day of the month.
        int eomday = t.end_of_month().day().as_number(); // Last day of the month.

        int month2; // The other month to interpolate from (one before or after)
        double dt;
        if ( day < eomday/2. )
        { // We're in the early part of the month and interpolate to the previous month
            month2 = month-1;
            dt     = eomday/2+day;
            if(month2==0)
                month2=12;
        }
        else
        { // We're in the late part of the month and interpolate to the next month
            month2 = month+1;
            dt     = eomday/2 + eomday - day;
            if(month2==13)
                month2=1;
        }

        // Now calculate snow thickness for the current day as an inexact temporal interpolation
        // I just assume all months are as long as the current one ... but the error incurred is small
        std::vector<double> lon = M_mesh.meanLon();
        std::vector<double> lat = M_mesh.meanLat();

        for ( int i=0; i<M_num_elements; i++ )
        {
            const double pi = std::acos(-1);

            double x = (90 - lat[i]) * std::cos(lon[i]*pi/180.);
            double y = (90 - lat[i]) * std::sin(lon[i]*pi/180.);

            M_snow_thick[i] = 1e-2*dt/eomday*std::max(0.,
                                                      coeffs[month-1][0] + coeffs[month-1][1]*x + coeffs[month-1][2]*y + coeffs[month-1][3]*x*y + coeffs[month-1][4]*x*x + coeffs[month-1][5]*y*y)
                              + 1e-2*(eomday-dt)/eomday*std::max(0.,
                                                                 coeffs[month2-1][0] + coeffs[month2-1][1]*x + coeffs[month2-1][2]*y + coeffs[month2-1][3]*x*y + coeffs[month2-1][4]*x*x + coeffs[month2-1][5]*y*y);
        }

    }//warrenClimatology
}
