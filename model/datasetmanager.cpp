/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   datasetmanager.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
* @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Mon Aug 24 11:02:45 2015
 */

#include <datasetmanager.hpp>

namespace Nextsim
{

void
DataSetManager::initOptAndParam()
{
    const boost::unordered_map<const std::string, setup::BathymetryType> str2bathymetry = boost::assign::map_list_of
        ("constant", setup::BathymetryType::CONSTANT)
        ("etopo", setup::BathymetryType::ETOPO);
        //! \param M_bathymetry_type (string) Option on the type of bathymetry (constant or ETOPO)
    M_bathymetry_type = OptionHandler::getOptionFromMap(vm, "setup.bathymetry-type", str2bathymetry);
    LOG(DEBUG) <<"BathymetryType= "<< (int) M_bathymetry_type <<"\n";

    //! Sets various time steps (init, thermo, output, mooring, restart) and options on data assimilation and restarts
    if (vm["simul.time_init"].as<std::string>() == "")
        throw std::runtime_error("Please provide simul.time_init option (start time)\n");
    else
        time_init = Nextsim::stringToDatenum(vm["simul.time_init"].as<std::string>()); //! \param time_init (string) Time at which the simulation is started

    M_spinup_duration = vm["simul.spinup_duration"].as<double>(); //! \param M_spinup_duration (double) duration of spinup of atmosphere/ocean forcing.
    LOG(DEBUG) << "M_spinup_duration " << M_spinup_duration << "\n";

    M_ensemble_member = vm["statevector.ensemble_member"].as<int>();

    M_use_nesting= vm["nesting.use_nesting"].as<bool>(); //! \param M_use_nesting (boolean) Option on the use of nested model meshes
    if (M_use_nesting)
    {
        M_use_ocean_nesting = vm["nesting.use_ocean_nesting"].as<bool>();
        M_nest_outer_mesh   = vm["nesting.outer_mesh"].as<std::string>();
        M_nest_inner_mesh   = vm["nesting.inner_mesh"].as<std::string>();
        M_nest_method       = vm["nesting.method"].as<std::string>();
        M_nudge_function    = vm["nesting.nudge_function"].as<std::string>();
        M_nudge_timescale   = vm["nesting.nudge_timescale"].as<double>();
        M_nudge_lengthscale = vm["nesting.nudge_lengthscale"].as<double>();
        M_nest_dynamic_vars = vm["nesting.nest_dynamic_vars"].as<bool>();
    }

    //! Sets options on the atmospheric and ocean forcing, initialization of ice, type of dynamics, bathymetry and on the use of nested meshes
    const boost::unordered_map<const std::string, setup::AtmosphereType> str2atmosphere = boost::assign::map_list_of
        ("constant", setup::AtmosphereType::CONSTANT)
        ("generic_ps", setup::AtmosphereType::GENERIC_PS)
        ("asr", setup::AtmosphereType::ASR)
        ("erai", setup::AtmosphereType::ERAi)
        ("era5", setup::AtmosphereType::ERA5)
        ("ec2", setup::AtmosphereType::EC2)
        ("ec_erai", setup::AtmosphereType::EC_ERAi)
        ("cfsr", setup::AtmosphereType::CFSR)
        ("cfsr_hi", setup::AtmosphereType::CFSR_HI)
        ("ec2_arome", setup::AtmosphereType::EC2_AROME)
        ("ec2_arome_ensemble", setup::AtmosphereType::EC2_AROME_ENSEMBLE);
    M_atmosphere_type = OptionHandler::getOptionFromMap(vm, "setup.atmosphere-type", str2atmosphere);
        //! \param M_atmosphere_type (enum) Option on the type of atm. forcing (constant, forecast or reanalyses)
    LOG(DEBUG)<<"AtmosphereType= "<< (int)M_atmosphere_type <<"\n";

    const boost::unordered_map<const std::string, setup::OceanType> str2ocean = boost::assign::map_list_of
        ("constant", setup::OceanType::CONSTANT)
        ("topaz", setup::OceanType::TOPAZR)
        ("topaz_atrest", setup::OceanType::TOPAZR_atrest)
        ("topaz_forecast", setup::OceanType::TOPAZF)
        ("topaz_altimeter", setup::OceanType::TOPAZR_ALTIMETER)
        ("coupled", setup::OceanType::COUPLED)
        ("glorys12", setup::OceanType::GLORYS12R);
    M_ocean_type = OptionHandler::getOptionFromMap(vm, "setup.ocean-type", str2ocean);
        //! \param M_ocean_type (enum) Option on the type of ocean forcing (constant or Topaz options)
    LOG(DEBUG) <<"OCEANTYPE= "<< (int)M_ocean_type <<"\n";
}

//------------------------------------------------------------------------------------------------------
//! Initializes the external data used by the model
//! * Datasets objects
//! * ExternalData objects
//! Called by the init() function.
void
DataSetManager::initExternalData()
{
    //! - 1) Init the datasets
    this->initDatasets();

    //! - 2) populates the forcing variables.
    LOG(DEBUG) <<"Initialize forcingAtmosphere\n";
    this->forcingAtmosphere();

    LOG(DEBUG) <<"Initialize forcingOcean\n";
    this->forcingOcean();

#ifdef OASIS
    LOG(DEBUG) <<"Initialize forcingWaves\n";
    if (M_couple_waves)
        this->forcingWaves();
#endif

    //! - 3) Initializes the bathymetry using the initBathymetry() function,
    LOG(DEBUG) <<"Initialize bathymetry\n";
    this->initBathymetry();

    //! - 4) Check the external data objects
    //       TODO add the nodes
    if(M_external_data_elements.size() != M_external_data_elements_names.size())
        throw std::runtime_error(
                "M_external_data_elements and M_external_data_elements_names should be the same size");
    for(int i=0; i<M_external_data_elements.size(); i++)
    {
        // check all the forcings on the elements are initialised
        std::string const msg = "ExternalData object "
            + M_external_data_elements_names[i] + " is not initialized";
        if(!M_external_data_elements[i]->isInitialized())
            throw std::logic_error(msg);
    }

    //! - 5) Initialise coupler interface (this can be considered "external data" too)
#ifdef OASIS
    LOG(DEBUG) <<"Initialize OASIS coupler\n";
    this->initOASIS();
#endif

}//initExternalData

//------------------------------------------------------------------------------------------------------
//! Set the 'reload flag' for datasets after remeshing
//! Called by the regrid() and initVariables() functions.
void
DataSetManager::setReloadDatasets()
{
    // The coupled system needs special handling
    if ( M_ocean_type != setup::OceanType::COUPLED )
    {
        M_ocean_nodes_dataset.loaded=false;
        M_ocean_elements_dataset.loaded=false;
        M_ocean_nodes_dataset.grid.loaded=false;
        M_ocean_elements_dataset.grid.loaded=false;
    }
#if defined OASIS
    if (M_couple_waves)
    {
        M_wave_elements_dataset.loaded=false;
        M_wave_elements_dataset.grid.loaded=false;
        if(M_recv_wave_stress)
        {
            M_wave_nodes_dataset.loaded=false;
            M_wave_nodes_dataset.grid.loaded=false;
        }
    }
#endif

    // reload the dataset
    M_atmosphere_nodes_dataset.loaded=false;
    M_atmosphere_elements_dataset.loaded=false;
    M_atmosphere_bis_elements_dataset.loaded=false;
    M_bathymetry_elements_dataset.loaded=false;


    // reload the grid
    M_atmosphere_nodes_dataset.grid.loaded=false;
    M_atmosphere_elements_dataset.grid.loaded=false;
    M_atmosphere_bis_elements_dataset.grid.loaded=false;
    M_bathymetry_elements_dataset.grid.loaded=false;


    // --------------------------------------------------------------
    // interpolation of the dataset
    M_atmosphere_nodes_dataset.interpolated=false;
    M_atmosphere_elements_dataset.interpolated=false;
    M_atmosphere_bis_elements_dataset.interpolated=false;
    M_ocean_nodes_dataset.interpolated=false;
    M_ocean_elements_dataset.interpolated=false;
#ifdef OASIS
    if(M_couple_waves)
    {
        if(M_recv_wave_stress)
            M_wave_nodes_dataset.interpolated = false;
        M_wave_elements_dataset.interpolated=false;
    }
#endif
    M_ice_topaz_elements_dataset.interpolated=false;
    M_ice_piomas_elements_dataset.interpolated=false;
    M_ice_amsre_elements_dataset.interpolated=false;
    M_ice_osisaf_elements_dataset.interpolated=false;
    M_ice_osisaf_type_elements_dataset.interpolated=false;
    M_ice_amsr2_elements_dataset.interpolated=false;
    M_ice_cs2_smos_elements_dataset.interpolated=false;
    M_ice_smos_elements_dataset.interpolated=false;
    M_bathymetry_elements_dataset.interpolated=false;

}//setReloadDatasets


//------------------------------------------------------------------------------------------------------
//! Initializes the Datasets used by the model
//! * Datasets objects
//! * ExternalData objects
//! Called by the initExternalData() function.
void
DataSetManager::initDatasets()
{
    //! - 1) Initializes the atmospheric forcing dataset
    switch(M_atmosphere_type){
        case setup::AtmosphereType::CONSTANT:
            break;

        case setup::AtmosphereType::GENERIC_PS:
            M_atmosphere_nodes_dataset=DataSet("generic_ps_atm_nodes");
            M_atmosphere_elements_dataset=DataSet("generic_ps_atm_elements");
            break;

        case setup::AtmosphereType::ASR:
            M_atmosphere_nodes_dataset=DataSet("asr_nodes");
            M_atmosphere_elements_dataset=DataSet("asr_elements");
            break;

        case setup::AtmosphereType::ERAi:
            M_atmosphere_nodes_dataset=DataSet("ERAi_nodes");
            M_atmosphere_elements_dataset=DataSet("ERAi_elements");
            break;

        case setup::AtmosphereType::ERA5:
            M_atmosphere_nodes_dataset=DataSet("ERA5_nodes");
            M_atmosphere_elements_dataset=DataSet("ERA5_elements");
            break;

        case setup::AtmosphereType::EC2:
            M_atmosphere_nodes_dataset=DataSet("ec2_nodes");
            M_atmosphere_elements_dataset=DataSet("ec2_elements");
            break;

        case setup::AtmosphereType::EC_ERAi:
            M_atmosphere_nodes_dataset=DataSet("ec_nodes");
            M_atmosphere_elements_dataset=DataSet("ec_elements");
            M_atmosphere_bis_elements_dataset=DataSet("ERAi_elements");
            break;

        case setup::AtmosphereType::CFSR:
            M_atmosphere_nodes_dataset=DataSet("cfsr_nodes");
            M_atmosphere_elements_dataset=DataSet("cfsr_elements");
            break;

        case setup::AtmosphereType::CFSR_HI:
            M_atmosphere_nodes_dataset=DataSet("cfsr_nodes_hi");
            M_atmosphere_elements_dataset=DataSet("cfsr_elements");
            break;

        case setup::AtmosphereType::EC2_AROME:
            M_atmosphere_nodes_dataset=DataSet("ec2_arome_nodes");
            M_atmosphere_elements_dataset=DataSet("ec2_arome_elements");
            break;

        case setup::AtmosphereType::EC2_AROME_ENSEMBLE:
            M_atmosphere_nodes_dataset=DataSet("ec2_arome_ensemble_nodes");
            M_atmosphere_elements_dataset=DataSet("ec2_arome_ensemble_elements");
            break;

        default:
            std::cout << "invalid atmospheric forcing"<<"\n";throw std::logic_error("invalid atmospheric forcing");
    }

    //! - 2) Initializes the oceanic forcing dataset
    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
        break;

        case setup::OceanType::TOPAZR:
        case setup::OceanType::TOPAZR_atrest:
            M_ocean_nodes_dataset=DataSet("topaz_nodes");
            M_ocean_elements_dataset=DataSet("topaz_elements");
            break;

        case setup::OceanType::TOPAZR_ALTIMETER:
            M_ocean_nodes_dataset=DataSet("ocean_currents_nodes");
            M_ocean_elements_dataset=DataSet("topaz_elements");
            break;

        case setup::OceanType::TOPAZF:
            M_ocean_nodes_dataset=DataSet("topaz_forecast_nodes");
            M_ocean_elements_dataset=DataSet("topaz_forecast_elements");
            break;
#ifdef OASIS
        case setup::OceanType::COUPLED:
            M_ocean_nodes_dataset=DataSet("ocean_cpl_nodes");
            M_ocean_elements_dataset=DataSet("ocean_cpl_elements");
            break;
#endif
        case setup::OceanType::GLORYS12R:
            M_ocean_nodes_dataset=DataSet("glorys12_nodes");
            M_ocean_elements_dataset=DataSet("glorys12_elements");
            break;
        default:
            std::cout << "invalid ocean forcing"<<"\n";throw std::logic_error("invalid ocean forcing");
    }

    //! - 3) Initializes the nesting datasets if needed
    if (M_use_nesting)
    {
        M_nesting_nodes_dataset=DataSet("nesting_nodes");
        M_nesting_ocean_elements_dataset=DataSet("nesting_ocean_elements");
        M_nesting_ice_elements_dataset=DataSet("nesting_ice_elements");
        M_nesting_dynamics_elements_dataset=DataSet("nesting_dynamics_elements");
        M_nesting_distance_elements_dataset=DataSet("nesting_distance_elements");
        M_nesting_distance_nodes_dataset=DataSet("nesting_distance_nodes");
    }

#ifdef OASIS
    if (M_couple_waves)
    {
        if(M_recv_wave_stress)
        {
            M_wave_nodes_dataset = DataSet("wave_cpl_nodes");
            M_datasets_regrid.push_back(&M_wave_nodes_dataset);
        }
        M_wave_elements_dataset = DataSet("wave_cpl_elements");
        M_datasets_regrid.push_back(&M_wave_elements_dataset);
    }
#endif

    //! - 4) Initializes the ice-init datasets
    //       TODO these probably don't need to be global variables
    //            - in fact they are probably taking up a significant
    //              amount of memory
    M_ice_topaz_elements_dataset=DataSet("ice_topaz_elements");
    M_ice_icesat_elements_dataset=DataSet("ice_icesat_elements");
    M_ice_piomas_elements_dataset=DataSet("ice_piomas_elements");
    M_ice_amsre_elements_dataset=DataSet("ice_amsre_elements");
    M_ice_osisaf_elements_dataset=DataSet("ice_osisaf_elements");
    M_ice_osisaf_type_elements_dataset=DataSet("ice_osisaf_type_elements");
    M_ice_amsr2_elements_dataset=DataSet("ice_amsr2_elements");
    M_ice_nic_elements_dataset=DataSet("ice_nic_elements");
    M_ice_nic_weekly_elements_dataset=DataSet("ice_nic_weekly_elements");
    M_ice_cs2_smos_elements_dataset=DataSet("ice_cs2_smos_elements");
    M_ice_smos_elements_dataset=DataSet("ice_smos_elements");

    // datasets that need to be re-interpolated after regridding
    // - not needed if only used at initialisation, or if not interpolated onto
    // mesh (eg wave datasets are interpolated onto a fixed grid)
    M_datasets_regrid.push_back(&M_atmosphere_nodes_dataset);
    M_datasets_regrid.push_back(&M_atmosphere_elements_dataset);
    M_datasets_regrid.push_back(&M_atmosphere_bis_elements_dataset);
    M_datasets_regrid.push_back(&M_ocean_nodes_dataset);
    M_datasets_regrid.push_back(&M_ocean_elements_dataset);
}//initDatasets

//------------------------------------------------------------------------------------------------------
//! Loads and checks on the loading of various datasets.
//! * for external data to be interpolated onto the mesh elements use
//!   RX = M_mesh.bCoordX(), RY = M_mesh.bCoordY()
//! * for external data to be interpolated onto the mesh nodes use
//!   RX = M_mesh.coordX(), RY = M_mesh.coordY()
//! * for external data to be interpolated onto the WIM elements use
//!   RX = M_wim.getX(), RY = M_wim.getY()
//! * NB we don't rotate RX, RY yet since rotation angle is not always defined at this point
//! Called by checkReloadMainDatasets(), and all the ice initialisation and assimilation routines.
void
DataSetManager::checkReloadDatasets(external_data_vec const& ext_data_vec,
        double const CRtime, std::vector<double> &RX, std::vector<double> &RY)
{
    M_timer->tick("checkReloadDatasets");
    if ( ext_data_vec.size()==0 )
    {
        LOG(DEBUG) <<"checkReloadDatasets - nothing to do\n";
        return;
    }

    //loop over ext_data_vec and call check and reload for each:
    int i = 0;
    for ( auto it = ext_data_vec.begin(); it != ext_data_vec.end(); ++it, ++i )
    {
        LOG(DEBUG) <<"checkReloadDatasets for variable "
            << i << ": " << (*it)->getVariableName()
            << " of dataset " << (*it)->getDatasetName() << "\n";

        M_timer->tick((*it)->getDatasetName());
        std::string msg = "checkReloadDatasets: ExternalData object "
                + (*it)->getDatasetName() + " is not initialised yet";
        if(!(*it)->isInitialized())
            throw std::runtime_error(msg);
#ifdef OASIS
        (*it)->check_and_reload(RX, RY, CRtime, M_comm, pcpt*time_step, cpl_time_step);
#else
        (*it)->check_and_reload(RX, RY, CRtime);
#endif
        M_timer->tock((*it)->getDatasetName());
    }
    M_timer->tock("checkReloadDatasets");
}//checkReloadDatasets

//------------------------------------------------------------------------------------------------------
//! Loads and checks on the loading of the time-dependant datasets.
//! * In practice this is done by looping over the ExternalData objects in
//!   M_external_data_elements and M_external_data_nodes and checking if the corresponding Dataset
//!   needs to be reloaded and/or reinterpolated
//! Called by init() and step()
void
DataSetManager::checkReloadMainDatasets(double const CRtime)
{
    // check the time-dependant ExternalData objects to see if they need to be reloaded
    // - mesh elements
    M_timer->tick("bCoord");
    auto RX = M_mesh.bCoordX();
    auto RY = M_mesh.bCoordY();
    M_timer->tock("bCoord");
    LOG(DEBUG) <<"checkReloadDatasets (time-dependant elements)\n";
    this->checkReloadDatasets(M_external_data_elements, CRtime, RX, RY);

    // - mesh nodes
    M_timer->tick("Coord");
    RX = M_mesh.coordX();
    RY = M_mesh.coordY();
    M_timer->tock("Coord");
    LOG(DEBUG) <<"checkReloadDatasets (time-dependant nodes)\n";
    this->checkReloadDatasets(M_external_data_nodes, CRtime, RX, RY);
}//checkReloadMainDatasets

//------------------------------------------------------------------------------------------------------
//! Sets the physical variables relevant to the atmosphere according to the chosen atmospheric forcing data (CONSTANT, ASR, ERAi, ...)
//! Called by the initExternalData() function.
void
DataSetManager::forcingAtmosphere()
{
    double air_temperature_correction=vm["forecast.air_temperature_correction"].as<double>();

    switch (M_atmosphere_type)
    {
        case setup::AtmosphereType::CONSTANT:
            M_wind=ExternalData(
                vm["ideal_simul.constant_wind_u"].as<double>(),
                vm["ideal_simul.constant_wind_v"].as<double>(),
                time_init, M_spinup_duration);

            M_tair=ExternalData(vm["ideal_simul.constant_tair"].as<double>());
            M_mixrat=ExternalData(vm["ideal_simul.constant_mixrat"].as<double>());
            M_mslp=ExternalData(vm["ideal_simul.constant_mslp"].as<double>());
            M_Qsw_in=ExternalData(vm["ideal_simul.constant_Qsw_in"].as<double>());
            if(!vm["thermo.use_parameterised_long_wave_radiation"].as<bool>())
                M_Qlw_in=ExternalData(vm["ideal_simul.constant_Qlw_in"].as<double>());
            else
                throw std::runtime_error("parameterised long wave radiation not implemented for setup.atmosphere-type=constant. Use thermo.use_parameterised_long_wave_radiation=false");
            M_snowfr=ExternalData(vm["ideal_simul.constant_snowfr"].as<double>());
            M_precip=ExternalData(vm["ideal_simul.constant_precip"].as<double>());
            M_dair=ExternalData(vm["ideal_simul.constant_dair"].as<double>());
        break;

        case setup::AtmosphereType::GENERIC_PS:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, M_spinup_duration);

            M_tair=ExternalData(&M_atmosphere_elements_dataset, M_mesh, 0, false,
                    time_init, 0, air_temperature_correction);
            M_dair=ExternalData(&M_atmosphere_elements_dataset, M_mesh, 1, false,
                    time_init, 0, air_temperature_correction);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            M_snowfall=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
        break;

        case setup::AtmosphereType::ASR:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, M_spinup_duration);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_mixrat=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            if(!vm["thermo.use_parameterised_long_wave_radiation"].as<bool>())
                M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            else
                throw std::runtime_error("parameterised long wave radiation not implemented for setup.atmosphere-type=asr. Use thermo.use_parameterised_long_wave_radiation=false");
            M_snowfr=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
        break;

        case setup::AtmosphereType::ERAi:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0,true ,
                time_init, M_spinup_duration);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            if(vm["thermo.use_parameterised_long_wave_radiation"].as<bool>())
                M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            else
                throw std::runtime_error("long wave radiation not implemented for setup.atmosphere-type=erai. Use thermo.use_parameterised_long_wave_radiation=true");
            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_snowfall=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
        break;

        case setup::AtmosphereType::ERA5:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0,true ,
                time_init, M_spinup_duration);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            if(!vm["thermo.use_parameterised_long_wave_radiation"].as<bool>())
                M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            else
                throw std::runtime_error("parameterised long wave radiation not implemented for setup.atmosphere-type=ERA5. Use thermo.use_parameterised_long_wave_radiation=false");
            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_snowfall=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
        break;

        case setup::AtmosphereType::EC2:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, M_spinup_duration);

            M_tair=ExternalData(&M_atmosphere_elements_dataset, M_mesh, 0, false,
                    time_init, 0, air_temperature_correction);
            M_dair=ExternalData(&M_atmosphere_elements_dataset, M_mesh, 1, false,
                    time_init, 0, air_temperature_correction);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            if(!vm["thermo.use_parameterised_long_wave_radiation"].as<bool>())
                M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            else
                M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
        break;

        case setup::AtmosphereType::EC_ERAi:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, M_spinup_duration);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            if(vm["thermo.use_parameterised_long_wave_radiation"].as<bool>())
                M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            else
                throw std::runtime_error("long wave radiation not implemented for setup.atmosphere-type=ec_erai. Use thermo.use_parameterised_long_wave_radiation=true");
            M_Qsw_in=ExternalData(&M_atmosphere_bis_elements_dataset,M_mesh,3,false,time_init);
            M_precip=ExternalData(&M_atmosphere_bis_elements_dataset,M_mesh,5,false,time_init);
        break;

        case setup::AtmosphereType::CFSR_HI:
        case setup::AtmosphereType::CFSR:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, M_spinup_duration);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_sphuma=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            if(!vm["thermo.use_parameterised_long_wave_radiation"].as<bool>())
                M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            else
                throw std::runtime_error("parameterised long wave radiation not implemented for setup.atmosphere-type=cfsr/cfsr_hi");
            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_snowfr=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
        break;

        case setup::AtmosphereType::EC2_AROME:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, M_spinup_duration);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_sphuma=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            M_snowfall=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
        case setup::AtmosphereType::EC2_AROME_ENSEMBLE:
            M_wind=ExternalData( &M_atmosphere_nodes_dataset, M_mesh, 0 ,true ,
                time_init, M_spinup_duration, 0, M_ensemble_member);
            M_tair=ExternalData(&M_atmosphere_elements_dataset, M_mesh, 0, false,
                    time_init, 0, 0, M_ensemble_member);
            M_sphuma=ExternalData(&M_atmosphere_elements_dataset, M_mesh, 1, false,
                    time_init, 0, 0, M_ensemble_member);
            M_mslp=ExternalData(&M_atmosphere_elements_dataset, M_mesh, 2, false,
                    time_init, 0, 0, M_ensemble_member);
            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset, M_mesh, 3, false,
                    time_init, 0, 0, M_ensemble_member);
            M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset, M_mesh, 4, false,
                    time_init, 0, 0, M_ensemble_member);
            M_snowfall=ExternalData(&M_atmosphere_elements_dataset, M_mesh, 5, false,
                    time_init, 0, 0, M_ensemble_member);
            M_precip=ExternalData(&M_atmosphere_elements_dataset, M_mesh, 6, false,
                    time_init, 0, 0, M_ensemble_member);
        break;

        default:
            std::cout << "invalid atmospheric forcing"<<"\n";
            throw std::logic_error("invalid atmospheric forcing");
    }

    // add the external data objects to M_external_data_nodes or M_external_data_elements
    // for looping
    // - these are common to all the forcings
    // - check if they are initialised here
    M_external_data_nodes.push_back(&M_wind);
    M_external_data_nodes_names.push_back("M_wind");
    if(!M_wind.isInitialized())
        throw std::logic_error("M_wind is not initialised");

    int i = M_external_data_elements.size();
    M_external_data_elements_names.push_back("M_tair");
    M_external_data_elements.push_back(&M_tair);
    M_external_data_elements_names.push_back("M_mslp");
    M_external_data_elements.push_back(&M_mslp);
    M_external_data_elements_names.push_back("M_precip");
    M_external_data_elements.push_back(&M_precip);

    // specific error for M_Qsw_in
    M_external_data_elements_names.push_back("M_Qsw_in");
    M_external_data_elements.push_back(&M_Qsw_in);
    if(!M_Qsw_in.isInitialized())
        throw std::logic_error(
                "The function approxSW is not yet implemented, you need to initialize M_Qsw_in");

    // either need the long wave input, or cloud cover to parameterise it
    if(M_Qlw_in.isInitialized())
    {
        M_external_data_elements_names.push_back("M_Qlw_in");
        M_external_data_elements.push_back(&M_Qlw_in);
    }
    else if(M_tcc.isInitialized())
    {
        M_external_data_elements_names.push_back("M_tcc");
        M_external_data_elements.push_back(&M_tcc);
    }
    else
        throw std::runtime_error("forcingAtmosphere: One of M_Qlw_in or M_tcc should be initialised");

    // - snowfall can come from M_snowfall, M_snowfr*M_precip, or just M_precip (if M_tair<0)
    if(M_snowfall.isInitialized())
    {
        M_external_data_elements_names.push_back("M_snowfall");
        M_external_data_elements.push_back(&M_snowfall);
    }
    else if (M_snowfr.isInitialized())
    {
        M_external_data_elements_names.push_back("M_snowfr");
        M_external_data_elements.push_back(&M_snowfr);
    }

    if(M_sphuma.isInitialized())
    {
        // have specific humidity from the forcing
        M_external_data_elements_names.push_back("M_sphuma");
        M_external_data_elements.push_back(&M_sphuma);
    }
    else if(M_mixrat.isInitialized())
    {
        // have mixing ratio (simple relationship to specific humidity) from the forcing
        M_external_data_elements_names.push_back("M_mixrat");
        M_external_data_elements.push_back(&M_mixrat);
    }
    else if(M_dair.isInitialized())
    {
        // need to estimate the specific humidity from the dew point
        M_external_data_elements_names.push_back("M_dair");
        M_external_data_elements.push_back(&M_dair);
    }
    else
        throw std::runtime_error("forcingAtmosphere: One of M_sphuma, M_mixrat or M_dair should be initialised");

}//forcingAtmosphere

//------------------------------------------------------------------------------------------------------
//! Sets the physical variables relevant to the ocean according to the chosen ocean state and data (CONSTANT, TOPAZR, ...)
//! Called by the initExternalData() function.
void
DataSetManager::forcingOcean()//(double const& u, double const& v)
{

    bool use_ocean_nesting = false;
    if(M_use_nesting)
    {
        if(M_use_ocean_nesting)
        {
            use_ocean_nesting = true;
            M_ocean_temp=ExternalData(&M_nesting_ocean_elements_dataset, M_mesh, 0,false,time_init);
            M_ocean_salt=ExternalData(&M_nesting_ocean_elements_dataset, M_mesh, 1,false,time_init);
        }
    }

    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
            M_ocean=ExternalData(
                vm["ideal_simul.constant_ocean_u"].as<double>(),
                vm["ideal_simul.constant_ocean_v"].as<double>(),
                time_init, M_spinup_duration);

            M_ssh=ExternalData(vm["ideal_simul.constant_ssh"].as<double>(),
                time_init, M_spinup_duration);

            if (!use_ocean_nesting)
            {
                M_ocean_temp=ExternalData(physical::ocean_freezing_temp);
                M_ocean_salt=ExternalData(physical::ocean_freezing_temp/physical::mu);
            }

            M_mld=ExternalData(vm["ideal_simul.constant_mld"].as<double>());
            break;

        case setup::OceanType::TOPAZR:
        case setup::OceanType::TOPAZF:
        case setup::OceanType::GLORYS12R:
            M_ocean=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 0, true,
                time_init, M_spinup_duration);

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
                time_init, M_spinup_duration);

            if (!use_ocean_nesting)
            {
                M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
                M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
            }

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            break;
#ifdef OASIS
        case setup::OceanType::COUPLED:
            M_ocean=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 0, true,
                time_init, M_spinup_duration);

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
                time_init, M_spinup_duration);

            M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
            M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
            M_qsrml=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            if ( vm["coupler.rcv_first_layer_depth"].as<bool>() )
                M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 3,false,time_init);
            break;
#endif
        case setup::OceanType::TOPAZR_ALTIMETER:
            M_ocean=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 0, true,
                time_init, M_spinup_duration);

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
                time_init, M_spinup_duration);

            if (!use_ocean_nesting)
            {
                M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
                M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
            }

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            break;

        case setup::OceanType::TOPAZR_atrest:
            M_ocean=ExternalData(
                vm["ideal_simul.constant_ocean_u"].as<double>(),
                vm["ideal_simul.constant_ocean_v"].as<double>(),
                time_init, M_spinup_duration);

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
                time_init, M_spinup_duration);

            if (!use_ocean_nesting)
            {
                M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
                M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
            }

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            // SYL: there was a capping of the mld at minimum vm["ideal_simul.constant_mld"].as<double>()
            // but Einar said it is not necessary, so it is not implemented
            break;

        default:
            std::cout << "invalid ocean forcing"<<"\n";
            throw std::logic_error("invalid ocean forcing");
    }

    // add the external data objects to M_external_data_nodes or M_external_data_elements
    // for looping
    // - nodes
    M_external_data_nodes.push_back(&M_ocean);
    M_external_data_nodes_names.push_back("M_ocean");
    M_external_data_nodes.push_back(&M_ssh);
    M_external_data_nodes_names.push_back("M_ssh");

    // - elements
    M_external_data_elements_names.push_back("M_ocean_temp");
    M_external_data_elements.push_back(&M_ocean_temp);
    M_external_data_elements_names.push_back("M_ocean_salt");
    M_external_data_elements.push_back(&M_ocean_salt);
    if (M_mld.isInitialized())
    {
        M_external_data_elements_names.push_back("M_mld");
        M_external_data_elements.push_back(&M_mld);
    }
    if ( M_qsrml.isInitialized() )
    {
        M_external_data_elements.push_back(&M_qsrml);
        M_external_data_elements_names.push_back("M_qsrml");
    }
}//forcingOcean

// -----------------------------------------------------------------------------------------------------------
//! Initializes the bathymetry according to the type used (constant or ETOPO).
//! Called by the step() function.
void
DataSetManager::initBathymetry()//(double const& u, double const& v)
{
    // This always needs to be done, regardless of which bathymetry type we
    // have as the object M_bathymetry_elements_dataset must be initialised. But
    // if we use CONSTANT then we don't put any data into the object.
    M_bathymetry_elements_dataset=DataSet("etopo_elements");

    switch (M_bathymetry_type)
    {
        case setup::BathymetryType::CONSTANT:
            M_element_depth=ExternalData(vm["ideal_simul.constant_bathymetry"].as<double>());
            M_external_data_elements.push_back(&M_element_depth);
            break;
        case setup::BathymetryType::ETOPO:
            M_element_depth=ExternalData(&M_bathymetry_elements_dataset,M_mesh,0,false,time_init);
            M_external_data_elements.push_back(&M_element_depth);
            break;
        default:
            std::cout << "invalid bathymetry"<<"\n";
            throw std::logic_error("invalid bathymetry");
    }

    M_external_data_elements_names.push_back("M_element_depth");
    M_datasets_regrid.push_back(&M_bathymetry_elements_dataset);//this needs to be reloaded if we are regridding
}//initBathymetry

//------------------------------------------------------------------------------------------------------
//! Sets the physical variables relevant to the couple wave model
//! Called by the initExternalData() function.
#ifdef OASIS
void
DataSetManager::forcingWaves()//(double const& u, double const& v)
{

    M_wlbk = ExternalData(&M_wave_elements_dataset, M_mesh, 0, false,
                time_init, 0);
    M_external_data_elements.push_back(&M_wlbk);
    M_external_data_elements_names.push_back("M_wlbk");

    if(M_recv_wave_stress)
    {
        M_tau_wi = ExternalData(&M_wave_nodes_dataset, M_mesh, 0, true,
                    time_init, M_spinup_duration);
        M_external_data_nodes.push_back(&M_tau_wi);
        M_external_data_nodes_names.push_back("M_tauwi");
    }
    //M_str_var = ExternalData(&M_wave_elements_dataset, M_mesh, 0, false,
    //            time_init, 0);
    //M_external_data_elements.push_back(&M_str_var);
    //M_external_data_elements_names.push_back("M_str_var");
    //M_tm02 = ExternalData(&M_wave_elements_dataset, M_mesh, 1, false,
    //            time_init, 0);
    //M_external_data_elements.push_back(&M_tm02);
    //M_external_data_elements_names.push_back("M_tm02");
}
#endif






}
