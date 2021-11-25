/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   externaldata.cpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#include <externaldata.hpp>
#include <date.hpp>

#include <string>
#include <iostream>
#include <fstream>

extern "C"
{
#include <mapx.h>
}
#if defined OASIS
#include<oasis_cpp_interface.h>
#endif


/**
 * @class ExternalData
 * @brief Manage the external data
 *
 * @see
 *
*/

namespace Nextsim
{

ExternalData::ExternalData( ):
M_initialized(false)
{}

ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId,
        bool is_vector, double StartingTime )//5 args
	:
    M_is_constant( false ),
    M_dataset( dataset ),
    M_VariableId( VariableId ),
    M_bias_correction( 0. ),
    M_is_vector( is_vector ),
    M_current_time( 0. ),
    M_StartingTime( StartingTime ),
    M_SpinUpDuration( 0. ),
    M_initialized(true),
    M_log_level(Environment::logLevel()),
    M_log_all(Environment::logAll()),
    M_comm(Environment::comm())
{
    M_datasetname = (boost::format( "%1%...%2%" )
                    % M_dataset->grid.prefix
                    % M_dataset->grid.postfix
                    ).str();
    fcoeff.resize(2);
}


ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId, bool is_vector,
        double StartingTime, double SpinUpDuration )//6 args
	:
    ExternalData(dataset, mesh, VariableId, is_vector, StartingTime )
    {
        M_SpinUpDuration= SpinUpDuration ;
    }

ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId, bool is_vector,
        double StartingTime, double SpinUpDuration, double bias_correction )//7 args
    :
    ExternalData(dataset, mesh, VariableId, is_vector, StartingTime, SpinUpDuration )
    {
        M_bias_correction= bias_correction ;
    }

ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId, bool is_vector,
        double StartingTime, double SpinUpDuration, double bias_correction,
        int const& ensemble_member )//8 args
    :
    ExternalData(dataset, mesh, VariableId, is_vector,
            StartingTime, SpinUpDuration, bias_correction )
    {
        if(ensemble_member<=0)
            throw std::runtime_error(
                    "option statevector.ensemble_member<=0; should be >= 1");
        M_ensemble_member = ensemble_member ;
    }

ExternalData::ExternalData( double ConstantValue )
	:
    M_is_constant( true ),
    M_constant_value( ConstantValue ),
    M_bias_correction( 0. ),
    M_is_vector( false ),
    M_current_time( 0. ),
    M_StartingTime( 0. ),
    M_SpinUpDuration( 0. ),
    M_initialized(true),
    M_log_level(Environment::logLevel()),
    M_log_all(Environment::logAll()),
    M_comm(Environment::comm())
    {}

ExternalData::ExternalData( double ConstantValue, double ConstantValuebis )
    :
    ExternalData( ConstantValue )
    {
        M_constant_valuebis= ConstantValuebis ;
        M_is_vector= true ;
    }

ExternalData::ExternalData( double ConstantValue, double StartingTime, double SpinUpDuration )
    :
    ExternalData( ConstantValue )
    {
        M_StartingTime= StartingTime ;
        M_SpinUpDuration= SpinUpDuration ;
    }

ExternalData::ExternalData( double ConstantValue, double ConstantValuebis,
        double StartingTime, double SpinUpDuration )
    :
    ExternalData( ConstantValue, ConstantValuebis )
    {
        M_StartingTime= StartingTime ;
        M_SpinUpDuration= SpinUpDuration ;
    }

ExternalData::~ExternalData()
{
	this->clear();
}

// add ensemble_member for ensemble run, otherwise, add ensemble_member in M_dataset in finiteelement.cpp
void ExternalData::check_and_reload(std::vector<double> const& RX_in,
            std::vector<double> const& RY_in, const double current_time, const int ensemble_member)
#ifdef OASIS
{
    Communicator comm;
    this->check_and_reload(RX_in, RY_in, current_time, comm, -1., -1.);
}

void ExternalData::check_and_reload(std::vector<double> const& RX_in,
            std::vector<double> const& RY_in, const double current_time, 
            Communicator comm, const int cpl_time, const int cpl_dt)
#endif
{
    M_target_size = RX_in.size();
    M_ensemble_member = ensemble_member;
    // Don't need to do nothing more for a constant dataset
    if (M_is_constant) return;

    M_current_time = current_time;

    double current_time_tmp=M_current_time;

    M_factor=1.;
    if((M_current_time-M_StartingTime)<M_SpinUpDuration)
        M_factor=(M_current_time-M_StartingTime)/M_SpinUpDuration;

    bool to_be_reloaded=false;

    if(M_dataset->grid.dataset_frequency=="constant")
        to_be_reloaded=!M_dataset->loaded;
    else if(M_dataset->grid.dataset_frequency=="nearest_daily")
        to_be_reloaded=(
                Nextsim::datenumToString(current_time, "%Y%m%d") != Nextsim::datenumToString(M_dataset->ftime_range[0], "%Y%m%d")
                || !M_dataset->loaded);
#ifdef OASIS
    else if(M_dataset->grid.dataset_frequency=="coupled")
        to_be_reloaded=((cpl_time < M_dataset->itime_range[0] ) || (M_dataset->itime_range[1] <= cpl_time) || !M_dataset->loaded );
#endif
    else
        to_be_reloaded=((current_time_tmp < M_dataset->ftime_range[0]) || (M_dataset->ftime_range[1] < current_time_tmp) || !M_dataset->loaded);

    if (to_be_reloaded)
    {
#ifdef OASIS
        // We call oasis_get every time step, but only actually receive data at coupling times
        if (M_dataset->coupled)
        {
            if(!M_dataset->grid.loaded)
            {
                // ---------------------------------
                // Define the mapping and rotation_angle
                mapx_class *mapNextsim;
                std::string configfileNextsim = (boost::format( "%1%/%2%" )
                                          % Environment::nextsimMeshDir().string()
                                          % Environment::vm()["mesh.mppfile"].as<std::string>()
                                          ).str();

                std::vector<char> strNextsim(configfileNextsim.begin(), configfileNextsim.end());
                strNextsim.push_back('\0');
                mapNextsim = init_mapx(&strNextsim[0]);

                if(M_dataset->grid.mpp_file!="")
                {
                    mapx_class *map;
                    std::string configfile = (boost::format( "%1%/%2%" )
                                          % Environment::nextsimMeshDir().string()
                                          % M_dataset->grid.mpp_file
                                          ).str();

                    std::vector<char> str(configfile.begin(), configfile.end());
                    str.push_back('\0');
                    map = init_mapx(&str[0]);
                    M_dataset->rotation_angle = -(mapNextsim->rotation-map->rotation)*PI/180.;

                    close_mapx(map);
                }
                else
                {
                    M_dataset->rotation_angle=0.;
                }

                // ---------------------------------
                // Projection of the mesh positions into the coordinate system of the data before the interpolation
                // (either the lat,lon projection or a polar stereographic projection with another rotation angle (for ASR))
                // we should need to that also for the TOPAZ native grid, so that we could use a gridtomesh, now we use the latlon of the TOPAZ grid

                std::vector<double> RX,RY;//size set in convertTargetXY
                M_dataset->convertTargetXY(&(M_dataset->grid), RX_in, RY_in, RX, RY,mapNextsim);

                // ---------------------------------
                // Load grid if unloaded
                // This would probably be more efficient with R(X|Y)_(max|min) ... but I didn't manage to get that to work
                M_dataset->loadGrid(mapNextsim, &(M_dataset->grid), M_StartingTime, M_current_time); //, RX_min, RX_max, RY_min, RY_max);

                // closing maps
                close_mapx(mapNextsim);
            }

            this->receiveCouplingData(M_dataset, cpl_time, comm);
            this->transformData(M_dataset);
            M_dataset->interpolated = false;
            M_dataset->itime_range[0] = cpl_time;
            M_dataset->itime_range[1] = cpl_time + cpl_dt;
            M_dataset->ftime_range[0] = M_current_time;
            M_dataset->ftime_range[1] = M_current_time + double(cpl_dt)*86400.;
        }
        else {
#endif
            LOG(DEBUG) << "Load " << M_datasetname << "\n";   
            this->loadDataset(M_dataset, RX_in, RY_in);       
            this->transformData(M_dataset);
#ifdef ENSEMBLE 
            //Todo: use float variable to save sources, since it's no need to have high precision perturbations.   
            // variables on nodes and elements are perturbed at the same time by default, but perturbation on  the elements can be turned off.
            // Load generic_atm_....nc is loaded twice for a given time, one for variables saved on elements, one for nodes. We only active perturbation after loaded atm node, since perturbation is independent on physical data.
            // todo, if it needs to avoid unexpected perturbation due to remesh process.  
            LOG(DEBUG) << "adding perturbations to loaded data\n"; 
            M_comm.barrier();
            if (strcmp (M_dataset->name.c_str(), "topaz_forecast_elements") == 0 || \
                strcmp (M_dataset->name.c_str(), "asr_nodes") == 0 || \
                strcmp (M_dataset->name.c_str(), "generic_atm_nodes") == 0 || \
                strcmp (M_dataset->name.c_str(), "asr_elements") == 0 || \
                strcmp (M_dataset->name.c_str(), "generic_atm_elements") == 0 )    
            {
                ensemble perturbation; 
                int M_full  = 1024; //std::pow(2,std::ceil(std::log2(M_dataset->grid.dimension_y_count_netcdf)));
                int N_full  = 1024; //std::pow(2,std::ceil(std::log2(M_dataset->grid.dimension_x_count_netcdf)));
                int MN_full = M_full*N_full;
                int y_start = M_dataset->grid.dimension_y_start;
                int x_start = M_dataset->grid.dimension_x_start;
                int y_count = M_dataset->grid.dimension_y_count;
                int x_count = M_dataset->grid.dimension_x_count; 
                std::vector<double> synforc1(1);
                std::vector<double> synforc2(1);
                synforc1.resize(MN_full);
                synforc2.resize(MN_full);

                LOG(DEBUG) << "### MN_FULL: " << MN_full <<" = " << M_full << "x" <<N_full<<"\n"; 
                LOG(DEBUG) << "### M_dataset_name: " << M_dataset->name << "\n";
                LOG(DEBUG) << "### M_current_time: " << M_current_time  << " = "<<datenumToString(M_current_time)<<"\n";

                // apply ocean perturbation
                // loaded_data[0]: previous, loaded_data[1]: current
                // int opr =1: +  //for most of the perturbed variables, =2: *// for variables using lognormal format, refered to rand_update() in mod_random_forcing.F90
                if (strcmp (M_dataset->name.c_str(), "topaz_forecast_elements") == 0 )   
                { 
                    M_dataset->N_ocean = floor(M_current_time - M_StartingTime + 0.5); // update file after 12:00
                    for (int it = 0; it<2; it++) //fstep
                    {
                        if (M_comm.rank() == 0) {  
                            // Perturbations_XX_perturbation_fileID.nc are linked from offline, generated by module/offline perturbations
                            // file series starts from N_ocean=0
                            // grid->dirname
                            // boost::filesystem::path const statevector_restart_input_path ( Environment::vm()["statevector.restart_path"].as<std::string>() );
                            std::string filename = Environment::vm()["statevector.restart_path"].as<std::string>() + "/Perturbations/OceanPerturbations_mem" + 
                                                std::to_string(M_ensemble_member) +"_series" + std::to_string(M_dataset->N_ocean + it) + ".nc";
                            LOG(DEBUG)<<"topaz load: "<<filename<<"\n";
                            netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
                            netCDF::NcVar data;
                            data = dataFile.getVar("sst");     
                            data.getVar(&synforc1[0]);
                            data = dataFile.getVar("sss");
                            data.getVar(&synforc2[0]);   
                        }
                        M_comm.barrier();
                        boost::mpi::broadcast(M_comm, &synforc1[0], MN_full, 0); 
                        boost::mpi::broadcast(M_comm, &synforc2[0], MN_full, 0); 
                        LOG(DEBUG) << "### Add previous/current perturbations to fields\n";  //{sst,sss,mld,conc,thick, snow_thick}
                        // todo, merge this part into addPerturbation 
                        
                        int final_MN = M_dataset->grid.reduced_nodes_ind.size();
                        LOG(DEBUG) << "### final_MN: "<< final_MN << " variable size: " << M_dataset->variables[0].loaded_data[it].size() <<"\n"; 
                        for(int n = 0; n < final_MN; n++)
                        {  //because of using part of ocean grid, i,j are recalculated. 
                            int j = M_dataset->grid.reduced_nodes_ind[n]/x_count;
                            int i = M_dataset->grid.reduced_nodes_ind[n] - j*x_count;
                            i = x_start + (y_start + j)*N_full + i;
                            M_dataset->variables[0].loaded_data[it][n] += synforc1[i]*( 1- M_dataset->variables[3].loaded_data[it][n]); //limit sst perturbation to ice covered area, variables[3] is conc, variables[0] is sst
                            M_dataset->variables[1].loaded_data[it][n] += synforc2[i]; 
                            if (i >= MN_full) {LOG(DEBUG)<<"maxmium synforc exceeded\n";}
                        }
                    }
                }                
                
                // wind velocity perturbation 
                if (strcmp (M_dataset->name.c_str(), "asr_nodes") == 0 || \
                    strcmp (M_dataset->name.c_str(), "generic_atm_nodes") == 0 )
                {   
                    for (int it = 0; it<2; it++)
                    {
                        if (M_comm.rank() == 0) {  
                            M_dataset->N_wind = floor((M_current_time - M_StartingTime)*4); // update files every 6 hours. unit is a day.
                            std::string filename = Environment::vm()["statevector.restart_path"].as<std::string>() + "/Perturbations/AtmospherePerturbations_mem" + std::to_string(M_ensemble_member) +"_series" + std::to_string(M_dataset->N_wind+it) + ".nc";                     
                            LOG(DEBUG)<<"wind node load: "<<filename<<"\n";
                            netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
                            netCDF::NcVar data;     
                            data = dataFile.getVar("uwind");   
                            data.getVar( &synforc1[0]);
                            data = dataFile.getVar("vwind");
                            data.getVar( &synforc2[0]);     
                        }
                        M_comm.barrier();
                        boost::mpi::broadcast(M_comm, &synforc1[0], MN_full, 0); 
                        boost::mpi::broadcast(M_comm, &synforc2[0], MN_full, 0); 
                        LOG(DEBUG) << "### Add previous/current perturbations to fields\n";  
                        perturbation.addPerturbation(M_dataset->variables[0].loaded_data[it], synforc1, M_full,N_full, x_start, y_start, x_count, y_count, 1);  
                        perturbation.addPerturbation(M_dataset->variables[1].loaded_data[it], synforc2, M_full,N_full, x_start, y_start, x_count, y_count, 1); 
                    }
                }
                // atmosphere scalar variable perturbation
                else if (strcmp (M_dataset->name.c_str(), "asr_elements") == 0 || \
                        strcmp (M_dataset->name.c_str(), "generic_atm_elements") == 0 )
                {   // index in M_dataset->variables[0] indicates to snowfall defined in dataset.cpp,generic_atm_elements,{ tair, dair, mslp, Qsw_in, Qlw_in, snowfall, precip }
                    for (int it = 0; it<2; it++)
                    {
                        if (M_comm.rank() == 0) {  
                            M_dataset->N_wind = floor((M_current_time - M_StartingTime)*4); // update files every 6 hours. unit is a day.
                            
                            std::string filename = Environment::vm()["statevector.restart_path"].as<std::string>() + "/Perturbations/AtmospherePerturbations_mem" + std::to_string(M_ensemble_member) +"_series" + std::to_string(M_dataset->N_wind+it) + ".nc";
                            LOG(DEBUG)<<"wind element load: "<<filename<<"\n";
                            netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
                            netCDF::NcVar data;
                            data = dataFile.getVar("Qlw_in");    //longwave downwelling radiation rate                      
                            data.getVar( &synforc1[0]);
                            data = dataFile.getVar("snowfall");  //snowfall rate
                            data.getVar( &synforc2[0]);
                        }
                        M_comm.barrier();
                        boost::mpi::broadcast(M_comm, &synforc1[0], MN_full, 0); 
                        boost::mpi::broadcast(M_comm, &synforc2[0], MN_full, 0); 
                        LOG(DEBUG) << "### Add previous/current perturbations to fields\n";  
                        // int opr =1: +  //for most of the perturbed variables, =2: *// for variables using lognormal format, refered to rand_update() in mod_random_forcing.F90
                        perturbation.addPerturbation(M_dataset->variables[4].loaded_data[it], synforc1, M_full,N_full, x_start, y_start, x_count, y_count, 1);  
                        perturbation.addPerturbation(M_dataset->variables[5].loaded_data[it], synforc2, M_full,N_full, x_start, y_start, x_count, y_count, 2); 
                    }
                    M_dataset->N_wind++; 
                }
                M_comm.barrier();
            }
#endif
        //interpolate the loaded external data onto the model grid
            M_dataset->interpolated = false;            
#ifdef OASIS
        }
#endif
    }
#ifdef OASIS
    if ( M_dataset->calc_nodal_weights )
    {
        LOG(DEBUG) << "set nodal weights for " << M_dataset->name << "\n";
        M_dataset->setNodalWeights(RX_in, RY_in);
        M_dataset->calc_nodal_weights = false;
        M_dataset->interpolated = false;
    }
#endif

    if (!M_dataset->interpolated)
    {
        LOG(DEBUG) << "Interpolate " << M_datasetname << "\n";
        this->interpolateDataset(M_dataset, RX_in, RY_in);
        LOG(DEBUG) << "Done\n";
    }
}

typename ExternalData::size_type
ExternalData::size()
{
    if ( M_is_vector )
        return static_cast<size_type>( 2*M_target_size );
    else
        return static_cast<size_type>( M_target_size );
}

typename ExternalData::value_type
ExternalData::operator [] (const size_type i)
{
    return static_cast<value_type>( get(i) );
}

typename ExternalData::value_type
ExternalData::get(const size_type i)
{
    value_type value;
    size_type i_tmp;
    int VariableId_tmp;
    if(M_is_constant)
    {
        // for the moment same value is given to all the components
        if(!M_is_vector)
        {
            ASSERT(i < M_target_size,
                    "invalid index for scalar (no dataset)");
            value = M_constant_value;
        }
        else
        {
            ASSERT(i < 2*M_target_size,
                    "invalid index for vector (no dataset)");

            if(i < M_target_size)
                value = M_constant_value;
            else
                value = M_constant_valuebis;
        }
    }
    else
    {
        // determine whether to use linear interpolation in time
        // or constant/step-wise interpolation
        bool interp_linear_time = (M_dataset->grid.dataset_frequency!="constant"
#ifdef OASIS
                && M_dataset->grid.dataset_frequency!="coupled"
#endif
                && M_dataset->grid.dataset_frequency!="nearest_daily");
        bool interp_const_wave  = (M_dataset->grid.waveOptions.wave_dataset);
#if defined(WAVES)
        interp_const_wave = (interp_const_wave
                && (Environment::vm()["wimsetup.wave-time-interp-option"].as<std::string>()=="step"));
#endif
                //&& (M_dataset->grid.waveOptions.time_interp_option=="step"));
        interp_linear_time  = (interp_linear_time && !interp_const_wave);

        if(interp_linear_time)
        {
            fdt = std::abs(M_dataset->ftime_range[1]-M_dataset->ftime_range[0]);
            fcoeff[0] = std::abs(M_current_time-M_dataset->ftime_range[1])/fdt;
            fcoeff[1] = std::abs(M_current_time-M_dataset->ftime_range[0])/fdt;

            //std::cout <<"LINEAR COEFF 1= "<< fcoeff[0] <<"\n";
            //std::cout <<"LINEAR COEFF 2= "<< fcoeff[1] <<"\n";

            if(!M_is_vector)
            {
                ASSERT(i < M_target_size, "invalid index for scalar id = "
                        + std::to_string(M_VariableId)+" ("+M_datasetname+"): i = "
                        + std::to_string(i)+"; M_target_size = "
                        + std::to_string(M_target_size));
                i_tmp=i;
                VariableId_tmp=M_VariableId;
            }
            else
            {
                ASSERT(i < 2*M_target_size, "invalid index for vector id = "
                        +std::to_string(M_VariableId)+" ("+M_datasetname+")");
                if(i < M_target_size)
                {
                    i_tmp=i;
                    VariableId_tmp=M_dataset->vectorial_variables[M_VariableId].components_Id[0];
                }
                else
                {
                    i_tmp=i-M_target_size;
                    VariableId_tmp=M_dataset->vectorial_variables[M_VariableId].components_Id[1];
                }
            }

            //time interpolation
            value =  M_factor*
                (fcoeff[0]*M_dataset->variables[VariableId_tmp].interpolated_data[0][i_tmp] +
                 fcoeff[1]*M_dataset->variables[VariableId_tmp].interpolated_data[1][i_tmp]);
        }//interp linear in time
        else
        {
            //step-function or constant in time
            if(!M_is_vector)
            {
                ASSERT(i < M_target_size, "invalid index for scalar id = "
                        +std::to_string(M_VariableId)+" ("+M_datasetname+")");
                i_tmp=i;
                VariableId_tmp=M_VariableId;
            }
            else
            {
                ASSERT(i < 2*M_target_size, "invalid index for vector id = "
                        +std::to_string(M_VariableId)+" ("+M_datasetname+")");

                if(i < M_target_size)
                {
                    // x component
                    i_tmp=i;
                    VariableId_tmp=M_dataset->vectorial_variables[M_VariableId].components_Id[0];

                }
                else
                {
                    // y component
                    i_tmp=i-M_target_size;
                    VariableId_tmp=M_dataset->vectorial_variables[M_VariableId].components_Id[1];
                }
            }
            value =  M_factor*M_dataset->variables[VariableId_tmp].interpolated_data[0][i_tmp];
        }
    }

	return static_cast<value_type>( value + M_bias_correction );
}
typename std::vector<double>
ExternalData::getVector()
{
    std::vector<double> vector_tmp;

    if(M_initialized)
    {
        int size_vector = M_target_size;
        if(M_is_vector)
            size_vector *= 2;

        vector_tmp.resize(size_vector);

        for (int i=0; i<size_vector; ++i)
        {
            vector_tmp[i] = (double) get(i);
        }
    }
    return vector_tmp;
}

#ifdef OASIS
void
ExternalData::receiveCouplingData(Dataset *dataset, int cpl_time, Communicator comm)
{
        // ierror = OASIS3::get_2d(var_id[1], pcpt*time_step, &field2_recv[0], M_cpl_out.M_ncols, M_cpl_out.M_nrows);
        LOG(DEBUG) << "receiveCouplingData at cpl_time " << cpl_time << "\n";
        for(int j=0; j<dataset->variables.size(); ++j)
        {
            int M_full  = dataset->grid.dimension_y_count_netcdf;
            int N_full  = dataset->grid.dimension_x_count_netcdf;
            int MN_full = M_full*N_full;
            std::vector<double> data_in_tmp(MN_full);

            if ( comm.rank() == 0 )
                int ierror = OASIS3::get_2d(dataset->M_cpl_id[j], cpl_time, &data_in_tmp[0], N_full, M_full);

            // TODO: This is not very efficient but it will do for now, premature optimisations, and all that!
            boost::mpi::broadcast(comm, &data_in_tmp[0], MN_full, 0);

            int y_start = dataset->grid.dimension_y_start;
            int x_start = dataset->grid.dimension_x_start;
            int M       = dataset->grid.dimension_y_count;
            int N       = dataset->grid.dimension_x_count;
            int MN      = M*N;
            int final_MN = MN;

            if(dataset->grid.reduced_nodes_ind.size()!=0)
            {
                if((dataset->grid.dimension_y.cyclic) || (dataset->grid.dimension_x.cyclic))
                    throw std::runtime_error("Using reduced grid and cyclic grid at the same time is not yet implemented");

                final_MN=dataset->grid.reduced_nodes_ind.size();
            }

            dataset->variables[j].loaded_data.resize(1);
            dataset->variables[j].loaded_data[0].resize(final_MN);
            int reduced_i;
            for (int i=0; i<(final_MN); ++i)
            {
                reduced_i=i;
                if(dataset->grid.reduced_nodes_ind.size()!=0)
                    reduced_i=dataset->grid.reduced_nodes_ind[i];

                //now add to loaded_data
                double tmp_data_i = data_in_tmp[reduced_i];
                tmp_data_i *= dataset->variables[j].a;
                tmp_data_i += dataset->variables[j].b;
                dataset->variables[j].loaded_data[0][i] = tmp_data_i;
            }
        }

        dataset->nb_forcing_step=1;
        dataset->loaded=true;
}
#endif

void
//ExternalData::loadDataset(Dataset *dataset, GmshMesh const& mesh)//(double const& u, double const& v)
ExternalData::loadDataset(Dataset *dataset, std::vector<double> const& RX_in,
        std::vector<double> const& RY_in)//(double const& u, double const& v)
{
    // ---------------------------------
    // Define the mapping and rotation_angle
    mapx_class *mapNextsim;
    std::string configfileNextsim = Environment::nextsimMppfile();
    std::vector<char> strNextsim(configfileNextsim.begin(), configfileNextsim.end());
    strNextsim.push_back('\0');
    mapNextsim = init_mapx(&strNextsim[0]);

    if(dataset->grid.mpp_file!="")
    {
        mapx_class *map;
        std::string configfile = (boost::format( "%1%/%2%" )
                          % Environment::nextsimMeshDir().string()
                          % dataset->grid.mpp_file
                          ).str();

        std::vector<char> str(configfile.begin(), configfile.end());
        str.push_back('\0');
        map = init_mapx(&str[0]);
        dataset->rotation_angle = -(mapNextsim->rotation-map->rotation)*PI/180.;

        close_mapx(map);
    }
    else
    {
        dataset->rotation_angle=0.;
    }
    // ---------------------------------
    // Load grid if unloaded
    LOG(DEBUG) << dataset->grid.loaded << "\n";

    if(!dataset->grid.loaded)
    {
        bool is_topaz_fc = (dataset->grid.dataset_frequency=="daily_forecast");//topaz forecast
        bool is_ec_fc = ((dataset->grid.prefix).find("start") != std::string::npos);//ec_[nodes,elements],ec2_[nodes,elements]
        bool true_forecast = (Environment::vm()["forecast.true_forecast"].as<bool>());
        double init_time = M_StartingTime;
        if((is_ec_fc||is_topaz_fc)&&true_forecast)
        {
            if (is_ec_fc)
            {
                // use forecast.time_init_atm_fc option to get init_time
                std::string tmpstr = (Environment::vm()["forecast.time_init_atm_fc"].as<std::string>());
                if(tmpstr!="")
                    init_time = Nextsim::stringToDatenum(tmpstr);
            }
            else
            {
                // use forecast.time_init_ocean_fc option to get init_time
                std::string tmpstr = (Environment::vm()["forecast.time_init_ocean_fc"].as<std::string>());
                if(tmpstr!="")
                    init_time = Nextsim::stringToDatenum(tmpstr);
            }
        }
        //only need init_time to get grid
        dataset->loadGrid(mapNextsim, &(dataset->grid), init_time, init_time, RX_in, RY_in);  // load a portion of grid based on sea ice area RX_in,RY_in for efficiency
        // dataset->loadGrid(mapNextsim, &(dataset->grid), init_time, init_time);  //, RX_in, RY_in);
    }

    // closing maps
    close_mapx(mapNextsim);

    // ---------------------------------
    std::vector<double> XTIME(1);
    std::vector<size_t> index_start(1);
    std::vector<size_t> index_count(1);

    // size of the data
    int M        = dataset->grid.dimension_y_count;
    int N        = dataset->grid.dimension_x_count;
    int MN       = M*N;
    int final_MN = MN;

    if(dataset->grid.reduced_nodes_ind.size()!=0)
    {
        if((dataset->grid.dimension_y.cyclic) || (dataset->grid.dimension_x.cyclic))
            throw std::runtime_error("Using reduced grid and cyclic grid at the same time is not yet implemented");

    	final_MN=dataset->grid.reduced_nodes_ind.size();
    }

    // initialization of the tmp data vector
    std::vector<double> data_in_tmp(MN);

    // netcdf objects
    // Attributes (scaling and offset)
    netCDF::NcVarAtt att;
    double scale_factor;
    double add_offset;

    std::vector<netCDF::NcVar> NcVars(dataset->variables.size());

    // ---------- Automatic identification of the file and time index

    // Define variables for this scope
    double ftime, time_prev, time_next;
    std::string filename;
    std::string filename_prev="";
    std::string filename_next="";
    std::vector<int> file_jump;

    int index = 0;
    int index_prev=-1;
    int index_next=-1;

    std::vector<std::string> filename_fstep;
    std::vector<int> index_fstep;

    std::string f_timestr;
    bool is_topaz_fc = (dataset->grid.dataset_frequency=="daily_forecast");//topaz forecast
    bool is_ec_fc = ((dataset->grid.prefix).find("start") != std::string::npos);//ec_[nodes,elements],ec2_[nodes,elements]
    bool true_forecast = (Environment::vm()["forecast.true_forecast"].as<bool>());
    double init_time = M_StartingTime;
        // - for forecasts, filename depends on start time

    // Filename depends on the date for time varying data
    if(dataset->grid.dataset_frequency!="constant"
        && dataset->grid.dataset_frequency!="nearest_daily")
    {
        ftime = M_current_time-dataset->averaging_period/2.;
        file_jump ={-1,0,1};

        LOG(DEBUG)<<"LOAD DATASET TIMES:\n";
        LOG(DEBUG)<<"init_time = "<<init_time<<" = "<<datenumToString(init_time)<<"\n";
        LOG(DEBUG)<<"M_current_time = "<<M_current_time<<" = "<<datenumToString(M_current_time)<<"\n";
        LOG(DEBUG)<<"ftime = "<<ftime<<" = "<<datenumToString(ftime)<<"\n";
        if((is_ec_fc||is_topaz_fc)&&true_forecast)
        {
            // when using forcing from ECMWF or topaz forecasts, we select the file based on the StartingTime
            if (is_ec_fc)
            {
                // - one file for all records
                // - ftime not used (only init_time)
                file_jump ={0};
                std::string tmpstr = (Environment::vm()["forecast.time_init_atm_fc"].as<std::string>());
                if(tmpstr!="")
                    init_time = Nextsim::stringToDatenum(tmpstr);
            }
            else
            {
                std::string tmpstr = (Environment::vm()["forecast.time_init_ocean_fc"].as<std::string>());
                if(tmpstr!="")
                    init_time = Nextsim::stringToDatenum(tmpstr);
            }
        }//forecasts

        for (int jump: file_jump)
        {
            if(is_ec_fc||is_topaz_fc)
            {
                double inittime = init_time;
                if(!true_forecast)
                    // * if(!true_forecast), take the forecast that started at the start of
                    //   the "current day" (ftime+jump)
                    // * also can't have init_time before start of
                    //   the "current day" (ftime+jump)
                    // NB jump is in days for these datasets
                    inittime = std::floor(ftime+jump);
                filename = dataset->getFilename(&(dataset->grid),inittime,ftime+jump);
            }
            else
                filename = dataset->getFilename(&(dataset->grid),init_time,ftime,jump);

            LOG(DEBUG)<<"FILENAME (JUMP = " <<jump<< ") = "<< filename
                << ". Exists? " << boost::filesystem::exists(filename)
                <<"\n";
            if ( ! boost::filesystem::exists(filename) )
                continue;

            index_start.resize(1);
            index_count.resize(1);

            // Open the netcdf file
            netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

            netCDF::NcVar FVTIME = dataFile.getVar(dataset->time.name);
            netCDF::NcDim timeDim = dataFile.getDim(dataset->time.name);

            // Read first two or the only time coordinates
            index_start[0]=0;
            index_count[0]=std::min(2, int(timeDim.getSize()));
            XTIME.resize(index_count[0]);

            FVTIME.getVar(index_start, index_count, &XTIME[0]);

            double const t_ref = stringToDatenum(dataset->grid.reference_date);
            // If we have two time indices we calculate a delta t and derive the correct index
            if ( index_count[0] == 2 )
            {
                double const t0 = (XTIME[0]*dataset->time.a+dataset->time.b)/24.0 + t_ref;
                double const t1 = (XTIME[1]*dataset->time.a+dataset->time.b)/24.0 + t_ref;
                double const dt = t1 - t0;
                int ntimes = timeDim.getSize();
                if(is_ec_fc && !true_forecast)
                    ntimes = std::round(1/dt);//only use 1st day of file

                // This is a double because we're most likely in between integer indices
                double const indx = (M_current_time - t0)/dt;

                /* There is an extremely small chance that this procedure puts
                 * us in the wrong interval. It can only happen if the rounding
                 * error on indx results in indx_floor being smaller than the
                 * correct value. But we know that M_current_time >
                 * M_dataset->ftime_range[1], so this is almost certainly not
                 * going to happen. */
                int const indx_floor = std::floor(indx);
                int const indx_ceil  = indx_floor + 1;

                // The index below the current time
                if ( indx_floor >= 0 && indx_floor < ntimes )
                {
                    index_start[0] = indx_floor;
                    index_count[0] = 1;
                    XTIME.resize(index_count[0]);

                    FVTIME.getVar(index_start, index_count, &XTIME[0]);

                    time_prev = (XTIME[0]*dataset->time.a+dataset->time.b)/24.0 + t_ref;
                    index_prev = indx_floor;
                    filename_prev = filename;
                }

                // The index above the current time
                if ( indx_ceil >= 0 && indx_ceil < ntimes )
                {
                    index_start[0] = indx_ceil;
                    index_count[0] = 1;
                    XTIME.resize(index_count[0]);

                    FVTIME.getVar(index_start, index_count, &XTIME[0]);

                    time_next = (XTIME[0]*dataset->time.a+dataset->time.b)/24.0 + t_ref;
                    index_next = indx_ceil;
                    filename_next = filename;
                }
            }
            // If there's only one time step per file
            else if ( index_count[0] == 1 )
            {
                double const f = (XTIME[0]*dataset->time.a+dataset->time.b)/24.0 + t_ref;
                if(f>M_current_time && index_next==-1)
                {
                    time_next=f;
                    index_next = 0;
                    filename_next = filename;
                }
                if(f<=M_current_time)
                {
                    time_prev=f;
                    index_prev = 0;
                    filename_prev = filename;
                }
            }
            // WTF?
            else
            {
                throw std::runtime_error("ExternalData::loadData: Empty time variable in dataset " + dataset->name + "\n");
            }
        }//loop over jump

        if(filename_prev!="")
        {
            filename_fstep.push_back(filename_prev);
            index_fstep.push_back(index_prev);
        }
        else
        {
            std::stringstream msg;
            msg << "Dataset " << dataset->name
                << ": file for before " << datenumToString(M_current_time)
                << " not found";
            throw std::runtime_error(msg.str());
        }


        if(filename_next!="")
        {
            filename_fstep.push_back(filename_next);
            index_fstep.push_back(index_next);
        }
        else
        {
            std::stringstream msg;
            msg << "Dataset " << dataset->name
                << ": file for after " << datenumToString(M_current_time)
                << " not found";
            throw std::runtime_error(msg.str());
        }

        dataset->ftime_range = {time_prev,time_next};
    }//not nearest_daily or constant
    else
    {
        if(dataset->grid.dataset_frequency=="nearest_daily")
        {
            ftime = M_current_time;
            double f= std::floor(ftime);
            f_timestr = datenumToString(f, "%Y%m%d");
            dataset->ftime_range = {f+.5};
        }
        else
            f_timestr ="";

        filename = (boost::format( "%1%/%2%/%3%%4%%5%" )
                    % Environment::nextsimDataDir().string()
                    % dataset->grid.dirname
                    % dataset->grid.prefix
                    % f_timestr
                    % dataset->grid.postfix
                    ).str();
        filename_fstep.push_back(filename);
        index_fstep.push_back(0);
    }

    // Initialise counters etc.
    int nb_forcing_step =filename_fstep.size();
    LOG(DEBUG)<<"Start loading data\n";
    for (int fstep=0; fstep < nb_forcing_step; ++fstep) // always need one step before and one after the target time
    {
        // Load each variable and copy its data into loaded_data
        for(int j=0; j<dataset->variables.size(); ++j)
        {
            filename=filename_fstep[fstep];
            index=index_fstep[fstep];

            // Replace the "prefix" if we have one variable per file
            if ( dataset->variables[j].filename_prefix != "" )
            {
                // extract the f_timestr for this fstep ( +1 for the directory / )
                std::string::size_type start = Environment::nextsimDataDir().string().size() + 1
                    + dataset->grid.dirname.size() + 1
                    + dataset->grid.prefix.size();

                std::string::size_type end = filename.find(dataset->grid.postfix, start);
                f_timestr = filename.substr(start, end-start);

                // Re-build filename with grid.prefix replaced with variables[j].filename_prefix
                filename = (boost::format( "%1%/%2%/%3%%4%%5%" )
                            % Environment::nextsimDataDir().string()
                            % dataset->grid.dirname
                            % dataset->variables[j].filename_prefix
                            % f_timestr
                            % dataset->grid.postfix
                            ).str();
            }

            LOG(DEBUG)<<"FILENAME= "<< filename <<"\n";
            if ( ! boost::filesystem::exists(filename) )
                throw std::runtime_error("File not found: " + filename);

            // Open the netcdf file
            netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

            LOG(DEBUG)<<"variables number:" << j  << "\n";
            if ((dataset->variables[j].wavDirOptions.isWavDir)
                    &&(!dataset->variables[j].wavDirOptions.xComponent))
            {
                // mwd is a scalar to be turned into a unit vector
                // then rotated to x-y grid coords
                // - do nothing now for y-component
                // (not used yet - only for storage later)
                dataset->variables[j].loaded_data[fstep]   =
                    dataset->variables[j-1].loaded_data[fstep];
                continue;
            }

            NcVars[j] = dataFile.getVar(dataset->variables[j].name);
            index_start.resize(dataset->variables[j].dimensions.size());
            index_count.resize(index_start.size());

            // here we find the start and count index for each dimensions
            for(int k=0; k<index_start.size(); ++k)
            {
                std::string const dimension_name=dataset->variables[j].dimensions[k].name;
                std::string const dn_lower=boost::algorithm::to_lower_copy(dimension_name);
                bool is_time = false;
                if(dn_lower.size()>=4)
                    // so far we have time, Time, time0 so this should work
                    is_time = dn_lower.substr(0, 4) == "time";

                // dimension_x case
                if ((dimension_name).find(dataset->grid.dimension_x.name) != std::string::npos)
                {
                    index_start[k] = dataset->grid.dimension_x_start;
                    index_count[k] = dataset->grid.dimension_x_count;
                }
                // dimension_y case
                else if ((dimension_name).find(dataset->grid.dimension_y.name) != std::string::npos)
                {
                    index_start[k] = dataset->grid.dimension_y_start;
                    index_count[k] = dataset->grid.dimension_y_count;
                }
                //time dimension
                // - for constant/nearest_daily use index_start = 0
                //   (as in other cases)
                else if (is_time
                        && dataset->grid.dataset_frequency!="constant"
                        && dataset->grid.dataset_frequency!="nearest_daily")
                {
                    LOG(DEBUG) << "Ordinary time dimension "
                        << k << "/" << dimension_name <<"\n";
                    index_start[k] = index;
                    index_count[k] = 1;
                }
                //depth dimension
                else if (dimension_name == "depth")
                {
                    index_start[k] = 0;
                    index_count[k] = 1;
                }
                // ensemble member case
                else if (dimension_name == "ensemble_member")
                {
                    LOG(DEBUG)<<"Using ensemble member = "
                        << M_ensemble_member << "\n";
                    index_start[k] = M_ensemble_member - 1;
                    index_count[k] = 1;
                }
                // other cases
                else
                {
                    auto tmpDim = dataFile.getDim(dimension_name);
                    index_start[k] = 0;
                    index_count[k] = tmpDim.getSize();
                }
                LOG(DEBUG) << "Dimension " << k << "(" 
                    << dimension_name << "), start/count = "
                    << index_start[k] << "/" << index_count[k]
                    << "\n";
            }

            // Reading the netcdf
            NcVars[j].getVar(index_start,index_count,&data_in_tmp[0]);

            //----------- Unit transformation ------------
            // scale factor and add offset are stored as variable attributes
            scale_factor=1.;
            try
            {
                att = NcVars[j].getAtt("scale_factor");
                att.getValues(&scale_factor);
            }
            catch(netCDF::exceptions::NcException& e)
            {}

            add_offset=0.;
            try
            {
                att = NcVars[j].getAtt("add_offset");
                att.getValues(&add_offset);
            }
            catch(netCDF::exceptions::NcException& e)
            {}

            // The factor and b are defined for each variables to get the right units for the neXtSIM model
            for(double& d : data_in_tmp)
                d=(d*scale_factor + add_offset)*dataset->variables[j].a+dataset->variables[j].b;


            // Check Nan and store the loaded data
            dataset->variables[j].loaded_data[fstep].resize(final_MN);

            double tmp_data_i;
            int reduced_i;
            for (int i=0; i<(final_MN); ++i)
            {
                reduced_i=i;
                if(dataset->grid.reduced_nodes_ind.size()!=0)
                    reduced_i=dataset->grid.reduced_nodes_ind[i];

                tmp_data_i=data_in_tmp[reduced_i];
                if ((!dataset->variables[j].wavDirOptions.isWavDir)
                        && (std::isnan(tmp_data_i)))
                {
                    tmp_data_i=0.;//keep as nan for now if mwd
                }

                //now add to loaded_data
                dataset->variables[j].loaded_data[fstep][i]=tmp_data_i;
            }
        }
    }

    dataset->nb_forcing_step=nb_forcing_step;
    dataset->loaded=true;
}//loadDataset


// Transformation of the vectorial variables from the coordinate system of the data to the polar stereographic projection used in the model
// Once we are in the polar stereographic projection, we can do spatial interpolation without bothering about the North Pole
void
ExternalData::transformData(Dataset *dataset)
{

    double tmp_data0, tmp_data1, new_tmp_data0, new_tmp_data1;
    double tmp_data0_deg, tmp_data1_deg;
    double lat_tmp, lon_tmp, lat_tmp_bis, lon_tmp_bis;
    double x_tmp, y_tmp, x_tmp_bis, y_tmp_bis;
    double speed, new_speed;
    int j0, j1;

    double R=mapx_Re_km*1000.; // Earth radius (m)
    double delta_t=1.; // 1 sec. This value needs to be small.
    double delta_r=10.; // 10m. This value needs to be small.
    double delta_x,delta_y,delta_lon,delta_lat;

    // ---------------------------------
    // Define the mapping and rotation_angle
    mapx_class *mapNextsim;
    std::string configfileNextsim = Environment::nextsimMppfile();
    std::vector<char> strNextsim(configfileNextsim.begin(), configfileNextsim.end());
    strNextsim.push_back('\0');
    mapNextsim = init_mapx(&strNextsim[0]);

    double cos_m_diff_angle, sin_m_diff_angle;
    cos_m_diff_angle=std::cos(-dataset->rotation_angle);
    sin_m_diff_angle=std::sin(-dataset->rotation_angle);

    // size of the data
    int M        = dataset->grid.dimension_y_count;
    int N        = dataset->grid.dimension_x_count;
    int MN       = M*N;
    int final_MN = MN;

    if(dataset->grid.reduced_nodes_ind.size()!=0)
    {
        if((dataset->grid.dimension_y.cyclic) || (dataset->grid.dimension_x.cyclic))
            throw std::runtime_error("Using reduced grid and cyclic grid at the same time is not yet implemented");

        final_MN=dataset->grid.reduced_nodes_ind.size();
    }


    LOG(DEBUG)<<"Start transforming the data\n";
    for (int fstep=0; fstep < dataset->nb_forcing_step; ++fstep) // always need one step before and one after the target time
    {
        for(int j=0; j<dataset->vectorial_variables.size(); ++j)
        {

            j0=dataset->vectorial_variables[j].components_Id[0];
            j1=dataset->vectorial_variables[j].components_Id[1];

            if (dataset->variables[j0].wavDirOptions.isWavDir)
            {
                LOG(DEBUG)<<"is wave direction\n";
                //change from wave-from to wave-to dirn
                double disp_factor  = 1.;
                if (dataset->variables[j0].wavDirOptions.waveFrom)
                    disp_factor = -1.;

                //get lon/lat's
                std::vector<double> LATtmp,LONtmp,Xtmp,Ytmp;
                dataset->getLatLonXYVectors(
                        LATtmp,LONtmp,Xtmp,Ytmp,mapNextsim);

                //now get data and rotate
                double mwd,uwave,vwave,delta_r_bis,lon_factor,lat_factor;
                double R_pol = 1.e3*mapNextsim->polar_radius;//m
                double R_eq  = 1.e3*mapNextsim->equatorial_radius;//m

#if 0
                double lon_check = 4.1589999;
                double lat_check = 76.560997;
                double xcheck = 0.;
                double ycheck = 0.;
                double Rcheck_min = 1000.e3;
                double mwd_check = 1000.;
                int Icheck = 0;

                forward_mapx(mapNextsim,lat_check,lon_check,
                    &xcheck,&ycheck);
#endif

                for (int i=0; i<final_MN; ++i)
                {
                    lon_tmp = LONtmp[i];
                    lat_tmp = LATtmp[i];
                    x_tmp   = Xtmp[i];
                    y_tmp   = Ytmp[i];
                    mwd     = dataset->variables[j0].loaded_data[fstep][i];
                    if ((std::isnan(mwd))||(lat_tmp>89.5))
                    {
                        // mwd ill-defined near north pole
                        // - not usually any waves there anyway
                        tmp_data0   = 0.;
                        tmp_data1   = 0.;
                    }
                    else
                    {
#if 1
                        //use shape of earth according to nextsim projection
                        lon_factor  = R_eq*std::cos((PI/180.)*lat_tmp);
                        lat_factor  = std::hypot( R_eq*std::sin((PI/180.)*lat_tmp),
                                                 R_pol*std::cos((PI/180.)*lat_tmp));
#else
                        // use spherical earth
                        // - difference is ~1 degree in the direction
                        lon_factor  = R*std::cos((PI/180.)*lat_tmp);
                        lat_factor  = R;
#endif

                        delta_x     = delta_r*disp_factor*std::sin( (PI/180.)*mwd );
                        delta_y     = delta_r*disp_factor*std::cos( (PI/180.)*mwd );
                        delta_lon   = (180./PI)*delta_x/lon_factor;
                        delta_lat   = (180./PI)*delta_y/lat_factor;


                        //new lon,lat after moving delta_r
                        lat_tmp_bis = lat_tmp+delta_lat;
                        lon_tmp_bis = lon_tmp+delta_lon;
                        lon_tmp_bis = lon_tmp_bis-360.*std::floor(lon_tmp_bis/(360.));

                        //(x,y) position after moving delta_r
                        forward_mapx(mapNextsim,lat_tmp_bis,lon_tmp_bis,
                                &x_tmp_bis,&y_tmp_bis);
                        delta_x     = x_tmp_bis-x_tmp;
                        delta_y     = y_tmp_bis-y_tmp;
                        delta_r_bis = std::hypot(delta_x,delta_y);
                        tmp_data0   = delta_x/delta_r_bis;
                        tmp_data1   = delta_y/delta_r_bis;

#if 0
                        double Rcheck  = std::hypot(x_tmp-xcheck,y_tmp-ycheck);
                        if (Rcheck<Rcheck_min)
                        {
                            Rcheck_min  = Rcheck;
                            Icheck      = i;
                            mwd_check   = mwd;
                        }
#endif

#if 0
                        if (i==83)
                        {
                            std::cout<<std::setprecision(15)
                                <<i<<" - checking rotation of mwd:\n";
#if 0
                            std::cout<<std::setprecision(15)
                                <<i<<" - R,a,b (km) = "<<R/1.e3<<","<<mapNextsim->equatorial_radius<<","<<mapNextsim->polar_radius<<"\n";
                            std::cout<<std::setprecision(15)
                                <<i<<" - dlon,dlat,dx,dy = "<<delta_lon<<","<<delta_lat<<","<<delta_x<<","<<delta_y<<"\n";
                            std::cout<<std::setprecision(15)
                                <<i<<" - x2,y2 = "<<x_tmp_bis<<","<<y_tmp_bis<<"\n";
#endif
                            std::cout<<std::setprecision(15)
                                <<i<<" - fstep,lon,lat,x,y,mwd = "<<fstep<<","<<lon_tmp<<","<<lat_tmp<<","<<x_tmp<<","<<y_tmp<<","<<mwd<<"\n";
                            std::cout<<std::setprecision(15)
                                <<i<<" - Hs,Tp,uwave,vwave ="
                                <<" "<<dataset->variables[0].loaded_data[fstep][i] //Hs
                                <<","<<dataset->variables[1].loaded_data[fstep][i] //Tp
                                <<","<<tmp_data0
                                <<","<<tmp_data1
                                <<"\n";
                        }
#endif
                    }//rotation of unit vector representing mwd

                    // modify original scalar variables
                    // - store unit vector representing wave direction
                    // as a unit vector in nextsim's x/y coordinates
                    dataset->variables[j0].loaded_data[fstep][i] = tmp_data0;
                    dataset->variables[j1].loaded_data[fstep][i] = tmp_data1;

                }//i loop
#if 0
                lon_tmp = LONtmp[Icheck];
                lat_tmp = LATtmp[Icheck];
                x_tmp   = Xtmp[Icheck];
                y_tmp   = Ytmp[Icheck];
                mwd     = mwd_check;
                //
                std::cout<<std::setprecision(15)
                    <<Icheck<<" - checking rotation of mwd near:\n";
                std::cout<<std::setprecision(15)
                    <<"time,lon,lat,x,y = "
                    <<M_current_time<<lon_check<<","<<lat_check<<","<<xcheck<<","<<ycheck<<"\n";
                std::cout<<std::setprecision(15)
                    <<Icheck<<" - fstep,lon,lat,x,y,mwd = "<<fstep<<","<<lon_tmp<<","<<lat_tmp<<","<<x_tmp<<","<<y_tmp<<","<<mwd<<"\n";
                std::cout<<std::setprecision(15)
                    <<Icheck<<" - Hs,Tp,uwave,vwave ="
                    <<" "<<dataset->variables[0].loaded_data[fstep][Icheck] //Hs
                    <<","<<dataset->variables[1].loaded_data[fstep][Icheck] //Tp
                    <<","<<dataset->variables[j0].loaded_data[fstep][Icheck]//uwave
                    <<","<<dataset->variables[j1].loaded_data[fstep][Icheck]//vwave
                    <<"\n";
#endif

                Xtmp.resize(0);
                Ytmp.resize(0);
                LONtmp.resize(0);
                LATtmp.resize(0);
            }//mwd option

            else if(dataset->vectorial_variables[j].east_west_oriented)
            // Method generalized for regular and irregular grids
            {
                //int M=dataset->grid.gridLAT.size();
                //int N=dataset->grid.gridLON.size();
                std::vector<double> LATtmp,LONtmp,Xtmp,Ytmp;
                dataset->getLatLonXYVectors(
                        LATtmp,LONtmp,Xtmp,Ytmp,mapNextsim);

                for (int i=0; i<LATtmp.size(); ++i)
                {
                    lat_tmp=LATtmp[i];
                    lon_tmp=LONtmp[i];

                    // velocity in the east (component 0) and north direction (component 1) in m/s
                    tmp_data0=dataset->variables[j0].loaded_data[fstep][i];
                    tmp_data1=dataset->variables[j1].loaded_data[fstep][i];
  
                    // velocity in the east (component 0) and north direction (component 1) in degree/s
                    if(lat_tmp<90.)
                    {
                        tmp_data0_deg=tmp_data0/(R*std::cos(lat_tmp*PI/180.))*180./PI;
                        tmp_data1_deg=tmp_data1/R*180./PI;
                    }
                    else
                    {
                        tmp_data0_deg=0.;
                        tmp_data1_deg=0.;
                    }
  
                    // position in lat, lon after delta_t
                    lat_tmp_bis=lat_tmp+tmp_data1_deg*delta_t;
                    lon_tmp_bis=lon_tmp+tmp_data0_deg*delta_t;
  
                    lon_tmp_bis = lon_tmp_bis-360.*std::floor(lon_tmp_bis/(360.));
  
                    // initial position x, y in meter
    		        forward_mapx(mapNextsim,lat_tmp,lon_tmp,&x_tmp,&y_tmp);
  
                    // position x, y after delta_t in meter
                    forward_mapx(mapNextsim,lat_tmp_bis,lon_tmp_bis,&x_tmp_bis,&y_tmp_bis);
  
                    // velocity in m/s
                    new_tmp_data0= (x_tmp_bis-x_tmp)/delta_t;
                    new_tmp_data1= (y_tmp_bis-y_tmp)/delta_t;
  
                    // normalisation
                    speed=std::hypot(tmp_data0,tmp_data1);
                    new_speed=std::hypot(new_tmp_data0,new_tmp_data1);
  
                    if(new_speed>0.)
                    {
                        new_tmp_data0=new_tmp_data0/new_speed*speed;
                        new_tmp_data1=new_tmp_data1/new_speed*speed;
                    }
  
                    dataset->variables[j0].loaded_data[fstep][i]=new_tmp_data0;
                    dataset->variables[j1].loaded_data[fstep][i]=new_tmp_data1;
                }

                // Treat the case of the North Pole where east-north components are ill-defined
                // We just take the mean of the velocity computed at a lower latitude in the polar stereo projection
                bool found_north_pole=false;
                int y_ind, y_ind_m1;
                if(dataset->grid.gridY[0]==90.)
                {
                    y_ind=0;
                    y_ind_m1=y_ind+1;
                    found_north_pole=true;
                }
                if(dataset->grid.gridY[M-1]==90.)
                {
                    y_ind=M-1;
                    y_ind_m1=y_ind-1;
                    found_north_pole=true;
                }
                if(found_north_pole)
                {
                    tmp_data0=0.;
                    tmp_data1=0.;
                    // loop over the nearest latitude that is not 90 to compute the mean value
                    for (int x_ind=0; x_ind<N; ++x_ind)
                    {
                        int i=y_ind_m1*N+x_ind;
                        tmp_data0=tmp_data0+dataset->variables[j0].loaded_data[fstep][i];
                        tmp_data1=tmp_data1+dataset->variables[j1].loaded_data[fstep][i];
                    }
                    tmp_data0=tmp_data0/N;
                    tmp_data1=tmp_data1/N;

                    // loop over the points defined at 90N
                    for (int x_ind=0; x_ind<N; ++x_ind)
                    {
                        int i=y_ind*N+x_ind;
                        dataset->variables[j0].loaded_data[fstep][i]=tmp_data0;
                        dataset->variables[j1].loaded_data[fstep][i]=tmp_data1;
                    }
                }
            }//east-west oriented & interpolating FromGridToMesh

            if(dataset->rotation_angle!=0.)
            {
                // rotate using cos_m_diff_angle and sin_m_diff_angle
                // - if dataset & nextsim use different stereographic projections
                for (int i=0; i<final_MN; ++i)
                {
                    tmp_data0=dataset->variables[j0].loaded_data[fstep][i];
                    tmp_data1=dataset->variables[j1].loaded_data[fstep][i];

                    new_tmp_data0= cos_m_diff_angle*tmp_data0+sin_m_diff_angle*tmp_data1;
                    new_tmp_data1=-sin_m_diff_angle*tmp_data0+cos_m_diff_angle*tmp_data1;

                    dataset->variables[j0].loaded_data[fstep][i]= new_tmp_data0;
                    dataset->variables[j1].loaded_data[fstep][i]= new_tmp_data1;
                }
            }//rotation due to different projections
#ifdef OASIS
            else if(dataset->grid.gridded_rotation_angle)
            {
                double rotation_angle = mapNextsim->rotation*PI/180.;
                // rotate using the rotation angle "Theta" read from the grid file
                for (int i=0; i<final_MN; ++i)
                {
                    tmp_data0=dataset->variables[j0].loaded_data[fstep][i];
                    tmp_data1=dataset->variables[j1].loaded_data[fstep][i];

                    cos_m_diff_angle=std::cos(rotation_angle - dataset->grid.gridTheta[i]);
                    sin_m_diff_angle=std::sin(rotation_angle - dataset->grid.gridTheta[i]);

                    new_tmp_data0= cos_m_diff_angle*tmp_data0+sin_m_diff_angle*tmp_data1;
                    new_tmp_data1=-sin_m_diff_angle*tmp_data0+cos_m_diff_angle*tmp_data1;

                    dataset->variables[j0].loaded_data[fstep][i]= new_tmp_data0;
                    dataset->variables[j1].loaded_data[fstep][i]= new_tmp_data1;
                }
            } // rotation due to distorted grid
#endif
        }//loop over vectorial variables
    }

    // closing maps
    close_mapx(mapNextsim);
}//transformData

void
ExternalData::interpolateDataset(Dataset *dataset, std::vector<double> const& RX_in,
        std::vector<double> const& RY_in)
{
    // ---------------------------------
    // Spatial interpolation
    LOG(DEBUG)<<"Spatial interpolation of the data\n";

    // size of the data
    int M  = dataset->grid.dimension_y_count;
    int N  = dataset->grid.dimension_x_count;

    int MN = M*N;

    int final_MN=MN;

    int cyclic_N=N;
    int cyclic_M=M;

    if(dataset->grid.reduced_nodes_ind.size()!=0)
    {
        final_MN=dataset->grid.reduced_nodes_ind.size();
    }
    else
    {
        double delta_y=dataset->grid.gridY[M-1]-dataset->grid.gridY[M-2];
        if(dataset->grid.dimension_y.cyclic && ((dataset->grid.dimension_y_start+M)==(0+dataset->grid.dimension_y_count_netcdf)))
        {
            cyclic_M=M+1;
            dataset->grid.gridY.push_back(dataset->grid.gridY[M-1]+delta_y);
        }

        double delta_x=dataset->grid.gridX[N-1]-dataset->grid.gridX[N-2];
        if(dataset->grid.dimension_x.cyclic && ((dataset->grid.dimension_x_start+N)==(0+dataset->grid.dimension_x_count_netcdf)))
        {
            cyclic_N=N+1;
            dataset->grid.gridX.push_back(dataset->grid.gridX[N-1]+delta_x);
        }

        final_MN=cyclic_M*cyclic_N;
    }

    // Collect all the data before the interpolation
    LOG(DEBUG) << "Collect the variables before interpolation:" <<"\n";
    std::vector<double> data_in(dataset->variables.size()*dataset->nb_forcing_step*final_MN);

    for (int fstep=0; fstep < dataset->nb_forcing_step; ++fstep)
    {
        for(int j=0; j<dataset->variables.size(); ++j)
        {
            // If one of the dimension is cyclic
            if((cyclic_M!=M) || (cyclic_N!=N))
            {
                for (int y_ind=0; y_ind<M; ++y_ind)
                {
                    for (int x_ind=0; x_ind<N; ++x_ind)
                    {
                        int i=y_ind*N+x_ind;
                        int cyclic_i=y_ind*cyclic_N+x_ind;
                        int ind = (dataset->variables.size()*dataset->nb_forcing_step)*cyclic_i+fstep*dataset->variables.size()+j;
                        data_in[ind] = dataset->variables[j].loaded_data[fstep][i];
                    }
                }

                if(cyclic_M!=M)
                {
                    for (int x_ind=0; x_ind<N; ++x_ind)
                    {
                        int i = 0*N+x_ind;
                        int cyclic_i = (cyclic_M-1)*cyclic_N+x_ind;
                        int ind = (dataset->variables.size()*dataset->nb_forcing_step)*cyclic_i+fstep*dataset->variables.size()+j;
                        data_in[ind] = dataset->variables[j].loaded_data[fstep][i];
                    }
                }
                if(cyclic_N!=N)
                {
                    for (int y_ind=0; y_ind<M; ++y_ind)
                    {
                        int i = y_ind*N+0;
                        int cyclic_i = y_ind*cyclic_N+(cyclic_N-1);
                        int ind = (dataset->variables.size()*dataset->nb_forcing_step)*cyclic_i+fstep*dataset->variables.size()+j;
                        data_in[ind] = dataset->variables[j].loaded_data[fstep][i];
                    }
                }
            }
            else // with no cyclic dimension, simply use the same indice i
                for (int i=0; i<final_MN; ++i)
                {
                    int ind = (dataset->variables.size()*dataset->nb_forcing_step)*i+fstep*dataset->variables.size()+j;
                    data_in[ind] = dataset->variables[j].loaded_data[fstep][i];
                }

        }
    }

    // ---------------------------------
    // Projection of the mesh positions into the coordinate system of the data before the interpolation
    // (either the lat,lon projection or a polar stereographic projection with another rotaion angle (for ASR))
    // we should need to that also for the TOPAZ native grid, so that we could use a gridtomesh, now we use the latlon of the TOPAZ grid

    // Define the mapping and rotation_angle
	mapx_class *mapNextsim;
	std::string configfileNextsim = Environment::nextsimMppfile();
	std::vector<char> strNextsim(configfileNextsim.begin(), configfileNextsim.end());
	strNextsim.push_back('\0');
	mapNextsim = init_mapx(&strNextsim[0]);

#if 0
    auto RX = mesh.coordX(dataset->rotation_angle);
    auto RY = mesh.coordY(dataset->rotation_angle);

    if(dataset->target_size==mesh.numTriangles())
    {
    	RX = mesh.bCoordX(dataset->rotation_angle);
        RY = mesh.bCoordY(dataset->rotation_angle);
    }

	if(dataset->grid.interpolation_in_latlon)
	{
		double lat, lon;

		for (int i=0; i<dataset->target_size; ++i)
		{
			inverse_mapx(mapNextsim,RX[i],RY[i],&lat,&lon);
			RY[i]=lat;
			RX[i]=lon;
			//tmp_latlon = XY2latLon(RX[i], RY[i], map, configfile);
			//RY[i]=tmp_latlon[0];
			//RX[i]=tmp_latlon[1];
		}
	}
#endif
    std::vector<double> RX,RY;//size set in convertTargetXY;
    dataset->convertTargetXY(&(dataset->grid), RX_in, RY_in, RX, RY,mapNextsim);

    // closing maps
    close_mapx(mapNextsim);


    LOG(DEBUG) << "Interpolation:" <<"\n";

    double* data_out;

    switch(dataset->grid.interpolation_method)
    {
        case InterpolationType::FromGridToMesh:
            InterpFromGridToMeshx(  data_out, &dataset->grid.gridX[0], dataset->grid.gridX.size(), &dataset->grid.gridY[0], dataset->grid.gridY.size(),
                                  &data_in[0], dataset->grid.gridY.size(), dataset->grid.gridX.size(),
                                  dataset->variables.size()*dataset->nb_forcing_step,
                                 &RX[0], &RY[0], M_target_size, 100000000., dataset->grid.interp_type); // We put an excessively high default value, so that it will most likely crashes when not finding data
        break;
        case InterpolationType::FromMeshToMesh2dx:
            InterpFromMeshToMesh2dx(&data_out,
                                &(dataset->grid.pfindex)[0],&dataset->grid.gridX[0],&dataset->grid.gridY[0],
                                        dataset->grid.gridX.size(),dataset->grid.pfnels,
                                        &data_in[0],
                                        dataset->grid.gridX.size(),dataset->variables.size()*dataset->nb_forcing_step,
                                        &RX[0], &RY[0], M_target_size,
                                        false);
        break;
#if defined OASIS
        case InterpolationType::ConservativeRemapping:
            ConservativeRemappingGridToMesh(data_out, data_in,
                                 dataset->variables.size()*dataset->nb_forcing_step,
                                 M_target_size,
                                 dataset->M_gridP, dataset->M_triangles, dataset->M_weights);
        break;
        case InterpolationType::FromMeshToMeshQuick:
            InterpFromMeshToMesh2dx_apply(&data_out,
                                 dataset->M_areacoord, dataset->M_vertex, dataset->M_it,
                                 dataset->grid.gridX.size(),dataset->grid.pfnels,
                                 &data_in[0],
                                 dataset->grid.gridX.size(),dataset->variables.size()*dataset->nb_forcing_step,
                                 &RX[0], &RY[0]);
        break;
#endif
        default:
            LOG(DEBUG) << "invalid interpolation type:" <<"\n";
            throw std::logic_error("invalid interpolation type");
    }

    // reset grid with cyclic dimension after interpolation
    if(cyclic_M!=M)
        dataset->grid.gridY.pop_back();

    if(cyclic_N!=N)
        dataset->grid.gridX.pop_back();

    // Redistribute all the data after the interpolation
    //std::cout << "Redistribution of the interpolated variables:" <<"\n";
    for (int fstep=0; fstep < dataset->nb_forcing_step; ++fstep)
    {
        for (int j=0; j<dataset->variables.size(); ++j)
        {
            dataset->variables[j].interpolated_data[fstep].resize(M_target_size);
            for (int i=0; i<M_target_size; ++i)
                dataset->variables[j].interpolated_data[fstep][i]=data_out[(dataset->variables.size()*dataset->nb_forcing_step)*i+fstep*dataset->variables.size()+j];
        }
    }

    xDelete<double>(data_out);
    dataset->interpolated=true;

}//interpolateDataset

void
ExternalData::clear()
{
    //xDelete<double>(fcoeff);
}

} // Nextsim
