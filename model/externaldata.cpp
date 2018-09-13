/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   externaldata.cpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#include <externaldata.hpp>
#include <date.hpp>
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

ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId, bool is_vector, double StartingTime )
	:
    M_is_constant( false ),
    M_dataset( dataset ),
    M_VariableId( VariableId ),
    M_bias_correction( 0. ),
    M_is_vector( is_vector ),
    M_current_time( 0. ),
    M_StartingTime( StartingTime ),
    M_SpinUpDuration( 0. ),
    M_initialized(true)
{
    M_datasetname = (boost::format( "%1%...%2%" )
                    % M_dataset->grid.prefix
                    % M_dataset->grid.postfix
                    ).str();

    fcoeff.resize(2);
}


ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId, bool is_vector, double StartingTime, double SpinUpDuration )
	:
    ExternalData(dataset, mesh, VariableId, is_vector, StartingTime )
    {
        M_SpinUpDuration= SpinUpDuration ;
    }

ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId, double bias_correction, bool is_vector, double StartingTime )
	:
    ExternalData(dataset, mesh, VariableId, is_vector, StartingTime )
    {
        M_bias_correction= bias_correction ;
    }

ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId, double bias_correction, bool is_vector, double StartingTime, double SpinUpDuration )
	:
    ExternalData(dataset, mesh, VariableId, is_vector, StartingTime, SpinUpDuration )
    {
        M_bias_correction= bias_correction ;
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
    M_initialized(true)
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

//void ExternalData::check_and_reload(GmshMesh const& mesh, const double current_time)
void ExternalData::check_and_reload(std::vector<double> const& RX_in,
            std::vector<double> const& RY_in, const double current_time)
#ifdef OASIS
{
    this->check_and_reload(RX_in, RY_in, current_time, -1);
}

void ExternalData::check_and_reload(std::vector<double> const& RX_in,
            std::vector<double> const& RY_in, const double current_time, const int cpl_time)
#endif
{
    M_current_time = current_time;
    M_target_size = RX_in.size();

    double current_time_tmp=M_current_time;

    M_factor=1.;
    if((M_current_time-M_StartingTime)<M_SpinUpDuration)
        M_factor=(M_current_time-M_StartingTime)/M_SpinUpDuration;

    if(!M_is_constant)
    {
        bool to_be_reloaded=false;

        if(M_dataset->grid.dataset_frequency=="constant")
            to_be_reloaded=!M_dataset->loaded;
        else if(M_dataset->grid.dataset_frequency=="nearest_daily")
            to_be_reloaded=(to_date_string_yd(current_time)!=to_date_string_yd(M_dataset->ftime_range[0]) || !M_dataset->loaded);
        else
            to_be_reloaded=((current_time_tmp < M_dataset->ftime_range[0]) || (M_dataset->ftime_range[1] < current_time_tmp) || !M_dataset->loaded);

#ifdef OASIS
        // We call oasis_get every time step, but only actually recieve data at coupling times
        if (M_dataset->coupled)
            this->recieveCouplingData(M_dataset, cpl_time);
#endif

        if (to_be_reloaded)
        {
            std::cout << "Load " << M_datasetname << "\n";
            this->loadDataset(M_dataset, RX_in, RY_in);
            this->transformData(M_dataset);
            std::cout << "Done\n";

            //need to interpolate again if reloading
            M_dataset->interpolated = false;
        }

        if (!M_dataset->interpolated)
        {
            std::cout << "Interpolate " << M_datasetname << "\n";
            this->interpolateDataset(M_dataset, RX_in, RY_in);
            std::cout << "Done\n";
        }
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
            ASSERT(i < M_target_size, "invalid index");
            value = M_constant_value;
        }
        else
        {
            ASSERT(i < 2*M_target_size, "invalid index");

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
                && M_dataset->grid.dataset_frequency!="nearest_daily");
        bool interp_const_wave  = ((M_dataset->grid.waveOptions.wave_dataset)
                && (Environment::vm()["wimsetup.wave-time-interp-option"].as<std::string>()=="step"));
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
                ASSERT(i < M_target_size, "invalid index");
                i_tmp=i;
                VariableId_tmp=M_VariableId;
            }
            else
            {
                ASSERT(i < 2*M_target_size, "invalid index");

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
                ASSERT(i < M_target_size, "invalid index");
                i_tmp=i;
                VariableId_tmp=M_VariableId;
            }
            else
            {
                ASSERT(i < 2*M_target_size, "invalid index");

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
ExternalData::recieveCouplingData(Dataset *dataset, int cpl_time)
{
        // ierror = OASIS3::get_2d(var_id[1], pcpt*time_step, &field2_recv[0], M_cpl_out.M_ncols, M_cpl_out.M_nrows);
        for(int j=0; j<dataset->vectorial_variables.size(); ++j)
            int ierror = OASIS3::get_2d(M_VariableId, cpl_time, &dataset->variables[j].loaded_data[0][0], dataset->grid.dimension_x_count, dataset->grid.dimension_y_count);
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
	std::string configfileNextsim = (boost::format( "%1%/%2%/%3%" )
                              % Environment::nextsimDir().string()
                              % "data"
                              % Environment::vm()["mesh.mppfile"].as<std::string>()
                              ).str();

	std::vector<char> strNextsim(configfileNextsim.begin(), configfileNextsim.end());
	strNextsim.push_back('\0');
	mapNextsim = init_mapx(&strNextsim[0]);

    if(dataset->grid.mpp_file!="")
    {
        mapx_class *map;
	    std::string configfile = (boost::format( "%1%/%2%/%3%" )
                              % Environment::nextsimDir().string()
                              % "data"
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
    // Projection of the mesh positions into the coordinate system of the data before the interpolation
    // (either the lat,lon projection or a polar stereographic projection with another rotation angle (for ASR))
    // we should need to that also for the TOPAZ native grid, so that we could use a gridtomesh, now we use the latlon of the TOPAZ grid

    std::vector<double> RX,RY;//size set in convertTargetXY
    this->convertTargetXY(dataset,RX_in,RY_in,RX,RY,mapNextsim);

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

    // closing maps
    close_mapx(mapNextsim);

    double RX_min=*std::min_element(RX.begin(),RX.end());
    double RX_max=*std::max_element(RX.begin(),RX.end());
    double RY_min=*std::min_element(RY.begin(),RY.end());
    double RY_max=*std::max_element(RY.begin(),RY.end());

    // ---------------------------------
    // Load grid if unloaded
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
                    init_time = Nextsim::from_date_time_string(tmpstr);
            }
            else
            {
                // use forecast.time_init_ocean_fc option to get init_time
                std::string tmpstr = (Environment::vm()["forecast.time_init_ocean_fc"].as<std::string>());
                if(tmpstr!="")
                    init_time = Nextsim::from_date_time_string(tmpstr);
            }
        }
        //only need init_time to get grid
        dataset->loadGrid(&(dataset->grid), init_time, init_time, RX_min, RX_max, RY_min, RY_max);
    }

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
    netCDF::NcDim tmpDim;

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
        bool is_topaz_fc = (dataset->grid.dataset_frequency=="daily_forecast");//topaz forecast
        bool is_ec_dataset = ((dataset->grid.prefix).find("start") != std::string::npos);//ec_[nodes,elements],ec2_[nodes,elements]
        bool true_forecast = (Environment::vm()["forecast.true_forecast"].as<bool>());

        ftime = M_current_time-dataset->averaging_period/2.;
        file_jump ={-1,0,1};

        std::cout<<"LOAD DATASET TIMES:\n";
        std::cout<<"init_time = "<<init_time<<" = "<<to_date_time_string(init_time)<<"\n";
        std::cout<<"M_current_time = "<<M_current_time<<" = "<<to_date_time_string(M_current_time)<<"\n";
        std::cout<<"ftime = "<<ftime<<" = "<<to_date_time_string(ftime)<<"\n";
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
                    init_time = Nextsim::from_date_time_string(tmpstr);
            }
            else
            {
                std::string tmpstr = (Environment::vm()["forecast.time_init_ocean_fc"].as<std::string>());
                if(tmpstr!="")
                    init_time = Nextsim::from_date_time_string(tmpstr);
            }
        }//forecasts

        for (auto jump_ptr = file_jump.begin() ; jump_ptr != file_jump.end(); ++jump_ptr)
        {
            int jump = *jump_ptr;//get jump as an integer
            std::cout<<"jump = "<<jump<<"\n";
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

            std::cout<<"FILENAME (JUMPS) = "<< filename <<"\n";
            if ( ! boost::filesystem::exists(filename) )
            {
                std::cout << "File not found: " + filename<< "; jump = "<<jump << "\n";
                continue;
                //throw std::runtime_error("File not found: " + filename);
            }
            else
                std::cout << "File found: " + filename<< "; jump = "<<jump << "\n";

            index_start.resize(1);
            index_count.resize(1);

            bool has_time_variable = true;

            try // we try because sometimes no time dimension is available in the netcdf
            {
                // Open the netcdf file
                netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

                // Set the time range XTIME
                netCDF::NcVar FVTIME = dataFile.getVar(dataset->time.name);
		        netCDF::NcDim timeDim = dataFile.getDim(dataset->time.name);
                index_start[0]=0;
                index_count[0]=timeDim.getSize();
                XTIME.resize(index_count[0]);
                FVTIME.getVar(index_start, index_count, &XTIME[0]);
            }
            catch(const std::exception& e) // if no time dimension is available in the netcdf, we define the time as
            {
                has_time_variable=false;
                index_start[0]=0;
                index_count[0]=1;
                XTIME.resize(index_count[0]);
                XTIME[0]=-1.;

                if(dataset->grid.dataset_frequency=="monthly" || dataset->grid.dataset_frequency=="yearly")
                    throw std::runtime_error("The case monthly and yearly when no time dimension is available is not implemented!");
            }

            double f;
            int nt = XTIME.size();
            std::cout<<"582: "<<nt<<"\n";
            if(is_ec_fc && (!true_forecast))
                nt = 4;// just use the first day of each file (1st 4 records, each 6 hours apart)

            for (int it=0; it < nt; ++it) // always need one step before and one after the target time
            {
                 if (!has_time_variable || ((dataset->name).find("ice_amsr2") != std::string::npos))
                     f = from_date_string((boost::format( "%1%-%2%-%3%" )
                                 % f_timestr.substr(0,4)
                                 % f_timestr.substr(4,2)
                                 % f_timestr.substr(6,2)).str())+0.5;
                 else
                     f = (XTIME[it]*dataset->time.a+dataset->time.b)/24.0
                          + from_date_string(dataset->grid.reference_date);

                 if(f>M_current_time && index_next==-1)
                 {
                     time_next=f;
                     index_next = it;
                     filename_next = filename;
                 }
                 if(f<=M_current_time)
                 { 
                     time_prev=f;
                     index_prev = it;
                     filename_prev = filename;
                 }
                 std::cout<<"608"<< f << ", " << M_current_time << ", " << filename_prev <<", " << filename_next <<"\n";
                 std::cout<<"609"<< time_prev << ", " << index_prev << ", " << time_next <<", " << index_next <<"\n";
            }
        }//loop over jump

        if(filename_prev!="")
        {
            filename_fstep.push_back(filename_prev);
            index_fstep.push_back(index_prev);
        }
        else
            throw std::runtime_error("Not found a file for before current_date!");


        if(filename_next!="")
        {
            filename_fstep.push_back(filename_next);
            index_fstep.push_back(index_next);
        }
        else
            throw std::runtime_error("Not found a file for after current_date!");

		dataset->ftime_range = {time_prev,time_next};
	}//not nearest_daily or constant
    else
    {
        if(dataset->grid.dataset_frequency=="nearest_daily")
        {
            ftime = M_current_time;
            f_timestr = to_date_string_yd(std::floor(ftime));
            double f=from_date_string((boost::format( "%1%-%2%-%3%" )
                        % f_timestr.substr(0,4)
                        % f_timestr.substr(4,2)
                        % f_timestr.substr(6,2)).str());
            dataset->ftime_range = {f+.5};
        }
        else
            f_timestr ="";
            filename = (boost::format( "%1%/%2%/%3%%4%%5%" )
                        % Environment::simdataDir().string()
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

    std::cout<<"Start loading data\n";
    for (int fstep=0; fstep < nb_forcing_step; ++fstep) // always need one step before and one after the target time
    {
        filename=filename_fstep[fstep];
        index=index_fstep[fstep];

        std::cout<<"FILENAME= "<< filename <<"\n";
        if ( ! boost::filesystem::exists(filename) )
            throw std::runtime_error("File not found: " + filename);

        // Open the netcdf file
        netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

        // Load each variable and copy its data into loaded_data
        for(int j=0; j<dataset->variables.size(); ++j)
        {
            std::cout<<"variables number:" << j  << "\n";
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
            index_count.resize(dataset->variables[j].dimensions.size());

            // here we find the start and count index for each dimensions
            for(int k=0; k<dataset->variables[j].dimensions.size(); ++k)
            {
                std::string dimension_name=dataset->variables[j].dimensions[k].name;

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
                // other cases
                else{
                    tmpDim = dataFile.getDim(dimension_name);

                    index_start[k] = 0;
                    index_count[k] = tmpDim.getSize();
                }
            }

            // time dimension
			if(dataset->variables[j].dimensions.size()>2
                    && dataset->grid.dataset_frequency!="constant"
                    && dataset->grid.dataset_frequency!="nearest_daily")
			{
            	index_start[0] = index;
            	index_count[0] = 1;
			}

            // depth dimension
			if(dataset->variables[j].dimensions.size()>3)
			{
            	index_start[1] = 0;
            	index_count[1] = 1;
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
}


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
    std::string configfileNextsim = (boost::format( "%1%/%2%/%3%" )
                              % Environment::nextsimDir().string()
                              % "data"
                              % Environment::vm()["mesh.mppfile"].as<std::string>()
                              ).str();

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


    std::cout<<"Start transforming the data\n";
    for (int fstep=0; fstep < dataset->nb_forcing_step; ++fstep) // always need one step before and one after the target time
    {
        for(int j=0; j<dataset->vectorial_variables.size(); ++j)
        {

            j0=dataset->vectorial_variables[j].components_Id[0];
            j1=dataset->vectorial_variables[j].components_Id[1];

            if (dataset->variables[j0].wavDirOptions.isWavDir)
            {
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

            else if(dataset->vectorial_variables[j].east_west_oriented
                    && dataset->grid.interpolation_method==InterpolationType::FromGridToMesh)
            {
                int M=dataset->grid.gridLAT.size();
                int N=dataset->grid.gridLON.size();

                for (int y_ind=0; y_ind<M; ++y_ind)
                {
                    for (int x_ind=0; x_ind<N; ++x_ind)
                    {
                        int i=y_ind*N+x_ind;

                        // lat lon of the data
                        lat_tmp=dataset->grid.gridLAT[y_ind];
                        lon_tmp=dataset->grid.gridLON[x_ind];

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
        }//loop over vectorial variables
    }

    // closing maps
    close_mapx(mapNextsim);
}


void
ExternalData::convertTargetXY(Dataset *dataset,
        std::vector<double> const& RX_in,  std::vector<double> const& RY_in,
        std::vector<double> & RX_out, std::vector<double> & RY_out,
        mapx_class *mapNextsim)//(double const& u, double const& v)
{
    RX_out.resize(M_target_size);
    RY_out.resize(M_target_size);

    if(dataset->grid.interpolation_in_latlon)
    {
        double lat, lon;
        for (int i=0; i<M_target_size; ++i)
        {
            //convert to lon,lat
			inverse_mapx(mapNextsim,RX_in[i],RY_in[i],&lat,&lon);
			RY_out[i]=lat;
            double bc_lon=dataset->grid.branch_cut_lon;
            bool close_on_right=false;
                //if true  make target lon >  bc_lon,<=bc_lon+180
                //if false make target lon >= bc_lon,< bc_lon+180
                //this shouldn't matter here though?
			RX_out[i]=dataset->thetaInRange(lon,bc_lon,close_on_right);
        }
    }
    else
    {
        double cos_rotangle = std::cos(dataset->rotation_angle);
        double sin_rotangle = std::sin(dataset->rotation_angle);
        //rotate to coord sys of dataset
        for (int i=0; i<M_target_size; ++i)
        {
			RX_out[i] =  cos_rotangle*RX_in[i]+sin_rotangle*RY_in[i];
			RY_out[i] = -sin_rotangle*RX_in[i]+cos_rotangle*RY_in[i];
        }
    }
}

void
ExternalData::interpolateDataset(Dataset *dataset, std::vector<double> const& RX_in,
        std::vector<double> const& RY_in)
{
    // ---------------------------------
    // Spatial interpolation
    std::cout<<"Spatial interpolation of the data\n";

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
    std::cout << "Collect the ivariables before interpolation:" <<"\n";
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
	std::string configfileNextsim = (boost::format( "%1%/%2%/%3%" )
                              % Environment::nextsimDir().string()
                              % "data"
                              % Environment::vm()["mesh.mppfile"].as<std::string>()
                              ).str();

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
    this->convertTargetXY(dataset,RX_in, RY_in, RX, RY,mapNextsim);

    // closing maps
    close_mapx(mapNextsim);


    std::cout << "Interpolation:" <<"\n";

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
        default:
            std::cout << "invalid interpolation type:" <<"\n";
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

}

void
ExternalData::clear()
{
    //xDelete<double>(fcoeff);
}

} // Nextsim
