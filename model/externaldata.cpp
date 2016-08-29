/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   externaldata.cpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#include <externaldata.hpp>
#include <date.hpp>
#include <dataset.hpp>
extern "C"
{
#include <mapx.h>
}
     

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

ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId, bool is_vector )
	:
    M_is_constant( false ),
    M_dataset( dataset ),
    M_VariableId( VariableId ),
    M_is_vector( is_vector ),
    M_current_time( 0. ),
    M_SpinUpStartingTime( 0. ),
    M_SpinUpDuration( 0. ),
    M_initialized(true)
{
    M_datasetname = (boost::format( "%1%...%2%" )
                    % M_dataset->prefix
                    % M_dataset->postfix
                    ).str();

    fcoeff.resize(2);
}


ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId, bool is_vector, double SpinUpStartingTime, double SpinUpDuration )
	:
    ExternalData(dataset, mesh, VariableId, is_vector )
    {
        M_SpinUpStartingTime= SpinUpStartingTime ;
        M_SpinUpDuration= SpinUpDuration ;
    }

ExternalData::ExternalData( double ConstantValue )
	:
    M_is_constant( true ),
    M_constant_value( ConstantValue ),
    M_is_vector( false ),
    M_current_time( 0. ),
    M_SpinUpStartingTime( 0. ),
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

ExternalData::ExternalData( double ConstantValue, double SpinUpStartingTime, double SpinUpDuration )
    :
    ExternalData( ConstantValue )
    {
        M_SpinUpStartingTime= SpinUpStartingTime ;
        M_SpinUpDuration= SpinUpDuration ;
    }

ExternalData::ExternalData( double ConstantValue, double ConstantValuebis, double SpinUpStartingTime, double SpinUpDuration )
    :
    ExternalData( ConstantValue, ConstantValuebis )
    {
        M_SpinUpStartingTime= SpinUpStartingTime ;
        M_SpinUpDuration= SpinUpDuration ;
    }

ExternalData::~ExternalData()
{
	this->clear();
}

void ExternalData::check_and_reload(GmshMesh const& mesh, const double current_time)
{
    M_current_time = current_time;
    
    double current_time_tmp=M_current_time;

    M_factor=1.;
    if((M_current_time-M_SpinUpStartingTime)<M_SpinUpDuration)
    {
        M_factor=(M_current_time-M_SpinUpStartingTime)/M_SpinUpDuration;
    }

    if(!M_is_constant)
    {
        bool to_be_reloaded=false;

        if(M_dataset->nb_timestep_day>0)
        {
            to_be_reloaded=((current_time_tmp < M_dataset->ftime_range[0]) || (M_dataset->ftime_range[1] < current_time_tmp) || !M_dataset->reloaded);
        }
        else
            to_be_reloaded=!M_dataset->reloaded;

        if (to_be_reloaded)
        {
            std::cout << "Load " << M_datasetname << "\n";
            loadDataset(M_dataset, mesh);
            std::cout << "Done\n";
        }
    }
}

typename ExternalData::value_type
ExternalData::operator [] (const size_type i)
{
    value_type value;
    size_type i_tmp;
    int VariableId_tmp;
    if(M_is_constant)
    {
        // for the moment same value is given to all the components
        value = M_constant_value;
    }
    else
    {
        if(M_dataset->nb_timestep_day>0)
        {
            fdt = std::abs(M_dataset->ftime_range[1]-M_dataset->ftime_range[0]);
            fcoeff[0] = std::abs(M_current_time-M_dataset->ftime_range[1])/fdt;
            fcoeff[1] = std::abs(M_current_time-M_dataset->ftime_range[0])/fdt;

            //std::cout <<"LINEAR COEFF 1= "<< fcoeff[0] <<"\n";
            //std::cout <<"LINEAR COEFF 2= "<< fcoeff[1] <<"\n";

            if(!M_is_vector)
            {
                ASSERT(i < M_dataset->target_size, "invalid index");
                i_tmp=i;
                VariableId_tmp=M_VariableId;
            }
            else
            {
                ASSERT(i < 2*M_dataset->target_size, "invalid index");

                if(i < M_dataset->target_size)
                {
                    i_tmp=i;
                    VariableId_tmp=M_dataset->vectorial_variables[M_VariableId].components_Id[0];
                }
                else
                {
                    i_tmp=i-M_dataset->target_size;
                    VariableId_tmp=M_dataset->vectorial_variables[M_VariableId].components_Id[1];
                }
            }
            value =  M_factor*
                (fcoeff[0]*M_dataset->variables[VariableId_tmp].data2[0][i_tmp] +
                 fcoeff[1]*M_dataset->variables[VariableId_tmp].data2[1][i_tmp]);
        }
        else
        {
            if(!M_is_vector)
            {
                ASSERT(i < M_dataset->target_size, "invalid index");
                i_tmp=i;
                VariableId_tmp=M_VariableId;
            }
            else
            {
                ASSERT(i < 2*M_dataset->target_size, "invalid index");

                if(i < M_dataset->target_size)
                {
                    i_tmp=i;
                    VariableId_tmp=M_dataset->vectorial_variables[M_VariableId].components_Id[0];

                }
                else
                {
                    i_tmp=i-M_dataset->target_size;
                    VariableId_tmp=M_dataset->vectorial_variables[M_VariableId].components_Id[1];
                }
            }
            value =  M_factor*M_dataset->variables[VariableId_tmp].data2[0][i_tmp];
        }
    }

	return static_cast<value_type>( value );
}

void
ExternalData::loadDataset(Dataset *dataset, GmshMesh const& mesh)//(double const& u, double const& v)
{
    // ---------------------------------
    // Define the mapping and rotation_angle
	mapx_class *mapNextsim;
	std::string configfileNextsim = (boost::format( "%1%/%2%/%3%" )
                              % Environment::nextsimDir().string()
                              % "data"
                              % "NpsNextsim.mpp"
                              ).str();

	std::vector<char> strNextsim(configfileNextsim.begin(), configfileNextsim.end());
	strNextsim.push_back('\0');
	mapNextsim = init_mapx(&strNextsim[0]);

	mapx_class *map;
    double rotation_angle, cos_m_diff_angle, sin_m_diff_angle;
    if(dataset->grid.mpp_file!="")
    {
	    std::string configfile = (boost::format( "%1%/%2%/%3%" )
                              % Environment::nextsimDir().string()
                              % dataset->grid.dirname
                              % dataset->grid.mpp_file
                              ).str();

	    std::vector<char> str(configfile.begin(), configfile.end());
	    str.push_back('\0');
	    map = init_mapx(&str[0]);
        rotation_angle = -(mapNextsim->rotation-map->rotation)*PI/180.;
        
        close_mapx(map);
    }
    else
    {
        rotation_angle=0.;
    }

    cos_m_diff_angle=std::cos(-rotation_angle);
    sin_m_diff_angle=std::sin(-rotation_angle);
    
    // ---------------------------------
    // Projection of the mesh positions into the coordinate system of the data before the interpolation
    // (either the lat,lon projection or a polar stereographic projection with another rotaion angle (for ASR))
    // we should need to that also for the TOPAZ native grid, so that we could use a gridtomesh, now we use the latlon of the TOPAZ grid

    auto RX = mesh.coordX(rotation_angle);
    auto RY = mesh.coordY(rotation_angle);

    if(dataset->target_size==mesh.numTriangles())
    {
    	RX = mesh.bcoordX(rotation_angle);
        RY = mesh.bcoordY(rotation_angle);
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
    
    double RX_min=*std::min_element(RX.begin(),RX.end());
    double RX_max=*std::max_element(RX.begin(),RX.end());
    double RY_min=*std::min_element(RY.begin(),RY.end());
    double RY_max=*std::max_element(RY.begin(),RY.end());

    // ---------------------------------
    // Load grid if unloaded
    if(!dataset->grid.loaded)
        dataset->loadGrid(&(dataset->grid), M_current_time, RX_min, RX_max, RY_min, RY_max);

    // ---------------------------------
    // Initialise counters etc.
	int nb_forcing_step =1;

	std::vector<double> XTIME(1);
	std::vector<size_t> index_start(1);
	std::vector<size_t> index_count(1);

	int index = 0;

	// interp_type for grid to mesh interpolation
	int interp_type = dataset->grid.interp_type;

    // create dataset->ftime_range for data sets which need to be interpolated in time
	if(dataset->nb_timestep_day>0)
	{
		double file_dt = 1./dataset->nb_timestep_day;
		
        double time_start, time_end;
        if(dataset->daily_mean)
        {
            time_start = std::floor((M_current_time-0.5)*dataset->nb_timestep_day)/dataset->nb_timestep_day+0.5;
    		time_end   = std::ceil ((M_current_time-0.5)*dataset->nb_timestep_day)/dataset->nb_timestep_day+0.5;
            std::cout << "time_start " << time_start << " " << to_date_string_yd(std::floor(time_start-0.5)) <<  "\n";
            std::cout << "time_end " << time_end     <<  " " << to_date_string_yd(std::floor(time_end-0.5)) <<"\n";
        }
        else
        {
            time_start = std::floor(M_current_time*dataset->nb_timestep_day)/dataset->nb_timestep_day;
    		time_end   = std::ceil (M_current_time*dataset->nb_timestep_day)/dataset->nb_timestep_day;
        }
            
		// We always need at least two time steps to interpolate between
		if (time_end == time_start)
		{
			time_end = time_start + (1./dataset->nb_timestep_day);
		}

		dataset->ftime_range.resize(0);
		for (double time_tmp=time_start; time_tmp<=time_end; time_tmp+=file_dt)
		{
			dataset->ftime_range.push_back(time_tmp);
		}
        
        std::cout << "dataset->ftime_range.size() " << dataset->ftime_range.size() << "\n";

		// for (int i=0; i<dataset->ftime_range.size(); ++i)
		// {
		// 	std::cout <<"TIMEVEC["<< i <<"]= "<< dataset->ftime_range[i] <<"\n";
		// }

		nb_forcing_step = dataset->ftime_range.size();
	}

    // Initialise variables for the fields
    std::vector<double> tmp_interpolated_field(dataset->target_size);

    int N_data =dataset->variables.size();
    int M  =dataset->grid.dimension_x_count;
    int N  = dataset->grid.dimension_y_count;

    int MN = M*N;

    int cyclic_N=N;
    int cyclic_M=M;
    double delta_y=dataset->grid.gridY[M-1]-dataset->grid.gridY[M-2];
    if(dataset->grid.dimension_y.cyclic)
    {
        cyclic_M=M+1;
        dataset->grid.gridY.push_back(dataset->grid.gridY[M-1]+delta_y);
    }

    double delta_x=dataset->grid.gridX[N-1]-dataset->grid.gridX[N-2];
    if(dataset->grid.dimension_x.cyclic)
    {
        cyclic_N=N+1;
        dataset->grid.gridX.push_back(dataset->grid.gridX[N-1]+delta_x);
    }

    int final_MN=cyclic_M*cyclic_N;

	if(dataset->grid.reduced_nodes_ind.size()!=0)
    {
        if((dataset->grid.dimension_y.cyclic) || (dataset->grid.dimension_x.cyclic))
            throw std::runtime_error("Using reduced grid and cyclic grid at the same time is not yet implemented");

    	final_MN=dataset->grid.reduced_nodes_ind.size();
    }

	// Memory leak:
    //double* data_in = new double[N_data*nb_forcing_step*final_MN];
    std::vector<double> data_in(N_data*nb_forcing_step*final_MN);

    std::vector<double> data_in_tmp(MN);
    //std::cout <<" \n";

    // Attributes (scaling and offset)
    netCDF::NcVarAtt att;
    double scale_factor;
    double add_offset;

    //std::cout <<"NB_FORCING_STEP= "<< nb_forcing_step <<"\n";

    // Read in data one time step at a time
    for (int fstep=0; fstep < nb_forcing_step; ++fstep)
    {
        // Define variables for this scope
        std::vector<netCDF::NcVar> NcVars(dataset->variables.size());
        double ftime;
        std::string filename;

        // Filename depends on the date for time varying data
		if(dataset->nb_timestep_day>0)
		{
            ftime = dataset->ftime_range[fstep];
         
            if(dataset->daily_mean)
                ftime = ftime-0.5;
            
            std::string f_timestr;
            if(dataset->grid.monthly_dataset)
                f_timestr = to_date_string_ym(std::floor(ftime));
            else
                f_timestr = to_date_string_yd(std::floor(ftime));
            
            std::cout <<"F_TIMESTR= "<< f_timestr <<"\n";

            filename = (boost::format( "%1%/%2%/%3%%4%%5%" )
                        % Environment::simdataDir().string()
                        % dataset->dirname
                        % dataset->prefix
                        % f_timestr
                        % dataset->postfix
                        ).str();
		}
        else
        {
            filename = (boost::format( "%1%/%2%/%3%%4%" )
                        % Environment::simdataDir().string()
                        % dataset->dirname
                        % dataset->prefix
                        % dataset->postfix
                        ).str();
        }

        std::cout<<"FILENAME= "<< filename <<"\n";
        if ( ! boost::filesystem::exists(filename) )
            throw std::runtime_error("File not found: " + filename);

        // change the reference_date if erai forcing according to the xxxx-01-01, where xxxx is the current year
        if ((dataset->name).find("ERAi") != std::string::npos)
        {
            dataset->reference_date = to_date_string_y(std::floor(ftime)) + "-01-01";
            //std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@DETECT ERAi: Dataset->reference_date= "<< dataset->reference_date <<"\n";
        }

        // Open the netcdf file
        netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

        netCDF::NcDim tmpDim;

        // Find the right time slice
        if (dataset->nb_timestep_day>1)
        {
            // Set the time range XTIME
            netCDF::NcVar FVTIME = dataFile.getVar(dataset->time.name);

            index_start.resize(1);
            index_count.resize(1);

            netCDF::NcDim timeDim = dataFile.getDim(dataset->time.name);

            index_start[0]=0;
            index_count[0]=timeDim.getSize();

            XTIME.resize(index_count[0]);

            FVTIME.getVar(index_start, index_count, &XTIME[0]);
            std::for_each(XTIME.begin(), XTIME.end(), [&](double& f){ f = f/24.0+from_date_string(dataset->reference_date); });

            auto it = std::find(XTIME.begin(), XTIME.end(), ftime);
            index = std::distance(XTIME.begin(),it);
            //std::cout <<"FIND "<< ftime <<" in index "<< index <<"\n";
        }

        // Load each variable and copy its data into data_in
        for(int j=0; j<dataset->variables.size(); ++j)
        {
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
			if(dataset->nb_timestep_day>0)
			{
            	index_start[0] = index;
            	index_count[0] = 1;
			}

            NcVars[j].getVar(index_start,index_count,&data_in_tmp[0]);

            // Need to multiply with scale factor and add offset - these are stored as variable attributes
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

            _printf_("For " << dataset->variables[j].name << " scale_factor is "  << scale_factor<<  " " <<  ", add_offset is " << add_offset << "\n");
            // Copy the data in data_in

            // If reduced_nodes is used
			if(dataset->grid.reduced_nodes_ind.size()!=0)
			{
            	for (int i=0; i<(final_MN); ++i)
                {
                    data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j]=
                        data_in_tmp[dataset->grid.reduced_nodes_ind[i]]*scale_factor + add_offset;
                    if(std::isnan(data_in_tmp[dataset->grid.reduced_nodes_ind[i]]*scale_factor + add_offset))
                    {
            			_printf_("found NaN at"  << data_in_tmp[dataset->grid.reduced_nodes_ind[i]]<<  " "<<  dataset->grid.reduced_nodes_ind[i] <<  ", default_value is used\n");
    				}
                }
			}
			else // if not reduced_node
            {
                // If one of the dimension is cyclic
                if((dataset->grid.dimension_y.cyclic) || (dataset->grid.dimension_x.cyclic))
                {
                    for (int y_ind=0; y_ind<M; ++y_ind)
                    {
                        for (int x_ind=0; x_ind<N; ++x_ind)
                        {
                            int i=y_ind*N+x_ind;
                            int cyclic_i=y_ind*cyclic_N+x_ind;

                            data_in[(dataset->variables.size()*nb_forcing_step)*cyclic_i+fstep*dataset->variables.size()+j]=data_in_tmp[i]*scale_factor + add_offset;
                        }
                    }
                }
                else // with no cyclic dimension, simply use the same indice i
                {
                    for (int i=0; i<(MN); ++i)
                    {
                        data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j]=data_in_tmp[i]*scale_factor + add_offset;
                        if(std::isnan(data_in_tmp[i]*scale_factor + add_offset))
                        {
                			_printf_("found NaN at"  << data_in_tmp[i] <<  " "<<  i <<  ", default_value is used\n");
        				}
                    }
                }
            }

            if(dataset->grid.dimension_y.cyclic)
                for (int x_ind=0; x_ind<N; ++x_ind)
                {
                    int i=0*N+x_ind;
                    int cyclic_i=(cyclic_M-1)*cyclic_N+x_ind;

                    data_in[(dataset->variables.size()*nb_forcing_step)*cyclic_i+fstep*dataset->variables.size()+j]=data_in_tmp[i]*scale_factor + add_offset;
                }

            if(dataset->grid.dimension_x.cyclic)
                for (int y_ind=0; y_ind<M; ++y_ind)
                {
                    int i=y_ind*N+0;
                    int cyclic_i=y_ind*cyclic_N+(cyclic_N-1);

                    data_in[(dataset->variables.size()*nb_forcing_step)*cyclic_i+fstep*dataset->variables.size()+j]=data_in_tmp[i]*scale_factor + add_offset;
                }

		}
    }

    // ---------------------------------
    // Transformation of the vectorial variables from the coordinate system of the data to the polar stereographic projection used in the model
    // Once we are in the polar stereographic projection, we can do spatial interpolation without bothering about the North Pole

    double tmp_data0, tmp_data1, new_tmp_data0, new_tmp_data1;
    double tmp_data0_deg, tmp_data1_deg;
    double lat_tmp, lon_tmp, lat_tmp_bis, lon_tmp_bis;
    double x_tmp, y_tmp, x_tmp_bis, y_tmp_bis;
    double speed, new_speed;
    int j0, j1;

    double R=mapx_Re_km*1000.; // Earth radius
    double delta_t=1.; // 1 sec. This value needs to be small.

    for (int fstep=0; fstep < nb_forcing_step; ++fstep)
    {
        for(int j=0; j<dataset->vectorial_variables.size(); ++j)
        {
            j0=dataset->vectorial_variables[j].components_Id[0];
            j1=dataset->vectorial_variables[j].components_Id[1];

            if(dataset->vectorial_variables[j].east_west_oriented && dataset->grid.interpolation_method==InterpolationType::FromGridToMesh)
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
                        tmp_data0=data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j0];
                        tmp_data1=data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j1];

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

                        data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j0]= new_tmp_data0;
                        data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j1]= new_tmp_data1;
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
                if(found_north_pole=true)
                {
                    tmp_data0=0.;
                    tmp_data1=0.;
                    // loop over the nearest latitude that is not 90 to compute the mean value
                    for (int x_ind=0; x_ind<N; ++x_ind)
                    {
                        int i=y_ind_m1*N+x_ind;
                        tmp_data0=tmp_data0+data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j0];
                        tmp_data1=tmp_data1+data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j1];
                    }
                    tmp_data0=tmp_data0/N;
                    tmp_data1=tmp_data1/N;

                    // loop over the points defined at 90N
                    for (int x_ind=0; x_ind<N; ++x_ind)
                    {
                        int i=y_ind*N+x_ind;
                        data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j0]=tmp_data0;
                        data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j1]=tmp_data1;
                    }
                }
            }

            if(rotation_angle!=0.)
            {
                for (int i=0; i<final_MN; ++i)
                {
                    tmp_data0=data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j0];
                    tmp_data1=data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j1];

                    new_tmp_data0= cos_m_diff_angle*tmp_data0+sin_m_diff_angle*tmp_data1;
                    new_tmp_data1=-sin_m_diff_angle*tmp_data0+cos_m_diff_angle*tmp_data1;

                    data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j0]= new_tmp_data0;
                    data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j1]= new_tmp_data1;
                }
            }
        }
    }

    // closing maps
    close_mapx(mapNextsim);

    // ---------------------------------
    // Spatial interpolation

    //std::cout <<"before interp " <<"\n";

    double* data_out;
    double tmp_data;

    switch(dataset->grid.interpolation_method)
    {
        case InterpolationType::FromGridToMesh:
            InterpFromGridToMeshx(  data_out, &dataset->grid.gridX[0], dataset->grid.gridX.size(), &dataset->grid.gridY[0], dataset->grid.gridY.size(),
                                  &data_in[0], dataset->grid.gridY.size(), dataset->grid.gridX.size(),
                                  dataset->variables.size()*nb_forcing_step,
                                 &RX[0], &RY[0], dataset->target_size, 100000000., interp_type); // We put an excessively high default value, so that it will most likely crashes when not finding data
        break;
        case InterpolationType::FromMeshToMesh2dx:
            InterpFromMeshToMesh2dx(&data_out,
                                dataset->grid.pfindex,&dataset->grid.gridX[0],&dataset->grid.gridY[0],
                                        dataset->grid.gridX.size(),dataset->grid.pfnels,
                                        &data_in[0],
                                        dataset->grid.gridX.size(),N_data*nb_forcing_step,
                                        &RX[0], &RY[0], dataset->target_size,
                                        false);
        break;
        default:
            std::cout << "invalid interpolation type:" <<"\n";
            throw std::logic_error("invalid interpolation type");
    }

    if(dataset->grid.dimension_y.cyclic)
        dataset->grid.gridY.pop_back();

    if(dataset->grid.dimension_x.cyclic)
        dataset->grid.gridX.pop_back();

    //std::cout <<"after interp " <<"\n";

    for (int fstep=0; fstep < nb_forcing_step; ++fstep)
    {
        for(int j=0; j<dataset->variables.size(); ++j)
        {
            for (int i=0; i<dataset->target_size; ++i)
            {
                tmp_data=data_out[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j];
                tmp_interpolated_field[i]=dataset->variables[j].a*tmp_data+dataset->variables[j].b;
            }

            dataset->variables[j].data2[fstep]=tmp_interpolated_field;
        }
    }

	xDelete<double>(data_out);

    dataset->reloaded=true;

    //std::cout <<"end load" <<"\n";
}

void
ExternalData::clear()
{
    //xDelete<double>(fcoeff);
}

} // Nextsim
