/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

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

ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId, bool is_vector, double StartingTime )
	:
    M_is_constant( false ),
    M_dataset( dataset ),
    M_VariableId( VariableId ),
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

ExternalData::ExternalData( double ConstantValue )
	:
    M_is_constant( true ),
    M_constant_value( ConstantValue ),
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

ExternalData::ExternalData( double ConstantValue, double ConstantValuebis, double StartingTime, double SpinUpDuration )
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

void ExternalData::check_and_reload(GmshMesh const& mesh, const double current_time)
{
    M_current_time = current_time;

    double current_time_tmp=M_current_time;

    M_factor=1.;
    if((M_current_time-M_StartingTime)<M_SpinUpDuration)
    {
        M_factor=(M_current_time-M_StartingTime)/M_SpinUpDuration;
    }

    if(!M_is_constant)
    {
        bool to_be_reloaded=false;

        if(M_dataset->grid.dataset_frequency=="constant")
            to_be_reloaded=!M_dataset->loaded;
        else if(M_dataset->grid.dataset_frequency=="nearest_daily")
            to_be_reloaded=(to_date_string_yd(current_time)!=to_date_string_yd(M_dataset->ftime_range[0]) || !M_dataset->loaded);
        else
            to_be_reloaded=((current_time_tmp < M_dataset->ftime_range[0]) || (M_dataset->ftime_range[1] < current_time_tmp) || !M_dataset->loaded);            

        if (to_be_reloaded)
        {
            std::cout << "Load " << M_datasetname << "\n";
            loadDataset(M_dataset, mesh);
            M_dataset->interpolated=false;
            std::cout << "Done\n";
        }
        
        if (!M_dataset->interpolated)
        {
            std::cout << "Interpolate " << M_datasetname << "\n";
            interpolateDataset(M_dataset, mesh);
            std::cout << "Done\n";
        }
    }
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
        value = M_constant_value;
    }
    else
    {
        if(M_dataset->grid.dataset_frequency!="constant" && M_dataset->grid.dataset_frequency!="nearest_daily")
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
                (fcoeff[0]*M_dataset->variables[VariableId_tmp].interpolated_data[0][i_tmp] +
                 fcoeff[1]*M_dataset->variables[VariableId_tmp].interpolated_data[1][i_tmp]);
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
            value =  M_factor*M_dataset->variables[VariableId_tmp].interpolated_data[0][i_tmp];
        }
    }

	return static_cast<value_type>( value );
}
typename std::vector<double>
ExternalData::getVector()
{
    std::vector<double> vector_tmp(1,0.);

    if(M_initialized)
    {
        int size_vector=M_dataset->target_size;
        if(M_is_vector)
            size_vector=2*M_dataset->target_size;

        vector_tmp.resize(size_vector);

        for (int i=0; i<size_vector; ++i)
        {
            vector_tmp[i]=(double) get(i);
        }
    }

	return vector_tmp;
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
                              % Environment::vm()["simul.proj_filename"].as<std::string>()
                              ).str();

	std::vector<char> strNextsim(configfileNextsim.begin(), configfileNextsim.end());
	strNextsim.push_back('\0');
	mapNextsim = init_mapx(&strNextsim[0]);

	mapx_class *map;
    double cos_m_diff_angle, sin_m_diff_angle;
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
        dataset->rotation_angle = -(mapNextsim->rotation-map->rotation)*PI/180.;

        close_mapx(map);
    }
    else
    {
        dataset->rotation_angle=0.;
    }

    cos_m_diff_angle=std::cos(-dataset->rotation_angle);
    sin_m_diff_angle=std::sin(-dataset->rotation_angle);

    // ---------------------------------
    // Projection of the mesh positions into the coordinate system of the data before the interpolation
    // (either the lat,lon projection or a polar stereographic projection with another rotaion angle (for ASR))
    // we should need to that also for the TOPAZ native grid, so that we could use a gridtomesh, now we use the latlon of the TOPAZ grid

    auto RX = mesh.coordX(dataset->rotation_angle);
    auto RY = mesh.coordY(dataset->rotation_angle);

    if(dataset->target_size==mesh.numTriangles())
    {
    	RX = mesh.bcoordX(dataset->rotation_angle);
        RY = mesh.bcoordY(dataset->rotation_angle);
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
	std::vector<double> XTIME(1);
	std::vector<size_t> index_start(1);
	std::vector<size_t> index_count(1);

    // size of the data
    int M  = dataset->grid.dimension_y_count;
    int N  = dataset->grid.dimension_x_count;

    int MN = M*N;
    
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
    std::string filename, filename_prev, filename_next;
    std::vector<int> file_jump;

	int index = 0;
    int index_prev=-1;
    int index_next=-1;

    std::vector<std::string> filename_fstep;
    std::vector<int> index_fstep;

    std::string f_timestr;

    // Filename depends on the date for time varying data
	if(dataset->grid.dataset_frequency!="constant" && dataset->grid.dataset_frequency!="nearest_daily")
	{
        // when using forcing from a forecast, we select the file based on the StartingTime
        if ((dataset->grid.prefix).find("start") != std::string::npos)
        {
            ftime = M_StartingTime;
            file_jump.push_back(0);
        }
        else // otherwise, we check for
        {
            ftime = M_current_time-dataset->averaging_period/2.;
            file_jump.push_back(-1);
            file_jump.push_back(0);
            file_jump.push_back(1);
        }

        for (std::vector<int>::iterator jump = file_jump.begin() ; jump != file_jump.end(); ++jump)
        {
            std::string myString;
            if(dataset->grid.dataset_frequency=="monthly")
            {
                f_timestr = to_date_string_ym(std::floor(ftime));

                myString = f_timestr.substr(4,2);
                std::cout <<"month= "<< myString <<"\n";
                int value_month = atoi(myString.c_str());
                myString = f_timestr.substr(0,4);
                std::cout <<"year= "<< myString <<"\n";
                int value_year = atoi(myString.c_str());

                std::cout <<"value_year= "<< value_year <<"\n";
                                std::cout <<"value_month= "<< value_month <<"\n";

                value_month+=*jump;
                if(value_month==13)
                {
                    value_month=1;
                    value_year++;
                }
                if(value_month==0)
                {
                    value_month=12;
                    value_year--;
                }
                f_timestr=(boost::format( "%1%%2%" ) % boost::io::group(std::setw(4), std::setfill('0'), value_year) % boost::io::group(std::setw(2), std::setfill('0'), value_month)).str();
            }
            else if(dataset->grid.dataset_frequency=="yearly")
            {
                f_timestr = to_date_string_y(std::floor(ftime));//yyyy
                int value_year = atoi(f_timestr.c_str());
                value_year+=*jump;
                f_timestr=(boost::format( "%1%" ) % boost::io::group(std::setw(4), std::setfill('0'), value_year)).str();
            }
            else
                f_timestr = to_date_string_yd(std::floor(ftime)+*jump);

            std::cout <<"F_TIMESTR= "<< f_timestr <<"\n";

            filename = (boost::format( "%1%/%2%/%3%%4%%5%" )
                        % Environment::simdataDir().string()
                        % dataset->grid.dirname
                        % dataset->grid.prefix
                        % f_timestr
                        % dataset->grid.postfix
                        ).str();

            std::cout<<"FILENAME= "<< filename <<"\n";
            if ( ! boost::filesystem::exists(filename) )
                continue;
                //throw std::runtime_error("File not found: " + filename);

            index_start.resize(1);
            index_count.resize(1);

            bool has_time_variable;

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
                
                has_time_variable=true;
            } 
            catch(const std::exception& e) // if no time dimension is available in the netcdf, we define the time as 
            {
                index_start[0]=0;
                index_count[0]=1;

                XTIME.resize(index_count[0]);
                
                XTIME[0]=-1.;
                
                has_time_variable=false;
                
                if(dataset->grid.dataset_frequency=="monthly" || dataset->grid.dataset_frequency=="yearly")
                    throw std::runtime_error("The case monthly and yearly when no time dimension is available is not implemented!");
            }
            
            double f;
            for (int it=0; it < XTIME.size(); ++it) // always need one step before and one after the target time
            {
                if (!has_time_variable || ((dataset->name).find("ice_amsr2") != std::string::npos))
                    f = from_date_string((boost::format( "%1%-%2%-%3%" ) % f_timestr.substr(0,4) % f_timestr.substr(4,2) % f_timestr.substr(6,2)).str())+0.5;
                else
                    f = (XTIME[it]*dataset->time.a+dataset->time.b)/24.0+from_date_string(dataset->grid.reference_date);

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
            }
        }

        filename_fstep.push_back(filename_prev);
        index_fstep.push_back(index_prev);

        filename_fstep.push_back(filename_next);
        index_fstep.push_back(index_next);

		dataset->ftime_range.resize(0);
		dataset->ftime_range.push_back(time_prev);
        dataset->ftime_range.push_back(time_next);
	}
    else 
    {
        if(dataset->grid.dataset_frequency=="nearest_daily")
        {
            ftime = M_current_time;
            f_timestr = to_date_string_yd(std::floor(ftime));
            
            double f=from_date_string((boost::format( "%1%-%2%-%3%" ) % f_timestr.substr(0,4) % f_timestr.substr(4,2) % f_timestr.substr(6,2)).str())+0.5;
    		
            dataset->ftime_range.resize(0);
            dataset->ftime_range.push_back(f);
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

    //std::cout<<"Start loading data\n";
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
			if(dataset->variables[j].dimensions.size()>2 && dataset->grid.dataset_frequency!="constant" && dataset->grid.dataset_frequency!="nearest_daily")
			{
            	index_start[0] = index;
            	index_count[0] = 1;
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
                
                if(std::isnan(tmp_data_i))
                    tmp_data_i=0.;
                
                dataset->variables[j].loaded_data[fstep][i]=tmp_data_i;
            }
        }
    

        // ---------------------------------
    
        //std::cout<<"Start transformation on the data\n";
    
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
            }

            if(dataset->rotation_angle!=0.)
            {
                for (int i=0; i<final_MN; ++i)
                {
                    tmp_data0=dataset->variables[j0].loaded_data[fstep][i];
                    tmp_data1=dataset->variables[j1].loaded_data[fstep][i];

                    new_tmp_data0= cos_m_diff_angle*tmp_data0+sin_m_diff_angle*tmp_data1;
                    new_tmp_data1=-sin_m_diff_angle*tmp_data0+cos_m_diff_angle*tmp_data1;

                    dataset->variables[j0].loaded_data[fstep][i]= new_tmp_data0;
                    dataset->variables[j1].loaded_data[fstep][i]= new_tmp_data1;
                }
            }
        }
    }

    // closing maps
    close_mapx(mapNextsim);

    dataset->nb_forcing_step=nb_forcing_step;
    dataset->loaded=true;
}   
    
void
ExternalData::interpolateDataset(Dataset *dataset, GmshMesh const& mesh)//(double const& u, double const& v)
{
    // ---------------------------------
    // Spatial interpolation
    //std::cout<<"Spatial interpolation of the data\n";
    
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
    //std::cout << "Collect the ivariables before interpolation:" <<"\n";
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
                            data_in[(dataset->variables.size()*dataset->nb_forcing_step)*cyclic_i+fstep*dataset->variables.size()+j]=dataset->variables[j].loaded_data[fstep][i];
                    }
                }
                
                if(cyclic_M!=M)
                {
                    for (int x_ind=0; x_ind<N; ++x_ind)
                    {
                        int i=0*N+x_ind;
                        int cyclic_i=(cyclic_M-1)*cyclic_N+x_ind;

                        data_in[(dataset->variables.size()*dataset->nb_forcing_step)*cyclic_i+fstep*dataset->variables.size()+j]=dataset->variables[j].loaded_data[fstep][i];
                    }
                }
                if(cyclic_N!=N)
                {
                    for (int y_ind=0; y_ind<M; ++y_ind)
                    {
                        int i=y_ind*N+0;
                        int cyclic_i=y_ind*cyclic_N+(cyclic_N-1);
                
                        data_in[(dataset->variables.size()*dataset->nb_forcing_step)*cyclic_i+fstep*dataset->variables.size()+j]=dataset->variables[j].loaded_data[fstep][i];
                    }
                }                
            }
            else // with no cyclic dimension, simply use the same indice i
                for (int i=0; i<final_MN; ++i)
                    data_in[(dataset->variables.size()*dataset->nb_forcing_step)*i+fstep*dataset->variables.size()+j]=dataset->variables[j].loaded_data[fstep][i];
            
        }
    }

    // ---------------------------------
    // Projection of the mesh positions into the coordinate system of the data before the interpolation
    // (either the lat,lon projection or a polar stereographic projection with another rotaion angle (for ASR))
    // we should need to that also for the TOPAZ native grid, so that we could use a gridtomesh, now we use the latlon of the TOPAZ grid

    //std::cout << "RX, RY:" <<"\n";
    // Define the mapping and rotation_angle
	mapx_class *mapNextsim;
	std::string configfileNextsim = (boost::format( "%1%/%2%/%3%" )
                              % Environment::nextsimDir().string()
                              % "data"
                              % Environment::vm()["simul.proj_filename"].as<std::string>()
                              ).str();

	std::vector<char> strNextsim(configfileNextsim.begin(), configfileNextsim.end());
	strNextsim.push_back('\0');
	mapNextsim = init_mapx(&strNextsim[0]);

    auto RX = mesh.coordX(dataset->rotation_angle);
    auto RY = mesh.coordY(dataset->rotation_angle);

    if(dataset->target_size==mesh.numTriangles())
    {
    	RX = mesh.bcoordX(dataset->rotation_angle);
        RY = mesh.bcoordY(dataset->rotation_angle);
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
    
    // closing maps
    close_mapx(mapNextsim);


    //std::cout << "Interpolation:" <<"\n";
    
    double* data_out;
    
    switch(dataset->grid.interpolation_method)
    {
        case InterpolationType::FromGridToMesh:
            InterpFromGridToMeshx(  data_out, &dataset->grid.gridX[0], dataset->grid.gridX.size(), &dataset->grid.gridY[0], dataset->grid.gridY.size(),
                                  &data_in[0], dataset->grid.gridY.size(), dataset->grid.gridX.size(),
                                  dataset->variables.size()*dataset->nb_forcing_step,
                                 &RX[0], &RY[0], dataset->target_size, 100000000., dataset->grid.interp_type); // We put an excessively high default value, so that it will most likely crashes when not finding data
        break;
        case InterpolationType::FromMeshToMesh2dx:
            InterpFromMeshToMesh2dx(&data_out,
                                dataset->grid.pfindex,&dataset->grid.gridX[0],&dataset->grid.gridY[0],
                                        dataset->grid.gridX.size(),dataset->grid.pfnels,
                                        &data_in[0],
                                        dataset->grid.gridX.size(),dataset->variables.size()*dataset->nb_forcing_step,
                                        &RX[0], &RY[0], dataset->target_size,
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
            dataset->variables[j].interpolated_data[fstep].resize(dataset->target_size);
            for (int i=0; i<dataset->target_size; ++i)
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
