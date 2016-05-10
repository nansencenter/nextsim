/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   externaldata.cpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#include <externaldata.hpp>
#include <date.hpp>

/**
 * @class ExternalData
 * @brief Manage the external data
 *
 * @see
 *
*/

namespace Nextsim
{

ExternalData::ExternalData( )
{}

ExternalData::ExternalData(Dataset * dataset, GmshMesh const& mesh, int VariableId )
	:
    M_dataset( dataset ),
    M_VariableId( VariableId ),
    M_is_initialized( false ),
    M_current_time( 0. )
{
    M_datasetname = (boost::format( "%1%...%2%" )
                    % M_dataset->prefix
                    % M_dataset->postfix
                    ).str();
    
    fcoeff.resize(2);
}

ExternalData::~ExternalData()
{
	this->clear();
}

void ExternalData::settime( const double current_time )
{
    M_current_time = current_time;
    M_is_initialized = true;
}

void ExternalData::check_and_reload(GmshMesh const& M_mesh)
{
    std::cout << "Before Assert" << "\n";
    
	ASSERT(M_is_initialized, "ExternalData not initialized");
	
    std::cout << "After Assert" << "\n";
    
    bool to_be_reloaded=false;
    
    if(M_dataset->nb_timestep_day>0) 
    {
        to_be_reloaded=((M_current_time < M_dataset->ftime_range[0]) || (M_dataset->ftime_range[1] < M_current_time) || !M_dataset->reloaded);
    }
    else    
        to_be_reloaded=!M_dataset->reloaded;
    
    if (to_be_reloaded)
    {
        std::cout << "Load " << M_datasetname << "\n";
        loadDataset(M_dataset, M_mesh);
        std::cout << "Done\n";
    }
    else
        std::cout << "Nothing to reload\n";
        
}

typename ExternalData::value_type
ExternalData::operator [] (const size_type i)
{
    ASSERT(i < M_dataset->target_size, "invalid index");
    
    value_type value;
    
	if(M_dataset->nb_timestep_day>0)
	{
        fdt = std::abs(M_dataset->ftime_range[1]-M_dataset->ftime_range[0]);
        fcoeff[0] = std::abs(M_current_time-M_dataset->ftime_range[1])/fdt;
        fcoeff[1] = std::abs(M_current_time-M_dataset->ftime_range[0])/fdt;
    
        std::cout <<"LINEAR COEFF 1= "<< fcoeff[0] <<"\n";
        std::cout <<"LINEAR COEFF 2= "<< fcoeff[1] <<"\n";
    
        value =  fcoeff[0]*M_dataset->variables[M_VariableId].data2[0][i] + 
                            fcoeff[1]*M_dataset->variables[M_VariableId].data2[1][i];
    }
    else
    {
        value =  -M_dataset->variables[M_VariableId].data2[0][i];
    }

	return static_cast<value_type>( value );
}

void
ExternalData::loadDataset(Dataset *dataset, GmshMesh const& M_mesh)//(double const& u, double const& v)
{

    // Load grid if unloaded
    if(!dataset->grid->loaded)
        loadGrid(dataset->grid);
        
    // Initialise counters etc.
	std::string current_timestr = "";
	int nb_forcing_step =1;

	std::vector<double> XTIME(1);
	std::vector<size_t> index_start(1);
	std::vector<size_t> index_count(1);

	int index = 0;

	// interp_type for grid to mesh interpolation
	int interp_type = dataset->grid->interp_type;

    // create dataset->ftime_range for data sets which need to be interpolated in time
	if(dataset->nb_timestep_day>0)
	{
		current_timestr = to_date_string_ym(M_current_time);

		double file_dt = 1./dataset->nb_timestep_day;
		double time_start = std::floor(M_current_time*dataset->nb_timestep_day)/dataset->nb_timestep_day;
		double time_end = std::ceil(M_current_time*dataset->nb_timestep_day)/dataset->nb_timestep_day;

		// We always need at least two time steps to interpolate between
		if (time_end == time_start)
		{
			time_end = time_start + (1./dataset->nb_timestep_day);
		}

		dataset->ftime_range.resize(0);
		for (double dt=time_start; dt<=time_end; dt+=file_dt)
		{
			dataset->ftime_range.push_back(dt);
		}

		for (int i=0; i<dataset->ftime_range.size(); ++i)
		{
			std::cout <<"TIMEVEC["<< i <<"]= "<< dataset->ftime_range[i] <<"\n";
		}

		nb_forcing_step = dataset->ftime_range.size();
	}

    // Initialise variables for the fields
    std::vector<double> tmp_interpolated_field(dataset->target_size);

    int N_data =dataset->variables.size();
    int M  =dataset->grid->dimension_y.end;
    int N  = dataset->grid->dimension_x.end;
    int MN = M*N;

	int reduced_MN=MN;
	if(dataset->grid->reduced_nodes_ind.size()!=0)
		reduced_MN=dataset->grid->reduced_nodes_ind.size();

	// Memory leak:
    //double* data_in = new double[N_data*nb_forcing_step*reduced_MN];
    std::vector<double> data_in(N_data*nb_forcing_step*reduced_MN);

    std::vector<double> data_in_tmp(MN);

    // Attributes (scaling and offset)
    netCDF::NcVarAtt att;
    double scale_factor;
    double add_offset;

    std::cout <<"NB_FORCING_STEP= "<< nb_forcing_step <<"\n";

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
            std::string f_timestr = to_date_string_ym(std::floor(ftime));
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

        // Open the netcdf file
        netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

        // Find the right time slice
        if (dataset->nb_timestep_day>0)
        {
            // Set the time range XTIME
            netCDF::NcVar FVTIME = dataFile.getVar(dataset->time.name);
            XTIME.resize(dataset->time.dimensions[0].end-dataset->time.dimensions[0].start);
            FVTIME.getVar(&XTIME[0]);
            std::for_each(XTIME.begin(), XTIME.end(), [&](double& f){ f = f/24.0+from_date_string(dataset->reference_date); });

            auto it = std::find(XTIME.begin(), XTIME.end(), ftime);
            index = std::distance(XTIME.begin(),it);
            std::cout <<"FIND "<< ftime <<" in index "<< index <<"\n";
        }

        for(int j=0; j<dataset->variables.size(); ++j)
        {
            NcVars[j] = dataFile.getVar(dataset->variables[j].name);

            index_start.resize(dataset->variables[j].dimensions.size());
            index_count.resize(dataset->variables[j].dimensions.size());

            for(int k=0; k<dataset->variables[j].dimensions.size(); ++k)
            {
                index_start[k] = dataset->variables[j].dimensions[k].start;
                index_count[k] = dataset->variables[j].dimensions[k].end-dataset->variables[j].dimensions[k].start;
            }

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

			if(dataset->grid->reduced_nodes_ind.size()!=0)
			{
            	for (int i=0; i<(reduced_MN); ++i)
                	data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j]=data_in_tmp[dataset->grid->reduced_nodes_ind[i]]*scale_factor + add_offset;
			}
			else
            	for (int i=0; i<(MN); ++i)
                	data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j]=data_in_tmp[i]*scale_factor + add_offset;

		}
    }

    double* data_out;
    double tmp_data;

    auto RX = M_mesh.coordX(dataset->grid->rotation_angle);
    auto RY = M_mesh.coordY(dataset->grid->rotation_angle);

    if(dataset->target_size==M_mesh.numTriangles())
    {
    	RX = M_mesh.bcoordX(dataset->grid->rotation_angle);
        RY = M_mesh.bcoordY(dataset->grid->rotation_angle);
    }

	if(dataset->grid->interpolation_in_latlon)
	{
		mapx_class *map;
		std::string configfile = (boost::format( "%1%/%2%/%3%" )
                                  % Environment::nextsimDir().string()
                                  % "data"
                                  % "NpsNextsim.mpp"
                                  ).str();

		std::vector<char> str(configfile.begin(), configfile.end());
		str.push_back('\0');
		map = init_mapx(&str[0]);

		double lat, lon;

		for (int i=0; i<dataset->target_size; ++i)
		{
			inverse_mapx(map,RX[i],RY[i],&lat,&lon);
			RY[i]=lat;
			RX[i]=lon;
			//tmp_latlon = XY2latLon(RX[i], RY[i], map, configfile);
			//RY[i]=tmp_latlon[0];
			//RX[i]=tmp_latlon[1];
		}

		close_mapx(map);
	}

    std::cout <<"before interp " <<"\n";

    switch(dataset->grid->interpolation_method)
    {
        case InterpolationType::FromGridToMesh:
            InterpFromGridToMeshx(  data_out, &dataset->grid->gridX[0], dataset->grid->gridX.size(), &dataset->grid->gridY[0], dataset->grid->gridY.size(),
                                  &data_in[0], dataset->grid->gridY.size(), dataset->grid->gridX.size(),
                                  dataset->variables.size()*nb_forcing_step,
                                 &RX[0], &RY[0], dataset->target_size, 1.0, interp_type);
        break;
        case InterpolationType::FromMeshToMesh2dx:
            InterpFromMeshToMesh2dx(&data_out,
                                dataset->grid->pfindex,&dataset->grid->gridX[0],&dataset->grid->gridY[0],
                                        dataset->grid->gridX.size(),dataset->grid->pfnels,
                                        &data_in[0],
                                        dataset->grid->gridX.size(),N_data*nb_forcing_step,
                                        &RX[0], &RY[0], dataset->target_size,
                                        false);
        break;
        default:
            std::cout << "invalid interpolation type:" <<"\n";
            throw std::logic_error("invalid interpolation type");
    }

std::cout <<"after interp " <<"\n";

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

    std::cout <<"end load" <<"\n";
}

void
ExternalData::loadGrid(Grid *grid)
{
    std::string current_timestr = to_date_string_ym(M_current_time);
    std::cout <<"TIMESTR= "<< current_timestr <<"\n";
    std::string filename = (boost::format( "%1%/%2%/%3%" )
                            % Environment::simdataDir().string()
                            % grid->dirname
                            % grid->filename
                            ).str();

    //switch (grid->latitude.dimensions.size())
    //{
    //    case 1:
	if(grid->latitude.dimensions.size()==1)
	{
		// read in coordinates
		std::vector<size_t> index_x_count(1);
		std::vector<size_t> index_y_count(1);

		std::vector<size_t> index_x_start(1);
		std::vector<size_t> index_y_start(1);

		index_y_start[0] = grid->dimension_y.start;
		index_y_count[0] = grid->dimension_y.end-grid->dimension_y.start;

		index_x_start[0] = grid->dimension_x.start;
		index_x_count[0] = grid->dimension_x.end-grid->dimension_x.start;

		std::vector<double> LAT(index_y_count[0]);
		std::vector<double> LON(index_x_count[0]);

		std::cout <<"GRID : READ NETCDF starts\n";
                if ( ! boost::filesystem::exists(filename) )
                    throw std::runtime_error("File not found: " + filename);
		netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
		netCDF::NcVar VLAT = dataFile.getVar(grid->latitude.name);
		netCDF::NcVar VLON = dataFile.getVar(grid->longitude.name);
		std::cout <<"GRID : READ NETCDF done\n";

		VLAT.getVar(index_y_start,index_y_count,&LAT[0]);
		VLON.getVar(index_x_start,index_x_count,&LON[0]);

		grid->gridY=LAT;
		grid->gridX=LON;
	}
	else
	{
//		break;
//    	case 2:
		// read in coordinates
		std::vector<size_t> index_px_count(2);
		std::vector<size_t> index_py_count(2);

		std::vector<size_t> index_px_start(2);
		std::vector<size_t> index_py_start(2);

		index_py_start[0] = grid->dimension_y.start;
		index_py_start[1] = grid->dimension_x.start;

		index_py_count[0] = grid->dimension_y.end-grid->dimension_y.start;
		index_py_count[1] = grid->dimension_x.end-grid->dimension_x.start;

		index_px_start[0] = grid->dimension_y.start;
		index_px_start[1] = grid->dimension_x.start;

		index_px_count[0] = grid->dimension_y.end-grid->dimension_y.start;
		index_px_count[1] = grid->dimension_x.end-grid->dimension_x.start;

		if(grid->interpolation_method==InterpolationType::FromGridToMesh)
		{
			index_py_count[1] = 1;
			index_px_count[0] = 1;
		}

		std::vector<double> XLAT(index_px_count[0]*index_px_count[1]);
		std::vector<double> XLON(index_px_count[0]*index_px_count[1]);
		std::vector<double> YLAT(index_py_count[0]*index_py_count[1]);
		std::vector<double> YLON(index_py_count[0]*index_py_count[1]);

		std::cout <<"GRID : READ NETCDF starts\n";
                if ( ! boost::filesystem::exists(filename) )
                    throw std::runtime_error("File not found: " + filename);
		netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
		netCDF::NcVar VLAT = dataFile.getVar(grid->latitude.name);
		netCDF::NcVar VLON = dataFile.getVar(grid->longitude.name);
		std::cout <<"GRID : READ NETCDF done\n";

		VLAT.getVar(index_px_start,index_px_count,&XLAT[0]);
		VLON.getVar(index_px_start,index_px_count,&XLON[0]);

		VLAT.getVar(index_py_start,index_py_count,&YLAT[0]);
		VLON.getVar(index_py_start,index_py_count,&YLON[0]);

		std::vector<double> X(index_px_count[0]*index_px_count[1]);
		std::vector<double> Y(index_py_count[0]*index_py_count[1]);

		mapx_class *map;
		std::string configfile = (boost::format( "%1%/%2%/%3%" )
                                  % Environment::nextsimDir().string()
                                  % grid->dirname
                                  % grid->mpp_file
                                  ).str();

		std::vector<char> str(configfile.begin(), configfile.end());
		str.push_back('\0');
		map = init_mapx(&str[0]);

	    double x;
	    double y;

		for (int i=0; i<index_px_count[0]; ++i)
		{
			for (int j=0; j<index_px_count[1]; ++j)
			{
			    forward_mapx(map,XLAT[index_px_count[1]*i+j],XLON[index_px_count[1]*i+j],&x,&y);
				X[index_px_count[1]*i+j]=x;
			}
		}

		for (int i=0; i<index_py_count[0]; ++i)
		{
			for (int j=0; j<index_py_count[1]; ++j)
			{
				forward_mapx(map,YLAT[index_py_count[1]*i+j],YLON[index_py_count[1]*i+j],&x,&y);
				Y[index_py_count[1]*i+j]=y;
			}
		}

		close_mapx(map);

		if(grid->interpolation_method==InterpolationType::FromMeshToMesh2dx)
		{
			if(grid->masking){
				netCDF::NcVar VMASK;

				VMASK = dataFile.getVar(grid->masking_variable.name);

				std::vector<double> data_in;

				std::vector<double> reduced_FX;
				std::vector<double> reduced_FY;
				std::vector<int> reduced_nodes_ind;

				std::vector<size_t> index_start(3,0);
				std::vector<size_t> index_count(3);

				index_start.resize(grid->masking_variable.dimensions.size());
				index_count.resize(grid->masking_variable.dimensions.size());

				for(int k=0; k<grid->masking_variable.dimensions.size(); ++k)
				{
					index_start[k] = grid->masking_variable.dimensions[k].start;
					index_count[k] = grid->masking_variable.dimensions[k].end-grid->masking_variable.dimensions[k].start;
				}
				index_start[0] = 0;
				index_count[0] = 1;

				if((index_px_count[0]!=index_count[grid->masking_variable.dimensions.size()-2]) || (index_px_count[1]!=index_count[grid->masking_variable.dimensions.size()-1]))
				{
                    std::cout << "index_px_count[0] = " << index_px_count[0] << " index_count[grid->masking_variable.dimensions.size()-2] = " << index_count[grid->masking_variable.dimensions.size()-2] <<"\n";
					std::cout << "index_px_count[1] = " << index_px_count[1] << " index_count[grid->masking_variable.dimensions.size()-1] = " << index_count[grid->masking_variable.dimensions.size()-1] <<"\n";
                    throw std::logic_error("Not the same dimension for the masking variable and the grid!!");
				}

				data_in.resize(index_px_count[0]*index_px_count[1]);
				VMASK.getVar(index_start,index_count,&data_in[0]);

				netCDF::NcVarAtt att;
				int FillValue;

				att = VMASK.getAtt("_FillValue");
				att.getValues(&FillValue);

				for (int i=0; i<index_px_count[0]; ++i)
				{
					for (int j=0; j<index_px_count[1]; ++j)
					{
						if (data_in[index_px_count[1]*i+j] != FillValue)
						{
							reduced_FX.push_back(X[index_px_count[1]*i+j]);
							reduced_FY.push_back(Y[index_px_count[1]*i+j]);
							reduced_nodes_ind.push_back(index_px_count[1]*i+j);
						}
					}
				}
				grid->gridX=reduced_FX;
				grid->gridY=reduced_FY;
				grid->reduced_nodes_ind=reduced_nodes_ind;
			}
			else // no masking of the Filled Value
			{
				grid->gridX=X;
				grid->gridY=Y;
			}

			std::cout <<"GRID : Triangulate starts\n";
			BamgTriangulatex(&grid->pfindex,&grid->pfnels,&grid->gridX[0],&grid->gridY[0],grid->gridX.size());
			std::cout <<"GRID : NUMTRIANGLES= "<< grid->pfnels <<"\n";
			std::cout <<"GRID : Triangulate done\n";
		}
		else
		{
			grid->gridX=X;
			grid->gridY=Y;
		}
        
	//	break;
	//
    //default:
    //   std::cout << "invalid ocean initialisation"<<"\n";
    //    throw std::logic_error("invalid ocean forcing");
	}
    
    grid->loaded=true;
}

void
ExternalData::clear()
{
    //xDelete<double>(fcoeff);
}

} // Nextsim
