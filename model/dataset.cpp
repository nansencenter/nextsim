/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   dataset.cpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#include <dataset.hpp>
#include <constants.hpp>

/**
 * @class DataSet
 * @brief Describe the datasets
 *
 * @see
 *
*/

namespace Nextsim
{

   
    DataSet::DataSet( )
    {}
   
   DataSet::DataSet(char *DatasetName, int target_size)
   {
//     Dataset *this;

     std::vector<std::vector<double>> data2_tmp;
     data2_tmp.resize(2);

     /*
      *	match projection name and initialize remaining parameters
      */
     if (strcmp (DatasetName, "asr_nodes") == 0)
     { 
     	// Definition of asr grid and datasets
         Dimension dimension_x={
             name:"x",
             cyclic:false
     	};

         Dimension dimension_y={
             name:"y",
             cyclic:false
     	};

         Dimension dimension_time={
             name:"time", // "Time"
             cyclic:false
     	};

         std::vector<Dimension> dimensions_latlon(2);
         dimensions_latlon[0] = dimension_y;
         dimensions_latlon[1] = dimension_x;

         std::vector<Dimension> dimensions(3);
         dimensions[0] = dimension_time;
         dimensions[1] = dimension_y;
         dimensions[2] = dimension_x;

         std::vector<Dimension> dimensions_time(1);
         dimensions_time[0] = dimension_time;

         Variable latitude={
             name: "XLAT",
             dimensions: dimensions_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp
     	};

         Variable longitude={
             name: "XLONG",
             dimensions: dimensions_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp
     	};

         Variable time={
             name: "time",
             dimensions: dimensions_time,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp
     	};

         // conversion factors: xnew = a*x + b
         Variable u={
             name: "U10", // U10M
             dimensions: dimensions,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v={
             name: "V10", // U10M
             dimensions: dimensions,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};
    
     	Grid M_grid={
     		interpolation_method: InterpolationType::FromGridToMesh,
     	    //interp_type : TriangleInterpEnum,  // slower
     	    interp_type : BilinearInterpEnum,
     	    //interp_type : NearestInterpEnum,

     		this->dirname="data",
     		filename:"asr30km.comb.2d.200803.nc",

     		latitude: latitude,
     		longitude: longitude,

     		dimension_x: dimension_x,
     		dimension_y: dimension_y,

     		mpp_file: "NpsASR.mpp",
      		interpolation_in_latlon: false,

             loaded: false,

     		masking: false
     	};

         std::vector<Variable> variables(2);
             variables[0] = u;
             variables[1] = v;
            
         std::vector<int> uv_tmp(2);
             uv_tmp[0] = 0;
             uv_tmp[1] = 1;
    
         Vectorial_Variable uv={
             components_Id: uv_tmp,
             east_west_oriented: false // if false, then we assume it is oriented following the input grid
         };
            
         std::vector<Vectorial_Variable> vectorial_variables(1);
         vectorial_variables[0] = uv;

             this->dirname= "data";
             this->prefix= "asr30km.comb.2d."; // "asr30km.comb.2D.";
             this->postfix=".nc";
             this->reference_date= "1901-01-01";

             this->variables= variables;
             this->vectorial_variables= vectorial_variables;
             this->target_size= target_size;
             this->grid= M_grid;

             this->reloaded=false;

             this->nb_timestep_day= 8;
             this->time= time;
}
     else if (strcmp (DatasetName, "asr_elements") == 0)
     { 
     	// Definition of asr grid and datasets
         Dimension dimension_x={
             name:"x",
             cyclic:false
     	};

         Dimension dimension_y={
             name:"y",
             cyclic:false
     	};

         Dimension dimension_time={
             name:"time", // "Time"
             cyclic:false
     	};

         std::vector<Dimension> dimensions_latlon(2);
         dimensions_latlon[0] = dimension_y;
         dimensions_latlon[1] = dimension_x;

         std::vector<Dimension> dimensions(3);
         dimensions[0] = dimension_time;
         dimensions[1] = dimension_y;
         dimensions[2] = dimension_x;

         std::vector<Dimension> dimensions_time(1);
         dimensions_time[0] = dimension_time;

         Variable latitude={
             name: "XLAT",
             dimensions: dimensions_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp
     	};

         Variable longitude={
             name: "XLONG",
             dimensions: dimensions_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp
     	};

         Variable time={
             name: "time",
             dimensions: dimensions_time,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp
     	};

         // conversion factors: xnew = a*x + b
         Variable u={
             name: "U10", // U10M
             dimensions: dimensions,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v={
             name: "V10", // U10M
             dimensions: dimensions,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};
    
     	Grid M_grid={
     		interpolation_method: InterpolationType::FromGridToMesh,
     	    //interp_type : TriangleInterpEnum,  // slower
     	    interp_type : BilinearInterpEnum,
     	    //interp_type : NearestInterpEnum,

     		this->dirname="data",
     		filename:"asr30km.comb.2d.200803.nc",

     		latitude: latitude,
     		longitude: longitude,

     		dimension_x: dimension_x,
     		dimension_y: dimension_y,

     		mpp_file: "NpsASR.mpp",
      		interpolation_in_latlon: false,

             loaded: false,

     		masking: false
     	};

         int nb_timestep_day=8;

         Variable tair={
             name:"T2",
             dimensions: dimensions,
             a:1.,
             b:-273.15,
             Units:"C",
             data2: data2_tmp
     	}; // T2M
         Variable mixrat={
             name:"Q2",
             dimensions: dimensions,
             a:1.,
             b:0.,
             Units:"",
             data2: data2_tmp
     	}; // Q2M
         Variable mslp={
             name:"SLP",
             dimensions: dimensions,
             a:1e2,
             b:0.,
             Units:"Pa",
             data2: data2_tmp
     	}; //PSFC, a=1.
         Variable Qsw_in={
             name:"SWDNB",
             dimensions: dimensions,
             a:1.,
             b:0.,
             Units:"W/m^2",
             data2: data2_tmp
     	};
         Variable Qlw_in={
             name:"LWDNB",
             dimensions: dimensions,
             a:1.,
             b:0.,
             Units:"W/m^2",
             data2: data2_tmp
     	};
         Variable snowfr={
             name:"SR",
             dimensions: dimensions,
             a:1.,
             b:0.,
             Units:"",
             data2: data2_tmp
     	};
         Variable precip={
             name:"RAINNC",
             dimensions: dimensions,
             a:physical::rhow/1000.*(nb_timestep_day)/(24.*3600),
             b:0.,
             Units:"kg/m^2/s",
             data2: data2_tmp
     	};

         std::vector<Variable> variables(7);
         variables[0] = tair;
         variables[1] = mixrat;
         variables[2] = mslp;
         variables[3] = Qsw_in;
         variables[4] = Qlw_in;
         variables[5] = snowfr;
         variables[6] = precip;
    
         std::vector<Vectorial_Variable> vectorial_variables(0);

             this->dirname="data";
             this->prefix="asr30km.comb.2d."; // "asr30km.comb.2D.";
             this->postfix=".nc";
             this->reference_date= "1901-01-01";

             this->variables= variables;
             this->vectorial_variables= vectorial_variables;
             this->target_size=target_size;
             this->grid= M_grid;

             this->reloaded=false;

             this->nb_timestep_day= 8;
             this->time= time;
     }
     else if (strcmp (DatasetName, "topaz_nodes") == 0)
     { 
     	// Definition of topaz grid and datasets
         Dimension dimension_x={
             name:"x",
             cyclic:false
     	};

         Dimension dimension_y={
             name:"y",
             cyclic:false
     	};

         Dimension dimension_time={
             name:"time", // "Time"
             cyclic:false
     	};

         Dimension dimension_depth={
             name:"depth", // "Time"
             cyclic:false
     	};

         std::vector<Dimension> dimensions_uv(4);
         dimensions_uv[0] = dimension_time;
         dimensions_uv[1] = dimension_depth;
         dimensions_uv[2] = dimension_y;
         dimensions_uv[3] = dimension_x;

         std::vector<Dimension> dimensions(3);
         dimensions[0] = dimension_time;
         dimensions[1] = dimension_y;
         dimensions[2] = dimension_x;

         std::vector<Dimension> dimensions_latlon(2);
         dimensions_latlon[0] = dimension_y;
         dimensions_latlon[1] = dimension_x;

         std::vector<Dimension> dimensions_time(1);
         dimensions_time[0] = dimension_time;

         Variable latitude={
             name: "latitude",
             dimensions: dimensions_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp};

         Variable longitude={
             name: "longitude",
             dimensions: dimensions_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp};

         Variable time={
             name: "time",
             dimensions: dimensions_time,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp};

         Variable u={
             name: "u",
             dimensions: dimensions_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v={
             name: "v",
             dimensions: dimensions_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

     	Variable ssh={
     		name: "ssh",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m/s",
     		data2: data2_tmp
     	};

     	Variable sst={
     		name: "temperature",
     		dimensions: dimensions_uv,
     		a: 1.,
     		b: 0.,
     		Units: "deg celsius",
     		data2: data2_tmp
     	};

     	Variable sss={
     		name: "salinity",
     		dimensions: dimensions_uv,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable mld={
     		name: "mlp",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable conc={
     		name: "fice",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable thick={
     		name: "hice",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable snow_thick={
     		name: "hsnow",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

         Grid M_grid={
             interpolation_method: InterpolationType::FromMeshToMesh2dx,
     		interp_type: -1,
             this->dirname= "data",
             filename: "TP4DAILY_200803_3m.nc",

             latitude: latitude,
             longitude: longitude,

             dimension_x: dimension_x,
             dimension_y: dimension_y,

             mpp_file: "NpsNextsim.mpp",
     		interpolation_in_latlon: false,

             loaded: false,

     		masking: true,
     		masking_variable: sss
     	};

         std::vector<Variable> variables(3);
         variables[0] = u;
         variables[1] = v;
         variables[2] = ssh;
    
         std::vector<int> uv_tmp(2);
             uv_tmp[0] = 0;
             uv_tmp[1] = 1;
    
         Vectorial_Variable uv={
             components_Id: uv_tmp,
             east_west_oriented: false // if false, then we assume it is oriented following the mpp_file defined for the grid
     	};
    
         std::vector<Vectorial_Variable> vectorial_variables(1);
         vectorial_variables[0] = uv;
    
             this->dirname= "data";
             this->prefix= "TP4DAILY_";
             this->postfix= "_30m.nc";
             this->reference_date= "1950-01-01";

             this->variables= variables;
             this->vectorial_variables= vectorial_variables;
             this->target_size= target_size;
             this->grid= M_grid;

             this->reloaded=false;

             this->nb_timestep_day= 1;
             this->time= time;
     }
     else if (strcmp (DatasetName, "topaz_elements") == 0)
     { 
     	// Definition of topaz grid and datasets
         Dimension dimension_x={
             name:"x",
             cyclic:false
     	};

         Dimension dimension_y={
             name:"y",
             cyclic:false
     	};

         Dimension dimension_time={
             name:"time", // "Time"
             cyclic:false
     	};

         Dimension dimension_depth={
             name:"depth", // "Time"
             cyclic:false
     	};

         std::vector<Dimension> dimensions_uv(4);
         dimensions_uv[0] = dimension_time;
         dimensions_uv[1] = dimension_depth;
         dimensions_uv[2] = dimension_y;
         dimensions_uv[3] = dimension_x;

         std::vector<Dimension> dimensions(3);
         dimensions[0] = dimension_time;
         dimensions[1] = dimension_y;
         dimensions[2] = dimension_x;

         std::vector<Dimension> dimensions_latlon(2);
         dimensions_latlon[0] = dimension_y;
         dimensions_latlon[1] = dimension_x;

         std::vector<Dimension> dimensions_time(1);
         dimensions_time[0] = dimension_time;

         Variable latitude={
             name: "latitude",
             dimensions: dimensions_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp};

         Variable longitude={
             name: "longitude",
             dimensions: dimensions_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp};

         Variable time={
             name: "time",
             dimensions: dimensions_time,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp};

         Variable u={
             name: "u",
             dimensions: dimensions_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v={
             name: "v",
             dimensions: dimensions_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

     	Variable ssh={
     		name: "ssh",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m/s",
     		data2: data2_tmp
     	};

     	Variable sst={
     		name: "temperature",
     		dimensions: dimensions_uv,
     		a: 1.,
     		b: 0.,
     		Units: "deg celsius",
     		data2: data2_tmp
     	};

     	Variable sss={
     		name: "salinity",
     		dimensions: dimensions_uv,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable mld={
     		name: "mlp",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable conc={
     		name: "fice",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable thick={
     		name: "hice",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable snow_thick={
     		name: "hsnow",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

         Grid M_grid={
             interpolation_method: InterpolationType::FromMeshToMesh2dx,
     		interp_type: -1,
             this->dirname= "data",
             filename: "TP4DAILY_200803_3m.nc",

             latitude: latitude,
             longitude: longitude,

             dimension_x: dimension_x,
             dimension_y: dimension_y,

             mpp_file: "NpsNextsim.mpp",
     		interpolation_in_latlon: false,

             loaded: false,

     		masking: true,
     		masking_variable: sss
     	};

         std::vector<Variable> variables(3);
         variables[0] = sst;
         variables[1] = sss;
         variables[2] = mld;
    
         std::vector<Vectorial_Variable> vectorial_variables(0);

             this->dirname= "data";
             this->prefix= "TP4DAILY_";
             this->postfix= "_3m.nc";
             this->reference_date= "1950-01-01";

             this->variables= variables;
             this->vectorial_variables= vectorial_variables;
             this->target_size= target_size;
             this->grid= M_grid;

             this->reloaded=false;

             this->nb_timestep_day= 1;
             this->time= time;
     }
     else if (strcmp (DatasetName, "ice_topaz_elements") == 0)
     { 
     	// Definition of topaz grid and datasets
         Dimension dimension_x={
             name:"x",
             cyclic:false
     	};

         Dimension dimension_y={
             name:"y",
             cyclic:false
     	};

         Dimension dimension_time={
             name:"time", // "Time"
             cyclic:false
     	};

         Dimension dimension_depth={
             name:"depth", // "Time"
             cyclic:false
     	};

         std::vector<Dimension> dimensions_uv(4);
         dimensions_uv[0] = dimension_time;
         dimensions_uv[1] = dimension_depth;
         dimensions_uv[2] = dimension_y;
         dimensions_uv[3] = dimension_x;

         std::vector<Dimension> dimensions(3);
         dimensions[0] = dimension_time;
         dimensions[1] = dimension_y;
         dimensions[2] = dimension_x;

         std::vector<Dimension> dimensions_latlon(2);
         dimensions_latlon[0] = dimension_y;
         dimensions_latlon[1] = dimension_x;

         std::vector<Dimension> dimensions_time(1);
         dimensions_time[0] = dimension_time;

         Variable latitude={
             name: "latitude",
             dimensions: dimensions_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp};

         Variable longitude={
             name: "longitude",
             dimensions: dimensions_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp};

         Variable time={
             name: "time",
             dimensions: dimensions_time,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp};

         Variable u={
             name: "u",
             dimensions: dimensions_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v={
             name: "v",
             dimensions: dimensions_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

     	Variable ssh={
     		name: "ssh",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m/s",
     		data2: data2_tmp
     	};

     	Variable sst={
     		name: "temperature",
     		dimensions: dimensions_uv,
     		a: 1.,
     		b: 0.,
     		Units: "deg celsius",
     		data2: data2_tmp
     	};

     	Variable sss={
     		name: "salinity",
     		dimensions: dimensions_uv,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable mld={
     		name: "mlp",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable conc={
     		name: "fice",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable thick={
     		name: "hice",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable snow_thick={
     		name: "hsnow",
     		dimensions: dimensions,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

         Grid M_grid={
             interpolation_method: InterpolationType::FromMeshToMesh2dx,
     		interp_type: -1,
             this->dirname= "data",
             filename: "TP4DAILY_200803_3m.nc",

             latitude: latitude,
             longitude: longitude,

             dimension_x: dimension_x,
             dimension_y: dimension_y,

             mpp_file: "NpsNextsim.mpp",
     		interpolation_in_latlon: false,

             loaded: false,

     		masking: true,
     		masking_variable: sss
     	};

         std::vector<Variable> variables(3);
         variables[0] = conc;
         variables[1] = thick;
         variables[2] = snow_thick;
    
         std::vector<Vectorial_Variable> vectorial_variables(0);
    
             this->dirname= "data";
             this->prefix= "TP4DAILY_";
             this->postfix= "_3m.nc";
             this->reference_date= "1950-01-01";

             this->variables= variables;
             this->vectorial_variables= vectorial_variables;
             this->target_size= target_size;
             this->grid= M_grid;

             this->reloaded=false;

     		this->nb_timestep_day= 1;
     	    this->time= time;
     }
     else if (strcmp (DatasetName, "etopo_elements") == 0)
     { 
     	// Definition of etopo grid and datasets
         Dimension dimension_x={
             name:"x",
             cyclic:false
     	};

         Dimension dimension_y={
             name:"y",
             cyclic:false
     	};

         std::vector<Dimension> dimensions_lon(1);
         dimensions_lon[0] = dimension_x;

         std::vector<Dimension> dimensions_lat(1);
         dimensions_lat[0] = dimension_y;

         std::vector<Dimension> dimensions(2);
         dimensions[0] = dimension_y;
         dimensions[1] = dimension_x;

         Variable latitude={
             name: "y",
             dimensions: dimensions_lat,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp
     	};

         Variable longitude={
             name: "x",
             dimensions: dimensions_lon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp
     	};

     	Grid M_grid={
     		interpolation_method: InterpolationType::FromGridToMesh,
     	    //interp_type : TriangleInterpEnum, // slower
     	    interp_type : BilinearInterpEnum,
     	    //interp_type : NearestInterpEnum,
     		this->dirname="data",
     		filename:"ETOPO1_Ice_g_gmt4.grd",

     		latitude: latitude,
     		longitude: longitude,

     		dimension_x: dimension_x,
     		dimension_y: dimension_y,

     		mpp_file: "",
     		interpolation_in_latlon: true,

             loaded: false,

     		masking: false
     	};

         Variable z={
             name:"z",
             dimensions: dimensions,
             a:-1.,
             b:0.,
             Units:"m",
             data2: data2_tmp
     	};

         std::vector<Variable> variables(1);
         variables[0] = z;
    
         std::vector<Vectorial_Variable> vectorial_variables(0);

         
             this->dirname="data";
             this->prefix="ETOPO1_Ice_g_gmt4";
             this->postfix=".grd";
             this->reference_date= "";

             this->variables= variables;
             this->vectorial_variables= vectorial_variables;
             this->target_size=target_size;
             this->grid= M_grid;

             this->reloaded=false;

             this->nb_timestep_day= 0;
     }
     else if (strcmp (DatasetName, "ERAi_elements") == 0)
     { 
     	// Definition of dimensions
        dimension_x={
             name:"lon",
             cyclic:true
     	};

        dimension_y={
             name:"lat",
             cyclic:false
     	};
    
        dimension_time={
             name:"time", // "Time"
             cyclic:false
     	};
        
        // Definition of the grid
         std::vector<Dimension> dimensions_lon(1);
         dimensions_lon[0] = this->dimension_x;

         std::vector<Dimension> dimensions_lat(1);
         dimensions_lat[0] = this->dimension_y;

         Variable latitude={
             name: "lat",
             dimensions: dimensions_lat,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp
     	};

         Variable longitude={
             name: "lon",
             dimensions: dimensions_lon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp
     	};

     	Grid M_grid={
     		interpolation_method: InterpolationType::FromGridToMesh,
     	    //interp_type : TriangleInterpEnum, // slower
     	    interp_type : BilinearInterpEnum,
     	    //interp_type : NearestInterpEnum,
     		this->dirname="data",
     		filename:"erai.6h.200803.nc",

     		latitude: latitude,
     		longitude: longitude,

     		dimension_x: dimension_x,
     		dimension_y: dimension_y,

     		mpp_file: "",
     		interpolation_in_latlon: true,

             loaded: false,

     		masking: false
     	};
        
        this->grid= M_grid;
        
        // Definition of the data
        
        std::vector<Dimension> dimensions_time(1);
        dimensions_time[0] = this->dimension_time;
   
        std::vector<Dimension> dimensions(3);
        dimensions[0] = this->dimension_time;
        dimensions[1] = this->dimension_y;
        dimensions[2] = this->dimension_x;
        
        Variable time={
            name: "time",
            dimensions: dimensions_time,
            a: 1.,
            b: 0.,
            Units: "hours",
            data2: data2_tmp
    	};
        
        int nb_timestep_day=4;

        Variable tair={
            name:"2T",
            dimensions: dimensions,
            a:1.,
            b:-273.15,
            Units:"C",
            data2: data2_tmp
    	}; // T2M
        Variable dair={
            name:"2D",
            dimensions: dimensions,
            a:1.,
            b:-273.15,
            Units:"C",
            data2: data2_tmp
    	}; // Q2M
        Variable mslp={
            name:"MSL",
            dimensions: dimensions,
            a:1.,
            b:0.,
            Units:"Pa",
            data2: data2_tmp
    	}; //PSFC, a=1.
        Variable Qsw_in={
            name:"SSRD",
            dimensions: dimensions,
            a:nb_timestep_day/(24.*3600),
            b:0.,
            Units:"W/m^2",
            data2: data2_tmp
    	};
        Variable tcc={
            name:"TCC",
            dimensions: dimensions,
            a:1.,
            b:0.,
            Units:"",
            data2: data2_tmp
    	};

        Variable precip={
            name:"TP",
            dimensions: dimensions,
            a:physical::rhow*(nb_timestep_day)/(24.*3600),
            b:0.,
            Units:"kg/m^2/s",
            data2: data2_tmp
    	};

        std::vector<Variable> variables(6);
        variables[0] = tair;
        variables[1] = dair;
        variables[2] = mslp;
        variables[3] = Qsw_in;
        variables[4] = tcc;
        variables[5] = precip;

         std::vector<Vectorial_Variable> vectorial_variables(0);

         this->dirname=  "data";
         this->prefix= "erai.6h.";
         this->postfix=".nc";
         this->reference_date="2008-01-01";

         this->variables= variables;
         this->vectorial_variables= vectorial_variables;
         this->target_size= target_size;
         
         this->reloaded= false;

         this->nb_timestep_day= 4;
         this->time= time;
     }
     else if (strcmp (DatasetName, "ERAi_nodes") == 0)
     { 
  	// Definition of dimensions
     dimension_x={
          name:"lon",
          cyclic:true
  	};

     dimension_y={
          name:"lat",
          cyclic:false
  	};
 
     dimension_time={
          name:"time", // "Time"
          cyclic:false};
     
     // Definition of the grid
      std::vector<Dimension> dimensions_lon(1);
      dimensions_lon[0] = this->dimension_x;

      std::vector<Dimension> dimensions_lat(1);
      dimensions_lat[0] = this->dimension_y;

      Variable latitude={
          name: "lat",
          dimensions: dimensions_lat,
          a: 1.,
          b: 0.,
          Units: "degree_north",
          data2: data2_tmp};

      Variable longitude={
          name: "lon",
          dimensions: dimensions_lon,
          a: 1.,
          b: 0.,
          Units: "degree_east",
          data2: data2_tmp};

  	Grid M_grid={
  		interpolation_method: InterpolationType::FromGridToMesh,
  	    //interp_type : TriangleInterpEnum, // slower
  	    interp_type : BilinearInterpEnum,
  	    //interp_type : NearestInterpEnum,
  		this->dirname="data",
  		filename:"erai.6h.200803.nc",

  		latitude: latitude,
  		longitude: longitude,

  		dimension_x: dimension_x,
  		dimension_y: dimension_y,

  		mpp_file: "",
  		interpolation_in_latlon: true,

          loaded: false,

  		masking: false
    };
     
     this->grid= M_grid;
     
     // Definition of the data
     
     std::vector<Dimension> dimensions_time(1);
     dimensions_time[0] = this->dimension_time;

     std::vector<Dimension> dimensions(3);
     dimensions[0] = this->dimension_time;
     dimensions[1] = this->dimension_y;
     dimensions[2] = this->dimension_x;
     
     Variable time={
         name: "time",
         dimensions: dimensions_time,
         a: 1.,
         b: 0.,
         Units: "hours",
         data2: data2_tmp
     };
     
      // conversion factors: xnew = a*x + b
      Variable u={
          name: "10U", // U10M
          dimensions: dimensions,
          a: 1.,
          b: 0.,
          Units: "m/s",
          data2: data2_tmp
      };

      Variable v={
          name: "10V", // U10M
          dimensions: dimensions,
          a: 1.,
          b: 0.,
          Units: "m/s",
          data2: data2_tmp
      };

      std::vector<Variable> variables(2);
      variables[0] = u;
      variables[1] = v;

      std::vector<int> uv_tmp(2);
          uv_tmp[0] = 0;
          uv_tmp[1] = 1;
 
      Vectorial_Variable uv={
          components_Id: uv_tmp,
          east_west_oriented: true
      };

      std::vector<Vectorial_Variable> vectorial_variables(1);
      vectorial_variables[0] = uv;

      this->dirname=  "data";
      this->prefix= "erai.6h.";
      this->postfix=".nc";
      this->reference_date="2008-01-01";

      this->variables= variables;
      this->vectorial_variables= vectorial_variables;
      this->target_size= target_size;
      
      this->reloaded= false;

      this->nb_timestep_day= 4;
      this->time= time;
      }
     else
       {
   	fprintf (stderr, "Dataset: unknown projection %s\n",DatasetName);
   	fprintf (stderr, "valid types are:\n");
   	fprintf (stderr, "asr_nodes\n");
   	fprintf (stderr, "asr_elements\n");
   	fprintf (stderr, "topaz_nodes\n");
   	fprintf (stderr, "topaz_elements\n");
   	fprintf (stderr, "ice_topaz_elements\n");
   	fprintf (stderr, "etopo_elements\n");
   	fprintf (stderr, "ERAi_nodes\n");
    fprintf (stderr, "ERAi_elements\n");
       //close_Dataset (this);
     }
     
     this->ftime_range.resize(2,0.);

   }

} // Nextsim
