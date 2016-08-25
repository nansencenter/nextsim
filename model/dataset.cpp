/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   dataset.cpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#include <dataset.hpp>
#include <date.hpp>
extern "C"
{
#include <mapx.h>
}
     
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

   DataSet::DataSet(char const *DatasetName, int target_size_tmp)
   {
//     Dataset *this;

       name = std::string(DatasetName);

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
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp
     	};

         Variable longitude={
             name: "XLONG",
             dimensions: dimensions_latlon,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp
     	};

         Variable time_tmp={
             name: "time",
             dimensions: dimensions_time,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp
     	};

         // conversion factors: xnew = a*x + b
         Variable u={
             name: "U10", // U10M
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v={
             name: "V10", // U10M
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

     	Grid grid_tmp={
     		interpolation_method: InterpolationType::FromGridToMesh,
     	    //interp_type : TriangleInterpEnum,  // slower
     	    interp_type : BilinearInterpEnum,
     	    //interp_type : NearestInterpEnum,

     		dirname="data",
     		//filename:"asr30km.comb.2d.200803.nc",
            prefix= "asr30km.comb.2d.", // "asr30km.comb.2D.";
            postfix=".nc",

     		latitude: latitude,
     		longitude: longitude,

     		dimension_x: dimension_x,
     		dimension_y: dimension_y,

     		mpp_file: "NpsASR.mpp",
      		interpolation_in_latlon: false,

            loaded: false,
            monthly_dataset:true,

     		masking: false
     	};

         std::vector<Variable> variables_tmp(2);
             variables_tmp[0] = u;
             variables_tmp[1] = v;

         std::vector<int> uv_tmp(2);
             uv_tmp[0] = 0;
             uv_tmp[1] = 1;

         Vectorial_Variable uv={
             components_Id: uv_tmp,
             east_west_oriented: false // if false, then we assume it is oriented following the input grid
         };

         std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
         vectorial_variables_tmp[0] = uv;

             dirname= "data";
             prefix= "asr30km.comb.2d."; // "asr30km.comb.2D.";
             postfix=".nc";
             reference_date= "1901-01-01";

             variables= variables_tmp;
             vectorial_variables= vectorial_variables_tmp;
             target_size= target_size_tmp;
             grid= grid_tmp;

             reloaded=false;

             nb_timestep_day= 8;
             daily_mean=false;
             time= time_tmp;
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
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp
     	};

         Variable longitude={
             name: "XLONG",
             dimensions: dimensions_latlon,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp
     	};

         Variable time_tmp={
             name: "time",
             dimensions: dimensions_time,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp
     	};

         // conversion factors: xnew = a*x + b
         Variable u={
             name: "U10", // U10M
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v={
             name: "V10", // U10M
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

     	Grid grid_tmp={
     		interpolation_method: InterpolationType::FromGridToMesh,
     	    //interp_type : TriangleInterpEnum,  // slower
     	    interp_type : BilinearInterpEnum,
     	    //interp_type : NearestInterpEnum,

     		dirname="data",
     		//filename:"asr30km.comb.2d.200803.nc",
            prefix="asr30km.comb.2d.", // "asr30km.comb.2D.";
            postfix=".nc",

     		latitude: latitude,
     		longitude: longitude,

     		dimension_x: dimension_x,
     		dimension_y: dimension_y,

     		mpp_file: "NpsASR.mpp",
      		interpolation_in_latlon: false,

            loaded: false,
            monthly_dataset:true,

     		masking: false
     	};

         nb_timestep_day=8;

         Variable tair={
             name:"T2",
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a:1.,
             b:-273.15,
             Units:"C",
             data2: data2_tmp
     	}; // T2M
         Variable mixrat={
             name:"Q2",
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a:1.,
             b:0.,
             Units:"",
             data2: data2_tmp
     	}; // Q2M
         Variable mslp={
             name:"SLP",
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a:1e2,
             b:0.,
             Units:"Pa",
             data2: data2_tmp
     	}; //PSFC, a=1.
         Variable Qsw_in={
             name:"SWDNB",
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a:1.,
             b:0.,
             Units:"W/m^2",
             data2: data2_tmp
     	};
         Variable Qlw_in={
             name:"LWDNB",
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a:1.,
             b:0.,
             Units:"W/m^2",
             data2: data2_tmp
     	};
         Variable snowfr={
             name:"SR",
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a:1.,
             b:0.,
             Units:"",
             data2: data2_tmp
     	};
         Variable precip={
             name:"RAINNC",
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a:physical::rhow/1000.*(nb_timestep_day)/(24.*3600),
             b:0.,
             Units:"kg/m^2/s",
             data2: data2_tmp
     	};

         std::vector<Variable> variables_tmp(7);
         variables_tmp[0] = tair;
         variables_tmp[1] = mixrat;
         variables_tmp[2] = mslp;
         variables_tmp[3] = Qsw_in;
         variables_tmp[4] = Qlw_in;
         variables_tmp[5] = snowfr;
         variables_tmp[6] = precip;

         std::vector<Vectorial_Variable> vectorial_variables_tmp(0);
         
         dirname="data";
         prefix="asr30km.comb.2d."; // "asr30km.comb.2D.";
         postfix=".nc";
         reference_date= "1901-01-01";

         variables= variables_tmp;
         vectorial_variables= vectorial_variables_tmp;
         target_size=target_size_tmp;
         grid= grid_tmp;

         reloaded=false;

         daily_mean=false;
         time= time_tmp;
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
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp};

         Variable longitude={
             name: "longitude",
             dimensions: dimensions_latlon,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp};

         Variable time_tmp={
             name: "time",
             dimensions: dimensions_time,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp};

         Variable u={
             name: "u",
             dimensions: dimensions_uv,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v={
             name: "v",
             dimensions: dimensions_uv,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

     	Variable ssh={
     		name: "ssh",
     		dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
     		a: 1.,
     		b: 0.,
     		Units: "m/s",
     		data2: data2_tmp
     	};

         Grid grid_tmp={
             interpolation_method: InterpolationType::FromMeshToMesh2dx,
     		interp_type: -1,
             dirname= "data",
             //filename: "TP4DAILY_200803_3m.nc",
             prefix= "TP4DAILY_",
             postfix= "_30m.nc",

             latitude: latitude,
             longitude: longitude,

             dimension_x: dimension_x,
             dimension_y: dimension_y,

             mpp_file: "NpsNextsim.mpp",
     		interpolation_in_latlon: false,

             loaded: false,
             monthly_dataset:true,

     		masking: true,
     		masking_variable: ssh
     	};

         std::vector<Variable> variables_tmp(3);
         variables_tmp[0] = u;
         variables_tmp[1] = v;
         variables_tmp[2] = ssh;

         std::vector<int> uv_tmp(2);
             uv_tmp[0] = 0;
             uv_tmp[1] = 1;

         Vectorial_Variable uv={
             components_Id: uv_tmp,
             east_west_oriented: false // if false, then we assume it is oriented following the mpp_file defined for the grid
     	};

         std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
         vectorial_variables_tmp[0] = uv;

         dirname= "data";
         prefix= "TP4DAILY_";
         postfix= "_30m.nc";
         reference_date= "1950-01-01";

         variables= variables_tmp;
         vectorial_variables= vectorial_variables_tmp;
         target_size= target_size_tmp;
         grid= grid_tmp;

         reloaded=false;

         nb_timestep_day= 1;
         daily_mean=true;
         time= time_tmp;
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
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp};

         Variable longitude={
             name: "longitude",
             dimensions: dimensions_latlon,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp};

         Variable time_tmp={
             name: "time",
             dimensions: dimensions_time,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp};

     	Variable sst={
     		name: "temperature",
     		dimensions: dimensions_uv,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
     		a: 1.,
     		b: 0.,
     		Units: "deg celsius",
     		data2: data2_tmp
     	};

     	Variable sss={
     		name: "salinity",
     		dimensions: dimensions_uv,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable mld={
     		name: "mlp",
     		dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

         Grid grid_tmp={
             interpolation_method: InterpolationType::FromMeshToMesh2dx,
     		interp_type: -1,
             dirname= "data",
             //filename: "TP4DAILY_200803_3m.nc",
             prefix= "TP4DAILY_",
             postfix= "_3m.nc",

             latitude: latitude,
             longitude: longitude,

             dimension_x: dimension_x,
             dimension_y: dimension_y,

             mpp_file: "NpsNextsim.mpp",
     		interpolation_in_latlon: false,

             loaded: false,
             monthly_dataset:true,

     		masking: true,
     		masking_variable: sss
     	};

         std::vector<Variable> variables_tmp(3);
         variables_tmp[0] = sst;
         variables_tmp[1] = sss;
         variables_tmp[2] = mld;

         std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

             dirname= "data";
             prefix= "TP4DAILY_";
             postfix= "_3m.nc";
             reference_date= "1950-01-01";

             variables= variables_tmp;
             vectorial_variables= vectorial_variables_tmp;
             target_size= target_size_tmp;
             grid= grid_tmp;

             reloaded=false;

             nb_timestep_day= 1;
             daily_mean=true;
             time= time_tmp;
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
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp};

         Variable longitude={
             name: "longitude",
             dimensions: dimensions_latlon,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp};

         Variable time_tmp={
             name: "time",
             dimensions: dimensions_time,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp};

     	Variable conc={
     		name: "fice",
     		dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable thick={
     		name: "hice",
     		dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable snow_thick={
     		name: "hsnow",
     		dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

         Grid grid_tmp={
             interpolation_method: InterpolationType::FromMeshToMesh2dx,
     		interp_type: -1,
             dirname= "data",
             //filename: "TP4DAILY_200803_3m.nc",
             prefix= "TP4DAILY_",
             postfix= "_3m.nc",

             latitude: latitude,
             longitude: longitude,

             dimension_x: dimension_x,
             dimension_y: dimension_y,

             mpp_file: "NpsNextsim.mpp",
     		interpolation_in_latlon: false,

             loaded: false,
             monthly_dataset:true,

     		masking: true,
     		masking_variable: conc
     	};

         std::vector<Variable> variables_tmp(3);
         variables_tmp[0] = conc;
         variables_tmp[1] = thick;
         variables_tmp[2] = snow_thick;

         std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

             dirname= "data";
             prefix= "TP4DAILY_";
             postfix= "_3m.nc";
             reference_date= "1950-01-01";

             variables= variables_tmp;
             vectorial_variables= vectorial_variables_tmp;
             target_size= target_size_tmp;
             grid= grid_tmp;

             reloaded=false;

     		nb_timestep_day= 1;
            daily_mean=true;
     	    time= time_tmp;
     }
     else if (strcmp (DatasetName, "asi_elements") == 0)
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
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp};

         Variable longitude={
             name: "longitude",
             dimensions: dimensions_latlon,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp};

         Variable time_tmp={
             name: "time",
             dimensions: dimensions_time,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1./3600,
             b: 0.,
             Units: "hours",
             data2: data2_tmp};

         Variable conc={
     		name: "icecon",
     		dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
     		a: 0.01,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

         Grid grid_tmp={
             interpolation_method: InterpolationType::FromMeshToMesh2dx,
     		interp_type: -1,
             dirname= "data",
             //filename: "TP4DAILY_200803_3m.nc",
             prefix= "asi-n6250-",
             postfix= "-v5i.nc",

             latitude: latitude,
             longitude: longitude,

             dimension_x: dimension_x,
             dimension_y: dimension_y,

             mpp_file: "NpsNextsim.mpp",
     		interpolation_in_latlon: false,

             loaded: false,
             monthly_dataset:false,

     		masking: true,
     		masking_variable: conc
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = conc;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        dirname= "data";
        prefix= "asi-n6250-";
        postfix= "-v5i.nc";
        reference_date= "2002-01-01";

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        target_size= target_size_tmp;
        grid= grid_tmp;

        reloaded=false;

        nb_timestep_day= 1;
        daily_mean=true;
     	time= time_tmp;
     }
     else if (strcmp (DatasetName, "arc_elements") == 0)
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
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp};

         Variable longitude={
             name: "longitude",
             dimensions: dimensions_latlon,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp};

         Variable time_tmp={
             name: "time",
             dimensions: dimensions_time,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 24.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp};

         Variable conc={
     		name: "sea_ice_concentration",
     		dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: 12500.,
            NaN_mask_defined: true,
            NaN_mask_value: 11500.,
     		a: 0.01,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

         Grid grid_tmp={
             interpolation_method: InterpolationType::FromMeshToMesh2dx,
     		interp_type: -1,
             dirname= "data",
             //filename: "TP4DAILY_200803_3m.nc",
             prefix= "LongitudeLatitudeGrid_3.125km_Arctic.nc",
             postfix= "",

             latitude: latitude,
             longitude: longitude,

             dimension_x: dimension_x,
             dimension_y: dimension_y,

             mpp_file: "NpsNextsim.mpp",
     		interpolation_in_latlon: false,

             loaded: false,
             monthly_dataset:false,

     		masking: true,
     		masking_variable: conc
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = conc;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        dirname= "data";
        prefix= "Arc_";
        postfix= "_res3.125_pyres.nc";
        reference_date= "0001-01-01";

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        target_size= target_size_tmp;
        grid= grid_tmp;

        reloaded=false;

        nb_timestep_day= 1;
        daily_mean=true;
     	time= time_tmp;
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
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp
     	};

         Variable longitude={
             name: "x",
             dimensions: dimensions_lon,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp
     	};

     	Grid grid_tmp={
     		interpolation_method: InterpolationType::FromGridToMesh,
     	    //interp_type : TriangleInterpEnum, // slower
     	    interp_type : BilinearInterpEnum,
     	    //interp_type : NearestInterpEnum,
     		dirname="data",
     		//filename:"ETOPO1_Ice_g_gmt4.grd",
            prefix="ETOPO1_Ice_g_gmt4.grd",
            postfix="",

     		latitude: latitude,
     		longitude: longitude,

     		dimension_x: dimension_x,
     		dimension_y: dimension_y,

     		mpp_file: "",
     		interpolation_in_latlon: true,

             loaded: false,
             monthly_dataset:false,

     		masking: false
     	};

         Variable z={
             name:"z",
             dimensions: dimensions,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a:-1.,
             b:0.,
             Units:"m",
             data2: data2_tmp
     	};

         std::vector<Variable> variables_tmp(1);
         variables_tmp[0] = z;

         std::vector<Vectorial_Variable> vectorial_variables_tmp(0);


         dirname="data";
         prefix="ETOPO1_Ice_g_gmt4";
         postfix=".grd";
         reference_date= "";
         
         variables= variables_tmp;
         vectorial_variables= vectorial_variables_tmp;
         target_size=target_size_tmp;
         
         grid= grid_tmp;
         reloaded=false;
         
         nb_timestep_day= 0;
         daily_mean=false;
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
         dimensions_lon[0] = dimension_x;

         std::vector<Dimension> dimensions_lat(1);
         dimensions_lat[0] = dimension_y;

         Variable latitude={
             name: "lat",
             dimensions: dimensions_lat,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp
     	};

         Variable longitude={
             name: "lon",
             dimensions: dimensions_lon,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp
     	};

     	Grid grid_tmp={
     		interpolation_method: InterpolationType::FromGridToMesh,
     	    //interp_type : TriangleInterpEnum, // slower
     	    interp_type : BilinearInterpEnum,
     	    //interp_type : NearestInterpEnum,
     		dirname="data",
     		//filename:"erai.6h.201304.nc",//"erai.6h.200803.nc",
            prefix= "erai.6h.",
            postfix=".nc",

     		latitude: latitude,
     		longitude: longitude,

     		dimension_x: dimension_x,
     		dimension_y: dimension_y,

     		mpp_file: "",
     		interpolation_in_latlon: true,

             loaded: false,
             monthly_dataset:true,

     		masking: false
     	};

        grid= grid_tmp;

        // Definition of the data

        std::vector<Dimension> dimensions_time(1);
        dimensions_time[0] = dimension_time;

        std::vector<Dimension> dimensions(3);
        dimensions[0] = dimension_time;
        dimensions[1] = dimension_y;
        dimensions[2] = dimension_x;

        Variable time_tmp={
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            a: 1.,
            b: 0.,
            Units: "hours",
            data2: data2_tmp
    	};

        nb_timestep_day=4;

        Variable tair={
            name:"2T",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            a:1.,
            b:-273.15,
            Units:"C",
            data2: data2_tmp
    	}; // T2M
        Variable dair={
            name:"2D",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            a:1.,
            b:-273.15,
            Units:"C",
            data2: data2_tmp
    	}; // Q2M
        Variable mslp={
            name:"MSL",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            a:1.,
            b:0.,
            Units:"Pa",
            data2: data2_tmp
    	}; //PSFC, a=1.
        Variable Qsw_in={
            name:"SSRD",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            a:nb_timestep_day/(24.*3600),
            b:0.,
            Units:"W/m^2",
            data2: data2_tmp
    	};
        Variable tcc={
            name:"TCC",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            a:1.,
            b:0.,
            Units:"",
            data2: data2_tmp
    	};

        Variable precip={
            name:"TP",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            a:physical::rhow*(nb_timestep_day)/(24.*3600),
            b:0.,
            Units:"kg/m^2/s",
            data2: data2_tmp
    	};

        std::vector<Variable> variables_tmp(6);
        variables_tmp[0] = tair;
        variables_tmp[1] = dair;
        variables_tmp[2] = mslp;
        variables_tmp[3] = Qsw_in;
        variables_tmp[4] = tcc;
        variables_tmp[5] = precip;

         std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

         dirname=  "data";
         prefix= "erai.6h.";
         postfix=".nc";
         reference_date="2013-01-01";//"2008-01-01";

         variables= variables_tmp;
         vectorial_variables= vectorial_variables_tmp;
         target_size= target_size_tmp;

         reloaded= false;

         daily_mean=false;
         
         time= time_tmp;
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
      dimensions_lon[0] = dimension_x;

      std::vector<Dimension> dimensions_lat(1);
      dimensions_lat[0] = dimension_y;

      Variable latitude={
          name: "lat",
          dimensions: dimensions_lat,
          land_mask_defined: false,
          land_mask_value: 0.,
          NaN_mask_defined: false,
          NaN_mask_value: 0.,
          a: 1.,
          b: 0.,
          Units: "degree_north",
          data2: data2_tmp};

      Variable longitude={
          name: "lon",
          dimensions: dimensions_lon,
          land_mask_defined: false,
          land_mask_value: 0.,
          NaN_mask_defined: false,
          NaN_mask_value: 0.,
          a: 1.,
          b: 0.,
          Units: "degree_east",
          data2: data2_tmp};

  	Grid grid_tmp={
  		interpolation_method: InterpolationType::FromGridToMesh,
  	    //interp_type : TriangleInterpEnum, // slower
  	    interp_type : BilinearInterpEnum,
  	    //interp_type : NearestInterpEnum,
  		dirname="data",
  		//filename:"erai.6h.201304.nc",//"erai.6h.200803.nc",
        prefix= "erai.6h.",
        postfix=".nc",

  		latitude: latitude,
  		longitude: longitude,

  		dimension_x: dimension_x,
  		dimension_y: dimension_y,

  		mpp_file: "",
  		interpolation_in_latlon: true,

        loaded: false,
        
        monthly_dataset:true,
        
        masking: false
    };

     grid= grid_tmp;

     // Definition of the data

     std::vector<Dimension> dimensions_time(1);
     dimensions_time[0] = dimension_time;

     std::vector<Dimension> dimensions(3);
     dimensions[0] = dimension_time;
     dimensions[1] = dimension_y;
     dimensions[2] = dimension_x;

     Variable time_tmp={
         name: "time",
         dimensions: dimensions_time,
         land_mask_defined: false,
         land_mask_value: 0.,
         NaN_mask_defined: false,
         NaN_mask_value: 0.,
         a: 1.,
         b: 0.,
         Units: "hours",
         data2: data2_tmp
     };

      // conversion factors: xnew = a*x + b
      Variable u={
          name: "10U", // U10M
          dimensions: dimensions,
          land_mask_defined: false,
          land_mask_value: 0.,
          NaN_mask_defined: false,
          NaN_mask_value: 0.,
          a: 1.,
          b: 0.,
          Units: "m/s",
          data2: data2_tmp
      };

      Variable v={
          name: "10V", // U10M
          dimensions: dimensions,
          land_mask_defined: false,
          land_mask_value: 0.,
          NaN_mask_defined: false,
          NaN_mask_value: 0.,
          a: 1.,
          b: 0.,
          Units: "m/s",
          data2: data2_tmp
      };

      std::vector<Variable> variables_tmp(2);
      variables_tmp[0] = u;
      variables_tmp[1] = v;

      std::vector<int> uv_tmp(2);
          uv_tmp[0] = 0;
          uv_tmp[1] = 1;

      Vectorial_Variable uv={
          components_Id: uv_tmp,
          east_west_oriented: true
      };

      std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
      vectorial_variables_tmp[0] = uv;

      dirname=  "data";
      prefix= "erai.6h.";
      postfix=".nc";
      reference_date="2013-01-01";//"2008-01-01";

      variables= variables_tmp;
      vectorial_variables= vectorial_variables_tmp;
      target_size= target_size_tmp;

      reloaded= false;

      nb_timestep_day= 4;
      daily_mean=false;
      time= time_tmp;
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
    fprintf (stderr, "asi_elements\n");
    fprintf (stderr, "arc_elements\n");
   	fprintf (stderr, "etopo_elements\n");
   	fprintf (stderr, "ERAi_nodes\n");
    fprintf (stderr, "ERAi_elements\n");
       //close_Dataset (this);
     }

     ftime_range.resize(2,0.);

   }

void
DataSet::loadGrid(Grid *grid_ptr, int current_time)
{
    // Attributes (scaling and offset)
    netCDF::NcVarAtt att;
    double scale_factor;
    double add_offset;
    
    //std::cout<<"---------------------fist loading ...\n";
    std::string current_timestr;
    if ( current_time > 0 )
    {
        if(grid_ptr->monthly_dataset)
            current_timestr = to_date_string_ym(current_time);
        else
            current_timestr = to_date_string_yd(current_time);
        
        std::cout <<"TIMESTR= "<< current_timestr <<"\n";
    }
    else
        current_timestr = "";


    // the ETOPO forcing filename does not contain the information on year and month
    if (grid_ptr->postfix=="")
    {
        current_timestr = "";
    }

    std::string filename = (boost::format( "%1%/%2%/%3%%4%%5%" )
                            % Environment::simdataDir().string()
                            % grid_ptr->dirname
                            % grid_ptr->prefix
                            % current_timestr
                            % grid_ptr->postfix
                            ).str();

    std::cout<<"GRID : FILENAME = "<< filename <<"\n";

	std::cout <<"GRID : READ NETCDF starts\n";
    if ( ! boost::filesystem::exists(filename) )
        throw std::runtime_error("File not found: " + filename);

	netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
    netCDF::NcDim tmpDim;

    tmpDim = dataFile.getDim(grid_ptr->dimension_y.name);
	grid_ptr->M  =  tmpDim.getSize();

    tmpDim = dataFile.getDim(grid_ptr->dimension_x.name);
	grid_ptr->N  =  tmpDim.getSize();

    //switch (grid_ptr->latitude.dimensions.size())
    //{
    // Here only two cases are considered, either the
    //    case 1:
	if(grid_ptr->latitude.dimensions.size()==1)
	{
		// read in coordinates
		std::vector<size_t> index_x_count(1);
		std::vector<size_t> index_y_count(1);

		std::vector<size_t> index_x_start(1);
		std::vector<size_t> index_y_start(1);

		index_y_start[0] = 0;
		index_y_count[0] = grid_ptr->M;

		index_x_start[0] = 0;
		index_x_count[0] = grid_ptr->N;

		std::vector<double> LAT(index_y_count[0]);
		std::vector<double> LON(index_x_count[0]);

		netCDF::NcVar VLAT = dataFile.getVar(grid_ptr->latitude.name);
		netCDF::NcVar VLON = dataFile.getVar(grid_ptr->longitude.name);
		std::cout <<"GRID : READ NETCDF done\n";

		VLAT.getVar(index_y_start,index_y_count,&LAT[0]);
		VLON.getVar(index_x_start,index_x_count,&LON[0]);

        // Need to multiply with scale factor and add offset - these are stored as variable attributes
        scale_factor=1.;
        try
        {
            att = VLAT.getAtt("scale_factor");
            att.getValues(&scale_factor);
        }
        catch(netCDF::exceptions::NcException& e)
        {}

        add_offset=0.;
        try
        {
            att = VLAT.getAtt("add_offset");
            att.getValues(&add_offset);
        }
        catch(netCDF::exceptions::NcException& e)
        {}

        for (int i=0; i<(index_y_count[0]); ++i) 
            LAT[i]=LAT[i]*scale_factor + add_offset;
                
        // Need to multiply with scale factor and add offset - these are stored as variable attributes
        scale_factor=1.;
        try
        {
            att = VLON.getAtt("scale_factor");
            att.getValues(&scale_factor);
        }
        catch(netCDF::exceptions::NcException& e)
        {}

        add_offset=0.;
        try
        {
            att = VLON.getAtt("add_offset");
            att.getValues(&add_offset);
        }
        catch(netCDF::exceptions::NcException& e)
        {}

        for (int i=0; i<(index_x_count[0]); ++i) 
            LON[i]=LON[i]*scale_factor + add_offset;
                

		grid_ptr->gridY=LAT;
		grid_ptr->gridX=LON;

        // Save lon and lat for possible output
        grid_ptr->gridLAT=LAT;
        grid_ptr->gridLON=LON;
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

		index_py_start[0] = 0;
		index_py_start[1] = 0;

		index_py_count[0] = grid_ptr->M;
		index_py_count[1] = grid_ptr->N;

		index_px_start[0] = 0;
		index_px_start[1] = 0;

		index_px_count[0] = grid_ptr->M;
		index_px_count[1] = grid_ptr->N;

		if(grid_ptr->interpolation_method==InterpolationType::FromGridToMesh)
		{
            // We the initial grid is actually regular, we can still use FromGridToMesh
            // by only taking the first line and column into account (only used for ASR so far)
			index_py_count[1] = 1;
			index_px_count[0] = 1;
		}

		std::vector<double> XLAT(index_px_count[0]*index_px_count[1]);
		std::vector<double> XLON(index_px_count[0]*index_px_count[1]);
		std::vector<double> YLAT(index_py_count[0]*index_py_count[1]);
		std::vector<double> YLON(index_py_count[0]*index_py_count[1]);

		netCDF::NcVar VLAT = dataFile.getVar(grid_ptr->latitude.name);
		netCDF::NcVar VLON = dataFile.getVar(grid_ptr->longitude.name);
		std::cout <<"GRID : READ NETCDF done\n";

        // Need to multiply with scale factor and add offset - these are stored as variable attributes
		VLAT.getVar(index_px_start,index_px_count,&XLAT[0]);
		VLON.getVar(index_px_start,index_px_count,&XLON[0]);

		VLAT.getVar(index_py_start,index_py_count,&YLAT[0]);
		VLON.getVar(index_py_start,index_py_count,&YLON[0]);

        // Apply the scale factor and offset if any
        scale_factor=1.;
        try
        {
            att = VLAT.getAtt("scale_factor");
            att.getValues(&scale_factor);
        }
        catch(netCDF::exceptions::NcException& e)
        {}

        add_offset=0.;
        try
        {
            att = VLAT.getAtt("add_offset");
            att.getValues(&add_offset);
        }
        catch(netCDF::exceptions::NcException& e)
        {}
        
        if(add_offset!=0. || scale_factor!=1.)
        {    
            for (int i=0; i<(index_px_count[0]*index_px_count[1]); ++i) 
            {
                XLON[i]=XLON[i]*scale_factor + add_offset;
                XLAT[i]=XLAT[i]*scale_factor + add_offset;
            }
        
            for (int i=0; i<(index_py_count[0]*index_py_count[1]); ++i) 
            {
                YLON[i]=YLON[i]*scale_factor + add_offset;
                YLAT[i]=YLAT[i]*scale_factor + add_offset;
            }
        }
        
        // projection

		std::vector<double> X(index_px_count[0]*index_px_count[1]);
		std::vector<double> Y(index_py_count[0]*index_py_count[1]);

		mapx_class *map;
		std::string configfile = (boost::format( "%1%/%2%/%3%" )
                                  % Environment::nextsimDir().string()
                                  % grid_ptr->dirname
                                  % grid_ptr->mpp_file
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

		if(grid_ptr->interpolation_method==InterpolationType::FromMeshToMesh2dx)
		{
			if(grid_ptr->masking){
				netCDF::NcVar VMASK;
                netCDF::NcDim tmpDim;

                // Open the datafile
				//VMASK = dataFile.getVar(grid_ptr->masking_variable.name);
                if ( current_time > 0 )
                {
                    if(grid_ptr->monthly_dataset)
                        current_timestr = to_date_string_ym(current_time);
                    else
                        current_timestr = to_date_string_yd(current_time);
        
                    std::cout <<"TIMESTR= "<< current_timestr <<"\n";
                }
                else
                    current_timestr = "";

                filename = (boost::format( "%1%/%2%/%3%%4%%5%" )
                                        % Environment::simdataDir().string()
                                        % dirname
                                        % prefix
                                        % current_timestr
                                        % postfix
                                        ).str();

                std::cout<<"GRID : FILENAME = "<< filename <<"\n";

            	std::cout <<"GRID : READ NETCDF starts\n";
                if ( ! boost::filesystem::exists(filename) )
                    throw std::runtime_error("File not found: " + filename);

            	netCDF::NcFile dataFile2(filename, netCDF::NcFile::read);

                // load the data
                VMASK = dataFile2.getVar(grid_ptr->masking_variable.name);
                
				std::vector<double> data_in;

				std::vector<double> reduced_FX;
				std::vector<double> reduced_FY;
				std::vector<int> reduced_nodes_ind;

				std::vector<size_t> index_start(3,0);
				std::vector<size_t> index_count(3);

				index_start.resize(grid_ptr->masking_variable.dimensions.size());
				index_count.resize(grid_ptr->masking_variable.dimensions.size());

				for(int k=0; k<grid_ptr->masking_variable.dimensions.size(); ++k)
				{
                    tmpDim = dataFile2.getDim(grid_ptr->masking_variable.dimensions[k].name);
                    index_start[k] = 0;
					index_count[k] = tmpDim.getSize();;
				}
				index_start[0] = 0;
				index_count[0] = 1;

				if((index_px_count[0]!=index_count[grid_ptr->masking_variable.dimensions.size()-2]) || (index_px_count[1]!=index_count[grid_ptr->masking_variable.dimensions.size()-1]))
				{
                    std::cout << "index_px_count[0] = " << index_px_count[0] << " index_count[grid_ptr->masking_variable.dimensions.size()-2] = " << index_count[grid_ptr->masking_variable.dimensions.size()-2] <<"\n";
					std::cout << "index_px_count[1] = " << index_px_count[1] << " index_count[grid_ptr->masking_variable.dimensions.size()-1] = " << index_count[grid_ptr->masking_variable.dimensions.size()-1] <<"\n";
                    throw std::logic_error("Not the same dimension for the masking variable and the grid!!");
				}

				data_in.resize(index_px_count[0]*index_px_count[1]);
				VMASK.getVar(index_start,index_count,&data_in[0]);

                // Read the attributes
				netCDF::NcVarAtt att;
				
                // Look for FillValue definition
                int FillValue;
                bool find_FillValue=true;
                try
                {
    				att = VMASK.getAtt("_FillValue");
    				att.getValues(&FillValue);
                }
                catch(netCDF::exceptions::NcException& e)
                {
                    find_FillValue=false;
                }
                
                // Look for missing_value definition
                int missing_value;
                bool find_missing_value=true;
                try
                {
    				att = VMASK.getAtt("missing_value");
    				att.getValues(&missing_value);
                }
                catch(netCDF::exceptions::NcException& e)
                {
                    find_missing_value=false;
                }
                
                bool find_land_mask     =grid_ptr->masking_variable.land_mask_defined;
                double land_mask_value  =grid_ptr->masking_variable.land_mask_value;
                bool find_NaN_mask      =grid_ptr->masking_variable.NaN_mask_defined;
                double NaN_mask_value   =grid_ptr->masking_variable.NaN_mask_value;

                double tmp_data;
				for (int i=0; i<index_px_count[0]; ++i)
				{
					for (int j=0; j<index_px_count[1]; ++j)
					{
                        tmp_data=data_in[index_px_count[1]*i+j];
						if (    (!find_FillValue        || (find_FillValue      && (tmp_data != FillValue))         )&&
                                (!find_missing_value    || (find_missing_value  && (tmp_data != missing_value))     )&&
                                (!find_land_mask        || (find_land_mask      && (tmp_data != land_mask_value))   )&&
                                (!find_NaN_mask         || (find_NaN_mask       && (tmp_data != NaN_mask_value))    )   )
						{
							reduced_FX.push_back(X[index_px_count[1]*i+j]);
							reduced_FY.push_back(Y[index_px_count[1]*i+j]);
							reduced_nodes_ind.push_back(index_px_count[1]*i+j);
						}
					}
				}
				grid_ptr->gridX=reduced_FX;
				grid_ptr->gridY=reduced_FY;
				grid_ptr->reduced_nodes_ind=reduced_nodes_ind;
			}
			else // no masking of the Filled Value
			{
				grid_ptr->gridX=X;
				grid_ptr->gridY=Y;
			}

			std::cout <<"GRID : Triangulate starts\n";
			BamgTriangulatex(&grid_ptr->pfindex,&grid_ptr->pfnels,&grid_ptr->gridX[0],&grid_ptr->gridY[0],grid_ptr->gridX.size());
			std::cout <<"GRID : NUMTRIANGLES= "<< grid_ptr->pfnels <<"\n";
			std::cout <<"GRID : Triangulate done\n";
		}
		else
		{
			grid_ptr->gridX=X;
			grid_ptr->gridY=Y;
		}

        // Save lon and lat for possible output
        grid_ptr->gridLAT=XLAT;
        grid_ptr->gridLON=YLON;

	//	break;
	//
    //default:
    //   std::cout << "invalid ocean initialisation"<<"\n";
    //    throw std::logic_error("invalid ocean forcing");
	}

    grid_ptr->loaded=true;
}

} // Nextsim
