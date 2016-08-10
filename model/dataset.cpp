/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   dataset.cpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#include <dataset.hpp>
     

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
#if 0
     if (strcmp (DatasetName, "asr_nodes") == 0)
     { 
     	// Definition of asr grid and datasets
         Dimension asr_dimension_x={
             name:"x",
             start:0,
             end:359,
             cyclic:false
     	};

         Dimension asr_dimension_y={
             name:"y",
             start:0,
             end:359,
             cyclic:false
     	};

         Dimension asr_dimension_time={
             name:"time", // "Time"
             start:0,
             end:247,
             cyclic:false
     	};

         std::vector<Dimension> dimensions_asr_latlon(2);
         dimensions_asr_latlon[0] = asr_dimension_y;
         dimensions_asr_latlon[1] = asr_dimension_x;

         std::vector<Dimension> dimensions_asr(3);
         dimensions_asr[0] = asr_dimension_time;
         dimensions_asr[1] = asr_dimension_y;
         dimensions_asr[2] = asr_dimension_x;

         std::vector<Dimension> dimensions_asr_time(1);
         dimensions_asr_time[0] = asr_dimension_time;

         Variable asr_latitude={
             name: "XLAT",
             dimensions: dimensions_asr_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp
     	};

         Variable asr_longitude={
             name: "XLONG",
             dimensions: dimensions_asr_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp
     	};

         Variable asr_time={
             name: "time",
             dimensions: dimensions_asr_time,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp
     	};

         // conversion factors: xnew = a*x + b
         Variable u10={
             name: "U10", // U10M
             dimensions: dimensions_asr,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v10={
             name: "V10", // U10M
             dimensions: dimensions_asr,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};
    
     	M_asr_grid={
     		interpolation_method: InterpolationType::FromGridToMesh,
     	    //interp_type : TriangleInterpEnum,  // slower
     	    interp_type : BilinearInterpEnum,
     	    //interp_type : NearestInterpEnum,

     		dirname:"data",
     		filename:"asr30km.comb.2d.200803.nc",

     		latitude: asr_latitude,
     		longitude: asr_longitude,

     		dimension_x: asr_dimension_x,
     		dimension_y: asr_dimension_y,

     		mpp_file: "NpsASR.mpp",
      		interpolation_in_latlon: false,

             loaded: false,

     		masking: false
     	};

         std::vector<Variable> variables_tmp0(2);
             variables_tmp0[0] = u10;
             variables_tmp0[1] = v10;
            
         std::vector<int> uv10_tmp0(2);
             uv10_tmp0[0] = 0;
             uv10_tmp0[1] = 1;
    
         Vectorial_Variable uv10={
             components_Id: uv10_tmp0,
             east_west_oriented: false // if false, then we assume it is oriented following the input grid
         };
            
         std::vector<Vectorial_Variable> vectorial_variables_tmp0(1);
             vectorial_variables_tmp0[0] = uv10;

         M_asr_nodes_dataset={
             dirname: "data",
             prefix: "asr30km.comb.2d.", // "asr30km.comb.2D.";
             postfix:".nc",
             reference_date: "1901-01-01",

             variables: variables_tmp0,
             vectorial_variables: vectorial_variables_tmp0,
             target_size: M_num_nodes,
             grid: &M_asr_grid,

             reloaded: false,

             nb_timestep_day: 8,
             time: asr_time,
             dimension_time: asr_dimension_time
     	};

         M_asr_nodes_dataset.ftime_range.resize(2,0.);
     }
     else if (strcmp (DatasetName, "asr_elements") == 0)
     { 
     	// Definition of asr grid and datasets
         Dimension asr_dimension_x={
             name:"x",
             start:0,
             end:359,
             cyclic:false
     	};

         Dimension asr_dimension_y={
             name:"y",
             start:0,
             end:359,
             cyclic:false
     	};

         Dimension asr_dimension_time={
             name:"time", // "Time"
             start:0,
             end:247,
             cyclic:false
     	};

         std::vector<Dimension> dimensions_asr_latlon(2);
         dimensions_asr_latlon[0] = asr_dimension_y;
         dimensions_asr_latlon[1] = asr_dimension_x;

         std::vector<Dimension> dimensions_asr(3);
         dimensions_asr[0] = asr_dimension_time;
         dimensions_asr[1] = asr_dimension_y;
         dimensions_asr[2] = asr_dimension_x;

         std::vector<Dimension> dimensions_asr_time(1);
         dimensions_asr_time[0] = asr_dimension_time;

         Variable asr_latitude={
             name: "XLAT",
             dimensions: dimensions_asr_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp
     	};

         Variable asr_longitude={
             name: "XLONG",
             dimensions: dimensions_asr_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp
     	};

         Variable asr_time={
             name: "time",
             dimensions: dimensions_asr_time,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp
     	};

         // conversion factors: xnew = a*x + b
         Variable u10={
             name: "U10", // U10M
             dimensions: dimensions_asr,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v10={
             name: "V10", // U10M
             dimensions: dimensions_asr,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};
    
     	M_asr_grid={
     		interpolation_method: InterpolationType::FromGridToMesh,
     	    //interp_type : TriangleInterpEnum,  // slower
     	    interp_type : BilinearInterpEnum,
     	    //interp_type : NearestInterpEnum,

     		dirname:"data",
     		filename:"asr30km.comb.2d.200803.nc",

     		latitude: asr_latitude,
     		longitude: asr_longitude,

     		dimension_x: asr_dimension_x,
     		dimension_y: asr_dimension_y,

     		mpp_file: "NpsASR.mpp",
      		interpolation_in_latlon: false,

             loaded: false,

     		masking: false
     	};

         int nb_timestep_day=8;

         Variable tair={
             name:"T2",
             dimensions: dimensions_asr,
             a:1.,
             b:-273.15,
             Units:"C",
             data2: data2_tmp
     	}; // T2M
         Variable mixrat={
             name:"Q2",
             dimensions: dimensions_asr,
             a:1.,
             b:0.,
             Units:"",
             data2: data2_tmp
     	}; // Q2M
         Variable mslp={
             name:"SLP",
             dimensions: dimensions_asr,
             a:1e2,
             b:0.,
             Units:"Pa",
             data2: data2_tmp
     	}; //PSFC, a=1.
         Variable Qsw_in={
             name:"SWDNB",
             dimensions: dimensions_asr,
             a:1.,
             b:0.,
             Units:"W/m^2",
             data2: data2_tmp
     	};
         Variable Qlw_in={
             name:"LWDNB",
             dimensions: dimensions_asr,
             a:1.,
             b:0.,
             Units:"W/m^2",
             data2: data2_tmp
     	};
         Variable snowfr={
             name:"SR",
             dimensions: dimensions_asr,
             a:1.,
             b:0.,
             Units:"",
             data2: data2_tmp
     	};
         Variable precip={
             name:"RAINNC",
             dimensions: dimensions_asr,
             a:nb_timestep_day/(24.*3600),
             b:0.,
             Units:"kg/m^2/s",
             data2: data2_tmp
     	};

         std::vector<Variable> variables_tmp1(7);
         variables_tmp1[0] = tair;
         variables_tmp1[1] = mixrat;
         variables_tmp1[2] = mslp;
         variables_tmp1[3] = Qsw_in;
         variables_tmp1[4] = Qlw_in;
         variables_tmp1[5] = snowfr;
         variables_tmp1[6] = precip;
    
         std::vector<Vectorial_Variable> vectorial_variables_tmp1(0);

         M_asr_elements_dataset={
             dirname:"data",
             prefix:"asr30km.comb.2d.", // "asr30km.comb.2D.";
             postfix:".nc",
             reference_date: "1901-01-01",

             variables: variables_tmp1,
             vectorial_variables: vectorial_variables_tmp1,
             target_size:M_num_elements,
             grid: &M_asr_grid,

             reloaded: false,

             nb_timestep_day: 8,
             time: asr_time,
             dimension_time: asr_dimension_time
     	};

         M_asr_elements_dataset.ftime_range.resize(2,0.);
     }
     else if (strcmp (DatasetName, "topaz_nodes") == 0)
     { 
     	// Definition of topaz grid and datasets
         Dimension topaz_dimension_x={
             name:"x",
             start:0,
             end:760,
             cyclic:false
     	};

         Dimension topaz_dimension_y={
             name:"y",
             start:0,
             end:1100,
             cyclic:false
     	};

         Dimension topaz_dimension_time={
             name:"time", // "Time"
             start:0,
             end:30,
             cyclic:false
     	};

         Dimension topaz_dimension_depth={
             name:"depth", // "Time"
             start:0,
             end:0,
             cyclic:false
     	};

         std::vector<Dimension> dimensions_topaz_uv(4);
         dimensions_topaz_uv[0] = topaz_dimension_time;
         dimensions_topaz_uv[1] = topaz_dimension_depth;
         dimensions_topaz_uv[2] = topaz_dimension_y;
         dimensions_topaz_uv[3] = topaz_dimension_x;

         std::vector<Dimension> dimensions_topaz(3);
         dimensions_topaz[0] = topaz_dimension_time;
         dimensions_topaz[1] = topaz_dimension_y;
         dimensions_topaz[2] = topaz_dimension_x;

         std::vector<Dimension> dimensions_topaz_latlon(2);
         dimensions_topaz_latlon[0] = topaz_dimension_y;
         dimensions_topaz_latlon[1] = topaz_dimension_x;

         std::vector<Dimension> dimensions_topaz_time(1);
         dimensions_topaz_time[0] = topaz_dimension_time;

         Variable topaz_latitude={
             name: "latitude",
             dimensions: dimensions_topaz_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp};

         Variable topaz_longitude={
             name: "longitude",
             dimensions: dimensions_topaz_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp};

         Variable topaz_time={
             name: "time",
             dimensions: dimensions_topaz_time,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp};

         Variable u={
             name: "u",
             dimensions: dimensions_topaz_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v={
             name: "v",
             dimensions: dimensions_topaz_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

     	Variable ssh={
     		name: "ssh",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m/s",
     		data2: data2_tmp
     	};

     	Variable sst={
     		name: "temperature",
     		dimensions: dimensions_topaz_uv,
     		a: 1.,
     		b: 0.,
     		Units: "deg celsius",
     		data2: data2_tmp
     	};

     	Variable sss={
     		name: "salinity",
     		dimensions: dimensions_topaz_uv,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable mld={
     		name: "mlp",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable conc={
     		name: "fice",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable thick={
     		name: "hice",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable snow_thick={
     		name: "hsnow",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

         M_topaz_grid={
             interpolation_method: InterpolationType::FromMeshToMesh2dx,
     		interp_type: -1,
             dirname: "data",
             filename: "TP4DAILY_200803_3m.nc",

             latitude: topaz_latitude,
             longitude: topaz_longitude,

             dimension_x: topaz_dimension_x,
             dimension_y: topaz_dimension_y,

             mpp_file: "NpsNextsim.mpp",
     		interpolation_in_latlon: false,

             loaded: false,

     		masking: true,
     		masking_variable: sss
     	};

         std::vector<Variable> variables_tmp2(3);
         variables_tmp2[0] = u;
         variables_tmp2[1] = v;
         variables_tmp2[2] = ssh;
    
         std::vector<int> uv_tmp0(2);
             uv_tmp0[0] = 0;
             uv_tmp0[1] = 1;
    
         Vectorial_Variable uv={
             components_Id: uv_tmp0,
             east_west_oriented: false // if false, then we assume it is oriented following the mpp_file defined for the grid
     	};
    
         std::vector<Vectorial_Variable> vectorial_variables_tmp2(1);
         vectorial_variables_tmp2[0] = uv;
    
         M_topaz_nodes_dataset={
             dirname: "data",
             prefix: "TP4DAILY_",
             postfix: "_30m.nc",
             reference_date: "1950-01-01",

             variables: variables_tmp2,
             vectorial_variables: vectorial_variables_tmp2,
             target_size: M_num_nodes,
             grid: &M_topaz_grid,

             reloaded: false,

             nb_timestep_day: 1,
             time: topaz_time,
             dimension_time: topaz_dimension_time
     	};

         M_topaz_nodes_dataset.ftime_range.resize(2,0.);
     }
     else if (strcmp (DatasetName, "topaz_elements") == 0)
     { 
     	// Definition of topaz grid and datasets
         Dimension topaz_dimension_x={
             name:"x",
             start:0,
             end:760,
             cyclic:false
     	};

         Dimension topaz_dimension_y={
             name:"y",
             start:0,
             end:1100,
             cyclic:false
     	};

         Dimension topaz_dimension_time={
             name:"time", // "Time"
             start:0,
             end:30,
             cyclic:false
     	};

         Dimension topaz_dimension_depth={
             name:"depth", // "Time"
             start:0,
             end:0,
             cyclic:false
     	};

         std::vector<Dimension> dimensions_topaz_uv(4);
         dimensions_topaz_uv[0] = topaz_dimension_time;
         dimensions_topaz_uv[1] = topaz_dimension_depth;
         dimensions_topaz_uv[2] = topaz_dimension_y;
         dimensions_topaz_uv[3] = topaz_dimension_x;

         std::vector<Dimension> dimensions_topaz(3);
         dimensions_topaz[0] = topaz_dimension_time;
         dimensions_topaz[1] = topaz_dimension_y;
         dimensions_topaz[2] = topaz_dimension_x;

         std::vector<Dimension> dimensions_topaz_latlon(2);
         dimensions_topaz_latlon[0] = topaz_dimension_y;
         dimensions_topaz_latlon[1] = topaz_dimension_x;

         std::vector<Dimension> dimensions_topaz_time(1);
         dimensions_topaz_time[0] = topaz_dimension_time;

         Variable topaz_latitude={
             name: "latitude",
             dimensions: dimensions_topaz_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp};

         Variable topaz_longitude={
             name: "longitude",
             dimensions: dimensions_topaz_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp};

         Variable topaz_time={
             name: "time",
             dimensions: dimensions_topaz_time,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp};

         Variable u={
             name: "u",
             dimensions: dimensions_topaz_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v={
             name: "v",
             dimensions: dimensions_topaz_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

     	Variable ssh={
     		name: "ssh",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m/s",
     		data2: data2_tmp
     	};

     	Variable sst={
     		name: "temperature",
     		dimensions: dimensions_topaz_uv,
     		a: 1.,
     		b: 0.,
     		Units: "deg celsius",
     		data2: data2_tmp
     	};

     	Variable sss={
     		name: "salinity",
     		dimensions: dimensions_topaz_uv,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable mld={
     		name: "mlp",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable conc={
     		name: "fice",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable thick={
     		name: "hice",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable snow_thick={
     		name: "hsnow",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

         M_topaz_grid={
             interpolation_method: InterpolationType::FromMeshToMesh2dx,
     		interp_type: -1,
             dirname: "data",
             filename: "TP4DAILY_200803_3m.nc",

             latitude: topaz_latitude,
             longitude: topaz_longitude,

             dimension_x: topaz_dimension_x,
             dimension_y: topaz_dimension_y,

             mpp_file: "NpsNextsim.mpp",
     		interpolation_in_latlon: false,

             loaded: false,

     		masking: true,
     		masking_variable: sss
     	};

         std::vector<Variable> variables_tmp3(3);
         variables_tmp3[0] = sst;
         variables_tmp3[1] = sss;
         variables_tmp3[2] = mld;
    
         std::vector<Vectorial_Variable> vectorial_variables_tmp3(0);

         M_topaz_elements_dataset={
             dirname: "data",
             prefix: "TP4DAILY_",
             postfix: "_3m.nc",
             reference_date: "1950-01-01",

             variables: variables_tmp3,
             vectorial_variables: vectorial_variables_tmp3,
             target_size: M_num_elements,
             grid: &M_topaz_grid,

             reloaded: false,

             nb_timestep_day: 1,
             time: topaz_time,
             dimension_time: topaz_dimension_time
     	};

         M_topaz_elements_dataset.ftime_range.resize(2,0.);
     }
     else if (strcmp (DatasetName, "ice_topaz_elements") == 0)
     { 
     	// Definition of topaz grid and datasets
         Dimension topaz_dimension_x={
             name:"x",
             start:0,
             end:760,
             cyclic:false
     	};

         Dimension topaz_dimension_y={
             name:"y",
             start:0,
             end:1100,
             cyclic:false
     	};

         Dimension topaz_dimension_time={
             name:"time", // "Time"
             start:0,
             end:30,
             cyclic:false
     	};

         Dimension topaz_dimension_depth={
             name:"depth", // "Time"
             start:0,
             end:0,
             cyclic:false
     	};

         std::vector<Dimension> dimensions_topaz_uv(4);
         dimensions_topaz_uv[0] = topaz_dimension_time;
         dimensions_topaz_uv[1] = topaz_dimension_depth;
         dimensions_topaz_uv[2] = topaz_dimension_y;
         dimensions_topaz_uv[3] = topaz_dimension_x;

         std::vector<Dimension> dimensions_topaz(3);
         dimensions_topaz[0] = topaz_dimension_time;
         dimensions_topaz[1] = topaz_dimension_y;
         dimensions_topaz[2] = topaz_dimension_x;

         std::vector<Dimension> dimensions_topaz_latlon(2);
         dimensions_topaz_latlon[0] = topaz_dimension_y;
         dimensions_topaz_latlon[1] = topaz_dimension_x;

         std::vector<Dimension> dimensions_topaz_time(1);
         dimensions_topaz_time[0] = topaz_dimension_time;

         Variable topaz_latitude={
             name: "latitude",
             dimensions: dimensions_topaz_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp};

         Variable topaz_longitude={
             name: "longitude",
             dimensions: dimensions_topaz_latlon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp};

         Variable topaz_time={
             name: "time",
             dimensions: dimensions_topaz_time,
             a: 1.,
             b: 0.,
             Units: "hours",
             data2: data2_tmp};

         Variable u={
             name: "u",
             dimensions: dimensions_topaz_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v={
             name: "v",
             dimensions: dimensions_topaz_uv,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

     	Variable ssh={
     		name: "ssh",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m/s",
     		data2: data2_tmp
     	};

     	Variable sst={
     		name: "temperature",
     		dimensions: dimensions_topaz_uv,
     		a: 1.,
     		b: 0.,
     		Units: "deg celsius",
     		data2: data2_tmp
     	};

     	Variable sss={
     		name: "salinity",
     		dimensions: dimensions_topaz_uv,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable mld={
     		name: "mlp",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable conc={
     		name: "fice",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "",
     		data2: data2_tmp
     	};

     	Variable thick={
     		name: "hice",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

     	Variable snow_thick={
     		name: "hsnow",
     		dimensions: dimensions_topaz,
     		a: 1.,
     		b: 0.,
     		Units: "m",
     		data2: data2_tmp
     	};

         M_topaz_grid={
             interpolation_method: InterpolationType::FromMeshToMesh2dx,
     		interp_type: -1,
             dirname: "data",
             filename: "TP4DAILY_200803_3m.nc",

             latitude: topaz_latitude,
             longitude: topaz_longitude,

             dimension_x: topaz_dimension_x,
             dimension_y: topaz_dimension_y,

             mpp_file: "NpsNextsim.mpp",
     		interpolation_in_latlon: false,

             loaded: false,

     		masking: true,
     		masking_variable: sss
     	};

         std::vector<Variable> variables_tmp4(3);
         variables_tmp4[0] = conc;
         variables_tmp4[1] = thick;
         variables_tmp4[2] = snow_thick;
    
         std::vector<Vectorial_Variable> vectorial_variables_tmp4(0);
    
         M_ice_topaz_elements_dataset={
             dirname: "data",
             prefix: "TP4DAILY_",
             postfix: "_3m.nc",
             reference_date: "1950-01-01",

             variables: variables_tmp4,
             vectorial_variables: vectorial_variables_tmp4,
             target_size: M_num_elements,
             grid: &M_topaz_grid,

             reloaded: false,

     		nb_timestep_day: 1,
     	    time: topaz_time,
     	    dimension_time: topaz_dimension_time
             };

         M_ice_topaz_elements_dataset.ftime_range.resize(2,0.);
     }
     else if (strcmp (DatasetName, "etopo_elements") == 0)
     { 
     	// Definition of etopo grid and datasets
         Dimension etopo_dimension_x={
             name:"x",
             start:0,
             end:21600,
             cyclic:false
     	};

         Dimension etopo_dimension_y={
             name:"y",
             start:0, approx 60 deg North
             end:10800,
             cyclic:false
     	};

         std::vector<Dimension> dimensions_etopo_lon(1);
         dimensions_etopo_lon[0] = etopo_dimension_x;

         std::vector<Dimension> dimensions_etopo_lat(1);
         dimensions_etopo_lat[0] = etopo_dimension_y;

         std::vector<Dimension> dimensions_etopo(2);
         dimensions_etopo[0] = etopo_dimension_y;
         dimensions_etopo[1] = etopo_dimension_x;

         Variable etopo_latitude={
             name: "y",
             dimensions: dimensions_etopo_lat,
             a: 1.,
             b: 0.,
             Units: "degree_north",
             data2: data2_tmp
     	};

         Variable etopo_longitude={
             name: "x",
             dimensions: dimensions_etopo_lon,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             data2: data2_tmp
     	};

     	M_etopo_grid={
     		interpolation_method: InterpolationType::FromGridToMesh,
     	    //interp_type : TriangleInterpEnum, // slower
     	    interp_type : BilinearInterpEnum,
     	    //interp_type : NearestInterpEnum,
     		dirname:"data",
     		filename:"ETOPO1_Ice_g_gmt4.grd",

     		latitude: etopo_latitude,
     		longitude: etopo_longitude,

     		dimension_x: etopo_dimension_x,
     		dimension_y: etopo_dimension_y,

     		mpp_file: "",
     		interpolation_in_latlon: true,

             loaded: false,

     		masking: false
     	};

         Variable z={
             name:"z",
             dimensions: dimensions_etopo,
             a:-1.,
             b:0.,
             Units:"m",
             data2: data2_tmp
     	};

         std::vector<Variable> variables_tmp5(1);
         variables_tmp5[0] = z;
    
         std::vector<Vectorial_Variable> vectorial_variables_tmp5(0);

         M_etopo_elements_dataset={
             dirname:"data",
             prefix:"ETOPO1_Ice_g_gmt4",
             postfix:".grd",
             reference_date: "",

             variables: variables_tmp5,
             vectorial_variables: vectorial_variables_tmp5,
             target_size:M_num_elements,
             grid: &M_etopo_grid,

             reloaded: false,

             nb_timestep_day: 0
     	};
     }
     else 
#endif         
         if (strcmp (DatasetName, "ERAi_nodes") == 0)
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
     		dirname:"data",
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
         Variable u10={
             name: "10U", // U10M
             dimensions: dimensions,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         Variable v10={
             name: "10V", // U10M
             dimensions: dimensions,
             a: 1.,
             b: 0.,
             Units: "m/s",
             data2: data2_tmp
     	};

         std::vector<Variable> variables(2);
         variables[0] = u10;
         variables[1] = v10;

         std::vector<int> uv_tmp0(2);
             uv_tmp0[0] = 0;
             uv_tmp0[1] = 1;
    
         Vectorial_Variable uv10={
             components_Id: uv_tmp0,
             east_west_oriented: true // if false, then we assume it is oriented following the input grid
     	};

         std::vector<Vectorial_Variable> vectorial_variables(1);
         vectorial_variables[0] = uv10;

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

         this->ftime_range.resize(2,0.);
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
       //close_Dataset (this);
     }

   }

} // Nextsim
