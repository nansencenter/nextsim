/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   dataset.cpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#include <dataset.hpp>
#include <date.hpp>
#include <constants.hpp>
#include <numeric>

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

DataSet::DataSet(char const *DatasetName)
    :
    M_log_level(Environment::logLevel()),
    M_log_all(Environment::logAll()),
    M_comm(Environment::comm())
{
    name = std::string(DatasetName);
    projfilename = Environment::vm()["mesh.mppfile"].as<std::string>();

    ftime_range.resize(2,0.);
#ifdef OASIS
    itime_range.resize(2,0.);
    calc_nodal_weights = false;
#endif


    std::vector<std::vector<double>> loaded_data_tmp;
    loaded_data_tmp.resize(2);

    std::vector<std::vector<double>> interpolated_data_tmp;
    interpolated_data_tmp.resize(2);

#ifdef OASIS
    // Extract the dirname, prefix, and postfix from coupler.exchange_grid_file (default is coupler/NEMO.nc).
    boost::filesystem::path const exchange_grid_file( Environment::vm()["coupler.exchange_grid_file"].as<std::string>() );
#endif

    /*
     *  match projection name and initialize remaining parameters
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
            name:"Time",
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
            filename_string: "", // All variables are in the same (grid) file
            name: "XLAT",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "XLONG",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "Time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // conversion factors: xnew = a*x + b
        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "U10M", // U10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "V10M", // V10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum,  // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,

            dirname:"",
            filename_mask: "asr30km.comb.2D.%Y%m.nc",
            gridfile: "",
            reference_date: "1901-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "NpsASR.mpp",
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"monthly",

            waveOptions: wavopt_none,

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

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
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
            name:"Time",
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
            filename_string: "", // All variables are in the same (grid) file
            name: "XLAT",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "XLONG",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "Time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum,  // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,

            dirname:"",
            filename_mask: "asr30km.comb.2D.%Y%m.nc",
            gridfile: "",
            reference_date: "1901-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "NpsASR.mpp",
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"monthly",

            waveOptions: wavopt_none,

            masking: false
        };


        // conversion factors: xnew = a*x + b
        Variable tair={
            filename_string: "", // All variables are in the same (grid) file
            name:"T2M",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // T2M
        Variable mixrat={
            filename_string: "", // All variables are in the same (grid) file
            name:"Q2M",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // Q2M
        Variable mslp={
            filename_string: "", // All variables are in the same (grid) file
            name:"PSFC",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"Pa",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; //PSFC, a=1.
        Variable Qsw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"SWDNB",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable Qlw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"LWDNB",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable snowfr={
            filename_string: "", // All variables are in the same (grid) file
            name:"SR",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable precip={
            filename_string: "", // All variables are in the same (grid) file
            name:"RAINNC",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:physical::rhow/1000./(3.*3600),
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
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

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "generic_ps_atm_nodes") == 0)
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
            name:"time",
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,//convert days to hours
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // conversion factors: xnew = a*x + b
        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "x_wind_10m", // U10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "y_wind_10m", // V10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum,  // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,

            dirname:"GENERIC_PS_ATM",
            filename_mask: "generic_ps_atm_%Y%m%d.nc",
            gridfile: "",
            reference_date: "1900-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "NpsNextsim.mpp",
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: false
        };

        std::vector<Variable> variables_tmp = { u, v };
        std::vector<int> uv_tmp(2);
        uv_tmp[0] = 0;
        uv_tmp[1] = 1;

        Vectorial_Variable uv={
            components_Id: uv_tmp,
            east_west_oriented: false // if false, then we assume it is oriented following the input grid
        };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
        vectorial_variables_tmp[0] = uv;

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "generic_ps_atm_elements") == 0)
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
            name:"time",
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,//convert days to hours
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum,  // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,

            dirname:"GENERIC_PS_ATM",
            filename_mask: "generic_ps_atm_%Y%m%d.nc",
            gridfile: "",
            reference_date: "1900-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "NpsNextsim.mpp",
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: false
        };


        // conversion factors: xnew = a*x + b
        Variable tair={
            filename_string: "", // All variables are in the same (grid) file
            name:"air_temperature_2m",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // T2M
        Variable dair={
            filename_string: "", // All variables are in the same (grid) file
            name:"dew_point_temperature_2m",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // Q2M
        Variable mslp={
            filename_string: "", // All variables are in the same (grid) file
            name:"atm_pressure",       
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"Pa",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; //PSFC, a=1.
        // conversion factors: xnew = a*x + b
        Variable Qsw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"instantaneous_downwelling_shortwave_radiation",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,//integrated for 1 hour - convert from total to rate
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable Qlw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"instantaneous_downwelling_longwave_radiation",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,//integrated for 1 hour - convert from total to rate
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable snowfall={
            filename_string: "", // All variables are in the same (grid) file
            name:"snowfall_rate",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,//integrated for 1 hour - convert from total to rate
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable precip={
            filename_string: "", // All variables are in the same (grid) file
            name:"total_precipitation_rate",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,//integrated for 1 hour - convert from total to rate
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        std::vector<Variable> variables_tmp
            = { tair, dair, mslp, Qsw_in, Qlw_in, snowfall, precip };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ec2_arome_nodes") == 0)
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
            name:"time",
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1./3600.,//convert seconds to hours
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // conversion factors: xnew = a*x + b
        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "x_wind_10m", // U10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "y_wind_10m", // V10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum,  // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,

            dirname:"EC2_AROME",
            filename_mask: "ec2_arome_blended_%Y%m%d.nc",
            gridfile: "",
            reference_date: "1970-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "NpsNextsim.mpp",
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: false
        };

        std::vector<Variable> variables_tmp = { u, v };
        std::vector<int> uv_tmp(2);
        uv_tmp[0] = 0;
        uv_tmp[1] = 1;

        Vectorial_Variable uv={
            components_Id: uv_tmp,
            east_west_oriented: false // if false, then we assume it is oriented following the input grid
        };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
        vectorial_variables_tmp[0] = uv;

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ec2_arome_elements") == 0)
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
            name:"time",
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1./3600.,//convert seconds to hours
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum,  // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,

            dirname:"EC2_AROME",
            filename_mask: "ec2_arome_blended_%Y%m%d.nc",
            gridfile: "",
            reference_date: "1970-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "NpsNextsim.mpp",
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: false
        };


        // conversion factors: xnew = a*x + b
        Variable tair={
            filename_string: "", // All variables are in the same (grid) file
            name:"air_temperature_2m",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // T2M
        Variable sphuma={
            filename_string: "", // All variables are in the same (grid) file
            name:"specific_humidity_2m",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // Q2M
        Variable mslp={
            filename_string: "", // All variables are in the same (grid) file
            name:"air_pressure_at_sea_level",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"Pa",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; //PSFC, a=1.
        // conversion factors: xnew = a*x + b
        Variable Qsw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./3600.,//integrated for 1 hour - convert from total to rate
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable Qlw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"integral_of_surface_downwelling_longwave_flux_in_air_wrt_time",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./3600.,//integrated for 1 hour - convert from total to rate
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable snowfall={
            filename_string: "", // All variables are in the same (grid) file
            name:"integral_of_snowfall_amount_wrt_time",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./3600.,//integrated for 1 hour - convert from total to rate
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable precip={
            filename_string: "", // All variables are in the same (grid) file
            name:"precipitation_amount_acc",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./3600.,//integrated for 1 hour - convert from total to rate
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        std::vector<Variable> variables_tmp
            = { tair, sphuma, mslp, Qsw_in, Qlw_in, snowfall, precip };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ec2_arome_ensemble_nodes") == 0)
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
            name:"time",
            cyclic:false
        };
        Dimension dimension_ens={
            name:"ensemble_member",
            cyclic:false
        };

        std::vector<Dimension> dimensions_latlon = {
            dimension_y, dimension_x};
        std::vector<Dimension> dimensions = {
            dimension_time, dimension_ens, dimension_y, dimension_x};
        std::vector<Dimension> dimensions_time(1);
        dimensions_time[0] = dimension_time;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1./3600.,//convert seconds to hours
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // conversion factors: xnew = a*x + b
        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "x_wind_10m", // U10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "y_wind_10m", // V10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum,  // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,

            dirname:"EC2_AROME_ENSEMBLE",
            filename_mask: "ec2_arome_blended_ensemble_%Y%m%d.nc",
            gridfile: "",
            reference_date: "1970-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "NpsNextsim.mpp",
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: false
        };

        std::vector<Variable> variables_tmp = { u, v };
        std::vector<int> uv_tmp(2);
        uv_tmp[0] = 0;
        uv_tmp[1] = 1;

        Vectorial_Variable uv={
            components_Id: uv_tmp,
            east_west_oriented: false // if false, then we assume it is oriented following the input grid
        };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
        vectorial_variables_tmp[0] = uv;

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ec2_arome_ensemble_elements") == 0)
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
            name:"time",
            cyclic:false
        };
        Dimension dimension_ens={
            name:"ensemble_member",
            cyclic:false
        };

        std::vector<Dimension> dimensions_latlon = {
            dimension_y, dimension_x};
        std::vector<Dimension> dimensions = {
            dimension_time, dimension_ens, dimension_y, dimension_x};
        std::vector<Dimension> dimensions_time(1);
        dimensions_time[0] = dimension_time;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1./3600.,//convert seconds to hours
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum,  // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,

            dirname:"EC2_AROME_ENSEMBLE",
            filename_mask: "ec2_arome_blended_ensemble_%Y%m%d.nc",
            gridfile: "",
            reference_date: "1970-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "NpsNextsim.mpp",
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: false
        };


        // conversion factors: xnew = a*x + b
        Variable tair={
            filename_string: "", // All variables are in the same (grid) file
            name:"air_temperature_2m",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // T2M
        Variable sphuma={
            filename_string: "", // All variables are in the same (grid) file
            name:"specific_humidity_2m",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // Q2M
        Variable mslp={
            filename_string: "", // All variables are in the same (grid) file
            name:"air_pressure_at_sea_level",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"Pa",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; //PSFC, a=1.
        // conversion factors: xnew = a*x + b
        Variable Qsw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./(3*3600.), //integrated for 3 hours - convert from total to rate
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable Qlw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"integral_of_surface_downwelling_longwave_flux_in_air_wrt_time",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./(3*3600.), //integrated for 3 hours - convert from total to rate
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable snowfall={
            filename_string: "", // All variables are in the same (grid) file
            name:"integral_of_snowfall_amount_wrt_time",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./(3*3600.), //integrated for 3 hours - convert from total to rate
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable precip={
            filename_string: "", // All variables are in the same (grid) file
            name:"precipitation_amount_acc",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./(3*3600.), //integrated for 3 hours - convert from total to rate
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        std::vector<Variable> variables_tmp = {
            tair, sphuma, mslp, Qsw_in, Qlw_in, snowfall, precip };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
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
            name:"time",
            cyclic:false
        };

        Dimension dimension_depth={
            name:"depth",
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "u",
            dimensions: dimensions_uv,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "v",
            dimensions: dimensions_uv,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable ssh={
            filename_string: "", // All variables are in the same (grid) file
            name: "ssh",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "TP4DAILY_%Y%m_30m.nc",
            gridfile: "",
            reference_date: "1950-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"monthly",

            waveOptions: wavopt_none,

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
            east_west_oriented: false
                // if false, then we assume it is oriented
                // following the mpp_file defined for the grid
        };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
        vectorial_variables_tmp[0] = uv;

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable sst={
            filename_string: "", // All variables are in the same (grid) file
            name: "temperature",
            dimensions: dimensions_uv,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "deg celsius",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable sss={
            filename_string: "", // All variables are in the same (grid) file
            name: "salinity",
            dimensions: dimensions_uv,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable mld={
            filename_string: "", // All variables are in the same (grid) file
            name: "mlp",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "TP4DAILY_%Y%m_3m.nc",
            gridfile: "",
            reference_date: "1950-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"monthly",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: sss
        };

        std::vector<Variable> variables_tmp(3);
        variables_tmp[0] = sst;
        variables_tmp[1] = sss;
        variables_tmp[2] = mld;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "glorys12_nodes") == 0)
    {
        // Definition of glorys12 grid and datasets
        Dimension dimension_x={
        name:"longitude",
        cyclic:false//true
        };

        Dimension dimension_y={
        name:"latitude",
        cyclic:false
        };

        Dimension dimension_time={
        name:"time",
        cyclic:false
        };

        Dimension dimension_depth={
        name:"depth",
        cyclic:false
        };

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        std::vector<Dimension> dimensions_time(1);
        dimensions_time[0] = dimension_time;

        std::vector<Dimension> dimensions_uv(4);
        dimensions_uv[0] = dimension_time;
        dimensions_uv[1] = dimension_depth;
        dimensions_uv[2] = dimension_y;
        dimensions_uv[3] = dimension_x;

        std::vector<Dimension> dimensions(3);
        dimensions[0] = dimension_time;
        dimensions[1] = dimension_y;
        dimensions[2] = dimension_x;

        //VARIABLES:
        Variable latitude={
        filename_string: "", // All variables are in the same (grid) file
        name: "latitude",
        dimensions: dimensions_lat,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units: "degree_north",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Variable longitude={
        filename_string: "", // All variables are in the same (grid) file
        name: "longitude",
        dimensions: dimensions_lon,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units: "degree_east",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Variable time_tmp={
        filename_string: "", // All variables are in the same (grid) file
        name: "time",
        dimensions: dimensions_time,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units: "hours",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Variable u={
        filename_string: "", // All variables are in the same (grid) file
        name: "uo",
        dimensions: dimensions_uv,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units: "m/s",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Variable v={
        filename_string: "", // All variables are in the same (grid) file
        name: "vo",
        dimensions: dimensions_uv,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units: "m/s",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Variable ssh={
        filename_string: "", // All variables are in the same (grid) file
        name: "zos",
        dimensions: dimensions,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units: "m",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
        interpolation_method: InterpolationType::FromGridToMesh,
        interp_type : BilinearInterpEnum,
        //interpolation_method: InterpolationType::FromMeshToMesh2dx,
        //interp_type: -1,
        dirname: "",
        filename_mask: "GLORYS12V1_%Y_30m.nc",
        gridfile: "",
        reference_date: "1950-01-01",

        latitude: latitude,
        longitude: longitude,

        dimension_x: dimension_x,
        dimension_y: dimension_y,

        mpp_file: "",
        interpolation_in_latlon: true,

        loaded: false,
        dataset_frequency:"yearly",

        waveOptions: wavopt_none,

        masking: false,
        //masking_variable:
        };

        //VAR
        std::vector<Variable> variables_tmp(3);
        variables_tmp[0] = u;
        variables_tmp[1] = v;
        variables_tmp[2] = ssh;

        //VECTOR
        std::vector<int> uv_tmp(2);
        uv_tmp[0] = 0;
        uv_tmp[1] = 1;

        Vectorial_Variable uv={
        components_Id: uv_tmp,
        east_west_oriented: true
        };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
        vectorial_variables_tmp[0] = uv;

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days lolo?
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "glorys12_elements") == 0)
    {
        // Definition of glorys12 grid and datasets
        Dimension dimension_x={
        name:"longitude",
        cyclic:false//true
        };

        Dimension dimension_y={
        name:"latitude",
        cyclic:false
        };

        Dimension dimension_time={
        name:"time",
        cyclic:false
        };

        Dimension dimension_depth={
        name:"depth",
        cyclic:false
        };

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        std::vector<Dimension> dimensions_time(1);
        dimensions_time[0] = dimension_time;

        std::vector<Dimension> dimensions_uv(4);
        dimensions_uv[0] = dimension_time;
        dimensions_uv[1] = dimension_depth;
        dimensions_uv[2] = dimension_y;
        dimensions_uv[3] = dimension_x;

        std::vector<Dimension> dimensions(3);
        dimensions[0] = dimension_time;
        dimensions[1] = dimension_y;
        dimensions[2] = dimension_x;

        //VARIABLES:
        Variable latitude={
        filename_string: "", // All variables are in the same (grid) file
        name: "latitude",
        dimensions: dimensions_lat,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units: "degree_north",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Variable longitude={
        filename_string: "", // All variables are in the same (grid) file
        name: "longitude",
        dimensions: dimensions_lon,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units: "degree_east",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Variable time_tmp={
        filename_string: "", // All variables are in the same (grid) file
        name: "time",
        dimensions: dimensions_time,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units: "hours",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Variable sst={
        filename_string: "", // All variables are in the same (grid) file
        name: "thetao",
        dimensions: dimensions_uv,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units:"C",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Variable sss={
        filename_string: "", // All variables are in the same (grid) file
        name: "so",
        dimensions: dimensions_uv,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units: "",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Variable mld={
        filename_string: "", // All variables are in the same (grid) file
        name: "mlotst",
        dimensions: dimensions,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        use_FillValue: true,
        use_missing_value: true,
        a: 1.,
        b: 0.,
        Units: "m",
        loaded_data: loaded_data_tmp,
        interpolated_data: interpolated_data_tmp,
        wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
        interpolation_method: InterpolationType::FromGridToMesh,
        interp_type : BilinearInterpEnum,
        //interpolation_method: InterpolationType::FromMeshToMesh2dx,
        //interp_type: -1,
        dirname: "",
        filename_mask: "GLORYS12V1_%Y_3m.nc",
        gridfile: "",
        reference_date: "1950-01-01",

        latitude: latitude,
        longitude: longitude,

        dimension_x: dimension_x,
        dimension_y: dimension_y,

        mpp_file: "",
        interpolation_in_latlon: true,

        loaded: false,
        dataset_frequency:"yearly",

        waveOptions: wavopt_none,

        masking: false,
        //masking_variable: sss
        };

        //VAR
        std::vector<Variable> variables_tmp(3);
        variables_tmp[0] = sst;
        variables_tmp[1] = sss;
        variables_tmp[2] = mld;

        //VECTOR
        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days lolo?
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
#ifdef OASIS
    else if (strcmp (DatasetName, "ocean_cpl_nodes") == 0)
    {
         // Definition of ocean coupling grid and datasets
        Dimension dimension_x={
            name:"y",
            cyclic:false
        };

        Dimension dimension_y={
            name:"x",
            cyclic:false
        };

        std::vector<Dimension> dimensions(2);
        dimensions[0] = dimension_y;
        dimensions[1] = dimension_x;

        // Variables received through OASIS
        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "I_Uocn",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "I_Vocn",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable ssh={
            filename_string: "", // All variables are in the same (grid) file
            name: "I_SSH",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // These are lat, lon, masking, and vector rotation variables in NEMO.nc
        Variable mask={
            filename_string: "", // All variables are in the same (grid) file
            name: "mask",
            dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable theta={
            filename_string: "", // All variables are in the same (grid) file
            name: "ptheta",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "radians",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "plat",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "plon",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMeshQuick,
            interp_type: -1,
            dirname: exchange_grid_file.parent_path().string(),
            filename_mask: exchange_grid_file.filename().string(),
            gridfile: "",
            reference_date: "1979-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"coupled",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: mask,

            gridded_rotation_angle: true,
            vector_rotation_variable: theta
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
            east_west_oriented: false
                // if false, then we assume it is oriented
                // following the mpp_file defined for the grid
        };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
        vectorial_variables_tmp[0] = uv;

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
        coupled = true;
    }
    else if (strcmp (DatasetName, "ocean_cpl_elements") == 0)
    {
         // Definition of ocean coupling grid and datasets
        Dimension dimension_x={
            name:"y",
            cyclic:false
        };

        Dimension dimension_y={
            name:"x",
            cyclic:false
        };

        std::vector<Dimension> dimensions(2);
        dimensions[0] = dimension_y;
        dimensions[1] = dimension_x;

        // Variables received through OASIS
        Variable sst={
            filename_string: "", // All variables are in the same (grid) file
            name: "I_SST",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "deg celsius",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable sss={
            filename_string: "", // All variables are in the same (grid) file
            name: "I_SSS",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable qsrml={
            filename_string: "", // All variables are in the same (grid) file
            name: "I_FrcQsr",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable mld={
            filename_string: "", // All variables are in the same (grid) file
            name: "I_MLD",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        // The masking, lon, and lat variables in NEMO.nc
        Variable mask={
            filename_string: "", // All variables are in the same (grid) file
            name: "mask",
            dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "plat",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "plon",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Grid grid_tmp={
            interpolation_method: InterpolationType::ConservativeRemapping,
            interp_type: -1,
            dirname: exchange_grid_file.parent_path().string(),
            filename_mask: exchange_grid_file.filename().string(),
            gridfile: "",
            reference_date: "1979-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"coupled",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: mask
        };

        std::vector<Variable> variables_tmp(3);
        variables_tmp[0] = sst;
        variables_tmp[1] = sss;
        variables_tmp[2] = qsrml;
        if ( Environment::vm()["coupler.rcv_first_layer_depth"].as<bool>() )
            variables_tmp.push_back(mld);

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
        coupled = true;
    }
    else if (strcmp (DatasetName, "wave_cpl_nodes") == 0)
    {
        Dimension dimension_x={
            name:"y",
            cyclic:false
        };

        Dimension dimension_y={
            name:"x",
            cyclic:false
        };

        std::vector<Dimension> dimensions = {dimension_y, dimension_x};

        // Variables received through OASIS
        Variable tau_wi_x={
            filename_string: "", // All variables are in the same (grid) file
            name: "I_tauwix",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: physical::rhow,//1., //WW3 WRS has unit m^2/s^2 TODO could be different in other models (option)
            b: 0.,
            Units: "Pa",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable tau_wi_y={
            filename_string: "", // All variables are in the same (grid) file
            name: "I_tauwiy",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: physical::rhow,//1., //WW3 WRS has unit m^2/s^2 TODO could be different in other models (option)
            b: 0.,
            Units: "Pa",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // These are lat, lon, masking, and vector rotation variables in NEMO.nc
        Variable mask={
            filename_string: "", // All variables are in the same (grid) file
            name: "mask",
            dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable theta={
            filename_string: "", // All variables are in the same (grid) file
            name: "ptheta",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "radians",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "plat",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "plon",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp = {
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Grid grid_tmp = {
            interpolation_method: InterpolationType::FromMeshToMeshQuick,
            interp_type: -1,
            dirname: exchange_grid_file.parent_path().string(),
            filename_mask: exchange_grid_file.filename().string(),
            gridfile: "",
            reference_date: "1979-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"coupled",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: mask,

            gridded_rotation_angle: true,
            vector_rotation_variable: theta
        };

        variables = {tau_wi_x, tau_wi_y};
        Vectorial_Variable tau_xy = {
            components_Id: {0, 1},
            east_west_oriented: true
                // if false, then we assume it is oriented
                // following the mpp_file defined for the grid
            };
        vectorial_variables = {tau_xy};

        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.; // days
        time= time_tmp;
        coupled = true;
    }
    else if (strcmp (DatasetName, "wave_cpl_elements") == 0)
    {
        Dimension dimension_x={
            name:"y",
            cyclic:false
        };

        Dimension dimension_y={
            name:"x",
            cyclic:false
        };

        std::vector<Dimension> dimensions = {dimension_y, dimension_x};

        // Variables received through OASIS
        //Variable str_var ={
        //    filename_string: "", // All variables are in the same (grid) file
        //    name: "I_str_var",
        //    dimensions: dimensions,
        //    land_mask_defined: false,
        //    land_mask_value: 0.,
        //    NaN_mask_defined: false,
        //    NaN_mask_value: 0.,
        //    use_FillValue: true,
        //    use_missing_value: true,
        //    a: 1.,
        //    b: 0.,
        //    Units: "",
        //    loaded_data: loaded_data_tmp,
        //    interpolated_data: interpolated_data_tmp,
        //    wavDirOptions: wavdiropt_none
        //};

        Variable wlbk ={
            filename_string: "", // All variables are in the same (grid) file
            name: "I_wlbk",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 2.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        //Variable Tm02 ={
        //    filename_preissue_241_oasis_ww3_rcv_dmaxfix: "", // All variables are in the same (grid) file
        //    name: "I_tm02",
        //    dimensions: dimensions,
        //    land_mask_defined: false,
        //    land_mask_value: 0.,
        //    NaN_mask_defined: false,
        //    NaN_mask_value: 0.,
        //    use_FillValue: true,
        //    use_missing_value: true,
        //    a: 1.,
        //    b: 0.,
        //    Units: "s",
        //    loaded_data: loaded_data_tmp,
        //    interpolated_data: interpolated_data_tmp,
        //    wavDirOptions: wavdiropt_none
        //};

        // These are lat, lon, masking, and vector rotation variables in NEMO.nc
        Variable mask={
            filename_string: "", // All variables are in the same (grid) file
            name: "mask",
            dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "plat",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "plon",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp = {
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Grid grid_tmp = {
            interpolation_method: InterpolationType::ConservativeRemapping,
            interp_type: -1,
            dirname: exchange_grid_file.parent_path().string(),
            filename_mask: exchange_grid_file.filename().string(),
            gridfile: "",
            reference_date: "1979-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"coupled",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: mask
        };

        variables = {wlbk};
        vectorial_variables = {};

        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.; // days
        time = time_tmp;
        coupled = true;
    }
#endif
    else if (strcmp (DatasetName, "nesting_ocean_elements") == 0)
    {
        // Definition of besting grid and datasets
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable sst={
            filename_string: "", // All variables are in the same (grid) file
            name: "sea_surface_temperature",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "deg celsius",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable sss={
            filename_string: "", // All variables are in the same (grid) file
            name: "sea_surface_salinity",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "deg celsius",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "nesting_"
               +Environment::vm()["nesting.outer_mesh"].as<std::string>()
               +"_%Y%m%d.nc",
            gridfile: "nesting_grid_"
                +Environment::vm()["nesting.inner_mesh"].as<std::string>()
                +".nc",
            reference_date: "1900-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: sst
        };

        std::vector<Variable> variables_tmp(2);
        variables_tmp[0] = sst;
        variables_tmp[1] = sss;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0; // hours
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "nesting_ice_elements") == 0)
    {
        if ( Environment::vm()["thermo.newice_type"].as<int>() == 4 ) {
            // Definition of nesting grid and datasets
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
            filename_string: "", // All variables are in the same (grid) file
                name: "latitude",
                dimensions: dimensions_latlon,
                land_mask_defined: false,
                land_mask_value: 0.,
                NaN_mask_defined: false,
                NaN_mask_value: 0.,
                use_FillValue: true,
                use_missing_value: true,
                a: 1.,
                b: 0.,
                Units: "degree_north",
                loaded_data: loaded_data_tmp,
                interpolated_data: interpolated_data_tmp,
                wavDirOptions: wavdiropt_none};

            Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
                name: "longitude",
                dimensions: dimensions_latlon,
                land_mask_defined: false,
                land_mask_value: 0.,
                NaN_mask_defined: false,
                NaN_mask_value: 0.,
                use_FillValue: true,
                use_missing_value: true,
                a: 1.,
                b: 0.,
                Units: "degree_east",
                loaded_data: loaded_data_tmp,
                interpolated_data: interpolated_data_tmp,
                wavDirOptions: wavdiropt_none};

            Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
                name: "time",
                dimensions: dimensions_time,
                land_mask_defined: false,
                land_mask_value: 0.,
                NaN_mask_defined: false,
                NaN_mask_value: 0.,
                use_FillValue: true,
                use_missing_value: true,
                a: 24.,
                b: 0.,
                Units: "hours",
                loaded_data: loaded_data_tmp,
                interpolated_data: interpolated_data_tmp,
                wavDirOptions: wavdiropt_none};

            Variable sit={
            filename_string: "", // All variables are in the same (grid) file
                name: "sea_ice_thickness",
                dimensions: dimensions,
                land_mask_defined: false,
                land_mask_value: 0.,
                NaN_mask_defined: false,
                NaN_mask_value: 0.,
                use_FillValue: true,
                use_missing_value: true,
                a: 1.,
                b: 0.,
                Units: "m",
                loaded_data: loaded_data_tmp,
                interpolated_data: interpolated_data_tmp,
                wavDirOptions: wavdiropt_none
            };

            Variable sic={
            filename_string: "", // All variables are in the same (grid) file
                name: "sea_ice_area_fraction",
                dimensions: dimensions,
                land_mask_defined: false,
                land_mask_value: 0.,
                NaN_mask_defined: false,
                NaN_mask_value: 0.,
                use_FillValue: true,
                use_missing_value: true,
                a: 1.,
                b: 0.,
                Units: "1",
                loaded_data: loaded_data_tmp,
                interpolated_data: interpolated_data_tmp,
                wavDirOptions: wavdiropt_none
            };

            Variable snt={
            filename_string: "", // All variables are in the same (grid) file
                name: "surface_snow_thickness",
                dimensions: dimensions,
                land_mask_defined: false,
                land_mask_value: 0.,
                NaN_mask_defined: false,
                NaN_mask_value: 0.,
                use_FillValue: true,
                use_missing_value: true,
                a: 1.,
                b: 0.,
                Units: "m",
                loaded_data: loaded_data_tmp,
                interpolated_data: interpolated_data_tmp,
                wavDirOptions: wavdiropt_none
            };

            Variable sit_young={
            filename_string: "", // All variables are in the same (grid) file
                name: "young_ice_thickness",
                dimensions: dimensions,
                land_mask_defined: false,
                land_mask_value: 0.,
                NaN_mask_defined: false,
                NaN_mask_value: 0.,
                use_FillValue: true,
                use_missing_value: true,
                a: 1.,
                b: 0.,
                Units: "m",
                loaded_data: loaded_data_tmp,
                interpolated_data: interpolated_data_tmp,
                wavDirOptions: wavdiropt_none
            };

            Variable sic_young={
            filename_string: "", // All variables are in the same (grid) file
                name: "young_ice_area_fraction",
                dimensions: dimensions,
                land_mask_defined: false,
                land_mask_value: 0.,
                NaN_mask_defined: false,
                NaN_mask_value: 0.,
                use_FillValue: true,
                use_missing_value: true,
                a: 1.,
                b: 0.,
                Units: "1",
                loaded_data: loaded_data_tmp,
                interpolated_data: interpolated_data_tmp,
                wavDirOptions: wavdiropt_none
            };

            Variable snt_young={
            filename_string: "", // All variables are in the same (grid) file
                name: "surface_snow_thickness_on_young_ice",
                dimensions: dimensions,
                land_mask_defined: false,
                land_mask_value: 0.,
                NaN_mask_defined: false,
                NaN_mask_value: 0.,
                use_FillValue: true,
                use_missing_value: true,
                a: 1.,
                b: 0.,
                Units: "m",
                loaded_data: loaded_data_tmp,
                interpolated_data: interpolated_data_tmp,
                wavDirOptions: wavdiropt_none
            };

            Grid grid_tmp={
                interpolation_method: InterpolationType::FromMeshToMesh2dx,
                interp_type: -1,
                dirname: "",
                filename_mask: "nesting_"
                   +Environment::vm()["nesting.outer_mesh"].as<std::string>()
                   +"_%Y%m%d.nc",
                gridfile: "nesting_grid_"
                    +Environment::vm()["nesting.inner_mesh"].as<std::string>()
                    +".nc",
                reference_date: "1900-01-01",

                latitude: latitude,
                longitude: longitude,

                dimension_x: dimension_x,
                dimension_y: dimension_y,

                mpp_file: projfilename,
                interpolation_in_latlon: false,

                loaded: false,
                dataset_frequency:"daily",

                waveOptions: wavopt_none,

                masking: true,
                masking_variable: sit
            };

            std::vector<Variable> variables_tmp(6);
            variables_tmp[0] = sit;
            variables_tmp[1] = sic;
            variables_tmp[2] = snt;
            variables_tmp[3] = sit_young;
            variables_tmp[4] = sic_young;
            variables_tmp[5] = snt_young;

            std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

            variables= variables_tmp;
            vectorial_variables= vectorial_variables_tmp;
            grid= grid_tmp;

            loaded=false;
            interpolated=false;

            averaging_period=0; // hours
            time= time_tmp;
    #ifdef OASIS
            coupled = false;
    #endif
        }
        else{
          // Definition of nesting grid and datasets
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
            filename_string: "", // All variables are in the same (grid) file
              name: "latitude",
              dimensions: dimensions_latlon,
              land_mask_defined: false,
              land_mask_value: 0.,
              NaN_mask_defined: false,
              NaN_mask_value: 0.,
              use_FillValue: true,
              use_missing_value: true,
              a: 1.,
              b: 0.,
              Units: "degree_north",
              loaded_data: loaded_data_tmp,
              interpolated_data: interpolated_data_tmp,
              wavDirOptions: wavdiropt_none};

          Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
              name: "longitude",
              dimensions: dimensions_latlon,
              land_mask_defined: false,
              land_mask_value: 0.,
              NaN_mask_defined: false,
              NaN_mask_value: 0.,
              use_FillValue: true,
              use_missing_value: true,
              a: 1.,
              b: 0.,
              Units: "degree_east",
              loaded_data: loaded_data_tmp,
              interpolated_data: interpolated_data_tmp,
              wavDirOptions: wavdiropt_none};

          Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
              name: "time",
              dimensions: dimensions_time,
              land_mask_defined: false,
              land_mask_value: 0.,
              NaN_mask_defined: false,
              NaN_mask_value: 0.,
              use_FillValue: true,
              use_missing_value: true,
              a: 24.,
              b: 0.,
              Units: "hours",
              loaded_data: loaded_data_tmp,
              interpolated_data: interpolated_data_tmp,
              wavDirOptions: wavdiropt_none};

          Variable sit={
            filename_string: "", // All variables are in the same (grid) file
              name: "sea_ice_thickness",
              dimensions: dimensions,
              land_mask_defined: false,
              land_mask_value: 0.,
              NaN_mask_defined: false,
              NaN_mask_value: 0.,
              use_FillValue: true,
              use_missing_value: true,
              a: 1.,
              b: 0.,
              Units: "m",
              loaded_data: loaded_data_tmp,
              interpolated_data: interpolated_data_tmp,
              wavDirOptions: wavdiropt_none
          };

          Variable sic={
            filename_string: "", // All variables are in the same (grid) file
              name: "sea_ice_area_fraction",
              dimensions: dimensions,
              land_mask_defined: false,
              land_mask_value: 0.,
              NaN_mask_defined: false,
              NaN_mask_value: 0.,
              use_FillValue: true,
              use_missing_value: true,
              a: 1.,
              b: 0.,
              Units: "1",
              loaded_data: loaded_data_tmp,
              interpolated_data: interpolated_data_tmp,
              wavDirOptions: wavdiropt_none
          };

          Variable snt={
            filename_string: "", // All variables are in the same (grid) file
              name: "surface_snow_thickness",
              dimensions: dimensions,
              land_mask_defined: false,
              land_mask_value: 0.,
              NaN_mask_defined: false,
              NaN_mask_value: 0.,
              use_FillValue: true,
              use_missing_value: true,
              a: 1.,
              b: 0.,
              Units: "m",
              loaded_data: loaded_data_tmp,
              interpolated_data: interpolated_data_tmp,
              wavDirOptions: wavdiropt_none
          };

          Grid grid_tmp={
              interpolation_method: InterpolationType::FromMeshToMesh2dx,
              interp_type: -1,
              dirname: "",
              filename_mask: "nesting_"
                   +Environment::vm()["nesting.outer_mesh"].as<std::string>()
                   +"_%Y%m%d.nc",
              gridfile: "nesting_grid_"
                    +Environment::vm()["nesting.inner_mesh"].as<std::string>()
                    +".nc",
              reference_date: "1900-01-01",

              latitude: latitude,
              longitude: longitude,

              dimension_x: dimension_x,
              dimension_y: dimension_y,

              mpp_file: projfilename,
              interpolation_in_latlon: false,

              loaded: false,
              dataset_frequency:"daily",

              waveOptions: wavopt_none,

              masking: true,
              masking_variable: sit
          };

          std::vector<Variable> variables_tmp(3);
          variables_tmp[0] = sit;
          variables_tmp[1] = sic;
          variables_tmp[2] = snt;

          std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

          variables= variables_tmp;
          vectorial_variables= vectorial_variables_tmp;
          grid= grid_tmp;

          loaded=false;
          interpolated=false;

          averaging_period=0; // hours
          time= time_tmp;
  #ifdef OASIS
          coupled = false;
  #endif
        }

    }
    else if (strcmp (DatasetName, "nesting_dynamics_elements") == 0)
    {
        // Definition of nesting grid and datasets
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable sigma1={
            filename_string: "", // All variables are in the same (grid) file
            name: "stress_tensor_first_component",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "1",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable sigma2={
            filename_string: "", // All variables are in the same (grid) file
            name: "stress_tensor_second_component",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "1",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable sigma3={
            filename_string: "", // All variables are in the same (grid) file
            name: "stress_tensor_third_component",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "1",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable damage={
            filename_string: "", // All variables are in the same (grid) file
            name: "sea_ice_damage",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "1",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable ridge_ratio={
            filename_string: "", // All variables are in the same (grid) file
            name: "ridge_ratio",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "1",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            //interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type: BilinearInterpEnum,
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "nesting_"
                   +Environment::vm()["nesting.outer_mesh"].as<std::string>()
                   +"_%Y%m%d.nc",
            gridfile: "nesting_grid_"
                    +Environment::vm()["nesting.inner_mesh"].as<std::string>()
                    +".nc",
            reference_date: "1900-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: damage
        };

        std::vector<Variable> variables_tmp(5);
        variables_tmp[0] = sigma1;
        variables_tmp[1] = sigma2;
        variables_tmp[2] = sigma3;
        variables_tmp[3] = damage;
        variables_tmp[4] = ridge_ratio;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0; // hours
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "nesting_nodes") == 0)
    {
        // Definition of nesting grid and datasets
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable siu={
            filename_string: "", // All variables are in the same (grid) file
            name: "sea_ice_x_velocity",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m s-1",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable siv={
            filename_string: "", // All variables are in the same (grid) file
            name: "sea_ice_y_velocity",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m s-1",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            //interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type: BilinearInterpEnum,
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "nesting_"
                   +Environment::vm()["nesting.outer_mesh"].as<std::string>()
                   +"_%Y%m%d.nc",
            gridfile: "nesting_grid_"
                    +Environment::vm()["nesting.inner_mesh"].as<std::string>()
                    +".nc",
            reference_date: "1900-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: siu
        };

        std::vector<Variable> variables_tmp(2);
        variables_tmp[0] = siu;
        variables_tmp[1] = siv;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0; // hours
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "nesting_distance_nodes") == 0)
    {
        // Definition of nesting grid and datasets
        Dimension dimension_x={
            name:"x",
            cyclic:false
        };

        Dimension dimension_y={
            name:"y",
            cyclic:false
        };

        std::vector<Dimension> dimensions(2);
        dimensions[0] = dimension_y;
        dimensions[1] = dimension_x;


        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable nesting_distance={
            filename_string: "", // All variables are in the same (grid) file
            name: "nesting_distance",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "nesting_grid_"
                   +Environment::vm()["nesting.inner_mesh"].as<std::string>()
                   +".nc",
            gridfile: "nesting_grid_"
                    +Environment::vm()["nesting.inner_mesh"].as<std::string>()
                    +".nc",
            reference_date: "",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"constant",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: nesting_distance
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = nesting_distance;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0; // hours
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "nesting_distance_elements") == 0)
    {
        // Definition of nesting grid and datasets
        Dimension dimension_x={
            name:"x",
            cyclic:false
        };

        Dimension dimension_y={
            name:"y",
            cyclic:false
        };

        std::vector<Dimension> dimensions(2);
        dimensions[0] = dimension_y;
        dimensions[1] = dimension_x;


        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable nesting_distance={
            filename_string: "", // All variables are in the same (grid) file
            name: "nesting_distance",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "nesting_grid_"
                   +Environment::vm()["nesting.inner_mesh"].as<std::string>()
                   +".nc",
            gridfile: "nesting_grid_"
                    +Environment::vm()["nesting.inner_mesh"].as<std::string>()
                    +".nc",
            reference_date: "",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"constant",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: nesting_distance
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = nesting_distance;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0; // hours
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "topaz_forecast_nodes") == 0)
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "u",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "v",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable ssh={
            filename_string: "", // All variables are in the same (grid) file
            name: "ssh",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            //interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interpolation_method: InterpolationType::FromGridToMesh,
            interp_type: -1,
            dirname: "TOPAZ4RC_daily",
            filename_mask: "topaz_forecast_%Y%m%d.nc",
            gridfile: "",
            reference_date: "1950-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: u
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

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "topaz_forecast_elements") == 0)
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable sst={
            filename_string: "", // All variables are in the same (grid) file
            name: "temperature",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "deg celsius",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable sss={
            filename_string: "", // All variables are in the same (grid) file
            name: "salinity",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable mld={
            filename_string: "", // All variables are in the same (grid) file
            name: "mlp",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable conc={
            filename_string: "", // All variables are in the same (grid) file
            name: "fice",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable thick={
            filename_string: "", // All variables are in the same (grid) file
            name: "hice",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable snow_thick={
            filename_string: "", // All variables are in the same (grid) file
            name: "hsnow",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            //interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interpolation_method: InterpolationType::FromGridToMesh,
            interp_type: -1,
            dirname: "TOPAZ4RC_daily",
            filename_mask: "topaz_forecast_%Y%m%d.nc",
            gridfile: "",
            reference_date: "1950-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable:sss
        };

        std::vector<Variable> variables_tmp(6);
        variables_tmp[0] = sst;
        variables_tmp[1] = sss;
        variables_tmp[2] = mld;
        variables_tmp[3] = conc;
        variables_tmp[4] = thick;
        variables_tmp[5] = snow_thick;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable conc={
            filename_string: "", // All variables are in the same (grid) file
            name: "fice",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable thick={
            filename_string: "", // All variables are in the same (grid) file
            name: "hice",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable snow_thick={
            filename_string: "", // All variables are in the same (grid) file
            name: "hsnow",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "TP4DAILY_%Y%m_3m.nc",
            gridfile: "",
            reference_date: "1950-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"monthly",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: conc
        };

        std::vector<Variable> variables_tmp(3);
        variables_tmp[0] = conc;
        variables_tmp[1] = thick;
        variables_tmp[2] = snow_thick;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ice_piomas_elements") == 0)
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable conc={
            filename_string: "", // All variables are in the same (grid) file
            name: "area",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable thick={
            filename_string: "", // All variables are in the same (grid) file
            name: "heff",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable snow_thick={
            filename_string: "", // All variables are in the same (grid) file
            name: "snow",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 1/physical::rhos,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "PIOMAS_%Y.nc",
            gridfile: "",
            reference_date: "1950-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"yearly",
            //monthly_dataset:true,

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: conc
        };

        std::vector<Variable> variables_tmp(3);
        variables_tmp[0] = conc;
        variables_tmp[1] = thick;
        variables_tmp[2] = snow_thick;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=365./12; // days
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ocean_currents_nodes") == 0)
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "U",
            dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: -10000.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "V",
            dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: -10000.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable ssh={
            filename_string: "", // All variables are in the same (grid) file
            name: "SSH",
            dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: -10000.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            //interp_type : NearestInterpEnum,

            dirname:"",
            filename_mask: "current_%Y.nc",
            gridfile: "",
            reference_date: "1950-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"yearly",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: u
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
            east_west_oriented: false
        };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
        vectorial_variables_tmp[0] = uv;

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=365./12; // days
        time= time_tmp;
    }
    else if (strcmp (DatasetName, "ice_amsre_elements") == 0)
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
             name: "longitude",
             dimensions: dimensions_latlon,
             land_mask_defined: false,
             land_mask_value: 0.,
             NaN_mask_defined: false,
             NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
             a: 1.,
             b: 0.,
             Units: "degree_east",
             loaded_data: loaded_data_tmp,
             interpolated_data: interpolated_data_tmp,
             wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1./3600,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable conc={
            filename_string: "", // All variables are in the same (grid) file
            name: "icecon",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 0.01,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "asi-n6250-%Y%m%d-v5i.nc",
            gridfile: "",
            reference_date: "2002-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: conc
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = conc;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ice_osisaf_elements") == 0)
    {
        // Definition of topaz grid and datasets
        Dimension dimension_x={
            name:"xc",
            cyclic:false
        };

        Dimension dimension_y={
            name:"yc",
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
            filename_string: "", // All variables are in the same (grid) file
            name: "lat",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lon",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1./3600,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable conc={
            filename_string: "", // All variables are in the same (grid) file
            name: "ice_conc",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: true,
            NaN_mask_value: -999.,
            use_FillValue: true,
            use_missing_value: true,
            a: 0.01,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable confidence={
            filename_string: "", // All variables are in the same (grid) file
            name: "confidence_level",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "ice_conc_nh_polstere-100_multi_%Y%m%d1200.nc",
            gridfile: "",
            reference_date: "1978-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"nearest_daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: conc
        };

        std::vector<Variable> variables_tmp(2);
        variables_tmp[0] = conc;
        variables_tmp[1] = confidence;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
     }
     else if (strcmp (DatasetName, "ice_osisaf_type_elements") == 0)
     {
        // Definition of topaz grid and datasets
        Dimension dimension_x={
            name:"xc",
            cyclic:false
        };

        Dimension dimension_y={
            name:"yc",
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
            filename_string: "", // All variables are in the same (grid) file
            name: "lat",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lon",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1./3600,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable type={
            filename_string: "", // All variables are in the same (grid) file
            name: "ice_type",
            dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: -1.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "ice_type_nh_polstere-100_multi_%Y%m%d1200.nc",
            gridfile: "",
            reference_date: "1978-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"nearest_daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: type
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = type;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
     }
     else if (strcmp (DatasetName, "ice_smos_elements") == 0)
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable thickness={
            filename_string: "", // All variables are in the same (grid) file
            name: "sea_ice_thickness",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "SMOS_Icethickness_v3.1_north_%Y%m%d.nc",
            gridfile: "",
            reference_date: "2010-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"nearest_daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: thickness
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = thickness;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ice_cs2_smos_elements") == 0)
    {
        // Definition of topaz grid and datasets
        Dimension dimension_x={
            name:"xc",
            cyclic:false
        };

        Dimension dimension_y={
            name:"yc",
            cyclic:false
        };

        Dimension dimension_time={
            name:"time",
            cyclic:false
        };

        std::vector<Dimension> dimensions_latlon = {dimension_y, dimension_x};
        std::vector<Dimension> dimensions = {dimension_time, dimension_y, dimension_x};

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lat",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lon",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable conc={
            filename_string: "", // All variables are in the same (grid) file
            name: "sea_ice_concentration",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 0.01,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable thickness={
            filename_string: "", // All variables are in the same (grid) file
            name: "analysis_sea_ice_thickness",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "cs2_smos_ice_thickness_%Y%m%d-fv2p2.nc",
            gridfile: "",
            reference_date: "",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"nearest_daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: conc
        };

        std::vector<Variable> variables_tmp(2);
        variables_tmp[0] = conc;
        variables_tmp[1] = thickness;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ice_amsr2_elements") == 0)
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable conc={
            filename_string: "", // All variables are in the same (grid) file
            name: "sea_ice_concentration",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 0.01,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable mask={
            filename_string: "", // All variables are in the same (grid) file
            name: "land",
            dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "Arc_%Y%m%d_res3.125_pyres.nc",
            gridfile: "LongitudeLatitudeGrid_3.125km_Arctic.nc",//older AMSR2 files don't have lon,lat inside them
            reference_date: "0001-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"nearest_daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: mask
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = conc;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

         averaging_period=1.; // days
         time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ice_nic_weekly_elements") == 0)
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable conc={
            filename_string: "", // All variables are in the same (grid) file
            name: "sea_ice_concentration",
            dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: 12500.,
            NaN_mask_defined: true,
            NaN_mask_value: 11500.,
            use_FillValue: false,
            use_missing_value: false,
            a: 0.01,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable mask={
            filename_string: "", // All variables are in the same (grid) file
            name: "land",
            dimensions: dimensions_latlon,
            land_mask_defined: true,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "NIC_weekly_%Y%m%d_res3.125_pyres.nc",
            gridfile: "",
            reference_date: "0001-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"nearest_daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: mask
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = conc;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

         averaging_period=1.; // days
         time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
     }
    else if (strcmp (DatasetName, "ice_nic_elements") == 0)
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
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,
            b: 12., // to center the time on the middle of the day
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable conc={
            filename_string: "", // All variables are in the same (grid) file
            name: "sea_ice_concentration",
            dimensions: dimensions,
            land_mask_defined: true,
            land_mask_value: 12500.,
            NaN_mask_defined: true,
            NaN_mask_value: 11500.,
            use_FillValue: false,
            use_missing_value: false,
            a: 0.01,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable mask={
            filename_string: "", // All variables are in the same (grid) file
            name: "land",
            dimensions: dimensions_latlon,
            land_mask_defined: true,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "NIC_%Y%m%d_res3.125_pyres.nc",
            gridfile: "",
            reference_date: "0001-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"nearest_daily",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: mask
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = conc;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

         averaging_period=1.; // days
         time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
     }
     else if (strcmp (DatasetName, "ice_icesat_elements") == 0)
     {
        // Definition of the icesat grid and datasets
         Dimension dimension_x={
             name:"x",
             cyclic:false
        };

         Dimension dimension_y={
             name:"y",
             cyclic:false
        };

        std::vector<Dimension> dimensions_latlon(2);
        dimensions_latlon[0] = dimension_y;
        dimensions_latlon[1] = dimension_x;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "Lat",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "Lon",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable thickness={
            filename_string: "", // All variables are in the same (grid) file
     		name: "Th",
     		dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 0.01,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "icesat_icethk_ON06_filled.nc",
            gridfile: "",
            reference_date: "",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"constant",

            waveOptions: wavopt_none,

            masking: false
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = thickness;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.; // days
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ice_cice_elements") == 0)
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
            filename_string: "", // All variables are in the same (grid) file
            name: "lat",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lon",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable conc={
            filename_string: "", // All variables are in the same (grid) file
            name: "aice_d",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable thick={
            filename_string: "", // All variables are in the same (grid) file
            name: "hi_d",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable snow_thick={
            filename_string: "", // All variables are in the same (grid) file
            name: "hs_d",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "CICE_%Y%m.nc",
            gridfile: "",
            reference_date: "1958-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"monthly",

            waveOptions: wavopt_none,

            masking: true,
            masking_variable: conc
        };

        std::vector<Variable> variables_tmp(3);
        variables_tmp[0] = conc;
        variables_tmp[1] = thick;
        variables_tmp[2] = snow_thick;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=1.; // days
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ice_nemo_elements") == 0)
    {
        // Definition of the icesat grid and datasets
        Dimension dimension_x={
            name:"x",
            cyclic:false
        };

        Dimension dimension_y={
            name:"y",
            cyclic:false
        };

        std::vector<Dimension> dimensions_latlon(2);
        dimensions_latlon[0] = dimension_y;
        dimensions_latlon[1] = dimension_x;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "nav_lat",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "nav_lon",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable thickness={
            filename_string: "", // All variables are in the same (grid) file
            name: "hicif",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable concentration={
            filename_string: "", // All variables are in the same (grid) file
            name: "frld",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: -1.,
            b: 1.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable snow={
            filename_string: "", // All variables are in the same (grid) file
            name: "hsnif",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: false,
            use_missing_value: false,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            interp_type: -1,
            dirname: "",
            filename_mask: "NEMO_icemod.nc",
            gridfile: "",
            reference_date: "",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"constant",

            waveOptions: wavopt_none,

            masking: false
        };

        std::vector<Variable> variables_tmp(3);
        variables_tmp[0] = concentration;
        variables_tmp[1] = thickness;
        variables_tmp[2] = snow;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.; // days
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "dist2coast_elements") == 0)
    {
        // Definition of dist2coast grid and datasets
        Dimension dimension_x={
            name:"lon",
            cyclic:true
        };

        Dimension dimension_y={
            name:"lat",
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
            filename_string: "", // All variables are in the same (grid) file
            name: "lat",
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lon",
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh, //NB can't use _FillValue with this
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"",
            filename_mask: "dist2coast_4deg.nc",
            //filename_mask: "dist2coast_1deg.nc",
            gridfile: "",
            reference_date: "",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,
            dataset_frequency:"constant",

            waveOptions: wavopt_none,

            masking: false
        };

        Variable dist={
            filename_string: "", // All variables are in the same (grid) file
            name:"dist",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1000.,
            b:0.,
            Units:"m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = dist;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "etopo_elements") == 0)
    {
        // Definition of etopo grid and datasets
        Dimension dimension_x={
            name:"lon", // for ETOPO_Arctic_1arcmin.nc
            //name:"x", // for ETOPO1_Ice_g_gmt4.grd
            cyclic:true
        };

        Dimension dimension_y={
            name:"lat", // for ETOPO_Arctic_1arcmin.nc
            //name:"y", // for ETOPO1_Ice_g_gmt4.grd
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
            filename_string: "", // All variables are in the same (grid) file
            name: "lat",  // for ETOPO_Arctic_1arcmin.nc
            //name: "y", // for ETOPO1_Ice_g_gmt4.grd
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lon", // for ETOPO_Arctic_1arcmin.nc
            //name: "x", // for ETOPO1_Ice_g_gmt4.grd
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // Extract the dirname, prefix, and postfix from setup.bathymetry-file (default is ETOPO_Arctic_2arcmin.nc).
        boost::filesystem::path const topo_file ( Environment::vm()["setup.bathymetry-file"].as<std::string>() );

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:topo_file.parent_path().string(),
            filename_mask:topo_file.filename().string(),
            gridfile: "",
            reference_date: "",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,
            dataset_frequency:"constant",

            waveOptions: wavopt_none,

            masking: false
        };

        Variable z={
            filename_string: "", // All variables are in the same (grid) file
            name:"z",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:-1.,
            b:0.,
            Units:"m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        std::vector<Variable> variables_tmp(1);
        variables_tmp[0] = z;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "cfsr_elements") == 0)
    {
        // Definition of dimensions
        Dimension dimension_x={
            name:"lon",
            cyclic:false
        };

        Dimension dimension_y={
            name:"lat",
            cyclic:false
        };

        Dimension dimension_time={
            name:"time0", // "Time"
            cyclic:false
        };

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lat",
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lon",
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"",
            filename_mask: "cfsr.6h.%Y%m.nc",
            gridfile: "",
            reference_date:"1901-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,
            dataset_frequency:"monthly",

            waveOptions: wavopt_none,

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
            filename_string: "", // All variables are in the same (grid) file
            name: "time0",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable tair={
            filename_string: "", // All variables are in the same (grid) file
            name:"TMP_L103",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // T2M
        Variable sphuma={
            filename_string: "", // All variables are in the same (grid) file
            name:"SPF_H_L103",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"kg/kg",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // Q2M
        Variable mslp={
            filename_string: "", // All variables are in the same (grid) file
            name:"PRES_L1",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"Pa",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; //PSFC, a=1.
        Variable Qsw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"DSWRF_L1_Avg_1",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };
        Variable Qlw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"DLWRF_L1_Avg_1",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable precip={
            filename_string: "", // All variables are in the same (grid) file
            name:"A_PCP_L1_Accum_1",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./(6.*3600),
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable snowfr={
            filename_string: "", // All variables are in the same (grid) file
            name:"CSNOW_L1_Avg_1",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        std::vector<Variable> variables_tmp(7);
        variables_tmp[0] = tair;
        variables_tmp[1] = sphuma;
        variables_tmp[2] = mslp;
        variables_tmp[3] = Qsw_in;
        variables_tmp[4] = Qlw_in;
        variables_tmp[5] = precip;
        variables_tmp[6] = snowfr;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "cfsr_nodes") == 0)
    {
        // Definition of dimensions
        Dimension dimension_x={
            name:"lon",
            cyclic:false
        };

        Dimension dimension_y={
            name:"lat",
            cyclic:false
        };

        Dimension dimension_time={
            name:"time0", // "Time"
            cyclic:false};

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lat",
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lon",
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"",
            filename_mask: "cfsr.6h.%Y%m.nc",
            gridfile: "",
            reference_date:"1901-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,

            dataset_frequency:"monthly",

            waveOptions: wavopt_none,

            masking: false
        }   ;

        grid= grid_tmp;

        // Definition of the data

        std::vector<Dimension> dimensions_time(1);
        dimensions_time[0] = dimension_time;

        std::vector<Dimension> dimensions(3);
        dimensions[0] = dimension_time;
        dimensions[1] = dimension_y;
        dimensions[2] = dimension_x;

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time0",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // conversion factors: xnew = a*x + b
        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "U_GRD_L103", // U10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "V_GRD_L103", // V10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
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

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "cfsr_nodes_hi") == 0)
    {
        // Definition of dimensions
        Dimension dimension_x={
            name:"lon",
            cyclic:false
        };

        Dimension dimension_y={
            name:"lat",
            cyclic:false
        };

        Dimension dimension_time={
            name:"time0", // "Time"
            cyclic:false};

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lat",
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lon",
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"",
            filename_mask: "cfsr_h.sh.%Y%m.nc",
            gridfile: "",
            reference_date:"1901-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,

            dataset_frequency:"monthly",

            waveOptions: wavopt_none,

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
            filename_string: "", // All variables are in the same (grid) file
            name: "time0",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // conversion factors: xnew = a*x + b
        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "U_GRD_L103", // U10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "V_GRD_L103", // U10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
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

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ERAi_elements") == 0)
    {
        // Definition of dimensions
        Dimension dimension_x={
             name:"longitude",
             cyclic:true
        };

        Dimension dimension_y={
             name:"latitude",
             cyclic:false
        };

        Dimension dimension_time={
             name:"time",
             cyclic:false
        };

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"",
            filename_mask: "erai.6h.%Y%m.nc",
            gridfile: "",
            reference_date:"1900-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,
            dataset_frequency:"monthly",

            waveOptions: wavopt_none,

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
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable tair={
            filename_string: "", // All variables are in the same (grid) file
            name:"t2m",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // T2M
        Variable dair={
            filename_string: "", // All variables are in the same (grid) file
            name:"d2m",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // D2M
        Variable mslp={
            filename_string: "", // All variables are in the same (grid) file
            name:"msl",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"Pa",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; //PSFC, a=1.

        Variable Qsw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"ssrd",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./(6.*3600),
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable tcc={
            filename_string: "", // All variables are in the same (grid) file
            name:"tcc",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable precip={
            filename_string: "", // All variables are in the same (grid) file
            name:"tp",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:physical::rhow/(6.*3600),
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable snowfall={
            filename_string: "", // All variables are in the same (grid) file
            name:"sf",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:physical::rhow/(6.*3600),
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp
        };

        std::vector<Variable> variables_tmp(7);
        variables_tmp[0] = tair;
        variables_tmp[1] = dair;
        variables_tmp[2] = mslp;
        variables_tmp[3] = Qsw_in;
        variables_tmp[4] = tcc;
        variables_tmp[5] = precip;
        variables_tmp[6] = snowfall;

         std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

         variables= variables_tmp;
         vectorial_variables= vectorial_variables_tmp;

         loaded=false;
         interpolated=false;

         averaging_period=0.;
         time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ERAi_nodes") == 0)
    {
        // Definition of dimensions
        Dimension dimension_x={
            name:"longitude",
            cyclic:true
        };

        Dimension dimension_y={
            name:"latitude",
            cyclic:false
        };

        Dimension dimension_time={
            name:"time",
            cyclic:false};

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"",
            filename_mask: "erai.6h.%Y%m.nc",
            gridfile: "",
            reference_date:"1900-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,

            dataset_frequency:"monthly",

            waveOptions: wavopt_none,

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
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // conversion factors: xnew = a*x + b
        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "u10", //U10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "v10", //V10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
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

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ERA5_elements") == 0)
    {
        // Definition of dimensions
        Dimension dimension_x={
             name:"longitude",
             cyclic:true
        };

        Dimension dimension_y={
             name:"latitude",
             cyclic:false
        };

        Dimension dimension_time={
             name:"time",
             cyclic:false
        };

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"",
            filename_mask: "ERA5_${VARSTRING}_y%Y.nc",
            gridfile: "",
            reference_date:"1900-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,
            dataset_frequency:"yearly",

            waveOptions: wavopt_none,

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
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable tair={
            filename_string: "t2m",
            name:"t2m",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // T2M
        Variable dair={
            filename_string: "d2m",
            name:"d2m",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // D2M
        Variable mslp={
            filename_string: "msl",
            name:"msl",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"Pa",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; //PSFC, a=1.

        Variable Qsw_in={
            filename_string: "msdwswrf",
            name:"msdwswrf",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable Qlw_in={
            filename_string: "msdwlwrf",
            name:"msdwlwrf",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable precip={
            filename_string: "mtpr",
            name:"mtpr",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable snowfall={
            filename_string: "msr",
            name:"msr",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp
        };

        std::vector<Variable> variables_tmp(7);
        variables_tmp[0] = tair;
        variables_tmp[1] = dair;
        variables_tmp[2] = mslp;
        variables_tmp[3] = Qsw_in;
        variables_tmp[4] = Qlw_in;
        variables_tmp[5] = precip;
        variables_tmp[6] = snowfall;

         std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

         variables= variables_tmp;
         vectorial_variables= vectorial_variables_tmp;

         loaded=false;
         interpolated=false;

         averaging_period=0.;
         time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ERA5_nodes") == 0)
    {
        // Definition of dimensions
        Dimension dimension_x={
            name:"longitude",
            cyclic:true
        };

        Dimension dimension_y={
            name:"latitude",
            cyclic:false
        };

        Dimension dimension_time={
            name:"time",
            cyclic:false};

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"",
            filename_mask: "ERA5_${VARSTRING}_y%Y.nc",
            gridfile: "",
            reference_date:"1900-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,

            dataset_frequency:"yearly",

            waveOptions: wavopt_none,

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
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // conversion factors: xnew = a*x + b
        Variable u={
            filename_string: "u10",
            name: "u10", //U10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "v10",
            name: "v10", //V10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
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

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ec2_elements") == 0)
    {
        // Definition of dimensions
        Dimension dimension_x={
            name:"lon",
            cyclic:true
        };

        Dimension dimension_y={
            name:"lat",
            cyclic:false
        };

        Dimension dimension_time={
            name:"time", // "Time"
            cyclic:false
        };

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lat",
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lon",
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"ECMWF_forecast_arctic",
            filename_mask: "ecmwf_forecast_%Y%m%d.nc",
            gridfile: "",
            reference_date:"1950-01-01",//"2008-01-01";

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: wavopt_none,

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
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable tair={
            filename_string: "", // All variables are in the same (grid) file
            name:"2T",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // T2M
        Variable dair={
            filename_string: "", // All variables are in the same (grid) file
            name:"2D",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:-273.15,
            Units:"C",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; // Q2M
        Variable mslp={
            filename_string: "", // All variables are in the same (grid) file
            name:"MSL",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"Pa",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        }; //PSFC, a=1.

        Variable Qsw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"SSRD",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./(Environment::vm()["forecast.ec2_time_res_hours"].as<double>()*3600),
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable Qlw_in={
            filename_string: "", // All variables are in the same (grid) file
            name:"STRD",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1./(Environment::vm()["forecast.ec2_time_res_hours"].as<double>()*3600),
            b:0.,
            Units:"W/m^2",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable tcc={
            filename_string: "", // All variables are in the same (grid) file
            name:"TCC",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:1.,
            b:0.,
            Units:"",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable precip={
            filename_string: "", // All variables are in the same (grid) file
            name:"TP",
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a:physical::rhow/(Environment::vm()["forecast.ec2_time_res_hours"].as<double>()*3600),
            b:0.,
            Units:"kg/m^2/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        std::vector<Variable> variables_tmp(6);
        variables_tmp[0] = tair;
        variables_tmp[1] = dair;
        variables_tmp[2] = mslp;
        variables_tmp[3] = precip;
        variables_tmp[4] = Qsw_in;
        if (Environment::vm()["thermo.use_parameterised_long_wave_radiation"].as<bool>())
            variables_tmp[5] = tcc;
        else
            variables_tmp[5] = Qlw_in;

        std::vector<Vectorial_Variable> vectorial_variables_tmp(0);

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        loaded=false;
        interpolated=false;

         averaging_period=0.;
         time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ec2_nodes") == 0)
    {
        // Definition of dimensions
        Dimension dimension_x={
            name:"lon",
            cyclic:true
        };

        Dimension dimension_y={
            name:"lat",
            cyclic:false
        };

        Dimension dimension_time={
            name:"time", // "Time"
            cyclic:false};

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lat",
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "lon",
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"ECMWF_forecast_arctic",
            filename_mask: "ecmwf_forecast_%Y%m%d.nc",
            gridfile: "",
            reference_date:"1950-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,

            dataset_frequency:"daily",

            waveOptions: wavopt_none,

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
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        // conversion factors: xnew = a*x + b
        Variable u={
            filename_string: "", // All variables are in the same (grid) file
            name: "10U", // U10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable v={
            filename_string: "", // All variables are in the same (grid) file
            name: "10V", // U10M
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
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

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;
        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "ww3a_elements") == 0)
    {
        // Definition of WW3 Arctic analysed - grid and datasets
        Dimension dimension_x={
            name:"longitude",
            cyclic:false
        };

        Dimension dimension_y={
            name:"latitude",
            cyclic:false
        };

        Dimension dimension_time={
            name:"time", // "Time"
            cyclic:false
        };

        // Definition of the grid
        std::vector<Dimension> dimensions_latlon(2);
        dimensions_latlon[0] = dimension_y;
        dimensions_latlon[1] = dimension_x;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_latlon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        // Definition of the data

        std::vector<Dimension> dimensions(3);
        dimensions[0] = dimension_time;
        dimensions[1] = dimension_y;
        dimensions[2] = dimension_x;


        std::vector<Dimension> dimensions_time(1);
        dimensions_time[0] = dimension_time;

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 24.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable SWH={
            filename_string: "", // All variables are in the same (grid) file
            name: "hs", // significant height of wind and swell waves
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        Variable FP={
            filename_string: "", // All variables are in the same (grid) file
            name: "fp", // wave peak frequency
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "/s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        Variable MWDx={
            filename_string: "", // All variables are in the same (grid) file
            name: "dir", // wave mean direction (wave_from_direction)
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions:
                {
                    isWavDir:true,
                    waveFrom:true,
                    xComponent:true
                }
        };


        Variable MWDy={
            filename_string: "", // All variables are in the same (grid) file
            name: "dir", // wave mean direction (wave_from_direction)
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions:
                {
                    isWavDir:true,
                    waveFrom:true,
                    xComponent:false
                }
        };


        Variable FICE={
            filename_string: "", // All variables are in the same (grid) file
            name: "ice", // ice concentration
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        Grid grid_tmp={
            interpolation_method: InterpolationType::FromMeshToMesh2dx,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"",
            filename_mask: "SWARP_WW3_ARCTIC-12K_%Y%m%d.nc",
            gridfile: "",
            reference_date:"1990-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily",

            waveOptions: {
               wave_dataset:true,
               use_mwp:false,
               use_ice:true,
               time_interp_option:"step"
            },

            // wave mask is only for land, not ice
            // ie it is not time-dependant
            masking: true,
            masking_variable: SWH
        };


        std::vector<Variable> variables_tmp(5);
        variables_tmp[0] = SWH;
        variables_tmp[1] = FP;
        variables_tmp[2] = MWDx;
        variables_tmp[3] = MWDy;
        variables_tmp[4] = FICE;

        // make a vectorial variable which
        // is a unit vector representing mwd
        // - better for interpolation
        std::vector<int> uv_tmp(2);
            uv_tmp[0] = 2;
            uv_tmp[1] = 3;

        Vectorial_Variable uv={
            components_Id: uv_tmp,
            east_west_oriented: true
        };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
        vectorial_variables_tmp[0] = uv;

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;
        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;

        time= time_tmp;
#ifdef OASIS
        coupled = false;
#endif
    }
    else if (strcmp (DatasetName, "erai_waves_1deg_elements") == 0)
    {
        // ERAI 1deg waves - grid and datasets
        Dimension dimension_x={
            name:"longitude",
            cyclic:true //need to add extra column, so that 1st and last are the same
        };

        Dimension dimension_y={
            name:"latitude",
            cyclic:false
        };

        Dimension dimension_time={
            name:"time", // "Time"
            cyclic:false
        };

        // Definition of the grid
        std::vector<Dimension> dimensions_lon(1);
        dimensions_lon[0] = dimension_x;

        std::vector<Dimension> dimensions_lat(1);
        dimensions_lat[0] = dimension_y;

        std::vector<Dimension> dimensions_latlon(2);
        dimensions_latlon[0] = dimension_y;
        dimensions_latlon[1] = dimension_x;

        Variable latitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "latitude",
            dimensions: dimensions_lat,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable longitude={
            filename_string: "", // All variables are in the same (grid) file
            name: "longitude",
            dimensions: dimensions_lon,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_east",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };



        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"",
            filename_mask: "erai_waves_1deg_%Y.nc",
            gridfile: "",
            reference_date:"1900-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,
            dataset_frequency:"yearly",

            waveOptions: {
                wave_dataset:true,
                use_mwp:true,
                use_ice:false,
                time_interp_option:"step"
            },

            masking: false
        };

        // Definition of the data

        std::vector<Dimension> dimensions(3);
        dimensions[0] = dimension_time;
        dimensions[1] = dimension_y;
        dimensions[2] = dimension_x;


        std::vector<Dimension> dimensions_time(1);
        dimensions_time[0] = dimension_time;

        Variable time_tmp={
            filename_string: "", // All variables are in the same (grid) file
            name: "time",
            dimensions: dimensions_time,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "hours",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };

        Variable SWH={
            filename_string: "", // All variables are in the same (grid) file
            name: "swh", // significant height of wind and swell waves
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "m",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        Variable MWP={
            filename_string: "", // All variables are in the same (grid) file
            name: "mwp", // wave peak period
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "s",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none
        };


        Variable MWDx={
            filename_string: "", // All variables are in the same (grid) file
            name: "mwd", // wave mean direction (wave_from_direction)
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions:
            {
                isWavDir:true,
                waveFrom:true,
                xComponent:true
            }
        };


        Variable MWDy={
            filename_string: "", // All variables are in the same (grid) file
            name: "mwd", // wave mean direction (wave_from_direction)
            dimensions: dimensions,
            land_mask_defined: false,
            land_mask_value: 0.,
            NaN_mask_defined: false,
            NaN_mask_value: 0.,
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions:
                {
                    isWavDir:true,
                    waveFrom:true,
                    xComponent:false
                }
        };

        std::vector<Variable> variables_tmp(4);
        variables_tmp[0] = SWH;
        variables_tmp[1] = MWP;
        variables_tmp[2] = MWDx;
        variables_tmp[3] = MWDy;


        // make a vectorial variable which
        // is a unit vector representing mwd
        // - better for interpolation
        std::vector<int> uv_tmp(2);
            uv_tmp[0] = 2;
            uv_tmp[1] = 3;

        Vectorial_Variable uv={
            components_Id: uv_tmp,
            east_west_oriented: true
        };

        std::vector<Vectorial_Variable> vectorial_variables_tmp(1);
        vectorial_variables_tmp[0] = uv;

        variables= variables_tmp;
        vectorial_variables= vectorial_variables_tmp;

        grid= grid_tmp;

        loaded=false;
        interpolated=false;

        averaging_period=0.;

        time= time_tmp;
    }
    else
    {
        fprintf (stderr, "Dataset: unknown dataset %s\n",DatasetName);
        fprintf (stderr, "valid types are:\n");
        fprintf (stderr, "asr_nodes\n");
        fprintf (stderr, "asr_elements\n");
        fprintf (stderr, "ec2_nodes\n");
        fprintf (stderr, "ec2_elements\n");
        fprintf (stderr, "ec2_arome_nodes\n");
        fprintf (stderr, "ec2_arome_elements\n");
        fprintf (stderr, "ec2_arome_ensemble_nodes\n");
        fprintf (stderr, "ec2_arome_ensemble_elements\n");
        fprintf (stderr, "cfsr_nodes\n");
        fprintf (stderr, "cfsr_elements\n");
        fprintf (stderr, "cfsr_nodes_hi\n");
        fprintf (stderr, "topaz_nodes\n");
        fprintf (stderr, "nesting_ocean_elements\n");
        fprintf (stderr, "nesting_ice_elements\n");
        fprintf (stderr, "nesting_nodes\n");
        fprintf (stderr, "nesting_distance_nodes\n");
        fprintf (stderr, "nesting_distance_elements\n");
        fprintf (stderr, "topaz_elements\n");
        fprintf (stderr, "ocean_cpl_nodes\n");
        fprintf (stderr, "ocean_cpl_elements\n");
        fprintf (stderr, "ice_topaz_elements\n");
        fprintf (stderr, "ice_amsre_elements\n");
        fprintf (stderr, "ice_osisaf_elements\n");
        fprintf (stderr, "ice_osisaf_type_elements\n");
        fprintf (stderr, "ice_amsr2_elements\n");
        fprintf (stderr, "ice_nic_elements\n");
        fprintf (stderr, "ice_nic_weekly_elements\n");
        fprintf (stderr, "ice_piomas_elements\n");
        fprintf (stderr, "dist2coast_elements\n");
        fprintf (stderr, "etopo_elements\n");
        fprintf (stderr, "ERAi_nodes\n");
        fprintf (stderr, "ERAi_elements\n");
        fprintf (stderr, "ww3a_elements\n");
        fprintf (stderr, "erai_waves_1deg_elements\n");
        fprintf (stderr, "ice_cs2_smos_elements\n");
        fprintf (stderr, "ice_smos_elements\n");
        fprintf (stderr, "ocean_currents_nodes\n");
        fprintf (stderr, "wave_cpl_nodes\n");
        fprintf (stderr, "wave_cpl_elements\n");
        throw std::runtime_error(std::string("Unknown dataset: ")+DatasetName);

        //close_Dataset (this);
    }
}

std::string
DataSet::getFilename(double const& current_time) const
{
    if ( current_time < 0 )
        throw std::runtime_error("getFilename: current time < 0");

    std::string fmask = grid.filename_mask;
    boost::replace_all(fmask, "${VARSTRING}", variables[0].filename_string);
    
    std::string const filename = (boost::format( "%1%/%2%/%3%" )
            % Environment::nextsimDataDir().string()
            % grid.dirname
            % datenumToString(current_time, fmask)
            ).str();
    return filename;
}//getFilename


std::string
DataSet::getFilename(double const& current_time,
        int const& jump) const
{
    if ( current_time < 0 )
        throw std::runtime_error("getFilename: current time < 0");

    double ftime;
    this->shiftDates(current_time, jump, ftime);
    return this->getFilename(ftime);
}//getFilename


void
DataSet::shiftDates(double const& current_time, int const& jump, double& ftime) const
{
    ftime = current_time + jump;//daily files

    // can just quit now for these dataset_frequency options
    std::vector<std::string> namelist = {
        "daily", "nearest_daily", "constant", "coupled"};
    if(std::count(namelist.begin(), namelist.end(), grid.dataset_frequency) > 0)
        return;

    auto pdate = datenumToPosixTime(current_time).date();
    // don't care about the day so move day to 1 to avoid problems with
    // leap years etc
    auto ptime = getPosixTime(pdate.year(), pdate.month(), 1);
    if(grid.dataset_frequency=="monthly")
    {
        //jump is in months
        ftime = posixTimeToDatenum(ptime + boost::gregorian::months(jump));
        return;
    }
    if(grid.dataset_frequency=="yearly")
    {
        //jump is in years
        ftime = posixTimeToDatenum(ptime + boost::gregorian::years(jump));
        return;
    }

    throw std::runtime_error(
            "This option for grid_ptr->dataset_frequency is not implemented: "
            + grid.dataset_frequency);
}//shiftDates()


//! get longitude range from the netcdf file
//! * only needed if grid.interpolation_in_latlon = true
//! called by loadGrid
void
DataSet::getLonRange(double &lonmin, double &lonmax, netCDF::NcVar &VLON)
{
    // - check 1st element
    std::vector<size_t> index_lon_start = {0};
    std::vector<size_t> index_lon_count = {1};
    auto LON = this->getNcVarData(VLON, index_lon_start, index_lon_count);
    double lon0 = LON[0];

    // - check last element
    size_t Nlon = grid.dimension_x_count_netcdf;
    index_lon_start = {Nlon - 1};
    LON = this->getNcVarData(VLON, index_lon_start, index_lon_count);
    double lon1 = LON[0];

    // get range
    lonmin = std::min(lon0, lon1);
    lonmax = std::max(lon0, lon1);
}//getLonRange


//! get nc data as a vector, with scale_factor and offset already applied
//! called by getLonRange
//  TODO also use this in eg getLatLonRegularLatLon...
std::vector<double>
DataSet::getNcVarData(netCDF::NcVar &ncvar, std::vector<size_t> const& start, std::vector<size_t> const& count)
{
    // Attributes (scaling and offset)
    // - Need to multiply with scale factor and add offset
    // - these are stored as variable attributes
    netCDF::NcVarAtt att;
    double scale_factor;
    double add_offset;
    scale_factor = 1.;
    try
    {
        att = ncvar.getAtt("scale_factor");
        att.getValues(&scale_factor);
    }
    catch(netCDF::exceptions::NcException& e)
    {}

    add_offset = 0.;
    try
    {
        att = ncvar.getAtt("add_offset");
        att.getValues(&add_offset);
    }
    catch(netCDF::exceptions::NcException& e)
    {}

    // read the raw data
    ASSERT(count.size()>0, "count is empty");
    ASSERT(count.size() == start.size(), "count and start should be the same size");
    int Ndata = std::accumulate(count.begin(), count.end(),
            1, std::multiplies<double>());
    std::vector<double> data(Ndata);
    ncvar.getVar(start, count, &data[0]);

    // apply the scale factor and offset
    for (auto it=data.begin(); it!=data.end(); it++)
        *it = (*it)*scale_factor + add_offset;

    return data;
}//getNcVarData

void
DataSet::loadGrid(mapx_class *mapNextsim, Grid *grid_ptr, double current_time)
{
    // Empty RXY vector is needed for the coupling so that all the domain is considered.
    /* TODO: Check if we really need to do this, or if the coupling can be made
     * more efficient by only considering a sub-domain */
    std::vector<double> RXY(0);
    loadGrid(mapNextsim, grid_ptr, current_time, RXY, RXY);

}//loadGrid

void
DataSet::loadGrid(mapx_class *mapNextsim, Grid *grid_ptr, double current_time, std::vector<double> const& RX_in, std::vector<double> const& RY_in)
{
    /* We make the loaded domain a bit larger to avoid problems.
     * This is now done by adding a "halo" of one grid cell around the grid
     * cells needed. Replaces the old "expansion_factor" approach. */
    int halo_size = 1;

    // Attributes (scaling and offset)
    netCDF::NcVarAtt att;
    double scale_factor;
    double add_offset;

    std::string filename = grid_ptr->gridfile;
    if ( filename!="" )
        filename = (boost::format( "%1%/%2%" )
                % Environment::nextsimDataDir().string()
                % grid_ptr->gridfile
                ).str();
    else
        filename = this->getFilename(current_time);

    LOG(DEBUG)<<"GRID : FILENAME = "<< filename <<"\n";

    //std::cout <<"GRID : READ NETCDF starts\n";
    if ( ! boost::filesystem::exists(filename) )
        throw std::runtime_error("File not found: " + filename);

    netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
    netCDF::NcDim tmpDim;

    // We initially set the size of the data domain to the full dataset domain, it will be reduce to the model domain in a second step.
    tmpDim = dataFile.getDim(grid_ptr->dimension_y.name);
    grid_ptr->dimension_y_count_netcdf  =  tmpDim.getSize();
    grid_ptr->dimension_y_count  = grid_ptr->dimension_y_count_netcdf;
    grid_ptr->dimension_y_start = 0;

    tmpDim = dataFile.getDim(grid_ptr->dimension_x.name);
    grid_ptr->dimension_x_count_netcdf =  tmpDim.getSize();
    grid_ptr->dimension_x_count =  grid_ptr->dimension_x_count_netcdf;
    grid_ptr->dimension_x_start = 0;

    if(grid_ptr->interpolation_in_latlon)
    {
        ASSERT((grid_ptr->latitude.dimensions.size()==1)
                && (grid_ptr->longitude.dimensions.size()==1),
                "lon & lat should be dimensions if Dataset::grid.interpolation_in_latlon = true");
        ASSERT(grid_ptr->interpolation_method == InterpolationType::FromGridToMesh,
                "Dataset::grid.interpolation_method should be FromGridToMesh if Dataset::grid.interpolation_in_latlon = true");
        netCDF::NcVar VLAT = dataFile.getVar(grid_ptr->latitude.name);
        netCDF::NcVar VLON = dataFile.getVar(grid_ptr->longitude.name);

        // first we determine where the branch cut should be from
        // the range in longitude
        double lonmin, lonmax;
        this->getLonRange(lonmin, lonmax, VLON);
        grid_ptr->branch_cut_lon = lonmin;
        LOG(DEBUG)<<name<<": branch_cut_lon = "<<lonmin<<"\n";

        // TODO we could determine what cyclic should be here too, ie
        // grid_ptr.dimension_x.cyclic = (lonmin + 360. != lonmax);
        // although that wouldn't work for a reduced grid but you get the idea

        // We load the full grid
        std::vector<double> LAT(grid_ptr->dimension_y_count);
        std::vector<double> LON(grid_ptr->dimension_x_count);
        this->getLatLonRegularLatLon(&LAT[0],&LON[0],&VLAT,&VLON);

        // Get the proc speciffic boundaries
        double RX_min, RX_max, RY_min, RY_max;
        this->getMinMax(mapNextsim, grid_ptr, RX_in, RY_in, RX_min, RX_max, RY_min, RY_max);

        // Then, we determine the reduced dimension
        int tmp_start, tmp_end;
        this->findMinMaxIndices(LAT, RY_min, RY_max, tmp_start, tmp_end);

        // Add a halo
        this->addHalo(halo_size, tmp_start, tmp_end,
                grid_ptr->dimension_y_start, grid_ptr->dimension_y_count);

        // Do the same for x
        this->findMinMaxIndices(LON, RX_min, RX_max, tmp_start, tmp_end);
        this->addHalo(halo_size, tmp_start, tmp_end,
                grid_ptr->dimension_x_start, grid_ptr->dimension_x_count);

        // Then we load the reduced grid
        LAT.resize(grid_ptr->dimension_y_count);
        LON.resize(grid_ptr->dimension_x_count);
        getLatLonRegularLatLon(&LAT[0],&LON[0],&VLAT,&VLON);
        grid_ptr->gridY=LAT;
        grid_ptr->gridX=LON;

        // Save lon and lat for possible output
        grid_ptr->gridLAT=LAT;
        grid_ptr->gridLON=LON;

        LOG(DEBUG) <<"GRID : READ NETCDF done\n";
    }//end regular lat=lon
    else if(grid_ptr->interpolation_method==InterpolationType::FromGridToMesh)
    {
        // regular x,y grid
        // - interp from grid to mesh
        // - need grid.mpp_file to be correct .mpp file
        netCDF::NcVar VLAT = dataFile.getVar(grid_ptr->latitude.name);
        netCDF::NcVar VLON = dataFile.getVar(grid_ptr->longitude.name);

        // We load the full grid
        std::vector<double> X(grid_ptr->dimension_x_count);
        std::vector<double> Y(grid_ptr->dimension_y_count);

        this->getXYRegularXY(&X[0],&Y[0],&VLAT,&VLON);
        // Get the proc specific boundaries
        double RX_min, RX_max, RY_min, RY_max;
        this->getMinMax(mapNextsim, grid_ptr, RX_in, RY_in, RX_min, RX_max, RY_min, RY_max);

        // Then, we determine the reduced dimension
        int tmp_start, tmp_end;
        this->findMinMaxIndices(Y, RY_min, RY_max, tmp_start, tmp_end);

        // Add a halo
        this->addHalo(halo_size, tmp_start, tmp_end,
                grid_ptr->dimension_y_start, grid_ptr->dimension_y_count);

        // Do the same for x
        this->findMinMaxIndices(X, RX_min, RX_max, tmp_start, tmp_end);
        this->addHalo(halo_size, tmp_start, tmp_end,
                grid_ptr->dimension_x_start, grid_ptr->dimension_x_count);

        // Then we load the reduced grid
        Y.resize(grid_ptr->dimension_y_count);
        X.resize(grid_ptr->dimension_x_count);
        getXYRegularXY(&X[0],&Y[0],&VLAT,&VLON);
        grid_ptr->gridX=X;
        grid_ptr->gridY=Y;

        // LAT/LON are not regular, better not to save them except if needed.
        //grid_ptr->gridLAT=XLAT;
        //grid_ptr->gridLON=YLON;

        //std::cout <<"GRID : READ NETCDF done\n";

    }//end interpolation_method==InterpolationType::FromGridToMesh
    else if(grid_ptr->interpolation_method==InterpolationType::FromMeshToMesh2dx
         || grid_ptr->interpolation_method==InterpolationType::FromMeshToMeshQuick)
    {
        // interpolation_method==InterpolationType::FromMeshToMesh2dx
        // - most general method
        // - project to x,y plane with nextsim .mpp file and do interpolation in x,y space
        netCDF::NcVar VLAT = dataFile.getVar(grid_ptr->latitude.name);
        netCDF::NcVar VLON = dataFile.getVar(grid_ptr->longitude.name);

        // We load the full grid
        std::vector<double> LAT(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
        std::vector<double> LON(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);

        std::vector<double> X(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
        std::vector<double> Y(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);

        this->getXYLatLonFromLatLon(&X[0],&Y[0],&LAT[0],&LON[0],&VLAT,&VLON);

#if 0
        // Then, we determine the reduced dimension
        std::vector<int> tmp_x_start(0);
        std::vector<int> tmp_x_end(0);
        std::vector<int> tmp_y_start(0);
        std::vector<int> tmp_y_end(0);

        std::vector<int> tmp_x_start(grid_ptr->dimension_y_count,-1);
        std::vector<int> tmp_x_end(grid_ptr->dimension_y_count,-1);
        std::vector<int> tmp_y_start(grid_ptr->dimension_x_count,-1);
        std::vector<int> tmp_y_end(grid_ptr->dimension_x_count,-1);

        for (int i=0; i<grid_ptr->dimension_x_count; ++i)
        {
            for (int j=0; j<grid_ptr->dimension_y_count; ++j)
            {
                if((Y[grid_ptr->dimension_x_count*j+i]>=RY_min) && (Y[grid_ptr->dimension_x_count*j+i]<=RY_max))
                {
                    tmp_y_end[i]=j;
                    if(tmp_y_start[i]==-1)
                        tmp_y_start[i]=j;
                }
            }
        }

        for (int i=0; i<grid_ptr->dimension_y_count; ++i)
        {
            for (int j=0; j<grid_ptr->dimension_x_count; ++j)
            {
                if((X[grid_ptr->dimension_x_count*i+j]>=RX_min) && (X[grid_ptr->dimension_x_count*i+j]<=RX_max))
                {
                    tmp_x_end[i]=j;
                    if(tmp_x_start[i]==-1)
                        tmp_x_start[i]=j;
                }
            }
        }

        int tmp_start=*std::min_element(tmp_y_start.begin(),tmp_y_start.end());
        int tmp_end=*std::max_element(tmp_y_end.begin(),tmp_y_end.end());

        grid_ptr->dimension_y_start=tmp_start;
        grid_ptr->dimension_y_count=tmp_end-tmp_start+1;

        tmp_start=*std::min_element(tmp_x_start.begin(),tmp_x_start.end());
        tmp_end=*std::max_element(tmp_x_end.begin(),tmp_x_end.end());

        grid_ptr->dimension_x_start=tmp_start;
        grid_ptr->dimension_x_count=tmp_end-tmp_start+1;
#endif

        // we just store the indices of all the points included in [RY_min, RY_max]

        // Get the proc specific boundaries
        double RX_min, RX_max, RY_min, RY_max;
        this->getMinMax(mapNextsim, grid_ptr, RX_in, RY_in, RX_min, RX_max, RY_min, RY_max);

        std::vector<int> tmp_tmp_x_id(0);
        std::vector<int> tmp_tmp_y_id(0);
        for (int i=0; i<grid_ptr->dimension_x_count; ++i)
        {
            for (int j=0; j<grid_ptr->dimension_y_count; ++j)
            {
                if(
                    (Y[grid_ptr->dimension_x_count*j+i]>=RY_min) &&
                    (Y[grid_ptr->dimension_x_count*j+i]<=RY_max) &&
                    (X[grid_ptr->dimension_x_count*j+i]>=RX_min) &&
                    (X[grid_ptr->dimension_x_count*j+i]<=RX_max) )
                {
                    tmp_tmp_y_id.push_back(j);
                    tmp_tmp_x_id.push_back(i);
                }
            }
        }

        // tmp_start and tmp_end are the start and end without halos
        int tmp_start=*std::min_element(tmp_tmp_y_id.begin(),tmp_tmp_y_id.end());
        int tmp_end=*std::max_element(tmp_tmp_y_id.begin(),tmp_tmp_y_id.end());

        // Add a halo
        this->addHalo(halo_size, tmp_start, tmp_end,
                grid_ptr->dimension_y_start, grid_ptr->dimension_y_count);

        // Do the same for x
        tmp_start=*std::min_element(tmp_tmp_x_id.begin(),tmp_tmp_x_id.end());
        tmp_end=*std::max_element(tmp_tmp_x_id.begin(),tmp_tmp_x_id.end());
        this->addHalo(halo_size, tmp_start, tmp_end,
                grid_ptr->dimension_x_start, grid_ptr->dimension_x_count);

        // Resize and read
        LAT.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
        LON.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);

        X.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
        Y.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);

        // Then we load the reduced grid
        this->getXYLatLonFromLatLon(&X[0],&Y[0],&LAT[0],&LON[0],&VLAT,&VLON);
#if defined OASIS
        // Read in the gridded rotation angle, if requested
        std::vector<double> Theta;
        if(grid_ptr->gridded_rotation_angle)
        {
            netCDF::NcVar VTHETA;

            // Open the datafile
            std::string filename = this->getFilename(current_time);

            LOG(DEBUG)<<"GRID for masking: FILENAME = "<< filename <<"\n";

            //std::cout <<"GRID : READ NETCDF starts\n";
            if ( ! boost::filesystem::exists(filename) )
                throw std::runtime_error("File not found: " + filename);

            netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

            // load the data
            VTHETA = dataFile.getVar(grid_ptr->vector_rotation_variable.name);

            // here we find the start and count index for each dimensions
            int dims = VTHETA.getDimCount();
            if ( dims != grid_ptr->vector_rotation_variable.dimensions.size() )
                throw std::logic_error( "DataSet::loadGrid: Wrong number of dimensions: Got " + std::to_string(dims) +
                        " from " + filename + " and " + std::to_string(grid_ptr->masking_variable.dimensions.size()) +
                        " from the settings in dataset.cpp");

            std::vector<size_t> index_count(dims);
            std::vector<size_t> index_start(dims);

            for (int i=0; i<dims; ++i)
            {
                netCDF::NcDim tmpDim = VTHETA.getDim(i);
                std::string name = tmpDim.getName();
                if ( name == grid_ptr->dimension_x.name )
                {
                    index_start[i] = grid_ptr->dimension_x_start;
                    index_count[i] = grid_ptr->dimension_x_count;
                }
                else if ( name == grid_ptr->dimension_y.name )
                {
                    index_start[i] = grid_ptr->dimension_y_start;
                    index_count[i] = grid_ptr->dimension_y_count;
                }
                else // We take the first slice of the time and/or depth dimension
                {
                    index_start[i] = 0;
                    index_count[i] = 1;
                }
            }

            Theta.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
            VTHETA.getVar(index_start,index_count,&Theta[0]);
        }
#endif

        // Then we apply the masking if activated
        if(grid_ptr->masking)
        {
            netCDF::NcVar VMASK;
            netCDF::NcDim tmpDim;

            // Open the datafile
            std::string filename = this->getFilename(current_time);

            LOG(DEBUG)<<"GRID for masking: FILENAME = "<< filename <<"\n";

            //std::cout <<"GRID : READ NETCDF starts\n";
            if ( ! boost::filesystem::exists(filename) )
                throw std::runtime_error("File not found: " + filename);

            netCDF::NcFile dataFile2(filename, netCDF::NcFile::read);

            // load the data
            VMASK = dataFile2.getVar(grid_ptr->masking_variable.name);

            std::vector<double> data_in;

            std::vector<double> reduced_X;
            std::vector<double> reduced_Y;
            std::vector<double> reduced_LAT;
            std::vector<double> reduced_LON;
            std::vector<int> reduced_nodes_ind;
#if defined OASIS
            std::vector<double> reduced_Theta;
#endif

            // here we find the start and count index for each dimensions
            int dims = VMASK.getDimCount();
            if ( dims != grid_ptr->masking_variable.dimensions.size() )
                throw std::logic_error( "DataSet::getXYLatLonFromLatLon: Wrong number of dimensions: Got " + std::to_string(dims) +
                        " from " + filename + " and " + std::to_string(grid_ptr->masking_variable.dimensions.size()) +
                        " from the settings in dataset.cpp");

            std::vector<size_t> index_count(dims);
            std::vector<size_t> index_start(dims);

            for (int i=0; i<dims; ++i)
            {
                netCDF::NcDim tmpDim = VMASK.getDim(i);
                std::string name = tmpDim.getName();
                if ( name == grid_ptr->dimension_x.name )
                {
                    index_start[i] = grid_ptr->dimension_x_start;
                    index_count[i] = grid_ptr->dimension_x_count;
                }
                else if ( name == grid_ptr->dimension_y.name )
                {
                    index_start[i] = grid_ptr->dimension_y_start;
                    index_count[i] = grid_ptr->dimension_y_count;
                }
                else // We take the first slice of the time and/or depth dimension
                {
                    index_start[i] = 0;
                    index_count[i] = 1;
                }
            }

            data_in.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
            VMASK.getVar(index_start, index_count, &data_in[0]);

            // Read the attributes
            netCDF::NcVarAtt att;

            // Look for FillValue definition
            // TODO: Should this be double?
            int FillValue;
            bool find_FillValue=true;
            if(grid_ptr->masking_variable.use_FillValue)
            {
                try
                {
                    att = VMASK.getAtt("_FillValue");
                    att.getValues(&FillValue);
                }
                catch(netCDF::exceptions::NcException& e)
                {
                    find_FillValue=false;
                }
            }
            else
                find_FillValue=false;

            // Look for missing_value definition
            // TODO: Should this be double?
            int missing_value;
            bool find_missing_value=true;
            if(grid_ptr->masking_variable.use_missing_value)
            {
                try
                {
                    att = VMASK.getAtt("missing_value");
                    att.getValues(&missing_value);
                }
                catch(netCDF::exceptions::NcException& e)
                {
                    find_missing_value=false;
                }
            }
            else
                find_missing_value=false;

#if 0
            // Look for FillValue definition
            double valid_min;
            bool find_valid_min=true;
            try
            {
                att = VMASK.getAtt("valid_min");
                att.getValues(&valid_min);
            }
            catch(netCDF::exceptions::NcException& e)
            {
                find_valid_min=false;
            }

            // Look for missing_value definition
            double valid_max;
            bool find_valid_max=true;
            try
            {
                att = VMASK.getAtt("valid_max");
                att.getValues(&valid_max);
            }
            catch(netCDF::exceptions::NcException& e)
            {
                find_valid_max=false;
            }
#endif
            bool find_land_mask     =grid_ptr->masking_variable.land_mask_defined;
            double land_mask_value  =grid_ptr->masking_variable.land_mask_value;
            bool find_NaN_mask      =grid_ptr->masking_variable.NaN_mask_defined;
            double NaN_mask_value   =grid_ptr->masking_variable.NaN_mask_value;

            double tmp_data;

            for (int i=0; i<grid_ptr->dimension_y_count; ++i)
            {
                for (int j=0; j<grid_ptr->dimension_x_count; ++j)
                {
                    tmp_data=data_in[grid_ptr->dimension_x_count*i+j];
                    if (    (!find_FillValue        || (tmp_data != FillValue)          )&&
                            (!find_missing_value    || (tmp_data != missing_value)      )&&
                            //(!find_valid_min        || (tmp_data >= valid_min)         )&&
                            //(!find_valid_max        || (tmp_data <= valid_max)         )&&
                            (!find_land_mask        || (tmp_data != land_mask_value)    )&&
                            (!find_NaN_mask         || (tmp_data != NaN_mask_value)     )&&
                            (!std::isnan(tmp_data)                                      )   )
                    {
                        reduced_X.push_back(X[grid_ptr->dimension_x_count*i+j]);
                        reduced_Y.push_back(Y[grid_ptr->dimension_x_count*i+j]);
                        reduced_LAT.push_back(LAT[grid_ptr->dimension_x_count*i+j]);
                        reduced_LON.push_back(LON[grid_ptr->dimension_x_count*i+j]);
                        reduced_nodes_ind.push_back(grid_ptr->dimension_x_count*i+j);
#ifdef OASIS
                        if(grid_ptr->gridded_rotation_angle)
                            reduced_Theta.push_back(Theta[grid_ptr->dimension_x_count*i+j]);
#endif
                    }
                }
            }
            grid_ptr->gridX=reduced_X;
            grid_ptr->gridY=reduced_Y;
            grid_ptr->gridLAT=reduced_LAT;
            grid_ptr->gridLON=reduced_LON;
#ifdef OASIS
            grid_ptr->gridTheta=reduced_Theta;
#endif
            grid_ptr->reduced_nodes_ind=reduced_nodes_ind;
        }
        else // no masking of the Filled Value
        {
            grid_ptr->gridX=X;
            grid_ptr->gridY=Y;
            grid_ptr->gridLAT=LAT;
            grid_ptr->gridLON=LON;
#ifdef OASIS
            grid_ptr->gridTheta=Theta;
#endif
        }

        LOG(DEBUG) <<"GRID : Triangulate starts\n";
        int* pfindex;
        BamgTriangulatex(&pfindex,&grid_ptr->pfnels,&grid_ptr->gridX[0],
                &grid_ptr->gridY[0],grid_ptr->gridX.size());
        grid_ptr->pfindex.resize(3*(grid_ptr->pfnels));
        for(int i=0;i<grid_ptr->pfindex.size();i++)
            grid_ptr->pfindex[i] = pfindex[i];
        xDelete<int>(pfindex);
        LOG(DEBUG) <<"GRID : NUMTRIANGLES= "<< grid_ptr->pfnels <<"\n";
        LOG(DEBUG) <<"GRID : Triangulate done\n";

    }//interpolation_method==InterpolationType::FromMeshToMesh2dx
#ifdef OASIS
    else if (grid_ptr->interpolation_method==InterpolationType::ConservativeRemapping)
    {
        netCDF::NcVar VLAT = dataFile.getVar(grid_ptr->latitude.name);
        netCDF::NcVar VLON = dataFile.getVar(grid_ptr->longitude.name);

        // We load the full grid
        std::vector<double> LAT(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
        std::vector<double> LON(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);

        std::vector<double> X(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
        std::vector<double> Y(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);

        this->getXYLatLonFromLatLon(&X[0],&Y[0],&LAT[0],&LON[0],&VLAT,&VLON);

        // Read in the gridded rotation angle, if requested
        std::vector<double> Theta;
        if(grid_ptr->gridded_rotation_angle)
        {
            // load the data
            netCDF::NcVar VTHETA = dataFile.getVar(grid_ptr->vector_rotation_variable.name);

            // here we find the start and count index for each dimensions
            int dims = VTHETA.getDimCount();
            if ( dims != grid_ptr->vector_rotation_variable.dimensions.size() )
                throw std::logic_error( "DataSet::loadGrid: Wrong number of dimensions: Got " + std::to_string(dims) +
                        " from " + filename + " and " + std::to_string(grid_ptr->masking_variable.dimensions.size()) +
                        " from the settings in dataset.cpp");

            std::vector<size_t> index_count(dims);
            std::vector<size_t> index_start(dims);

            for (int i=0; i<dims; ++i)
            {
                netCDF::NcDim tmpDim = VTHETA.getDim(i);
                std::string name = tmpDim.getName();
                if ( name == grid_ptr->dimension_x.name )
                {
                    index_start[i] = grid.dimension_x_start;
                    index_count[i] = grid.dimension_x_count;
                }
                else if ( name == grid_ptr->dimension_y.name )
                {
                    index_start[i] = grid.dimension_y_start;
                    index_count[i] = grid.dimension_y_count;
                }
                else // We take the first slice of the time and/or depth dimension
                {
                    index_start[i] = 0;
                    index_count[i] = 1;
                }
            }

            Theta.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
            VTHETA.getVar(index_start,index_count,&Theta[0]);

        }
        grid_ptr->gridX=X;
        grid_ptr->gridY=Y;
        grid_ptr->gridLAT=LAT;
        grid_ptr->gridLON=LON;
        grid_ptr->gridTheta=Theta;
    }
#endif

    grid_ptr->loaded=true;
}//loadGrid

// Add a halo - keeping in mind that dime_count is the total dimension size at the moment.
void inline
DataSet::addHalo(int const halo_size, int const tmp_start, int const tmp_end, int& dim_start, int& dim_count)
{
        dim_start = std::max(0,tmp_start-halo_size);
        double const end = std::min(dim_count-1, tmp_end+halo_size);
        dim_count = end - dim_start + 1;
}

// Determine the reduced dimension
void inline
DataSet::findMinMaxIndices(std::vector<double>& XY, double const R_min, double const R_max, int& tmp_start, int& tmp_end)
{
        tmp_start = -1;
        tmp_end   = -1;
        for (int i=0; i<(XY.size()); ++i)
        {
            if(XY[i]>=R_min && XY[i]<=R_max)
            {
                tmp_end = i;
                if ( tmp_start == -1 )
                    tmp_start = i;
            }
        }
}

void
DataSet::getLatLonRegularLatLon(double* LAT, double* LON,netCDF::NcVar* VLAT_ptr,netCDF::NcVar* VLON_ptr)
{
    // Attributes (scaling and offset)
    netCDF::NcVarAtt att;
    double scale_factor;
    double add_offset;

    // read in coordinates
    std::vector<size_t> index_x_count(1);
    std::vector<size_t> index_y_count(1);

    std::vector<size_t> index_x_start(1);
    std::vector<size_t> index_y_start(1);

    index_y_start[0] = grid.dimension_y_start;
    index_y_count[0] = grid.dimension_y_count;

    index_x_start[0] = grid.dimension_x_start;
    index_x_count[0] = grid.dimension_x_count;
    VLAT_ptr->getVar(index_y_start,index_y_count,&LAT[0]);
    //std::cout<<"getLatLonRegularLatLon x "<<index_x_start[0]<<","<<index_x_count[0]<<"\n";
    VLON_ptr->getVar(index_x_start,index_x_count,&LON[0]);

    // Need to multiply with scale factor and add offset - these are stored as variable attributes
    scale_factor=1.;
    try
    {
        att = VLAT_ptr->getAtt("scale_factor");
        att.getValues(&scale_factor);
    }
    catch(netCDF::exceptions::NcException& e)
    {}

    add_offset=0.;
    try
    {
        att = VLAT_ptr->getAtt("add_offset");
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
        att = VLON_ptr->getAtt("scale_factor");
        att.getValues(&scale_factor);
    }
    catch(netCDF::exceptions::NcException& e)
    {}

    add_offset=0.;
    try
    {
        att = VLON_ptr->getAtt("add_offset");
        att.getValues(&add_offset);
    }
    catch(netCDF::exceptions::NcException& e)
    {}

    for (int i=0; i<(index_x_count[0]); ++i)
        LON[i]=LON[i]*scale_factor + add_offset;

}//getLatLonRegularLatLon

void
DataSet::getXYRegularXY(double* X, double* Y,netCDF::NcVar* VLAT_ptr,netCDF::NcVar* VLON_ptr)
{
    // Attributes (scaling and offset)
    netCDF::NcVarAtt att;
    double scale_factor;
    double add_offset;

    // read in coordinates
    std::vector<size_t> index_px_count(2);
    std::vector<size_t> index_py_count(2);

    std::vector<size_t> index_px_start(2);
    std::vector<size_t> index_py_start(2);

    // Since the initial grid is actually regular, we can still use FromGridToMesh
    // by only taking the first line and column into account (only used for ASR so far)
    index_py_start[0] = grid.dimension_y_start;
    index_py_start[1] = 0;

    index_py_count[0] = grid.dimension_y_count;
    index_py_count[1] = 1;

    index_px_start[0] = 0;
    index_px_start[1] = grid.dimension_x_start;

    index_px_count[0] = 1;
    index_px_count[1] = grid.dimension_x_count;

    std::vector<double> XLAT(index_px_count[0]*index_px_count[1]);
    std::vector<double> XLON(index_px_count[0]*index_px_count[1]);
    std::vector<double> YLAT(index_py_count[0]*index_py_count[1]);
    std::vector<double> YLON(index_py_count[0]*index_py_count[1]);

    //std::cout <<"GRID : READ NETCDF done\n";

    // Need to multiply with scale factor and add offset - these are stored as variable attributes
    VLAT_ptr->getVar(index_px_start,index_px_count,&XLAT[0]);
    VLON_ptr->getVar(index_px_start,index_px_count,&XLON[0]);

    VLAT_ptr->getVar(index_py_start,index_py_count,&YLAT[0]);
    VLON_ptr->getVar(index_py_start,index_py_count,&YLON[0]);

    // Apply the scale factor and offset if any
    scale_factor=1.;
    try
    {
        att = VLAT_ptr->getAtt("scale_factor");
        att.getValues(&scale_factor);
    }
    catch(netCDF::exceptions::NcException& e)
    {}

    add_offset=0.;
    try
    {
        att = VLAT_ptr->getAtt("add_offset");
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
    mapx_class *map;
    std::string configfile = (boost::format( "%1%/%2%" )
                          % Environment::nextsimMeshDir().string()
                          % grid.mpp_file
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

}//getXYRegularXY


void
DataSet::getLatLonXYVectors(std::vector<double> &LAT,std::vector<double> &LON,
                            std::vector<double> &X,std::vector<double> &Y,
                            mapx_class *mapNextsim)
{
    if((this->grid.latitude.dimensions.size()==1)
            && (this->grid.longitude.dimensions.size()==1))
    {
        //regular lon-lat grid
        // i=y_ind*N+x_ind;
        int M=this->grid.gridLAT.size();
        int N=this->grid.gridLON.size();
        int x_ind,y_ind;
        int MN = M*N;
        LAT.resize(MN);
        LON.resize(MN);
        X.resize(MN);
        Y.resize(MN);

        for (int i=0; i<MN; ++i)
        {
            // i=y_ind*N+x_ind;
            x_ind   = i%N;
            y_ind   = i/N;
            LON[i]  = this->grid.gridLON[x_ind];
            LAT[i]  = this->grid.gridLAT[y_ind];
            forward_mapx(mapNextsim,LAT[i],LON[i],&X[i],&Y[i]);
        }
    }
    else if(this->grid.interpolation_method==InterpolationType::FromGridToMesh)
    {
        //regular x-y grid
        int M=this->grid.gridX.size();
        int N=this->grid.gridY.size();
        int x_ind,y_ind;
        int MN = M*N;
        LAT.resize(MN);
        LON.resize(MN);
        X.resize(MN);
        Y.resize(MN);

        for (int i=0; i<MN; ++i)
        {
            // i=y_ind*N+x_ind;
            x_ind   = i%N;
            y_ind   = i/N;
            X[i] = this->grid.gridX[x_ind];
            Y[i] = this->grid.gridY[y_ind];
            inverse_mapx(mapNextsim,X[i],Y[i],
                    &LAT[i],&LON[i]);
        }
    }
    else
    {
        LON = this->grid.gridLON;
        LAT = this->grid.gridLAT;
        X   = this->grid.gridX;
        Y   = this->grid.gridY;
    }
}//getLatLonXYVectors


void
DataSet::getXYLatLonFromLatLon(double* X, double* Y, double* LAT, double* LON,netCDF::NcVar* VLAT_ptr,netCDF::NcVar* VLON_ptr)
{
    // Attributes (scaling and offset)
    netCDF::NcVarAtt att;
    double scale_factor;
    double add_offset;

    // read in coordinates
    std::vector<size_t> index_count(2);
    std::vector<size_t> index_start(2);

    index_start[0] = grid.dimension_y_start;
    index_start[1] = grid.dimension_x_start;

    index_count[0] = grid.dimension_y_count;
    index_count[1] = grid.dimension_x_count;

    // Need to multiply with scale factor and add offset - these are stored as variable attributes
    VLAT_ptr->getVar(index_start,index_count,&LAT[0]);
    VLON_ptr->getVar(index_start,index_count,&LON[0]);

    // Apply the scale factor and offset if any
    scale_factor=1.;
    try
    {
        att = VLAT_ptr->getAtt("scale_factor");
        att.getValues(&scale_factor);
    }
    catch(netCDF::exceptions::NcException& e)
    {}

    add_offset=0.;
    try
    {
        att = VLAT_ptr->getAtt("add_offset");
        att.getValues(&add_offset);
    }
    catch(netCDF::exceptions::NcException& e)
    {}

    if(add_offset!=0. || scale_factor!=1.)
    {
        for (int i=0; i<(index_count[0]*index_count[1]); ++i)
        {
            LON[i]=LON[i]*scale_factor + add_offset;
            LAT[i]=LAT[i]*scale_factor + add_offset;
        }
    }

    // projection
    mapx_class *map;
    std::string configfile = (boost::format( "%1%/%2%" )
                              % Environment::nextsimMeshDir().string()
                              % grid.mpp_file
                              ).str();

    std::vector<char> str(configfile.begin(), configfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    double x;
    double y;

    for (int i=0; i<index_count[0]; ++i)
    {
        for (int j=0; j<index_count[1]; ++j)
        {
            forward_mapx(map,LAT[index_count[1]*i+j],LON[index_count[1]*i+j],&x,&y);
            X[index_count[1]*i+j]=x;
            Y[index_count[1]*i+j]=y;
        }
    }

    close_mapx(map);
}//getXYLatLonFromLatLon


double
DataSet::thetaInRange(double const& th_, double const& th1, bool const& close_on_right)
{
    //if close_on_right: convert th_ to angle in (th1,th1+360]
    //else: convert th_ to angle in [th1,th1+360)
    double th2, dth, th;
    int njump;

    th2   = th1 + 360.;
    if (th_ < th1)
    {
        dth   = th1 - th_;
        njump = std::ceil(dth/360.);
        th    = th_ + njump*360.;
    }
    else if (th_ > th2)
    {
        dth   = th_ - th2;
        njump = std::ceil(dth/360.);
        th = th_ - njump*360.;
    }
#if 0
    else if (th_ == th2)
    {
        th = th1;
    }
#endif
    else
    {
        th = th_;
    }

    if (close_on_right && abs(th-th1)<1.e-12)
        th = th2;

    return th;
}

#if defined OASIS
void
DataSet::setElementWeights(std::vector<int> const &gridP, std::vector<std::vector<int>> const &triangles, std::vector<std::vector<double>> const &weights)
{
    M_gridP = gridP;
    M_triangles = triangles;
    M_weights = weights;
}

void
DataSet::setNodalWeights(const std::vector<double>& RX, const std::vector<double>& RY)
{
    // One call to set the node weights
    InterpFromMeshToMesh2dx_weights(
          M_areacoord, M_vertex, M_it,
          &grid.pfindex[0],&grid.gridX[0],&grid.gridY[0],
          grid.gridX.size(),grid.pfnels,
          grid.gridX.size(),
          &RX[0], &RY[0], RX.size());

    // And another one to set the element weights for a drop-in-the-bucket interpolation
    /* This one's not currently used and hence commented out - this leaves M_it empty with length = 0
     * InterpFromMeshToMesh2dx_weights(
     *       M_areacoord, M_vertex, M_it,
     *       &grid.pfindex[0],&grid.gridX[0],&grid.gridY[0],
     *       grid.gridX.size(),grid.pfnels,
     *       grid.pfnels,
     *       &RX[0], &RY[0], RX.size()); */
}
#endif

void
DataSet::convertTargetXY(Grid *grid_ptr,
        std::vector<double> const& RX_in,  std::vector<double> const& RY_in,
        std::vector<double> & RX_out, std::vector<double> & RY_out,
        mapx_class *mapNextsim)//(double const& u, double const& v)
{
    double const target_size = RX_in.size();
    assert(target_size==RY_in.size());

    RX_out.resize(target_size);
    RY_out.resize(target_size);

    if(grid_ptr->interpolation_in_latlon)
    {
        double lat, lon;
        for (int i=0; i<target_size; ++i)
        {
            //convert to lon,lat
            inverse_mapx(mapNextsim,RX_in[i],RY_in[i],&lat,&lon);
            RY_out[i]=lat;
            double bc_lon=grid_ptr->branch_cut_lon;
            bool close_on_right=false;
                //if true  make target lon >  bc_lon,<=bc_lon+180
                //if false make target lon >= bc_lon,< bc_lon+180
                //this shouldn't matter here though?
            RX_out[i]=this->thetaInRange(lon,bc_lon,close_on_right);
        }
    }
    else
    {
        double cos_rotangle = std::cos(rotation_angle);
        double sin_rotangle = std::sin(rotation_angle);
        //rotate to coord sys of dataset
        for (int i=0; i<target_size; ++i)
        {
            RX_out[i] =  cos_rotangle*RX_in[i]+sin_rotangle*RY_in[i];
            RY_out[i] = -sin_rotangle*RX_in[i]+cos_rotangle*RY_in[i];
        }
    }
}//convertTargetXY

void
DataSet::getMinMax(mapx_class *mapNextsim, Grid *grid_ptr,
        std::vector<double> const& RX_in, std::vector<double> const& RY_in,
        double& RX_min, double& RX_max, double& RY_min, double& RY_max)
{
    // ---------------------------------
    // Projection of the mesh positions into the coordinate system of the data before the interpolation
    // (either the lat,lon projection or a polar stereographic projection with another rotation angle (for ASR))
    // we should need to that also for the TOPAZ native grid, so that we could use a gridtomesh, now we use the latlon of the TOPAZ grid


    if ( RX_in.size() > 0 && RY_in.size() > 0 )
    {
        std::vector<double> RX,RY;//size set in convertTargetXY
        this->convertTargetXY(grid_ptr, RX_in,RY_in,RX,RY,mapNextsim);

        RX_min=*std::min_element(RX.begin(),RX.end());
        RX_max=*std::max_element(RX.begin(),RX.end());
        RY_min=*std::min_element(RY.begin(),RY.end());
        RY_max=*std::max_element(RY.begin(),RY.end());

        //std::cout <<"RX_min= "<< RX_min << "RX_max= "<< RX_max <<"RY_min= "<< RY_min <<"RY_max= "<< RY_max <<"\n";
    } else {
        RX_min = std::numeric_limits<double>::lowest();
        RX_max = std::numeric_limits<double>::max();
        RY_min = RX_min;
        RY_max = RX_max;
    }
}

} // Nextsim
