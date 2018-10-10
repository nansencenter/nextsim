/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   dataset.cpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#include <dataset.hpp>
#include <date.hpp>
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

DataSet::DataSet(char const *DatasetName)
{

    name = std::string(DatasetName);
    projfilename = Environment::vm()["mesh.mppfile"].as<std::string>();

    std::vector<std::vector<double>> loaded_data_tmp;
    loaded_data_tmp.resize(2);

    std::vector<std::vector<double>> interpolated_data_tmp;
    interpolated_data_tmp.resize(2);

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

                dirname:"data",
                prefix: "asr30km.comb.2D.",
                postfix:".nc",
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

            dirname:"data",
            prefix:"asr30km.comb.2D.",
            postfix:".nc",
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
            dirname: "data",
            prefix: "TP4DAILY_",
            postfix: "_30m.nc",
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
            dirname: "data",
            prefix: "TP4DAILY_",
            postfix: "_3m.nc",
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
            dirname: "data",
            prefix: "nesting_"+Environment::vm()["nesting.outer_mesh"].as<std::string>()+"_",
            postfix: ".nc",
            gridfile: "nesting_grid_"+Environment::vm()["nesting.inner_mesh"].as<std::string>()+".nc",
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

            Variable sit_thin={
                name: "thin_ice_thickness",
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

            Variable sic_thin={
                name: "thin_ice_area_fraction",
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

            Variable snt_thin={
                name: "surface_snow_thickness_on_thin_ice",
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
                dirname: "data",
                prefix: "nesting_"+Environment::vm()["nesting.outer_mesh"].as<std::string>()+"_",
                postfix: ".nc",
                gridfile: "nesting_grid_"+Environment::vm()["nesting.inner_mesh"].as<std::string>()+".nc",
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
            variables_tmp[3] = sit_thin;
            variables_tmp[4] = sic_thin;
            variables_tmp[5] = snt_thin;

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
              dirname: "data",
              prefix: "nesting_"+Environment::vm()["nesting.outer_mesh"].as<std::string>()+"_",
              postfix: ".nc",
              gridfile: "nesting_grid_"+Environment::vm()["nesting.inner_mesh"].as<std::string>()+".nc",
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
            dirname: "data",
            prefix: "nesting_"+Environment::vm()["nesting.outer_mesh"].as<std::string>()+"_",
            postfix: ".nc",
            gridfile: "nesting_grid_"+Environment::vm()["nesting.inner_mesh"].as<std::string>()+".nc",
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
            dirname: "data",
            prefix: "nesting_"+Environment::vm()["nesting.outer_mesh"].as<std::string>()+"_",
            postfix: ".nc",
            gridfile: "nesting_grid_"+Environment::vm()["nesting.inner_mesh"].as<std::string>()+".nc",
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
            dirname: "data",
            prefix: "nesting_grid_"+Environment::vm()["nesting.inner_mesh"].as<std::string>()+".nc",
            postfix: "",
            gridfile: "nesting_grid_"+Environment::vm()["nesting.inner_mesh"].as<std::string>()+".nc",
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
            dirname: "data",
            prefix: "nesting_grid_"+Environment::vm()["nesting.inner_mesh"].as<std::string>()+".nc",
            postfix: "",
            gridfile: "nesting_grid_"+Environment::vm()["nesting.inner_mesh"].as<std::string>()+".nc",
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
    else if (strcmp (DatasetName, "topaz_cpl_nodes") == 0)
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
            dirname: "data",
            prefix: "TP4DAILY_",
            postfix: "_30m.nc",
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
    else if (strcmp (DatasetName, "topaz_cpl_elements") == 0)
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
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
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
            dirname: "data",
            prefix: "TP4DAILY_",
            postfix: "_3m.nc",
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
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
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
            dirname: "data",
            prefix: "_dm-metno-MODEL-topaz4-ARC-b",
            postfix: "-fv02.0.nc",
            gridfile: "",
            reference_date: "1950-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily_forecast",

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
            use_FillValue: true,
            use_missing_value: true,
            a: 1.,
            b: 0.,
            Units: "degree_north",
            loaded_data: loaded_data_tmp,
            interpolated_data: interpolated_data_tmp,
            wavDirOptions: wavdiropt_none};

        Variable longitude={
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
            dirname: "data",
            prefix: "_dm-metno-MODEL-topaz4-ARC-b",
            postfix: "-fv02.0.nc",
            gridfile: "",
            reference_date: "1950-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: projfilename,
            interpolation_in_latlon: false,

            loaded: false,
            dataset_frequency:"daily_forecast",

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
            dirname: "data",
            prefix: "TP4DAILY_",
            postfix: "_3m.nc",
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
            dirname: "data",
            prefix: "PIOMAS_",
            postfix: ".nc",
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

                dirname:"data",
                prefix: "current_",
                postfix:".nc",
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
            dirname: "data",
            prefix: "asi-n6250-",
            postfix: "-v5i.nc",
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
            dirname: "data",
            prefix:"ice_conc_nh_polstere-100_multi_",
            postfix: "1200.nc",
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
            dirname: "data",
            prefix:"ice_type_nh_polstere-100_multi_",
            postfix: "1200.nc",
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
            dirname: "data",
            prefix:"SMOS_Icethickness_v2.1_north_",
            postfix: ".nc",
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

        std::vector<Dimension> dimensions_latlon(2);
        dimensions_latlon[0] = dimension_y;
        dimensions_latlon[1] = dimension_x;

        Variable latitude={
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

        Variable conc={
            name: "ice_concentration",
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

        Variable thickness={
            name: "analysis_thickness",
            dimensions: dimensions_latlon,
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
            dirname: "data",
            //prefix:"cs2_smos_ice_thickness_20130121_20130127.nc",
            //postfix: "",
            prefix:"cs2_smos_ice_thickness_",
            postfix: ".nc",
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
            dirname: "data",
            prefix: "Arc_",
            postfix: "_res3.125_pyres.nc",
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
            dirname: "data",
            prefix: "NIC_weekly_",
            postfix: "_res3.125_pyres.nc",
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
            dirname: "data",
            prefix: "NIC_",
            postfix: "_res3.125_pyres.nc",
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
            dirname: "data",
            prefix: "icesat_icethk_ON06",
            postfix: "_filled.nc",
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
            dirname:"data",
            prefix:"dist2coast_4deg.nc",
            //prefix:"dist2coast_1deg.nc",
            postfix:"",
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
            cyclic:false
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

        Grid grid_tmp={
            interpolation_method: InterpolationType::FromGridToMesh,
            //interp_type : TriangleInterpEnum, // slower
            interp_type : BilinearInterpEnum,
            //interp_type : NearestInterpEnum,
            dirname:"data",
            prefix:"ETOPO_Arctic_1arcmin.nc",
            //prefix:"ETOPO1_Ice_g_gmt4.grd",
            postfix:"",
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
            dirname:"data",
            prefix: "cfsr.6h.",
            postfix:".nc",
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
            dirname:"data",
            prefix: "cfsr.6h.",
            postfix:".nc",
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
            dirname:"data",
            prefix: "cfsr_h.3h.",
            postfix:".nc",
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
            dirname:"data",
            prefix: "erai.6h.",
            postfix:".nc",
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
        }; // Q2M

        Variable mslp={
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
            dirname:"data",
            prefix: "erai.6h.",
            postfix:".nc",
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
    else if (strcmp (DatasetName, "ec_elements") == 0)
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
            dirname:"data",
            prefix: "ec_start",
            postfix:".nc",
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
        Variable dair={
            name:"D2M",
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

        Variable tcc={
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

        std::vector<Variable> variables_tmp(4);
        variables_tmp[0] = tair;
        variables_tmp[1] = dair;
        variables_tmp[2] = mslp;
        variables_tmp[3] = tcc;

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
    else if (strcmp (DatasetName, "ec_nodes") == 0)
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
            dirname:"data",
            prefix: "ec_start",
            postfix:".nc",
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
            dirname:"data",
            prefix: "ec2_start",
            postfix:".nc",
            gridfile: "",
            reference_date:"1950-01-01",//"2008-01-01";

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,
            dataset_frequency:"daily_ec2_forecast",

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
            name:"SSRD",
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

        Variable Qlw_in={
            name:"STRD",
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
            name:"TP",
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

        std::vector<Variable> variables_tmp(7);
        variables_tmp[0] = tair;
        variables_tmp[1] = dair;
        variables_tmp[2] = mslp;
        variables_tmp[3] = Qsw_in;
        variables_tmp[4] = Qlw_in;
        variables_tmp[5] = tcc;
        variables_tmp[6] = precip;

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
            dirname:"data",
            prefix: "ec2_start",
            postfix:".nc",
            gridfile: "",
            reference_date:"1950-01-01",

            latitude: latitude,
            longitude: longitude,

            dimension_x: dimension_x,
            dimension_y: dimension_y,

            mpp_file: "",
            interpolation_in_latlon: true,

            loaded: false,

            dataset_frequency:"daily_ec2_forecast",

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
            dirname:"data",
            prefix: "SWARP_WW3_ARCTIC-12K_",
            postfix:".nc",
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
            dirname:"data",
            prefix: "erai_waves_1deg_",
            postfix:".nc",
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
        fprintf (stderr, "ec_nodes\n");
        fprintf (stderr, "ec_elements\n");
        fprintf (stderr, "ec2_nodes\n");
        fprintf (stderr, "ec2_elements\n");
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
        fprintf (stderr, "topaz_cpl_nodes\n");
        fprintf (stderr, "topaz_cpl_elements\n");
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

        //close_Dataset (this);
    }

    ftime_range.resize(2,0.);

}

void
DataSet::loadGrid(Grid *grid_ptr, double init_time, double current_time)
{
    loadGrid(grid_ptr, init_time, current_time, 0., 0., 0., 0.);
}

std::string
DataSet::getFilename(Grid *grid_ptr, double init_time, double current_time,int jump)
{

    bool is_topaz_fc = (grid_ptr->dataset_frequency=="daily_forecast");//topaz forecast
    bool is_ec_fc = ((grid_ptr->prefix).find("start") != std::string::npos);//ECMWF forecasts (ec_[nodes,elements],ec2_[nodes,elements])

    if ( current_time < 0 )
        throw std::runtime_error("getFilename: current time < 0");

    std::string current_timestr="";
    std::string filename="";
    if(is_ec_fc||is_topaz_fc)
    {
        // if current (shifted) time before init_time, reduce the init_time to find a file
        // NB jump is in days
        double inittime = std::floor(init_time);
        if(std::floor(current_time+jump)<inittime)
            inittime = std::floor(current_time + jump);
        std::string init_timestr= to_date_string_yd(inittime);//yyyymmdd

        if(is_ec_fc)
    {
            //get filename
            filename = (boost::format( "%1%/%2%%3%%4%" )
                    % Environment::nextsimDataDir().string()
                    % grid_ptr->prefix
                    % init_timestr
                    % grid_ptr->postfix
                    ).str();
            return filename;
    }
    else
    {

            // also need current time for filename
            current_timestr = to_date_string_yd(std::floor(current_time+jump));//yyyymmdd

            filename = (boost::format( "%1%/%2%%3%%4%%5%" )
                    % Environment::nextsimDataDir().string()
                    % current_timestr
                    % grid_ptr->prefix
                    % init_timestr
                    % grid_ptr->postfix
                    ).str();
            return filename;
    }
    }

    if(grid_ptr->dataset_frequency=="monthly")
    {
        //jump is in months
        current_timestr = to_date_string_ym(current_time);//yyyymm

        std::string myString = current_timestr.substr(4,2);
        std::cout <<"month= "<< myString <<"\n";
        int value_month = atoi(myString.c_str());
        myString = current_timestr.substr(0,4);
        std::cout <<"year= "<< myString <<"\n";
        int value_year = atoi(myString.c_str());

        std::cout <<"value_year= "<< value_year <<"\n";
                        std::cout <<"value_month= "<< value_month <<"\n";

        value_month+=jump;
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
        current_timestr=(boost::format( "%1%%2%" )
                % boost::io::group(std::setw(4), std::setfill('0'), value_year)
                % boost::io::group(std::setw(2), std::setfill('0'), value_month)).str();
    }
    else if(grid_ptr->dataset_frequency=="yearly")
    {
        //jump is in years
        current_timestr = to_date_string_y(std::floor(current_time));//yyyy
        int value_year = atoi(current_timestr.c_str());
        value_year+=jump;
        current_timestr=(boost::format( "%1%" )
                % boost::io::group(std::setw(4), std::setfill('0'), value_year)).str();
    }
    else if(grid_ptr->dataset_frequency=="daily"
            || grid_ptr->dataset_frequency=="nearest_daily")
        current_timestr = to_date_string_yd(std::floor(current_time+jump));//yyyymmdd
    else if(grid_ptr->dataset_frequency=="constant")
        current_timestr = "";
    else
        throw std::runtime_error(
                "This option for grid_ptr->dataset_frequency is not implemented: "
                + grid_ptr->dataset_frequency);

    filename = (boost::format( "%1%/%2%%3%%4%" )
                    % Environment::nextsimDataDir().string()
                    % grid_ptr->prefix
                    % current_timestr
                    % grid_ptr->postfix
                    ).str();

    return filename;
}


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
DataSet::loadGrid(Grid *grid_ptr, double init_time, double current_time, double RX_min, double RX_max, double RY_min, double RY_max)
{
    /* We make the loaded domain a bit larger to avoid problems.
     * This is now done by adding a "halo" of one grid cell around the grid
     * cells needed. Replaces the old "expansion_factor" approach. */
    int halo_size = 1;

    //std::cout <<"RX_min= "<< RX_min << "RX_max= "<< RX_max <<"RY_min= "<< RY_min <<"RY_max= "<< RY_max <<"\n";

    // Attributes (scaling and offset)
    netCDF::NcVarAtt att;
    double scale_factor;
    double add_offset;

    //std::string filename = getFilename(grid_ptr, init_time, current_time);
    std::string filename = grid_ptr->gridfile;
    if ( filename=="" )
        filename = getFilename(grid_ptr, init_time, init_time);
    else
        filename = (boost::format( "%1%/%2%" )
                % Environment::nextsimDataDir().string()
                % grid_ptr->gridfile
                ).str();

    std::cout<<"GRID : FILENAME = "<< filename <<"\n";

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
        //TODO implement LOG(DEBUG) here
        //if(Environment::comm().rank()==0)
        //    std::cout<<name<<": branch_cut_lon = "<<lonmin<<"\n";

        // TODO we could determine what cyclic should be here too, ie
        // grid_ptr.dimension_x.cyclic = (lonmin + 360. != lonmax);
        // although that wouldn't work for a reduced grid but you get the idea

        // We load the full grid
		std::vector<double> LAT(grid_ptr->dimension_y_count);
		std::vector<double> LON(grid_ptr->dimension_x_count);
        getLatLonRegularLatLon(&LAT[0],&LON[0],&VLAT,&VLON);
#if 0
        // Then, we determine the reduced dimension
        int tmp_start=-1;
        int tmp_end=-1;
        for (int i=0; i<(LAT.size()); ++i)
        {
            if(LAT[i]>=RY_min && LAT[i]<=RY_max)
            {
                tmp_end=i;
                if(tmp_start==-1)
                    tmp_start=i;
            }
        }

        grid_ptr->dimension_y_start = std::max(0,tmp_start-halo_size);
        grid_ptr->dimension_y_count = std::min(grid_ptr->dimension_y_count, tmp_end+1-tmp_start+2*halo_size);

        tmp_start=-1;
        tmp_end=-1;
        for (int i=0; i<(LON.size()); ++i)
        {
            if((LON[i]>=RX_min) && (LON[i]<=RX_max))
            {
                tmp_end=i;
                if(tmp_start==-1)
                    tmp_start=i;
            }
        }
        grid_ptr->dimension_x_start = std::max(0,tmp_start-halo_size);
        grid_ptr->dimension_x_count = std::min(grid_ptr->dimension_x_count, tmp_end+1-tmp_start+2*halo_size);

		LAT.resize(grid_ptr->dimension_y_count);
		LON.resize(grid_ptr->dimension_x_count);

        std::cout<<tmp_start<<","<<tmp_end<<","<<tmp_end-tmp_start+1<<"\n";
        // Then we load the reduced grid
        getLatLonRegularLatLon(&LAT[0],&LON[0],&VLAT,&VLON);
#endif
        grid_ptr->gridY=LAT;
        grid_ptr->gridX=LON;

        // Save lon and lat for possible output
        grid_ptr->gridLAT=LAT;
        grid_ptr->gridLON=LON;

		std::cout <<"GRID : READ NETCDF done\n";
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

        getXYRegularXY(&X[0],&Y[0],&VLAT,&VLON);
#if 0
        // Then, we determine the reduced dimension
        int tmp_start=-1;
        int tmp_end=-1;
        for (int i=0; i<(Y.size()); ++i)
        {
            if(Y[i]>=RY_min && Y[i]<=RY_max)
            {
                tmp_end=i;
                if(tmp_start==-1)
                    tmp_start=i;
            }
        }

        grid_ptr->dimension_y_start = std::max(0,tmp_start-halo_size);
        grid_ptr->dimension_y_count = std::min(grid_ptr->dimension_y_count, tmp_end+1-tmp_start+2*halo_size);

        tmp_start=-1;
        tmp_end=-1;
        for (int i=0; i<(X.size()); ++i)
        {
            if((X[i]>=RX_min) && (X[i]<=RX_max))
            {
                tmp_end=i;
                if(tmp_start==-1)
                    tmp_start=i;
            }
        }
        grid_ptr->dimension_x_start = std::max(0,tmp_start-halo_size);
        grid_ptr->dimension_x_count = std::min(grid_ptr->dimension_x_count, tmp_end+1-tmp_start+2*halo_size);

		Y.resize(grid_ptr->dimension_y_count);
		X.resize(grid_ptr->dimension_x_count);

        // Then we load the reduced grid
        getXYRegularXY(&X[0],&Y[0],&VLAT,&VLON);
#endif
        grid_ptr->gridX=X;
        grid_ptr->gridY=Y;

        // LAT/LON are not regular, better not to save them except if needed.
        //grid_ptr->gridLAT=XLAT;
        //grid_ptr->gridLON=YLON;

        //std::cout <<"GRID : READ NETCDF done\n";

	}//end interpolation_method==InterpolationType::FromGridToMesh
	else
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

        getXYLatLonFromLatLon(&X[0],&Y[0],&LAT[0],&LON[0],&VLAT,&VLON);


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

        // Add a halo - keeping in mind that grid_ptr->dimension_(x|y)_count is the total dimension
        // size at the moment.
        grid_ptr->dimension_y_start = std::max(0,tmp_start-halo_size);
        tmp_end = std::min(grid_ptr->dimension_y_count, tmp_end+halo_size);
        grid_ptr->dimension_y_count = tmp_end - grid_ptr->dimension_y_start;

        // Do the same for x
        tmp_start=*std::min_element(tmp_tmp_x_id.begin(),tmp_tmp_x_id.end());
        tmp_end=*std::max_element(tmp_tmp_x_id.begin(),tmp_tmp_x_id.end());

        grid_ptr->dimension_x_start = std::max(0,tmp_start-halo_size);
        tmp_end = std::min(grid_ptr->dimension_x_count, tmp_end+halo_size);
        grid_ptr->dimension_x_count = tmp_end - grid_ptr->dimension_x_start;

        // Resize and read
		LAT.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
		LON.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);

		X.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
		Y.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);

        // Then we load the reduced grid
        getXYLatLonFromLatLon(&X[0],&Y[0],&LAT[0],&LON[0],&VLAT,&VLON);

        // Then we apply the masking if activated
		if(grid_ptr->masking){
			netCDF::NcVar VMASK;
            netCDF::NcDim tmpDim;

            // Open the datafile
            filename = getFilename(grid_ptr, init_time, init_time);

            std::cout<<"GRID for masking: FILENAME = "<< filename <<"\n";

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


            // here we find the start and count index for each dimensions
            auto mask_dims = VMASK.getDims();
        	std::vector<size_t> index_count(mask_dims.size());
            std::vector<size_t> index_start(mask_dims.size());
            for(int k=0; k<mask_dims.size(); ++k)
            {
                auto ncdim = mask_dims[k];
                auto dimension_name = ncdim.getName();

                // dimension_x case
                if (dimension_name == grid_ptr->dimension_x.name)
                {
                    index_start[k] = grid_ptr->dimension_x_start;
                    index_count[k] = grid_ptr->dimension_x_count;
                }
                // dimension_y case
                else if (dimension_name == grid_ptr->dimension_y.name)
                {
                    index_start[k] = grid_ptr->dimension_y_start;
                    index_count[k] = grid_ptr->dimension_y_count;
                }
                // other cases: only time or depth
                // - just take 1st record as time/depth-varying masks are not implemented
                // - see below: data_in just has the size of the spatial dimensions
                else
                {
                    index_start[k] = 0;
                    index_count[k] = 1;
                }
            }

			data_in.resize(grid_ptr->dimension_y_count*grid_ptr->dimension_x_count);
			VMASK.getVar(index_start, index_count, &data_in[0]);

            // Read the attributes
			netCDF::NcVarAtt att;

            // Look for FillValue definition
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
					}
				}
			}
			grid_ptr->gridX=reduced_X;
			grid_ptr->gridY=reduced_Y;
			grid_ptr->gridLAT=reduced_LAT;
			grid_ptr->gridLON=reduced_LON;
			grid_ptr->reduced_nodes_ind=reduced_nodes_ind;
		}
		else // no masking of the Filled Value
		{
			grid_ptr->gridX=X;
			grid_ptr->gridY=Y;
			grid_ptr->gridLAT=LAT;
			grid_ptr->gridLON=LON;
		}

		std::cout <<"GRID : Triangulate starts\n";
        int* pfindex;
        BamgTriangulatex(&pfindex,&grid_ptr->pfnels,&grid_ptr->gridX[0],&grid_ptr->gridY[0],grid_ptr->gridX.size());
        grid_ptr->pfindex.resize(3*(grid_ptr->pfnels));
        for(int i=0;i<grid_ptr->pfindex.size();i++)
            grid_ptr->pfindex[i] = pfindex[i];
        xDelete<int>(pfindex);
		std::cout <<"GRID : NUMTRIANGLES= "<< grid_ptr->pfnels <<"\n";
		std::cout <<"GRID : Triangulate done\n";

	}//interpolation_method==InterpolationType::FromMeshToMesh2dx

    grid_ptr->loaded=true;
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

}

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

	// We the initial grid is actually regular, we can still use FromGridToMesh
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

}

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
}

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
}

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

} // Nextsim
