/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   dataset.hpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#ifndef __DataSet_H
#define __DataSet_H 1

#include <environment.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <Bamgx.h>
#include <InterpFromMeshToMesh2dCavities.h>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromGridToMeshx.h>
#include <InterpFromMeshToGridx.h>
#include <netcdf>
#include <BamgTriangulatex.h>
extern "C"
{
#include <mapx.h>
}

/**
 * @class DataSet
 * @brief Describe the dataset
 * The structure of this class is inspired from the structure of mapx.h
 *
 * @see
 *
 */

namespace Nextsim
{
    enum InterpolationType
    {
        None = -1,
        FromGridToMesh = 0,
        FromMeshToMesh2dx = 1,
        FromMeshToMesh2dCavities = 2,
    };

    typedef struct WaveOptions
    {
       bool wave_dataset;
       bool use_mwp;
       bool use_ice;
       std::string time_interp_option;
    } WaveOptions;

    typedef struct WaveDirOptions
    {
        bool isWavDir;
        bool waveFrom;
        bool xComponent;
    } WaveDirOptions;

class DataSet
{

public:

    typedef struct Dimension
    {
        // Information on the input data
        std::string name;   // name of the dimension in the input file
        bool cyclic;        // if cyclic, then the first value will also be used for interpolation after the last value
    } Dimesion;

    typedef struct Variable
    {
        // Information on the input data
        std::string name;   // name of the variable in the input file
        std::vector<Dimension> dimensions; // dimensions in the input file

        // Information on the fill values and land mask when not available through the netcdf Attributes
        bool land_mask_defined;
        double land_mask_value;
        bool NaN_mask_defined;
        double NaN_mask_value;

        // Information on the unit transform
        double a;           // scale_factor defined by us to have the data in the units system used by nextsim
        double b;           // add_offset defined by us to have the data in the units system used by nextsim
        std::string Units;  // units used in neXtSIM for this variable

        // Storage of the loaded and interpolated data
        std::vector<std::vector<double>> loaded_data;       // 2 vectors, one for the previous and one for the next data timestep
        std::vector<std::vector<double>> interpolated_data; // 2 vectors, one for the previous and one for the next data timestep

        WaveDirOptions wavDirOptions;
            // to determine if this is the mean wave direction
            // - and if so, wave-from or wave-to convention;
            // store as unit vector rotated onto nextsim grid
    } Variable;

    typedef struct Vectorial_Variable
    {
        std::vector<int> components_Id;
        bool east_west_oriented;
    } Vectorial_Variable;

    //so we don't have to change all the datasets
    //every time we change the wave dataset options
    // - default value of Grid.waveOptions.wave_dataset should be false
    WaveOptions wavopt_none  = {
        wave_dataset: false,
        use_mwp: false,
        use_ice: false,
        time_interp_option:"linear",
    };

    WaveDirOptions wavdiropt_none  = {
        isWavDir: false,
        waveFrom: true,
        xComponent: false
    };


    typedef struct Grid
    {
        InterpolationType interpolation_method;
		int interp_type;
        std::string dirname;
        std::string prefix;
        std::string postfix;
        std::string reference_date;

        Variable latitude;
        Variable longitude;

        Dimension dimension_x;
        Dimension dimension_y;

        std::string mpp_file;
		bool interpolation_in_latlon;
		double branch_cut_lon;//where the discontinuity in lon is (only for if interpolation_in_latlon=true)

        bool loaded;
        std::string dataset_frequency;
        std::string target_location;// "mesh_nodes", "mesh_elements", "wim_grid"

        WaveOptions waveOptions;

		bool masking;
		Variable masking_variable;
		std::vector<int> reduced_nodes_ind;

        int* pfindex;
        int pfnels;
        std::vector<double> gridX;
        std::vector<double> gridY;

        std::vector<double> gridLAT;
        std::vector<double> gridLON;

        int dimension_x_start;
        int dimension_x_count;
        int dimension_x_count_netcdf;
        int dimension_y_start;
        int dimension_y_count;
        int dimension_y_count_netcdf;
    } Grid;

    DataSet();

    DataSet(char const *DatasetName,int target_size);

public:
    Grid grid;

    std::vector<Variable> variables; // vector listing all the available variables, included the components of the vectorial variables
    std::vector<Vectorial_Variable> vectorial_variables; // vectors listing the vectorial variables
    int target_size;
    int nb_forcing_step;
    double rotation_angle;

    bool loaded;
    bool interpolated;

    double averaging_period; // 0. for snapshot, 1. for daily values, 365./12 for monthly values.
    Variable time;

    std::vector<double> ftime_range;

    void loadGrid(Grid *grid, int current_time);

    void loadGrid(Grid *grid, int current_time, double RX_min, double RX_max, double RY_min, double RY_max);

    void getLatLonRegularLatLon(double* LAT, double* LON,
                                  netCDF::NcVar* VLAT_ptr,netCDF::NcVar* VLON_ptr);

    void getXYRegularXY(double* X, double* Y,netCDF::NcVar* VLAT_ptr,netCDF::NcVar* VLON_ptr);

    void getXYLatLonFromLatLon(double* X, double* Y,double* LAT, double* LON, netCDF::NcVar* VLAT_ptr,netCDF::NcVar* VLON_ptr);
    double thetaInRange(double const& th_, double const& th1, bool const& close_on_right=false);

    void getLatLonXYVectors(std::vector<double> &LAT,std::vector<double> &LON,
        std::vector<double> &X,std::vector<double> &Y,mapx_class *mapNextsim);

    // name of the dataSet
    std::string name;
    std::string projfilename;
};

} // Nextsim

#endif // __DataSet_H
