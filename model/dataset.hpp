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
#include <assert.hpp>
#include <debug.hpp>
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
#if defined (OASIS) || defined (ENSEMBLE)
        ConservativeRemapping = 3
#endif
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
        std::string name;   //! name of the dimension in the input file
        bool cyclic;        //! if cyclic, then the first value will also be used for interpolation after the last value
    } Dimension;

    typedef struct Variable
    {
        // Information on the input data
        std::string name;   //! name of the variable in the input file
        std::vector<Dimension> dimensions; //! dimensions in the input file

        // Information on the fill values and land mask when not available through the netcdf Attributes
        bool land_mask_defined;
        double land_mask_value;
        bool NaN_mask_defined;
        double NaN_mask_value;

        bool use_FillValue;
        bool use_missing_value;

        // Information on the unit transform
        double a;           //! scale_factor defined by us to have the data in the units system used by nextsim
        double b;           //! add_offset defined by us to have the data in the units system used by nextsim
        std::string Units;  //! units used in neXtSIM for this variable

        // Storage of the loaded and interpolated data
        std::vector<std::vector<double>> loaded_data;       //! 2 vectors, one for the previous and one for the next data timestep
        std::vector<std::vector<double>> interpolated_data; //! 2 vectors, one for the previous and one for the next data timestep

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
        use_ice: false
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
        std::string gridfile;
        std::string reference_date;

        Variable latitude;
        Variable longitude;

        Dimension dimension_x;
        Dimension dimension_y;

        std::string mpp_file;
        bool interpolation_in_latlon;

        bool loaded;
        std::string dataset_frequency;

        WaveOptions waveOptions;

        bool masking;
        Variable masking_variable;

        // optional variables below here...
#if defined (OASIS) || defined (ENSEMBLE)
        bool gridded_rotation_angle;
        Variable vector_rotation_variable;

        std::vector<double> gridTheta;
#endif

		double branch_cut_lon;
                // where the discontinuity in lon is (only used if interpolation_in_latlon=true,
                // and is now determined automatically in loadGrid)
        std::vector<int> reduced_nodes_ind;

        std::vector<int> pfindex;
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

    DataSet(char const *DatasetName);

public:
    Grid grid;

    std::vector<Variable> variables; // vector listing all the available variables, included the components of the vectorial variables
    std::vector<Vectorial_Variable> vectorial_variables; // vectors listing the vectorial variables
    int nb_forcing_step;
    double rotation_angle;

    bool loaded;
    bool interpolated;

    double averaging_period; // 0. for snapshot, 1. for daily values, 365./12 for monthly values.
    Variable time;

    std::vector<double> ftime_range;
#ifdef OASIS
    std::vector<double> itime_range;
#endif
            
    std::string getFilename(Grid *grid, double init_time, double current_time, int jump=0); 

    void loadGrid(Grid *grid, double init_time, double current_time);

    void loadGrid(Grid *grid, double init_time, double current_time, double RX_min, double RX_max, double RY_min, double RY_max);

    std::vector<double> getNcVarData(netCDF::NcVar &ncvar, std::vector<size_t> const& start, std::vector<size_t> const& count);
    void getLonRange(double &lonmin, double &lonmax, netCDF::NcVar &VLON);
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

#if defined OASIS
    bool coupled;
    std::vector<int> M_cpl_id;

    void setWeights(std::vector<int> const &gridP, std::vector<std::vector<int>> const &triangles, std::vector<std::vector<double>> const &weights);

    std::vector<int> M_gridP;
    std::vector<std::vector<int>> M_triangles;
    std::vector<std::vector<double>> M_weights;
#endif

private:

    LogLevel M_log_level;
    bool M_log_all;
    Communicator M_comm;
};

} // Nextsim

#endif // __DataSet_H
