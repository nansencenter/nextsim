/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gridoutput.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Thu Aug  4 09:47:27 CEST 2016
 */

#ifndef __GridOutput_H
#define __GridOutput_H 1

#include <environment.hpp>
#include <assert.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/format.hpp>
#include <gmshmesh.hpp>
#include <Bamgx.h>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromMeshToGridx.h>
#include <BamgTriangulatex.h>
#include <netcdf>

/**
 * @class GridOutput
 * @brief Manage output of gridded data
 *
 * @see
 *
 */

namespace Nextsim
{
class GridOutput
{
public:

    //////////////////////////////////////////////////////////////////////
    // Misc. enumerations and type definitions
    //////////////////////////////////////////////////////////////////////

    enum fileLength
    {
        // How much data to put into each output file
        inf     = 0, // Never start a new file
        daily   = 1,
        weekly  = 2,
        monthly = 3,
        yearly  = 4
    };

    enum variableKind
    {
        nodal      =  0,
        elemental  =  1
    };

    typedef struct Grid
    {
        std::string gridFile;
        std::string dirname;
        std::string mpp_file;
        std::string dimNameX;
        std::string dimNameY;
        std::string latName;
        std::string lonName;

        bool loaded;
        std::vector<double> gridLAT;
        std::vector<double> gridLON;
        std::vector<double> gridX;
        std::vector<double> gridY;
    } Grid;

    //////////////////////////////////////////////////////////////////////
    // Output variable definitions
    //////////////////////////////////////////////////////////////////////
    /*
     * We need to define the variable ID and the associated name,
     * longName, stdName and units. Try to follow the CF table of
     * standard names whenever possible
     * (http://cfconventions.org/standard-names.html).
     */

    enum variableID
    {
        lsm         =  0,

        // Prognostic variables
        conc        =  1,
        thick       =  2,
        damage      =  3,
        snow        =  4,
        VT_x        =  5,
        VT_y        =  6,
        tsurf       =  7,
        sst         =  8,
        sss         =  9,
        tsurf_ice   = 10,
        t1          = 11,
        t2          = 12,
        h_thin      = 13,
        hs_thin     = 14,
        conc_thin   = 15,

        // Diagnostic variables
        Qa          = 100,
        Qsw         = 101,
        Qlw         = 102,
        Qsh         = 103,
        Qlh         = 104,
        Qo          = 105,
        delS        = 106
    };

    typedef struct Variable
    {
        Variable(){}

        Variable(variableID id, std::vector<double> dmesh, std::vector<double> dgrid)
            : varID(id), data_mesh(dmesh), data_grid(dgrid)
        {
            switch (varID)
            {
                case (variableID::lsm):
                    name     = "lsm";
                    longName = "Land Sea Mask";
                    stdName  = "land_sea_mask";
                    Units    = "1";
                    break;

                // Prognostic variables
                case (variableID::conc):
                    name     = "sic";
                    longName = "Sea Ice Concentration";
                    stdName  = "sea_ice_area_fraction";
                    Units    = "1";
                    break;
                case (variableID::thick):
                    name     = "sit";
                    longName = "Sea Ice Thickness";
                    stdName  = "sea_ice_thickness";
                    Units    = "m";
                    break;
                case (variableID::damage):
                    name     = "damage";
                    longName = "Sea Ice Damage";
                    stdName  = "sea_ice_damage";
                    Units    = "1";
                    break;
                case (variableID::snow):
                    name     = "snt";
                    longName = "Surface Snow Thickness";
                    stdName  = "surface_snow_thickness";
                    Units    = "m";
                    break;
                case (variableID::VT_x):
                    name     = "siu";
                    longName = "Sea Ice X Velocity";
                    stdName  = "sea_ice_x_velocity";
                    Units    = "m s-1";
                    break;
                case (variableID::VT_y):
                    name     = "siv";
                    longName = "Sea Ice Y Velocity";
                    stdName  = "sea_ice_y_velocity";
                    Units    = "m s-1";
                    break;
                case (variableID::tsurf):
                    name     = "ts";
                    longName = "Surface Temperature";
                    stdName  = "surface_temperature";
                    Units    = "degree_Celsius";
                    // CF cannonical units are K, but we can use C also
                    break;
                case (variableID::sst):
                    name     = "sst";
                    longName = "Sea Surface Temperature";
                    stdName  = "sea_surface_temperature";
                    Units    = "degree_Celsius";
                    // CF cannonical units are K, but we can use C also
                    break;
                case (variableID::sss):
                    name     = "sss";
                    longName = "Sea Surface Salinity";
                    stdName  = "sea_surface_salinity";
                    Units    = "1e-3";
                    break;
                case (variableID::tsurf_ice):
                    name     = "tsi";
                    longName = "Sea Ice Surface Temperature";
                    stdName  = "sea_ice_surface_temperature";
                    Units    = "degree_Celsius";
                    // CF cannonical units are K, but we can use C also
                    break;
                case (variableID::t1):
                    name     = "t1";
                    longName = "Ice Temperature 1";
                    stdName  = "ice_temperature_1";
                    Units    = "degree_Celsius";
                    // CF cannonical units are K, but we can use C also
                    break;
                case (variableID::t2):
                    name     = "t2";
                    longName = "Ice Temperature 2";
                    stdName  = "ice_temperature_2";
                    Units    = "degree_Celsius";
                    // CF cannonical units are K, but we can use C also
                    break;
                case (variableID::conc_thin):
                    name     = "sic_thin";
                    longName = "Thin Ice Concentration";
                    stdName  = "thin_ice_area_fraction";
                    Units    = "1";
                    break;
                case (variableID::h_thin):
                    name     = "sit_thin";
                    longName = "Thin Ice Thickness";
                    stdName  = "thin_ice_thickness";
                    Units    = "m";
                    break;
                case (variableID::hs_thin):
                    name     = "snt_thin";
                    longName = "Surface Snow Thickness on thin ice";
                    stdName  = "surface_snow_thickness_on_thin_ice";
                    Units    = "m";
                    break;

                // Diagnostic variables
                case (variableID::Qa):
                    name     = "hfs";
                    longName = "Surface Upward Heat Flux In Air";
                    stdName  = "surface_upward_heat_flux_in_air";
                    Units    = "W m-2";
                    break;
                case (variableID::Qsw):
                    name     = "rss";
                    longName = "Surface Net Upward Shortwave Flux";
                    stdName  = "surface_net_upward_shortwave_flux";
                    Units    = "W m-2";
                    break;
                case (variableID::Qlw):
                    name     = "rls";
                    longName = "Surface Net Upward Longwave Flux";
                    stdName  = "surface_net_upward_longwave_flux";
                    Units    = "W m-2";
                    break;
                case (variableID::Qsh):
                    name     = "hfss";
                    longName = "Surface Upward Sensible Heat Flux";
                    stdName  = "surface_upward_sensible_heat_flux";
                    Units    = "W m-2";
                    break;
                case (variableID::Qlh):
                    name     = "hfsl";
                    longName = "Surface Upward Latent Heat Flux";
                    stdName  = "surface_upward_latent_heat_flux";
                    Units    = "W m-2";
                    break;
                case (variableID::Qo):
                    name     = "hfos";
                    longName = "Surface Upward Heatflux In Ocean";
                    stdName  = "surface_upward_heatflux_in_ocean";
                    Units    = "W m-2";
                    break;
                case (variableID::delS):
                    name     = "sfo";
                    longName = "Downward Saltflux In Ocean";
                    stdName  = "downward_slatflux_in_ocean";
                    Units    = "kg m-2 s-1";
                    break;
            }
        }

        int varID;
        std::string name;
        std::string longName;
        std::string stdName;
        std::string Units;

        std::vector<double> data_mesh;
        std::vector<double> data_grid;

    } Variable;

    typedef struct Vectorial_Variable
    {
        std::vector<int> components_Id;
        bool east_west_oriented;
    } Vectorial_Variable;

    ///////////////////////////////////////////////////////////////////////
    // Constructors (and destructor)
    ///////////////////////////////////////////////////////////////////////
    GridOutput();

    GridOutput(int ncols, int nrows, double mooring_spacing, double xmin, double ymin, std::vector<Variable> variables, variableKind kind);

    GridOutput(Grid grid, std::vector<Variable> variables, variableKind kind);

    GridOutput(int ncols, int nrows, double mooring_spacing, double xmin, double ymin, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables);

    GridOutput(Grid grid, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables);

    GridOutput(int ncols, int nrows, double mooring_spacin, double xmin, double yming, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables);

    GridOutput(Grid grid, std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables);

    ~GridOutput();

    void updateGridMean(GmshMesh const& mesh);
    void resetGridMean();
    void resetMeshMean(GmshMesh const& mesh);
    std::string initNetCDF(std::string file_prefix, fileLength file_length, double current_time);
    void appendNetCDF(std::string filename, double timestamp);

    std::vector<int> getMask(GmshMesh const &mesh, variableKind kind);

    int M_ncols;
    int M_nrows;
    double M_mooring_spacing;
    int M_grid_size;
    Grid M_grid;
    std::vector<Variable> M_nodal_variables;
    std::vector<Variable> M_elemental_variables;
    std::vector<Vectorial_Variable> M_vectorial_variables;

    double M_miss_val = -1e+14; // Must be smaller than any expected result

private:

    GridOutput(std::vector<Variable> variables, variableKind kind);

    GridOutput(std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables);

    GridOutput(std::vector<Variable> nodal_variables, std::vector<Variable> elemental_variables, std::vector<Vectorial_Variable> vectorial_variables);

    void initRegularGrid(int ncols, int nrows, double mooring_spacing, double xmin, double ymin);

    void initArbitraryGrid(Grid grid);

    void updateGridMeanWorker(GmshMesh const& mesh, int mesh_size, std::vector<Variable>& variables, double miss_val);

    void rotateVectors(GmshMesh const& mesh, Vectorial_Variable const& vectorial_variable, std::vector<Variable>& variables);

    size_t M_nc_step;
};
} // Nextsim
#endif // __GridOutput_H
