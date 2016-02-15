/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   finiteelement.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Mon Aug 24 10:40:29 2015
 */

#ifndef __FiniteElement_HPP
#define __FiniteElement_HPP 1

#include <solverpetsc.hpp>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/version.hpp>
#include <boost/format.hpp>
#include <BamgConvertMeshx.h>
#include <BamgTriangulatex.h>
#include <Bamgx.h>
#include <InterpFromMeshToMesh2dCavities.h>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromGridToMeshx.h>
#include <gmshmesh.hpp>
#include <graphcsr.hpp>
#include "enums.hpp"
#include <netcdf>
#include <omp.h>

extern "C"
{
#include <mapx.h>
}


namespace Nextsim
{

struct Dimension
{
    std::string name;
    int start;
    int end;
};

struct Variable
{
    std::string name;
    std::vector<Dimension> dimensions;
    double a;
    double b;
    std::string Units;
    netCDF::NcVar NcVar;
};

struct Dataset
{
    int case_number;
    std::string dirname;
    std::string prefix;
    std::string postfix;
    std::string reference_date;
    int nb_timestep_day;
    Variable latitude;
    Variable longitude;
    Variable time;
    Dimension dimension_x;
    Dimension dimension_y;
    Dimension dimension_time;
    std::vector<Variable> variables;
    int target_size;
    std::string mpp_file;
    double rotation_angle;
    std::vector<double> ftime_range;
    int* pfindex;
    int pfnels;
    std::vector<double> gridX;
    std::vector<double> gridY;
};

class FiniteElement
{
public:

    typedef typename GmshMesh::point_type point_type;
    typedef typename GmshMesh::element_type element_type;

    typedef GmshMesh mesh_type;
    typedef SolverPetsc solver_type;
    typedef boost::shared_ptr<solver_type> solver_ptrtype;
    typedef MatrixPetsc matrix_type;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;
    typedef VectorPetsc vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    typedef GraphCSR graph_type;
    typedef boost::shared_ptr<graph_type> graph_ptrtype;

    typedef struct Dimension
    {
        std::string name;
        int start;
        int end;
    } Dimesion;

    typedef struct Variable
    {
        std::string name;
        std::vector<Dimension> dimensions;
        double a;
        double b;
        std::string Units;
        netCDF::NcVar NcVar;
    } Variable;

    typedef struct Grid
    {
        setup::InterpolationType interpolation_method;
        std::string dirname;
        std::string prefix;
        std::string postfix;

        Variable latitude;
        Variable longitude;

        Dimension dimension_x;
        Dimension dimension_y;

        std::string mpp_file;
        double rotation_angle;
        int* pfindex;
        int pfnels;
        std::vector<double> gridX;
        std::vector<double> gridY;
    } Grid;

    typedef struct Dataset
    {
        int case_number;
        std::string dirname;
        std::string prefix;
        std::string postfix;
        std::string reference_date;
        
        int nb_timestep_day;
        Variable time;
        Dimension dimension_time;
        
        std::vector<Variable> variables;
        int target_size;
        Grid *grid;
        std::vector<double> ftime_range;
    } Dataset;

    FiniteElement();

    mesh_type const& mesh() const {return M_mesh;}

    solver_ptrtype const& solverPtr() const {return M_solver;}
    matrix_ptrtype const& matrixPtr() const {return M_matrix;}
    vector_ptrtype const& rhsPtr() const {return M_vector;}
    vector_ptrtype const& solutionPtr() const {return M_solution;}

    solver_type const& solver() const {return *M_solver;}
    matrix_type const& matrix() const {return *M_matrix;}
    vector_type const& rhs() const {return *M_vector;}
    vector_type const& solution() const {return *M_solution;}

    void init();
    void createGMSHMesh(std::string const& geofilename);
    double jacobian(element_type const& element, mesh_type const& mesh) const;
    double jacobian(element_type const& element, mesh_type const& mesh,
                    std::vector<double> const& um, double factor = 1.) const;
    std::vector<double> sides(element_type const& element, mesh_type const& mesh) const;
    std::vector<double> minMaxSide(mesh_type const& mesh) const;
    //void movedMesh(std::vector<double> const& um, double factor = 0);
    double measure(element_type const& element, mesh_type const& mesh) const;
    double measure(element_type const& element, mesh_type const& mesh,
                   std::vector<double> const& um, double factor = 1.) const;
    std::vector<double> shapeCoeff(element_type const& element, mesh_type const& mesh) const;
    void regrid(bool step = true);
    void adaptMesh();

    void assemble(int cpt);
    void assembleSeq(int cpt);
    void solve();
    void run();
    void error();

    void thermo();
    void thermoDamage(double old_vol, int i);
    void slabOcean(double old_conc, double old_vol, double old_snow_vol, double evap, double Qio, double Qow, double Qdw, double Fdw, int i);
    void thermoOW(double &hi, double &hs, double &Qow, double &evap, double wspeed, double del_hi, int i);
    void thermoIce0(double &hi, double &hs, double &Qio, double &del_hi, double &Qai, double wspeed, int i);
    double getAlbedo( double hs, double alb_ice, double alb_sn, double I_0, int alb_scheme );
    double getAlbedoCCSM3( double tsurf, double hs, double alb_ice, double alb_sn );
    void bulkIce(double tsurf, double tair, double mslp, double wspeed, double mixrat, double dair, double Cd, double &Qout, double &dQoutdT, double &subl);
    void openWaterFlux(double &Qow, double &evap, double wspeed, double const ocean_albedo, double const drag_ocean_t, double const drag_ocean_q, int i);
    double calcSphumA(double mslp, double dair, double mixrat);
    void calcSphumI(double mslp, double tsurf, double &sphumi, double &dsphumidT);
    double calcSphumW(double mslp, double sst, double sss);

    Dataset M_asr_nodes_dataset;
    Dataset M_asr_elements_dataset;
    Dataset M_topaz_nodes_dataset;
    Dataset M_topaz_elements_dataset;

    Grid M_asr_grid;
    Grid M_topaz_grid;

    double minAngles(element_type const& element, mesh_type const& mesh) const;
    double minAngle(mesh_type const& mesh) const;

    double minAngle(mesh_type const& mesh, std::vector<double> const& um, double factor) const;

    bool flip(mesh_type const& mesh, std::vector<double> const& um, double factor) const;

    double resolution(mesh_type const& mesh) const;

    std::vector<double> hminVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh) const;
    std::vector<double> hmaxVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh) const;

    std::vector<double> latLon2XY(double const& lat, double const& lon, mapx_class* map, std::string const& configfile);
    double latLon2X(double const& lat, double const& lon, mapx_class* map, std::string const& configfile);
    double latLon2Y(double const& lat, double const& lon, mapx_class* map, std::string const& configfile);

    void initBamg();
    void initConstant();
    void forcing();
    void forcingAtmosphere(bool reload);//(double const& u, double const& v);
    void forcingOcean(bool reload);//(double const& u, double const& v);

    void initConcentration();
    void initThickness();
    void initDamage();
    void initSnowThickness();
    void initThermodynamics();
    void initSlabOcean();
    void initDrifter();
    void coriolis();
    void bathymetry();
    void timeInterpolation(int step);
    void nodesToElements(double const* depth, std::vector<double>& v);

    void PwlInterp2D();
    void importBamg(BamgMesh const* bamg_mesh);
    void initSimulation();
    void tensors();
    void cohesion();
    void updateVelocity();
    void scalingVelocity();
    void update();
    void updateSeq();
    void exportResults(int step, bool export_mesh = true);

private:
    po::variables_map vm;
    mesh_type M_mesh;
    graph_type M_graph;
    mesh_type M_mesh_init;
    mesh_type M_mesh_previous;
    solver_ptrtype M_solver;
    matrix_ptrtype M_matrix;
    //matrix_ptrtype M_mass;
    vector_ptrtype M_vector;
    vector_ptrtype M_solution;
    vector_ptrtype M_exact;

    std::vector<point_type> M_nodes;
    std::vector<element_type> M_edges;
    std::vector<element_type> M_elements;

    int M_num_nodes;
    int M_num_elements;

    std::vector<int> M_boundary_flags;
    std::vector<int> M_dirichlet_flags;
    std::vector<int> M_dirichlet_nodes;
    std::vector<int> M_neumann_flags;
    std::vector<int> M_neumann_nodes;

    boost::mpi::timer chrono, chrono_tot;

    setup::AtmosphereType M_atmosphere_type;
    setup::OceanType M_ocean_type;

    setup::IceCategoryType M_ice_type;
    setup::ConcentrationType M_conc_type;
    setup::ThicknessType M_thick_type;
    setup::DamageType M_damage_type;
    setup::SnowThicknessType M_snow_thick_type;
    setup::DrifterType M_drifter_type;
    setup::DomainType M_domain_type;
    setup::MeshType M_mesh_type;

    std::string M_mesh_filename;

    int M_flag_fix;

    std::vector<double> M_vector_reduction;
    std::vector<bool> M_valid_conc;


    std::vector<double> M_surface;
    std::vector<double> M_sigma;
    std::vector<double> M_divergence_rate;
    //std::vector<double> M_UT;
    std::vector<double> M_UM;
    std::vector<double> M_VT;
    std::vector<double> M_VTM;
    std::vector<double> M_VTMM;

    std::vector<double> M_bathy_depth;

    std::vector<double> M_hminVertices;
    std::vector<double> M_hmaxVertices;

    std::vector<double> M_element_depth;
    std::vector<double> M_Vair_factor;
    std::vector<double> M_Voce_factor;
    std::vector<double> M_basal_factor;
    std::vector<double> M_water_elements;
    std::vector<double> M_ssh;
    std::vector<double> M_h_thin;
    std::vector<double> M_hs_thin;
    std::vector<double> M_h_ridged_thin_ice;
    std::vector<double> M_h_ridged_thick_ice;

    std::vector<double> M_fcor;

    std::vector<double> M_Dunit;
    std::vector<double> M_Dunit_comp;
    std::vector<double> M_Mass;
    std::vector<std::vector<double>> M_shape_coeff;
    std::vector<std::vector<double>> M_B0T;
    std::vector<std::vector<double>> M_B0T_Dunit_B0T;
    std::vector<std::vector<double>> M_B0T_Dunit_comp_B0T;
    std::vector<double> M_random_number;
    std::vector<double> M_Cohesion;
    std::vector<double> M_Compressive_strength;
    std::vector<double> M_time_relaxation_damage;

private:

    double nu0;
    double young;
    double rhoi;
    double rhos;
    double days_in_sec;
    double time_init;
    double output_time_step;
    double time_step;
    double duration;
    double spinup_duration;
    double divergence_min;
    double compression_factor;
    double exponent_compression_factor;
    double ocean_turning_angle_rad;
    double ridging_exponent;
    double quad_drag_coef_air;
    double quad_drag_coef_water;
    double time_relaxation_damage;
    double deltaT_relaxation_damage;
    //double water_depth;
    double ssh_coef;
    double Vair_coef;
    double Voce_coef;

    double basal_k2;
    double basal_drag_coef_air;
    double basal_u_0;
    double basal_Cb;

    double h_thin_max;
    double c_thin_max;

    double compr_strength;
    double tract_coef;
    double scale_coef;
    double alea_factor;
    double cfix;
    double C_fix;
    double C_alea;
    double tan_phi;
    double ridge_h;
    double current_time;
    bool M_reuse_prec;
    bool M_regrid;

private:

    std::vector<std::vector<double>> M_uoce2;
    std::vector<std::vector<double>> M_voce2;
    std::vector<std::vector<double>> M_ssh2;

    std::vector<std::vector<double>> M_uair2;
    std::vector<std::vector<double>> M_vair2;

    std::vector<std::vector<double>> M_tair2;
    std::vector<std::vector<double>> M_mixrat2;
    std::vector<std::vector<double>> M_dair2;
    std::vector<std::vector<double>> M_mslp2;
    std::vector<std::vector<double>> M_Qsw_in2;
    std::vector<std::vector<double>> M_Qlw_in2;
    std::vector<std::vector<double>> M_precip2;
    std::vector<std::vector<double>> M_snowfr2;

    std::vector<std::vector<double>> M_sst2;
    std::vector<std::vector<double>> M_sss2;
    std::vector<std::vector<double>> M_mld2;


private:

    // bamgopt_ptrtype bamgopt;
    // bamgmesh_ptrtype bamgmesh;
    // bamggeom_ptrtype bamggeom;

    // bamgmesh_ptrtype bamgmeshout;
    // bamggeom_ptrtype bamggeomout;

    BamgOpts *bamgopt;
    BamgMesh *bamgmesh;
    BamgGeom *bamggeom;

    // BamgMesh *bamgmeshout;
    // BamgGeom *bamggeomout;

    BamgOpts *bamgopt_previous;
    BamgMesh *bamgmesh_previous;
    BamgGeom *bamggeom_previous;

    //Options *options;

private:

    // Thermodynamic and dynamic forcing
    // Atmosphere
    std::vector<double> M_wind;         // Surface wind [m/s]
    std::vector<double> M_tair;         // 2 m temperature [C]
    std::vector<double> M_mixrat;       // Mixing ratio
    std::vector<double> M_dair;         // 2 m dew point [C]
    std::vector<double> M_mslp;         // Atmospheric pressure [Pa]
    std::vector<double> M_Qsw_in;       // Incoming short-wave radiation [W/m2]
    std::vector<double> M_Qlw_in;       // Incoming long-wave radiation [W/m2]
    std::vector<double> M_precip;       // Total precipitation [m]
    std::vector<double> M_snowfr;       // Fraction of precipitation that is snow
    // Ocean
    std::vector<double> M_ocean;        // "Geostrophic" ocean currents [m/s]
    std::vector<double> M_ocean_temp;   // Ocean temperature in top layer [C]
    std::vector<double> M_ocean_salt;   // Ocean salinity in top layer [C]
    std::vector<double> M_mld;          // Mixed-layer depth [m]

    // Drifters
    std::vector<double> M_drifter;

    // Prognostic variables
    std::vector<double> M_conc;         // Ice concentration
    std::vector<double> M_thick;        // Effective ice thickness [m]
    std::vector<double> M_damage;       // Ice damage
    std::vector<double> M_snow_thick;   // Effective snow thickness [m]
    std::vector<double> M_tsurf;        // Ice surface temperature [C]
    std::vector<double> M_sst;          // Sea-surface temperature [C]
    std::vector<double> M_sss;          // Sea-surface salinity [psu]

private:
    void constantAtmosphere();
    void constantOcean();

    void constantConc();
    void constantThick();
    void constantDamage();
    void constantSnowThick();
    void equallySpacedDrifter();

    void topazConc();
    void topazThick();
    void topazSnowThick();

    void asrAtmosphere(bool reload);//(double const& u, double const& v);
    void topazOcean(bool reload);//(double const& u, double const& v);

    void loadDataset(Dataset *dataset);//(double const& u, double const& v);
    void loadGrid(Grid *grid);
};
} // Nextsim
#endif
