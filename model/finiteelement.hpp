/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   finiteelement.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Mon Aug 24 10:40:29 2015
 */

#ifndef __FiniteElement_HPP
#define __FiniteElement_HPP 1

#include <solverpetsc.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/version.hpp>
#include <boost/format.hpp>
#include <BamgConvertMeshx.h>
#include <Bamgx.h>
#include <InterpFromMeshToMesh2dCavities.h>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromGridToMeshx.h>
#include <InterpFromMeshToGridx.h>
#include <gmshmesh.hpp>
#include <graphcsr.hpp>
#include <wimdiscr.hpp>
#include "enums.hpp"
#include <debug.hpp>
#include <omp.h>
#include <externaldata.hpp>

extern "C"
{
#include <mapx.h>
}


namespace Nextsim
{

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

    typedef ExternalData external_data;
    typedef ExternalData::Dataset Dataset;
    typedef ExternalData::Grid Grid;
    typedef ExternalData::Dimension Dimension;
    typedef ExternalData::Variable Variable;
    typedef boost::ptr_vector<external_data> externaldata_ptr_vector;

    typedef Wim::WimDiscr<double> wim_type;

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

    void initMesh(setup::DomainType const& domain_type, setup::MeshType const& mesh_type);
    void initDatasets();
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
    void solve();
    void run();
    void init();
    void step(bool &is_running, int &pcpt);
    void finalise();
    void error();

    void thermo();
    void thermoIce0(int i, double wspeed, double sphuma, double conc, double voli, double vols,
            double &hi, double &hs, double &hi_old, double &Qio, double &del_hi, double &Tsurf);

    Dataset M_asr_nodes_dataset;
    Dataset M_asr_elements_dataset;
    Dataset M_topaz_nodes_dataset;
    Dataset M_topaz_elements_dataset;
	Dataset M_ice_topaz_elements_dataset;
    Dataset M_etopo_elements_dataset;


    Grid M_asr_grid;
    Grid M_topaz_grid;
	Grid M_etopo_grid;

    double minAngles(element_type const& element, mesh_type const& mesh) const;
    double minAngle(mesh_type const& mesh) const;

    double minAngle(mesh_type const& mesh, std::vector<double> const& um, double factor) const;

    bool flip(mesh_type const& mesh, std::vector<double> const& um, double factor) const;

    double resolution(mesh_type const& mesh) const;

    std::vector<double> hminVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh) const;
    std::vector<double> hmaxVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh) const;

    void initBamg();
    void initConstant();
    void forcing();
    void forcingAtmosphere();
    void forcingOcean();
	void bathymetry();

    void initIce();
    void initThermodynamics();
    void initSlabOcean();
    void initDrifter();
    void coriolis();
    void initNFloes();
    void nodesToElements(double const* depth, std::vector<double>& v);

    void PwlInterp2D();
    void importBamg(BamgMesh const* bamg_mesh);
    void initVariables();
    void initModelState();
    void tensors();
    void cohesion();
    void updateVelocity();
    void scalingVelocity();
    void update();
    void exportResults(int step, bool export_mesh = true);

    void writeRestart(int pcpt, int step);
    void readRestart(int &pcpt, int step);

    void nextsimToWim(bool step);
    void wimToNextsim(bool step);
    void clear();

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
    setup::IceType M_ice_type;
    setup::BathymetryType M_bathymetry_type;

    setup::IceCategoryType M_ice_cat_type;
    setup::DrifterType M_drifter_type;
    setup::DomainType M_domain_type;
    setup::MeshType M_mesh_type;

    LogLevel M_log_level;

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

    //std::vector<double> M_element_depth;
    external_data M_element_depth;
    std::vector<double> M_Vair_factor;
    std::vector<double> M_Voce_factor;
    std::vector<double> M_basal_factor;
    std::vector<double> M_water_elements;
    std::vector<double> M_h_thin;
    std::vector<double> M_hs_thin;
    std::vector<double> M_h_ridged_thin_ice;
    std::vector<double> M_h_ridged_thick_ice;

    //externaldata_ptr_vector M_external_data;
    std::vector<external_data*> M_external_data;

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

    // variables needed for coupling with wim
    wim_type wim;
    std::vector<double> M_nfloes;
    std::vector<double> M_dfloe;
    std::vector<double> M_tau;

    std::vector<double> M_icec_grid;
    std::vector<double> M_iceh_grid;
    std::vector<double> M_nfloes_grid;

    std::vector<double> M_taux_grid;
    std::vector<double> M_tauy_grid;

private:

    double nu0;
    double young;
    double rhoi;
    double rhos;
    double days_in_sec;
    double time_init;
    double output_time_step;
    double mooring_output_time_step;
    double restart_time_step;
    double time_step;
    double duration;
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
    bool M_run_wim;

    bool M_use_restart;
    bool M_write_restart;

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
    external_data M_wind;         // Surface wind [m/s]
    external_data M_tair;         // 2 m temperature [C]
    external_data M_mixrat;       // Mixing ratio
    external_data M_mslp;         // Atmospheric pressure [Pa]
    external_data M_Qsw_in;       // Incoming short-wave radiation [W/m2]
    external_data M_Qlw_in;       // Incoming long-wave radiation [W/m2]
    external_data M_precip;       // Total precipitation [m]
    external_data M_snowfr;       // Fraction of precipitation that is snow
    external_data M_dair;         // 2 m dew point [C]

    // Ocean
    external_data M_ocean;        // "Geostrophic" ocean currents [m/s]
    external_data M_ssh;          // Sea surface elevation [m]

    external_data M_ocean_temp;   // Ocean temperature in top layer [C]
    external_data M_ocean_salt;   // Ocean salinity in top layer [C]
    external_data M_mld;          // Mixed-layer depth [m]

    // Drifters
    boost::unordered_map<int, std::array<double,2>> M_drifter; // Drifters are kept in an unordered map containing number and coordinates
    std::fstream M_iabp_file;             // The file we read the IABP buoy data from
    std::fstream M_drifters_out;    // The file we write our simulated drifter positions into

    // Prognostic ice variables
    std::vector<double> M_conc;         // Ice concentration
    std::vector<double> M_thick;        // Effective ice thickness [m]
    std::vector<double> M_damage;       // Ice damage
    std::vector<double> M_snow_thick;   // Effective snow thickness [m]

	// Prognostic slab ocean variables
    std::vector<double> M_sst;          // Sea-surface temperature [C]
    std::vector<double> M_sss;          // Sea-surface salinity [psu]

	// Non-prognostic variables used to speed up the convergence of a non-linear equation in thermodynamics
    std::vector<double> M_tsurf;        // Ice surface temperature [C]
    std::vector<double> M_tsurf_thin;   // Ice surface temperature of thin ice [C]

private:
    // Variables for the moorings

    bool M_use_moorings;
    int M_grid_size, M_ncols, M_nrows;

    std::vector<double> M_conc_mean;    // Mean concentration (on the mesh)
    std::vector<double> M_thick_mean;   // Mean ice thickness (on the mesh)
    std::vector<double> M_snow_thick_mean;  // Mean snow thickness (on the mesh)
    std::vector<double> M_VT_mean;      // Mean velocity (on the mesh)

    std::vector<double> M_conc_grid;    // Mean concentration (on the grid)
    std::vector<double> M_thick_grid;   // Mean ice thickness (on the grid)
    std::vector<double> M_snow_thick_grid;  // Mean snow thickness (on the grid)
    std::vector<double> M_VT_grid;      // Mean velocity (on the grid)

private:
    void constantIce();
    void targetIce();
    void topazIce();

    void equallySpacedDrifter();
    void outputDrifter(std::fstream &iabp_out);
    void initIABPDrifter();
    void updateIABPDrifter();

    void updateMeans();
    int initMoorings(int &ncols, int &nrows);
    void updateMoorings(int grid_size, int ncols, int nrows);
    void exportMoorings(int grid_size);

#if 0
    void topazIce();
	void etopoBathymetry();
#endif
};
} // Nextsim
#endif
