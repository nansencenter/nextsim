/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

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
#include <BamgTriangulatex.h>
#include <Bamgx.h>
#include <InterpFromMeshToMesh2dCavities.h>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromGridToMeshx.h>
#include <gmshmesh.hpp>
#include <gmshmeshseq.hpp>
#include <graphcsr.hpp>
#include <graphcsrmpi.hpp>
#include <externaldata.hpp>
#include <gridoutput.hpp>
#include <dataset.hpp>
#include <drifters.hpp>
#include "enums.hpp"
#include <debug.hpp>
#include <omp.h>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_01.hpp>

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
    typedef typename GmshMesh::bimap_type bimap_type;

    typedef GmshMesh mesh_type;
    typedef GmshMeshSeq mesh_type_root;
    typedef SolverPetsc solver_type;
    typedef boost::shared_ptr<solver_type> solver_ptrtype;
    typedef MatrixPetsc matrix_type;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;
    typedef VectorPetsc vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    typedef GraphCSR graph_type;
    typedef boost::shared_ptr<graph_type> graph_ptrtype;

    typedef GraphCSRMPI graphmpi_type;
    typedef boost::shared_ptr<graphmpi_type> graphmpi_ptrtype;

    // typedef ExternalData external_data;
    // typedef ExternalData::Dataset Dataset;
    // typedef ExternalData::Grid Grid;
    // typedef ExternalData::Dimension Dimension;
    // typedef ExternalData::Variable Variable;
    // typedef ExternalData::Vectorial_Variable Vectorial_Variable;
    // typedef boost::ptr_vector<external_data> externaldata_ptr_vector;

    typedef DataSet Dataset;
    typedef ExternalData external_data;
    typedef typename std::vector<external_data*> external_data_vec ;
    typedef typename std::vector<Dataset*> Dataset_vec ;

    typedef boost::ptr_vector<external_data> externaldata_ptr_vector;

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

    void initMesh();
    void initForcings();
    void createGMSHMesh(std::string const& geofilename);

    double jacobian(element_type const& element, mesh_type const& mesh) const;
    double jacobian(element_type const& element, mesh_type const& mesh,
                    std::vector<double> const& um, double factor = 1.) const;

    double jacobian(element_type const& element, mesh_type_root const& mesh) const;
    double jacobian(element_type const& element, mesh_type_root const& mesh,
                    std::vector<double> const& um, double factor = 1.) const;

    std::vector<double> sides(element_type const& element, mesh_type const& mesh) const;
    std::vector<double> sides(element_type const& element, mesh_type const& mesh,
                              std::vector<double> const& um, double factor) const;

    std::vector<double> sides(element_type const& element, mesh_type_root const& mesh) const;
    std::vector<double> sides(element_type const& element, mesh_type_root const& mesh,
                              std::vector<double> const& um, double factor) const;

    std::vector<double> minMaxSide(mesh_type_root const& mesh) const;

    template<typename FEMeshType>
    double measure(element_type const& element, FEMeshType const& mesh) const;

    template<typename FEMeshType>
    double measure(element_type const& element, FEMeshType const& mesh,
                   std::vector<double> const& um, double factor = 1.) const;

    std::vector<double> shapeCoeff(element_type const& element, mesh_type const& mesh) const;

    std::vector<double> shapeCoeff(element_type const& element, mesh_type_root const& mesh) const;

    void regrid(bool step = true);
    void adaptMesh();
    void updateNodeIds();
    void updateBoundaryFlags();

    void gatherSizes();
    void gatherFieldsElement(std::vector<double>& interp_in_elements);
    void scatterFieldsElement(double* interp_elt_out);

    //void gatherUM(std::vector<double>& um);
    void gatherNodalField(std::vector<double> const& field_local, std::vector<double>& field_root);
    void scatterNodalField(std::vector<double> const& field_root, std::vector<double>& field_local);

    // other interfaces
    void gatherNodalField(std::vector<double> const& field1_local, std::vector<double> const& field2_local,
                          std::vector<double>& field1_root, std::vector<double>& field2_root);

    void gatherElementField(std::vector<double> const& field_local, std::vector<double>& field_root, int nb_fields = 1);
    void scatterElementField(std::vector<double> const& field_root, std::vector<double>& field_local, int nb_fields = 1);

    void gatherFieldsNode(std::vector<double>& interp_in_elements, std::vector<int> const& rmap_nodes, std::vector<int> sizes_nodes);
    void scatterFieldsNode(double* interp_nd_out);

    void interpFields(std::vector<int> const& rmap_nodes, std::vector<int> sizes_nodes);

    void assemble(int pcpt);
    void solve();
    void init();
    void step();
    void run();

    void nestingIce();
    void nestingDynamics();
    void thermo(double dt);
    void thermoIce0(int i, double dt, double wspeed, double sphuma, double conc, double voli, double vols, double Qlw_in, double Qsw_in, double mld, double snowfall,
                    double &hi, double &hs, double &hi_old, double &Qio, double &del_hi, double &Tsurf,
                    double &Qai, double &Qsw, double &Qlw, double &Qsh, double &Qlh);
    void thermoWinton(int i, double dt, double wspeed, double sphuma, double conc, double voli, double vols,
                      double Qlw_in, double Qsw_in, double mld, double snowfall,
                      double &hi, double &hs, double &hi_old, double &Qio, double &del_hi, double &Tsurf, double &T1, double &T2,
                      double &Qai, double &Qsw, double &Qlw, double &Qsh, double &Qlh);
    double albedo(int alb_scheme, double Tsurf, double hs, double alb_sn, double alb_ice, double I_0);
    void atmFluxBulk(int i, double Tsurf, double sphuma, double drag_ice_t, double Qsw, double Qlw_in, double wspeed,
                     double &Qai, double &dQaidT, double &subl,
                     double &Qsh, double &Qlh, double &Qlw);
    double iceOceanHeatflux(int cpt, double sst, double tbot, double mld, double dt);

    void checkReloadDatasets(external_data_vec const& ext_data_vec,
                             double const& CRtime, std::string const& printout);

    Dataset M_atmosphere_nodes_dataset;
    Dataset M_atmosphere_elements_dataset;
    Dataset M_atmosphere_bis_elements_dataset;
    Dataset M_ocean_nodes_dataset;
    Dataset M_ocean_elements_dataset;
    Dataset M_bathymetry_elements_dataset;

    Dataset M_ice_topaz_elements_dataset;
    Dataset M_ice_icesat_elements_dataset;
    Dataset M_ice_piomas_elements_dataset;
    Dataset M_ice_amsre_elements_dataset;
    Dataset M_ice_osisaf_elements_dataset;
    Dataset M_ice_osisaf_type_elements_dataset;
    Dataset M_ice_amsr2_elements_dataset;
    Dataset M_ice_nic_elements_dataset;
    Dataset M_ice_nic_weekly_elements_dataset;
    Dataset M_ice_cs2_smos_elements_dataset;
    Dataset M_ice_smos_elements_dataset;

    // Datasets for nesting from outer domain with coarse resolution
    Dataset M_nesting_ocean_elements_dataset;
    Dataset M_nesting_nodes_dataset;
    Dataset M_nesting_ice_elements_dataset;
    Dataset M_nesting_distance_nodes_dataset;
    Dataset M_nesting_distance_elements_dataset;
    Dataset M_nesting_dynamics_elements_dataset;

    template<typename FEMeshType>
    double minAngles(element_type const& element, FEMeshType const& mesh) const;
    template<typename FEMeshType>
    double minAngles(element_type const& element, FEMeshType const& mesh,
                     std::vector<double> const& um, double factor) const;

    template<typename FEMeshType>
    double minAngle(FEMeshType const& mesh) const;
    template<typename FEMeshType>
    double minAngle(FEMeshType const& mesh, std::vector<double> const& um, double factor, bool root = false) const;

    template<typename FEMeshType>
    bool flip(FEMeshType const& mesh, std::vector<double> const& um, double factor) const;

    double resolution(mesh_type_root const& mesh) const;

    std::vector<double> hminVertices(mesh_type_root const& mesh, BamgMesh const* bamg_mesh) const;
    std::vector<double> hmaxVertices(mesh_type_root const& mesh, BamgMesh const* bamg_mesh) const;

    void initBamg();
    void initOptAndParam();
    void forcing();
    void forcingAtmosphere();
    void forcingOcean();
    void forcingNesting();
	void initBathymetry();

    void assimilateIce();
    void assimilateSlabOcean();
    void initIce();
    void checkConsistency();
    void initThermodynamics();
    void initSlabOcean();
    void initDrifters();
    void updateDrifterPosition();

    void calcCoriolis();
    //void timeInterpolation(int step);
    void nodesToElements(double const* depth, std::vector<double>& v);

    void PwlInterp2D();
    void importBamg(BamgMesh const* bamg_mesh);
    void createGraph();//(BamgMesh const* bamg_mesh);
    void assignVariables();
    void initVariables();
    void initModelState();
    void DataAssimilation();
    void FETensors();
    void calcCohesion();
    void updateVelocity();
    void updateFreeDriftVelocity();
    void speedScaling(std::vector<double>& speed_scaling);
    void scalingVelocity();
    void update();

    void checkOutputs(bool const& at_init_time);
    void exportResults(bool const& export_mesh,
            bool const& export_fields, bool const& apply_displacement);
    void exportResults(std::string const& name_str, bool const& export_mesh,
            bool const& export_fields, bool const& apply_displacement);
    void exportResults(std::vector<std::string> const& filenames, bool const& export_mesh,
            bool const& export_fields, bool const& apply_displacement);

    bool writingRestart();
    void writeRestart();
    void writeRestart(std::string const& name_string);
    int readRestart(int step);
    int readRestart(std::string step);
    void partitionMeshRestart();
    void collectRootRestart(std::vector<double>& interp_elt_out, std::vector<double>& interp_nd_out,
            std::vector<std::vector<double>*> &data,
            std::vector<int> &num_components);

    void rootMeshProcessing();

    void rootMeshRenumbering();

    void distributedMeshProcessing(bool start = false);

    void interpVertices();

    void bcMarkedNodes();

    void finalise();

public:
    std::string gitRevision();
    std::string system(std::string const& command);
    std::string getEnv(std::string const& envname);
    void writeLogFile();

private:
    void advect(std::vector<double> const& interp_elt_in, std::vector<double>& interp_elt_out);
    void advectRoot(std::vector<double> const& interp_elt_in, std::vector<double>& interp_elt_out);
    void diffuse(std::vector<double>& variable_elt, double diffusivity_parameters, double dx);

    void collectVariables(std::vector<double>& interp_elt_in_local, bool ghosts = false);
    void redistributeVariables(std::vector<double> const& out_elt_values, bool check_conc = false);

    // IO
    void collectVariablesIO(std::vector<double>& interp_elt_in_local, bool ghosts, bool thin_ice);
    void gatherFieldsElementIO(std::vector<double>& interp_in_elements, bool thin_ice);

    std::vector<std::string> getRestartVariableNames();
    void getVariablesIO(
            std::vector<std::vector<double>*> &data,
            std::vector<int> &num_components,
            std::vector<std::string> const &names);
    void redistributeVariablesIO(std::vector<double> const& out_elt_values,
            std::vector<std::vector<double>*> &data,
            std::vector<int> const &num_components);
    void scatterFieldsElementIO(std::vector<double> const& out_elt_values,
            std::vector<std::vector<double>*> &data,
            std::vector<int> const &num_components);

    void scatterElementConnectivity();

private:
    po::variables_map vm;
    mesh_type M_mesh;
    graph_type M_graph;
    graphmpi_type M_graphmpi;
    mesh_type M_mesh_init;
    mesh_type M_mesh_previous;
    solver_ptrtype M_solver;
    matrix_ptrtype M_matrix;
    vector_ptrtype M_vector;
    vector_ptrtype M_solution;
    vector_ptrtype M_exact;

    std::map<int, point_type > M_nodes;
    //std::vector<point_type> M_nodes;
    std::vector<element_type> M_edges;
    std::vector<element_type> M_elements;

    int M_num_nodes;
    int M_num_elements;

    int M_ndof;
    int M_local_ndof;
    int M_local_ndof_ghost;
    int M_local_nelements;
    int M_rank;
    Communicator M_comm;

    int M_nb_var_element;
    int M_nb_var_node;

    int M_prv_local_ndof;
    int M_prv_num_nodes;
    int M_prv_num_elements;
    int M_prv_global_num_nodes;
    int M_prv_global_num_elements;

    int pcpt;
    int niter;
    int mesh_adapt_step;
    bool had_remeshed;
    double minang;
    double M_res_root_mesh;

    std::vector<int> M_boundary_flags;
    std::vector<int> M_dirichlet_flags;
    std::vector<int> M_dirichlet_nodes;
    std::vector<int> M_neumann_flags;
    std::vector<int> M_neumann_nodes;

    boost::mpi::timer chrono, chrono_tot;
    std::map<std::string,std::pair<boost::mpi::timer,double> > timer;

    setup::AtmosphereType M_atmosphere_type;
    setup::OceanType M_ocean_type;
    setup::IceType M_ice_type;
    setup::BathymetryType M_bathymetry_type;
    setup::BasalStressType M_basal_stress_type;
    setup::ThermoType M_thermo_type;
    setup::DynamicsType M_dynamics_type;

    setup::IceCategoryType M_ice_cat_type;
    setup::MeshType M_mesh_type;
    mesh::Partitioner M_partitioner;
    mesh::PartitionSpace M_partition_space;

    std::string M_mesh_filename;
    std::string M_mesh_fileformat;

    int M_flag_fix;
    // local masks
    std::vector<bool> M_mask;
    std::vector<bool> M_mask_dirichlet;

    // global masks on root
    std::vector<bool> M_mask_root;
    std::vector<bool> M_mask_dirichlet_root;

    // interpolation method
    std::vector<int> M_interp_method;

    // diffusivity parameters
    std::vector<double> M_diffusivity_parameters;

    std::vector<double> M_surface;
    std::vector<double> M_sigma;
    std::vector<double> M_UM;
    std::vector<double> M_UT;
    std::vector<double> M_VT;
    std::vector<double> M_VTM;
    std::vector<double> M_VTMM;

    // needed for the drifters
    std::vector<double> M_UT_root;
    std::vector<double> M_conc_root;

    std::vector<double> M_bathy_depth;

    std::vector<double> M_hminVertices;
    std::vector<double> M_hmaxVertices;

    std::vector<double> M_Vair_factor;
    std::vector<double> M_Voce_factor;
    std::vector<double> M_basal_factor;
    std::vector<double> M_water_elements;

    std::vector<double> M_h_thin;
    std::vector<double> M_conc_thin;
    std::vector<double> M_hs_thin;

    std::vector<double> M_ridge_ratio;
    std::vector<double> M_h_ridged_thin_ice;
    std::vector<double> M_h_ridged_thick_ice;

    external_data_vec M_external_data;
    external_data_vec M_external_data_tmp;
    Dataset_vec M_datasets_regrid;
    //std::vector<external_data*> M_external_data;

    std::vector<double> M_fcor;

    std::vector<int> M_sizes_nodes;
    std::vector<int> M_sizes_nodes_with_ghost;
    std::vector<int> M_sizes_elements;
    std::vector<int> M_sizes_elements_with_ghost;
    std::vector<int> M_id_nodes;
    std::vector<int> M_rmap_nodes;

    std::vector<int> M_id_elements;
    std::vector<int> M_rmap_elements;
    //std::vector<double> M_speed_scaling;
    std::vector<double> M_element_connectivity;

    std::vector<double> M_Dunit;
    //std::vector<double> M_Dunit_comp;
    std::vector<double> M_Mass;
    std::vector<double> M_Diag;
    std::vector<std::vector<double>> M_shape_coeff;
    std::vector<std::vector<double>> M_B0T;
    std::vector<std::vector<double>> M_B0T_Dunit_B0T;
    //std::vector<std::vector<double>> M_B0T_Dunit_comp_B0T;
    std::vector<double> M_random_number;
    std::vector<double> M_Cohesion;
    std::vector<double> M_Compressive_strength;
    std::vector<double> M_time_relaxation_damage;

    // =============================================================================
    // variables needed for nesting
    bool M_use_nesting;
    bool M_use_ocean_nesting;
    std::string M_nest_outer_mesh;
    std::string M_nest_inner_mesh;
    std::string M_nest_method;
    std::string M_nudge_function;
    double M_nudge_timescale;
    double M_nudge_lengthscale;
    bool M_nest_dynamic_vars;
    // =============================================================================

    LogLevel M_log_level;

private:

    double nu0;
    double young;
    double rhoi;
    double rhos;
    double days_in_sec;
    double time_init;
    double output_time_step;
    double ptime_step;
    double mooring_output_time_step;
    double mooring_time_factor;
    double drifter_output_time_step;
    double restart_time_step;
    double time_step;
    double thermo_timestep;
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

    double basal_k2;
    double basal_drag_coef_air;
    double basal_u_0;
    double basal_Cb;

    double h_thin_min;
    double h_thin_max;

    double compr_strength;
    double tract_coef;
    double scale_coef;
    double alea_factor;
    double cfix;
    double C_fix;
    double C_alea;
    double tan_phi;
    double ridge_h;
    double M_current_time;
    bool M_reuse_prec;
    bool M_regrid;
    int M_nb_regrid;

    bool M_use_restart;
    bool M_write_restart;
    bool M_use_assimilation;

    std::string M_export_path;

private:

    BamgOpts *bamgopt;
    BamgMesh *bamgmesh;
    BamgGeom *bamggeom;

private: // only on root process (rank 0)

    mesh_type_root M_mesh_root;
    mesh_type_root M_mesh_init_root;
    mesh_type_root M_mesh_previous_root;

    //std::vector<double> M_UM_root;
    std::vector<double> M_surface_root;
    std::vector<int> M_dirichlet_flags_root;
    std::vector<int> M_neumann_flags_root;

    std::vector<int> M_dirichlet_nodes_root;
    std::vector<int> M_neumann_nodes_root;

    std::vector<double> M_random_number_root;

    std::vector<std::vector<double>> M_B0T_root;

    BamgMesh *bamgmesh_root;
    BamgGeom *bamggeom_root;

    BamgOpts *bamgopt_previous;
    BamgMesh *bamgmesh_previous;
    BamgGeom *bamggeom_previous;



private:

    // Thermodynamic and dynamic forcing
    // Atmosphere
    external_data M_wind;         // Surface wind [m/s]
    external_data M_tair;         // 2 m temperature [C]
    external_data M_mixrat;       // Mixing ratio
    external_data M_mslp;         // Atmospheric pressure [Pa]
    external_data M_Qsw_in;       // Incoming short-wave radiation [W/m2]
    external_data M_Qlw_in;       // Incoming long-wave radiation [W/m2]
    external_data M_tcc;          // Incoming long-wave radiation [W/m2]
    external_data M_precip;       // Total precipitation [m]
    external_data M_snowfall;     // Snowfall rate [kg/m^2/s]
    external_data M_snowfr;       // Fraction of precipitation that is snow
    external_data M_dair;         // 2 m dew point [C]
    external_data M_sphuma;       // Speciffic humidity of the atmosphere [kg/kg]

    // Ocean
    external_data M_ocean;        // "Geostrophic" ocean currents [m/s]
    external_data M_ssh;          // Sea surface elevation [m]

    external_data M_ocean_temp;   // Ocean temperature in top layer [C]
    external_data M_ocean_salt;   // Ocean salinity in top layer [C]
    external_data M_mld;          // Mixed-layer depth [m]

    // Nesting
    external_data M_nesting_dist_elements; // Distance to the nearest open boundaries
    external_data M_nesting_dist_nodes; // Distance to the nearest open boundaries
    external_data M_ice_conc; // sea_ice_area_fraction from the outer domain
    external_data M_ice_thick; // sea_ice_thickness from the outer domain
    external_data M_ice_snow_thick; // surface_snow_thickness from the outer domain
    external_data M_ice_h_thin ; // thin_ice_thickness from the outer domain
    external_data M_ice_conc_thin ; // thin_ice_area_fraction from the outer domain
    external_data M_ice_hs_thin ; // surface_snow_thickness_on_thin_ice from the outer domain
    external_data M_nesting_damage; // damage from the outer domain
    external_data M_nesting_ridge_ratio; // ridge_ratio from the outer domain
    external_data M_nesting_VT1; // X-velocity from the outer domain
    external_data M_nesting_VT2; // Y-velocity from the outer domain
    external_data M_nesting_sigma1; // 1st component stress tensor from the outer domain
    external_data M_nesting_sigma2; // 2nd component stress tensor from the outer domain
    external_data M_nesting_sigma3; // 3rd component stress tensor from the outer domain

    // Bathymetry
    external_data M_element_depth;

    // Drifters
    bool M_use_drifters;
    bool M_use_iabp_drifters;
    boost::unordered_map<int, std::array<double,2>> M_iabp_drifters; // Drifters are kept in an unordered map containing number and coordinates
    std::fstream M_iabp_file;   // The file we read the IABP buoy data from
    std::fstream M_iabp_out;    // The file we write our simulated drifter positions into

    // Drifters on a grid
    double M_equallyspaced_drifters_output_time_step;
    bool M_use_equallyspaced_drifters;
    Drifters M_equallyspaced_drifters;

    // Drifters as in the RGPS data
    double M_rgps_drifters_output_time_step;
    bool M_use_rgps_drifters;
    Drifters M_rgps_drifters;

    // drifters for the OSISAF emulation
    bool M_use_osisaf_drifters;
    std::vector<Drifters> M_osisaf_drifters;

    // Element variable
    std::vector<double> M_element_age;         // Age of the element (model time since its last adaptation)

    // Prognostic ice variables
    std::vector<double> M_conc;         // Ice concentration
    std::vector<double> M_thick;        // Effective ice thickness [m]
    std::vector<double> M_damage;       // Ice damage
    std::vector<double> M_snow_thick;   // Effective snow thickness [m]

	// Prognostic slab ocean variables
    std::vector<double> M_sst;          // Sea-surface temperature [C]
    std::vector<double> M_sss;          // Sea-surface salinity [psu]

    // Non-prognostic variables used to speed up the convergence of a non-linear equation in thermodynamics
    std::vector<std::vector<double>> M_tice;    // Ice temperature - 0 for surface and higher ordinals for layers in the ice
    std::vector<double> M_tsurf_thin;   // Ice surface temperature of thin ice [C]


private:
    // Variables for the moorings

    std::vector<double> M_conc_mean;    // Mean concentration (on the mesh)
    std::vector<double> M_thick_mean;   // Mean ice thickness (on the mesh)
    std::vector<double> M_snow_thick_mean;  // Mean snow thickness (on the mesh)
    std::vector<double> M_VT_mean;      // Mean velocity (on the mesh)

    std::vector<double> M_conc_grid;    // Mean concentration (on the grid)
    std::vector<double> M_thick_grid;   // Mean ice thickness (on the grid)
    std::vector<double> M_snow_thick_grid;  // Mean snow thickness (on the grid)
    std::vector<double> M_VT_grid;      // Mean velocity (on the grid)

private:
    // Variables for the moorings

    bool M_use_moorings;
    bool M_moorings_snapshot;
    bool M_moorings_parallel_output;
    std::string M_moorings_file;
    GridOutput::fileLength M_moorings_file_length;
    GridOutput M_moorings;

private:

    void constantIce();
    void targetIce();
    void binaryIce();
    void topazIce();
    void topazIceOsisafIcesat();
    void piomasIce();
    void topazForecastIce();
    void topazForecastAmsr2Ice();
    void topazForecastAmsr2OsisafIce();
    void topazForecastAmsr2OsisafNicIce(bool use_weekly_nic);
    void assimilate_topazForecastAmsr2OsisafIce();
    void assimilate_topazForecastAmsr2OsisafNicIce(bool use_weekly_nic);
    void concBinsNic(double &thin_conc_obs_min,double &thin_conc_obs_max,double ci,bool use_weekly_nic);
    void cs2SmosIce();
    void cs2SmosAmsr2Ice();
    void warrenClimatology();
    void smosIce();

    void topazAmsreIce();
    void topazAmsr2Ice();

    void initOSISAFDrifters();
    void initRGPSDrifters();
    void updateRGPSDrifters();
    void equallySpacedDrifter();
    void outputDrifter(std::fstream &iabp_out);
    void initIABPDrifter();
    void updateIABPDrifter();

    //void updateMeans(GridOutput &means);
    void updateMeans(GridOutput& means, double time_factor);
    void initMoorings();
    void updateMoorings();
    void mooringsAppendNetcdf(double const &output_time);
    void checkFields();

private:

    // Diagnostic variables
    std::vector<double> D_Qa; // Heat loss to atmosphere [W/m2]
    std::vector<double> D_Qsw; // Total short wave at surface [W/m2]
    std::vector<double> D_Qlw; // Total long wave at surface [W/m2]
    std::vector<double> D_Qsh; // Total sensible heat flux at surface [W/m2]
    std::vector<double> D_Qlh; // Total latent heat flux at surface [W/m2]
    std::vector<double> D_Qo; // Heat loss from ocean [W/m2]
    std::vector<double> D_delS; // Salt flux to ocean

};
} // Nextsim
#endif
