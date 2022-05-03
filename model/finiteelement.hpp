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

#include "version.hpp"
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/version.hpp>
#include <boost/format.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <BamgConvertMeshx.h>
#include <BamgTriangulatex.h>
#include <Bamgx.h>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromGridToMeshx.h>
#include <gmshmesh.hpp>
#include <gmshmeshseq.hpp>
#include <graphcsr.hpp>
#include <graphcsrmpi.hpp>
#include <externaldata.hpp>
#include <gridoutput.hpp>
#include <dataset.hpp>
#include <model_variable.hpp>
#include <drifters.hpp>
#include "enums.hpp"
#include <debug.hpp>
#include <omp.h>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_01.hpp>
#if defined OASIS
#include<oasis_cpp_interface.h>
#endif
#ifdef AEROBULK
#include "aerobulk.hpp"
#endif
#include "timer.hpp"

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

    FiniteElement(Communicator const& comm = Environment::comm());

    // FiniteElement(Communicator const& comm = Environment::comm());

    mesh_type const& mesh() const {return M_mesh;}

    void initMesh();
    void initExternalData();
    void initDatasets();
    void createGMSHMesh(std::string const& geofilename);

    double jacobian(std::vector<std::vector<double>> const& vertices) const;

    template<typename FEMeshType>
    double jacobian(element_type const& element, FEMeshType const& mesh) const
    { return this->jacobian(mesh.vertices(element.indices)); }

    template<typename FEMeshType>
    double jacobian(element_type const& element, FEMeshType const& mesh,
                    std::vector<double> const& um, double factor = 1.) const
    { return this->jacobian(mesh.vertices(element.indices, um, factor)); }

    std::vector<double> sides(element_type const& element, mesh_type const& mesh) const;
    std::vector<double> sides(element_type const& element, mesh_type const& mesh,
                              std::vector<double> const& um, double factor = 1.) const;

    std::vector<double> sides(element_type const& element, mesh_type_root const& mesh) const;
    std::vector<double> sides(element_type const& element, mesh_type_root const& mesh,
                              std::vector<double> const& um, double factor = 1.) const;

    std::vector<double> minMaxSide(mesh_type_root const& mesh) const;

    template<typename FEMeshType>
    double measure(element_type const& element, FEMeshType const& mesh) const;

    template<typename FEMeshType>
    double measure(element_type const& element, FEMeshType const& mesh,
                   std::vector<double> const& um, double factor = 1.) const;

    std::vector<double> shapeCoeff(element_type const& element) const;
    template<typename FEMeshType>
    std::vector<double> surface(FEMeshType const& mesh);
    template<typename FEMeshType>
    std::vector<double> surface(FEMeshType const& mesh,
            std::vector<double> const& um, double const& factor=1);

    bool checkRegridding();
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

    void init();
    void step();
    void run();

    inline void updateSigmaEVP(double const dte, double const e, double const Pstar, double const C, double const delta_min);
    inline void updateSigmaMEVP(double const dte, double const e, double const Pstar, double const C, double const delta_min, double const alpha);
    void explicitSolve();

    void nestingIce();
    void nestingDynamics();
    void thermo(int dt);
    inline void thermoIce0(const double dt, const double conc, const double voli, const double vols, const double mld, const double snowfall,
            const double Qia, const double dQiadT, const double subl, const double Tbot,
            double &Qio, double &hi, double &hs, double &hi_old, double &del_hi, double &del_hs_mlt, double &mlt_hi_top, double &mlt_hi_bot, double &del_hi_s2i, double &Tsurf);
    inline void thermoWinton(const double dt, const double I_0, const double conc, const double voli, const double vols, const double mld, const double snowfall,
            double const Qia, double const dQiadT, const double Qsw, const double subl, const double Tbot,
            double &Qio, double &hi, double &hs, double &hi_old, double &del_hi, double &del_hs_mlt, double &mlt_hi_top, double &mlt_hi_bot, double &del_hi_s2i,
            double &Tsurf, double &T1, double &T2);
    void OWBulkFluxes(std::vector<double>& Qow, std::vector<double>& Qlw, std::vector<double>& Qsw,
                 std::vector<double>& Qlh, std::vector<double>& Qsh, std::vector<double>& evap, ModelVariable& tau);
    void IABulkFluxes(
            const std::vector<double>& Tsurf, const std::vector<double>& snow_thick,
            const std::vector<double>& conc, std::vector<double>& Qia,
            std::vector<double>& Qlw, std::vector<double>& Qsw,
            std::vector<double>& Qlh, std::vector<double>& Qsh,
            std::vector<double>& subl, std::vector<double>& dQiadT,
            std::vector<double>& alb_tot);
    inline double albedo(const double Tsurf, const double hs,
        int alb_scheme, double alb_ice, double alb_sn, double I_0);
    inline std::pair<double,double> specificHumidity(schemes::specificHumidity scheme, const int i, double temp = -999.);
    inline double iceOceanHeatflux(const int cpt, const double sst, const double tbot, const double mld, const double dt);
    inline double incomingLongwave(const int i);
    inline double freezingPoint(const double sss);
    inline double windSpeedElement(const int i);

    void checkReloadDatasets(external_data_vec const& ext_data_vec,
                    double const CRtime, std::vector<double> &RX, std::vector<double> &RY);
    void checkReloadMainDatasets(double const CRtime);

    Dataset M_atmosphere_nodes_dataset;
    Dataset M_atmosphere_elements_dataset;
    Dataset M_atmosphere_bis_elements_dataset;
    Dataset M_ocean_nodes_dataset;
    Dataset M_ocean_elements_dataset;
    Dataset M_bathymetry_elements_dataset;
#ifdef OASIS
    Dataset M_wave_nodes_dataset;
    Dataset M_wave_elements_dataset;
#endif

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
    void initFETensors();
    template<typename enum_type>
    enum_type getOptionFromMap(std::string const &opt_name,
        boost::unordered_map<const std::string, enum_type> map);
    void forcing();
    void forcingAtmosphere();
    void forcingOcean();
#ifdef OASIS
    void forcingWaves();
#endif
    void forcingNesting();
    void initBathymetry();

    void assimilateIce();
    void assimilateSlabOcean();
    void initIce();
    void checkConsistency();
    void initSlabOcean();

    void calcCoriolis();
    //void timeInterpolation(int step);
    void nodesToElements(double const* depth, std::vector<double>& v);

    void PwlInterp2D();
    void importBamg(BamgMesh const* bamg_mesh);
    void createGraph();//(BamgMesh const* bamg_mesh);
    void assignVariables();
    void initVariables();
    void calcAuxiliaryVariables();
    void initModelVariables();
    void sortPrognosticVars();
    void initModelState();
    void DataAssimilation();

    void calcCohesion();
    void updateFreeDriftVelocity();
    void speedScaling(std::vector<double>& speed_scaling);
    void update(std::vector<double> const & UM_P);
    void updateSigmaDamage(double const dt);

    void updateGhosts(std::vector<double>& mesh_nodal_vec);
    void initUpdateGhosts();
    int globalNumToprocId(int global_num);

#ifdef OASIS
    bool M_couple_waves;
    bool M_recv_wave_stress;
    // FSD related functions
    void initFsd();
    void redistributeFSD();
    void updateFSD();
    std::vector<double> computeWaveBreakingProb();
    double computeLateralAreaFSD(const int cpt);
    double computeLeadFractionFSD(const int cpt);
    void weldingRoach(const int cpt, double ddt);
    void redistributeThermoFSD(const int i,double ddt, double lat_melt_rate, double young_ice_growth, double old_conc, double old_conc_young) ;
    double lateralMeltFSD(const int i,double ddt) ;
#endif

    void checkOutputs(bool const& at_init_time);
    void exportResults(bool const& export_mesh,
            bool const& export_fields, bool const& apply_displacement);
    void exportResults(std::string const& name_str, bool const& export_mesh,
            bool const& export_fields, bool const& apply_displacement);
    void exportResults(std::vector<std::string> const& filenames, bool const& export_mesh,
            bool const& export_fields, bool const& apply_displacement);
    void updateIceDiagnostics();

    void writeRestart();
    void writeRestart(std::string const& name_string);
    void readRestart(std::string const& name_string);
    void partitionMeshRestart();
    void collectNodesRestart(std::vector<double>& interp_nd_out);
    void collectElementsRestart(std::vector<double>& interp_elt_out,
            std::vector<std::vector<double>*> &data_elements_root);

    void rootMeshProcessing();

    void rootMeshRenumbering();

    void distributedMeshProcessing(bool start = false);

    void interpVertices();

    void bcMarkedNodes();

    void finalise(std::string current_time_system);

public:
    std::string system(std::string const& command);
    std::string getEnv(std::string const& envname);
    void writeLogFile();

private:
    void advect(std::vector<double> const& interp_elt_in, std::vector<double>& interp_elt_out);
    void advectRoot(std::vector<double> const& interp_elt_in, std::vector<double>& interp_elt_out);
    void diffuse(std::vector<double>& variable_elt, double diffusivity_parameters, double dx);

    void collectVariables(std::vector<double>& interp_elt_in_local, bool ghosts);
    void redistributeVariables(std::vector<double> const& out_elt_values, bool const& apply_maxima);

    // IO
    void collectVariablesIO(std::vector<double>& elt_values_local,
            std::vector<ModelVariable*> const& vars_elements,
            std::vector<ExternalData*> const& ext_data_elements,
            bool const& ghosts);
    void gatherFieldsElementIO(std::vector<double>& elt_values_root,
            std::vector<ModelVariable*> const& vars_elements,
            std::vector<ExternalData*> const& ext_data_elements);
    void gatherFieldsElementIO(std::vector<double>& elt_values_root,
            std::vector<ModelVariable*> const& vars_elements)
    {
        std::vector<ExternalData*> ext_data_elements = {};// add a place-holder
        this->gatherFieldsElementIO(elt_values_root, vars_elements, ext_data_elements);
    }

    void redistributeVariablesIO(std::vector<double> const& out_elt_values,
            std::vector<ModelVariable*> &vars_elements);
    void scatterFieldsElementIO(std::vector<double> const& interp_elt_out,
        std::vector<ModelVariable*> &vars_elements);

    void scatterElementConnectivity();

private:
    po::variables_map vm;
    mesh_type M_mesh;
    graph_type M_graph;
    graphmpi_type M_graphmpi;
    mesh_type M_mesh_init;
    mesh_type M_mesh_previous;

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
    Timer M_timer;

    setup::AtmosphereType M_atmosphere_type;
    setup::OceanType M_ocean_type;
    setup::IceType M_ice_type;
    setup::BathymetryType M_bathymetry_type;
    setup::BasalStressType M_basal_stress_type;
    setup::ThermoType M_thermo_type;
    setup::DynamicsType M_dynamics_type;
    int M_ensemble_member;

#ifdef AEROBULK
    aerobulk::algorithm M_ocean_bulk_formula;
#endif

    setup::FreezingPointType M_freezingpoint_type;
    setup::OceanHeatfluxScheme M_Qio_type;
    setup::IceCategoryType M_ice_cat_type;
    setup::MeshType M_mesh_type;
    mesh::Partitioner M_partitioner;
    mesh::PartitionSpace M_partition_space;
    //fsd related
    setup::WeldingType M_welding_type    ;
    setup::BreakupType M_breakup_type    ;
    setup::FSDType M_fsd_type    ;

    bool M_flooding;

    std::string M_mesh_basename;
    std::string M_mesh_filename;
    std::string M_partitioned_mesh_filename;
    std::string M_mesh_fileformat;

    int M_flag_fix;
    // local masks
    std::vector<bool> M_mask;
    std::vector<bool> M_mask_dirichlet;

    // global masks on root
    std::vector<bool> M_mask_root;
    std::vector<bool> M_mask_dirichlet_root;

    // diffusivity parameters
    std::vector<double> M_diffusivity_parameters;

    std::vector<double> M_surface;
    std::vector<double> M_delta_x;

    std::vector<double> M_UM;
    std::vector<double> M_UT;
    std::vector<double> M_VT;


    std::vector<double> M_hminVertices;
    std::vector<double> M_hmaxVertices;

    std::vector<double> M_Vair_factor;
    std::vector<double> M_Voce_factor;
    std::vector<double> M_basal_factor;
    std::vector<double> M_water_elements;


#ifdef OASIS
    ExternalData M_tau_wi;
    ExternalData M_wlbk;
//    ExternalData M_str_var;
//    ExternalData M_tm02;
#endif

    external_data_vec M_external_data_elements;
    std::vector<std::string> M_external_data_elements_names;//list of names for debugging and exporting
    external_data_vec M_external_data_nodes;
    std::vector<std::string> M_external_data_nodes_names;//list of names for debugging and exporting
    Dataset_vec M_datasets_regrid;

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
    std::vector<std::vector<double>> M_shape_coeff;
    std::vector<std::vector<double>> M_B0T;
    std::vector<double> M_Cohesion;
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
    bool M_log_all;


private:

    double nu0;
    double young;
    double rhoi;
    double rhos;
    double const days_in_sec  = 86400.;
    double time_init;
    int output_time_step;
    int ptime_step;
    int mooring_output_time_step;
    double mooring_time_factor;
    int restart_time_step;
    int time_step;
    double dtime_step;
    double duration;
    double divergence_min;
    double compression_factor;
    double exponent_compression_factor;
    double exponent_cohesion;
    double ocean_turning_angle_rad;
    double compaction_param;
    double undamaged_time_relaxation_sigma;
    double exponent_relaxation_sigma;
    double quad_drag_coef_air;
    double quad_drag_coef_water;
    double lin_drag_coef_air;
    double lin_drag_coef_water;
    double time_relaxation_damage;
    double deltaT_relaxation_damage;

    double basal_k2;
    double basal_drag_coef_air;
    double basal_u_0;
    double basal_Cb;

    double h_young_min;
    double h_young_max;
    double M_ks;
    double M_ocean_albedo;

    double compr_strength;
    double tract_coef;
    double scale_coef;
    double alea_factor;
    double C_lab;
    double C_fix;
    double C_alea;
    double tan_phi;
    double ridge_h;
    double M_current_time;
    bool M_reuse_prec;
    bool M_regrid;
    int M_nb_regrid;

    bool M_use_assimilation;

    bool M_use_restart;
    bool M_check_restart;
    bool M_write_restart_interval;
    bool M_write_restart_end;
    bool M_write_restart_start;

    double M_spinup_duration;

    std::string M_export_path;

private: // update solution from explicit solver
    std::vector<std::vector<int>> M_extract_local_index;
    std::vector<int> M_recipients_proc_id;
    std::vector<int> M_local_ghosts_proc_id;
    std::vector<std::vector<int>> M_local_ghosts_local_index;

private:

    BamgOpts *bamgopt;
    BamgMesh *bamgmesh;
    BamgGeom *bamggeom;

private: // only on root process (rank 0)

    mesh_type_root M_mesh_root;
    mesh_type_root M_mesh_init_root;
    mesh_type_root M_mesh_previous_root;

    std::vector<int> M_connectivity_root;
    std::vector<int> M_dirichlet_flags_root;
    std::vector<int> M_neumann_flags_root;

    std::vector<int> M_dirichlet_nodes_root;
    std::vector<int> M_neumann_nodes_root;

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

    external_data M_qsrml;        // Fraction of short wave radiation absorbed by the mixed layer

    // Nesting
    external_data M_nesting_dist_elements; // Distance to the nearest open boundaries
    external_data M_nesting_dist_nodes; // Distance to the nearest open boundaries
    external_data M_nesting_conc; // sea_ice_area_fraction from the outer domain
    external_data M_nesting_thick; // sea_ice_thickness from the outer domain
    external_data M_nesting_snow_thick; // surface_snow_thickness from the outer domain
    external_data M_nesting_h_young ; // young_ice_thickness from the outer domain
    external_data M_nesting_conc_young ; // young_ice_area_fraction from the outer domain
    external_data M_nesting_hs_young ; // surface_snow_thickness_on_young_ice from the outer domain
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
    std::vector<Drifters> M_drifters;// vector of all the Drifters objects (including IABP ones)
    std::vector<int> M_osisaf_drifters_indices;// indices of OSISAF drifters in M_drifters

    // Element variable
    std::vector<double> M_element_age;         // Age of the element (model time since its last adaptation)

    // vectors of pointers to variables (for looping)
    std::vector<ModelVariable*> M_variables_elt;
    std::vector<ModelVariable*> M_prognostic_variables_elt;//for restart, regrid
    std::vector<ModelVariable*> M_export_variables_elt;

    // other vectors related to export/restart
    std::vector<std::string> M_restart_names_elt;
    std::vector<std::string> M_export_names_elt;

    // Prognostic variables
    ModelVariable M_conc;               // Ice concentration
    ModelVariable M_thick;              // Effective ice thickness [m]
    ModelVariable M_damage;             // Ice damage
#ifdef OASIS
    ModelVariable M_cum_damage;         // Ice cumulated damage
    ModelVariable M_cum_wave_damage;    // Ice cumulated damage due to wave
#endif
    ModelVariable M_snow_thick;         // Effective snow thickness [m]
    ModelVariable M_ridge_ratio;
    std::vector<ModelVariable> M_tice;  // Ice temperature - 0 for surface and higher ordinals for layers in the ice
    std::vector<ModelVariable> M_sigma; // Internal stress tensor
    ModelVariable M_sst;                // Sea-surface (slab ocean) temperature [C]
    ModelVariable M_sss;                // Sea-surface (slab ocean) salinity [psu]
    ModelVariable M_tsurf_young;        // Ice surface temperature of young ice [C]
    ModelVariable M_h_young;
    ModelVariable M_hs_young;
    ModelVariable M_conc_young;
    ModelVariable M_random_number;
    ModelVariable M_fyi_fraction;
    ModelVariable M_age_det;
    ModelVariable M_age;
    ModelVariable M_conc_upd;           // Ice concentration update by assimilation

#ifdef OASIS
    // Following variables are related to floe size distribution
    std::vector<ModelVariable> M_conc_fsd;
    //std::vector<ModelVariable> M_conc_fsd_thick;
    //std::vector<ModelVariable> M_conc_fsd_young ;
    std::vector<ModelVariable> M_conc_mech_fsd;
    int M_num_fsd_bins;
    std::vector<double> M_fsd_bin_widths;
    double M_fsd_bin_cst_width;
    double M_fsd_min_floe_size;
    std::vector<double> M_fsd_bin_centres;
    std::vector<double> M_fsd_bin_low_limits;
    std::vector<double> M_fsd_bin_up_limits;
    double M_fsd_unbroken_floe_size     ;
    // Non-circularity of floes
    double M_floe_shape                              ;
    // Lettie's variables
    std::vector<double> M_floe_area_up             ;
    std::vector<double> M_floe_area_low            ;
    std::vector<double> M_floe_area_centered       ;
    std::vector<double> M_floe_area_binwidth       ;

    std::vector<double> M_fsd_area_scaled_up       ;
    std::vector<double> M_fsd_area_scaled_low      ;
    std::vector<double> M_fsd_area_scaled_centered ;
    std::vector<double> M_fsd_area_scaled_binwidth ;
    std::vector<double> M_fsd_area_lims            ;
    std::vector<double> M_fsd_area_lims_scaled     ;

    std::vector<std::vector<int> > M_alpha_fsd_merge ;
    // In namelist
    bool   M_distinguish_mech_fsd             ;
    bool   M_debug_fsd                        ;
    int    M_fsd_damage_type                  ;
    double M_floes_flex_strength              ;
    double M_floes_flex_young                 ;
    double M_welding_kappa                    ;
    bool   M_fsd_welding_use_scaled_area      ;
    double M_dmax_c_threshold                 ;
    double M_breakup_thick_min                ;
    bool   M_breakup_in_dt                    ;
    bool   M_breakup_cell_average_thickness   ;
    // Horvat et Tziperman (2015) lead fraction, lat. surf and lead width
    // double M_lead_width    ;
#endif

    // Diagnostic variables
    ModelVariable D_conc; //total concentration
    ModelVariable D_thick; // total thickness [m]
    ModelVariable D_snow_thick;// total snow thickness [m]
    ModelVariable D_tsurf; //mean surface temperature (thick + young ice + slab ocean) [deg C]
    std::vector<ModelVariable> D_sigma; //principal stresses [Pa]
    ModelVariable D_divergence; // ice velocity divergence [1/s]
    ModelVariable D_Qa; // Heat loss to atmosphere [W/m2]
    ModelVariable D_Qsw; // Total short wave at surface [W/m2]
    ModelVariable D_Qlw; // Total long wave at surface [W/m2]
    ModelVariable D_Qsh; // Total sensible heat flux at surface [W/m2]
    ModelVariable D_Qlh; // Total latent heat flux at surface [W/m2]
    ModelVariable D_Qo; // Heat loss from ocean [W/m2]
    ModelVariable D_Qnosun; // Non-solar heat loss from ocean [W/m2]
    ModelVariable D_Qsw_ocean; // SW flux out of the ocean [W/m2]
    ModelVariable D_Qassim; // flux from assim [W/m2]
    ModelVariable D_vice_melt; // ice volume (/element_area) melted/formed [m/day]
    ModelVariable D_del_vi_young; // young ice volume (/element_area) melted/formed [m/day]
    ModelVariable D_del_hi; // ice growth/melt rate [m/sdat
    ModelVariable D_del_hi_young; // young ice growth/melt rate [m/day]
    ModelVariable D_newice; // ice volume (/element_area) formed in open water [m/day]
    ModelVariable D_mlt_top; // ice volume (/element_area) melted at top [m/day]
    ModelVariable D_mlt_bot; // ice volume (/element_area) melted at bottom [m/day]
    ModelVariable D_snow2ice; // ice volume (/element_area) melted at bottom [m/day]
    ModelVariable D_delS; // Salt flux to ocean
    ModelVariable D_fwflux; // Fresh-water flux at ocean surface [kg/m2/s]
    ModelVariable D_fwflux_ice; // Fresh-water flux at ocean surface due to ice processes [kg/m2/s]
    ModelVariable D_brine; // Brine release into the ocean [kg/m2/s]
    ModelVariable D_dmax; //max floe size [m]
    ModelVariable D_dmean; //mean floe size [m]
    ModelVariable D_tau_ow; // Ocean atmosphere drag coefficient - still needs to be multiplied with the wind [Pa/s/m] (for the coupled ice-ocean system)
    ModelVariable D_evap; // Evaporation out of the ocean [kg/m2/s]
    ModelVariable D_rain; // Rain into the ocean [kg/m2/s]
    ModelVariable D_albedo; // total albedo

    // Temporary variables
    std::vector<double> D_tau_w; // Ice-ocean drag [Pa]
    std::vector<double> D_tau_a; // Ice-atmosphere drag [Pa]

private:
    // Variables for the moorings

    std::vector<double> M_conc_mean;    // Mean concentration (on the mesh)
    std::vector<double> M_thick_mean;   // Mean ice thickness (on the mesh)
    std::vector<double> M_snow_thick_mean;  // Mean snow thickness (on the mesh)
    std::vector<double> M_VT_mean;      // Mean velocity (on the mesh)
    std::vector<double> M_fyi_fraction_mean;  // Fraction of the first year ice (FYI) (on the mesh)
    std::vector<double> M_age_det_mean;       // Ice age observable from space (area weighted) [timestep] (on the mesh)
    std::vector<double> M_age_mean;           // Effective ice age [timestep] (on the mesh)

    std::vector<double> M_conc_grid;    // Mean concentration (on the grid)
    std::vector<double> M_thick_grid;   // Mean ice thickness (on the grid)
    std::vector<double> M_snow_thick_grid;  // Mean snow thickness (on the grid)
    std::vector<double> M_VT_grid;      // Mean velocity (on the grid)
    std::vector<double> M_fyi_fraction_grid;  // Fraction of the first year ice (FYI) (on the grid)
    std::vector<double> M_age_det_grid;       // Ice age observable from space (area weighted) [timestep] (on the grid)
    std::vector<double> M_age_grid;           // Effective ice age [timestep] (on the grid)

private:
    // Variables for the moorings

    bool M_use_moorings;
    bool M_moorings_snapshot;
    bool M_moorings_parallel_output;
    std::string M_moorings_file;
    GridOutput::fileLength M_moorings_file_length;
    GridOutput M_moorings;
    bool M_moorings_false_easting;
    double M_moorings_averaging_period;

#ifdef OASIS
    // Coupling with OASIS
    GridOutput M_cpl_out;
    std::vector<int> var_id_snd;
    std::vector<int> var_id_rcv;

    std::vector<std::string> var_snd;
    std::vector<std::string> var_rcv;

    int cpl_time_step;
    void initOASIS();
    void setCplId_rcv(DataSet &dataset);
    void setCplId_snd(std::vector<GridOutput::Variable> &cpl_var);
#endif

private:

    //ice-init functions
    void constantIce();
    void topazIce();
    void topazIceOsisafIcesat();
    void piomasIce();
    void nemoIce();
    void ciceIce();
    void topazForecastIce();
    void topazForecastAmsr2Ice();
    void topazForecastAmsr2OsisafIce();
    void topazForecastAmsr2OsisafNicIce(bool use_weekly_nic);
    void concBinsNic(double &young_conc_obs_min,double &young_conc_obs_max,double ci,bool use_weekly_nic);
    void cs2SmosIce();
    void cs2SmosAmsr2Ice();
    void smosIce();

    //no ice-type option to activate these
    void topazAmsreIce();
    void topazAmsr2Ice();
    void amsr2ConstThickIce();

    void warrenClimatology();
    void assimilate_topazForecastAmsr2OsisafIce();
    void assimilate_topazForecastAmsr2OsisafNicIce(bool use_weekly_nic);

    //drifter functions
    void checkMoveDrifters();
    void checkUpdateDrifters();
    void instantiateDrifters();
    void synchroniseOsisafDrifters();

    //void updateMeans(GridOutput &means);
    void updateMeans(GridOutput& means, double time_factor);
    void initMoorings();
    void updateMoorings();
    void mooringsAppendNetcdf(double const &output_time);
    void checkFields();
    void checkFieldsFast();
    void checkVelocityFields();

};
} // Nextsim
#endif
