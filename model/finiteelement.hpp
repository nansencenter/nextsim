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
#include <Bamgx.h>
#include <InterpFromMeshToMesh2dCavities.h>
#include <InterpFromMeshToMesh2dx.h>
#include <InterpFromGridToMeshx.h>
#include <InterpFromMeshToGridx.h>
#include <gmshmesh.hpp>
#include <graphcsr.hpp>
#if defined (WAVES)
#include <wimdiscr.hpp>
#endif
#if defined OASIS
#include<oasis_cpp_interface.h>
#endif
#include "enums.hpp"
#include <debug.hpp>
#include <omp.h>
#include <externaldata.hpp>
#include <gridoutput.hpp>
#include <dataset.hpp>
#include <drifters.hpp>

#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

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

    typedef DataSet Dataset;
    typedef ExternalData external_data;
    typedef typename std::vector<external_data*> external_data_vec ;
    typedef typename std::vector<Dataset*> Dataset_vec ;

    typedef boost::ptr_vector<external_data> externaldata_ptr_vector;

#if defined (WAVES)
    typedef Wim::WimDiscr<double> wim_type;
#endif

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

    void initMesh(setup::MeshType const& mesh_type);
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
    void nodeMaxConc();
    void solve();
    void run();
    int init();
    void step(int &pcpt);
    void finalise();
    void error();

    void thermo();
    void thermoIce0(int i, double wspeed, double sphuma, double conc, double voli, double vols, double Qlw_in, double Qsw_in, double mld, double snowfall,
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
    Dataset M_ice_cs2_smos_elements_dataset;
    Dataset M_ice_smos_elements_dataset;
#if defined (WAVES)
    Dataset M_wave_elements_dataset;
#endif
    double minAngles(element_type const& element, mesh_type const& mesh) const;
    double minAngle(mesh_type const& mesh) const;

    double minAngle(mesh_type const& mesh, std::vector<double> const& um, double factor) const;

    bool flip(mesh_type const& mesh, std::vector<double> const& um, double factor) const;

    double resolution(mesh_type const& mesh) const;

    std::vector<double> hminVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh) const;
    std::vector<double> hmaxVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh) const;

    std::vector<double> AllMinAngle(mesh_type const& mesh, std::vector<double> const& um, double factor) const;

    void initBamg();
    void initConstant();
    void forcing();
    void forcingAtmosphere();
    void forcingOcean();

#if defined (WAVES)
    void forcingWave();
    WaveOptions wim_forcing_options;
    bool        wim_ideal_forcing;
#endif

	void bathymetry();
    void checkReloadDatasets(external_data_vec const& ext_data_vec,
        double const& CRtime, std::string const& printout);
    void assimilateIce();
    void assimilateSlabOcean();
    void initIce();
    void initThermodynamics();
    void initSlabOcean();
    void initDrifter();
    void coriolis();
    void nodesToElements(double const* depth, std::vector<double>& v);

    void PwlInterp2D();
    void importBamg(BamgMesh const* bamg_mesh);
    void initVariables();
    void initModelState();
    void DataAssimilation();
    void tensors();
    void cohesion();
    void updateVelocity();
    void updateFreeDriftVelocity();
    void update();
    void exportInitMesh();
    void exportResults(int step,
            bool export_mesh = true, bool export_fields = true, bool apply_displacement = true);
    void exportResults(std::string const name_str,
            bool export_mesh = true, bool export_fields = true, bool apply_displacement = true);
    void exportResults(std::vector<std::string> const &filenames,
            bool export_mesh = true, bool export_fields = true, bool apply_displacement = true);

    void writeRestart(int pcpt, int step);
    int readRestart(int step);

#if defined (WAVES)
    void initWim(int const pcpt);
    void initWimVariables();
    void nextsimToWim();
    void wimToNextsim();
#if 0
    std::vector<double> FiniteElements::rotatedWimElementsX(double const& rotangle) const;
    std::vector<double> FiniteElements::rotatedWimElementsY(double const& rotangle) const;
#endif
#endif

    std::string gitRevision();
    std::string system(std::string const& command);
    std::string getEnv(std::string const& envname);
    void writeLogFile();

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
    setup::WaveType M_wave_type;
    setup::BathymetryType M_bathymetry_type;
    setup::BasalStressType M_basal_stress_type;
    setup::ThermoType M_thermo_type;
    setup::DynamicsType M_dynamics_type;

    setup::IceCategoryType M_ice_cat_type;
    setup::MeshType M_mesh_type;

    LogLevel M_log_level;

    std::string M_mesh_filename;

    int M_flag_fix;
    std::vector<bool> M_mask ;
    std::vector<bool> M_mask_dirichlet ;

    int mesh_adapt_step;
    bool had_remeshed;

    std::vector<double> M_surface;
    std::vector<double> M_sigma;
    std::vector<double> M_UM;
    std::vector<double> M_UT;
    std::vector<double> M_VT;
    std::vector<double> M_VTM;
    std::vector<double> M_VTMM;

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

    external_data_vec M_external_data;
    external_data_vec M_external_data_tmp;
    Dataset_vec M_datasets_regrid;

    std::vector<double> M_fcor;

    std::vector<double> M_Dunit;
    std::vector<double> M_Mass;
    std::vector<double> M_Diag;
    std::vector<std::vector<double>> M_shape_coeff;
    std::vector<std::vector<double>> M_B0T;
    std::vector<std::vector<double>> M_B0T_Dunit_B0T;
    std::vector<double> M_random_number;
    std::vector<double> M_Cohesion;
    std::vector<double> M_Compressive_strength;
    std::vector<double> M_time_relaxation_damage;

    // =============================================================================
    // variables needed for coupling with wim
#if defined (WAVES)
    wim_type wim;
    std::vector<double> M_nfloes;
    std::vector<double> M_dfloe;

    std::vector<double> M_SWH_grid;
    std::vector<double> M_MWD_grid;
    std::vector<double> M_MWP_grid;

    wim_type::WimGrid wim_grid;
    double xmin_wim,xmax_wim;
    double ymin_wim,ymax_wim;
    int num_elements_wim_grid;
    int wim_cpt;//no of times WIM has been called
    int steps_since_last_wim_call;//no of time steps since WIM was last called

    std::vector<double> M_stokes_drift;

    std::vector<double> M_icec_grid;
    std::vector<double> M_iceh_grid;
    std::vector<double> M_nfloes_grid;

    std::vector<double> M_taux_grid;
    std::vector<double> M_tauy_grid;
    std::vector<double> M_stokes_drift_x_grid;
    std::vector<double> M_stokes_drift_y_grid;
    bool M_export_stokes_drift_mesh;
#endif
    std::vector<double> M_tau;//this can just be set to zero if not using WIM
    // =============================================================================

private:

    double nu0;
    double young;
    double rhoi;
    double rhos;
    double days_in_sec;
    double time_init;
    double output_time_step;
    double mooring_output_time_step;
    double mooring_time_factor;
    double restart_time_step;
    double time_step;
    double duration;
    double ocean_turning_angle_rad;
    double ridging_exponent;
    double quad_drag_coef_air;
    double quad_drag_coef_water;
    double time_relaxation_damage;
    double deltaT_relaxation_damage;
    //double water_depth;

    double basal_k2;
    double basal_u_0;
    double basal_Cb;

    double h_thin_max;
    double h_thin_min;

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
    int M_nb_regrid;
#if defined (WAVES)
    bool M_run_wim;
    bool M_use_wim;
    bool M_interp_fsd;
#endif

    bool M_use_restart;
    bool M_write_restart;
    bool M_use_assimilation;

    std::string M_export_path;

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
    external_data M_tcc;          // Incoming long-wave radiation [W/m2]
    external_data M_precip;       // Total precipitation rate [kg/m^2/s]
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

    // Wave
    external_data M_SWH;	      // Significant wave height [m]
    external_data M_MWD;	      // Mean wave direction (deg)
    external_data M_MWP;          // Peak wave period (s)
    external_data M_fice_waves;   // Waves masked if ice used in external wave model 
                                  // - due to inconsistent ice masks,
                                  // there could be attenuation in the open ocean

    // Bathymetry
    external_data M_element_depth;

    // IABP-like drifters
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
    // std::vector<double> M_tsurf;        // Ice surface temperature [C]
    std::vector<std::vector<double>> M_tice;    // Ice temperature - 0 for surface and higher ordinals for layers in the ice
    std::vector<double> M_tsurf_thin;   // Ice surface temperature of thin ice [C]

private:
    // Variables for the moorings

    bool M_use_moorings;
    bool M_moorings_snapshot;
    std::string M_moorings_file;
    GridOutput::fileLength M_moorings_file_length;
    GridOutput M_moorings;

#ifdef OASIS
    // Coupling with OASIS
    GridOutput M_cpl_out;
    std::vector<int> var_id;
    int cpl_time_step;
    void initOASIS();
#endif

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
    void topazForecastAmsr2OsisafNicIce();
    void assimilate_topazForecastAmsr2OsisafIce();
    void assimilate_topazForecastAmsr2OsisafNicIce();
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

    void updateMeans(GridOutput &means, double time_factor);
    void initMoorings();

    void redistributeVariables(double* interp_elt_out,int nb_var, bool check_conc_max);
    int collectVariables(double** interp_elt_in_ptr, int** interp_elt_method, double** diffusivity_parameters, int num_elements);
    void advect (double** interp_elt_out_ptr,double* interp_elt_in, int* interp_method, int nb_var);
    void diffuse(double* variable_elt, double diffusivity_parameters, double dx);

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
