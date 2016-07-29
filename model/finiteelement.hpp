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
#include <gmshmeshseq.hpp>
#include <graphcsr.hpp>
#include <graphcsrmpi.hpp>
#include "enums.hpp"
#include <netcdf>
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

    typedef struct Dimension
    {
        std::string name;
        int start;
        int end;
    } Dimension;

    typedef struct Variable
    {
        std::string name;
        std::vector<Dimension> dimensions;
        double a;
        double b;
        std::string Units;
        std::vector<std::vector<double>> data2;
    } Variable;

    typedef struct Grid
    {
        setup::InterpolationType interpolation_method;
		int interp_type;
        std::string dirname;
        std::string filename;

        Variable latitude;
        Variable longitude;

        Dimension dimension_x;
        Dimension dimension_y;

        std::string mpp_file;
        double rotation_angle;
		bool interpolation_in_latlon;

		bool masking;
		Variable masking_variable;
		std::vector<int> reduced_nodes_ind;

        int* pfindex;
        int pfnels;
        std::vector<double> gridX;
        std::vector<double> gridY;

        std::vector<double> gridLAT;
        std::vector<double> gridLON;
    } Grid;

    typedef struct Dataset
    {
        std::string dirname;
        std::string prefix;
        std::string postfix;
        std::string reference_date;

        std::vector<Variable> variables;
        int target_size;
        Grid *grid;

        int nb_timestep_day;
        Variable time;
        Dimension dimension_time;

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

    void initMesh();
    void initDatasets();
    void createGMSHMesh(std::string const& geofilename);

    double jacobian(element_type const& element, mesh_type const& mesh) const;

    double jacobian(element_type const& element, mesh_type_root const& mesh) const;

    double jacobian(element_type const& element, mesh_type const& mesh,
                    std::vector<double> const& um, double factor = 1.) const;

    double jacobian(element_type const& element, mesh_type_root const& mesh,
                    std::vector<double> const& um, double factor = 1.) const;

    std::vector<double> sides(element_type const& element, mesh_type const& mesh) const;

    std::vector<double> sides(element_type const& element, mesh_type_root const& mesh) const;

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

    void gatherSizes();
    void gatherFieldsElement(std::vector<double>& interp_in_elements);
    void scatterFieldsElement(double* interp_elt_out);
    void interpFieldsElement();

    void gatherUM(std::vector<double>& um);
    void interpFieldsNode(bimap_type const& rmap_nodes, std::vector<int> sizes_nodes);

    void assemble(int cpt);
    void solve();
    void run();
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

    template<typename FEMeshType>
    double minAngles(element_type const& element, FEMeshType const& mesh) const;

    template<typename FEMeshType>
    double minAngle(FEMeshType const& mesh) const;

    template<typename FEMeshType>
    double minAngle(FEMeshType const& mesh, std::vector<double> const& um, double factor) const;

    template<typename FEMeshType>
    bool flip(FEMeshType const& mesh, std::vector<double> const& um, double factor) const;

    double resolution(mesh_type_root const& mesh) const;

    std::vector<double> hminVertices(mesh_type_root const& mesh, BamgMesh const* bamg_mesh) const;
    std::vector<double> hmaxVertices(mesh_type_root const& mesh, BamgMesh const* bamg_mesh) const;

    void initBamg();
    void initConstant();
    void forcing();
    void forcingAtmosphere(bool reload);
    void forcingOcean(bool reload);
	void bathymetry(bool reload);

    void initIce();
    void initThermodynamics();
    void initSlabOcean();
    void initDrifter();
    void coriolis();
    void timeInterpolation(int step);
    void nodesToElements(double const* depth, std::vector<double>& v);

    void PwlInterp2D();
    void importBamg(BamgMesh const* bamg_mesh);
    void createGraph();//(BamgMesh const* bamg_mesh);
    void assignVariables();
    void initVariables();
    void initModelState();
    void tensors();
    void tensorsOnRoot();
    void cohesion();
    void updateVelocity();
    void speedScaling(std::vector<double>& speed_scaling);
    void scalingVelocity();
    void update();
    void updateOnRoot();
    void exportResults(int step, bool export_mesh = true);

    void writeRestart(int pcpt, int step);
    void readRestart(int &pcpt, int step);

    void rootMeshProcessing();

    void rootMeshRenumbering();

    void distributedMeshProcessing(bool start = false);

    void interpVertices();

    void bcMarkedNodes();

    void clear();

private:
    po::variables_map vm;
    mesh_type M_mesh;
    graph_type M_graph;
    graphmpi_type M_graphmpi;
    mesh_type M_mesh_init;
    mesh_type M_mesh_previous;
    solver_ptrtype M_solver;
    matrix_ptrtype M_matrix;
    //matrix_ptrtype M_mass;
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

    setup::IceCategoryType M_ice_cat_type;
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

    std::vector<int> M_sizes_nodes;
    std::vector<int> M_sizes_nodes_with_ghost;
    std::vector<int> M_sizes_elements;
    std::vector<int> M_sizes_elements_with_ghost;
    std::vector<int> M_id_nodes;
    //std::vector<double> M_speed_scaling;

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

    LogLevel M_log_level;

private:

    double nu0;
    double young;
    double rhoi;
    double rhos;
    double days_in_sec;
    double time_init;
    double output_time_step;
    double restart_time_step;
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

    BamgOpts *bamgopt;
    BamgMesh *bamgmesh;
    BamgGeom *bamggeom;

    // BamgOpts *bamgopt_previous;
    // BamgMesh *bamgmesh_previous;
    // BamgGeom *bamggeom_previous;

    //Options *options;

private: // only on root process (rank 0)

    mesh_type_root M_mesh_root;
    mesh_type_root M_mesh_init_root;
    mesh_type_root M_mesh_previous_root;

    std::vector<double> M_UM_root;
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
    boost::unordered_map<int, std::array<double,2>> M_drifter; // Drifters are kept in an unordered map containing number and coordinates
    std::fstream M_iabp_file;             // The file we read the IABP buoy data from

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
    void constantAtmosphere();
    void constantOcean();
    void constantIce();
	void constantBathymetry();

    void equallySpacedDrifter();
    void outputDrifter(std::fstream &iabp_out);
    void initIABPDrifter();
    void updateIABPDrifter();

    void asrAtmosphere(bool reload);
    void topazOcean(bool reload);
    void topazIce();
	void etopoBathymetry(bool reload);

    void loadDataset(Dataset *dataset);//(double const& u, double const& v);
    void loadGrid(Grid *grid);
};
} // Nextsim
#endif
