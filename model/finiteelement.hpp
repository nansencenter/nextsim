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

#include <datasetmanager.hpp>
#include <variablemanager.hpp>
#include <optionhandler.hpp>
#include <meshhandler.hpp>
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

class FiniteElement : public VariableManager, public DataSetManager
{
public:

    typedef typename GmshMesh::bimap_type bimap_type;

    typedef boost::shared_ptr<graphmpi_type> graphmpi_ptrtype;

    // typedef ExternalData external_data;
    // typedef ExternalData::Dataset Dataset;
    // typedef ExternalData::Grid Grid;
    // typedef ExternalData::Dimension Dimension;
    // typedef ExternalData::Variable Variable;
    // typedef ExternalData::Vectorial_Variable Vectorial_Variable;
    // typedef boost::ptr_vector<external_data> externaldata_ptr_vector;

    FiniteElement();

    // FiniteElement(Communicator const& comm = Environment::comm());

    mesh_type const& mesh() const {return M_mesh;}

    std::vector<double> shapeCoeff(element_type const& element) const;

    bool checkRegridding();
    void regrid(bool step = true);

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

    template<typename FEMeshType>
    bool flip(FEMeshType const& mesh, std::vector<double> const& um, double factor) const;

    void initOptAndParam();
    void initFETensors();
    void forcingNesting();

    void assimilateIce();
    void assimilateSlabOcean();
    void initIce();
    void checkConsistency();
    void initSlabOcean();

    //void timeInterpolation(int step);
    void nodesToElements(double const* depth, std::vector<double>& v);

    void PwlInterp2D();
    void calcAuxiliaryVariables();
    void initModelState();
    void DataAssimilation();

    void calcCohesion();
    void updateFreeDriftVelocity();
    void speedScaling(std::vector<double>& speed_scaling);
    void update(std::vector<double> const & UM_P);
    void updateSigmaDamage(double const dt);

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

private:
    po::variables_map vm;
    mesh_type M_mesh_init;

    std::vector<element_type> M_edges;

    int M_rank;

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

    std::vector<int> M_boundary_flags;

    boost::mpi::timer chrono, chrono_tot;
    Timer::timer M_timer;

    setup::IceType M_ice_type;
    setup::BasalStressType M_basal_stress_type;
    setup::ThermoType M_thermo_type;
    setup::DynamicsType M_dynamics_type;

#ifdef AEROBULK
    aerobulk::algorithm M_ocean_bulk_formula;
#endif

    setup::FreezingPointType M_freezingpoint_type;
    setup::OceanHeatfluxScheme M_Qio_type;
    setup::IceCategoryType M_ice_cat_type;
    //fsd related
    setup::WeldingType M_welding_type    ;
    setup::BreakupType M_breakup_type    ;
    setup::FSDType M_fsd_type    ;

    bool M_flooding;

    // diffusivity parameters
    std::vector<double> M_diffusivity_parameters;

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

    std::vector<double> M_fcor;

    std::vector<double> M_Dunit;
    std::vector<std::vector<double>> M_shape_coeff;
    std::vector<std::vector<double>> M_B0T;

    // =============================================================================
    // variables needed for nesting
    bool M_use_nesting;
    std::string M_nudge_function;
    double M_nudge_timescale;
    double M_nudge_lengthscale;
    bool M_nest_dynamic_vars;
    // =============================================================================

private:

    double nu0;
    double young;
    double rhoi;
    double rhos;
    double const days_in_sec  = 86400.;
    int output_time_step;
    double time_init;
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

private: // only on root process (rank 0)

    std::vector<int> M_connectivity_root;

    std::vector<std::vector<double>> M_B0T_root;

private:
    // Drifters
    std::vector<Drifters> M_drifters;// vector of all the Drifters objects (including IABP ones)
    std::vector<int> M_osisaf_drifters_indices;// indices of OSISAF drifters in M_drifters

    // Element variable
    std::vector<double> M_element_age;         // Age of the element (model time since its last adaptation)

private:
    // Variables for the moorings

    std::vector<double> M_conc_mean;    // Mean concentration (on the mesh)
    std::vector<double> M_thick_mean;   // Mean ice thickness (on the mesh)
    std::vector<double> M_snow_thick_mean;  // Mean snow thickness (on the mesh)
    std::vector<double> M_VT_mean;      // Mean velocity (on the mesh)
    std::vector<double> M_fyi_fraction_mean;  // Fraction of the first year ice (FYI) (on the mesh)
    std::vector<double> M_age_det_mean;       // Ice age observable from space (area weighted) [timestep] (on the mesh)
    std::vector<double> M_age_mean;           // Effective ice age [timestep] (on the mesh)
    std::vector<double> M_conc_myi_mean;  // Mean concentration of multiyear ice (MYI) (on the mesh)
    std::vector<double> M_thick_myi_mean;  // Mean thickness of multiyear ice (MYI) (on the mesh)
    std::vector<double> M_freeze_days_mean;  // Mean number of consecutive days freezing has been occurring (on the mesh)
    std::vector<double> M_freeze_onset_mean;  // 1 if freezeing has been occurring (on the mesh)

    std::vector<double> M_conc_grid;    // Mean concentration (on the grid)
    std::vector<double> M_thick_grid;   // Mean ice thickness (on the grid)
    std::vector<double> M_snow_thick_grid;  // Mean snow thickness (on the grid)
    std::vector<double> M_VT_grid;      // Mean velocity (on the grid)
    std::vector<double> M_fyi_fraction_grid;  // Fraction of the first year ice (FYI) (on the grid)
    std::vector<double> M_age_det_grid;       // Ice age observable from space (area weighted) [timestep] (on the grid)
    std::vector<double> M_age_grid;           // Effective ice age [timestep] (on the grid)
    // NDGB: Not sure the lines below are needed, what are these grid useful?
    std::vector<double> M_conc_myi_grid;  // Mean concentration of multiyear ice (MYI) (on the grid)
    std::vector<double> M_thick_myi_grid;  // Mean thickness of multiyear ice (MYI) (on the grid)
    std::vector<double> M_freeze_days_grid;  // Mean number of consecutive days freezing has been occurring (on the grid)
    std::vector<double> M_freeze_onset_grid;  // 1 if freezing has been occurring (on the grid)

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
    void glorys12Ice();

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
