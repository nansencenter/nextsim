/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   finiteelement.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug 24 10:40:29 2015
 */

#ifndef __FiniteElement_HPP
#define __FiniteElement_HPP 1

#include <solverpetsc.hpp>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <boost/version.hpp>
#include <boost/format.hpp>
#include <BamgConvertMeshx.h>
#include <BamgTriangulatex.h>
#include <Bamgx.h>
#include <InterpFromMeshToMesh2dx.h>
#include <pwl_interp_2d_scattered.hpp>
#include <gmshmesh.hpp>
#include "enums.hpp"

namespace Nextsim
{
class FiniteElement
{
public:

    typedef typename GmshMesh::point_type point_type;
    typedef typename GmshMesh::element_type element_type;

    typedef GmshMesh mesh_type;
    typedef MatrixPetsc matrix_type;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;
    typedef VectorPetsc vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

    FiniteElement();

    mesh_type const& mesh() const {return M_mesh;}

    matrix_ptrtype const& matrixPtr() const {return M_matrix;}
    vector_ptrtype const& rhsPtr() const {return M_vector;}
    vector_ptrtype const& solutionPtr() const {return M_solution;}

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

    void assemble();
    void solve();
    void run();
    void error();

    double minAngles(element_type const& element, mesh_type const& mesh) const;
    double minAngle(mesh_type const& mesh) const;

    double minAngle(mesh_type const& mesh, std::vector<double> const& um, double factor) const;

    bool flip(mesh_type const& mesh, std::vector<double> const& um, double factor) const;

    double resolution(mesh_type const& mesh) const;

    std::vector<double> hminVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh) const;
    std::vector<double> hmaxVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh) const;

    void computeFactors(int cpt);
    void initBamg();
    void initConstant();
    void forcing();
    void forcingWind(double const& u, double const& v);
    void forcingOcean(double const& u, double const& v);
    void forcingThermo(double const& u, double const& v);

    void initConcentration();
    void initThickness();
    void initDamage();
    void initSnowThickness();
    void initThermodynamics();
    void initDrifter();
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
    void exportResults(int step);

private:
    po::variables_map vm;
    mesh_type M_mesh;
    mesh_type M_mesh_init;
    mesh_type M_mesh_previous;
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

    std::vector<int> M_dirichlet_flags;
    std::vector<int> M_neumann_flags;

    boost::timer chrono;

    forcing::WindType M_wind_type;
    forcing::OceanType M_ocean_type;
    forcing::ThermoType M_thermo_type;

    forcing::ConcentrationType M_conc_type;
    forcing::ThicknessType M_thick_type;
    forcing::DamageType M_damage_type;
    forcing::SnowThicknessType M_snow_thick_type;
    forcing::DrifterType M_drifter_type;

    std::vector<double> M_sigma;
    std::vector<double> M_surface;
    std::vector<double> M_divergence_rate;
    //std::vector<double> M_UT;
    std::vector<double> M_UM;
    std::vector<double> M_VT;
    std::vector<double> M_VTM;
    std::vector<double> M_VTMM;
    std::vector<double> M_Vcor;

    std::vector<double> M_bathy_depth;
    std::vector<double> M_element_depth;
    std::vector<double> M_Vair_factor;
    std::vector<double> M_Voce_factor;
    std::vector<double> M_basal_factor;
    std::vector<double> M_water_elements;
    std::vector<double> M_ssh;
    std::vector<double> M_element_ssh;
    std::vector<double> M_h_thin;
    std::vector<double> M_hs_thin;
    std::vector<double> M_h_ridged_thin_ice;
    std::vector<double> M_h_ridged_thick_ice;

    std::vector<double> M_norm_Voce_ice;
    std::vector<double> M_norm_Vair_ice;
    std::vector<double> M_norm_Vice;
    std::vector<double> M_fcor;

    std::vector<double> M_Dunit;
    std::vector<double> M_Dunit_comp;
    std::vector<double> M_Mass;
    // std::vector<std::vector<double>> M_B0T;
    // std::vector<std::vector<double>> M_B0T_Dunit_B0T;
    // std::vector<std::vector<double>> M_B0T_Dunit_comp_B0T;
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
    double ssh_coef;
    double time_relaxation_damage;
    double water_depth;
    double critical_h;
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

    Options *options;

private:
    std::vector<double> M_wind;
    std::vector<double> M_ocean;
    std::vector<double> M_thermo;

    std::vector<double> M_conc;
    std::vector<double> M_thick;
    std::vector<double> M_damage;
    std::vector<double> M_snow_thick;
    std::vector<double> M_drifter;

private:
    void constantWind(double const& u, double const& v);
    void constantOcean(double const& u, double const& v);
    void constantThermo(double const& u, double const& v);

    void constantConc();
    void constantThick();
    void constantDamage();
    void constantSnowThick();
    void equallySpacedDrifter();

};
} // Nextsim
#endif
