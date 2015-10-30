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
#include <boost/version.hpp>
#include <boost/format.hpp>
#include <BamgConvertMeshx.h>
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
    std::vector<double> sides(element_type const& element, mesh_type const& mesh) const;
    std::vector<double> minMaxSide(mesh_type const& mesh) const;
    void movedMesh(std::vector<double> const& um, double factor = 0);
    double measure(element_type const& element, mesh_type const& mesh) const;
    std::vector<double> shapeCoeff(element_type const& element, mesh_type const& mesh) const;

    void assemble();
    void solve();
    void run();
    void error();

    double minAngles(element_type const& element, mesh_type const& mesh) const;
    double minAngle(mesh_type const& mesh) const;

    double resolution(mesh_type const& mesh) const;

    std::vector<double> hminVertices(mesh_type const& mesh, BamgMesh const* bamgmesh) const;
    std::vector<double> hmaxVertices(mesh_type const& mesh, BamgMesh const* bamgmesh) const;

    void performSimulation();
    void initialConditions();
    void initBamg();
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

    void PwlInterp2D();
    void importBamg();

private:
    po::variables_map vm;
    mesh_type M_mesh;
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

    std::vector<int> dirichlet_flags;
    boost::timer chrono;

    forcing::WindType M_wind;
    forcing::OceanType M_ocean;
    forcing::ThermoType M_thermo;

    forcing::ConcentrationType M_conc;
    forcing::ThicknessType M_thick;
    forcing::DamageType M_damage;
    forcing::SnowThicknessType M_snow_thick;
    forcing::DrifterType M_drifter;

    std::vector<double> M_sigma;
    std::vector<double> M_surface;
    std::vector<double> M_divergence_rate;
    std::vector<double> M_UT;
    std::vector<double> M_UM;
    std::vector<double> M_VT;
    std::vector<double> M_VTM;
    std::vector<double> M_VTMM;

private:

    // bamgopt_ptrtype bamgopt;
    // bamgmesh_ptrtype bamgmesh;
    // bamggeom_ptrtype bamggeom;

    // bamgmesh_ptrtype bamgmeshout;
    // bamggeom_ptrtype bamggeomout;

    BamgOpts *bamgopt;
    BamgMesh *bamgmesh;
    BamgGeom *bamggeom;

    BamgMesh *bamgmeshout;
    BamgGeom *bamggeomout;

private:
    std::vector<double> wind;
    std::vector<double> ocean;
    std::vector<double> thermo;

    std::vector<double> conc;
    std::vector<double> thick;
    std::vector<double> damage;
    std::vector<double> snow_thick;
    std::vector<double> drifter;

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
