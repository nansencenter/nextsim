/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   metric.hpp
 * @author Fabien Salmon <fabien.salmon@inria.fr>
 * @date   Thu Aug  4 09:47:27 CEST 2024
 */

#ifndef __metric_H
#define __metric_H

#include <vector>
#include <math.h>
#include <numeric>
#include <limits>
#include <algorithm>
#include <gmshmesh.hpp>
#include <gmshmeshseq.hpp>

/**
 * @class metric
 * @brief Compute the metric for anisotropic remeshing with mmg
 *
 * @see
 *
 */

namespace Nextsim
{
class Metric
{
public:

    double quality_threshold = 0.; // Limit of quality that triggers the remeshing process
    int p = 2; // p of norm L^p
    int Nst = 1; // Number of points 
    bool is_vtk_written = false; // Set to true to write a vtk file containing the mesh and the metric. It is true when the debug option for MMG is activated
    bool isotropic = true; // Set to true to have an isotropic metric (Salmon's metric), set to false to have an anisotropic metric (Loseille's metric)
    double scale_factor_min = 1.; // Scale factor for the minimum size of the elements (scale_factor * size of the smallest initial element)
    double scale_factor_max = 1.; // Scale factor for the maximum size of the elements (scale_factor * size of the greatest initial element)
    int process;

    // Edge length in the non-euclidean geometry defined by the metric
    double metric_square_edge_length(int id1, int id2, std::vector<double> const pt1, std::vector<double> const pt2) const;

    // Areas of the triangles in the non-euclidean geometry
    template<typename FEMeshType>
    std::vector<double> metric_triangle_areas(FEMeshType const& mesh) const;

    // Quality function for every triangle
    template<typename FEMeshType>
    std::vector<double> triangle_quality(FEMeshType const& mesh) const;

    // Read the node coordinates according to the type of mesh
    std::vector<double> nodes_coordinates(GmshMesh const& mesh, int const n_triangle, int const i_index) const
    {return mesh.nodes().find(mesh.triangles()[n_triangle].indices[i_index])->second.coords;}
    std::vector<double> nodes_coordinates(GmshMeshSeq const& mesh, int const n_triangle, int const i_index) const
    {return mesh.nodes()[mesh.triangles()[n_triangle].indices[i_index]-1].coords;}

    // Return a list of triangles that do not respect the quality tolerance
    std::vector<int> check_triangle_quality(GmshMesh const& mesh) const;
    std::vector<int> check_triangle_quality(GmshMeshSeq const& mesh) const;

    // Metric solution of the minimization problem
    void compute_optimal_metric(GmshMesh const& mesh, std::vector<double> const& field, double hmin, double hmax); 
    void compute_optimal_metric(GmshMeshSeq const& mesh, std::vector<double> const& field, double hmin, double hmax);

    double error_min_max(std::vector<double> const& areas_vertex, std::vector<int> const& order_number,
                         std::vector<double> const& sum_gamma, int nb_vertices, double d_min, double d_max, double alpha, double z1);
    double error_max(std::vector<double> const& areas_vertex, std::vector<int> const& order_number,
                     std::vector<double> const& sum_gamma, int nb_vertices, double d_max, double alpha, double z2);
    void compute_optimal_bounded_isotropic_metric(GmshMesh const& mesh, std::vector<double> const& field, double hmin, double hmax, Communicator const & comm);
    void compute_optimal_bounded_isotropic_metric(GmshMeshSeq const& mesh, std::vector<double> const& field, double hmin, double hmax, Communicator const & comm);

    // Compute the areas of each triangle of the mesh
    template<typename FEMeshType>
    std::vector<double> compute_triangle_euclidean_areas(FEMeshType const& mesh) const ; 

    // Compute the Hessian matrix at each vertex for a chosen field
    template<typename FEMeshType>
    std::vector<std::vector<double>> compute_hessian_field(FEMeshType const& mesh, std::vector<double> const& field) const ;

    // Compute the determinant of a matrix (Hessian) associated to each vertex
    std::vector<double> compute_determinant_field(std::vector<std::vector<double>> const& matrix) const ;

    // Compute an average over several triangles with weights associated to their areas
    double area_weighted_average(std::vector<double> const& scalar_field, std::vector<double> const& area, std::vector<int> const& list_triangles) const;

    // Compute the gradient of a field over a triangle
    std::vector<double> compute_grad(std::vector<double> const& scalar_field, std::vector<std::vector<double>> const& points_coord, double const area) const;

    // Compute the Hessian of a field over a triangle
    std::vector<double> compute_hessian(std::vector<std::vector<double>> const& gradient_field, std::vector<std::vector<double>> const& points_coord, double const area) const;

    // Diagonalize the Hessian and replace the eigen values by their absolute value to make the Hessian symmetric positive definite
    std::vector<std::vector<double>> diagonalize_hessian_positive(std::vector<std::vector<double>> const& matrix) const;

    // Compute the list of triangles containing each vertex
    template<typename FEMeshType>
    std::vector<std::vector<int>> compute_list_triangles(FEMeshType const& mesh) const;

    // Write a vtk file with the mesh nodes and the metric field
    template<typename FEMeshType>
    void write_metric_vtk(FEMeshType const& mesh) const;

    std::vector<std::vector<double>> components; // Components of the metric at each vertex

    // Constructors
    Metric(){}

    Metric(int i1, int i2)
    :
        p(i1),
        Nst(i2)
    {}

    Metric(int i1) // Norm L^1 by default
    :
        Nst(i1)
    {}

    // Destructor
    ~Metric(){}

};
}//Nextsim
#endif
