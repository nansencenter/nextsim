/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   interpolation.hpp
 * @author Fabien Salmon <fabien.salmon@inria.fr>
 * @date   Thu Aug  4 09:47:27 CEST 2024
 */

#ifndef __iterpolation_H
#define __interpolation_H

#include <vector>
#include <InterpFromMeshToMesh2dx.h>
#include <gmshmeshseq.hpp>

/**
 * @brief Different functions allowing interpolations between two meshes
 *
 * @see
 *
 */

namespace Nextsim
{

void ConservativeRemappingWithWeights(GmshMeshSeq const& mesh, std::vector<double> &gridX, std::vector<double> &gridY, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY,
                                  std::vector<int> &gridP, std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights);

// Remapping from mesh to mesh.
// In this case we want to both calculate weights and apply them in the same step
void ConservativeRemappingFromMeshToMesh(double* &interp_out, std::vector<double> &interp_in, int nb_var, GmshMeshSeq const& mesh_old, GmshMeshSeq const& mesh_new);

// Apply weights for a mesh-to-grid remapping
void ConservativeRemappingMeshToGrid(double* &interp_out, std::vector<double> &interp_in, int const nb_var, int grid_size, double miss_val,
                                     std::vector<int> &gridP, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY,
                                     std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights, int num_corners);

// Recursive function to check the current triangle
void checkTriangle(GmshMeshSeq const& mesh, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY, int current_number, // inputs
                   std::vector<std::vector<int>> const& list_triangles, std::vector<std::vector<int>> const& list_neighbours, // inputs
                   std::vector<int> &number, std::vector<double> &weight);  // outputs

// Check if we've already visited this triangle
bool visited(int current_triangle, std::vector<int> const &triangles);

// Check if points are inside polygon
bool checkIfInside(std::vector<double> const &vertx, std::vector<double> const &verty, double testx, double testy);

// Calculate the area
double area(std::vector<std::pair<double,double>> &points);

// Check for intersection and add intersecting points to the list
// https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
bool checkIfIntersecting(double X, double Y, double Xprev, double Yprev, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY, // inputs
                         std::vector<std::pair<double,double>> &points);

// Sort the points in a clock wise manner around the centre of the cloud of points
void sortClockwise(std::vector<std::pair<double,double>> &points);

// The custom sorting function
// https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order <- the second code exmple
bool CWSort(std::pair<double,double> p1, std::pair<double,double> p2);

}
#endif
