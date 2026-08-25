/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   interpolation.hpp
 * @author Fabien Salmon <fabien.salmon@inria.fr>
 * @date   Thu Aug  4 09:47:27 CEST 2024
 */

#ifndef __interpolation_H
#define __interpolation_H

#include <vector>
#include <gmshmeshseq.hpp>
#include <MemOps.h>

/**
 * @brief Different functions allowing interpolations between two meshes
 *
 * @see
 *
 */

namespace Nextsim
{

// Check if a point is inside a triangle
int inside(std::vector<double> const& points, double xp, double yp);

// Compute the area of a triangle
double triangle_area(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &c);

// Function to interpolate fields on nodes
int InterpFromMeshToMesh2d(double** interp_out, int* background_triangles, double* background_x, double* background_y,
                           int n_background_nodes, int n_background_triangles, double* interp_in, int source_size, int nb_var,
                           double* x, double* y, int n_nodes, bool isdefault=false, double defaultvalue=1e-24);

void ConservativeRemappingWithWeights(GmshMeshSeq const& mesh, std::vector<double> &gridX, std::vector<double> &gridY, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY,
                                  std::vector<int> &gridP, std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights);

// Remapping from mesh to mesh.
// In this case we want to both calculate weights and apply them in the same step
void ConservativeRemappingFromMeshToMesh(double* &interp_out, std::vector<double>& interp_in, int nb_var,
                                         std::vector<int>& indexTr, std::vector<double>& coordX, std::vector<double>& coordY,
                                         const std::vector<int>& new_indexTr, const std::vector<double>& new_coordX, const std::vector<double>& new_coordY, const std::vector<double>& UM);

// Apply weights for a mesh-to-grid remapping
void ConservativeRemappingMeshToGrid(double* &interp_out, std::vector<double> &interp_in, int const nb_var, int grid_size, double miss_val,
                                     std::vector<int> &gridP, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY,
                                     std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights, int num_corners);

// Recursive function to check the current triangle
void checkTriangle(std::vector<int> &indexTr, std::vector<double> &coordX, std::vector<double> &coordY, 
                   std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY, int current_number, // inputs
                   std::vector<std::vector<int>> const& list_triangles, std::vector<std::vector<int>> const& list_neighbours, // inputs
                   std::vector<int> &number, std::vector<double> &weight);  // outputs

// Check if we've already visited this triangle
bool visited(int current_triangle, std::vector<int> const &triangles);

// Check if points are inside polygon
bool checkIfInside(const std::vector<double>& vertx, const std::vector<double>& verty, double testx, double testy, bool inclusive);

// Calculate the area
double area(std::vector<std::pair<double,double>>& points);

// Check for intersection and add intersecting points to the list
// https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
bool checkIfIntersecting(double X, double Y, double Xprev, double Yprev, const std::vector<double>& gridCornerX, const std::vector<double>& gridCornerY, // inputs
                         std::vector<std::pair<double,double>>& points); // sideffects

// Sort the points in a clock wise manner around the centre of the cloud of points
void sortClockwise(std::vector<std::pair<double,double>>& points);

// The custom sorting function
// https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order <- the second code exmple
bool CWSort(std::pair<double,double> p1, std::pair<double,double> p2);

}
#endif
