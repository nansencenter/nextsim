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
void ConservativeRemappingFromMeshToMesh(double* &interp_out, std::vector<double> &interp_in, int nb_var, std::vector<int> &indexTr, std::vector<double> &coordX, std::vector<double> &coordY,
                                         std::vector<int> &new_indexTr, std::vector<double> &new_coordX, std::vector<double> &new_coordY);

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

double compute_area(double a0, double a1, double b0, double b1, double c0, double c1);
std::vector<double> interpolationNode(int* indexTria, double* coordX, double* coordY, int numTriangles, int numNodes, int M_nb_var_node, double x, double y, double* data_in);
void interpolationNodes(double** interp_out, int* indexTria_prv, double* coordX_prv, double* coordY_prv, int numTriangles_prv,                                        
                        int numNodes_prv, int M_nb_var_node, double* coordX, double* coordY, int numNodes, double* data_in);

std::vector<std::vector<std::vector<int>>> grid_size_triangles(double* coordX, double* coordY, int* indexTria, int numTriangles, double minX, double minY, double maxX, double maxY, int GRID_SIZE);
void find_bounding_box(double* coordX, double* coordY, int numNodes, double* xmin, double* xmax, double* ymin, double* ymax);
std::vector<double> find_segment_intersection(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, int* is_intersected);
void interpolationElements(double** interp_out, int* indexTria_prv, double* coordX_prv, double* coordY_prv, int numTriangles_prv,
                           int numNodes_prv, int nb_var_element, double* coordX, double* coordY, int* indexTria, int numTriangles, double* data_in);
double cross_product(std::vector<double> v1, std::vector<double> v2);
int is_inside(std::vector<double> Point, std::vector<std::vector<double>> Triangle);
std::vector<std::vector<double>> sort_points(std::vector<std::vector<double>> list_intersection_points);
double intersection_area(std::vector<std::vector<double>> Tria1, std::vector<std::vector<double>> Tria2);
}
#endif
