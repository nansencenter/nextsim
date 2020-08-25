/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/*!\file:  InterpFromMeshToMesh2dx.h
 * \brief header file for Bamg module
 */

#ifndef _CONSERVATIVEREMAPPING_H
#define _CONSERVATIVEREMAPPING_H

#include "Mesh.h"
#include <stdexcept>
#include <algorithm>
#include <fstream>

// Calculate the weights
void ConservativeRemappingWeights(BamgMesh* bamgmesh, std::vector<double> &gridX, std::vector<double> &gridY, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY,
        std::vector<int> &gridP, std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights);

// Apply weights for a mesh-to-grid remapping
void ConservativeRemappingMeshToGrid(double* &interp_out, std::vector<double> &interp_in, int const nb_var, int grid_size, double miss_val,
        std::vector<int> &gridP, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY,
        std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights);

// Apply weights going from grid to mesh
void ConservativeRemappingGridToMesh(double* &interp_out, std::vector<double> &interp_in, int const nb_var, int const numElements, std::vector<int> &gridP, std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights);

// Recursive function to check the current triangle
inline void checkTriangle(BamgMesh* bamgmesh, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY, int current_number, // inputs
		std::vector<int> &number, std::vector<double> &weight);  // outputs

// Check if we've already visited this triangle
inline bool visited(int current_triangle, std::vector<int> const &triangles);

// Check if points are inside polygon
inline bool checkIfInside(std::vector<double> const &vertx, std::vector<double> const &verty, double testx, double testy);

// Check for intersection and add intersecting points to the list
inline bool checkIfIntersecting(double X, double Y, double Xnext, double Ynext, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY, // inputs
		std::vector<std::pair<double,double>> &points); // side-effect

// Calculate the area
inline double area(std::vector<std::pair<double,double>> &points);

// Sort the points in a clock wise manner
inline void sortClockwise(std::vector<std::pair<double,double>> &points);

// The custom sorting function
inline bool CWSort(std::pair<double,double> p1, std::pair<double,double> p2);

#endif
