/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   ConservativeRemapping.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Tue 24 Apr 10:58:29 CEST 2018
 */

#include "./InterpFromMeshToMesh2dx.h"
#include "ConservativeRemapping.hpp"

using namespace bamg;

// Calculate the weights
void ConservativeRemappingWeights(BamgMesh* bamgmesh, std::vector<double> &gridX, std::vector<double> &gridY, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY,
        std::vector<int> &gridP, std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights)
{
    // ---------- Initialisation ---------- //
    // Copy the triangle information
    int numTriangles = bamgmesh->TrianglesSize[0];
    std::vector<int> indexTr(3*numTriangles);
    std::vector<double> elnum(numTriangles);
    for (int tr=0; tr<numTriangles; ++tr)
    {
        // NB: Maintain bamg numbering for the call to InterpFromMeshToMesh (i.e. the first triangle is numbered 1 so we don't do Triangles[...]-1)
        indexTr[3*tr  ] = bamgmesh->Triangles[4*tr];
        indexTr[3*tr+1] = bamgmesh->Triangles[4*tr+1];
        indexTr[3*tr+2] = bamgmesh->Triangles[4*tr+2];

        // NB: Here we use the C/C++ numbering (i.e. the first triangle is numbered 0)
        elnum[tr] = tr;
    }

    // Copy the node information
    int numNodes     = bamgmesh->VerticesSize[0];
    std::vector<double> coordX(numNodes);
    std::vector<double> coordY(numNodes);
    for (int id=0; id<numNodes; ++id)
    {
        coordX[id] = bamgmesh->Vertices[3*id];
        coordY[id] = bamgmesh->Vertices[3*id+1];
    }

    // Copy the grid information
    int grid_size = gridX.size();
    assert(grid_size==gridY.size());

    // Reset the size of gridP, triangles, and weights
    // Doing this and using push_back is probably sub-optimal, but that's an
    // optimisation for another day
    gridP.resize(0);
    triangles.resize(0);
    weights.resize(0);

    // Find which element each P-point hits - use -1 for land points
    double* elnum_out;
    InterpFromMeshToMesh2dx(&elnum_out,
                &indexTr[0],&coordX[0],&coordY[0],
                numNodes,numTriangles,
                &elnum[0],
                numTriangles,1,
                &gridX[0],&gridY[0],grid_size,
                true, -1);

    // Calculate weights
    for (int ppoint=0; ppoint<grid_size; ++ppoint)
    {
        // Carefully take the right integer value for element number
        int i_elnum_out = std::round(elnum_out[ppoint]);

        // Don't do anything for land points
        if ( i_elnum_out >= 0 )
        {
            // Save the ppoints
            gridP.push_back(ppoint);

            // coordinates and size of the grid cell
            std::vector<double> cornerX(4);
            std::vector<double> cornerY(4);
            for (int corner=0; corner<4; ++corner)
            {
                cornerX[corner] = gridCornerX[4*ppoint+corner];
                cornerY[corner] = gridCornerY[4*ppoint+corner];
            }

            // get the numbers for the elements that contribute and the weight of the contribution
            std::vector<int> local_triangles;
            std::vector<double> local_weights;
            // When using a "split" mesh reserving 10 covers about 95% of cases, while reserving 20 covers just about 100%
            // local_triangles.reserve(10);
            // local_weights.reserve(10);

            // Call the recursive function (this is our work horse here)
            checkTriangle(bamgmesh, cornerX, cornerY, i_elnum_out, local_triangles, local_weights);

            // Save the weights and triangle numbers
            triangles.push_back(local_triangles);
            weights.push_back(local_weights);
        }
    }
    xDelete<double>(elnum_out);
}

// Incremental remapping using ConservativeRemappingWeights and ConservativeRemappingMeshToGrid
void IncrementalRemapping(double* &interp_out, std::vector<double> &interp_in, int const nb_var, BamgMesh* bamgmesh, std::vector<double> &UM)
{
    // Copy the node information
    int numNodes = bamgmesh->VerticesSize[0];
    assert(UM.size() == 2*numNodes);

    // Copy the triangle information of the mesh
    // Keep the nomenclature for a grid (even if it's a mesh)
    int grid_size = bamgmesh->TrianglesSize[0];

    // Initialise gridP, triangles, and weights
    std::vector<int> gridP(grid_size);
    std::vector<std::vector<int>> triangles(grid_size);
    std::vector<std::vector<double>> weights(grid_size);

    // Calculate weights
    for (int ppoint=0; ppoint<grid_size; ++ppoint)
    {
        // Save the ppoints
        gridP[ppoint] = ppoint;

        // vertices of the moved mesh in a format checkTriangle understands
        std::vector<double> cornerX(3);
        std::vector<double> cornerY(3);
        for (int i=0; i<3; ++i)
        {
            int id = bamgmesh->Triangles[4*ppoint+i] - 1; // Here we need C/C++ numbering
            cornerX[i] = bamgmesh->Vertices[3*id];
            cornerY[i] = bamgmesh->Vertices[3*id+1];

            /* Apply UM only if it's sufficiently large. The intersection and
             * point-in-polygon algorithms don't neccesarily agree when the
             * displacement is close to the rounding error. We also set UM to
             * zero only if it's large enough. That way, no information is lost
             * here. */
            if ( std::hypot(UM[id],UM[id+numNodes]) > 1e-3 )
            {
                cornerX[i] += UM[id];
                cornerY[i] += UM[id+numNodes];
                UM[id] = 0.;
                UM[id+numNodes] = 0.;
            }
        }

        // Call the recursive function (this is our work horse here)
        // Here element number is the same as ppoint
        checkTriangle(bamgmesh, cornerX, cornerY, ppoint, triangles[ppoint], weights[ppoint]);
        assert(weights[ppoint][0] > 0.);
    }

    // Now we apply the weights with GridToMesh
    ConservativeRemappingGridToMesh(interp_out, interp_in, nb_var, grid_size, gridP, triangles, weights);
}

// Apply weights for a mesh-to-grid remapping
void ConservativeRemappingMeshToGrid(double* &interp_out, std::vector<double> &interp_in, int const nb_var, int grid_size, double miss_val,
        std::vector<int> &gridP, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY,
        std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights,
        int num_corners)  // Default for num_corners is 4
{
    // Initialise interp_out
    interp_out = xNew<double>(nb_var*grid_size);
    for (int i=0; i<nb_var*grid_size; ++i)
        interp_out[i] = miss_val;

    // Apply the weights
    for (int i=0; i<gridP.size(); ++i)
    {
        // Calculate cell area
        int ppoint = gridP[i];
        std::vector<std::pair<double,double>> points(num_corners);
        for (int corner=0; corner<num_corners; ++corner)
            points[corner] = std::make_pair(gridCornerX[num_corners*ppoint+corner],gridCornerY[num_corners*ppoint+corner]);

        double r_cell_area = 1./area(points);

        // walk through all the variables
        for (int var=0; var<nb_var; ++var)
        {
            // ... and contributing elements
            interp_out[gridP[i]*nb_var+var] = 0;
            for (int tr=0; tr<triangles[i].size(); ++tr)
                interp_out[gridP[i]*nb_var+var] += interp_in[triangles[i][tr]*nb_var+var]*weights[i][tr];

            interp_out[gridP[i]*nb_var+var] *= r_cell_area;
        }
    }
}

// Apply weights going from grid to mesh
void ConservativeRemappingGridToMesh(double* &interp_out, std::vector<double> &interp_in, int const nb_var, int const numElements, std::vector<int> &gridP, std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights)
{
    assert(weights.size()==triangles.size());

    interp_out = xNew<double>(nb_var*numElements);

    // Initialise interp_out and triangle_area
    std::vector<double> triangle_area(numElements,0.);
    for (int tr=0; tr<numElements; ++tr)
        for (int var=0; var<nb_var; ++var)
            interp_out[tr*nb_var+var] = 0;

    // Walk through wet grid cells
    for (int i=0; i<gridP.size(); ++i)
    {
        // Apply the weights: walk through the contributing elements
        for (int j=0; j<triangles[i].size(); ++j)
        {
            int const tr = triangles[i][j];
            triangle_area[tr] += weights[i][j];
            // ... and all the variables
            for (int var=0; var<nb_var; ++var)
                interp_out[tr*nb_var+var] += interp_in[gridP[i]*nb_var+var]*weights[i][j];
        }
    }

    // Post-process interp_out
    for (int tr=0; tr<numElements; ++tr)
    {
        double const r_triangle_area = 1./triangle_area[tr];
        for (int var=0; var<nb_var; ++var)
            interp_out[tr*nb_var+var] *= r_triangle_area;
    }
}

// Remapping from mesh to mesh.
// In this case we want to both calculate weights and apply them in the same step
// Drop-in-replacement for InterpFromMeshToMesh2dCavities
void ConservativeRemappingMeshToMesh(double* &interp_out, std::vector<double> &interp_in, int nb_var,
      BamgMesh* bamgmesh_old, BamgMesh* bamgmesh_new)
{
    // We start off the same as ConservativeRemappingWeights - only here the new mesh replaces the grid

    // ---------- Initialisation ---------- //
    // Copy the triangle information of the _old mesh
    int numTriangles = bamgmesh_old->TrianglesSize[0];
    std::vector<int> indexTr(3*numTriangles);
    std::vector<double> elnum(numTriangles);
    for (int tr=0; tr<numTriangles; ++tr)
    {
        // NB: Maintain bamg numbering for the call to InterpFromMeshToMesh (i.e. the first triangle is numbered 1 so we don't do Triangles[...]-1)
        indexTr[3*tr  ] = bamgmesh_old->Triangles[4*tr];
        indexTr[3*tr+1] = bamgmesh_old->Triangles[4*tr+1];
        indexTr[3*tr+2] = bamgmesh_old->Triangles[4*tr+2];

        // NB: Here we use the C/C++ numbering (i.e. the first triangle is numbered 0)
        elnum[tr] = tr;
    }

    // Copy the node information
    int numNodes     = bamgmesh_old->VerticesSize[0];
    std::vector<double> coordX(numNodes);
    std::vector<double> coordY(numNodes);
    for (int id=0; id<numNodes; ++id)
    {
        coordX[id] = bamgmesh_old->Vertices[3*id];
        coordY[id] = bamgmesh_old->Vertices[3*id+1];
    }

    // Copy the triangle information of the _new mesh and calculate the barycentre
    // Keep the nomenclature for a grid (even if it's a mesh)
    int grid_size = bamgmesh_new->TrianglesSize[0];
    std::vector<double> gridX(grid_size);
    std::vector<double> gridY(grid_size);
    std::vector<double> gridCornerX(3*grid_size);
    std::vector<double> gridCornerY(3*grid_size);
    for (int tr=0; tr<grid_size; ++tr)
    {
        gridX[tr] = 0; // barycentre
        gridY[tr] = 0;
        for (int i=0; i<3; ++i)
        {
            // grid corner == vertice
            int id = bamgmesh_new->Triangles[4*tr+i] - 1; // Here we need C/C++ numbering
            gridCornerX[3*tr+i] = bamgmesh_new->Vertices[3*id];
            gridCornerY[3*tr+i] = bamgmesh_new->Vertices[3*id+1];

            gridX[tr] += gridCornerX[3*tr+i];
            gridY[tr] += gridCornerY[3*tr+i];
        }
        gridX[tr] /= 3;
        gridY[tr] /= 3;
    }

    // Initialise gridP, triangles, and weights
    std::vector<int> gridP(grid_size);
    std::vector<std::vector<int>> triangles(grid_size);
    std::vector<std::vector<double>> weights(grid_size);

    // Find which element of the _old mesh each P-point of the _new mesh hits (no need for defaults)
    /* TODO: Virtually all of the time is spent in InterpFromMeshToMesh2dx. We
     * could rewrite the current function so that we use the PreviousNumbering
     * information to find out which elements are new and only call
     * InterpFromMeshToMesh2dx for those. That should speed things up quite a
     * lot (probably an order of magnitude). But as the current function is
     * already fast enough that remains is an optimisation for another day.
     */
    double* elnum_out;
    InterpFromMeshToMesh2dx(&elnum_out,
                &indexTr[0],&coordX[0],&coordY[0],
                numNodes,numTriangles,
                &elnum[0],
                numTriangles,1,
                &gridX[0],&gridY[0],grid_size,
                true, -1);

    // Calculate weights
    for (int ppoint=0; ppoint<grid_size; ++ppoint)
    {
        assert( elnum_out[ppoint] >= 0. );

        // Carefully take the right integer value for element number
        int i_elnum_out = std::round(elnum_out[ppoint]);

        // Save the ppoint - for compatibility with ConservativeRemappingMeshToGrid
        gridP[ppoint] = ppoint;

        // initialise local vectors for the numbers for the elements that contribute and the weight of the contribution
        std::vector<int> local_triangles;
        std::vector<double> local_weights;

        // Check if the new and the old elements are the same
        // Take care to use the previous numbering
        std::vector<int> nodes_new(3);
        for ( int i=0; i<3; ++i )
        {
            int id = bamgmesh_new->Triangles[4*ppoint+i] - 1; // Here we need C/C++ numbering
/*  The previous numbering of the nodes on the boundaries given by bamg is wrong.
 *  Bamg adds these nodes at the beginning of the list of nodes.
 *  In our case, these nodes are always the same as the boundary is neither
 *  moving nor adapted. */
            if ( id > bamgmesh_new->VerticesOnGeomVertexSize[0] )
                nodes_new[i] = bamgmesh_new->PreviousNumbering[id] - 1; // Here we need C/C++ numbering
            else
                nodes_new[i] = id;
        }

        std::vector<int> nodes_old(3);
        for ( int i=0; i<3; ++i )
            nodes_old[i] = bamgmesh_old->Triangles[4*i_elnum_out+i] - 1; // Here we need C/C++ numbering

        // If the old and the new are the same we just save this information. Otherwise we must do a checkTriangle
        std::sort(nodes_old.begin(), nodes_old.end());
        std::sort(nodes_new.begin(), nodes_new.end());
        if ( nodes_old == nodes_new )
        {
            std::vector<std::pair<double,double>> points(3);
            for (int corner=0; corner<3; ++corner)
                points[corner] = std::make_pair(gridCornerX[3*ppoint+corner], gridCornerY[3*ppoint+corner]);

            local_triangles.push_back(i_elnum_out);
            local_weights.push_back(area(points));
            triangles[ppoint] = local_triangles;
            weights[ppoint] = local_weights;
        }
        else
        {
            // vertices of the old mesh in a format checkTriangle understands
            std::vector<double> cornerX(3);
            std::vector<double> cornerY(3);
            for (int corner=0; corner<3; ++corner)
            {
                cornerX[corner] = gridCornerX[3*ppoint+corner];
                cornerY[corner] = gridCornerY[3*ppoint+corner];
            }

            // Call the recursive function (this is our work horse here)
            checkTriangle(bamgmesh_old, cornerX, cornerY, i_elnum_out, local_triangles, local_weights);

            // Save the weights and triangle numbers
            triangles[ppoint] = local_triangles;
            weights[ppoint] = local_weights;
        }
    }
    xDelete<double>(elnum_out);

    // Now we apply the weights with meshToGrid - specify num_corners as 3
    ConservativeRemappingMeshToGrid(interp_out, interp_in, nb_var, grid_size, std::nan(""),
        gridP, gridCornerX, gridCornerY, triangles, weights, 3);
}

// Recursive function to check the current triangle
inline void checkTriangle(BamgMesh* bamgmesh, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY, int current_triangle, // inputs
		std::vector<int> &triangles, std::vector<double> &weights)  // outputs
{
    /*
     * We must plant the flag here, before any recursion can happen. We must
     * also push_back on weights and remember where to write the final weight
     * because some of the recursive calles may do a push_back of their own.
     */
    triangles.push_back(current_triangle);
    weights.push_back(0.);
    int loc = weights.size()-1;

    // This is a list of the points we use to calculate the area
    std::vector<std::pair<double,double>> points;

    int num_corners = gridCornerX.size();
    assert(num_corners = gridCornerY.size());

    // ---------- Now we do the three checks ---------- //

    // 1: check which of the nodes are inside the cell
    std::vector<bool> inCell(3); // Is the node inside the cell
    std::vector<int> nodeID(3);  // ID of the node
    std::vector<double> X(3);    // Node coordinates
    std::vector<double> Y(3);
    for (int i=0; i<3; ++i)
    {
        // ID and coordinates of the node - remember them for later
        nodeID[i] = bamgmesh->Triangles[4*current_triangle+i] - 1; // Here we need C/C++ numbering
        X[i] = bamgmesh->Vertices[3*nodeID[i]];
        Y[i] = bamgmesh->Vertices[3*nodeID[i]+1];

        // If we're inside we note the point and call self for the surrounding triangles
        inCell[i] = checkIfInside(gridCornerX, gridCornerY, X[i], Y[i]);
        if ( inCell[i] )
        {
            points.push_back(std::make_pair(X[i],Y[i]));

            int num_elements = bamgmesh->NodalElementConnectivitySize[1];
            for (int j=0; j<num_elements; j++)
            {
                int elt_num = bamgmesh->NodalElementConnectivity[num_elements*nodeID[i]+j] - 1; // Here we need C/C++ numbering
                // Negative elt_num means there are no more elements belonging to this node
                if ( elt_num < 0 ) continue;

                if ( ! visited(elt_num, triangles) )
                    checkTriangle(bamgmesh, gridCornerX, gridCornerY, elt_num, triangles, weights);
            }
        }
    }

    // If all the nodes are inside the cell we just calculate the triangle's area and return
    if ( inCell[0] && inCell[1] && inCell[2] )
    {
        weights[loc] = area(points);
        return;
    }

    // 2: Check and record which of the grid cell corners are inside the triangle
    int counter = 0;
    for (int i=0; i<num_corners; ++i)
    {
        if ( checkIfInside(X, Y, gridCornerX[i], gridCornerY[i]) )
        {
            points.push_back(std::make_pair(gridCornerX[i],gridCornerY[i]));
            ++counter;
        }
    }

    // If all the the cell is within the triangle we just calculate the cell's area and return
    if ( counter == num_corners )
    {
        weights[loc] = area(points);
        return;
    }

    // 3: Look for intersections between the triangle and cell
    int prev = 2;
    for (int i=0; i<3; prev=i++)
    {
        // If both the current and previous node are in the cell then we can skip this round
        if ( inCell[i] && inCell[prev] )
            continue;

        // Check for intersection, note the point and call self for the adjacent triangle
        if ( checkIfIntersecting(X[i], Y[i], X[prev], Y[prev], gridCornerX, gridCornerY, points) )
        {
            // We don't know which of the (up to) three connected triangles is the right one so we must search
            for (int j=0; j<3; j++)
            {
                int elt_num = bamgmesh->ElementConnectivity[3*current_triangle+j] - 1; // Here we need C/C++ numbering
                // Negative elt_num means there are no more elements belonging to this node
                if ( elt_num < 0 ) continue;

                /*
                 * Find the triangle adjacent to the current one by comparing
                 * the node IDs of the current edge with the node IDs of the
                 * selected triangle
                 */
                int counter = 0;
                for (int k=0; k<3; ++k)
                {
                    int myID = bamgmesh->Triangles[4*elt_num+k] - 1; // Here we need C/C++ numbering
                    if ( myID==nodeID[i] || myID==nodeID[prev] )
                        ++counter;
                }

                // Sanity check - the two triangles must share at least one node and no more than two!
                assert(counter>=1);
                assert(counter<=2);

                if ( counter==2 && !visited(elt_num, triangles) )
                    checkTriangle(bamgmesh, gridCornerX, gridCornerY, elt_num, triangles, weights);
            }
        }
    }

    // We've checked everything. Finish by calculating the area
    weights[loc] = area(points);

}

// Check if we've already visited this triangle
inline bool visited(int current_triangle, std::vector<int> const &triangles)
{
    for (auto it=triangles.begin(); it!=triangles.end(); ++it)
        if ( current_triangle == *it )
            return true;

    return false;
}

/*
 * Check if points are inside polygon
 * Short and works for both triangles and quadrangles (or any higher order polygon).
 *
 * from https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
 *
 * Copyright (c) 1970-2003, Wm. Randolph Franklin
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 *  1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
 *  2. Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
 *  3. The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
inline bool checkIfInside(std::vector<double> const &vertx, std::vector<double> const &verty, double testx, double testy)
{
    // Initilisation and sanity check
    bool inside = false;
    int nvert = vertx.size();
    assert(nvert==verty.size());

    int i;
    int j=nvert-1;
    for (i=0; i<nvert; j=i++)
        if ( ((verty[i]>testy) != (verty[j]>testy)) &&
            (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
                inside = !inside;

    return inside;
}

// Check for intersection and add intersecting points to the list
// https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
inline bool checkIfIntersecting(double X, double Y, double Xprev, double Yprev, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY, // inputs
		std::vector<std::pair<double,double>> &points) // side-effect
{
    // Initialise
    int num_corners = gridCornerX.size();
    assert(num_corners = gridCornerY.size());

    bool ret_val = false;
    double s1_x = X - Xprev;
    double s1_y = Y - Yprev;

    // Loop over the grid
    int prev=num_corners-1;
    for (int i=0; i<num_corners; prev=i++)
    {
        double s2_x = gridCornerX[i] - gridCornerX[prev];
        double s2_y = gridCornerY[i] - gridCornerY[prev];

        double det = -s2_x * s1_y + s1_x * s2_y;
        if ( det == 0 )
            continue; // The lines are parallel!

        double rdet = 1./det;
        double s = (-s1_y * (Xprev - gridCornerX[prev]) + s1_x * (Yprev - gridCornerY[prev])) * rdet;
        double t = ( s2_x * (Yprev - gridCornerY[prev]) - s2_y * (Xprev - gridCornerX[prev])) * rdet;

        /*
         * Here we assume that the case of overlaping points is an
         * intersection. It will result in double counting in some cases, but
         * not doing it would result in us not catching all the points all the
         * time.
         */
        if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
        {
            // Intersection detected
            points.push_back(std::make_pair(Xprev + (t * s1_x), Yprev + (t * s1_y)));
            ret_val = true;
        }
    }
    return ret_val;
}

// Calculate the area
inline double area(std::vector<std::pair<double,double>> &points)
{
    // Just a quick check
    if ( points.size() < 3 )
        return 0.;

    // First sort the points
    sortClockwise(points);

    // Calculate value of shoelace formula
    double area = 0.;
    int n = points.size();
    int j = n-1;
    for (int i=0; i < n; j=i++)
        area += (points[j].first + points[i].first) * (points[j].second - points[i].second);

    // Return half the absolute value
    return std::abs(area) * 0.5;
}

// Sort the points in a clock wise manner around the centre of the cloud of points
inline void sortClockwise(std::vector<std::pair<double,double>> &points)
{
    // Calculate the centre point
    int n = points.size();
    double centreX = 0.;
    double centreY = 0.;
    for (auto it=points.begin(); it!=points.end(); ++it)
    {
        centreX += it->first;
        centreY += it->second;
    }
    double rn = 1./(double)n;
    centreX *= rn;
    centreY *= rn;

    // Subtract the centre
    for (auto it=points.begin(); it!=points.end(); ++it)
    {
        it->first  -= centreX;
        it->second -= centreY;
    }

    // Sort with a custom function
    std::sort(points.begin(), points.end(), CWSort);

    // Erase duplicates - these may occur when we double count nodes on grid corners (i.e. costal points)
    //points.erase( std::unique( points.begin(), points.end() ), points.end() );

    // Add the centre again
    for (auto it=points.begin(); it!=points.end(); ++it)
    {
        it->first  += centreX;
        it->second += centreY;
    }
}

// The custom sorting function
// https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order <- the second code exmple
inline bool CWSort(std::pair<double,double> p1, std::pair<double,double> p2)
{
    // Computes the quadrant for p1 and p2 (0-3):
    //     ^
    //   1 | 0
    //  ---+-->
    //   2 | 3

    const int dax = (p1.first > 0) ? 1 : 0;
    const int day = (p1.second > 0) ? 1 : 0;
    const int qa = (1 - dax) + (1 - day) + ((dax & (1 - day)) << 1);

    /* The previous computes the following:

       const int qa =
       (  p1.first
        ? ( p1.second
            ? 0 : 3)
        : (p1.second)
            ? 1 : 2)); */

    const int dbx = (p2.first > 0) ? 1 : 0;
    const int dby = (p2.second > 0) ? 1 : 0;
    const int qb = (1 - dbx) + (1 - dby) + ((dbx & (1 - dby)) << 1);

    if (qa == qb)
        return p2.first*p1.second < p2.second*p1.first;
    else
        return qa < qb;
}
