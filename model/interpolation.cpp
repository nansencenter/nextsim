/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   interpolation.cpp
 * @author Fabien Salmon <fabien.salmon@inria.fr>
 * @date   Thu Aug  4 09:47:27 CEST 2024
 */

#include "interpolation.hpp"
#include <vector>

/**
 * @brief Different functions allowing interpolations between two meshes
 * The functions come from bamg and are adapted to neXtSIM
 *
 * @see
 *
 */

namespace Nextsim
{

void
ConservativeRemappingWithWeights(GmshMeshSeq const& mesh, std::vector<double> &gridX, std::vector<double> &gridY, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY,
                                 std::vector<int> &gridP, std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights)
{
    // ---------- Initialisation ---------- //
    // Copy the triangle information
    int numTriangles = mesh.numTriangles();
    std::vector<int> indexTr = mesh.indexTr();
    std::vector<double> elnum(numTriangles);
    for (int tr = 0; tr < numTriangles; ++tr)
    {
        elnum[tr] = tr;
    }

    // Copy the node information
    int numNodes     = mesh.numNodes();
    std::vector<double> coordX = mesh.coordX();
    std::vector<double> coordY = mesh.coordY();

    // Copy the list of triangles containing each point
    std::vector<std::vector<int>> list_triangles(mesh.numNodes());
    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        for (int k = 0; k <= 2; k++)
        {
            list_triangles[mesh.triangles()[n].indices[k]-1].push_back(n);
        }
    }

    // List of triangles associated with each edge
    bool exit_loop;
    std::vector<std::vector<int>> list_triangles_edges(3*mesh.numTriangles());
    int n_nodes1, n_nodes2;
    for (int i = 0; i < mesh.numTriangles(); i++)
    {
        for (int l = 0; l < 3; l++)
        {
            n_nodes1 = mesh.triangles()[i].indices[l]-1;
            n_nodes2 = mesh.triangles()[i].indices[(l+1)%3]-1;
            exit_loop = false;

            for (int j = 0; j < list_triangles[n_nodes1].size() && !exit_loop; j++)
            {
                for (int k = 0; k < list_triangles[n_nodes2].size(); k++)
                {
                    if (list_triangles[n_nodes1][j] == list_triangles[n_nodes2][k] && list_triangles[n_nodes1][j] != i)
                    {
                        list_triangles_edges[3*i+l].push_back(i);
                        list_triangles_edges[3*i+l].push_back(list_triangles[n_nodes1][j]);
                        exit_loop = true;
                        break;
                    }
                }
            }
        }
    }

    // Copy the list of element neighbours of each element
    int n_elt1, n_elt2, index;
    std::vector<std::vector<int>> list_neighbours(mesh.numTriangles(),std::vector<int>(3,-1));
    for (int i = 0; i < list_triangles_edges.size(); i++)
    {
        // Check whether the edge is on a boundary
        if (list_triangles_edges[i].empty()) {continue;}

        n_elt1 = list_triangles_edges[i][0];
        n_elt2 = list_triangles_edges[i][1];

        // Check whether the edge has already been treated (occurs once for each edge)
        if (std::find(list_neighbours[n_elt1].begin(), list_neighbours[n_elt1].end(), n_elt2) != list_neighbours[n_elt1].end()) {continue;}

        index = std::distance(list_neighbours[n_elt1].begin(), std::min_element(list_neighbours[n_elt1].begin(), list_neighbours[n_elt1].end()));
        list_neighbours[n_elt1][index] = n_elt2;
        index = std::distance(list_neighbours[n_elt2].begin(), std::min_element(list_neighbours[n_elt2].begin(), list_neighbours[n_elt2].end()));
        list_neighbours[n_elt2][index] = n_elt1;
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

            // Call the recursive function (this is our work horse here)
            checkTriangle(indexTr, coordX, coordY, cornerX, cornerY, i_elnum_out, list_triangles, list_neighbours, local_triangles, local_weights);

            // Save the weights and triangle numbers
            triangles.push_back(local_triangles);
            weights.push_back(local_weights);
        }
    }

    xDelete<double>(elnum_out);

}//ConservativeRemappingWithWeights

// Remapping from mesh to mesh.
// In this case we want to both calculate weights and apply them in the same step
// Drop-in-replacement for InterpFromMeshToMesh2dCavities
void ConservativeRemappingFromMeshToMesh(double* &interp_out, std::vector<double> &interp_in, int nb_var,
                                         std::vector<int> &indexTr, std::vector<double> &coordX, std::vector<double> &coordY,
                                         std::vector<int> &new_indexTr, std::vector<double> &new_coordX, std::vector<double> &new_coordY)
{
    // We start off the same as ConservativeRemappingWeights - only here the new mesh replaces the grid

    // ---------- Initialisation ---------- //
    // Copy the triangle information of the _old mesh
    int numTriangles = indexTr.size()/3;
    std::vector<double> elnum(numTriangles);
    for (int tr = 0; tr < numTriangles; ++tr)
    {
        elnum[tr] = tr;
    }

    int numNodes     = coordX.size();

    // Copy the list of triangles containing each point
    std::vector<std::vector<int>> list_triangles(numNodes);
    for (int n = 0; n < numTriangles; n++)
    {
        for (int k = 0; k <= 2; k++)
        {
            list_triangles[indexTr[3*n+k]-1].push_back(n);
        }
    }

    // List of triangles associated with each edge
    bool exit_loop;
    std::vector<std::vector<int>> list_triangles_edges(3*numTriangles);
    int n_nodes1, n_nodes2;
    for (int i = 0; i < numTriangles; i++)
    {
        for (int l = 0; l < 3; l++)
        {
            n_nodes1 = indexTr[3*i+l]-1;
            n_nodes2 = indexTr[3*i+(l+1)%3]-1;
            exit_loop = false;

            for (int j = 0; j < list_triangles[n_nodes1].size() && !exit_loop; j++)
            {
                for (int k = 0; k < list_triangles[n_nodes2].size(); k++)
                {
                    if (list_triangles[n_nodes1][j] == list_triangles[n_nodes2][k] && list_triangles[n_nodes1][j] != i)
                    {
                        list_triangles_edges[3*i+l].push_back(i);
                        list_triangles_edges[3*i+l].push_back(list_triangles[n_nodes1][j]);
                        exit_loop = true;
                        break;
                    }
                }
            }
        }
    }

    // Copy the list of element neighbours of each element
    int n_elt1, n_elt2;
    std::vector<int> index(numTriangles,0);
    std::vector<std::vector<int>> list_neighbours(numTriangles,std::vector<int>(3,-1));
    for (int i = 0; i < list_triangles_edges.size(); i++)
    {
        // Check whether the edge is on a boundary
        if (list_triangles_edges[i].empty()) {continue;}

        n_elt1 = list_triangles_edges[i][0];
        n_elt2 = list_triangles_edges[i][1];

        // Check whether the edge has already been treated (occurs once for each edge)
        if (std::find(list_neighbours[n_elt1].begin(), list_neighbours[n_elt1].end(), n_elt2) != list_neighbours[n_elt1].end()) {continue;}

        list_neighbours[n_elt1][index[n_elt1]++] = n_elt2;
        list_neighbours[n_elt2][index[n_elt2]++] = n_elt1;
    }

    // Copy the triangle information of the _new mesh and calculate the barycentre
    // Keep the nomenclature for a grid (even if it's a mesh)
    int grid_size = new_indexTr.size()/3;
    std::vector<double> gridX(grid_size);
    std::vector<double> gridY(grid_size);
    std::vector<double> gridCornerX(3*grid_size);
    std::vector<double> gridCornerY(3*grid_size);
    for (int tr = 0; tr < grid_size; ++tr)
    {
        gridX[tr] = 0.; // barycentre
        gridY[tr] = 0.;
        for (int i=0; i<3; ++i)
        {
            // grid corner == vertice
            gridCornerX[3*tr+i] = new_coordX[new_indexTr[3*tr+i]-1];
            gridCornerY[3*tr+i] = new_coordY[new_indexTr[3*tr+i]-1];

            gridX[tr] += gridCornerX[3*tr+i];
            gridY[tr] += gridCornerY[3*tr+i];
        }
        gridX[tr] /= 3.;
        gridY[tr] /= 3.;
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

        // vertices of the old mesh in a format checkTriangle understands
        std::vector<double> cornerX(3);
        std::vector<double> cornerY(3);
        for (int corner=0; corner<3; ++corner)
        {
            cornerX[corner] = gridCornerX[3*ppoint+corner];
            cornerY[corner] = gridCornerY[3*ppoint+corner];
        }

        // Call the recursive function (this is our work horse here)
        checkTriangle(indexTr, coordX, coordY, cornerX, cornerY, i_elnum_out, list_triangles, list_neighbours, local_triangles, local_weights);

        // Save the weights and triangle numbers
        triangles[ppoint] = local_triangles;
        weights[ppoint] = local_weights;
    }
    xDelete<double>(elnum_out);

    // Now we apply the weights with meshToGrid - specify num_corners as 3
    ConservativeRemappingMeshToGrid(interp_out, interp_in, nb_var, grid_size, std::nan(""),
                                    gridP, gridCornerX, gridCornerY, triangles, weights, 3);

}//ConservativeRemappingFromMeshToMesh

// Apply weights for a mesh-to-grid remapping
void 
ConservativeRemappingMeshToGrid(double* &interp_out, std::vector<double> &interp_in, int const nb_var, int grid_size, double miss_val,
                                std::vector<int> &gridP, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY,
                                std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &weights, int num_corners)
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
            /*if (var == 2) {
                interp_out[gridP[i]*nb_var+var] = interp_in[triangles[i][0]*nb_var+var];
                continue;
            }*/
            interp_out[gridP[i]*nb_var+var] = 0.;

            for (int tr=0; tr<triangles[i].size(); ++tr) 
                interp_out[gridP[i]*nb_var+var] += interp_in[triangles[i][tr]*nb_var+var]*weights[i][tr];

            interp_out[gridP[i]*nb_var+var] *= r_cell_area;
        }
    }
} // ConservativeRemappingMeshToGrid

// Recursive function to check the current triangle
void 
checkTriangle(std::vector<int> &indexTr, std::vector<double> &coordX, std::vector<double> &coordY, //inputs
              std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY, int current_triangle, //inputs
              std::vector<std::vector<int>> const& list_triangles, std::vector<std::vector<int>> const& list_neighbours, // inputs
              std::vector<int> &triangles, std::vector<double> &weights)  // outputs
{
    /*
     * We must plant the flag here, before any recursion can happen. We must
     * also push_back on weights and remember where to write the final weight
     * because some of the recursive calls may do a push_back of their own.
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
        nodeID[i] = indexTr[3*current_triangle+i] - 1; // mesh.triangles()[current_triangle].indices[i] - 1;
        X[i] = coordX[nodeID[i]]; //mesh.nodes()[nodeID[i]].coords[0];
        Y[i] = coordY[nodeID[i]]; //mesh.nodes()[nodeID[i]].coords[1];

        // If we're inside we note the point and call self for the surrounding triangles
        inCell[i] = checkIfInside(gridCornerX, gridCornerY, X[i], Y[i]);
        if ( inCell[i] )
        {
            points.push_back(std::make_pair(X[i],Y[i]));

            for (int j = 0; j < list_triangles[nodeID[i]].size(); j++)
            {
                int elt_num = list_triangles[nodeID[i]][j];

                if ( ! visited(elt_num, triangles) )
                    checkTriangle(indexTr, coordX, coordY, gridCornerX, gridCornerY, elt_num, list_triangles, list_neighbours, triangles, weights);
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
                int elt_num = list_neighbours[current_triangle][j]; 

                // Negative elt_num means there are no more elements belonging to this node
                if ( elt_num < 0 ) break;

                /*
                 * Find the triangle adjacent to the current one by comparing
                 * the node IDs of the current edge with the node IDs of the
                 * selected triangle
                 */
                int counter = 0;
                for (int k=0; k<3; ++k)
                {
                    int myID = indexTr[3*elt_num+k]-1; //mesh.triangles()[elt_num].indices[k] - 1;
                    if ( myID==nodeID[i] || myID==nodeID[prev] )
                        ++counter;
                }

                // Sanity check - the two triangles must share at least one node and no more than two!
                assert(counter>=1);
                assert(counter<=2);

                if ( counter == 2 && !visited(elt_num, triangles) )
                    checkTriangle(indexTr, coordX, coordY, gridCornerX, gridCornerY, elt_num, list_triangles, list_neighbours, triangles, weights);
            }
        }
    }

    // We've checked everything. Finish by calculating the area
    weights[loc] = area(points);

}//checkTriangle

// Check if we've already visited this triangle
bool 
visited(int current_triangle, std::vector<int> const &triangles)
{
    for (auto it = triangles.begin(); it != triangles.end(); ++it)
        if ( current_triangle == *it )
            return true;

    return false;
} // visited

// Check if points are inside polygon
bool checkIfInside(std::vector<double> const &vertx, std::vector<double> const &verty, double testx, double testy)
{
    // Initilisation and sanity check
    bool inside = false;
    double EPSILON = 1.e-10;
    int nvert = vertx.size();
    assert(nvert==verty.size());

    double dx = std::max(fabs(vertx[1] - vertx[0]), fabs(vertx[2] - vertx[0]));

    // Check if the point is on a vertex of the triangle
    for (int i = 0; i < nvert; ++i) {
        if (fabs(testx - vertx[i]) < EPSILON*dx && fabs(testy - verty[i]) < EPSILON*dx) {
            return true;
        }
    }

    // Check if the point is on an edge of the triangle
    for (int i = 0, j = nvert - 1; i < nvert; j = i++) {
        if (fabs((testx - vertx[i]) * (verty[j] - verty[i]) - (testy - verty[i]) * (vertx[j] - vertx[i])) < EPSILON*dx*dx &&
            (testx >= std::min(vertx[i], vertx[j]) && testx <= std::max(vertx[i], vertx[j])) &&
            (testy >= std::min(verty[i], verty[j]) && testy <= std::max(verty[i], verty[j]))) {
            return true;
        }
    }

    // Ray-casting
    int j = nvert-1;
    for (int i = 0; i<nvert; j=i++)
        if ( ((verty[i]>testy) != (verty[j]>testy)) &&
            (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
                inside = !inside;

    return inside;
} //checkIfInside

// Check for intersection and add intersecting points to the list
bool 
checkIfIntersecting(double X, double Y, double Xprev, double Yprev, std::vector<double> const &gridCornerX, std::vector<double> const &gridCornerY, // inputs
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
} // checkIfIntersecting

// Calculate the area
double 
area(std::vector<std::pair<double,double>> &points)
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
} // area

// Sort the points in a clock wise manner around the centre of the cloud of points
void 
sortClockwise(std::vector<std::pair<double,double>> &points)
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
} // sortClockwise

// The custom sorting function
bool CWSort(std::pair<double,double> p1, std::pair<double,double> p2)
{
    // Computes the quadrant for p1 and p2 (0-3):
    //     ^
    //   1 | 0
    //  ---+-->
    //   2 | 3

    const int dax = (p1.first > 0) ? 1 : 0;
    const int day = (p1.second > 0) ? 1 : 0;
    const int qa = (1 - dax) + (1 - day) + ((dax & (1 - day)) << 1);

    const int dbx = (p2.first > 0) ? 1 : 0;
    const int dby = (p2.second > 0) ? 1 : 0;
    const int qb = (1 - dbx) + (1 - dby) + ((dbx & (1 - dby)) << 1);

    if (qa == qb)
        return p2.first*p1.second < p2.second*p1.first;
    else
        return qa < qb;
} // CWSort

}// nextsim
