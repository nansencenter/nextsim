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
 * The functions come from bamg
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

    delete elnum_out;

}//ConservativeRemappingWithWeights

// Remapping from mesh to mesh.
// In this case we want to both calculate weights and apply them in the same step
// Drop-in-replacement for InterpFromMeshToMesh2dCavities
void ConservativeRemappingFromMeshToMesh(double* &interp_out, std::vector<double> &interp_in, int nb_var,
//                                         GmshMeshSeq const& mesh_old, GmshMeshSeq const& mesh_new)
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
            gridCornerX[3*tr+i] = new_coordX[new_indexTr[3*tr+i]-1];// mesh_new.nodes()[mesh_new.triangles()[tr].indices[i]-1].coords[0];
            gridCornerY[3*tr+i] = new_coordY[new_indexTr[3*tr+i]-1];// mesh_new.nodes()[mesh_new.triangles()[tr].indices[i]-1].coords[1];

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
        assert( elnum_out[ppoint] >= -0.1 );

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
    delete elnum_out;

    // Now we apply the weights with meshToGrid - specify num_corners as 3
    ConservativeRemappingMeshToGrid(interp_out, interp_in, nb_var, grid_size, std::nan(""),
                                    gridP, gridCornerX, gridCornerY, triangles, weights, 3);

}//ConservativeRemappingMeshToMesh

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
int nb = 0;
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
            if (var == 2) {
                interp_out[gridP[i]*nb_var+var] = interp_in[triangles[i][0]*nb_var+var];
                continue;
            }
            interp_out[gridP[i]*nb_var+var] = 0.;

double test = 0.;
for (int tr=0; tr<triangles[i].size(); ++tr) test += weights[i][tr];
if (var == 0 && fabs(test*r_cell_area-1) > 1.e-8) nb++;

            for (int tr=0; tr<triangles[i].size(); ++tr) 
                interp_out[gridP[i]*nb_var+var] += interp_in[triangles[i][tr]*nb_var+var]*weights[i][tr];

            interp_out[gridP[i]*nb_var+var] *= r_cell_area;
        }
    }
std::cout << "nobre faux " << nb << std::endl;
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

double compute_area(double a0, double a1, double b0, double b1, double c0, double c1)
{return 0.5*((b0 - a0)*(c1 - a1) - (b1 - a1)*(c0 - a0));}

// Find the closest triangle or the containing triangle
std::vector<double> interpolationNode(int* indexTria, double* coordX, double* coordY, int numTriangles, int numNodes, int M_nb_var_node, double x, double y, double* data_in)
{
    int k,i;
    int isinside = 0;
    std::vector<double> contrib(3);
    double closestDist = 1.e20;
    int closestTria;

    // Loop over the triangles
    for (k = 0; k < numTriangles; k++)
    {
        // Compute the area of the triangle
        double area = compute_area( coordX[indexTria[3*k]-1], coordY[indexTria[3*k]-1], coordX[indexTria[3*k+1]-1], coordY[indexTria[3*k+1]-1], coordX[indexTria[3*k+2]-1], coordY[indexTria[3*k+2]-1]);
        for (i = 0; i < 3; i++)
        {
            contrib[i] = compute_area( x, y, coordX[indexTria[3*k+(i+1)%3]-1], coordY[indexTria[3*k+(i+1)%3]-1], coordX[indexTria[3*k+(i+2)%3]-1], coordY[indexTria[3*k+(i+2)%3]-1] ) / area;
        }

        // Check whether the point is inside the triangle
        if ( std::min(contrib[0], std::min(contrib[1], contrib[2])) > -1.e-12 ) isinside = 1;

        // Distance from center of mass
        double dist[2];
        dist[0] = x;
        dist[1] = y;
        
        for (i = 0; i < 3; i++)
        {
            dist[0] -= coordX[indexTria[3*k+i]-1] /3.;
            dist[1] -= coordY[indexTria[3*k+i]-1] /3.;
        }

        double norm = sqrt(dist[0]*dist[0] + dist[1]*dist[1]);

        if (norm < closestDist)
        {
            closestDist = norm;
            closestTria = k;
        }
        
        if (isinside) break;
    }

    std::vector<double> data_out(M_nb_var_node,0.);
    for (k = 0; k < M_nb_var_node; k++)
    {
        for (i = 0; i < 3; i++)
        {
            data_out[k] += contrib[i] * data_in[indexTria[3*closestTria+i]-1 + k*numNodes];
        }
    }

    return data_out;
} // interpolationNode

void interpolationNodes(double** interp_out, int* indexTria_prv, double* coordX_prv, double* coordY_prv, int numTriangles_prv, 
                        int numNodes_prv, int M_nb_var_node, double* coordX, double* coordY, int numNodes, double* data_in)
{
    int k, i;
    double* data_interp=NULL;

    data_interp = xNew<double>(numNodes*M_nb_var_node);

    for (k = 0; k < numNodes; k++)
    {
        std::vector<double> data_out = interpolationNode(&indexTria_prv[0], &coordX_prv[0], &coordY_prv[0], numTriangles_prv, numNodes_prv, M_nb_var_node, coordX[k], coordY[k], &data_in[0]);
        for (i = 0; i < M_nb_var_node; i++) data_interp[k+i*numNodes] = data_out[i];
    }

    *interp_out = data_interp;
} //interpolationNodes

void find_bounding_box(double* coordX, double* coordY, int numNodes, double* xmin, double* xmax, double* ymin, double* ymax) 
{
    *xmax = -1.e20;
    *xmin = 1.e20;
    *ymax = -1.e20;
    *ymin = 1.e20;

    for (int k = 0; k < numNodes; k++) 
    {
        if (coordX[k] < *xmin) *xmin = coordX[k];
        if (coordX[k] > *xmax) *xmax = coordX[k];
        if (coordY[k] < *ymin) *ymin = coordY[k];
        if (coordY[k] > *ymax) *ymax = coordY[k];
    }
}

 // Divide the domain as a Cartesian grid and store the triangles in each cell of this grid.
std::vector<std::vector<std::vector<int>>> grid_size_triangles(double* coordX, double* coordY, int* indexTria, int numTriangles, double minX, double minY, double maxX, double maxY, int GRID_SIZE, double EPSILON) 
{
    int i, j, k;
    int idx, idy, maxx, maxy, minx, miny;
    int list_idx[6], list_idy[6];

    double cellWidth = (maxX - minX) / GRID_SIZE;
    double cellHeight = (maxY - minY) / GRID_SIZE;

    // Initialize list_triangles
    std::vector<std::vector<std::vector<int>>> list_triangles;
    list_triangles.resize(GRID_SIZE);

    for (i = 0; i < GRID_SIZE; i++)
    {
        list_triangles[i].resize(GRID_SIZE);
        for (j = 0; j < GRID_SIZE; j++)
        {
            list_triangles[i][j].push_back(0);
        }
    }

    for (k = 0; k < numTriangles; k++) {

        // List of cells containing at least a part of the triangle
        for (i = 0; i < 3; i++) {
          list_idx[i] = (int)((coordX[indexTria[3*k+i]-1] - minX) / cellWidth + EPSILON);
          list_idy[i] = (int)((coordY[indexTria[3*k+i]-1] - minY) / cellHeight + EPSILON);
          list_idx[3+i] = (int)((coordX[indexTria[3*k+i]-1] - minX) / cellWidth - EPSILON);
          list_idy[3+i] = (int)((coordY[indexTria[3*k+i]-1] - minY) / cellHeight - EPSILON);
        }

        minx = *std::min_element(list_idx, list_idx+6);
        if (minx < 0) minx = 0;
        maxx = *std::max_element(list_idx, list_idx+6);
        if (maxx >= GRID_SIZE) maxx = GRID_SIZE-1;
        miny = *std::min_element(list_idy, list_idy+6);
        if (miny < 0) miny = 0;
        maxy = *std::max_element(list_idy, list_idy+6);
        if (maxy >= GRID_SIZE) maxy = GRID_SIZE-1;

        for (idx = minx; idx <= maxx; idx++) 
        {
            for (idy = miny; idy <= maxy; idy++) 
            {
                list_triangles[idx][idy][0]++;
                list_triangles[idx][idy].push_back(k);
            }
        }
    }

    return list_triangles;
}

void interpolationElements(double** interp_out, int* indexTria_prv, double* coordX_prv, double* coordY_prv, int numTriangles_prv,
                           int numNodes_prv, int nb_var_element, double* coordX, double* coordY, int* indexTria, int numTriangles, double* data_in)
{
    int i, j, k, l, n, p;
    int GRID_SIZE = sqrt(numTriangles);
    int list_idx[6], list_idy[6];
    double EPSILON = 1.e-12;
    double area;

    double xmin, xmax, ymin, ymax, minx, maxx, miny, maxy;
    find_bounding_box(&coordX_prv[0], &coordY_prv[0], numNodes_prv, &xmin, &xmax, &ymin, &ymax);

    double cellWidth = (xmax - xmin) / GRID_SIZE;
    double cellHeight = (ymax - ymin) / GRID_SIZE;

    std::vector<std::vector<std::vector<int>>> list_triangles = grid_size_triangles(&coordX_prv[0], &coordY_prv[0], &indexTria_prv[0], numTriangles_prv, xmin, ymin, xmax, ymax, GRID_SIZE, EPSILON);
    double* data_interp=NULL;

    data_interp = xNew<double>(numTriangles*nb_var_element);
    for (k = 0; k < numTriangles*nb_var_element; k++) data_interp[k] = 0.;

//std::cout << "MERDE " << intersection_area(Trian1, Trian2) << std::endl;
 // One file per process
// std::ofstream vtkfile {"Triangle.vtk"};
//vtkfile << "# vtk DataFile Version 2.0\nmesh_metric\nASCII\nDATASET POLYDATA\nPOINTS 10 float\n";

    for (k = 0; k < numTriangles; k++)
    { 
//std::cout << (1.*k)/numTriangles <<std::endl;
        std::vector<std::vector<double>> Triangle(3, std::vector<double>(2));
        for (i = 0; i < 3; i++)
        {
            Triangle[i][0] = coordX[indexTria[3*k+i]-1];
            Triangle[i][1] = coordY[indexTria[3*k+i]-1];
//if (k == 1830) std::cout << Triangle[i][0] << " " << Triangle[i][1] << std::endl;
        }

        // Compute triangle area
        area = compute_area( coordX[indexTria[3*k]-1], coordY[indexTria[3*k]-1], coordX[indexTria[3*k+1]-1], coordY[indexTria[3*k+1]-1], coordX[indexTria[3*k+2]-1], coordY[indexTria[3*k+2]-1]);

        // List of cells containing at least a part of the triangle
        for (i = 0; i < 3; i++) {
            list_idx[i] = (int)((coordX[indexTria[3*k+i]-1] - xmin) / cellWidth + EPSILON);
            list_idy[i] = (int)((coordY[indexTria[3*k+i]-1] - ymin) / cellHeight + EPSILON);
            list_idx[3+i] = (int)((coordX[indexTria[3*k+i]-1] - xmin) / cellWidth - EPSILON);
            list_idy[3+i] = (int)((coordY[indexTria[3*k+i]-1] - ymin) / cellHeight - EPSILON);
        }

        minx = *std::min_element(list_idx, list_idx+6);
        if (minx < 0) minx = 0;
        maxx = *std::max_element(list_idx, list_idx+6);
        if (maxx >= GRID_SIZE) maxx = GRID_SIZE-1;
        miny = *std::min_element(list_idy, list_idy+6);
        if (miny < 0) miny = 0;
        maxy = *std::max_element(list_idy, list_idy+6);
        if (maxy >= GRID_SIZE) maxy = GRID_SIZE-1;

        // Loop over all the triangles that intersect triangle k potentially
        std::vector<int> triangle_color(numTriangles_prv,0);

        for (i = minx; i <= maxx; i++)
        {
            for (j = miny; j <= maxy; j++)
            {
                for (n = 1; n <= list_triangles[i][j][0]; n++)
                {
                    // Triangle index
                    p = list_triangles[i][j][n];

                    if (triangle_color[p]) continue; // Already treated

                    triangle_color[p] = 1;
                    std::vector<std::vector<double>> Triangle2(3, std::vector<double>(2));
                    for (l = 0; l < 3; l++)
                    {
                        Triangle2[l][0] = coordX_prv[indexTria_prv[3*p+l]-1];
                        Triangle2[l][1] = coordY_prv[indexTria_prv[3*p+l]-1];
//if (k == 1830) vtkfile << Triangle2[l][0] << " " << Triangle2[l][1] << " 0.\n";
//if (k == 1830) std::cout << Triangle2[l][0] << " " << Triangle2[l][1] << std::endl;
                    }
            
                    // Compute the intersection area
                    double area_intersection = intersection_area(Triangle, Triangle2);
                    
                    if (area_intersection < 1.e-12) continue;
//std::cout << k << " " << area << " " << area_intersection << " " << p << " " << Triangle2[0][0] << " " << Triangle2[0][1] << " " << Triangle2[1][0] << " " << Triangle2[1][1] << " " << Triangle2[2][0] << " " << Triangle2[2][1] << std::endl;
//std::cout << Triangle[0][0] << " " << Triangle[0][1] << " " << Triangle[1][0] << " " << Triangle[1][1] << " " << Triangle[2][0] << " " << Triangle[2][1] << std::endl;

//std::cout << k << " " << area << " " << area_intersection << std::endl;//" " << compute_area(Triangle2[0][0], Triangle2[0][1], Triangle2[1][0], Triangle2[1][1], Triangle2[2][0], Triangle2[2][1]) << std::endl;

                    // Add the contribution of this triangle
                    for (l = 0; l < nb_var_element; l++)
                    {
                        data_interp[k*nb_var_element + l] += area_intersection / area * data_in[p*nb_var_element + l];
                    }

                }
            }
        }
//if (k == 1830) {    vtkfile << "VERTICES " << k << " " << 2*k <<"\n";
//vtkfile.close();
//}
    }

    *interp_out = data_interp;
} //interpolationElements

// Intersection between two segments
std::vector<double> find_segment_intersection(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, int* is_intersected)
{
    std::vector<double> intersection(2);

    double determinant = (b[0]-a[0]) * (d[1]-c[1]) - (b[1] - a[1]) * (d[0] - c[0]);

    // Check whether the two segments are parallel
    if (fabs(determinant) < 1.e-12)
    {
        *is_intersected = 0; 
        return intersection;        
    }

    double t = ((c[0] - a[0]) * (d[1] - c[1]) - (c[1] - a[1]) * (d[0] - c[0])) / determinant;
    double u = ((c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0])) / determinant;
    
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1)
    {
        *is_intersected = 1;
    }
    else
    {
        *is_intersected = 0;
    }

    intersection[0] = a[0] + t * (b[0] - a[0]);
    intersection[1] = a[1] + t * (b[1] - a[1]);
//std::cout << *is_intersected << intersection[0] << " " << intersection[1] << " " << a[0] << " " << a[1] << " " << b[0] << " " << b[1] << " " << c[0] << " " << c[1] << " " << d[0] << " " << d[1] << " " << t <<std::endl;
    return intersection;
} // find_segment_intersection

// Cross product
double cross_product(std::vector<double> v1, std::vector<double> v2)
{
    return v1[0] * v2[1] - v1[1] * v2[0];
}

// Check whether a point is inside a triangle
int is_inside(std::vector<double> Point, std::vector<std::vector<double>> Triangle)
{

    // Check the sign of the cross products
    std::vector<std::vector<double>> vectors(6);
    for (int i = 0; i < 3; i++) 
    {
        vectors[2*i].resize(2);
        vectors[2*i][0] = Triangle[(i+1)%3][0] - Triangle[i][0];
        vectors[2*i][1] = Triangle[(i+1)%3][1] - Triangle[i][1];
        vectors[2*i+1].resize(2);
        vectors[2*i+1][0] = Point[0] - Triangle[i][0];
        vectors[2*i+1][1] = Point[1] - Triangle[i][1];
    }

    double cross1 = cross_product(vectors[0], vectors[1]);
    double cross2 = cross_product(vectors[2], vectors[3]);
    double cross3 = cross_product(vectors[4], vectors[5]);

    if ((cross1 >=0 && cross2 >= 0 && cross3 >=0) || (cross1 <=0 && cross2 <= 0 && cross3 <=0))
    {
        return 1;
    }
    return 0;

} // is_inside

std::vector<std::vector<double>> sort_points(std::vector<std::vector<double>> list_intersection_points)
{
    std::vector<double> center(2);
    int k, i;
    int point_number = list_intersection_points.size();
    
    // Compute the barycenter of the points
    for (k = 0; k < point_number; k++)
    {
        center[0] += list_intersection_points[k][0];
        center[1] += list_intersection_points[k][1];
    }

    center[0] /= point_number;
    center[1] /= point_number;

    // Compute the angle with the x-axis
    std::vector<double> angle(point_number);
    for (k = 0; k < point_number; k++)
    {
        angle[k] = atan2(list_intersection_points[k][1] - center[1], list_intersection_points[k][0] - center[0]);
    }

    std::vector<std::vector<double>> list_sorted_points(point_number);
    for (k = 0; k < point_number; k++)
    {
        double minimum = 10.;
        int i_min;
        for (i = 0; i < point_number; i++)
        {
            if (angle[i] < minimum)
            {
                minimum = angle[i];
                i_min = i;
            }
        }
        angle[i_min] = 11.;
        list_sorted_points[k] = list_intersection_points[i_min];
    }

    return list_sorted_points;

} // sort_points

int aligned(std::vector<std::vector<double>> list_intersection_points)
{

    for (int k = 1; k < list_intersection_points.size()-1; k++)
    {
        std::vector<double> v1(2), v2(2);
        for (int i = 0; i < 2; i++)
        {
            v1[i] = list_intersection_points[k-1][i] - list_intersection_points[k][i];
            v2[i] = list_intersection_points[k+1][i] - list_intersection_points[k][i];
        }
        if (fabs(cross_product(v1, v2)) > 1.e-8) return 0;
    }

    return 1;
}

// Compute the intersection area between two triangles
double intersection_area(std::vector<std::vector<double>> Tria1, std::vector<std::vector<double>> Tria2)
{
    int k1, k2;
    int is_intersected;
    std::vector<double> intersection(2);
    std::vector<std::vector<double>> list_intersection_points;

    // Find the intersection points between the triangles
    for (k1 = 0; k1 < 3; k1++)
    {
        for (k2 = 0; k2 < 3; k2++)
        {
            intersection = find_segment_intersection(Tria1[k1], Tria1[(k1+1)%3], Tria2[k2], Tria2[(k2+1)%3], &is_intersected);
            if (is_intersected) {
list_intersection_points.push_back(intersection);
//if (k == 1830) std::cout <<"GGG " << intersection[0] << " " << intersection[1] << std::endl;
}
        }
    }

    // Find the points of the triangles that are inside the other one
    for (k1 = 0; k1 < 3; k1++)
    {
        if (is_inside(Tria1[k1], Tria2)) list_intersection_points.push_back(Tria1[k1]);
        if (is_inside(Tria2[k1], Tria1)) list_intersection_points.push_back(Tria2[k1]);
    }

    int number_points = list_intersection_points.size();
//std::cout <<"GGG "<< number_points << std::endl;
    if (number_points < 3) return 0.;
    if (aligned(list_intersection_points)) return 0.;

    // Sort the inside and intersection points to form a polygon
    std::vector<std::vector<double>> list_points_polygon = sort_points(list_intersection_points);

    // Compute the area of the polygon
    double area = 0.;
    for (k1 = 0; k1 < number_points; k1++)
    {
//if (k == 1830) std::cout << list_points_polygon[k1][0] << " " << list_points_polygon[k1][1] << std::endl;
        area += list_points_polygon[k1][0] * list_points_polygon[(k1+1)%number_points][1] - list_points_polygon[k1][1] * list_points_polygon[(k1+1)%number_points][0];
    }

    return 0.5*fabs(area);

} // intersection_area

}// nextsim
