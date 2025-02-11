/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   metric.cpp
 * @author Fabien Salmon <fabien.salmon@inria.fr>
 * @date   Thu Aug  4 09:47:27 CEST 2024
 */

#include "metric.hpp"
#include <vector>

/**
 * @class metric
 * @brief Compute the metric for anisotropic remeshing with mmg
 *
 * @see
 *
 */

namespace Nextsim
{
// Member functions
// Edge length in the non-euclidean geometry defined by the metric
double 
Metric::metric_square_edge_length(int id1, int id2, std::vector<double> const pt1, std::vector<double> const pt2) const
{
    // Compute the average metric on the edge
    std::vector<double> average_metric(3);
    for (int i = 0; i<= 2; i++)
    {
        average_metric[i] = 0.5*(this->components[id1][i] + this->components[id2][i]);
    }

    // Compute non-euclidean length with the average metric
    std::vector<double> u(2);
    u[0] = pt2[0] - pt1[0];
    u[1] = pt2[1] - pt1[1];

    double res = u[0] * (average_metric[0] * u[0] + average_metric[1] * u[1]) + u[1] * (average_metric[1] * u[0] + average_metric[2] * u[1]);
    
    return res;
} // metric_edge_length

// Areas of the triangles in the non-euclidean geometry
template<typename FEMeshType>
std::vector<double>
Metric::metric_triangle_areas(FEMeshType const& mesh) const
{
    std::vector<double> barycenter(2);
    std::vector<std::vector<double>> points(3, std::vector<double>(2, 0.));
    std::vector<double> inverse_distance(3);
    std::vector<double> average_metric(3);
    std::vector<int> ids(3);
    std::vector<double> areas = compute_triangle_euclidean_areas(mesh);
    std::vector<double> res(mesh.numTriangles());

    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        // points of the triangle
        for (int i = 0; i <= 2; i++)
        {
            ids[i] = mesh.triangles()[n].indices[i] - 1;
            points[i] = nodes_coordinates(mesh, n, i);
        }

        // Barycenter coordinates
        barycenter[0] = (points[0][0] + points[1][0] + points[2][0])/3.;
        barycenter[1] = (points[0][1] + points[1][1] + points[2][1])/3.;

        // Distance of the points from the barycenter
        for (int i = 0; i <= 2; i++)
        {   
            inverse_distance[i] = 1./std::max(sqrt(pow(points[i][0]-barycenter[0],2)+pow(points[i][1]-barycenter[1],2)), std::numeric_limits<double>::epsilon());
        }

        // Metric at the barycenter as the average of the metrics at the three vertices
        average_metric = {0.,0.,0.};
        for (int i = 0; i <= 2; i++)
        {
            for (int j = 0; j <= 2; j++)
            {
                average_metric[i] += inverse_distance[j] * this->components[ids[j]][i];
            }
            average_metric[i] /= std::accumulate(std::cbegin(inverse_distance), std::cend(inverse_distance), 0.);
        }

        // Compute the determinant of the average_metric
        double det = average_metric[0]*average_metric[2] - average_metric[1]*average_metric[1];

        res[n] = areas[n] * sqrt(det);
    }

    return res;
} // metric_triangle_areas

// Quality function
template<typename FEMeshType>
std::vector<double> 
Metric::triangle_quality(FEMeshType const& mesh) const
{
    double const coef = 4.*sqrt(3); // coef to have a quality equal to one for unit elements
    std::vector<double> non_euclidean_areas = metric_triangle_areas(mesh);
    std::vector<double> edge_length_square(3);
    std::vector<std::vector<double>> points(3, std::vector<double>(2, 0));
    std::vector<int> ids(3);
    std::vector<double> res(mesh.numTriangles());

    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        // If the area is null, the quality of the triangle is bad so a great result is returned
        if (non_euclidean_areas[n] < std::numeric_limits<double>::epsilon()) 
        {
            res[n] = 1.e10;
            continue;
        }

        for (int i = 0; i <= 2; i++)
        {
            ids[i] = mesh.triangles()[n].indices[i] - 1;
            points[i] = nodes_coordinates(mesh, n, i);
        }

        for (int i = 0; i <= 2; i++)
        {
            edge_length_square[i] = metric_square_edge_length(ids[(i+1)%3], ids[i], points[(i+1)%3], points[i]);
        }

        res[n] = coef * non_euclidean_areas[n] / std::accumulate(std::cbegin(edge_length_square), std::cend(edge_length_square), 0.);

    }

    return res;

} // triangle_quality

// Return a list of triangles that do not respect the quality tolerance
std::vector<int> 
Metric::check_triangle_quality(GmshMesh const& mesh) const
{
    std::vector<int> res {};
    std::vector<double> quality = triangle_quality(mesh);

    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        if (quality[n] < this->quality_threshold)
        {
            res.push_back(n);
        }
    }

    return res;

} // check_triangle_quality

std::vector<int>
Metric::check_triangle_quality(GmshMeshSeq const& mesh) const
{
    std::vector<int> res {};
    std::vector<double> quality = triangle_quality(mesh);

    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        if (quality[n] < this->quality_threshold)
        {
            res.push_back(n);
        }
    }

    return res;

} // check_triangle_quality

// Metric Loseille solution of the minimization problem
void
Metric::compute_optimal_metric(GmshMesh const& mesh, std::vector<double> const& field, double hmin, double hmax)
{

    // First step: compute triangle areas
    std::vector<double> areas_triangle = compute_triangle_euclidean_areas(mesh);

    // Second step: compute the area associated to each vertex. 
    // This area corresponds to the sum of one third of the area of each triangle contaning the vertex.
    std::vector<double> areas_vertex(mesh.numNodes(),0.);
    int id;
    
    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        for (int i = 0; i <= 2; i++)
        {
            id = mesh.triangles()[n].indices[i] - 1;
            areas_vertex[id] += areas_triangle[n]/3.;
        }
    }

    // Third step: compute the Hessian and its determinant
    std::vector<std::vector<double>> hessian_components = compute_hessian_field(mesh, field);
    std::vector<double> hessian_determinant = compute_determinant_field(hessian_components);

    // Fourth step: compute the integral in the relation giving the optimal metric field
    double integral = 0.;
    for (int n = 0; n < mesh.numNodes(); n++)
    {
        integral += pow(hessian_determinant[n], this->p/(2.*this->p+3.)) * areas_vertex[n];
    }

    // Last step: compute the optimized metric tensor field
    std::vector<std::vector<double>> res(mesh.numNodes(), std::vector<double>(3));
    for (int n = 0; n < mesh.numNodes(); n++)
    {
        for (int i = 0; i <= 2; i++)
        {
            res[n][i] = pow(integral / this->Nst, -2./3.) * pow(hessian_determinant[n], -1/(2.*this->p+3.)) * hessian_components[n][i];
        }

        // Apply a threshold to avoid too small of too big elements
        res[n][0] = std::min(std::max(res[n][0], 1./(hmax*hmax)), 1./(hmin*hmin));
        res[n][2] = std::min(std::max(res[n][2], 1./(hmax*hmax)), 1./(hmin*hmin));
    }
    
    this->components = res;

    // Write a vtk file containing the mesh and the metric
    if (this->is_vtk_written) write_metric_vtk(mesh);

} // compute_optimal_metric

void
Metric::compute_optimal_metric(GmshMeshSeq const& mesh, std::vector<double> const& field, double hmin, double hmax)
{
    // First step: compute triangle areas
    std::vector<double> areas_triangle = compute_triangle_euclidean_areas(mesh);

    // Second step: compute the area associated to each vertex. 
    // This area corresponds to the sum of one third of the area of each triangle contaning the vertex.
    std::vector<double> areas_vertex(mesh.numNodes(),0.);
    int id;

    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        for (int i = 0; i <= 2; i++)
        {
            id = mesh.triangles()[n].indices[i] - 1;
            areas_vertex[id] += areas_triangle[n]/3.;
        }
    }

    // Third step: compute the Hessian and its determinant
    std::vector<std::vector<double>> hessian_components = compute_hessian_field(mesh, field);
    std::vector<double> hessian_determinant = compute_determinant_field(hessian_components);

    // Fourth step: compute the integral in the relation giving the optimal metric field
    double integral = 0.;
    for (int n = 0; n < mesh.numNodes(); n++)
    {
        integral += pow(hessian_determinant[n], this->p/(2.*this->p+2.)) * areas_vertex[n];
    }

    // Last step: compute the optimized metric tensor field
    std::vector<std::vector<double>> res(mesh.numNodes(), std::vector<double>(3));
    for (int n = 0; n < mesh.numNodes(); n++)
    {
        for (int i = 0; i <= 2; i++)
        {
            res[n][i] = this->Nst / integral * pow(hessian_determinant[n], -1/(2.*this->p+2.)) * hessian_components[n][i];
        }

        // Apply a threshold to avoid too small of too big elements
        res[n][0] = std::min(std::max(res[n][0], 1./(hmax*hmax)), 1./(hmin*hmin));
        res[n][2] = std::min(std::max(res[n][2], 1./(hmax*hmax)), 1./(hmin*hmin));
    }

    this->components = res;

    // Write a vtk file containing the mesh and the metric
    if (this->is_vtk_written) write_metric_vtk(mesh);

} // compute_optimal_metric

// Metric Salmon solution of the minimization problem
double 
Metric::error_min_max(std::vector<double> const& areas_vertex, std::vector<int> const& order_number, 
                      std::vector<double> const& sum_gamma, int nb_vertices, double d_min, double d_max, double alpha, double z1)
{
    double z2 = z1 * pow(d_max/d_min, 1./alpha);

    double area_A = 0.;
    double area_B = 0.;
    double integral = 0.;

    for (int k = 0; k < nb_vertices; k++) 
    {
        int n = order_number[k];
        if (sum_gamma[n] < z1)
        {
            area_A += areas_vertex[n];
        }
        else if (sum_gamma[n] < z2)
        {
            integral += pow(sum_gamma[n],alpha) * areas_vertex[n];
        }
        else
        {
            area_B += areas_vertex[n];
        }
    }

    return (this->Nst - area_A*d_min - area_B*d_max) * pow(z1, alpha) / integral - d_min;

} // error_min_max

double
Metric::error_max(std::vector<double> const& areas_vertex, std::vector<int> const& order_number, 
                  std::vector<double> const& sum_gamma, int nb_vertices, double d_max, double alpha, double z2)
{

    double area_B = 0.;
    double integral = 0.;

    for (int k = 0; k < nb_vertices; k++)
    {
        int n = order_number[k];
        if (sum_gamma[n] < z2)
        {
            integral += pow(sum_gamma[n],alpha) * areas_vertex[n];
        }
        else
        {
            area_B += areas_vertex[n];
        }
    }

    return (this->Nst - area_B*d_max) * pow(z2, alpha) / integral - d_max;

} // error_max

void
Metric::compute_optimal_bounded_isotropic_metric(GmshMesh const& mesh, std::vector<double> const& field, double hmin, double hmax, Communicator const & comm)
{
    double d_max = 1./(hmin*hmin);
    double d_min = 1./(hmax*hmax);
    double alpha = 1./(1.+this->p);

    // First step: compute triangle areas
    std::vector<double> areas_triangle = compute_triangle_euclidean_areas(mesh);

    // Second step: compute the area associated to each vertex. 
    // This area corresponds to the sum of one third of the area of each triangle contaning the vertex.
    std::vector<double> areas_vertex(mesh.numNodes(),0.);
    int id;

    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        for (int i = 0; i <= 2; i++)
        {
            id = mesh.triangles()[n].indices[i] - 1;
            areas_vertex[id] += areas_triangle[n]/3.;
        }
    }

    // Third step: compute the Hessian
    std::vector<std::vector<double>> hessian_components = compute_hessian_field(mesh, field);

    // Fourth step: compute the second term in the integral error and order it
    std::vector<double> sum_gamma(mesh.numNodes());
    for (int n = 0; n < mesh.numNodes(); n++) 
    {
        sum_gamma[n] = pow(hessian_components[n][0] + hessian_components[n][2], this->p);
    }

    std::vector<int> sizes(comm.size());
    boost::mpi::all_gather(comm, mesh.numNodes(), sizes);

    int total_size = std::accumulate(sizes.begin(), sizes.end(), 0);

    std::vector<double> global_sum_gamma(total_size);
    boost::mpi::gatherv(comm, sum_gamma, &global_sum_gamma[0], sizes, 0);
    boost::mpi::broadcast(comm, &global_sum_gamma[0], total_size, 0);

    std::vector<int> order_number(total_size);
    std::iota(order_number.begin(), order_number.end(), 0); 
    std::sort(order_number.begin(), order_number.end(), [&global_sum_gamma](int i, int j) {
              return global_sum_gamma[i] < global_sum_gamma[j];}
             );

    std::vector<double> global_areas_vertex(total_size);
    boost::mpi::gatherv(comm, areas_vertex, &global_areas_vertex[0], sizes, 0);
    boost::mpi::broadcast(comm, &global_areas_vertex[0], total_size, 0);

    // Fifth step: find z1 and z2 by dichotomy 
    // NOTE: the method might be improved but error_min_max and error_max are constant piecewise functions
    double z1,z2;
    int n = 0;
    int N_ITER = 100;
    auto a = *std::min_element(global_sum_gamma.begin(), global_sum_gamma.end());
    auto b = *std::max_element(std::begin(global_sum_gamma), std::end(global_sum_gamma));
    double TOL = 1.e-15;
    double c;
    if (error_min_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_min, d_max, alpha, a) *
        error_min_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_min, d_max, alpha, b) > 0.)
    {
        // d_min is smaller than the smallest element using Loseille solution
        a += (b-a)/10000.;
        if (error_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_max, alpha, a) *
            error_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_max, alpha, b) > 0.)
        {
            // d_max is greater than the greatest element using Loseille solution
            z1 = 0.;
            z2 = 1.e40;
        }
        else
        {
            while (fabs(error_min_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_min, d_max, alpha, a) - 
                        error_min_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_min, d_max, alpha, b)) > TOL && 
                   n < N_ITER)
            {
                double c = (a+b) / 2.;
                if (error_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_max, alpha, a) *
                    error_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_max, alpha, c) < 0.)
                {
                    b = c;
                }
                else 
                {
                    a = c;
                }
                n += 1;
            }

            z2 = (a+b) / 2.;
            z1 = 0.;
        }
    }
    else
    {
        while (fabs(error_min_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_min, d_max, alpha, a) - 
                    error_min_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_min, d_max, alpha, b)) > TOL &&
               n < N_ITER)
        {
            double c = (a+b) / 2.;
            if (error_min_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_min, d_max, alpha, a) *
                error_min_max(global_areas_vertex, order_number, global_sum_gamma, total_size, d_min, d_max, alpha, c) < 0.)
            {
                b = c;
            }
            else
            {
                a = c;
            }
            n += 1;
        }

        z1 = (a+b) / 2.;
        z2 = z1 * pow(d_max/d_min, 1./alpha);
    }

    // Sixth step: compute the eigenvalue that minimize the integral error
    double area_A = 0.;
    double area_B = 0.;
    double integral = 0.;
    double error = 0.;
    double constraint = 0.;
    std::vector<double> lambda(total_size);
    for (int k = 0; k < total_size; k++) 
    {
        n = order_number[k];
        if (global_sum_gamma[n] < z1)
        {
            area_A += global_areas_vertex[n];
        }
        else if (global_sum_gamma[n] < z2)
        {
            integral += pow(global_sum_gamma[n], alpha) * global_areas_vertex[n];
        }
        else
        {
            area_B += global_areas_vertex[n];
        }
    }

    std::vector<int> displs(comm.size(), 0);
    for (int k = 1; k < comm.size(); ++k) {
        displs[k] = displs[k - 1] + sizes[k - 1];
    }
    for (int k = 0; k < total_size; k++)
    {
        n = order_number[k];
        if (global_sum_gamma[n] < z1)
        {
            lambda[n] = d_min;
        }
        else if (global_sum_gamma[n] < z2)
        {
            lambda[n] = (this->Nst - area_A*d_min -area_B*d_max) * pow(global_sum_gamma[n], alpha) / integral;
        }
        else
        {
            lambda[n] = d_max;
        }
        //constraint += lambda[n]*global_areas_vertex[n];
        //error += sum_gamma[n]*pow(lambda[n],-this->p)*areas_vertex[n];
    }

    //double result = 0.;
    //boost::mpi::reduce(comm, constraint, result, std::plus<double>(), 0);
    //if (comm.rank() == 0) std::cout << "Both should be equal " << constraint << " " << this->Nst << std::endl;

    // Last step: compute the metric based on the eigenvalue lambda
    std::vector<double> v1(2);
    std::vector<double> v2(2);
    std::vector<std::vector<double>> res(mesh.numNodes(), std::vector<double>(3));
    for (n = 0; n < mesh.numNodes(); n++)
    {
        // Eigen vectors of the Hessian matrix
        a = hessian_components[n][0];
        b = hessian_components[n][1];
        c = hessian_components[n][2];

        // Check whether the hessian is diagonal
        if (fabs(b) < 1e-12)
        {
            v1[0] = 1.;
            v1[1] = 0.;
            v2[1] = 1.;
            v2[0] = 0.;
        }
        else
        {
            // Eigen vectors
            double delta = sqrt((a-c)*(a-c) + 4*b*b);
            double m = (-a+c+delta)/(2*b);
            double inv_det = 1./sqrt(1+m*m);
            v1[0] = inv_det;
            v1[1] = m*inv_det;
            v2[0] = -m*inv_det;
            v2[1] = inv_det;
        }

        int l = n + displs[comm.rank()];
        res[n][0] = lambda[l]*pow(v1[0],2) + lambda[l]*pow(v1[1],2);
        res[n][2] = lambda[l]*pow(v2[0],2) + lambda[l]*pow(v2[1],2);
        res[n][1] = lambda[l]*v1[0]*v2[0] + lambda[l]*v1[1]*v2[1];

    }

    this->components = res;

    // Write a vtk file containing the mesh and the metric
    if (this->is_vtk_written) write_metric_vtk(mesh);

} // compute_optimal_bounded_metric

void
Metric::compute_optimal_bounded_isotropic_metric(GmshMeshSeq const& mesh, std::vector<double> const& field, double hmin, double hmax, Communicator const & comm)
{
    double d_max = 1./(hmin*hmin);
    double d_min = 1./(hmax*hmax);
    double alpha = 1./(1.+this->p);

    // First step: compute triangle areas
    std::vector<double> areas_triangle = compute_triangle_euclidean_areas(mesh);

    // Second step: compute the area associated to each vertex. 
    // This area corresponds to the sum of one third of the area of each triangle contaning the vertex.
    std::vector<double> areas_vertex(mesh.numNodes(),0.);
    int id;

    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        for (int i = 0; i <= 2; i++)
        {
            id = mesh.triangles()[n].indices[i] - 1;
            areas_vertex[id] += areas_triangle[n]/3.;
        }
    }

    // Third step: compute the Hessian
    std::vector<std::vector<double>> hessian_components = compute_hessian_field(mesh, field);

    // Fourth step: compute the second term in the integral error and order it
    std::vector<double> sum_gamma(mesh.numNodes());
    for (int n = 0; n < mesh.numNodes(); n++)
    {
        sum_gamma[n] = pow(hessian_components[n][0] + hessian_components[n][2], this->p);
    }

    std::vector<int> order_number(mesh.numNodes());
    std::iota(order_number.begin(), order_number.end(), 0);
    std::sort(order_number.begin(), order_number.end(), [&sum_gamma](int i, int j) {
              return sum_gamma[i] < sum_gamma[j];}
             );

    // Fifth step: find z1 and z2 by dichotomy 
    // NOTE: the method might be improved but error_min_max and error_max are constant piecewise functions
    double z1,z2;
    int n = 0;
    int N_ITER = 100;
    double TOL = 1.e-15;
    auto a = *std::min_element(sum_gamma.begin(), sum_gamma.end());
    auto b = *std::max_element(std::begin(sum_gamma), std::end(sum_gamma));
    double c;
    if (error_min_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_min, d_max, alpha, a) *
        error_min_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_min, d_max, alpha, b) > 0.)
    {
        // d_min is smaller than the smallest element using Loseille solution
        a += (b-a)/10000.;
        if (error_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_max, alpha, a) * 
            error_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_max, alpha, b) > 0.)
        {
            // d_max is greater than the greatest element using Loseille solution
            z1 = 0.;
            z2 = 0.;
        }
        else
        {
            while (fabs(error_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_max, alpha, a) -
                        error_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_max, alpha, b)) > TOL &&
                   n < N_ITER)
            {
                double c = (a+b) / 2.;
                if (error_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_max, alpha, a) *
                    error_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_max, alpha, c) < 0.)
                {
                    b = c;
                }
                else
                {
                    a = c;
                }
                n += 1;
            }

            z2 = (a+b) / 2.;
            z1 = 0.;
        }
    }
    else
    {
        while (fabs(error_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_max, alpha, a) -
                    error_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_max, alpha, b)) > TOL &&
               n < N_ITER)
        {
            double c = (a+b) / 2.;
            if (error_min_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_min, d_max, alpha, a) *
                error_min_max(areas_vertex, order_number, sum_gamma, mesh.numNodes(), d_min, d_max, alpha, c) < 0.)
            {
                b = c;
            }
            else
            {
                a = c;
            }
            n += 1;
        }

        z1 = (a+b) / 2.;
        z2 = z1 * pow(d_max/d_min, 1./alpha);
    }

    // Sixth step: compute the eigenvalue that minimize the integral error
    double area_A = 0.;
    double area_B = 0.;
    double integral = 0.;
    double error = 0.;
    std::vector<double> lambda(mesh.numNodes());
    for (int k = 0; k < mesh.numNodes(); k++)
    {
        n = order_number[k];
        if (sum_gamma[n] < z1)
        {
            area_A += areas_vertex[n];
        }
        else if (sum_gamma[n] < z2)
        {
            integral += pow(sum_gamma[n], alpha) * areas_vertex[n];
        }
        else
        {
            area_B += areas_vertex[n];
        }
    }

    for (int k = 0; k < mesh.numNodes(); k++)
    {
        n = order_number[k];
        if (sum_gamma[n] < z1)
        {
            lambda[n] = d_min;
        }
        else if (sum_gamma[n] < z2)
        {
            lambda[n] = (this->Nst - area_A*d_min -area_B*d_max) * pow(sum_gamma[n], alpha) / integral;
        }
        else
        {
            lambda[n] = d_max;
        }

        //error += sum_gamma[n]*pow(lambda[n],-this->p)*areas_vertex[n];
    }

    // Last step: compute the metric based on the eigenvalue lambda
    std::vector<double> v1(2);
    std::vector<double> v2(2);
    std::vector<std::vector<double>> res(mesh.numNodes(), std::vector<double>(3));
    for (n = 0; n < mesh.numNodes(); n++)
    {
        // Eigen vectors of the Hessian matrix
        a = hessian_components[n][0];
        b = hessian_components[n][1];
        c = hessian_components[n][2];

        // Check whether the hessian is diagonal
        if (fabs(b) < 1e-12)
        {
            v1[0] = 1.;
            v1[1] = 0.;
            v2[1] = 1.;
            v2[0] = 0.;
        }
        else
        {
            // Eigen vectors
            double delta = sqrt((a-c)*(a-c) + 4*b*b);
            double m = (-a+c+delta)/(2*b);
            double inv_det = 1./sqrt(1+m*m);
            v1[0] = inv_det;
            v1[1] = m*inv_det;
            v2[0] = -m*inv_det;
            v2[1] = inv_det;
        }

        res[n][0] = lambda[n]*pow(v1[0],2) + lambda[n]*pow(v1[1],2);
        res[n][2] = lambda[n]*pow(v2[0],2) + lambda[n]*pow(v2[1],2);
        res[n][1] = lambda[n]*v1[0]*v2[0] + lambda[n]*v1[1]*v2[1];

    }

    this->components = res;

    // Write a vtk file containing the mesh and the metric
    if (this->is_vtk_written) write_metric_vtk(mesh);

} // compute_optimal_bounded_metric

// Write the mesh vertices and the metric as a VTK file for debug
template<typename FEMeshType>
void
Metric::write_metric_vtk(FEMeshType const& mesh) const
{
    // One file per process
    std::stringstream name_file {};
    name_file << "mesh_metric_";
    name_file << this->process;
    name_file << ".vtk";
    std::string file_name;
    name_file >> file_name;
    std::ofstream vtkfile {file_name};

    vtkfile << "# vtk DataFile Version 2.0\nmesh_metric\nASCII\nDATASET POLYDATA\nPOINTS " << mesh.numNodes() << " float\n";

    std::vector<double> coordX = mesh.coordX();
    std::vector<double> coordY = mesh.coordY();

    for (int n = 0; n < mesh.numNodes(); n++)
    {
        vtkfile << coordX[n] << " " << coordY[n] << " 0.\n";
    }

    vtkfile << "VERTICES " << mesh.numNodes() << " " << 2*mesh.numNodes() <<"\n";

    for (int n = 0; n < mesh.numNodes(); n++)
    {
        vtkfile << "1 " << n << "\n";
    }

    vtkfile << "POINT_DATA " << mesh.numNodes() << "\n";
    for (int j = 0; j <= 2; j++)
    {
        vtkfile << "SCALARS metric" << j << " float\nLOOKUP_TABLE default\n";
        for (int n = 0; n < mesh.numNodes(); n++)
        {
            vtkfile << this->components[n][j] <<"\n";
        }
    }

    vtkfile.close();

} // write_metric_vtk

// Compute the areas of each triangle of the mesh (Euclidean calculus)
template<typename FEMeshType>
std::vector<double>
Metric::compute_triangle_euclidean_areas(FEMeshType const& mesh) const
{
    std::vector<std::vector<double>> points(3, std::vector<double>(2, 0));
    std::vector<double> u(2);
    std::vector<double> v(2);
    std::vector<double> area(mesh.numTriangles());

    for (int n = 0; n < mesh.numTriangles(); n++)
    {

        for (int i = 0; i <= 2; i++)
        {
            points[i] = nodes_coordinates(mesh, n, i);
        }

        for (int i = 0; i <= 1; i++)
        {
            u[i] = points[1][i] - points[0][i];
            v[i] = points[2][i] - points[0][i];
        }

        area[n] = 0.5 * abs(u[0] * v[1] - v[0] * u[1]);
    }

    return area;
} // compute_triangle_euclidean_areas

// Compute an average of a scalar field over several triangles with weights associated to their areas
double
Metric::area_weighted_average(std::vector<double> const& scalar_field, std::vector<double> const& area, std::vector<int> const& list_triangles) const
{
    double sum_field = 0.;
    double sum_area = 0.;

    for (auto it = list_triangles.begin(); it != list_triangles.end(); ++it)
    {
        sum_field += area[*it] * scalar_field[*it];
        sum_area += area[*it];
    }

    return sum_field/sum_area;
} // area_weighted_average

// Compute the gradient of a field over a triangle defined by the coordinates of its points
std::vector<double> 
Metric::compute_grad(std::vector<double> const& scalar_field, std::vector<std::vector<double>> const& points_coord, double const area) const
{
    double dx;
    double dy;

    std::vector<double> res {0,0};

    for (int i = 0; i<= 2; i++)
    {
        dx = points_coord[(i+1)%3][0] - points_coord[(i+2)%3][0];
        dy = points_coord[(i+1)%3][1] - points_coord[(i+2)%3][1];
        res[0] += scalar_field[i] * dy;
        res[1] -= scalar_field[i] * dx;
    }

    res[0] /= 2*area;
    res[1] /= 2*area;

    return res;
} // compute_grad

// Compute the Hessian of a field over a triangle (same procedure as for the gradient but using the gradient field as input)
std::vector<double> 
Metric::compute_hessian(std::vector<std::vector<double>> const& gradient_field, std::vector<std::vector<double>> const& points_coord, double const area) const
{
    double dx;
    double dy;

    std::vector<double> hessian(4, 0);

    for (int i = 0; i <= 2; i++)
    {
        dx = points_coord[(i+1)%3][0] - points_coord[(i+2)%3][0];
        dy = points_coord[(i+1)%3][1] - points_coord[(i+2)%3][1];

        hessian[0] += gradient_field[i][0] * dy;
        hessian[1] -= gradient_field[i][0] * dx;
        hessian[2] += gradient_field[i][1] * dy;
        hessian[3] -= gradient_field[i][1] * dx;
    }

    // To enforce the symmetry of the Hessian, the non-diagonal terms are equal to 1/2 (h_xy + h_yx)
    std::vector<double> res(3);
    res[0] = hessian[0]/(2.*area);
    res[1] = (hessian[1] + hessian[2])/(4.*area);
    res[2] = hessian[3]/(2.*area);

    return res;
}   

// Compute the list of triangles containing each vertex
template<typename FEMeshType>
std::vector<std::vector<int>> 
Metric::compute_list_triangles(FEMeshType const& mesh) const
{
    std::vector<std::vector<int>> list_triangles;
    list_triangles.resize(mesh.numNodes());

    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        for (int k = 0; k <= 2; k++)
        {
            list_triangles[mesh.triangles()[n].indices[k]-1].push_back(n); // The node indices begin with 1
        }
    }

    return list_triangles;
} // compute_list_triangles

// Compute the Hessian matrix at each vertex for a chosen field
template<typename FEMeshType>
std::vector<std::vector<double>>
Metric::compute_hessian_field(FEMeshType const& mesh, std::vector<double> const& field) const
{

    // Preliminary: compute the list of triangles associated to each vertex on which the average must be achieved
    std::vector<std::vector<int>> list_triangles = compute_list_triangles(mesh);

    // Then compute the Euclidean area of each triangle
    std::vector<double> areas = compute_triangle_euclidean_areas(mesh);

    // Check whether the field is on vertices or triangles. In the second case, make interpolations to have a node field
    std::vector<double> field_copy;
    if (field.size() == mesh.numTriangles())
    {
        std::vector<double> field_nodes(mesh.numNodes());
        for (int n = 0; n < mesh.numNodes(); n++)
        {
            field_nodes[n] = area_weighted_average(field, areas, list_triangles[n]);
        }
        field_copy = field_nodes;
    }
    else
    {
        field_copy = field;
    }

    // Firs step: compute the gradient of the field on each triangle
    std::vector<std::vector<double>> points(3, std::vector<double>(2, 0));
    std::vector<double> sub_field(3);
    std::vector<double> gradient_triangle(2);
    std::vector<double> gradient_triangle_x(mesh.numTriangles());
    std::vector<double> gradient_triangle_y(mesh.numTriangles());
    
    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        // Coordinates of the three points of the triangle and field at the three points
        for (int i = 0; i <= 2; i++)
        {  
            points[i] = nodes_coordinates(mesh, n, i); 
            sub_field[i] = field_copy[mesh.triangles()[n].indices[i]-1];
        }
        
        gradient_triangle = compute_grad(sub_field, points, areas[n]);
        gradient_triangle_x[n] = gradient_triangle[0];
        gradient_triangle_y[n] = gradient_triangle[1];
    }

    // Second step: area weighted average of the gradient field at each vertex
    std::vector<std::vector<double>> average_gradient(mesh.numNodes(), std::vector<double>(2, 0));

    for (int n = 0; n < mesh.numNodes(); n++)
    {
        average_gradient[n][0] = area_weighted_average(gradient_triangle_x, areas, list_triangles[n]);
        average_gradient[n][1] = area_weighted_average(gradient_triangle_y, areas, list_triangles[n]);
    }
 
    // Third step: compute the Hessian of the field on each triangle

    std::vector<double> Hessian_triangle(3);
    std::vector<std::vector<double>> sub_gradient_field(3,std::vector<double>(2,0));
    std::vector<double> Hessian_triangle_xx(mesh.numTriangles());
    std::vector<double> Hessian_triangle_xy(mesh.numTriangles());
    std::vector<double> Hessian_triangle_yy(mesh.numTriangles());

    for (int n = 0; n < mesh.numTriangles(); n++)
    {
        // Coordinates of the three points of the triangle
        for (int i = 0; i <= 2; i++)
        {
            points[i] = nodes_coordinates(mesh, n, i);
            sub_gradient_field[i] = average_gradient[mesh.triangles()[n].indices[i]-1];
        }

        Hessian_triangle = compute_hessian(sub_gradient_field, points, areas[n]);
        Hessian_triangle_xx[n] = Hessian_triangle[0];
        Hessian_triangle_xy[n] = Hessian_triangle[1];
        Hessian_triangle_yy[n] = Hessian_triangle[2];
    }

    // Fourth step: area weighted average of the Hessian at each vertex

    std::vector<std::vector<double>> average_hessian(mesh.numNodes(), std::vector<double>(3, 0));
    for (int n = 0; n < mesh.numNodes(); n++)
    {
        average_hessian[n][0] = area_weighted_average(Hessian_triangle_xx, areas, list_triangles[n]);
        average_hessian[n][1] = area_weighted_average(Hessian_triangle_xy, areas, list_triangles[n]);
        average_hessian[n][2] = area_weighted_average(Hessian_triangle_yy, areas, list_triangles[n]);
    }

    // Fifth step: diagonalize the Hessian at each vertex and replace the eigenvalues by their absolute value
    std::vector<std::vector<double>> hessian = diagonalize_hessian_positive(average_hessian);

    return hessian; 
} // compute_hessian_field

// Diagonalize the Hessian matrices and replace the eigenvalues by their absolute value
std::vector<std::vector<double>>
Metric::diagonalize_hessian_positive(std::vector<std::vector<double>> const& matrix) const
{
    double lambda1, lambda2, delta, a, b, c, m, inv_det;
    std::vector<double> v1(2), v2(2);
    std::vector<std::vector<double>> res(matrix.size(), std::vector<double>(3, 0));

    for (int n = 0; n < matrix.size(); n++)
    {
        // Coefficients of the matrix to make the code readable
        a = matrix[n][0];
        b = matrix[n][1];
        c = matrix[n][2];

        // First check whether the matrix is already diagonal
        if (b < std::numeric_limits<double>::epsilon())
        {
            // A threshold is simultaneously applied to avoid tiny eigenvalues
            res[n][0] = std::max(abs(a), std::numeric_limits<double>::epsilon());
            res[n][2] = std::max(abs(c), std::numeric_limits<double>::epsilon());

            continue;
        }
        
        // Eigen values
        delta = pow((a-c)*(a-c) + 4*b*b, 0.5);
        lambda1 = 0.5*(a+c+delta);
        lambda2 = 0.5*(a+c-delta);

        // Eigen vectors
        m = (-a+c+delta)/(2*b);
        inv_det = 1./pow(1+m*m,0.5);
        v1[0] = inv_det;
        v1[1] = m*inv_det;
        v2[0] = -m*inv_det;
        v2[1] = inv_det;

        // SPD hessian matrix - A threshold is simultaneously applied to avoid tiny eigenvalues
        lambda1 = std::max(abs(lambda1), std::numeric_limits<double>::epsilon());
        lambda2 = std::max(abs(lambda2), std::numeric_limits<double>::epsilon());
        res[n][0] = v1[0]*v1[0]*lambda1 + v1[1]*v1[1]*lambda2;
        res[n][1] = v1[0]*v2[0]*lambda1 + v1[1]*v2[1]*lambda2;
        res[n][2] = v2[0]*v2[0]*lambda1 + v2[1]*v2[1]*lambda2;
    }

    return res;

} // diagonalize_hessian_positive

// Compute the determinant of a matrix associated to each vertex
std::vector<double>
Metric::compute_determinant_field(std::vector<std::vector<double>> const& matrix) const
{
    std::vector<double> determinant(matrix.size()); 
    
    for (int i = 0; i < determinant.size(); i++)
    {
        determinant[i] = matrix[i][0] * matrix[i][2] - matrix[i][1] * matrix[i][1];
    }
    
    return determinant;
} // compute_determinant_field


} // Nextsim
