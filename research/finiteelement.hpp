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
#include <gmshmesh.hpp>


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
    double measure(element_type const& element) const;
    void assemble();
    void solve();
    void run();


private:
    po::variables_map vm;
    mesh_type M_mesh;
    matrix_ptrtype M_matrix;
    vector_ptrtype M_vector;
    vector_ptrtype M_solution;

    std::map<int,std::vector<double> > M_local_contrib;
    std::map<int,std::vector<double> > M_local_rows;
    std::map<int,std::vector<double> > M_local_cols;

    std::map<int, point_type > M_nodes;
    std::map<int, element_type > M_lines;
    std::map<int, element_type > M_elements;

    int M_num_nodes;
    int M_num_elements;

    std::vector<int> dirichlet_flags;
    boost::timer chrono;
};
} // Nextsim
#endif
