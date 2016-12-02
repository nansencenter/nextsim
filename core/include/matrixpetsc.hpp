/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   matrixpetsc.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Wed Jul  1 16:33:49 2015
 */

#ifndef __MatrixPetsc_H
#define __MatrixPetsc_H 1

#include <environment.hpp>
#include <graphcsr.hpp>
#include <graphcsrmpi.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <assert.hpp>


extern "C"
{
#include <petscmat.h>
}

# include <petscversion.h>
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
#   include <petscsles.h>
# else
#   include <petscksp.h>
# endif


/**
 * @class MatrixPetsc
 * @brief Wrapper for petsc matrices
 *
 * @see
 *
 */

namespace Nextsim
{

class VectorPetsc;

class MatrixPetsc
{

public:

    typedef GraphCSR graph_type;
    typedef GraphCSRMPI graphmpi_type;
    typedef std::size_t size_type;
    typedef double value_type;

    MatrixPetsc( Communicator const& comm = Environment::comm() );

    ~MatrixPetsc();

    Communicator const& comm() const { return M_comm; }

    Mat mat () const;

    Mat& mat ();

    void init( const size_type m, const size_type n, const size_type nnz );
    void init( const size_type m, const size_type n, graph_type const& graph );

    void init( const size_type m,
               const size_type n,
               const size_type m_l,
               const size_type n_l,
               graphmpi_type const& graph );

    void initLocalToGlobalMapping(graphmpi_type const& graph);

    void zero();

    void close() const;

    void setValue(size_type const& i, size_type const& j, value_type const& value);

    void addValue(size_type const& i, size_type const& j, value_type const& value);

    void setMatrix(int* rows, int nrows, int* cols, int ncols, value_type* data);

    void addMatrix(int* rows, int nrows, int* cols, int ncols, value_type* data);

    value_type l1Norm() const;

    value_type linftyNorm() const;

    value_type operator () (size_type const& i, size_type const& j) const;

    MatrixPetsc& operator = ( MatrixPetsc const& M );

    size_type size1() const;

    size_type size2() const;

    void printScreen() const;

    void printMatlab(std::string const& filename = "NULL") const;

    void printBinary(std::string const& filename = "NULL") const;

    void loadBinary(std::string const& filename = "NULL");

    void clear();

    bool isInitialized() const { return M_is_initialized; }

    bool closed() const;

    bool isSymmetric() const;

    void diagonal( VectorPetsc& out ) const;

    void on(std::vector<int> const& flags, VectorPetsc& rhs);
    //void on(std::vector<int> const& flags, VectorPetsc& rhs, graphmpi_type const& graph);

    value_type energy(VectorPetsc& u) const;

    size_type rowStart() const;

    size_type rowStop() const;

    graphmpi_type const& graph() const
    {
        return M_graph;
    }

    void setGraph( graphmpi_type const& graph )
    {
        M_graph = graph;
    }


private:

    Mat M_mat;
    Communicator M_comm;
    bool M_is_initialized;
    graphmpi_type M_graph;

};

} // Nextsim
#endif // __MatrixPetsc_H
