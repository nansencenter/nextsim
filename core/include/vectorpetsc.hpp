/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   vectorpetsc.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Jul  6 11:43:08 2015
 */

#ifndef __VectorPetsc_H
#define __VectorPetsc_H 1

#include <environment.hpp>
#include <matrixpetsc.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <assert.hpp>
#include <petsc.hpp>

extern "C"
{
#include <petscmat.h>
}

/**
 * @class VectorPetsc
 * @brief Wrapper for petsc vectors
 *
 * @see
 *
 */

namespace Nextsim
{

//class MatrixPetsc;

class VectorPetsc
{

public:

	typedef std::size_t size_type;
    typedef double value_type;

    VectorPetsc( Communicator const& comm = Environment::comm() );

	~VectorPetsc();

	Communicator const& comm() const { return M_comm; }

	Vec vec () const;

    Vec& vec ();

    void init( const size_type n, bool fast = false );

    void zero();

    void resize(const size_type n, bool fast = false);

	void close();

	void set(const value_type& value);

	void set(size_type i, const value_type& value);

	void setConstant(value_type value);

	void setOnes();

	void setVector(int* i, int n, value_type* v);

	void add(size_type i, const value_type& value);

	void addVector(int* i, int n, value_type* v);

    void addVector(VectorPetsc const& v_in, MatrixPetsc const& A_in);

	void add(const value_type& v_in);

	void add(const VectorPetsc& v);

	void add(const value_type& a_in, const VectorPetsc& v_in);

	void scale(value_type factor_in);

	value_type l1Norm() const;

	value_type l2Norm() const;

	value_type linftyNorm() const;

	value_type min() const;

	value_type max() const;

	value_type sum() const;

	value_type operator() (const size_type i) const;

	value_type& operator() (const size_type i);

    VectorPetsc& operator += (const VectorPetsc& v_in);

    VectorPetsc& operator -= (const VectorPetsc& v_in);

	size_type size() const;

    void printScreen() const;

	void printMatlab(std::string const& filename = "NULL") const;

    void printBinary(std::string const& filename = "NULL") const;

    void loadBinary(std::string const& filename = "NULL");

	void clear();

    std::vector<value_type> container();
    
    void container(double *values);

	bool isInitialized() const { return M_is_initialized; }

	bool closed() const { return M_is_closed; }


private:

	Vec M_vec;
	Communicator M_comm;
	bool M_is_initialized;
	bool M_is_closed;
};

} // Nextsim
#endif // __VectorPetsc_H
