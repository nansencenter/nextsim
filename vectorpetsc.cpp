/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   vectorpetsc.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Jul  6 15:56:26 2015
 */

#include <vectorpetsc.hpp>

/**
 * @class VectorPetsc
 * @brief Wrapper for petsc vectors
 *
 * @see
 *
 */

namespace Nextsim
{
VectorPetsc::VectorPetsc( const size_type n, Communicator const& comm )
	:
	M_comm( comm )
{
	int ierr=0;
	int petsc_n=static_cast<int>( n );

	//ierr = VecCreateSeq ( PETSC_COMM_SELF, petsc_n, &M_vec );
	ierr = VecCreateSeq ( M_comm, petsc_n, &M_vec );
	CHKERRABORT( M_comm,ierr );

	ierr = VecSetFromOptions ( M_vec );
	CHKERRABORT( M_comm,ierr );

	M_is_initialized = true;
}

VectorPetsc::~VectorPetsc()
{
	this->clear();
}

Vec
VectorPetsc::vec() const
{
	ASSERT(M_vec != 0, "invalid petsc vector");
	return M_vec;
}

Vec&
VectorPetsc::vec()
{
	ASSERT(M_vec != 0, "invalid petsc vector");
	return M_vec;
}

void
VectorPetsc::zero()
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");

	int ierr=0;

	PetscScalar z=0.;
	this->close();

	ierr = VecSet ( M_vec, z );
	CHKERRABORT( M_comm,ierr );
}

void
VectorPetsc::close()
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");

	int ierr=0;

	ierr = VecAssemblyBegin( M_vec );
	CHKERRABORT( M_comm,ierr );
	ierr = VecAssemblyEnd( M_vec );
	CHKERRABORT( M_comm,ierr );

	M_is_closed = true;
}

void
VectorPetsc::set(const value_type& value)
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");

	int ierr=0;
	PetscScalar petsc_value = static_cast<PetscScalar>( value );

	ierr = VecSet ( M_vec, petsc_value );
	CHKERRABORT( M_comm,ierr );
}

void
VectorPetsc::set(size_type i, const value_type& value)
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");
	ASSERT(i < this->size(), "invalid index");

	int ierr=0;
	int i_val = static_cast<int>( i );
	PetscScalar petsc_value = static_cast<PetscScalar>( value );

	ierr = VecSetValues ( M_vec, 1, &i_val, &petsc_value, INSERT_VALUES );
	CHKERRABORT( M_comm,ierr );
}

void
VectorPetsc::setVector(int* i, int n, value_type* v)
{
#if PETSC_VERSION_LESS_THAN(3,5,3)
	if ( n == 0 ) return;
#endif
	int ierr=0;
	ierr = VecSetValues ( M_vec, n, i, v, INSERT_VALUES );
	CHKERRABORT( M_comm,ierr );
}

void
VectorPetsc::setConstant(value_type value)
{
	this->set(value);
}

void
VectorPetsc::setOnes()
{
	this->setConstant(1.);
}

void
VectorPetsc::add(size_type i, const value_type& value)
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");
	ASSERT(i < this->size(), "invalid index");

	int ierr=0;
	int i_val = static_cast<int>( i );
	PetscScalar petsc_value = static_cast<PetscScalar>( value );

	ierr = VecSetValues ( M_vec, 1, &i_val, &petsc_value, ADD_VALUES );
	CHKERRABORT( M_comm,ierr );
}

void
VectorPetsc::addVector(int* i, int n, value_type* v)
{
#if PETSC_VERSION_LESS_THAN(3,5,3)
	if ( n == 0 ) return;
#endif
	int ierr=0;
	ierr = VecSetValues ( M_vec, n, i, v, ADD_VALUES );
	CHKERRABORT( M_comm,ierr );
}

void
VectorPetsc::addVector(VectorPetsc const& v_in, MatrixPetsc const& A_in)
{
    this->close();
    A_in.close();

    int ierr = 0;

    ierr = MatMultAdd( A_in.mat(), v_in.vec(), M_vec, M_vec );
    CHKERRABORT( M_comm,ierr );

}

void
VectorPetsc::add(const value_type& v_in)
{
	int ierr=0;
	PetscScalar* values;
	const PetscScalar v = static_cast<PetscScalar>( v_in );
	const int n = static_cast<int>( this->size() );
	//const int fli = static_cast<int>( this->firstLocalIndex() );

	for ( int i = 0; i < n; i++ )
	{
		ierr = VecGetArray ( M_vec, &values );
		CHKERRABORT( M_comm,ierr );

		//int ig = fli + i;

		PetscScalar value = ( values[i] + v );

		ierr = VecRestoreArray ( M_vec, &values );
		CHKERRABORT( M_comm,ierr );

		ierr = VecSetValues ( M_vec, 1, &i, &value, INSERT_VALUES );
		CHKERRABORT( M_comm,ierr );
	}
 }

void
VectorPetsc::add(const VectorPetsc& v)
{
	this->add(1., v);
}

void
VectorPetsc::add(const value_type& a_in, const VectorPetsc& v_in)
{
	int ierr = 0;
	PetscScalar a = static_cast<PetscScalar>( a_in );

	ASSERT(this->size() == v_in.size(), "invalid vector sizes");

	// 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

	ierr = VecAXPY( &a, v->M_vec, M_vec );
	CHKERRABORT( M_comm,ierr );

	// 2.3.x & later style
#else

	ierr = VecAXPY( M_vec, a, v_in.M_vec );
	CHKERRABORT( M_comm,ierr );

#endif
}

void
VectorPetsc::scale(value_type factor_in)
{
	int ierr = 0;
	PetscScalar factor = static_cast<PetscScalar>( factor_in );

	// 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

	ierr = VecScale( &factor, M_vec );
	CHKERRABORT( M_comm,ierr );

	// 2.3.x & later style
#else

	ierr = VecScale( M_vec, factor );
	CHKERRABORT( M_comm,ierr );

#endif
}

typename VectorPetsc::value_type
VectorPetsc::l1Norm() const
{
	if ( !this->closed() )
		const_cast<VectorPetsc*>( this )->close();

	int ierr=0;
	double value=0.;

	ierr = VecNorm ( M_vec, NORM_1, &value );
	CHKERRABORT( M_comm,ierr );

	//return static_cast<value_type>( value );
	return value;
}

typename VectorPetsc::value_type
VectorPetsc::l2Norm() const
{
	if ( !this->closed() )
		const_cast<VectorPetsc*>( this )->close();

	int ierr=0;
	double value=0.;

	ierr = VecNorm ( M_vec, NORM_2, &value );
	CHKERRABORT( M_comm,ierr );

	//return static_cast<Real>( value );
	return value;
}

typename VectorPetsc::value_type
VectorPetsc::linftyNorm() const
{
	if ( !this->closed() )
		const_cast<VectorPetsc*>( this )->close();

	int ierr=0;
	double value=0.;

	ierr = VecNorm ( M_vec, NORM_INFINITY, &value );
	CHKERRABORT( M_comm,ierr );

	//return static_cast<Real>( value );
	return value;
}

typename VectorPetsc::value_type
VectorPetsc::min() const
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");

	int index=0, ierr=0;
	PetscReal min=0.;

	ierr = VecMin ( M_vec, &index, &min );
	CHKERRABORT( M_comm,ierr );

	// this return value is correct: VecMin returns a PetscReal
	return static_cast<value_type>( min );
}

typename VectorPetsc::value_type
VectorPetsc::max() const
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");

	int index=0, ierr=0;
	PetscReal max=0.;

	ierr = VecMax ( M_vec, &index, &max );
	CHKERRABORT( M_comm,ierr );

	// this return value is correct: VecMax returns a PetscReal
	return static_cast<value_type>( max );
}

typename VectorPetsc::value_type
VectorPetsc::sum() const
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");

	if ( !this->closed() )
		const_cast<VectorPetsc*>( this )->close();

	int ierr=0;
	double value=0.;

	ierr = VecSum ( M_vec, &value );
	CHKERRABORT( M_comm,ierr );

	//return static_cast<value_type>( value );
	return value;
}

typename VectorPetsc::value_type
VectorPetsc::operator () (const size_type i) const
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");
	ASSERT(i < this->size(), "invalid index");

	int ierr=0;
	PetscScalar *values, value=0.;

	ierr = VecGetArray( M_vec, &values );
	CHKERRABORT( M_comm,ierr );

	value = values[i];

	ierr = VecRestoreArray ( M_vec, &values );
	CHKERRABORT( M_comm,ierr );

	return static_cast<value_type>( value );
}

typename VectorPetsc::value_type&
VectorPetsc::operator () (const size_type i)
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");
	ASSERT(i < this->size(), "invalid index");

	int ierr=0;
	PetscScalar *values;

	ierr = VecGetArray( M_vec, &values );
	CHKERRABORT( M_comm,ierr );

	PetscScalar& value = values[i];

	ierr = VecRestoreArray ( M_vec, &values );
	CHKERRABORT( M_comm,ierr );

	return static_cast<value_type&>( value );
}

VectorPetsc&
VectorPetsc::operator += (const VectorPetsc& v_in)
{
    ASSERT(M_is_closed, "VectorPetsc is not closed");

    this->add( 1., v_in );

    return *this;
}

VectorPetsc&
VectorPetsc::operator -= (const VectorPetsc& v_in)
{
    ASSERT(M_is_closed, "VectorPetsc is not closed");

    this->add( -1., v_in );

    return *this;
}

typename VectorPetsc::size_type
VectorPetsc::size() const
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");

	if ( !M_is_initialized )
		return 0;

	int petsc_size=0;
	int ierr = VecGetSize( M_vec, &petsc_size );
	CHKERRABORT( M_comm,ierr );
	return static_cast<size_type>( petsc_size );
}

void
VectorPetsc::printMatlab(std::string const& filename) const
{
	ASSERT(M_is_initialized, "VectorPetsc not initialized");

	if ( !this->closed() )
		const_cast<VectorPetsc*>( this )->close();

	PetscObjectSetName((PetscObject)M_vec,boost::filesystem::path("out_"+filename).stem().string().c_str());

	int ierr=0;

	PetscViewer petsc_viewer;

	ierr = PetscViewerCreate ( M_comm, &petsc_viewer );
	CHKERRABORT( M_comm,ierr );

	if ( filename != "NULL" )
	{
		ierr = PetscViewerASCIIOpen( M_comm,
		                             filename.c_str(),
		                             &petsc_viewer );
		CHKERRABORT( M_comm,ierr );

		ierr = PetscViewerSetFormat ( petsc_viewer,
		                              PETSC_VIEWER_ASCII_MATLAB );
		CHKERRABORT( M_comm,ierr );

		ierr = VecView ( const_cast<Vec>( M_vec ), petsc_viewer );
		CHKERRABORT( M_comm,ierr );
	}
	else
	{
		ierr = PetscViewerSetFormat ( PETSC_VIEWER_STDOUT_WORLD,
		                              PETSC_VIEWER_ASCII_MATLAB );
		CHKERRABORT( M_comm,ierr );

		ierr = VecView ( const_cast<Vec>( M_vec ), PETSC_VIEWER_STDOUT_WORLD );
		CHKERRABORT( M_comm,ierr );
	}

	ierr = PetscViewerDestroy ( &petsc_viewer );
	CHKERRABORT( M_comm,ierr );
}

void
VectorPetsc::clear()
{
	int ierr=0;

	if (M_is_initialized)
	{
		ierr = VecDestroy(&M_vec);
		CHKERRABORT( M_comm,ierr );
	}

	M_is_initialized = true;
	M_is_closed = true;
}

} // Nextsim
