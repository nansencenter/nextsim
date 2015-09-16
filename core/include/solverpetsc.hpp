/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   solverpetsc.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Tue Jul  7 12:11:40 2015
 */

#ifndef __SolverPetsc_H
#define __SolverPetsc_H 1

#include <matrixpetsc.hpp>
#include <vectorpetsc.hpp>
#include <parameter.hpp>
#include <enums.hpp>
#include <assert.hpp>
#include <petsc.hpp>

extern "C" {
# include <petscversion.h>
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
#   include <petscsles.h>
# else
#   include <petscksp.h>
# endif
}

namespace Nextsim
{
class SolverPetsc
{

public:

	typedef std::size_t size_type;
    typedef double value_type;

    typedef MatrixPetsc matrix_type;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;
    typedef VectorPetsc vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;


	SolverPetsc( Communicator const& comm = Environment::comm() );

	~SolverPetsc();

	Communicator const& comm() const { return M_comm; }

    void init();

	void solveLinearSystem(matrix_ptrtype const& matrix,
                           vector_ptrtype& solution,
                           vector_ptrtype const& rhs);

    void setSolverType(std::string const& st);
    void setPreconditionerType(std::string const& pct);
    void setMatSolverPackageType(std::string const& mspackt);
    void PetscPCFactorSetMatSolverPackage(PC & pc, MatSolverPackageType mspackt);
    void setReusePrec(bool reuse);

	void clear();

	bool initialized() const { return M_is_initialized; }
	value_type rTolerance() const { return M_rtolerance; }
	value_type dTolerance() const { return M_dtolerance; }
	value_type aTolerance() const { return M_atolerance; }
	size_type maxIterations() const { return M_maxit; }
	SolverType solverType() const { return M_solver_type; }
    PreconditionerType preconditionerType() const { return M_preconditioner_type; }
    MatSolverPackageType matSolverPackageType() const { return M_matSolverPackage_type; }
    bool reusePrec() const { return M_reuse_prec;}

    std::string PetscConvertKSPReasonToString(KSPConvergedReason reason);
    size_type nIterations() const { return M_iteration; }
    value_type residual() const { return M_residual; }
    std::string convergedReason() const { return M_reason; }

	KSP ksp();

	PC pc();

    BOOST_PARAMETER_MEMBER_FUNCTION( ( void ),
                                     solve,
                                     tag,
                                     ( required
                                       ( matrix, (matrix_ptrtype) )
                                       ( solution, ( vector_ptrtype ) )
                                       ( rhs, ( vector_ptrtype ) )
                                       )
                                     ( optional
                                       ( ksp,( std::string ),M_ksp_type/*"preonly"*/ )
                                       ( pc,( std::string ),M_pc_type/*"cholesky"*/ )
                                       ( pcfactormatsolverpackage,( std::string ), M_pcfactormatsolverpackage_type/*"cholmod"*/ )
                                       ( rtolerance,( double ), M_rtolerance/*1e-13*/ )
                                       ( atolerance,( double ), M_atolerance/*1e-50*/ )
                                       ( dtolerance,( double ), M_dtolerance/*1e5*/ )
                                       ( maxit,( size_type ), M_maxit/*1000*/ )
                                       ( reuse_prec,( bool ), M_reuse_prec )
                                       ) )
    {
        this->setSolverType(ksp);
        this->setPreconditionerType(pc);
        this->setMatSolverPackageType(pcfactormatsolverpackage);
        this->setReusePrec(reuse_prec);

        setTolerances(_rtolerance=rtolerance,
                      _dtolerance=dtolerance,
                      _atolerance=atolerance,
                      _maxit=maxit);

        this->solveLinearSystem(matrix, solution, rhs);;
    }



	BOOST_PARAMETER_MEMBER_FUNCTION( ( void ),
	                                 setTolerances,
	                                 tag,
	                                 ( required
	                                   ( rtolerance,( double ) )
	                                   )
	                                 ( optional
	                                   ( maxit,( size_type ), 1000 )
	                                   ( atolerance,( double ), 1e-50 )
	                                   ( dtolerance,( double ), 1e5 )
	                                   ) )
	{
		M_rtolerance = rtolerance;
		M_dtolerance = dtolerance;
		M_atolerance = atolerance;
		M_maxit = maxit;
	}



private:

	Communicator M_comm;
	bool M_is_initialized;

	// Enum stating which type of iterative solver to use
	SolverType M_solver_type;

	// Enum statitng with type of preconditioner to use
	PreconditionerType M_preconditioner_type;

    // Enum the software that is used to perform the factorization
    MatSolverPackageType M_matSolverPackage_type;

	// Preconditioner context
	PC M_pc;

	// Krylov subspace context
	KSP M_ksp;

	// relative tolerance
	value_type M_rtolerance;

	// divergence tolerance
	value_type M_dtolerance;

	// absolute tolerance
	value_type M_atolerance;

	// maximum number of iterations
	size_type M_maxit;

    // number of iterations required for ksp solver
    size_type M_iteration;

    // norm of the final residual
    value_type M_residual;

    // reason ksp iteration was stopped
    std::string M_reason;

    std::string M_ksp_type;
    std::string M_pc_type;
    std::string M_pcfactormatsolverpackage_type;
    bool M_reuse_prec;

private:

	void setPetscSolverType();
	void setPetscPreconditionerType();
    void setInitialized(bool init);
};

} // Nextsim

#endif // __SolverPetsc_H
