/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   enums.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Tue Jul  7 13:08:30 2015
 */

#include <string>

namespace Nextsim
{

enum SolverType
{
     CG=0,
     CGN,
     CGS,
     CR,
     QMR,
     TCQMR,
     TFQMR,
     BICG,
     BICGSTAB,
     MINRES,
     GMRES,
     FGMRES,
     LSQR,
     JACOBI,
     SOR_FORWARD,
     SOR_BACKWARD,
     SSOR,
     RICHARDSON,
     CHEBYSHEV,
     PREONLY,
     GCR,
     INVALID_SOLVER
};

enum PreconditionerType
{
	 IDENTITY_PRECOND=0,
     JACOBI_PRECOND,
     BLOCK_JACOBI_PRECOND,
     SOR_PRECOND,
     SSOR_PRECOND,
     EISENSTAT_PRECOND,
     ASM_PRECOND,
     GASM_PRECOND,
     CHOLESKY_PRECOND,
     ICC_PRECOND,
     ILU_PRECOND,
     LU_PRECOND,
     AMG_PRECOND,
     USER_PRECOND,
     SHELL_PRECOND,
     FIELDSPLIT_PRECOND,
     LSC_PRECOND,
     LSC2_PRECOND,
     FEELPP_BLOCKNS_PRECOND,
     FEELPP_BLOCKMS_PRECOND,
     ML_PRECOND,
     GAMG_PRECOND,
     BOOMERAMG_PRECOND,
     REDUNDANT_PRECOND,
     NONE_PRECOND,
     INVALID_PRECONDITIONER
};

enum MatSolverPackageType
{
	 MATSOLVER_SPOOLES=0,
	 MATSOLVER_SUPERLU,
	 MATSOLVER_SUPERLU_DIST,
	 MATSOLVER_UMFPACK,
	 MATSOLVER_ESSL,
	 MATSOLVER_LUSOL,
	 MATSOLVER_MUMPS,
	 MATSOLVER_MKL_PARDISO,
	 MATSOLVER_PASTIX,
	 MATSOLVER_DSCPACK,
	 MATSOLVER_MATLAB,
	 MATSOLVER_PETSC,
	 MATSOLVER_PLAPACK,
	 MATSOLVER_BAS,
	 MATSOLVER_BOOMERAMG,
	 MATSOLVER_EUCLID,
	 MATSOLVER_PILUT,
     MATSOLVER_CHOLMOD
};

SolverType
kspTypeConvertStrToEnum( std::string const& type );

PreconditionerType
pcTypeConvertStrToEnum( std::string const& type );

MatSolverPackageType
matSolverPackageConvertStrToEnum( std::string const& type );

} // Nextsim
