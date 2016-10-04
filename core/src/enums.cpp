/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

/**
 * @file   enums.cpp
 * @author Abdoulaye Samake <abdoulaye.samake.nersc.no>
 * @date   Mon Sep 14 13:10:44 2015
 */

#include <enums.hpp>

namespace Nextsim
{
SolverType
kspTypeConvertStrToEnum( std::string const& type )
{
    /**/ if ( type=="cg" )         return SolverType::CG;
    else if ( type=="cr" )         return SolverType::CR;
    else if ( type=="cgs" )        return SolverType::CGS;
    else if ( type=="bicg" )       return SolverType::BICG;
    else if ( type=="tcqmr" )      return SolverType::TCQMR;
    else if ( type=="tfqmr" )      return SolverType::TFQMR;
    else if ( type=="lsqr" )       return SolverType::LSQR;
    else if ( type=="bicgstab" )   return SolverType::BICGSTAB;
    else if ( type=="minres" )     return SolverType::MINRES;
    else if ( type=="gmres" )      return SolverType::GMRES;
    else if ( type=="fgmres" )     return SolverType::FGMRES;
    else if ( type=="richardson" ) return SolverType::RICHARDSON;
    else if ( type=="chebyshev" )  return SolverType::CHEBYSHEV;
    else if ( type=="preonly" )    return SolverType::PREONLY;
    else if ( type=="gcr" )        return SolverType::GCR;
    else                           return SolverType::GMRES;
}

PreconditionerType
pcTypeConvertStrToEnum( std::string const& type )
{
	/**/ if ( type=="lu" )           return PreconditionerType::LU_PRECOND;
	else if ( type=="ilu" )          return PreconditionerType::ILU_PRECOND;
	else if ( type=="id" )           return PreconditionerType::IDENTITY_PRECOND;
    else if ( type=="cholesky" )     return PreconditionerType::CHOLESKY_PRECOND;
    else if ( type=="icc" )          return PreconditionerType::ICC_PRECOND;
    else if ( type=="asm" )          return PreconditionerType::ASM_PRECOND;
    else if ( type=="gasm" )         return PreconditionerType::GASM_PRECOND;
    else if ( type=="jacobi" )       return PreconditionerType::JACOBI_PRECOND;
    else if ( type=="block_jacobi" ) return PreconditionerType::BLOCK_JACOBI_PRECOND;
    else if ( type=="bjacobi" )      return PreconditionerType::BLOCK_JACOBI_PRECOND;
    else if ( type=="sor" )          return PreconditionerType::SOR_PRECOND;
    else if ( type=="eisenstat" )    return PreconditionerType::EISENSTAT_PRECOND;
    else if ( type=="shell" )        return PreconditionerType::SHELL_PRECOND;
    else if ( type=="fieldsplit" )   return PreconditionerType::FIELDSPLIT_PRECOND;
    else if ( type=="lsc" )          return PreconditionerType::LSC_PRECOND;
    else if ( type=="lsc2" )         return PreconditionerType::LSC2_PRECOND;
    else if ( type=="blockns" )      return PreconditionerType::FEELPP_BLOCKNS_PRECOND;
    else if ( type=="blockms" )      return PreconditionerType::FEELPP_BLOCKMS_PRECOND;
    else if ( type=="ml" )           return PreconditionerType::ML_PRECOND;
    else if ( type=="gamg" )         return PreconditionerType::GAMG_PRECOND;
    else if ( type=="boomeramg" )    return PreconditionerType::BOOMERAMG_PRECOND;
    else if ( type=="redundant" )    return PreconditionerType::REDUNDANT_PRECOND;
    else if ( type=="none" )         return PreconditionerType::NONE_PRECOND;
    else                             return PreconditionerType::LU_PRECOND;
}

MatSolverPackageType
matSolverPackageConvertStrToEnum( std::string const& type )
{
    /**/ if ( type =="spooles" )     return MatSolverPackageType::MATSOLVER_SPOOLES;
    else if ( type=="superlu" )      return MatSolverPackageType::MATSOLVER_SUPERLU;
    else if ( type=="superlu-dist" ) return MatSolverPackageType::MATSOLVER_SUPERLU_DIST;
    else if ( type=="umfpack" )      return MatSolverPackageType::MATSOLVER_UMFPACK;
    else if ( type=="cholmod" )      return MatSolverPackageType::MATSOLVER_CHOLMOD;
    else if ( type=="essl" )         return MatSolverPackageType::MATSOLVER_ESSL;
    else if ( type=="lusol" )        return MatSolverPackageType::MATSOLVER_LUSOL;
    else if ( type=="mumps" )        return MatSolverPackageType::MATSOLVER_MUMPS;
    else if ( type=="mkl_pardiso" )  return MatSolverPackageType::MATSOLVER_MKL_PARDISO;
    else if ( type=="pastix" )       return MatSolverPackageType::MATSOLVER_PASTIX;
    else if ( type=="dscpack" )      return MatSolverPackageType::MATSOLVER_DSCPACK;
    else if ( type=="matlab" )       return MatSolverPackageType::MATSOLVER_MATLAB;
    else if ( type=="petsc" )        return MatSolverPackageType::MATSOLVER_PETSC;
    else if ( type=="plapack" )      return MatSolverPackageType::MATSOLVER_PLAPACK;
    else if ( type=="bas" )          return MatSolverPackageType::MATSOLVER_BAS;
    else if ( type=="boomeramg" )    return MatSolverPackageType::MATSOLVER_BOOMERAMG;
    else if ( type=="euclid" )       return MatSolverPackageType::MATSOLVER_EUCLID;
    else if ( type=="pilut" )        return MatSolverPackageType::MATSOLVER_PILUT;
    else                             return MatSolverPackageType::MATSOLVER_PETSC;
}


} // Nextsim
