/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   solverpetsc.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Tue Jul  7 15:20:55 2015
 */

#include <solverpetsc.hpp>

namespace Nextsim
{
SolverPetsc::SolverPetsc(Communicator const& comm)
	:
	M_comm(comm),
	M_solver_type(PREONLY),
	M_preconditioner_type(LU_PRECOND),
    M_matSolverPackage_type(MATSOLVER_UMFPACK),
    M_rtolerance(1e-13),
    M_dtolerance(1e5),
    M_atolerance(1e-50),
    M_maxit(1000),
	M_is_initialized(false)
{}

SolverPetsc::~SolverPetsc()
{
	this->clear();
}

void
SolverPetsc::init()
{
    if ( !this->initialized() )
    {
        M_is_initialized = true;

        int ierr=0;

        // Create the linear solver context
        ierr = KSPCreate ( M_comm, &M_ksp );
        CHKERRABORT( M_comm,ierr );

        // Create the preconditioner context
        ierr = KSPGetPC ( M_ksp, &M_pc );
        CHKERRABORT( M_comm,ierr );

        // Have the Krylov subspace method use our good initial guess rather than 0
        bool useInitialGuessNonZero = false; //boption(_name="ksp-use-initial-guess-nonzero", _prefix=this->prefix() );
        ierr = KSPSetInitialGuessNonzero ( M_ksp, (useInitialGuessNonZero)?PETSC_TRUE:PETSC_FALSE );
        CHKERRABORT( M_comm,ierr );

        // Set user-specified  solver and preconditioner types
        this->setPetscSolverType();
        //this->setPetscConstantNullSpace();

        ierr = KSPSetFromOptions ( M_ksp );
        CHKERRABORT( M_comm,ierr );

#if PETSC_VERSION_LESS_THAN(3,0,0)
        KSPType ksp_type;
#else
#if PETSC_VERSION_LESS_THAN(3,4,0)
        const KSPType ksp_type;
#else
        KSPType ksp_type;
#endif
#endif

        ierr = KSPGetType ( M_ksp, &ksp_type );
        CHKERRABORT( M_comm,ierr );


        ierr = KSPSetResidualHistory( M_ksp,
                                      PETSC_NULL,   // pointer to the array which holds the history
                                      PETSC_DECIDE, // size of the array holding the history
                                      PETSC_TRUE ); // Whether or not to reset the history for each solve.
        CHKERRABORT( M_comm,ierr );

        this->setPetscPreconditionerType();
        PetscPCFactorSetMatSolverPackage( M_pc,this->matSolverPackageType() );

        KSPMonitorSet( M_ksp,KSPMonitorDefault,PETSC_NULL,PETSC_NULL );
    }
}

void
SolverPetsc::solve(MatrixPetsc const& matrix,
                   VectorPetsc& solution,
                   VectorPetsc const& rhs,
                   const value_type tolerance,
                   const size_type maxit)
{
    this->init();

    matrix.close();
    solution.close();
    //rhs.close();

    int ierr=0;
    int its=0;
    PetscReal final_resid=0.;

#if PETSC_VERSION_LESS_THAN(3,5,0)
    ierr = KSPSetOperators( M_ksp, matrix.mat(), matrix.mat(), DIFFERENT_NONZERO_PATTERN);
#else
    bool same_preconditioner=true;
    ierr = KSPSetReusePreconditioner( M_ksp, (same_preconditioner)? PETSC_TRUE : PETSC_FALSE );
    CHKERRABORT( M_comm,ierr );
    ierr = KSPSetOperators( M_ksp, matrix.mat(), matrix.mat() );
#endif
    CHKERRABORT( M_comm,ierr );

    ierr = KSPSetTolerances ( M_ksp,
                              this->rTolerance(),
                              //1e-15,
                              this->aTolerance(),
                              this->dTolerance(),
                              this->maxIterations() );
    CHKERRABORT( M_comm,ierr );

#if PETSC_VERSION_LESS_THAN(3,5,0)
    KSPDefaultConvergedSetUIRNorm( M_ksp );
#else
    KSPConvergedDefaultSetUIRNorm( M_ksp );
#endif

    // Solve the linear system
    ierr = KSPSolve ( M_ksp, rhs.vec(), solution.vec() );
    CHKERRABORT( M_comm,ierr );

    // Get the number of iterations required for convergence
    ierr = KSPGetIterationNumber ( M_ksp, &its );
    M_iteration = its;

    std::cout << "NITER = " << M_iteration << "\n";
    CHKERRABORT( M_comm,ierr );

    // Get the norm of the final residual to return to the user.
    //ierr = KSPGetResidualNorm ( M_ksp, &final_resid );
    ierr = KSPGetResidualNorm ( M_ksp, &M_residual );
    std::cout << "RESIDUAL = " << M_residual << "\n";
    CHKERRABORT( M_comm,ierr );

    KSPConvergedReason reason;
    KSPGetConvergedReason( M_ksp,&reason );
    M_reason = PetscConvertKSPReasonToString(reason);

    ierr = KSPView( M_ksp, PETSC_VIEWER_STDOUT_WORLD );

    std::cout<< "[solverpetsc] reason = " << PetscConvertKSPReasonToString(reason) <<"\n";

    if (reason == KSP_DIVERGED_INDEFINITE_PC)
    {
        std::cout << "[solverpetsc] Divergence because of indefinite preconditioner \n";
        std::cout << "[solverpetsc] Run the executable again but with '-pc_factor_shift_type POSITIVE_DEFINITE' option \n";
    }
    else if (reason < 0)
    {
        std::cout <<"[solverpetsc] Other kind of divergence: this should not happen \n";
    }
}

void
SolverPetsc::setSolverType(const SolverType st)
{
	M_solver_type = st;
}

void SolverPetsc::setPreconditionerType(const PreconditionerType pct)
{
    M_preconditioner_type = pct;
}

void
SolverPetsc::setMatSolverPackageType(const MatSolverPackageType mspackt)
{
    M_matSolverPackage_type = mspackt;
}

void
SolverPetsc::clear()
{
    PetscBool pinit;
    PetscInitialized( &pinit );

    if ( pinit && this->initialized() )
    {
        this->setInitialized( false );

        int ierr=0;

        ASSERT(M_ksp != 0, "invalid ksp");
        ierr = KSPDestroy( &M_ksp );
        CHKERRABORT( M_comm,ierr );
    }
}

void
SolverPetsc::setInitialized( bool init )
{
    M_is_initialized = init;
}


void
SolverPetsc::setPetscSolverType()
{
	int ierr = 0;

	switch ( this->solverType() )
	{
	case CG:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPCG );
		CHKERRABORT( M_comm,ierr );
		return;

	case CR:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPCR );
		CHKERRABORT( M_comm,ierr );
		return;

	case CGS:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPCGS );
		CHKERRABORT( M_comm,ierr );
		return;

	case BICG:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPBICG );
		CHKERRABORT( M_comm,ierr );
		return;

	case TCQMR:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPTCQMR );
		CHKERRABORT( M_comm,ierr );
		return;

	case TFQMR:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPTFQMR );
		CHKERRABORT( M_comm,ierr );
		return;

	case LSQR:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPLSQR );
		CHKERRABORT( M_comm,ierr );
		return;

	case BICGSTAB:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPBCGS );
		CHKERRABORT( M_comm,ierr );
		return;

	case MINRES:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPMINRES );
		CHKERRABORT( M_comm,ierr );
		return;

	case GMRES:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPGMRES );
		CHKERRABORT( M_comm,ierr );
		return;

	case FGMRES:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPFGMRES );
		CHKERRABORT( M_comm,ierr );
		return;

	case RICHARDSON:
		ierr = KSPSetType ( M_ksp, ( char* ) KSPRICHARDSON );
		CHKERRABORT( M_comm,ierr );
		return;
	case PREONLY :
		ierr = KSPSetType ( M_ksp, ( char* ) KSPPREONLY );
		CHKERRABORT( M_comm,ierr );
		return;

	case GCR :
		ierr = KSPSetType ( M_ksp, ( char* ) KSPGCR );
		CHKERRABORT( M_comm,ierr );
		return;

	default:
		std::cerr << "ERROR:  Unsupported PETSC Solver: "
		          << this->solverType()               << std::endl
		          << "Continuing with PETSC defaults" << std::endl;

	}
}

void
SolverPetsc::setPetscPreconditionerType()
{
	int ierr = 0;
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
	ierr = PCFactorSetMatSolverPackage( M_pc,MATSOLVERUMFPACK );

	if ( ierr )
	{
		ierr = PCFactorSetMatSolverPackage( M_pc,MATSOLVERSUPERLU );

		if ( ierr )
		{
			ierr = PCFactorSetMatSolverPackage( M_pc,MATSOLVERPETSC );
		}
	}

#elif (PETSC_VERSION_MAJOR >= 3)
	ierr = PCFactorSetMatSolverPackage( M_pc,MAT_SOLVER_UMFPACK );

	if ( ierr )
	{
		ierr = PCFactorSetMatSolverPackage( M_pc,MAT_SOLVER_SUPERLU );

		if ( ierr )
		{
			ierr = PCFactorSetMatSolverPackage( M_pc,MAT_SOLVER_PETSC );
		}
	}

#endif

    switch ( this->preconditionerType() )
    {
    case IDENTITY_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCNONE );
        CHKERRABORT( M_comm,ierr );
        return;

    case CHOLESKY_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCCHOLESKY );
        CHKERRABORT( M_comm,ierr );
        return;

    case ICC_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCICC );
        CHKERRABORT( M_comm,ierr );
        return;

    case ILU_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCILU );
        CHKERRABORT( M_comm,ierr );

    case LU_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCLU );
        CHKERRABORT( M_comm,ierr );
        return;

    case ASM_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCASM );
        CHKERRABORT( M_comm,ierr );
        return;

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
    case GASM_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCGASM );
        CHKERRABORT( M_comm,ierr );
        return;
#endif

    case JACOBI_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCJACOBI );
        CHKERRABORT( M_comm,ierr );
        return;

    case BLOCK_JACOBI_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCBJACOBI );
        CHKERRABORT( M_comm,ierr );
        return;

    case SOR_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCSOR );
        CHKERRABORT( M_comm,ierr );
        return;

    case EISENSTAT_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCEISENSTAT );
        CHKERRABORT( M_comm,ierr );
        return;

    default:
        std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
                  << this->preconditionerType()       << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }
}

void
SolverPetsc::PetscPCFactorSetMatSolverPackage(PC & pc, MatSolverPackageType mspackt)
{
    int ierr = 0;

    switch ( mspackt )
    {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
#if PETSC_VERSION_LESS_THAN(3,4,0)
    case MATSOLVER_SPOOLES :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERSPOOLES );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;
#endif

    case MATSOLVER_SUPERLU :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERSUPERLU );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_SUPERLU_DIST :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERSUPERLU_DIST );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_UMFPACK :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERUMFPACK );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_ESSL :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERESSL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_LUSOL :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERLUSOL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_MUMPS :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERMUMPS );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,5,0)
    case MATSOLVER_MKL_PARDISO :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERMKL_PARDISO );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;
#endif
    case MATSOLVER_PASTIX :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERPASTIX );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_MATLAB :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERMATLAB );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_PETSC :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERPETSC );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

#if PETSC_VERSION_LESS_THAN(3,4,0)
    case MATSOLVER_PLAPACK :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERPLAPACK );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;
#endif
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 1)

    case MATSOLVER_BAS :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERBAS );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;
#endif

#else // PETSC < 3.2

    case MATSOLVER_SPOOLES :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_SPOOLES );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_SUPERLU :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_SUPERLU );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_SUPERLU_DIST :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_SUPERLU_DIST );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_UMFPACK :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_UMFPACK );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_ESSL :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_ESSL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_LUSOL :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_LUSOL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_MUMPS :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_MUMPS );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_PASTIX :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_PASTIX );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_DSCPACK :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_DSCPACK );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_MATLAB :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_MATLAB );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_PETSC :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_PETSC );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case MATSOLVER_PLAPACK :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_PLAPACK );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 1)

    case MATSOLVER_BAS :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_BAS );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;
#endif
#endif

    default:
        std::cerr << "ERROR:  Unsupported PETSC mat solver package: "
                  << mspackt               << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }

    const MatSolverPackage t;
    ierr = PCFactorGetMatSolverPackage( pc, &t );
}


std::string
SolverPetsc::PetscConvertKSPReasonToString(KSPConvergedReason reason)
{
    switch ( reason )
    {
        /* converged */
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3, 2, 0 )
    case KSP_CONVERGED_RTOL_NORMAL     : return "CONVERGED_RTOL_NORMAL";
#endif
    case KSP_CONVERGED_RTOL            : return "CONVERGED_RTOL";
    case KSP_CONVERGED_ATOL            : return "CONVERGED_ATOL";
    case KSP_CONVERGED_ITS             : return "CONVERGED_ITS";
    case KSP_CONVERGED_CG_NEG_CURVE    : return "CONVERGED_CG_NEG_CURVE";
    case KSP_CONVERGED_CG_CONSTRAINED  : return "CONVERGED_CG_CONSTRAINED";
    case KSP_CONVERGED_STEP_LENGTH     : return "CONVERGED_STEP_LENGTH";
    case KSP_CONVERGED_HAPPY_BREAKDOWN : return "CONVERGED_HAPPY_BREAKDOWN";
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3, 2, 0 )
    case KSP_CONVERGED_ATOL_NORMAL     : return "CONVERGED_ATOL_NORMAL";
#endif

        /* diverged */
    case KSP_DIVERGED_NULL           : return "DIVERGED_NULL";
    case KSP_DIVERGED_ITS            : return "DIVERGED_ITS";
    case KSP_DIVERGED_DTOL           : return "DIVERGED_DTOL";
    case KSP_DIVERGED_BREAKDOWN      : return "DIVERGED_BREAKDOWN";
    case KSP_DIVERGED_BREAKDOWN_BICG : return "DIVERGED_BREAKDOWN_BICG";
    case KSP_DIVERGED_NONSYMMETRIC   : return "DIVERGED_NONSYMMETRIC";
    case KSP_DIVERGED_INDEFINITE_PC  : return "DIVERGED_INDEFINITE_PC";
#if PETSC_VERSION_LESS_THAN(3,4,0)
    case KSP_DIVERGED_NAN            : return "DIVERGED_NAN";
#endif
    case KSP_DIVERGED_INDEFINITE_MAT : return "DIVERGED_INDEFINITE_MAT";

    case KSP_CONVERGED_ITERATING : return "CONVERGED_ITERATING";

    default: return "INDEFINE_KSP_REASON";

    }
}

} // Nextsim
