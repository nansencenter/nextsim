/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

//#include <environment.hpp>
//#include <matrixpetsc.hpp>
//#include <vectorpetsc.hpp>
#include <solverpetsc.hpp>

int main(int argc, char** argv )
{
    using namespace Nextsim;

    Environment env(argc, argv);

    //std::cout<<"STATUS= "<< boost::mpi::environment::finalized() <<"\n";
    //Communicator comm;
    //Communicator comm = Environment::comm();

    //MatrixPetsc testmat(10,5,15,Environment::comm());

#if 1
    MatrixPetsc testmat(10,10,5);

    std::cout<<"SIZE1= "<< testmat.size1() <<"\n";
    std::cout<<"SIZE2= "<< testmat.size2() <<"\n";

    for (int i = 0; i < testmat.size1(); ++i)
        testmat.setValue(i,i,1.);

    VectorPetsc testvec(10);
    testvec.setOnes();


#if 0
    int rows[3];
    rows[0] = 0; rows[1] = 4; rows[2] = 9;
    int nrows = 3;

    int cols[4];
    cols[0] = 0; cols[1] = 2; cols[2] = 4, cols[3] = 8;
    int ncols = 4;
    double data[12];

    for (int i=0; i < 12; ++i)
        data[i] = i+1;

    testmat.addMatrix(rows, nrows, cols, ncols, data);
#endif

    //Environment::logMemoryUsage("toto");

    // int row=2;
    // int nrow=1;
    // testmat.addMatrix(&row, nrow, cols, ncols, data);

    testmat.close();

    // for (int i = 0; i < testmat.size1(); ++i)
    //     std::cout<<"MAT["<< i <<"," << i <<"]= "<< testmat(i,i) <<"\n";

    //testmat.printScreen();
    //testmat.zero();
    //testmat.printMatlab("test.m");

    //std::cout<<"CLOSED= "<< testmat.closed() <<"\n";

    //std::cout<<"SYMM= "<< testmat.isSymmetric() <<"\n";
    //std::cout<<"NORM= "<< testmat.linftyNorm() <<"\n";
#endif

    //testmat.clear();
    //Environment::logMemoryUsage("PASS");

#if 0
    VectorPetsc testvec(10);
    testvec.setOnes();

    //testvec.printMatlab("vec.m");

    //testvec.clear();

    std::cout<<"SIZE= "<< testvec.size() <<"\n";
    std::cout<<"NORM= "<< testvec.linftyNorm() <<"\n";

    //testvec.zero();

    //testvec.set(1,3.);

    testvec.printMatlab("vec.m");

    std::cout<<"MAX= "<< testvec.min() <<"\n";

    //testvec.add(2,5.);

    //testvec.scale(2.);

    // for (int i = 0; i < testvec.size(); ++i)
    //     std::cout<<"VEC["<< i <<"]= "<< testvec(i) <<"\n";

    testvec.printMatlab("vec.m");

    VectorPetsc testvec2(10);
    //testvec2.setConstant(2.);

    //testvec.add(2.,testvec2);

    //testvec.close();

    testvec2.addVector(testvec, testmat);

    //testvec -= testvec2;

    testvec2.printMatlab("vec2.m");
#endif

#if 0
    SolverPetsc ksp;

    VectorPetsc rhs(10);
    rhs.setOnes();
    rhs.printMatlab("rhs.m");

    VectorPetsc solution(10);

    boost::timer chrono;
    //std::map<std::string, std::pair<boost::mpi::timer, double> > timers;
    chrono.restart();
    ksp.solve(testmat, solution, rhs,1e-12,100);
    std::cout<<"TMER= " << chrono.elapsed() <<"s\n";

    solution.printMatlab("solution.m");

#endif
    //testsolver.init();

    // std::cout<<"solverType= "<< testsolver.solverType() <<"\n";
    // std::cout<<"PreconditionerType= "<< testsolver.preconditionerType() <<"\n";
    // std::cout<<"MatSolverPackageType= "<< testsolver.matSolverPackageType() <<"\n";


    Mat A;
    PetscInt M=79, N=16;
    PetscInt       i,j,Ii,J,Istart,Iend, n=4, m=4;
    PetscViewer viewer;
    PetscScalar v;
    int ierr = 0;
    //std::string file = "spd-real-int64-float64";
    PetscLogEvent MATRIX_GENERATE,MATRIX_READ;

#if 0
    // /* Read the matrix again as a sequential matrix */
    // ierr = PetscViewerBinaryOpen(Environment::comm(),"spd-real-int64-float64",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    // ierr = MatCreate(Environment::comm(),&A);CHKERRQ(ierr);
    // ierr = MatSetType(A,MATSEQAIJ);CHKERRQ(ierr);
    // ierr = MatLoad(A,viewer);CHKERRQ(ierr);
    // ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

    /* PART 1:  Generate matrix, then write it in binary format */

    ierr = PetscLogEventRegister("Generate Matrix",0,&MATRIX_GENERATE);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(MATRIX_GENERATE,0,0,0,0);CHKERRQ(ierr);

    /* Generate matrix */
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
    for (Ii=Istart; Ii<Iend; Ii++) {
        v = -1.0; i = Ii/n; j = Ii - i*n;
        if (i>0)   {J = Ii - n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
        if (i<m-1) {J = Ii + n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
        if (j>0)   {J = Ii - 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
        if (j<n-1) {J = Ii + 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
        v = 4.0; ierr = MatSetValues(A,1,&Ii,1,&Ii,&v,ADD_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"writing matrix in binary to matrix.dat ...\n");CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"matrix.dat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = MatView(A,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(MATRIX_GENERATE,0,0,0,0);CHKERRQ(ierr);




    /* Read the matrix again as a sequential matrix */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"reading matrix in binary from matrix.dat ...\n");CHKERRQ(ierr);
    Mat B;
    ierr = PetscViewerBinaryOpen(Environment::comm(),"matrix.dat",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = MatCreate(Environment::comm(),&B);CHKERRQ(ierr);
    ierr = MatSetType(B,MATSEQAIJ);CHKERRQ(ierr);
    ierr = MatLoad(B,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

#endif

    //testmat.printScreen();

    std::string file = "matrix.dat";
    std::cout<<"writing matrix in binary to matrix.dat ...\n";
    // ierr = PetscViewerBinaryOpen(Environment::comm(),file.c_str(),FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    // ierr = MatView(testmat.mat(),viewer);CHKERRQ(ierr);
    // ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    testmat.printBinary(file);
    testmat.clear();

    ierr = PetscPrintf(PETSC_COMM_WORLD,"reading matrix in binary from matrix.dat ...\n");CHKERRQ(ierr);
#if 0
    Mat B;
    ierr = PetscViewerBinaryOpen(Environment::comm(),file.c_str(),FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = MatCreate(Environment::comm(),&B);CHKERRQ(ierr);
    ierr = MatSetType(B,MATSEQAIJ);CHKERRQ(ierr);
    ierr = MatLoad(B,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = MatDestroy(&B);CHKERRQ(ierr);
#endif

    MatrixPetsc B;
    B.loadBinary(file);
    B.printScreen();

    std::cout<<"B.size1= "<< B.size1() <<"\n";
    std::cout<<"B.size2= "<< B.size2() <<"\n";

    file = "vector.dat";
    std::cout<<"writing vector in binary to vector.dat ...\n";
    testvec.printBinary(file);
    testvec.clear();
    ierr = PetscPrintf(PETSC_COMM_WORLD,"reading vector in binary from vector.dat ...\n");CHKERRQ(ierr);

#if 0
    Vec C;
    ierr = PetscViewerBinaryOpen(Environment::comm(),file.c_str(),FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecCreate(Environment::comm(),&C);CHKERRQ(ierr);
    ierr = VecSetType(C,VECSEQ);CHKERRQ(ierr);
    ierr = VecLoad(C,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = VecView(C,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecDestroy(&C);CHKERRQ(ierr);
#endif

    VectorPetsc C;
    C.loadBinary(file);
    C.printScreen();



}
