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

    testmat.printScreen();
    //testmat.zero();
    testmat.printMatlab("test.m");

    //std::cout<<"CLOSED= "<< testmat.closed() <<"\n";

    std::cout<<"SYMM= "<< testmat.isSymmetric() <<"\n";
    std::cout<<"NORM= "<< testmat.linftyNorm() <<"\n";
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

#if 1
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



}
