/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#include <solverpetsc.hpp>
#include <gmshmesh.hpp>
#include <boost/program_options.hpp>
#include <boost/version.hpp>
//#include <src/wrappers/BamgConvertMesh/BamgConvertMesh.h>

#include <src/c/main/globals.h>
//#include <src/c/bamg/Mesh.h>
#include <src/c/modules/Bamgx/Bamgx.h>
#include <src/c/modules/BamgConvertMeshx/BamgConvertMeshx.h>

#include <netcdf>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

using namespace netCDF;
using namespace netCDF::exceptions;

// We are reading 2D data, a 6 x 12 grid.
static const int NX = 6;
static const int NY = 12;

// Return this in event of a problem.
static const int NC_ERR = 2;

namespace Nextsim
{
    po::options_description descrOptions();
}

int main(int argc, char** argv )
{
    using namespace Nextsim;

    Environment env(argc, argv);

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

#if 0
    Environment::logMemoryUsage("toto");

    int row=2;
    int nrow=1;
    testmat.addMatrix(&row, nrow, cols, ncols, data);
#endif

    testmat.close();

#if 0
    for (int i = 0; i < testmat.size1(); ++i)
        std::cout<<"MAT["<< i <<"," << i <<"]= "<< testmat(i,i) <<"\n";

    testmat.printScreen();
    testmat.zero();
    testmat.printMatlab("test.m");

    std::cout<<"CLOSED= "<< testmat.closed() <<"\n";

    std::cout<<"SYMM= "<< testmat.isSymmetric() <<"\n";
    std::cout<<"NORM= "<< testmat.linftyNorm() <<"\n";
#endif

#endif

#if 0
    testmat.clear();
    Environment::logMemoryUsage("PASS");
#endif


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
    PetscLogEvent MATRIX_GENERATE,MATRIX_READ;

#if 1
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
    //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"writing matrix in binary to matrix.dat ...\n");CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"../data/matrix.dat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = MatView(A,viewer);CHKERRQ(ierr);
    ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(MATRIX_GENERATE,0,0,0,0);CHKERRQ(ierr);

#else

    /* Read the matrix again as a sequential matrix */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"reading matrix in binary from matrix.dat ...\n");CHKERRQ(ierr);
    Mat B;
    ierr = PetscViewerBinaryOpen(Environment::comm(),"../data/matrix.dat",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = MatCreate(Environment::comm(),&B);CHKERRQ(ierr);
    ierr = MatSetType(B,MATSEQAIJ);CHKERRQ(ierr);
    ierr = MatLoad(B,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

#endif

    std::cout<<"reading matrix in binary from matrix.dat ...\n";
    std::string file = "../data/matrix.dat";
    MatrixPetsc B;
    B.loadBinary(file);
    B.printScreen();

    std::cout<<"B.size1= "<< B.size1() <<"\n";
    std::cout<<"B.size2= "<< B.size2() <<"\n";

    file = "../data/vector.dat";
    std::cout<<"writing vector in binary to vector.dat ...\n";
    testvec.printBinary(file);
    testvec.clear();

    std::cout<<"reading vector in binary from vector.dat ...\n";
    VectorPetsc C;
    C.loadBinary(file);
    C.printScreen();


    //namespace po = boost::program_options;
    po::options_description desc = descrOptions();
    po::variables_map vm;

    try
    {
        //po::store(po::parse_command_line(argc, argv, desc),vm);

        po::store(po::command_line_parser(argc, argv)
                  .options(desc)
                  .style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise)
                  //.style(po::command_line_style::allow_long | po::command_line_style::long_allow_adjacent | po::command_line_style::long_allow_next)
                  .allow_unregistered()
                  .run(),
                  vm);

        if ( vm.count("help")  )
        {
            std::cout<< "BOOST VERSION= "<< BOOST_LIB_VERSION <<"\n";
            std::cout << "Basic Command Line Parameter Application" <<"\n"
                      << desc << "\n";
            //return SUCCESS;
            return 0;
        }

        if ( vm.count( "config-file" ) || vm.count( "config-files" ) )
        {
#if 0
            if ( fs::exists( vm["config-file"].as<std::string>() ) )
            {
                std::ifstream ifs( vm["config-file"].as<std::string>().c_str() );
                po::store( parse_config_file( ifs, desc, true ), vm );
                po::notify( vm );
            }
            else
            {
                std::cout << "Cannot found " << "config-file `" << vm["config-file"].as<std::string>() <<"`\n";
                return 1;
            }
#endif

            if ( vm.count( "config-file" ) )
            {
                if ( fs::exists( vm["config-file"].as<std::string>() ) )
                {
                    std::cout << "Reading " << vm["config-file"].as<std::string>() << "...\n";
                    std::ifstream ifs( vm["config-file"].as<std::string>().c_str() );
                    po::store( parse_config_file( ifs, desc, true ), vm );
                    //po::notify( vm );
                }
                else
                {
                    std::cout << "Cannot found " << "config-file `" << vm["config-file"].as<std::string>() <<"`\n";
                    return 1;
                }
            }

            if ( vm.count( "config-files" ) )
            {
                std::vector<std::string> configFiles = vm["config-files"].as<std::vector<std::string> >();

                // reverse order (priorty for the last)
                //std::reverse(configFiles.begin(),configFiles.end());

                for ( std::string cfgfile : configFiles )
                {
                    if ( fs::exists( cfgfile ) )
                    {
                        std::cout << "Reading " << cfgfile << "...\n";
                        std::ifstream ifs( cfgfile.c_str() );
                        po::store( parse_config_file( ifs, desc, true ), vm );
                    }
                    else
                    {
                        std::cout << "Cannot found " << "config-file `" << cfgfile <<"`\n";
                        return 1;
                    }
                }
            }

        }

        po::notify(vm);

    }

#if 0
    catch(po::error& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl <<"\n";
        std::cerr << desc <<"\n";
        //return ERROR_IN_COMMAND_LINE;
        return -1;
    }

    catch(po::duplicate_option_error const& e)
    {
        return -1;
    }
#endif

    // catches program_options exceptions

    catch (po::multiple_occurrences const& e)
    {
        std::cout << "Command line or config file option parsing error: " << e.what() << "\n"
                  << "  o faulty option: " << e.get_option_name() << "\n"
                  << "Error: the .cfg file or some options may not have been read properly\n";

        return -2;
    }

    catch (po::ambiguous_option const& e)
    {
        std::cout << "Command line or config file option parsing error: " << e.what() << "\n"
                  << "  o faulty option: " << e.get_option_name() << "\n"
                  << "  o possible alternatives: " ;
        std::for_each( e.alternatives().begin(), e.alternatives().end(), []( std::string const& s )
                       {
                           std::cout << s << " ";
                       } );
        std::cout << "\n"
                  << "Error: the .cfg file or some options may not have been read properly\n";

        return -3;
    }

    catch ( std::exception& e )
    {
        std::cout << "Application option parsing: unknown option:" << e.what() << " (the .cfg file or some options may not have been read properly)\n";
    }


    catch ( ... )
    {
        std::cout << "Application option parsing: unknown exception triggered  (the .cfg file or some options may not have been read properly)\n";
    }

    // std::cout<<"nx= "<< vm["nx"]. as<int>() <<"\n";
    // std::cout<<"ny= "<< vm["ny"]. as<int>() <<"\n";
    // std::cout<<"dx= "<< vm["dx"]. as<double>() <<"\n";
    // std::cout<<"dy= "<< vm["dy"]. as<double>() <<"\n";

#if 0
    BamgOpts *bamgopt = NULL;
    BamgMesh *bamgmesh = NULL;
    BamgGeom *bamggeom = NULL;

    bamgopt=new BamgOpts();
    bamggeom=new BamgGeom();
    bamgmesh=new BamgMesh();

    bamgopt->anisomax=1;

    std::cout<<"Anisomax= "<< bamgopt->anisomax <<"\n";

    //bamgopt->Check();

    delete bamggeom;
    delete bamgmesh;
    delete bamgopt;
#endif

    try
    {
        // This is the array we will read.
        int dataIn[NX][NY];

        // Open the file for read access
        NcFile dataFile("../data/simple_xy.nc", NcFile::read);

        // Retrieve the variable named "data"
        NcVar data=dataFile.getVar("data");
        if(data.isNull()) return NC_ERR;
        data.getVar(dataIn);

        // Check the values.
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
            {
                //std::cout<<"simple_xy.nc= "<< dataIn[i][j] <<"\n";
                if (dataIn[i][j] != i * NY + j)
                    return NC_ERR;
            }

        // The netCDF file is automatically closed by the NcFile destructor
        //cout << "*** SUCCESS reading example file simple_xy.nc!" << endl;

        //return 0;
    }
    catch(NcException& e)
    {
        e.what();
        cout<<"FAILURE*************************************"<<endl;
        return NC_ERR;
    }

    GmshMesh mesh;

    std::cout <<"VERSION= "<< mesh.version() <<"\n";


    mesh.readFromFile("bigarctic10km.msh");

    auto nodes = mesh.nodes();

    for (auto it=nodes.begin(), end=nodes.end(); it!=end; ++it)
    {
        auto coords = it->second.coords;

        if (it->first < 10)
            std::cout<< "Nodes "<< it->first << ": coords= ("<< coords[0] << ","<< coords[1] << ","<< coords[2] <<")\n";
    }



    auto elements = mesh.elements();

    for (auto it=elements.begin(), end=elements.end(); it!=end; ++it)
    {
        //auto coords = it->second.coords;

        if (it->first < 10)
        {
            std::cout<< "Elements : "<< it->first <<"\n"
                     << "           number= " << it->second.number <<"\n"
                     << "             type= " << it->second.type <<"\n"
                     << "         physical= " << it->second.physical <<"\n"
                     << "       elementary= " << it->second.elementary <<"\n"
                     << "      numVertices= " << it->second.numVertices <<"\n";
                //<< "          indices= (" << indices[] <<"\n"

        }
    }

    std::cout<<"NumNodes   = "<< mesh.numNodes() <<"\n";
    std::cout<<"NumElements= "<< mesh.numElements() <<"\n";

    int nods = mesh.numNodes();
    int nels = mesh.numElements();

    int *index;
    double *x;
    double *y;

    index = new int[nels];
    x = new double[nods];
    y = new double[nods];

    // int cpt = 0;

    // for (auto it=elements.begin(), end=elements.end(); it!=end; ++it)
    // {
    //     index[cpt] = it->first;
    //     ++cpt;
    // }

    int cpt = 0;

    for (auto it=nodes.begin(), end=nodes.end(); it!=end; ++it)
    {
        auto coords = it->second.coords;

        index[cpt] = it->first-1;

        x[cpt] = coords[0];
        y[cpt] = coords[1];


        //if (cpt < 10)
        //    std::cout<<"I= "<< index[cpt] << " X= "<< x[cpt] << " Y= "<< y[cpt] <<"\n";

        //if (cpt < 10)
        //    std::cout<<"I= "<< index[cpt] <<"\n";

        ++cpt;
    }

    //int Bmesh = BamgConvertMeshx;


    BamgOpts *bamgopt = NULL;
    BamgMesh *bamgmesh = NULL;
    BamgGeom *bamggeom = NULL;

    bamgopt=new BamgOpts();
    bamggeom=new BamgGeom();
    bamgmesh=new BamgMesh();

    //int test = BamgConvertMeshx(bamgmesh,bamggeom,index,x,y,nods,nels);

    //bamg::Mesh Th(index,x,y,nods,nels);

    bamg::Mesh Th(x,y,nods);

    std::cout<< "BAMG OBJECTS*******************\n";
    std::cout<< "MAX_NVERTEX= " << Th.maxnbv <<"\n";
    std::cout<< "MAX_NELEMENT= " << Th.maxnbt <<"\n";
    //maxnbt

    Th.Gh.WriteGeometry(bamggeom,bamgopt);
    Th.WriteMesh(bamgmesh,bamgopt);

    //std::cout<<"TRIAN= "<< bamgmesh->TrianglesSize[1] <<"\n";

    delete bamggeom;
    delete bamgmesh;
    delete bamgopt;

}
