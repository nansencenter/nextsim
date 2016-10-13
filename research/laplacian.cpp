/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   finiteelement.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug 24 11:02:45 2015
 */

#include <laplacian.hpp>

#define GMSH_EXECUTABLE gmsh

namespace Nextsim
{
Laplacian::Laplacian()
    :
    vm(Environment::vm()),
    M_mesh(),
    M_matrix(),
    M_vector()
{}

void Laplacian::init()
{
    // std::cout <<"GMSH VERSION= "<< M_mesh.version() <<"\n";
    // M_mesh.setOrdering("bamg");
    // M_mesh.readFromFile("bigarctic10km.msh");

    createGMSHMesh("hypercube.geo");
    M_mesh.setOrdering("gmsh");
    M_mesh.readFromFile("hypercube.msh");

 #if 0
    auto nodes = M_mesh.nodes();
    for (auto it=nodes.begin(), end=nodes.end(); it!=end; ++it)
    {
        auto coords = it->second.coords;

        if (it->first < 10)
            std::cout<< "Nodes "<< it->first
                     << ": coords= ("<< coords[0] << ","<< coords[1] << ","<< coords[2] <<")\n";
    }


    //auto elements = M_mesh.elements();
    auto elements = M_mesh.triangles();
    int nprint = 0;
    for (auto it=elements.begin(), end=elements.end(); it!=end; ++it)
    {
        //if (it->first < 10)
        if (nprint < 10)
        {
            std::cout<< "Elements : "<< it->first <<"\n"
                     << "           number= " << it->second.number <<"\n"
                     << "             type= " << it->second.type <<"\n"
                     << "         physical= " << it->second.physical <<"\n"
                     << "       elementary= " << it->second.elementary <<"\n"
                     << "      numVertices= " << it->second.numVertices <<"\n"
                     << "          indices= (" << it->second.indices[0] << ","
                     << it->second.indices[1] << "," << it->second.indices[2] <<")\n";

        }

        ++nprint;
    }

  #endif

    M_lines = M_mesh.lines();
    int first = 0;
    for (auto it=M_lines.begin(), end=M_lines.end(); it!=end; ++it)
    {
        //if (it->second.physical==158)
        if (it->second.physical==1)
        {
            // ++first;
            // if (first==1)
            //     continue;

            dirichlet_flags.push_back(it->second.indices[0]-1);
        }
    }

    //std::cout<<"DIRICHLET_FLAGS= "<< dirichlet_flags.size() <<"\n";

    std::cout<<"NumNodes     = "<< M_mesh.numNodes() <<"\n";
    std::cout<<"NumElements  = "<< M_mesh.numElements() <<"\n";
    std::cout<<"NumTriangles = "<< M_mesh.numTriangles() <<"\n";
    std::cout<<"NumLines     = "<< M_mesh.numLines() <<"\n";

    // std::string str = "../data/Nps.mpp";
    // M_mesh.project(str);

    // M_mesh.writeTofile("../data/arctic10km.msh");

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_num_nodes = M_mesh.numNodes();

    M_matrix = matrix_ptrtype(new matrix_type(M_num_nodes,M_num_nodes,10));
    M_matrix->zero();

    // M_mass = matrix_ptrtype(new matrix_type(M_num_nodes,M_num_nodes,10));
    // M_mass->zero();

    M_vector = vector_ptrtype(new vector_type(M_num_nodes));
    M_vector->zero();
    //M_vector->setOnes();

    M_solution = vector_ptrtype(new vector_type(M_num_nodes));
    M_solution->zero();

}

void Laplacian::createGMSHMesh(std::string const& geofilename)
{
    std::string gmshgeofile = Environment::nextsimDir().string() + "/mesh/" + geofilename;

    if (fs::exists(gmshgeofile))
    {
        //std::cout<<"NOT FOUND " << fs::absolute( gmshgeofile ).string() <<"\n";
        std::ostringstream gmshstr;
        gmshstr << BOOST_PP_STRINGIZE( GMSH_EXECUTABLE )
                << " -" << 2 << " -part " << 1 << " -clmax " << vm["hsize"].as<double>() << " " << gmshgeofile;

        std::cout << "[Gmsh::generate] execute '" <<  gmshstr.str() << "'\n";
        auto err = ::system( gmshstr.str().c_str() );
    }
    else
    {
        std::cout << "Cannot found " << gmshgeofile <<"\n";
    }
}

double Laplacian::measure(element_type const& element) const
{
    std::vector<double> vertex_0 = M_nodes.find(element.indices[0])->second.coords;
    std::vector<double> vertex_1 = M_nodes.find(element.indices[1])->second.coords;
    std::vector<double> vertex_2 = M_nodes.find(element.indices[2])->second.coords;

    double area = (vertex_1[0]-vertex_0[0])*(vertex_2[1]-vertex_0[1]);
    area -= (vertex_2[0]-vertex_0[0])*(vertex_1[1]-vertex_0[1]);

    return (1./2)*std::abs(area);
}

void Laplacian::assemble()
{
    chrono.restart();
    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        // if (cpt <10)
        //     std::cout<<"***Element "<< it->first <<"\n";

        double area = measure(it->second);
        std::vector<double> x(3);
        std::vector<double> y(3);
        std::vector<double> data(9);
        //std::vector<double> mass_data(9);

        std::vector<int> rcindices(3);
        for (int s=0; s<rcindices.size(); ++s)
            rcindices[s] = it->second.indices[s]-1;

        double m_jk = 0;
        //double mass_jk = 0;

        std::vector<double> fvdata(3);
        double f_j = 0;
        double x_b = 0;
        double y_b = 0;

        for (int i=0; i<3; ++i)
        {
            x[i] = M_nodes.find(it->second.indices[i])->second.coords[0];
            y[i] = M_nodes.find(it->second.indices[i])->second.coords[1];

            x_b += x[i];
            y_b += y[i];
        }
        x_b = x_b/3.0;
        y_b = y_b/3.0;

        int lc = 0;
        for (int j=0; j<3; ++j)
        {
            f_j = 0;
            // x-axis
            int jp1 = (j+1)%3;
            int jp2 = (j+2)%3;

            for (int k=0; k<3; ++k)
            {
                // y-axis
                int kp1 = (k+1)%3;
                int kp2 = (k+2)%3;

                m_jk = (y[jp1]-y[jp2])*(y[kp1]-y[kp2])+(x[jp1]-x[jp2])*(x[kp1]-x[kp2]);
                m_jk = m_jk/(4.0*area);

                //add mass matrix contribution
                //m_jk += ((j == k) ? 2.0 : 1.0)*area/12.0;
                //mass_jk = ((j == k) ? 2.0 : 1.0)*area/12.0;

                data[lc] = m_jk;

                //mass_data[lc] = mass_jk;

                ++lc;

                // if (cpt <1)
                //     std::cout<<"DATA["<< j << ","<< k <<"]= "<< m_jk <<"\n";

                //compute right-hand side contribution
                //f_j += ((j == k) ? 2.0 : 1.0)*2.0*PI*PI*std::sin(PI*x[k])*std::sin(PI*y[k]);
            }

            //f_j = f_j*area/12.0;
            f_j = 2*PI*PI*std::sin(PI*x_b)*std::sin(PI*y_b);
            //f_j = 2*PI*PI*std::sin(PI*x[j])*std::sin(PI*y[j]);
            f_j = f_j*area/3.0;
            fvdata[j] = f_j;
        }

        M_matrix->addMatrix(&rcindices[0], rcindices.size(),
                            &rcindices[0], rcindices.size(), &data[0]);

        // M_mass->addMatrix(&rcindices[0], rcindices.size(),
        //                   &rcindices[0], rcindices.size(), &mass_data[0]);


        M_vector->addVector(&rcindices[0], rcindices.size(), &fvdata[0]);

        ++cpt;
    }

    M_matrix->close();

    std::cout<<"TIMER ASSEMBLY= " << chrono.elapsed() <<"s\n";

    //M_mass->close();

    // apply homogeneous dirichlet boundary conditions
    chrono.restart();
    M_matrix->on(dirichlet_flags,*M_vector);
    std::cout<<"TIMER DBCA= " << chrono.elapsed() <<"s\n";

    std::cout<<"[PETSC MATRIX] CLOSED      = "<< M_matrix->closed() <<"\n";
    std::cout<<"[PETSC MATRIX] SIZE        = "<< M_matrix->size1() << " " << M_matrix->size2() <<"\n";
    std::cout<<"[PETSC MATRIX] SYMMETRIC   = "<< M_matrix->isSymmetric() <<"\n";
    std::cout<<"[PETSC MATRIX] NORM        = "<< M_matrix->linftyNorm() <<"\n";

#if 1
    //M_matrix->printScreen();
    M_matrix->printMatlab("stiffness.m");
    M_vector->printMatlab("rhs.m");

    vector_type v(M_num_nodes);
    M_matrix->diagonal(v);
    v.printMatlab("diag.m");
#endif

}

void Laplacian::solve()
{
    SolverPetsc ksp;

    chrono.restart();
    //ksp.solve(M_matrix, M_solution, M_vector);

    ksp.solve(_matrix=M_matrix,
              _solution=M_solution,
              _rhs=M_vector,
              _ksp="preonly",
              _pc="cholesky",
              _pcfactormatsolverpackage="cholmod",
              _reuse_prec=true
              );

    std::cout<<"TIMER SOLUTION= " << chrono.elapsed() <<"s\n";
    M_solution->printMatlab("solution.m");

    Environment::logMemoryUsage("");
}

void Laplacian::run()
{

    this->init();
    this->assemble();
    this->solve();

    // test bamg usage
    std::vector<int> index(3*M_num_elements);
    std::vector<double> x(M_num_nodes);
    std::vector<double> y(M_num_nodes);

    vector_type mindex(3*M_num_elements);
    vector_type mx(3*M_num_elements);
    vector_type my(3*M_num_elements);
    vector_type mc(3*M_num_elements);

    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        index[3*cpt] = it->second.indices[0];//it->first;
        index[3*cpt+1] = it->second.indices[1];
        index[3*cpt+2] = it->second.indices[2];

  #if 1
        mindex(3*cpt) = it->second.indices[0];//it->first;
        mindex(3*cpt+1) = it->second.indices[1];
        mindex(3*cpt+2) = it->second.indices[2];

        double sum = 0;
        for (int j=0; j<3; ++j)
        {
            sum += M_solution->operator()(it->second.indices[j]-1);
        }
        mc(3*cpt) = sum;
        mc(3*cpt+1) = sum;
        mc(3*cpt+2) = sum;

        for (int i=0; i<3; ++i)
        {
            mx(3*cpt+i) = M_nodes.find(it->second.indices[i])->second.coords[0];
            my(3*cpt+i) = M_nodes.find(it->second.indices[i])->second.coords[1];
        }
#endif

        ++cpt;
    }

    //vector_type M_exact(M_num_nodes);
    M_exact = vector_ptrtype(new vector_type(M_num_nodes));
    std::vector<double> data(M_num_nodes);
    double exact = 0;
    cpt = 0;
    for (auto it=M_nodes.begin(), end=M_nodes.end(); it!=end; ++it)
    {
        x[cpt] = it->second.coords[0];
        y[cpt] = it->second.coords[1];

        exact = std::sin(PI*x[cpt])*std::sin(PI*y[cpt]);
        M_exact->operator()(cpt) = exact;

        data[cpt] = exact;

        //std::cout<<"data["<< cpt <<"]= "<< data[cpt] <<"\n";

        ++cpt;
    }

    mindex.printMatlab("mindex.m");
    mx.printMatlab("mx.m");
    my.printMatlab("my.m");
    mc.printMatlab("mc.m");

    M_exact->printMatlab("exact.m");

    BamgOpts *bamgopt = NULL;
    BamgMesh *bamgmesh = NULL;
    BamgGeom *bamggeom = NULL;

    bamgopt=new BamgOpts();
    bamggeom=new BamgGeom();
    bamgmesh=new BamgMesh();

    //bamgopt->anisomax=1;
    //std::cout<<"Anisomax= "<< bamgopt->anisomax <<"\n";
    //bamgopt->Check();

    int test = BamgConvertMeshx(bamgmesh,bamggeom,&index[0],&x[0],&y[0],M_num_nodes,M_num_elements);

    delete bamggeom;
    delete bamgmesh;
    delete bamgopt;

#if 1
    std::vector<double> x_interp(M_num_nodes);
    std::vector<double> y_interp(M_num_nodes);
    //std::vector<double> data_interp(M_num_nodes);
    double* data_interp;

    x_interp = x;
    y_interp = y;

    Options* options;
    options = new Options();
    //options->AddOption("default");
    // std::cout<<"DEFAULT= "<< options->GetOption("default") <<"\n";
    // double defaultvalue;
    // options->Get(&defaultvalue,"default");

    //std::cout<<"VALUE= "<< defaultvalue <<"\n";

    //std::vector<double> data(M_num_nodes);
    //std::iota(data.begin(),data.end(),0.);

    InterpFromMeshToMesh2dx(&data_interp,&index[0],&x[0],&y[0],M_num_nodes,M_num_elements,
                            &data[0],M_num_nodes,1,&x_interp[0],&y_interp[0],M_num_nodes,
                            options);

    //std::cout<<"First= "<< data_interp[0] <<"\n";

    // for (int i=0; i<M_num_nodes; ++i)
    //     std::cout<<"data_interp["<< i <<"]= "<< data_interp[i] <<"\n";

#endif

    delete options;

    M_exact->scale(-1.);
    M_exact->add(*M_solution);

    //std::cout<<"L2  = "<< M_exact.l2Norm() <<"\n";
    //std::cout<<"LINF= "<< M_exact.linftyNorm() <<"\n";
    //std::cout<<"||u-uh||_L2  = "<< std::sqrt(M_mass->energy(*M_exact)) <<"\n";
    //std::cout<<"||u-uh||_H1  = "<< std::sqrt(M_mass->energy(*M_exact)+M_matrix->energy(*M_exact)) <<"\n";

    this->error();

}

void Laplacian::error()
{
    double l2_error = 0;
    double sh1_error = 0;

    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        double area = measure(it->second);
        std::vector<double> x(3);
        std::vector<double> y(3);

        double l2_contrib = 0;
        double sh1_contrib = 0;
        double entry_contrib = 0;

        for (int i=0; i<3; ++i)
        {
            x[i] = M_nodes.find(it->second.indices[i])->second.coords[0];
            y[i] = M_nodes.find(it->second.indices[i])->second.coords[1];
        }

        int lc = 0;
        for (int j=0; j<3; ++j)
        {
            l2_contrib = 0;
            sh1_contrib = 0;
            // x-axis
            int jp1 = (j+1)%3;
            int jp2 = (j+2)%3;

            for (int k=0; k<3; ++k)
            {
                // y-axis
                int kp1 = (k+1)%3;
                int kp2 = (k+2)%3;

                // semi h1 error
                entry_contrib = (y[jp1]-y[jp2])*(y[kp1]-y[kp2])+(x[jp1]-x[jp2])*(x[kp1]-x[kp2]);
                sh1_contrib += entry_contrib*M_exact->operator()(it->second.indices[k]-1)/(4.0*area);

                // l2 error
                l2_contrib += M_exact->operator()(it->second.indices[k]-1)*((j == k) ? 2.0 : 1.0)*area/12.0;
            }

            l2_error += l2_contrib*M_exact->operator()(it->second.indices[j]-1);
            sh1_error += sh1_contrib*M_exact->operator()(it->second.indices[j]-1);
        }
    }

    std::cout<<"||u-uh||_L2  = "<< std::sqrt(l2_error) <<"\n";
    std::cout<<"||u-uh||_H1  = "<< std::sqrt(l2_error+sh1_error) <<"\n";
}

} // Nextsim
