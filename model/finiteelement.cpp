/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   finiteelement.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug 24 11:02:45 2015
 */

#include <finiteelement.hpp>
#include <constants.hpp>

#include <boost/format.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#define GMSH_EXECUTABLE gmsh

namespace Nextsim
{
FiniteElement::FiniteElement()
    :
    vm(Environment::vm()),
    M_mesh(),
    M_matrix(),
    M_vector()
{}

void
FiniteElement::init()
{
    // std::cout <<"GMSH VERSION= "<< M_mesh.version() <<"\n";
    M_mesh.setOrdering("bamg");
    M_mesh.readFromFile("bigarctic10km.msh");

    // createGMSHMesh("hypercube.geo");
    // M_mesh.setOrdering("gmsh");
    // M_mesh.readFromFile("hypercube.msh");

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
    //int first = 0;
    for (auto it=M_lines.begin(), end=M_lines.end(); it!=end; ++it)
    {
        //if (it->second.physical==158)
        //if (it->second.physical==161)
        //if (it->second.physical==1)
        if (it->second.physical==161)
        {
            //++first;
            // if (first==1)
            //     continue;

            dirichlet_flags.push_back(it->second.indices[0]-1);
            dirichlet_flags.push_back(it->second.indices[0]-1+M_mesh.numNodes());
        }
    }

    // for (int flag=0; flag < 10; ++flag)
    // {
    //     std::cout<<"flag["<< flag <<"]= "<< dirichlet_flags[flag] <<"\n";
    // }

    std::cout<<"DIRICHLET_FLAGS= "<< dirichlet_flags.size() <<"\n";

    std::cout<<"NumNodes     = "<< M_mesh.numNodes() <<"\n";
    std::cout<<"NumElements  = "<< M_mesh.numElements() <<"\n";
    std::cout<<"NumTriangles = "<< M_mesh.numTriangles() <<"\n";
    std::cout<<"NumLines     = "<< M_mesh.numLines() <<"\n";

    // std::string str = "../data/Nps.mpp";
    // M_mesh.project(str);

    M_mesh.stereographicProjection();
    M_mesh.writeTofile("arctic10km.msh");

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_num_nodes = M_mesh.numNodes();

    M_matrix = matrix_ptrtype(new matrix_type(2*M_num_nodes,2*M_num_nodes,20));
    M_matrix->zero();

    // M_mass = matrix_ptrtype(new matrix_type(M_num_nodes,M_num_nodes,10));
    // M_mass->zero();

    M_vector = vector_ptrtype(new vector_type(2*M_num_nodes));
    M_vector->zero();
    //M_vector->setOnes();

    M_solution = vector_ptrtype(new vector_type(2*M_num_nodes));
    M_solution->zero();

}

void
FiniteElement::createGMSHMesh(std::string const& geofilename)
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

double
FiniteElement::jacobian(element_type const& element) const
{
    std::vector<double> vertex_0 = M_nodes.find(element.indices[0])->second.coords;
    std::vector<double> vertex_1 = M_nodes.find(element.indices[1])->second.coords;
    std::vector<double> vertex_2 = M_nodes.find(element.indices[2])->second.coords;

    double jac = (vertex_1[0]-vertex_0[0])*(vertex_2[1]-vertex_0[1]);
    jac -= (vertex_2[0]-vertex_0[0])*(vertex_1[1]-vertex_0[1]);

    return jac;
}

double
FiniteElement::measure(element_type const& element) const
{
    return (1./2)*std::abs(jacobian(element));
}

std::vector<double>
FiniteElement::shapeCoeff(element_type const& element) const
{
    std::vector<double> x(3);
    std::vector<double> y(3);

    for (int i=0; i<3; ++i)
    {
        x[i] = M_nodes.find(element.indices[i])->second.coords[0];
        y[i] = M_nodes.find(element.indices[i])->second.coords[1];
    }

    std::vector<double> coeff(6);
    double jac = jacobian(element);

    for (int k=0; k<6; ++k)
    {
        int kp1 = (k+1)%3;
        int kp2 = (k+2)%3;

        if (k<3)
        {
            coeff[k] = (y[kp1]-y[kp2])/jac;
        }
        else
        {
            coeff[k] = (x[kp2]-x[kp1])/jac;
        }
    }

    return coeff;
}


void
FiniteElement::assemble()
{
    chrono.restart();
    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        // if (cpt <10)
        //     std::cout<<"***Element "<< it->first <<"\n";

        if (cpt==0)
        {
            for (int tt=0; tt<6; ++tt)
                std::cout<<"shape["<< tt <<"]= "<< this->shapeCoeff(it->second)[tt] <<"\n";
        }

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


        //M_vector->addVector(&rcindices[0], rcindices.size(), &fvdata[0]);

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

#if 0
    //M_matrix->printScreen();
    M_matrix->printMatlab("stiffness.m");
    M_vector->printMatlab("rhs.m");

    vector_type v(M_num_nodes);
    M_matrix->diagonal(v);
    v.printMatlab("diag.m");
#endif

}

void
FiniteElement::solve()
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

void
FiniteElement::run()
{

    this->init();

#if 1
    //this->assemble();
    this->initialdata();
    this->solve();

    // test bamg usage
    std::vector<int> index(3*M_num_elements);
    std::vector<double> x(M_num_nodes);
    std::vector<double> y(M_num_nodes);

    //vector_type mindex(3*M_num_elements);
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
        // mindex(3*cpt) = it->second.indices[0];//it->first;
        // mindex(3*cpt+1) = it->second.indices[1];
        // mindex(3*cpt+2) = it->second.indices[2];

        double sum = 0;
        for (int j=0; j<3; ++j)
        {
            sum += M_solution->operator()(it->second.indices[j]-1+M_num_nodes);
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

    mx.printMatlab("mx.m");
    my.printMatlab("my.m");
    mc.printMatlab("mc.m");

#if 0
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

    //mindex.printMatlab("mindex.m");
    mx.printMatlab("mx.m");
    my.printMatlab("my.m");
    mc.printMatlab("mc.m");

    // M_exact->printMatlab("exact.m");

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

    std::vector<double> x_interp(M_num_nodes);
    std::vector<double> y_interp(M_num_nodes);
    //std::vector<double> data_interp(M_num_nodes);
    double* data_interp;

    x_interp = x;
    y_interp = y;

    Options* options;
    options = new Options();

    InterpFromMeshToMesh2dx(&data_interp,&index[0],&x[0],&y[0],M_num_nodes,M_num_elements,
                            &data[0],M_num_nodes,1,&x_interp[0],&y_interp[0],M_num_nodes,
                            options);

    //std::cout<<"First= "<< data_interp[0] <<"\n";

    // for (int i=0; i<M_num_nodes; ++i)
    //     std::cout<<"data_interp["<< i <<"]= "<< data_interp[i] <<"\n";

    delete options;

    // M_exact->scale(-1.);
    // M_exact->add(*M_solution);
    // this->error();


    double date_time = 10.;
    //boost::gregorian::date dt = boost::date_time::parse_date<boost::gregorian::date>( "1899-12-30", boost::date_time::ymd_order_iso );
    //dt = boost::gregorian::date_duration( static_cast<long>( floor(date_time) ) );
    //std::cout<< dt <<"\n";
    //from_date_string( )
    //std::string str = boost::gregorian::to_date_string( date_time );

    //boost::gregorian::date dt = boost::date_time::parse_date<boost::gregorian::date>( "1899-12-30", boost::date_time::ymd_order_iso );
    //boost::gregorian::date dt = boost::date_time::parse_date<boost::gregorian::date>( "1400-01-01", boost::date_time::ymd_order_iso );
    //dt += boost::gregorian::date_duration( static_cast<long>( floor(date_time) ) );

 #if 0
    std::string value = "1900-01-30";
    boost::gregorian::date epoch = boost::date_time::parse_date<boost::gregorian::date>( "1899-12-30", boost::date_time::ymd_order_iso);
    boost::gregorian::date dt = boost::date_time::parse_date<boost::gregorian::date>( value, boost::date_time::ymd_order_iso);

    std::cout<<"time 1 = "<< epoch <<"\n";
    std::cout<<"time 2 = "<< dt <<"\n";

    boost::gregorian::date_duration diff = dt - epoch;
    std::cout<<"diff = "<< diff.days()*24 <<"\n";
#endif

    using namespace boost::gregorian;
    boost::gregorian::date epoch = from_string("1899-12-30");
    boost::gregorian::date dt = from_string("1900-01-30");

    std::cout<<"epoch = "<< epoch <<"\n";
    std::cout<<"dt    = "<< dt <<"\n";
#endif

    this->initialdata();

#endif


}

void
FiniteElement::error()
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

void
FiniteElement::performSimulation()
{
    bool simul_running=true;
    int n_regrid = 0;
    int regrid_counter = 0;
    int cpt=0;
}

void
FiniteElement::initialdata()
{
    int sys_size = 2*M_num_nodes;
    //std::vector<double> VT(2*M_num_nodes,0);
    double VT = 1;
    double Vcor = 1;
    double Vair = 1;
    double Voce = 1;

    std::vector<double> VTold(2*M_num_nodes,0);
    std::vector<double> VToldold(2*M_num_nodes,0);
    std::vector<double> UT(2*M_num_nodes,0);
    std::vector<double> UM(2*M_num_nodes,0);
    //std::vector<double> Vair(2*M_num_nodes,0);
    //std::vector<double> Voce(2*M_num_nodes,0);
    double Vair_factor  = 1.;
    double Voce_factor  = 1.;
    double basal_factor = 1.;
    double g_ssh = 0.;
    double fcor = 0.;
    //surface=mesh.element.surf;
    std::vector<double> sigma(3,0);
    double divergence_rate = 0;
    double h_ridged_thin_ice = 0;
    double h_ridged_thick_ice = 0;
    double damage = 0;

    bool VTold_exists = false;
    bool VToldold_exists = false;

    if (vm["simul_in.newice_type"].as<std::string>() == "Thin-ice-category" )
    {
        std::vector<double> h_thin(M_num_elements,0);
        std::vector<double> hs_thin(M_num_elements,0);
    }

    // unit elasticity tensor
    std::vector<double> Dunit(9,0);
    std::vector<double> Dunit_comp(9,0);
    std::vector<double> Mass(9);
    std::vector<double> B0T(18,0);

    double nu0 = vm["simul_in.nu0"].as<double>();
    double young = vm["simul_in.young"].as<double>();
    double rhoi = physical::rhoi;
    double rhos = physical::rhos;
    double timestep = vm["simul_in.timestep"].as<double>();
    double ridging_exponent = vm["simul_in.ridging_exponent"].as<double>();
    double ocean_turning_angle_rad = (PI/180.)*vm["simul_in.oceanic_turning_angle"].as<double>();

    double days_in_sec = 24.*3600.;
    double divergence_min = (1./days_in_sec)*vm["simul_in.divergence_min"].as<double>();
    double compression_factor = vm["simul_in.compression_factor"].as<double>();
    double exponent_compression_factor = vm["simul_in.exponent_compression_factor"].as<double>();

    double concentr = vm["simul_in.init_concentration"].as<double>();
    //std::vector<double> concentration(M_num_elements,concentr);
    double h = 0.;
    double h_thin = 0.;

    double hs = 0.;
    double hs_thin = 0.;

    double coef_V, coef_Voce, coef_Vair, coef_basal, coef_X, coef_Y, coef_C;
    double coef = 0;
    double coef_P = 0.;
    double mass_e = 0;
    double surface_e = 0;
    double g_ssh_e_x = 0.;
    double g_ssh_e = 0.;
    double g_ssh_e_y = 0.;

    std::vector<double> sigma_P(3,0); /* temporary variable for the resistance to the compression */
    std::vector<double> B0Tj_Dunit(6,0);
    std::vector<double> B0Tj_Dunit_B0Ti(4,0);
    std::vector<double> B0T_Dunit_B0T(36,0);
    std::vector<double> B0Tj_Dunit_comp(6,0);
    std::vector<double> B0Tj_Dunit_comp_B0Ti(4,0);
    std::vector<double> B0T_Dunit_comp_B0T(36,0);
    std::vector<double> B0Tj_sigma_h(2,0);
    double B0Tj_Dunit_tmp0, B0Tj_Dunit_tmp1;
    double B0Tj_Dunit_B0Ti_tmp0, B0Tj_Dunit_B0Ti_tmp1, B0Tj_Dunit_B0Ti_tmp2, B0Tj_Dunit_B0Ti_tmp3;
    double B0Tj_Dunit_comp_tmp0, B0Tj_Dunit_comp_tmp1;
    double B0Tj_Dunit_comp_B0Ti_tmp0, B0Tj_Dunit_comp_B0Ti_tmp1, B0Tj_Dunit_comp_B0Ti_tmp2, B0Tj_Dunit_comp_B0Ti_tmp3;
    double mloc = 0;


    for (int k=0; k<6; k+=3)
    {
        for (int kk=0; kk<2; ++kk )
        {
            Dunit[k+kk] = (1-((k+kk)%2)*(1-nu0))/(1-std::pow(nu0,2.));
            Dunit_comp[k+kk] = 1.;
        }
    }
    Dunit[8] = (1-nu0)/(2.*(1-std::pow(nu0,2.)));

    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            Mass[3*i+j] = ((i == j) ? 2.0 : 1.0)/12.0;
            //std::cout<< std::left << std::setw(12) << Mass[3*i+j] <<"  ";
        }

        //std::cout<<"\n";
    }

    chrono.restart();
    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        // if (cpt <10)
        //     std::cout<<"***Element "<< it->first <<"\n";

        // if (cpt==0)
        // {
        //     for (int tt=0; tt<6; ++tt)
        //         std::cout<<"shape["<< tt <<"]= "<< this->shapeCoeff(it->second)[tt] <<"\n";
        // }

        for (int i=0; i<18; ++i)
        {
            if (i < 3)
            {
                B0T[2*i] = this->shapeCoeff(it->second)[i];
                B0T[12+2*i] = this->shapeCoeff(it->second)[i+3];
                B0T[13+2*i] = this->shapeCoeff(it->second)[i];
            }
            else if (i < 6)
            {
                B0T[2*i+1] = this->shapeCoeff(it->second)[i];
            }

        }

        //std::cout<<"\n";
        // if (cpt == 0)
        //     for (int i=0; i<3; ++i)
        //     {
        //         for (int j=0; j<6; ++j)
        //         {
        //             std::cout<< std::left << std::setw(12) << B0T[6*i+j] <<"  ";
        //         }
        //         std::cout<<"\n";
        //     }


        coef = young*(1-damage)*h*std::exp(ridging_exponent*(1-concentr));
        coef_P = 0.;
        if(divergence_rate < 0.)
        {
            coef_P = compression_factor*std::pow(h,exponent_compression_factor)*std::exp(ridging_exponent*(1-concentr));
            coef_P = coef_P/(std::abs(divergence_rate)+divergence_min);
        }


        /* Compute the value that only depends on the element */
        mass_e = rhoi*h + rhos*hs;
        mass_e = (concentr > 0.) ? (mass_e/concentr):0.;
        surface_e = this->measure(it->second);

        // /* compute the x and y derivative of g*ssh */
        g_ssh_e_x = 0.;
        g_ssh_e_y = 0.;
        for(int i=0; i<3; i++)
        {
            g_ssh_e = g_ssh;   /* g*ssh at the node k of the element e */
            g_ssh_e_x += this->shapeCoeff(it->second)[i]*g_ssh_e; /* x derivative of g*ssh */
            g_ssh_e_y += this->shapeCoeff(it->second)[i+3]*g_ssh_e; /* y derivative of g*ssh */
        }

        coef_C    = mass_e*fcor;              /* for the Coriolis term */
        coef_V    = mass_e/timestep;             /* for the inertial term */
        coef_X    = - mass_e*g_ssh_e_x;              /* for the ocean slope */
        coef_Y    = - mass_e*g_ssh_e_y;              /* for the ocean slope */
        coef_Vair = Vair_factor;             /* for the wind stress */
        coef_Voce = Voce_factor;             /* for the ocean stress */
        coef_basal = basal_factor;            /* for the basal stress */


        for(int j=0; j<3; j++)
        {

            /* The rigidity matrix that will be multiplied by E and by the surface
             * is given by the product B0'*matrix.Dunit*B0
             * This product is computed for the indices:
             * 2*i  ,2*j   -> B0Tj_Dunit_B0Ti[0]
             * 2*i  ,2*j+1 -> B0Tj_Dunit_B0Ti[1]
             * 2*i+1,2*j   -> B0Tj_Dunit_B0Ti[2]
             * 2*i+1,2*j+1 -> B0Tj_Dunit_B0Ti[3] */

            /* new version without assembling zero */
            for(int i=0; i<3; i++)
            {
                /* product of the first line of B0T' and the matrix Dunit */
                B0Tj_Dunit_tmp0 = 0.;
                B0Tj_Dunit_tmp1 = 0.;

                B0Tj_Dunit_comp_tmp0 = 0.;
                B0Tj_Dunit_comp_tmp1 = 0.;

                for(int kk=0; kk<3; kk++)
                {
                    B0Tj_Dunit_tmp0 += B0T[kk*6+2*j]*Dunit[3*i+kk];
                    B0Tj_Dunit_tmp1 += B0T[kk*6+2*j+1]*Dunit[3*i+kk];

                    B0Tj_Dunit_comp_tmp0 += B0T[kk*6+2*j]*Dunit_comp[3*i+kk];
                    B0Tj_Dunit_comp_tmp1 += B0T[kk*6+2*j+1]*Dunit_comp[3*i+kk];
                }

                B0Tj_Dunit[2*i] = B0Tj_Dunit_tmp0;
                B0Tj_Dunit[2*i+1] = B0Tj_Dunit_tmp1;

                B0Tj_Dunit_comp[2*i] = B0Tj_Dunit_comp_tmp0;
                B0Tj_Dunit_comp[2*i+1] = B0Tj_Dunit_comp_tmp1;
            }

            for(int i=0; i<3; i++)
            {
                /* The rigidity matrix */
                /* scalar product of B0Ti_Dunit and the first column of B0T */
                B0Tj_Dunit_B0Ti_tmp0 = 0.;
                B0Tj_Dunit_B0Ti_tmp1 = 0.;
                B0Tj_Dunit_B0Ti_tmp2 = 0.;
                B0Tj_Dunit_B0Ti_tmp3 = 0.;

                /* scalar product of B0Ti_Dunit_comp and the first column of B0T */
                B0Tj_Dunit_comp_B0Ti_tmp0 = 0.;
                B0Tj_Dunit_comp_B0Ti_tmp1 = 0.;
                B0Tj_Dunit_comp_B0Ti_tmp2 = 0.;
                B0Tj_Dunit_comp_B0Ti_tmp3 = 0.;

                for(int kk=0; kk<3; kk++)
                {
                    B0Tj_Dunit_B0Ti_tmp0 += B0Tj_Dunit[2*kk]*B0T[kk*6+2*i];
                    B0Tj_Dunit_B0Ti_tmp1 += B0Tj_Dunit[2*kk]*B0T[kk*6+2*i+1];
                    B0Tj_Dunit_B0Ti_tmp2 += B0Tj_Dunit[2*kk+1]*B0T[kk*6+2*i];
                    B0Tj_Dunit_B0Ti_tmp3 += B0Tj_Dunit[2*kk+1]*B0T[kk*6+2*i+1];

                    B0Tj_Dunit_comp_B0Ti_tmp0 += B0Tj_Dunit_comp[2*kk]*B0T[kk*6+2*i];
                    B0Tj_Dunit_comp_B0Ti_tmp1 += B0Tj_Dunit_comp[2*kk]*B0T[kk*6+2*i+1];
                    B0Tj_Dunit_comp_B0Ti_tmp2 += B0Tj_Dunit_comp[2*kk+1]*B0T[kk*6+2*i];
                    B0Tj_Dunit_comp_B0Ti_tmp3 += B0Tj_Dunit_comp[2*kk+1]*B0T[kk*6+2*i+1];
                }

                B0T_Dunit_B0T[(2*i)*6+2*j] = B0Tj_Dunit_B0Ti_tmp0;
                B0T_Dunit_B0T[(2*i+1)*6+2*j] = B0Tj_Dunit_B0Ti_tmp1;
                B0T_Dunit_B0T[(2*i)*6+2*j+1] = B0Tj_Dunit_B0Ti_tmp2;
                B0T_Dunit_B0T[(2*i+1)*6+2*j+1] = B0Tj_Dunit_B0Ti_tmp3;

                B0T_Dunit_comp_B0T[(2*i)*6+2*j] = B0Tj_Dunit_comp_B0Ti_tmp0;
                B0T_Dunit_comp_B0T[(2*i+1)*6+2*j] = B0Tj_Dunit_comp_B0Ti_tmp1;
                B0T_Dunit_comp_B0T[(2*i)*6+2*j+1] = B0Tj_Dunit_comp_B0Ti_tmp2;
                B0T_Dunit_comp_B0T[(2*i+1)*6+2*j+1] = B0Tj_Dunit_comp_B0Ti_tmp3;
            }
        }

        /*
         * B0T_Dunit_B0T should be symmetric but is not exactly after the calcultation here above
         * because the sequence of operation is not the same for the component i,j and j,i.
         * We force the matrix to be symmetric by copying the upper part onto the lower part
         */
        for(int i=1; i<6; i++)
        {
            for(int j=0; j<i; j++)
            {
                B0T_Dunit_B0T[i*6+j] = B0T_Dunit_B0T[j*6+i];
                B0T_Dunit_comp_B0T[i*6+j] = B0T_Dunit_comp_B0T[j*6+i];
            }
        }

        /* Loop over the 6 by 6 components of the finite element intergrale
         * this is done smartly by looping over j=0:2 and i=0:2
         * col = (mwIndex)it[2*j]-1  , row = (mwIndex)it[2*i]-1;
         * col  , row   -> UU component
         * col  , row+1 -> VU component
         * col+1, row   -> VV component
         * col+1, row+1 -> UV component */

        std::vector<double> data(36);
        std::vector<double> fvdata(6);
        //int lc = 0;

        for(int j=0; j<3; j++)
        {
            /* Column corresponding to indice j (we also assemble terms in col+1) */
            //col = (mwIndex)it[2*j]-1; /* -1 to use the indice convention of C */

            for(int i=0; i<3; i++)
            {
                /* Row corresponding to indice i (we also assemble terms in row+1) */
                //row = (mwIndex)it[2*i]-1; /* -1 to use the indice convention of C */

                B0Tj_Dunit_B0Ti[0] = B0T_Dunit_B0T[(2*i)*6+2*j];
                B0Tj_Dunit_B0Ti[1] = B0T_Dunit_B0T[(2*i+1)*6+2*j];
                B0Tj_Dunit_B0Ti[2] = B0T_Dunit_B0T[(2*i)*6+2*j+1];
                B0Tj_Dunit_B0Ti[3] = B0T_Dunit_B0T[(2*i+1)*6+2*j+1];

                B0Tj_Dunit_comp_B0Ti[0] = B0T_Dunit_comp_B0T[(2*i)*6+2*j];
                B0Tj_Dunit_comp_B0Ti[1] = B0T_Dunit_comp_B0T[(2*i+1)*6+2*j];
                B0Tj_Dunit_comp_B0Ti[2] = B0T_Dunit_comp_B0T[(2*i)*6+2*j+1];
                B0Tj_Dunit_comp_B0Ti[3] = B0T_Dunit_comp_B0T[(2*i+1)*6+2*j+1];

                /* Select the nodal weight values from M_loc */
                mloc = Mass[3*j+i];

                B0Tj_sigma_h[0] = 0.;
                B0Tj_sigma_h[1] = 0.;

                for(int k=0; k<3; k++)
                {
                    B0Tj_sigma_h[0] += B0T[k*6+2*i]*(sigma[k]*h+sigma_P[k]);
                    B0Tj_sigma_h[1] += B0T[k*6+2*i+1]*(sigma[k]*h+sigma_P[k]);
                }

                /* ---------- UU component */
                double duu = surface_e*( mloc*(coef_Vair+coef_Voce*std::cos(ocean_turning_angle_rad)+coef_V+coef_basal)+B0Tj_Dunit_B0Ti[0]*coef*timestep+B0Tj_Dunit_comp_B0Ti[0]*coef_P);

                /* ---------- VU component */
                double dvu = surface_e*(+B0Tj_Dunit_B0Ti[1]*coef*timestep+B0Tj_Dunit_comp_B0Ti[1]*coef_P);

                /* ---------- UV component */
                double duv = surface_e*(+B0Tj_Dunit_B0Ti[2]*coef*timestep+B0Tj_Dunit_comp_B0Ti[2]*coef_P);

                /* ---------- VV component */
                double dvv = surface_e*( mloc*(coef_Vair+coef_Voce*std::cos(ocean_turning_angle_rad)+coef_V+coef_basal)+B0Tj_Dunit_B0Ti[3]*coef*timestep+B0Tj_Dunit_comp_B0Ti[3]*coef_P);

                // if (cpt ==1)
                // {
                //     std::cout<<"duu= "<< duu <<"\n";
                //     std::cout<<"dvu= "<< dvu <<"\n";
                //     std::cout<<"duv= "<< duv <<"\n";
                //     std::cout<<"dvv= "<< dvv <<"\n";
                // }

                data[12*i+2*j] = duu;
                data[12*i+2*j+1] = duv;
                data[(12*i+2*j)+6] = dvu;
                data[(12*i+2*j+1)+6] = dvv;

            }

            double fuu = surface_e*( mloc*( coef_Vair*Vair+coef_Voce*std::cos(ocean_turning_angle_rad)*Voce+coef_X+coef_V*VT) - B0Tj_sigma_h[0]/3);
            double fvv = surface_e*( mloc*( coef_Voce*std::sin(ocean_turning_angle_rad)*(Voce-VT)-coef_C*Vcor) );

            fvdata[2*j] = fuu;
            fvdata[2*j+1] = fvv;
        }


#if 1
        // if (cpt == 0)
        //     for (int i=0; i<6; ++i)
        //     {
        //         for (int j=0; j<6; ++j)
        //         {
        //             std::cout<< std::left << std::setw(12) << data[6*i+j] <<"  ";
        //         }
        //         std::cout<<"\n";
        //     }


        std::vector<int> rcindices(6);
        for (int s=0; s<3; ++s)
        {
            rcindices[2*s] = it->second.indices[s]-1;
            //rcindices[2*s+1] = it->second.indices[s];
            rcindices[2*s+1] = it->second.indices[s]-1+M_num_nodes;
        }

        // if (cpt == 0)
        //     for(int k=0; k<6; k++)
        //         std::cout<<"INDEX["<< rcindices[k] <<"]\n";

        M_matrix->addMatrix(&rcindices[0], rcindices.size(),
                            &rcindices[0], rcindices.size(), &data[0]);

        M_vector->addVector(&rcindices[0], rcindices.size(), &fvdata[0]);

#endif

        ++cpt;
    }

    M_matrix->close();

    std::cout<<"TIMER ASSEMBLY= " << chrono.elapsed() <<"s\n";

    chrono.restart();
    M_matrix->on(dirichlet_flags,*M_vector);
    std::cout<<"TIMER DBCA= " << chrono.elapsed() <<"s\n";

    std::cout<<"[PETSC MATRIX] CLOSED      = "<< M_matrix->closed() <<"\n";
    std::cout<<"[PETSC MATRIX] SIZE        = "<< M_matrix->size1() << " " << M_matrix->size2() <<"\n";
    std::cout<<"[PETSC MATRIX] SYMMETRIC   = "<< M_matrix->isSymmetric() <<"\n";
    std::cout<<"[PETSC MATRIX] NORM        = "<< M_matrix->linftyNorm() <<"\n";

    M_matrix->printMatlab("stiffness.m");
    M_vector->printMatlab("rhs.m");
}

} // Nextsim
