/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   finiteelement.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug 24 11:02:45 2015
 */

#include <finiteelement.hpp>
#include <constants.hpp>
#include <date.hpp>

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
    std::cout <<"GMSH VERSION= "<< M_mesh.version() <<"\n";
    M_mesh.setOrdering("bamg");
    M_mesh.readFromFile("bigarctic10km.msh");

    // createGMSHMesh("hypercube.geo");
    // //M_mesh.setOrdering("gmsh");
    // M_mesh.readFromFile("hypercube.msh");

    this->initBamg();

    //BamgConvertMeshx(bamgmesh,bamggeom,&index[0],&x[0],&y[0],M_num_nodes,M_num_elements);
    BamgConvertMeshx(
                     bamgmesh,bamggeom,
                     &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                     M_mesh.numNodes(), M_mesh.numTriangles()
                     );


    auto h = this->minMaxSide(M_mesh);
    bamgopt->hmin = h[0];
    bamgopt->hmax = h[1];

    double minang = this->minAngle(M_mesh);

    double resol = this->resolution(M_mesh);

    std::cout<<"Lenght min = "<< h[0] <<"\n";
    std::cout<<"Lenght max = "<< h[1] <<"\n";
    std::cout<<"Angle min  = "<< minang <<"\n";
    std::cout<<"Resolution = "<< resol <<"\n";

    if (minang < vm["simul_in.regrid_angle"].as<double>())
    {
        std::cout<<"invalid regridding angle: should be smaller than the minimal angle in the intial grid\n";
        throw std::logic_error("invalid regridding angle: should be smaller than the minimal angle in the intial grid");
    }

    auto hminVertices = this->hminVertices(M_mesh, bamgmesh);
    auto hmaxVertices = this->hmaxVertices(M_mesh, bamgmesh);

    bamgopt->hminVertices = &hminVertices[0];
    bamgopt->hmaxVertices = &hmaxVertices[0];

    if (!vm["simul_in.use_simul_out"].as<bool>())
    {
        //if (vm["simul_in.regrid"].as<std::string>() == "bamg")
        if (vm["simul_in.regrid"].as<std::string>() == "No-regridding")
        {
            if((bamgopt->hminVerticesSize[0] == 0) && (bamgopt->hminVerticesSize[1] == 0))
            {
                //bamgopt->KeepVertices=0;
                Bamgx(bamgmeshout,bamggeomout,bamgmesh,bamggeom,bamgopt);
            }
            else
            {
                bamgopt->KeepVertices=0;
                Bamgx(bamgmeshout,bamggeomout,bamgmesh,bamggeom,bamgopt);
            }
        }
    }

    this->importBamg();

    M_edges = M_mesh.edges();
    //int first = 0;
    for (auto it=M_edges.begin(), end=M_edges.end(); it!=end; ++it)
    {
        //if (it->second.physical==158)
        //if (it->second.physical==161)
        //if (it->second.physical==1)

        if (it->physical==0)
        {
            //++first;
            // if (first==1)
            //     continue;

            dirichlet_flags.push_back(it->indices[0]-1);
            dirichlet_flags.push_back(it->indices[0]-1+M_mesh.numNodes());
        }
    }

    // for (int flag=0; flag < 183; ++flag)
    // {
    //     std::cout<<"flag["<< flag <<"]= "<< dirichlet_flags[flag] <<"\n";
    // }

    std::cout<<"DIRICHLET_FLAGS= "<< dirichlet_flags.size() <<"\n";
    std::cout<<"NumNodes     = "<< M_mesh.numNodes() <<"\n";
    std::cout<<"NumTriangles = "<< M_mesh.numTriangles() <<"\n";
    std::cout<<"NumEdges     = "<< M_mesh.numEdges() <<"\n";

    M_mesh.stereographicProjection();
    //M_mesh.writeTofile("arctic10km.msh");

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_num_nodes = M_mesh.numNodes();

    M_matrix = matrix_ptrtype(new matrix_type(2*M_num_nodes,2*M_num_nodes,22));
    M_matrix->zero();

    M_vector = vector_ptrtype(new vector_type(2*M_num_nodes));
    M_vector->zero();
    //M_vector->setOnes();

    M_solution = vector_ptrtype(new vector_type(2*M_num_nodes));
    M_solution->zero();
}

void
FiniteElement::initBamg()
{
    bamgopt = new BamgOpts();

    bamgopt->Crack             = 0;
    bamgopt->anisomax          = 1e30;
    bamgopt->coeff             = 1;
    bamgopt->cutoff            = 1e-5;
    //bamgopt->err               = 0.01;
    bamgopt->errg              = 0.1;
    bamgopt->field             = NULL;
    bamgopt->gradation         = 1.5;
    bamgopt->Hessiantype       = 0;
    bamgopt->hmin              = 1e-100;
    bamgopt->hmax              = 1e100;
    bamgopt->hminVertices      = NULL;
    bamgopt->hmaxVertices      = NULL;
    bamgopt->hVertices         = NULL;
    bamgopt->KeepVertices      = 1;
    bamgopt->MaxCornerAngle    = 10;
    bamgopt->maxnbv            = 1e7;
    bamgopt->maxsubdiv         = 10;
    bamgopt->metric            = NULL;
    bamgopt->Metrictype        = 0;
    bamgopt->nbjacobi          = 1;
    bamgopt->nbsmooth          = 3;
    bamgopt->omega             = 1.8;
    bamgopt->power             = 1.;
    bamgopt->splitcorners      = 0; //the Devil!  Changed to 0, original 1 Phil
    bamgopt->geometricalmetric = 0;
    bamgopt->random            = true;
    bamgopt->verbose           = vm["verbose"].as<int>();

    bamggeom = new BamgGeom();
    bamgmesh = new BamgMesh();

    bamggeomout = new BamgGeom();
    bamgmeshout = new BamgMesh();

    bamgopt->Check();
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
FiniteElement::jacobian(element_type const& element, mesh_type const& mesh) const
{
    std::vector<double> vertex_0 = mesh.nodes()[element.indices[0]-1].coords;
    std::vector<double> vertex_1 = mesh.nodes()[element.indices[1]-1].coords;
    std::vector<double> vertex_2 = mesh.nodes()[element.indices[2]-1].coords;

    double jac = (vertex_1[0]-vertex_0[0])*(vertex_2[1]-vertex_0[1]);
    jac -= (vertex_2[0]-vertex_0[0])*(vertex_1[1]-vertex_0[1]);

    return  jac;
}

std::vector<double>
FiniteElement::sides(element_type const& element, mesh_type const& mesh) const
{
    std::vector<double> vertex_0 = mesh.nodes()[element.indices[0]-1].coords;
    std::vector<double> vertex_1 = mesh.nodes()[element.indices[1]-1].coords;
    std::vector<double> vertex_2 = mesh.nodes()[element.indices[2]-1].coords;

    // std::vector<double> vertex_0 = dispnodes[element.indices[0]-1].coords;
    // std::vector<double> vertex_1 = dispnodes[element.indices[1]-1].coords;
    // std::vector<double> vertex_2 = dispnodes[element.indices[2]-1].coords;

    std::vector<double> side(3);

    side[0] = std::sqrt(std::pow(vertex_1[0]-vertex_0[0],2.) + std::pow(vertex_1[1]-vertex_0[1],2.));
    side[1] = std::sqrt(std::pow(vertex_2[0]-vertex_1[0],2.) + std::pow(vertex_2[1]-vertex_1[1],2.));
    side[2] = std::sqrt(std::pow(vertex_2[0]-vertex_0[0],2.) + std::pow(vertex_2[1]-vertex_0[1],2.));

    return side;
}

std::vector<double>
FiniteElement::minMaxSide(mesh_type const& mesh) const
{
    std::vector<double> minmax(2);
    std::vector<double> all_min_side(mesh.numTriangles());
    std::vector<double> all_max_side(mesh.numTriangles());

    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        auto side = this->sides(*it,mesh);
        all_min_side[cpt] = *std::min_element(side.begin(),side.end());
        all_max_side[cpt] = *std::max_element(side.begin(),side.end());
        ++cpt;
    }

    minmax[0] = *std::min_element(all_min_side.begin(),all_min_side.end());
    minmax[1] = *std::max_element(all_max_side.begin(),all_max_side.end());

    //minmax[0] = std::accumulate(all_min_side.begin(),all_min_side.end(),0.)/(all_min_side.size());
    //minmax[1] = std::accumulate(all_max_side.begin(),all_max_side.end(),0.)/(all_max_side.size());

    return minmax;
}

double
FiniteElement::minAngles(element_type const& element, mesh_type const& mesh) const
{
    std::vector<double> side = this->sides(element,mesh);
    std::sort(side.begin(),side.end());
    double minang = std::acos( (std::pow(side[1],2.) + std::pow(side[2],2.) - std::pow(side[0],2.) )/(2*side[1]*side[2]) );
    minang = minang*45.0/std::atan(1.0);

    return minang;
}

double
FiniteElement::minAngle(mesh_type const& mesh) const
{
    std::vector<double> all_min_angle(mesh.numTriangles());

    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        all_min_angle[cpt] = this->minAngles(*it,mesh);
        ++cpt;
    }

    return *std::min_element(all_min_angle.begin(),all_min_angle.end());;
}

double
FiniteElement::resolution(mesh_type const& mesh) const
{
    std::vector<double> all_min_measure(mesh.numTriangles());

    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        all_min_measure[cpt] = this->measure(*it,mesh);
        ++cpt;
    }

    double resol = std::accumulate(all_min_measure.begin(),all_min_measure.end(),0.)/(all_min_measure.size());
    resol = std::pow(resol,0.5);

    return resol;
}


std::vector<double>
FiniteElement::hminVertices(mesh_type const& mesh, BamgMesh const* bamgmesh) const
{
    // std::cout<<"Connectivity[0]= "<< (bamgmesh->NodalElementConnectivitySize[0]) <<"\n";
    // std::cout<<"Connectivity[1]= "<< (bamgmesh->NodalElementConnectivitySize[1]) <<"\n";
    // std::cout<<"Connectivity[1]= "<< (mesh.numTriangles()) <<"\n";

    std::vector<double> hmin(bamgmesh->NodalElementConnectivitySize[0]);

    for (int i=0; i<bamgmesh->NodalElementConnectivitySize[0]; ++i)
    {
        std::vector<double> measure(bamgmesh->NodalElementConnectivitySize[1]);
        int j = 0;
        for (j=0; j<bamgmesh->NodalElementConnectivitySize[1]; ++j)
        {
            int elt_num = bamgmesh->NodalElementConnectivity[bamgmesh->NodalElementConnectivitySize[1]*i+j]-1;

            if ((0 <= elt_num) && (elt_num < mesh.numTriangles()) && (elt_num != NAN))
            {
                measure[j] = this->measure(mesh.triangles()[elt_num],mesh);
            }
            else
            {
                break;
            }
        }

        measure.resize(j);
        hmin[i] = std::sqrt(std::sqrt(2.)*(*std::min_element(measure.begin(),measure.end()))*0.9);
    }

    return hmin;
}

std::vector<double>
FiniteElement::hmaxVertices(mesh_type const& mesh, BamgMesh const* bamgmesh) const
{
    std::vector<double> hmax = this->hminVertices(mesh,bamgmesh);

    std::for_each(hmax.begin(), hmax.end(), [&](double& f){ f = 1.2*f; });

    return hmax;
}

double
FiniteElement::measure(element_type const& element, mesh_type const& mesh) const
{
    return (1./2)*std::abs(jacobian(element,mesh));
}

std::vector<double>
FiniteElement::shapeCoeff(element_type const& element, mesh_type const& mesh) const
{
    std::vector<double> x(3);
    std::vector<double> y(3);

    for (int i=0; i<3; ++i)
    {
        x[i] = mesh.nodes()[element.indices[i]-1].coords[0];
        y[i] = mesh.nodes()[element.indices[i]-1].coords[1];
    }

    std::vector<double> coeff(6);
    double jac = jacobian(element,mesh);

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

    std::vector<double> um(2*M_mesh.numNodes(),1.);
    for (int al=0; al<um.size();++al)
    {
        um[al] = static_cast <double> (std::rand()) / static_cast <double> (RAND_MAX);
    }

    chrono.restart();
    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        for (int i=0; i<18; ++i)
        {
            if (i < 3)
            {
                B0T[2*i] = this->shapeCoeff(*it,M_mesh)[i];
                B0T[12+2*i] = this->shapeCoeff(*it,M_mesh)[i+3];
                B0T[13+2*i] = this->shapeCoeff(*it,M_mesh)[i];
            }
            else if (i < 6)
            {
                B0T[2*i+1] = this->shapeCoeff(*it,M_mesh)[i];
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
        surface_e = this->measure(*it,M_mesh);

        // /* compute the x and y derivative of g*ssh */
        g_ssh_e_x = 0.;
        g_ssh_e_y = 0.;
        for(int i=0; i<3; i++)
        {
            g_ssh_e = g_ssh;   /* g*ssh at the node k of the element e */
            g_ssh_e_x += this->shapeCoeff(*it,M_mesh)[i]*g_ssh_e; /* x derivative of g*ssh */
            g_ssh_e_y += this->shapeCoeff(*it,M_mesh)[i+3]*g_ssh_e; /* y derivative of g*ssh */
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
            rcindices[2*s] = it->indices[s]-1;
            rcindices[2*s+1] = it->indices[s]-1+M_num_nodes;
        }

        // if (cpt == 0)
        //     for(int k=0; k<6; k++)
        //         std::cout<<"INDEX["<< rcindices[k] <<"]\n";

        M_matrix->addMatrix(&rcindices[0], rcindices.size(),
                            &rcindices[0], rcindices.size(), &data[0]);

        M_vector->addVector(&rcindices[0], rcindices.size(), &fvdata[0]);

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

    //M_matrix->printMatlab("stiffness.m");
    //M_vector->printMatlab("rhs.m");
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
              _ksp=vm["solver.ksp-type"].as<std::string>()/*"preonly"*/,
              _pc=vm["solver.pc-type"].as<std::string>()/*"cholesky"*/,
              _pcfactormatsolverpackage=vm["solver.mat-package-type"].as<std::string>()/*"cholmod"*/,
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

    double days_in_sec = 24.0*3600.0;
    int pcpt = 0;
    int niter = 0;
    double time_init = dateStr2Num(vm["simul_in.time_init"].as<std::string>());
    double output_time_step =  days_in_sec/vm["simul_in.output_per_day"].as<int>();
    double current_time = time_init /*+ pcpt*time_step/(24*3600.0)*/;

    double time_step = vm["simul_in.timestep"].as<double>();
    double duration = (vm["simul_in.duration"].as<double>())*days_in_sec;

    std::cout<<"TIMESTEP= "<< time_step <<"\n";
    std::cout<<"DURATION= "<< duration <<"\n";

    //gregorian::date epoch = date_time::parse_date<gregorian::date>( "01-01-1900", date_time::ymd_order_dmy);
    gregorian::date epoch = date_time::parse_date<gregorian::date>(
                                                                   vm["simul_in.time_init"].as<std::string>(),
                                                                   //date_time::ymd_order_dmy
                                                                   date_time::ymd_order_iso
                                                                   );

    std::string time_init_ym = to_iso_string(epoch).substr(0,6);
    std::string init_topaz_file = (boost::format( "TP4DAILY_%1%_3m.nc" ) % time_init_ym ).str();
    std::cout<<"INIT_TOPAZ_FILE "<< init_topaz_file <<"\n";

    std::string init_mit_file = (boost::format( "MITgcm_%1%_3m.nc" ) % time_init_ym ).str();
    std::cout<<"INIT_MIT_FILE "<< init_mit_file <<"\n";

    // Interpolation of the bathymetry
    if (vm["simul_in.Lemieux_basal_k2"].as<double>() > 0 )
    {
        std::vector<double> data_in(M_mesh_previous.numNodes(),1.);

        for (int i=0; i<data_in.size(); ++i)
            data_in[i] = static_cast <double> (std::rand()) / static_cast <double> (RAND_MAX);

        double *data_out;

        Options* options;
        options = new Options();

        InterpFromMeshToMesh2dx(&data_out,
                                &M_mesh_previous.indexTr()[0],&M_mesh_previous.coordX()[0],&M_mesh_previous.coordY()[0],
                                M_mesh_previous.numNodes(),M_mesh_previous.numTriangles(),
                                &data_in[0],
                                M_mesh_previous.numNodes(),1,
                                &M_mesh.coordX()[0],&M_mesh.coordY()[0],M_mesh.numNodes(),
                                options);
#if 0
        cout << "\n";
        cout << "     K      Xi(K)       Yi(K)       Zi(K)       Z(X,Y)\n";
        cout << "\n";
        //for (int k = 0; k < bamgmeshout->VerticesSize[0]; k++ )
        for (int k = 0; k < M_mesh.numNodes(); k++ )
        {
            //ze = xyi[0+k*2] + 2.0 * xyi[1+k*2];
            cout << "  " << setw(4) << k
                 << "  " << setw(10) << M_mesh.coordX()[k]
                 << "  " << setw(10) << M_mesh.coordY()[k]
                 << "  " << setw(10) << data_out[k] << "\n";
            //<< "  " << setw(10) << data_in[k] << "\n";
        }
#endif
    }

    this->forcingWind(0.,0.);
    this->forcingOcean(0.,0.);
    this->forcingThermo(0.,0.);

    this->initConcentration();
    this->initThickness();
    this->initDamage();

    for (int i=0; i<M_num_elements; ++i)
    {
        if ((conc[i] < 0) || (thick[i] < 0) )
        {
            conc[i] = 0;
            thick[i] = 0;
        }
    }

    // main loop for nextsim program

    bool is_running = true;
    while (is_running)
    {
        is_running = ((pcpt+1)*time_step) < duration;
        //is_running = false;

        current_time = time_init + pcpt*time_step/(24*3600.0);

        // step 0: preparation
        // remeshing and remapping of the prognostic variables

        if (vm["simul_in.regrid"].as<std::string>() == "No-regridding")
        {
            //std::cout<<"YES\n";
            //minAngleMesh(&M_mesh.coordX()[0],&M_mesh.coordY()[0],M_mesh.coordY().size());
        }

        //if (pcpt <30)
        //std::cout<<"CURRENT_TIME= "<< current_time <<"\n";

        ++pcpt;
    }

    std::cout<<"pcpt= "<< pcpt <<"\n";
    this->assemble();
    this->solve();

    //this->initSnowThickness();
    //this->performSimulation();
}

void
FiniteElement::error()
{
    double l2_error = 0;
    double sh1_error = 0;

    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        double area = measure(*it,M_mesh);
        std::vector<double> x(3);
        std::vector<double> y(3);

        double l2_contrib = 0;
        double sh1_contrib = 0;
        double entry_contrib = 0;

        for (int i=0; i<3; ++i)
        {
            // x[i] = M_nodes.find(it->second.indices[i])->second.coords[0];
            // y[i] = M_nodes.find(it->second.indices[i])->second.coords[1];

            x[i] = M_nodes[it->indices[i]-1].coords[0];
            y[i] = M_nodes[it->indices[i]-1].coords[1];
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
                sh1_contrib += entry_contrib*M_exact->operator()(it->indices[k]-1)/(4.0*area);

                // l2 error
                l2_contrib += M_exact->operator()(it->indices[k]-1)*((j == k) ? 2.0 : 1.0)*area/12.0;
            }

            l2_error += l2_contrib*M_exact->operator()(it->indices[j]-1);
            sh1_error += sh1_contrib*M_exact->operator()(it->indices[j]-1);
        }
    }

    std::cout<<"||u-uh||_L2  = "<< std::sqrt(l2_error) <<"\n";
    std::cout<<"||u-uh||_H1  = "<< std::sqrt(l2_error+sh1_error) <<"\n";
}

void
FiniteElement::performSimulation()
{
}

void
FiniteElement::initialConditions()
{

}

void
FiniteElement::forcingWind(double const& u, double const& v)
{
    switch (M_wind)
    {
        case forcing::WindType::CONSTANT:
            this->constantWind(u,v);
            break;

        default:
            std::cout << "invalid wind forcing"<<"\n";
            throw std::logic_error("invalid wind forcing");
    }
}

void
FiniteElement::constantWind(double const& u, double const& v)
{
    if (wind.size() ==0)
        wind.resize(2*M_num_nodes);

    for (int i=0; i<M_num_nodes; ++i)
    {
        wind[i] = u;
        wind[i+M_num_nodes] = v;
    }
}

void
FiniteElement::forcingOcean(double const& u, double const& v)
{
    switch (M_ocean)
    {
        case forcing::OceanType::CONSTANT:
            this->constantOcean(u,v);
            break;

        default:
            std::cout << "invalid ocean forcing"<<"\n";
            throw std::logic_error("invalid ocean forcing");
    }
}

void
FiniteElement::constantOcean(double const& u, double const& v)
{
    if (ocean.size() ==0)
        ocean.resize(2*M_num_nodes);

    for (int i=0; i<M_num_nodes; ++i)
    {
        ocean[i] = u;
        ocean[i+M_num_nodes] = v;
    }
}

void
FiniteElement::forcingThermo(double const& u, double const& v)
{
    switch (M_thermo)
    {
        case forcing::ThermoType::CONSTANT:
            this->constantThermo(u,v);
            break;

        default:
            std::cout << "invalid thermo forcing"<<"\n";
            throw std::logic_error("invalid thermo forcing");
    }
}

void
FiniteElement::constantThermo(double const& u, double const& v)
{
    if (thermo.size() ==0)
        thermo.resize(2*M_num_nodes);

    for (int i=0; i<M_num_nodes; ++i)
    {
        thermo[i] = u;
        thermo[i+M_num_nodes] = v;
    }
}

void
FiniteElement::initConcentration()
{
    switch (M_conc)
    {
        case forcing::ConcentrationType::CONSTANT:
            this->constantConc();
            break;

        default:
            std::cout << "invalid initialization of concentration"<<"\n";
            throw std::logic_error("invalid initialization of concentration");
    }
}

void
FiniteElement::constantConc()
{
    if (conc.size() ==0)
        conc.resize(M_num_elements);

    std::fill(conc.begin(), conc.end(), vm["simul_in.init_concentration"].as<double>());
}

void
FiniteElement::initThickness()
{
    switch (M_thick)
    {
        case forcing::ThicknessType::CONSTANT:
            this->constantThick();
            break;

        default:
            std::cout << "invalid initialization of thickness"<<"\n";
            throw std::logic_error("invalid initialization of thickness");
    }
}

void
FiniteElement::constantThick()
{
    if (thick.size() ==0)
        thick.resize(M_num_elements);

    for (int i=0; i<M_num_elements; ++i)
    {
        thick[i] = (vm["simul_in.init_thickness"].as<bool>())*conc[i];
    }
}

void
FiniteElement::initDamage()
{
    switch (M_damage)
    {
        case forcing::DamageType::CONSTANT:
            this->constantDamage();
            break;

        default:
            std::cout << "invalid initialization of damage"<<"\n";
            throw std::logic_error("invalid initialization of damage");
    }
}

void
FiniteElement::constantDamage()
{
    if (damage.size() ==0)
        damage.resize(M_num_elements);

    std::fill(damage.begin(), damage.end(), 0.);

#if 0
    for (int i=0; i<M_num_elements; ++i)
    {
        damage[i] = 1.0 - conc[i];
    }
#endif
}

void
FiniteElement::initSnowThickness()
{
    switch (M_snow_thick)
    {
        case forcing::SnowThicknessType::CONSTANT:
            this->constantSnowThick();
            break;

        default:
            std::cout << "invalid initialization of snow thickness"<<"\n";
            throw std::logic_error("invalid initialization of snow thickness");
    }
}

void
FiniteElement::constantSnowThick()
{
    if (snow_thick.size() ==0)
        snow_thick.resize(M_num_elements);

    for (int i=0; i<M_num_elements; ++i)
    {
        snow_thick[i] = (vm["simul_in.init_snow_thickness"].as<double>())*conc[i];
    }
}


void
FiniteElement::initThermodynamics()
{
    if (damage.size() ==0)
        damage.resize(M_num_elements);

    std::fill(damage.begin(), damage.end(), 0.);

#if 0
    for (int i=0; i<M_num_elements; ++i)
    {
        damage[i] = 1.0 - conc[i];
    }
#endif
}

void
FiniteElement::initDrifter()
{
    switch (M_drifter)
    {
        case forcing::DrifterType::EQUALLYSPACED:
            this->equallySpacedDrifter();
            break;

        default:
            std::cout << "invalid initialization of drifter"<<"\n";
            throw std::logic_error("invalid initialization of drifter");
    }
}

void
FiniteElement::equallySpacedDrifter()
{
    if (drifter.size() ==0)
        drifter.resize(M_num_elements);

    std::fill(drifter.begin(), drifter.end(), 0.);
}

void
FiniteElement::importBamg()
{
    //mesh_type mesh;
    std::vector<point_type> mesh_nodes;
    std::vector<element_type> mesh_edges;
    std::vector<element_type> mesh_triangles;
    std::vector<double> coords(3,0);

    mesh_nodes.resize(bamgmeshout->VerticesSize[0]);

    for (int id=0; id<bamgmeshout->VerticesSize[0]; ++id)
    {
        coords[0] = bamgmeshout->Vertices[3*id];
        coords[1] = bamgmeshout->Vertices[3*id+1];

        mesh_nodes[id].id = id+1;
        mesh_nodes[id].coords = coords;
    }

    int type = 2;
    int physical = 0;
    int elementary = 0;
    int numVertices = 2;
    std::vector<int> edges(numVertices);

    for (int edg=0; edg<bamgmeshout->EdgesSize[0]; ++edg)
    {
        edges[0] = bamgmeshout->Edges[3*edg];
        edges[1] = bamgmeshout->Edges[3*edg+1];
        //std::cout<< "Edges= "<< bamgmeshout->Edges[3*edg+2] <<"\n";

        element_type gmshElt( edg,
                              type,
                              physical,
                              elementary,
                              numVertices,
                              edges );

        //mesh_edges.insert(std::make_pair(edg,gmshElt));
        mesh_edges.push_back(gmshElt);
    }

    numVertices = 3;
    std::vector<int> indices(numVertices);

    for (int tr=0; tr<bamgmeshout->TrianglesSize[0]; ++tr)
    {
        indices[0] = bamgmeshout->Triangles[4*tr];
        indices[1] = bamgmeshout->Triangles[4*tr+1];
        indices[2] = bamgmeshout->Triangles[4*tr+2];

        element_type gmshElt( tr,
                              type,
                              physical,
                              elementary,
                              numVertices,
                              indices );

        //mesh_triangles.insert(std::make_pair(tr,gmshElt));
        mesh_triangles.push_back(gmshElt);
    }

    M_mesh_previous = M_mesh;
    M_mesh = mesh_type(mesh_nodes,mesh_edges,mesh_triangles);
    M_mesh.writeTofile("out.msh");

    std::cout<<"\n";
    std::cout<<"INFO: Previous  NumNodes     = "<< M_mesh_previous.numNodes() <<"\n";
    std::cout<<"INFO: Previous  NumTriangles = "<< M_mesh_previous.numTriangles() <<"\n";
    std::cout<<"INFO: Previous  NumEdges     = "<< M_mesh_previous.numEdges() <<"\n";
    std::cout<<"\n";
    std::cout<<"INFO: Current  NumNodes      = "<< M_mesh.numNodes() <<"\n";
    std::cout<<"INFO: Current  NumTriangles  = "<< M_mesh.numTriangles() <<"\n";
    std::cout<<"INFO: Current  NumEdges      = "<< M_mesh.numEdges() <<"\n";
    std::cout<<"\n";
}

} // Nextsim
