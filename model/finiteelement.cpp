/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   finiteelement.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug 24 11:02:45 2015
 */

#include <finiteelement.hpp>
#include <constants.hpp>
#include <date.hpp>
#include <thin_ice_redistribute.hpp>

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

    M_mesh.stereographicProjection();
    // M_mesh.writeTofile("arctic10km.msh");

    // createGMSHMesh("hypercube.geo");
    // //M_mesh.setOrdering("gmsh");
    // M_mesh.readFromFile("hypercube.msh");

    this->initConstant();

    this->initBamg();

    BamgConvertMeshx(
                     bamgmesh,bamggeom,
                     &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                     M_mesh.numNodes(), M_mesh.numTriangles()
                     );

    _printf_("VerticesSize[0] = " << bamgmesh->VerticesSize[0] << "\n");
    _printf_("EdgesSize[0] = " << bamgmesh->EdgesSize[0] << "\n");
    _printf_("rianglesSize[0] = " << bamgmesh->TrianglesSize[0] << "\n");

    _printf_("SubDomainsSize[0] = " << bamgmesh->SubDomainsSize[0] << "\n");
    _printf_("SubDomainsFromGeomSize[0] = " << bamgmesh->SubDomainsFromGeomSize[0] << "\n");

    _printf_("VerticesOnGeomVertexSize[0] = " << bamgmesh->VerticesOnGeomVertexSize[0] << "\n");
    _printf_("VerticesOnGeomEdgeSize[0] = " << bamgmesh->VerticesOnGeomEdgeSize[0] << "\n");
    _printf_("EdgesOnGeomEdgeSize[0] = " << bamgmesh->EdgesOnGeomEdgeSize[0] << "\n");

    _printf_("IssmEdgesSize[0] = " << bamgmesh->IssmEdgesSize[0] << "\n");
    _printf_("IssmSegmentsSize[0] = " << bamgmesh->IssmSegmentsSize[0] << "\n");

    _printf_("ElementConnectivitySize[1] = " << bamgmesh->ElementConnectivitySize[1] << "\n");
    _printf_("NodalConnectivitySize[1] = " << bamgmesh->NodalConnectivitySize[1] << "\n");
    _printf_("NodalElementConnectivitySize[1] = " << bamgmesh->NodalElementConnectivitySize[1] << "\n");

    _printf_("bamggeom->VerticesSize[0] = " << bamggeom->VerticesSize[0] << "\n");
    _printf_("bamggeom->EdgesSize[0] = " << bamggeom->EdgesSize[0] << "\n");
    _printf_("bamggeom->RequiredVerticesSize[0] = " << bamggeom->RequiredVerticesSize[0] << "\n");



    /*for (int i=0; i<bamgmesh->VerticesOnGeomVertexSize[0]; ++i)
    {
        _printf_( bamgmesh->VerticesOnGeomVertex[i*2] << " " << bamgmesh->VerticesOnGeomVertex[i*2+1] <<"\n");
    }

    for (int i=0; i<bamggeom->RequiredVerticesSize[0]; ++i)
    {
        _printf_( bamggeom->RequiredVertices[i]  <<"\n");
    }*/

    for (auto it=M_mesh.edges().begin(), end=M_mesh.edges().end(); it!=end; ++it)
    {
        if (it->physical==161)
        {
            M_dirichlet_flags.push_back(it->indices[0]-1);
            M_dirichlet_flags.push_back(it->indices[1]-1);
        }
    }

    std::sort(M_dirichlet_flags.begin(), M_dirichlet_flags.end());
    M_dirichlet_flags.erase(std::unique( M_dirichlet_flags.begin(), M_dirichlet_flags.end() ), M_dirichlet_flags.end());

    importBamg(bamgmesh);

    M_mesh_init = M_mesh;

    M_edges = M_mesh.edges();

    auto h = this->minMaxSide(M_mesh);
    bamgopt->hmin = h[0];
    bamgopt->hmax = h[1];

    std::cout<<"HMIN= "<< h[0] <<"\n";
    std::cout<<"HMAX= "<< h[1] <<"\n";
    std::cout<<"RES = "<< this->resolution(M_mesh) <<"\n";

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_num_nodes = M_mesh.numNodes();

    M_matrix = matrix_ptrtype(new matrix_type());
    M_vector = vector_ptrtype(new vector_type());
    M_solution = vector_ptrtype(new vector_type());


    const boost::unordered_map<const std::string, forcing::WindType> str2wind = boost::assign::map_list_of
        ("constant", forcing::WindType::CONSTANT)
        ("asr", forcing::WindType::ASR);
    M_wind_type = str2wind.find(vm["wind-type"].as<std::string>())->second;

    //std::cout<<"WINDTYPE= "<< (int)M_wind_type <<"\n";

    const boost::unordered_map<const std::string, forcing::OceanType> str2ocean = boost::assign::map_list_of
        ("constant", forcing::OceanType::CONSTANT)
        ("topaz", forcing::OceanType::TOPAZR);
    M_ocean_type = str2ocean.find(vm["ocean-type"].as<std::string>())->second;

    //std::cout<<"OCEANTYPE= "<< (int)M_ocean_type <<"\n";

    const boost::unordered_map<const std::string, forcing::ConcentrationType> str2conc = boost::assign::map_list_of
        ("constant", forcing::ConcentrationType::CONSTANT)
        ("topaz", forcing::ConcentrationType::TOPAZ4);
    M_conc_type = str2conc.find(vm["concentration-type"].as<std::string>())->second;

    //std::cout<<"CONCTYPE= "<< (int)M_conc_type <<"\n";

    const boost::unordered_map<const std::string, forcing::ThicknessType> str2thick = boost::assign::map_list_of
        ("constant", forcing::ThicknessType::CONSTANT)
        ("topaz", forcing::ThicknessType::TOPAZ4);
    M_thick_type = str2thick.find(vm["thickness-type"].as<std::string>())->second;

    //std::cout<<"THICKTYPE= "<< (int)M_thick_type <<"\n";

    // init options for interpolation from mesh to mesh
    // options = new Options();
}

void
FiniteElement::initSimulation()
{

    M_matrix->init(2*M_num_nodes,2*M_num_nodes,22);
    M_vector->init(2*M_num_nodes);
    M_solution->init(2*M_num_nodes);

    M_VT.resize(2*M_num_nodes,0.);
    M_VTM.resize(2*M_num_nodes,0.);
    M_VTMM.resize(2*M_num_nodes,0.);

    M_wind.resize(2*M_num_nodes);
    M_ocean.resize(2*M_num_nodes);
    M_thermo.resize(2*M_num_nodes);
    M_bathy_depth.resize(M_mesh_init.numTriangles(),200.);

    M_vair.resize(2);
    M_voce.resize(2);
    M_vssh.resize(M_num_nodes);

    M_UM.resize(2*M_num_nodes,0.);

    M_element_depth.resize(M_num_elements,0.);

    M_h_thin.resize(M_num_elements,0.);
    M_hs_thin.resize(M_num_elements,0.);

    M_h_ridged_thin_ice.resize(M_num_elements,0.);
    M_h_ridged_thick_ice.resize(M_num_elements,0.);

    M_divergence_rate.resize(M_num_elements,0.);
    M_sigma.resize(3*M_num_elements,0.);

    M_random_number.resize(M_num_elements);
    for (int i=0; i<M_random_number.size(); ++i)
        M_random_number[i] = static_cast <double> (std::rand()) / static_cast <double> (RAND_MAX);


    M_conc.resize(M_num_elements);
    M_thick.resize(M_num_elements);
    M_damage.resize(M_num_elements);
    M_snow_thick.resize(M_num_elements);

    this->initConcentration();
    this->initThickness();
    this->initDamage();
    this->initSnowThickness();

    for (int i=0; i<M_num_elements; ++i)
    {
        if ((M_conc[i] <= 0.) || (M_thick[i] <= 0.) )
        {
            M_conc[i] = 0.;
            M_thick[i] = 0.;
        }
    }

    M_norm_Voce_ice.resize(M_num_elements);
    M_norm_Vair_ice.resize(M_num_elements);
    M_norm_Vice.resize(M_num_elements);
    M_element_ssh.resize(M_num_elements,0.);
    M_ssh.resize(M_num_nodes,0.);

    M_Vair_factor.resize(M_num_elements);
    M_Voce_factor.resize(M_num_elements);

    M_basal_factor.resize(M_num_elements);
    M_fcor.resize(M_num_elements);
    M_Vcor.resize(M_VT.size());

    M_ftime_wind_range.resize(2,0.);
    M_ftime_ocean_range.resize(2,0.);
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

    bamgopt_previous = new BamgOpts();
    bamggeom_previous = new BamgGeom();
    bamgmesh_previous = new BamgMesh();

    bamgopt->Check();
}

void
FiniteElement::initConstant()
{
    nu0 = vm["simul_in.nu0"].as<double>();
    young = vm["simul_in.young"].as<double>();
    rhoi = physical::rhoi;
    rhos = physical::rhos;

    days_in_sec = 24.0*3600.0;
    time_init = dateStr2Num(vm["simul_in.time_init"].as<std::string>());
    output_time_step =  days_in_sec/vm["simul_in.output_per_day"].as<int>();

    time_step = vm["simul_in.timestep"].as<double>();
    duration = (vm["simul_in.duration"].as<double>())*days_in_sec;
    spinup_duration = (vm["simul_in.spinup_duration"].as<double>())*days_in_sec;

    divergence_min = (1./days_in_sec)*vm["simul_in.divergence_min"].as<double>();
    compression_factor = vm["simul_in.compression_factor"].as<double>();
    exponent_compression_factor = vm["simul_in.exponent_compression_factor"].as<double>();
    ocean_turning_angle_rad = (PI/180.)*vm["simul_in.oceanic_turning_angle"].as<double>();
    ridging_exponent = vm["simul_in.ridging_exponent"].as<double>();

    quad_drag_coef_air = vm["simul_in.ASR_quad_drag_coef_air"].as<double>();
    quad_drag_coef_water = vm["simul_in.quad_drag_coef_water"].as<double>();

    basal_k2 = vm["simul_in.Lemieux_basal_k2"].as<double>();
    basal_drag_coef_air = vm["simul_in.Lemieux_drag_coef_air"].as<double>();
    basal_u_0 = vm["simul_in.Lemieux_basal_u_0"].as<double>();
    basal_Cb = vm["simul_in.Lemieux_basal_Cb"].as<double>();

    time_relaxation_damage = vm["simul_in.time_relaxation_damage"].as<double>();

    h_thin_max = vm["simul_in.h_thin_max"].as<double>();
    c_thin_max = vm["simul_in.c_thin_max"].as<double>();

    compr_strength = vm["simul_in.compr_strength"].as<double>();
    tract_coef = vm["simul_in.tract_coef"].as<double>();
    scale_coef = vm["simul_in.scale_coef"].as<double>();
    alea_factor = vm["simul_in.alea_factor"].as<double>();
    cfix = vm["simul_in.cfix"].as<double>();

    C_fix    = cfix*scale_coef;          // C_fix;...  : cohesion (mohr-coulomb) in MPa (40000 Pa)
    C_alea   = alea_factor*C_fix;        // C_alea;... : alea sur la cohesion (Pa)
    tan_phi = vm["simul_in.tan_phi"].as<double>();
    ridge_h = vm["simul_in.ridge_h"].as<double>();
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

double
FiniteElement::jacobian(element_type const& element, mesh_type const& mesh,
                        std::vector<double> const& um, double factor) const
{
    std::vector<double> vertex_0 = mesh.nodes()[element.indices[0]-1].coords;
    std::vector<double> vertex_1 = mesh.nodes()[element.indices[1]-1].coords;
    std::vector<double> vertex_2 = mesh.nodes()[element.indices[2]-1].coords;

    for (int i=0; i<2; ++i)
    {
        vertex_0[i] += factor*um[element.indices[0]-1+i*(M_num_nodes)];
        vertex_1[i] += factor*um[element.indices[1]-1+i*(M_num_nodes)];
        vertex_2[i] += factor*um[element.indices[2]-1+i*(M_num_nodes)];
    }

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

    // minmax[0] = *std::min_element(all_min_side.begin(),all_min_side.end());
    // minmax[1] = *std::max_element(all_max_side.begin(),all_max_side.end());

    minmax[0] = std::accumulate(all_min_side.begin(),all_min_side.end(),0.)/(all_min_side.size());
    minmax[1] = std::accumulate(all_max_side.begin(),all_max_side.end(),0.)/(all_max_side.size());

    return minmax;
}

double
FiniteElement::minAngles(element_type const& element, mesh_type const& mesh) const
{
    std::vector<double> side = this->sides(element,mesh);
    //std::for_each(side.begin(), side.end(), [&](double& f){ f = 1000.*f; });
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

    return *std::min_element(all_min_angle.begin(),all_min_angle.end());
}

double
FiniteElement::minAngle(mesh_type const& mesh, std::vector<double> const& um, double factor) const
{
    auto movedmesh = mesh;
    movedmesh.move(um,factor);

    std::vector<double> all_min_angle(movedmesh.numTriangles());

    int cpt = 0;
    for (auto it=movedmesh.triangles().begin(), end=movedmesh.triangles().end(); it!=end; ++it)
    {
        all_min_angle[cpt] = this->minAngles(*it,movedmesh);
        ++cpt;
    }

    return *std::min_element(all_min_angle.begin(),all_min_angle.end());;
}

bool
FiniteElement::flip(mesh_type const& mesh, std::vector<double> const& um, double factor) const
{
    auto movedmesh = mesh;
    movedmesh.move(um,factor);

    std::vector<double> area(movedmesh.numTriangles());

    int cpt = 0;
    for (auto it=movedmesh.triangles().begin(), end=movedmesh.triangles().end(); it!=end; ++it)
    {
        area[cpt] = this->jacobian(*it,movedmesh);
        ++cpt;
    }

    double minarea = *std::min_element(area.begin(),area.end());
    double maxarea = *std::max_element(area.begin(),area.end());

    return ((minarea <= 0.) && (maxarea >= 0.));
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
FiniteElement::hminVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh) const
{
    std::vector<double> hmin(bamg_mesh->NodalElementConnectivitySize[0]);

    for (int i=0; i<bamg_mesh->NodalElementConnectivitySize[0]; ++i)
    {
        std::vector<double> measure(bamg_mesh->NodalElementConnectivitySize[1]);
        int j = 0;
        for (j=0; j<bamg_mesh->NodalElementConnectivitySize[1]; ++j)
        {
            int elt_num = bamg_mesh->NodalElementConnectivity[bamg_mesh->NodalElementConnectivitySize[1]*i+j]-1;

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
FiniteElement::hmaxVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh) const
{
    std::vector<double> hmax = this->hminVertices(mesh,bamg_mesh);

    std::for_each(hmax.begin(), hmax.end(), [&](double& f){ f = 1.2*f; });

    return hmax;
}

double
FiniteElement::measure(element_type const& element, mesh_type const& mesh) const
{
    return (1./2)*std::abs(jacobian(element,mesh));
}

double
FiniteElement::measure(element_type const& element, mesh_type const& mesh,
                       std::vector<double> const& um, double factor) const
{
    return (1./2)*std::abs(jacobian(element,mesh,um,factor));
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
FiniteElement::regrid(bool step)
{
    double displacement_factor = 2.;
    bool flip = true;
    int substep = 0;

    double* hmin_vertices;
    double* hmax_vertices;

    std::vector<double> hmin_vertices_first;
    std::vector<double> hmax_vertices_first;

    if (step)
    {
        while (flip)
        {
            ++substep;
            displacement_factor /= 2.;
            flip = this->flip(M_mesh,M_UM,displacement_factor);

            if (substep > 1)
                std::cout<<"FLIP DETECTED "<< substep-1 <<"\n";
        }

        std::cout<<"displacement_factor= "<< displacement_factor <<"\n";

        M_mesh.move(M_UM,displacement_factor);

#if 0
        chrono.restart();
        InterpFromMeshToMesh2dx(&hmin_vertices,
                                &M_mesh_init.indexTr()[0],&M_mesh_init.coordX()[0],&M_mesh_init.coordY()[0],
                                M_mesh_init.numNodes(),M_mesh_init.numTriangles(),
                                bamgopt->hminVertices,
                                M_mesh_init.numNodes(),1,
                                &M_mesh.coordX()[0],&M_mesh.coordY()[0],M_mesh.numNodes(),
                                false);

        InterpFromMeshToMesh2dx(&hmax_vertices,
                                &M_mesh_init.indexTr()[0],&M_mesh_init.coordX()[0],&M_mesh_init.coordY()[0],
                                M_mesh_init.numNodes(),M_mesh_init.numTriangles(),
                                bamgopt->hmaxVertices,
                                M_mesh_init.numNodes(),1,
                                &M_mesh.coordX()[0],&M_mesh.coordY()[0],M_mesh.numNodes(),
                                false);

        std::cout<<"TIMER INTERPOLATION= " << chrono.elapsed() <<"s\n";

#endif

        //auto hminVertices = this->hminVertices(M_mesh_init, bamgmesh);
        //auto hmaxVertices = this->hmaxVertices(M_mesh_init, bamgmesh);

#if 0
        cout << "\n";
        cout << "     K      Xi(K)       Yi(K)       Zi(K)       Z(X,Y)\n";
        cout << "\n";
        //for (int k = 0; k < bamgmeshout->VerticesSize[0]; k++ )
        for (int k = 0; k < M_mesh_init.numNodes(); k++ )
        {
            //ze = xyi[0+k*2] + 2.0 * xyi[1+k*2];
            cout << "  " << setw(4) << k
                 << "  " << setw(10) << M_mesh_init.coordX()[k]
                 << "  " << setw(10) << M_mesh_init.coordY()[k]
                //<< "  " << setw(10) << bamgopt->hminVertices[k]
                //<< "  " << setw(10) << bamgopt->hmaxVertices[k] << "\n";
                 << "  " << setw(10) << hmin_vertices[k]
                 << "  " << setw(10) << hmax_vertices[k] << "\n";
            //<< "  " << setw(10) << data_in[k] << "\n";
        }
#endif
        // if (bamgopt->KeepVertices!=0)
        //     bamgopt->KeepVertices=0;
    }
#if 0
    else
    {
        hmin_vertices_first = this->hminVertices(M_mesh, bamgmesh);
        hmax_vertices_first = this->hmaxVertices(M_mesh, bamgmesh);
    }

    bamgopt->hminVertices = new double[M_mesh.numNodes()];
    bamgopt->hmaxVertices = new double[M_mesh.numNodes()];

    for (int i=0; i<M_mesh.numNodes();++i)
    {
        bamgopt->hminVertices[i] = (step) ? hmin_vertices[i] : hmin_vertices_first[i];
        bamgopt->hmaxVertices[i] = (step) ? hmax_vertices[i] : hmax_vertices_first[i];
    }
#endif

    BamgConvertMeshx(
                     bamgmesh_previous,bamggeom_previous,
                     &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                     M_mesh.numNodes(), M_mesh.numTriangles()
                     );

    std::vector<int> vec;
    int fnd = 0;
    int snd = 0;

#if 1
    for (int edg=0; edg<bamgmesh_previous->EdgesSize[0]; ++edg)
    {
        fnd = bamgmesh_previous->Edges[3*edg]-1;
        snd = bamgmesh_previous->Edges[3*edg+1]-1;

        if ((std::binary_search(M_dirichlet_flags.begin(),M_dirichlet_flags.end(),fnd))
            || (std::binary_search(M_dirichlet_flags.begin(),M_dirichlet_flags.end(),snd)))
        {
            bamggeom_previous->Edges[3*edg+2] = 161;
            bamgmesh_previous->Edges[3*edg+2] = 161;

        }
    }

    //std::sort(vec.begin(), vec.end());
    // for (int i=0; i<vec.size(); ++i)
    // {
    //     //std::cout<<"VEC["<< i <<"]= "<< vec[i] <<"\n";
    //     //std::cout<<"INIT["<< vec[i] << "]= ("<< M_mesh.coordX()[vec[i]] <<","<< M_mesh.coordY()[vec[i]] <<")\n";
    // }
#endif

    *bamgopt_previous = *bamgopt;
    //Bamgx(bamgmesh,bamggeom_previous,bamgmesh_previous,bamggeom_previous,bamgopt_previous);
    Bamgx(bamgmesh,bamggeom,bamgmesh_previous,bamggeom_previous,bamgopt_previous);
    this->importBamg(bamgmesh);

    // update dirichlet nodes
    M_boundary_flags.resize(0);
    M_dirichlet_flags.resize(0);
    for (int edg=0; edg<bamgmesh->EdgesSize[0]; ++edg)
    {
        M_boundary_flags.push_back(bamgmesh->Edges[3*edg]-1);
        M_boundary_flags.push_back(bamgmesh->Edges[3*edg+1]-1);

        if (bamgmesh->Edges[3*edg+2] == 161)
        {
            M_dirichlet_flags.push_back(bamgmesh->Edges[3*edg]-1);
            M_dirichlet_flags.push_back(bamgmesh->Edges[3*edg+1]-1);

            //std::cout<<"NODES["<< edg <<"]= "<< bamgmesh->Edges[3*edg]-1 << " and "<< bamgmesh->Edges[3*edg+1]-1 <<"\n";
        }
    }

    std::sort(M_dirichlet_flags.begin(), M_dirichlet_flags.end());
    M_dirichlet_flags.erase( std::unique(M_dirichlet_flags.begin(), M_dirichlet_flags.end() ), M_dirichlet_flags.end());

    std::sort(M_boundary_flags.begin(), M_boundary_flags.end());
    M_boundary_flags.erase( std::unique(M_boundary_flags.begin(), M_boundary_flags.end() ), M_boundary_flags.end());

    M_neumann_flags.resize(0);
    std::set_difference(M_boundary_flags.begin(), M_boundary_flags.end(),
                        M_dirichlet_flags.begin(), M_dirichlet_flags.end(),
                        std::back_inserter(M_neumann_flags));

    // for (const int& edg : M_dirichlet_flags)
    // {
    //     std::cout<<"AFTER["<< edg << "]= ("<< M_mesh.coordX()[edg] <<","<< M_mesh.coordY()[edg] <<")\n";
    // }

    M_dirichlet_nodes.resize(2*(M_dirichlet_flags.size()));
    for (int i=0; i<M_dirichlet_flags.size(); ++i)
    {
        M_dirichlet_nodes[2*i] = M_dirichlet_flags[i];
        M_dirichlet_nodes[2*i+1] = M_dirichlet_flags[i]+M_num_nodes;
    }


    M_neumann_nodes.resize(2*(M_neumann_flags.size()));
    for (int i=0; i<M_neumann_flags.size(); ++i)
    {
        M_neumann_nodes[2*i] = M_neumann_flags[i];
        M_neumann_nodes[2*i+1] = M_neumann_flags[i]+M_num_nodes;
    }

    if (step)
    {
        // NODAL INTERPOLATION
        int prv_num_nodes = M_mesh_previous.numNodes();
        double* interp_in;
        interp_in = new double[6*prv_num_nodes];

        double* interp_out;

        std::cout<<"NODAL: Interp starts\n";
        chrono.restart();

        for (int i=0; i<prv_num_nodes; ++i)
        {
            // VT
            interp_in[6*i] = M_VT[i];
            interp_in[6*i+1] = M_VT[i+prv_num_nodes];

            // VTM
            interp_in[6*i+2] = M_VTM[i];
            interp_in[6*i+3] = M_VTM[i+prv_num_nodes];

            // VTMM
            interp_in[6*i+4] = M_VTMM[i];
            interp_in[6*i+5] = M_VTMM[i+prv_num_nodes];
        }

        InterpFromMeshToMesh2dx(&interp_out,
                                &M_mesh_previous.indexTr()[0],&M_mesh_previous.coordX()[0],&M_mesh_previous.coordY()[0],
                                M_mesh_previous.numNodes(),M_mesh_previous.numTriangles(),
                                interp_in,
                                M_mesh_previous.numNodes(),6,
                                &M_mesh.coordX()[0],&M_mesh.coordY()[0],M_mesh.numNodes(),
                                false);

        M_VT.resize(2*M_num_nodes,0.);
        M_VTM.resize(2*M_num_nodes,0.);
        M_VTMM.resize(2*M_num_nodes,0.);
        M_UM.resize(2*M_num_nodes,0.);

        for (int i=0; i<M_num_nodes; ++i)
        {
            // VT
            M_VT[i] = interp_out[6*i];
            M_VT[i+M_num_nodes] = interp_out[6*i+1];

            // VTM
            M_VTM[i] = interp_out[6*i+2];
            M_VTM[i+M_num_nodes] = interp_out[6*i+3];

            // VTMM
            M_VTMM[i] = interp_out[6*i+4];
            M_VTMM[i+M_num_nodes] = interp_out[6*i+5];

            // UM
            M_UM[i] = 0.;
            M_UM[i+M_num_nodes] = 0.;
        }

        std::cout<<"NODAL: Interp done\n";
        std::cout<<"NODAL: TIMER INTERPOLATION= " << chrono.elapsed() <<"s\n";

        // ELEMENT INTERPOLATION
        int prv_num_elements = M_mesh_previous.numTriangles();
        double* interp_elt_in;
        interp_elt_in = new double[11*prv_num_elements];

        double* interp_elt_out;

        std::cout<<"ELEMENT: Interp starts\n";
        chrono.restart();

        for (int i=0; i<prv_num_elements; ++i)
        {
            // concentration
            interp_elt_in[11*i] = M_conc[i];

            // thickness
            interp_elt_in[11*i+1] = M_thick[i];

            // snow thickness
            interp_elt_in[11*i+2] = M_snow_thick[i];

            // integrated_stress1
            interp_elt_in[11*i+3] = M_sigma[3*i]*M_thick[i];

            // integrated_stress2
            interp_elt_in[11*i+4] = M_sigma[3*i+1]*M_thick[i];

            // integrated_stress3
            interp_elt_in[11*i+5] = M_sigma[3*i+2]*M_thick[i];

            // compliance
            interp_elt_in[11*i+6] = 1./(1.-M_damage[i]);

            // divergence_rate
            interp_elt_in[11*i+7] = M_divergence_rate[i];

            // h_ridged_thin_ice
            interp_elt_in[11*i+8] = M_h_ridged_thin_ice[i];

            // h_ridged_thick_ice
            interp_elt_in[11*i+9] = M_h_ridged_thick_ice[i];

            // random_number
            interp_elt_in[11*i+10] = M_random_number[i];

        }

        double* surface_previous = new double[prv_num_elements];
        double* surface = new double[M_num_elements];

        int cpt = 0;
        for (auto it=M_mesh_previous.triangles().begin(), end=M_mesh_previous.triangles().end(); it!=end; ++it)
        {
            surface_previous[cpt] = this->measure(*it,M_mesh_previous);
            ++cpt;
        }

        cpt = 0;
        for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
        {
            surface[cpt] = this->measure(*it,M_mesh);
            ++cpt;
        }

        //#if 0
        InterpFromMeshToMesh2dCavities(&interp_elt_out,interp_elt_in,11,
             surface_previous, surface, bamgmesh_previous, bamgmesh);
        //#endif

        #if 0
        InterpFromMeshToMesh2dx(&interp_elt_out,
                                &M_mesh_previous.indexTr()[0],&M_mesh_previous.coordX()[0],&M_mesh_previous.coordY()[0],
                                M_mesh_previous.numNodes(),M_mesh_previous.numTriangles(),
                                interp_elt_in,
                                M_mesh_previous.numTriangles(),11,
                                &M_mesh.bCoordX()[0],&M_mesh.bCoordY()[0],M_mesh.numTriangles(),
                                false);
        #endif

        M_conc.resize(M_num_elements,0.);
        M_thick.resize(M_num_elements,0.);
        M_snow_thick.resize(M_num_elements,0.);
        M_sigma.resize(3*M_num_elements,0.);
        M_damage.resize(M_num_elements,0.);

        M_divergence_rate.resize(M_num_elements,0.);
        M_h_ridged_thin_ice.resize(M_num_elements,0.);
        M_h_ridged_thick_ice.resize(M_num_elements,0.);

        M_random_number.resize(M_num_elements);

        for (int i=0; i<M_num_elements; ++i)
        {
            // concentration
            M_conc[i] = std::max(0., std::min(1.,interp_elt_out[11*i]));

            // thickness
            M_thick[i] = std::max(0., std::min(1.,interp_elt_out[11*i+1]));

            // snow thickness
            M_snow_thick[i] = std::max(0., std::min(1.,interp_elt_out[11*i+2]));

            if (M_thick[i] != 0.)
            {
                // integrated_stress1
                M_sigma[3*i] = interp_elt_out[11*i+3]/M_thick[i];

                // integrated_stress2
                M_sigma[3*i+1] = interp_elt_out[11*i+4]/M_thick[i];

                // integrated_stress3
                M_sigma[3*i+2] = interp_elt_out[11*i+5]/M_thick[i];
            }

            // compliance
            if (interp_elt_out[11*i+6] != 0.)
            {
                M_damage[i] = std::max(0., std::min(1.,1./(1.-interp_elt_out[11*i+6])));
            }

            // divergence_rate
            M_divergence_rate[i] = interp_elt_out[11*i+7];

            // h_ridged_thin_ice
            M_h_ridged_thin_ice[i] = interp_elt_out[11*i+8];

            // h_ridged_thick_ice
            M_h_ridged_thick_ice[i] = interp_elt_out[11*i+9];

            // random_number
            M_random_number[i] = interp_elt_in[11*i+10];
        }

        std::cout<<"ELEMENT: Interp done\n";
        std::cout<<"ELEMENT: TIMER INTERPOLATION= " << chrono.elapsed() <<"s\n";
    }

    M_matrix->init(2*M_num_nodes,2*M_num_nodes,22);
    M_vector->resize(2*M_num_nodes);
    M_solution->resize(2*M_num_nodes);

    M_wind.resize(2*M_num_nodes,0.);
    M_ocean.resize(2*M_num_nodes,0.);
    M_thermo.resize(2*M_num_nodes,0.);

    //M_vair.resize(2*M_num_nodes,0.);

    M_element_depth.resize(M_num_elements,0.);
    M_h_thin.resize(M_num_elements,0.);
    M_hs_thin.resize(M_num_elements,0.);

    M_norm_Voce_ice.resize(M_num_elements);
    M_norm_Vair_ice.resize(M_num_elements);
    M_norm_Vice.resize(M_num_elements);
    M_element_ssh.resize(M_num_elements,0.);
    M_ssh.resize(M_num_nodes,0.);
    M_vssh.resize(M_num_nodes);

    M_Vair_factor.resize(M_num_elements);
    M_Voce_factor.resize(M_num_elements);

    M_basal_factor.resize(M_num_elements);
    M_fcor.resize(M_num_elements);
    M_Vcor.resize(M_VT.size());

    M_Cohesion.resize(M_num_elements);
    M_Compressive_strength.resize(M_num_elements);

}

void
FiniteElement::assemble()
{
    double coef_V, coef_Voce, coef_Vair, coef_basal, coef_X, coef_Y, coef_C;
    double coef = 0;
    double coef_P = 0.;
    double mass_e = 0;
    double surface_e = 0;
    double g_ssh_e_x = 0.;
    double g_ssh_e = 0.;
    double g_ssh_e_y = 0.;

    std::vector<double> B0T(18,0);
    std::vector<double> B0Tj_Dunit(6,0);
    std::vector<double> B0Tj_Dunit_B0Ti(4,0);
    std::vector<double> B0T_Dunit_B0T(36,0);
    std::vector<double> B0Tj_Dunit_comp(6,0);
    std::vector<double> B0Tj_Dunit_comp_B0Ti(4,0);
    std::vector<double> B0T_Dunit_comp_B0T(36,0);

    std::vector<double> sigma_P(3,0); /* temporary variable for the resistance to the compression */
    std::vector<double> B0Tj_sigma_h(2,0);

    double B0Tj_Dunit_tmp0, B0Tj_Dunit_tmp1;
    double B0Tj_Dunit_B0Ti_tmp0, B0Tj_Dunit_B0Ti_tmp1, B0Tj_Dunit_B0Ti_tmp2, B0Tj_Dunit_B0Ti_tmp3;
    double B0Tj_Dunit_comp_tmp0, B0Tj_Dunit_comp_tmp1;
    double B0Tj_Dunit_comp_B0Ti_tmp0, B0Tj_Dunit_comp_B0Ti_tmp1, B0Tj_Dunit_comp_B0Ti_tmp2, B0Tj_Dunit_comp_B0Ti_tmp3;
    double mloc = 0;

    // std::vector<double> B0Tj_Dunit_B0Ti(4,0);
    // std::vector<double> B0Tj_Dunit_comp_B0Ti(4,0);
    // std::vector<double> B0T_Dunit_B0T(36,0);
    // std::vector<double> B0T_Dunit_comp_B0T(36,0);

    std::cout<<"Assembling starts\n";
    chrono.restart();
    std::vector<double> data(36);
    std::vector<double> fvdata(6);
    std::vector<int> rcindices(6);

    double duu, dvu, duv, dvv, fuu, fvv;
    int index_u, index_v;


    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        std::vector<double> shapecoeff = this->shapeCoeff(*it,M_mesh);

        for (int i=0; i<18; ++i)
        {
            if (i < 3)
            {
                B0T[2*i] = shapecoeff[i];
                B0T[12+2*i] = shapecoeff[i+3];
                B0T[13+2*i] = shapecoeff[i];
            }
            else if (i < 6)
            {
                B0T[2*i+1] = shapecoeff[i];
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

        coef = young*(1-M_damage[cpt])*M_thick[cpt]*std::exp(ridging_exponent*(1-M_conc[cpt]));

        coef_P = 0.;
        if(M_divergence_rate[cpt] < 0.)
        {
            coef_P = compression_factor*std::pow(M_thick[cpt],exponent_compression_factor)*std::exp(ridging_exponent*(1-M_conc[cpt]));
            coef_P = coef_P/(std::abs(M_divergence_rate[cpt])+divergence_min);
            //std::cout<<"Coeff= "<< coef_P <<"\n";
        }

        /* Compute the value that only depends on the element */
        mass_e = rhoi*M_thick[cpt] + rhos*M_snow_thick[cpt];
        mass_e = (M_conc[cpt] > 0.) ? (mass_e/M_conc[cpt]):0.;
        surface_e = this->measure(*it,M_mesh);

        // /* compute the x and y derivative of g*ssh */
        g_ssh_e_x = 0.;
        g_ssh_e_y = 0.;
        for(int i=0; i<3; i++)
        {
            g_ssh_e = (vm["simul_in.gravity"].as<double>())*M_ssh[it->indices[i]-1] /*g_ssh*/;   /* g*ssh at the node k of the element e */
            g_ssh_e_x += shapecoeff[i]*g_ssh_e; /* x derivative of g*ssh */
            g_ssh_e_y += shapecoeff[i+3]*g_ssh_e; /* y derivative of g*ssh */
        }

        coef_C     = mass_e*M_fcor[cpt];              /* for the Coriolis term */
        coef_V     = mass_e/time_step;             /* for the inertial term */
        coef_X     = - mass_e*g_ssh_e_x;              /* for the ocean slope */
        coef_Y     = - mass_e*g_ssh_e_y;              /* for the ocean slope */
        coef_Vair  = M_Vair_factor[cpt];             /* for the wind stress */
        coef_Voce  = M_Voce_factor[cpt];             /* for the ocean stress */
        coef_basal = M_basal_factor[cpt];            /* for the basal stress */

        // std::cout<<"************************\n";
        // std::cout<<"Coef_X= "<< coef_X <<"\n";
        // std::cout<<"Coef_Y= "<< coef_Y <<"\n";

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
                    B0Tj_Dunit_tmp0 += B0T[kk*6+2*j]*M_Dunit[3*i+kk];
                    B0Tj_Dunit_tmp1 += B0T[kk*6+2*j+1]*M_Dunit[3*i+kk];

                    B0Tj_Dunit_comp_tmp0 += B0T[kk*6+2*j]*M_Dunit_comp[3*i+kk];
                    B0Tj_Dunit_comp_tmp1 += B0T[kk*6+2*j+1]*M_Dunit_comp[3*i+kk];
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

        // std::vector<double> data(36);
        // std::vector<double> fvdata(6);
        // std::vector<int> rcindices(6);

        for(int j=0; j<3; j++)
        {
            /* Column corresponding to indice j (we also assemble terms in col+1) */
            //col = (mwIndex)it[2*j]-1; /* -1 to use the indice convention of C */
            fuu=0.;
            fvv=0.;

            int index_u = it->indices[j]-1;
            int index_v = it->indices[j]-1+M_num_nodes;

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
                mloc = M_Mass[3*j+i];

                B0Tj_sigma_h[0] = 0.;
                B0Tj_sigma_h[1] = 0.;

                for(int k=0; k<3; k++)
                {
                    B0Tj_sigma_h[0] += B0T[k*6+2*i]*(M_sigma[3*cpt+k]*M_thick[k]+sigma_P[k]);
                    B0Tj_sigma_h[1] += B0T[k*6+2*i+1]*(M_sigma[3*cpt+k]*M_thick[k]+sigma_P[k]);
                }

                /* ---------- UU component */
                duu = surface_e*( mloc*(coef_Vair+coef_Voce*std::cos(ocean_turning_angle_rad)+coef_V+coef_basal)+B0Tj_Dunit_B0Ti[0]*coef*time_step+B0Tj_Dunit_comp_B0Ti[0]*coef_P);

                /* ---------- VU component */
                dvu = surface_e*(+B0Tj_Dunit_B0Ti[1]*coef*time_step+B0Tj_Dunit_comp_B0Ti[1]*coef_P);

                /* ---------- UV component */
                duv = surface_e*(+B0Tj_Dunit_B0Ti[2]*coef*time_step+B0Tj_Dunit_comp_B0Ti[2]*coef_P);

                /* ---------- VV component */
                dvv = surface_e*( mloc*(coef_Vair+coef_Voce*std::cos(ocean_turning_angle_rad)+coef_V+coef_basal)+B0Tj_Dunit_B0Ti[3]*coef*time_step+B0Tj_Dunit_comp_B0Ti[3]*coef_P);

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

#if 1
                fuu += surface_e*( mloc*( coef_Vair*M_wind[index_u]+coef_Voce*std::cos(ocean_turning_angle_rad)*M_ocean[index_u]+coef_X+coef_V*M_VT[index_u]) - B0Tj_sigma_h[0]/3);
                fuu += surface_e*( mloc*( -coef_Voce*std::sin(ocean_turning_angle_rad)*(M_ocean[index_u]-M_VT[index_u])-coef_C*M_Vcor[index_u]) );

                fvv += surface_e*( mloc*( coef_Vair*M_wind[index_v]+coef_Voce*std::cos(ocean_turning_angle_rad)*M_ocean[index_v]+coef_Y+coef_V*M_VT[index_v]) - B0Tj_sigma_h[1]/3);
                fvv += surface_e*( mloc*( -coef_Voce*std::sin(ocean_turning_angle_rad)*(M_ocean[index_v]-M_VT[index_v])-coef_C*M_Vcor[index_v]) );
#endif
            }

            if (cpt < 0)
            {
                std::cout<<"**********************************\n";
                std::cout<<"mloc        = "<< mloc <<"\n";
                std::cout<<"coef_Vair   = "<< coef_Vair <<"\n";
                std::cout<<"coef_Voce   = "<< coef_Voce <<"\n";
                std::cout<<"coef_X      = "<< coef_X <<"\n";
                std::cout<<"coef_Y      = "<< coef_Y <<"\n";
                std::cout<<"coef_V      = "<< coef_V <<"\n";
                std::cout<<"wind U      = "<< M_wind[it->indices[j]-1] <<"\n";
                std::cout<<"wind V      = "<< M_wind[it->indices[j]-1+M_num_nodes] <<"\n";
            }

            // std::cout<<"**********************************\n";
            // std::cout<<"fuu= "<< fuu <<"\n";
            // std::cout<<"fvv= "<< fvv <<"\n";

            fvdata[2*j] = fuu;
            fvdata[2*j+1] = fvv;

            rcindices[2*j] = index_u;
            rcindices[2*j+1] = index_v;
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

        // if (cpt == 0)
        //     for(int k=0; k<6; k++)
        //         std::cout<<"INDEX["<< rcindices[k] <<"]\n";

        M_matrix->addMatrix(&rcindices[0], rcindices.size(),
                            &rcindices[0], rcindices.size(), &data[0]);

        M_vector->addVector(&rcindices[0], rcindices.size(), &fvdata[0]);

        ++cpt;
    }

    M_matrix->close();
    std::cout<<"Assembling done\n";
    std::cout<<"TIMER ASSEMBLY= " << chrono.elapsed() <<"s\n";

    chrono.restart();
    M_matrix->on(M_dirichlet_nodes,*M_vector);
    std::cout<<"TIMER DBCA= " << chrono.elapsed() <<"s\n";

    std::cout<<"[PETSC MATRIX] CLOSED      = "<< M_matrix->closed() <<"\n";
    std::cout<<"[PETSC MATRIX] SIZE        = "<< M_matrix->size1() << " " << M_matrix->size2() <<"\n";
    std::cout<<"[PETSC MATRIX] SYMMETRIC   = "<< M_matrix->isSymmetric() <<"\n";
    std::cout<<"[PETSC MATRIX] NORM        = "<< M_matrix->linftyNorm() <<"\n";

    //M_matrix->printMatlab("stiffness.m");
    //M_vector->printMatlab("rhs.m");
}

void
FiniteElement::tensors()
{
    M_Dunit.resize(9,0);
    M_Dunit_comp.resize(9,0);
    M_Mass.resize(9,0);

    for (int k=0; k<6; k+=3)
    {
        for (int kk=0; kk<2; ++kk )
        {
            M_Dunit[k+kk] = (1-((k+kk)%2)*(1-nu0))/(1-std::pow(nu0,2.));
            M_Dunit_comp[k+kk] = 1.;
        }
    }
    M_Dunit[8] = (1-nu0)/(2.*(1-std::pow(nu0,2.)));

    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            M_Mass[3*i+j] = ((i == j) ? 2.0 : 1.0)/12.0;
            //std::cout<< std::left << std::setw(12) << Mass[3*i+j] <<"  ";
        }

        //std::cout<<"\n";
    }

#if 0

    M_B0T.resize(M_num_elements);
    M_B0T_Dunit_B0T.resize(M_num_elements);
    M_B0T_Dunit_comp_B0T.resize(M_num_elements);

    std::vector<double> B0T(18,0);
    std::vector<double> B0Tj_Dunit(6,0);
    std::vector<double> B0T_Dunit_B0T(36,0);
    std::vector<double> B0Tj_Dunit_comp(6,0);
    std::vector<double> B0T_Dunit_comp_B0T(36,0);

    double B0Tj_Dunit_tmp0, B0Tj_Dunit_tmp1;
    double B0Tj_Dunit_B0Ti_tmp0, B0Tj_Dunit_B0Ti_tmp1, B0Tj_Dunit_B0Ti_tmp2, B0Tj_Dunit_B0Ti_tmp3;
    double B0Tj_Dunit_comp_tmp0, B0Tj_Dunit_comp_tmp1;
    double B0Tj_Dunit_comp_B0Ti_tmp0, B0Tj_Dunit_comp_B0Ti_tmp1, B0Tj_Dunit_comp_B0Ti_tmp2, B0Tj_Dunit_comp_B0Ti_tmp3;

    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        std::vector<double> shapecoeff = this->shapeCoeff(*it,M_mesh);

        for (int i=0; i<18; ++i)
        {
            if (i < 3)
            {
                B0T[2*i] = shapecoeff[i];
                B0T[12+2*i] = shapecoeff[i+3];
                B0T[13+2*i] = shapecoeff[i];
            }
            else if (i < 6)
            {
                B0T[2*i+1] = shapecoeff[i];
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
                    B0Tj_Dunit_tmp0 += B0T[kk*6+2*j]*M_Dunit[3*i+kk];
                    B0Tj_Dunit_tmp1 += B0T[kk*6+2*j+1]*M_Dunit[3*i+kk];

                    B0Tj_Dunit_comp_tmp0 += B0T[kk*6+2*j]*M_Dunit_comp[3*i+kk];
                    B0Tj_Dunit_comp_tmp1 += B0T[kk*6+2*j+1]*M_Dunit_comp[3*i+kk];
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

        M_B0T[cpt]                = B0T;
        M_B0T_Dunit_B0T[cpt]      = B0T_Dunit_B0T;
        M_B0T_Dunit_comp_B0T[cpt] = B0T_Dunit_comp_B0T;

        ++cpt;
    }

#endif

}

void
FiniteElement::cohesion()
{
    for (int i=0; i<M_Cohesion.size(); ++i)
        M_Cohesion[i] = C_fix+C_alea*(M_random_number[i]-0.5);

    for (int i=0; i<M_Compressive_strength.size(); ++i)
        M_Compressive_strength[i] = compr_strength*scale_coef;
}

void
FiniteElement::update()
{
    std::vector<double> thickness_new = M_thick;
    std::vector<double> snow_thickness_new = M_snow_thick;
    std::vector<double> concentration_new = M_conc;
    std::vector<double> damage_new = M_damage;
    std::vector<double> thin_thickness_new = M_h_thin;
    std::vector<double> thin_snow_thickness_new = M_hs_thin;
    std::vector<double> h_ridged_thin_ice_new = M_h_ridged_thin_ice;
    std::vector<double> h_ridged_thick_ice_new = M_h_ridged_thick_ice;

    /* deformation, deformation rate and internal stress tensor and temporary variables */
    double epsilon_veloc_i;
    std::vector<double> epsilon_veloc(3);
    std::vector<double> sigma_pred(3);
    double sigma_dot_i;

    /* some variables used for the advection*/
    double surface, surface_new, ar, ar_new;
    double ice_surface, ice_volume, snow_volume;
    double thin_ice_surface, thin_ice_volume, thin_snow_volume;
    double ridging_thin_ice, ridging_thick_ice, ridging_snow_thin_ice;
    double ridged_thin_ice_volume, ridged_thick_ice_volume;

    /* invariant of the internal stress tensor and some variables used for the damaging process*/
    double sigma_s, sigma_n;
    double tract_max;
    double tmp, sigma_target;

    /* some variables used for the ice redistribution*/
    double tanalpha, rtanalpha, del_v, del_c, del_vs, new_v_thin;

    /* set constants for the ice redistribution */
	tanalpha  = h_thin_max/c_thin_max;
	rtanalpha = 1/tanalpha;

    std::vector<double> UM_P = M_UM;

    //std::cout<<"CURRENT: DIRICHLET NODES= "<< M_dirichlet_flags.size() <<"\n";
    //std::cout<<"CURRENT: NEUMANN   NODES= "<< M_neumann_flags.size() <<"\n";

    // std::cout<<"UM MIN B= "<< *std::min_element(M_UM.begin(),M_UM.end()) <<"\n";
    // std::cout<<"UM MAX B= "<< *std::max_element(M_UM.begin(),M_UM.end()) <<"\n";

    for (int nd=0; nd<M_UM.size(); ++nd)
    {
        M_UM[nd] = M_UM[nd] + time_step*M_VT[nd];
    }

    for (const int& nd : M_neumann_nodes)
    {
        M_UM[nd] = UM_P[nd];
    }

    // std::cout<<"UM MIN A= "<< *std::min_element(M_UM.begin(),M_UM.end()) <<"\n";
    // std::cout<<"UM MAX A= "<< *std::max_element(M_UM.begin(),M_UM.end()) <<"\n";

    std::cout<<"VT MIN= "<< *std::min_element(M_VT.begin(),M_VT.end()) <<"\n";
    std::cout<<"VT MAX= "<< *std::max_element(M_VT.begin(),M_VT.end()) <<"\n";

    // std::cout<<"Jacobian init   = "<< this->measure(*M_elements.begin(),M_mesh, UM_P) <<"\n";
    // std::cout<<"Jacobian current= "<< this->measure(*M_elements.begin(),M_mesh,M_UM) <<"\n";

    std::vector<double> B0T(18,0);
    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        std::vector<double> shapecoeff = this->shapeCoeff(*it,M_mesh);

        for (int i=0; i<18; ++i)
        {
            if (i < 3)
            {
                B0T[2*i] = shapecoeff[i];
                B0T[12+2*i] = shapecoeff[i+3];
                B0T[13+2*i] = shapecoeff[i];
            }
            else if (i < 6)
            {
                B0T[2*i+1] = shapecoeff[i];
            }
        }

        /*======================================================================
         * Diagnostic:
         * Elastic deformation and instantaneous deformation rate
         *======================================================================
         */

        /* Compute the elastic deformation and the instantaneous deformation rate */
        for(int i=0;i<3;i++)
        {
            epsilon_veloc_i = 0.0;
            for(int j=0;j<3;j++)
            {
                /* deformation */
                //col = (mwIndex)it[j]-1;
                epsilon_veloc_i += B0T[i*6 + 2*j]*M_VT[it->indices[j]-1]  ;
                epsilon_veloc_i += B0T[i*6 + 2*j + 1]*M_VT[it->indices[j]-1+M_num_nodes]  ;
            }

            epsilon_veloc[i] = epsilon_veloc_i;
        }

        M_divergence_rate[cpt]= (epsilon_veloc[0]+epsilon_veloc[1]);

        //std::cout<<"Rate= "<< M_divergence_rate[cpt] <<"\n";

        /*======================================================================
         * Update the internal stress
         *======================================================================
         */

        for(int i=0;i<3;i++)
        {
            sigma_dot_i = 0.0;
            for(int j=0;j<3;j++)
            {
                sigma_dot_i += std::exp(ridging_exponent*(1-M_conc[cpt]))*young*(1.-M_damage[cpt])*M_Dunit[i*3 + j]*epsilon_veloc[j];
            }

            //M_sigma[3*cpt+i] = M_sigma[3*cpt+i] + time_step*sigma_dot_i;
            M_sigma[3*cpt+i] += time_step*sigma_dot_i;
            sigma_pred[i]    = M_sigma[3*cpt+i] + time_step*sigma_dot_i;
        }

        /*======================================================================
         * Correct the internal stress and the damage
         *======================================================================
         */

        /* Compute the shear and normal stress, which are two invariants of the internal stress tensor */

        sigma_s=std::sqrt(std::pow((sigma_pred[0]-sigma_pred[1])/2.,2.)+std::pow(sigma_pred[2],2.));
        sigma_n=         (sigma_pred[0]+sigma_pred[1])/2.;

        //std::cout<<"sigma_n= "<< sigma_n <<"\n";

        /* minimum and maximum normal stress */
        tract_max=tract_coef*M_Cohesion[cpt]/tan_phi;

        /* Correction of the damage */

        if((sigma_n>tract_max) || (sigma_n<(-M_Compressive_strength[cpt])))
        {
            if(sigma_n>tract_max)
            {
                sigma_target=tract_max;
            }
            else
            {
                sigma_target=-M_Compressive_strength[cpt];
            }

            tmp=1.0-sigma_target/sigma_n*(1-M_damage[cpt]);

            if(tmp>damage_new[cpt])
            {
                damage_new[cpt]=tmp;
            }
        }

        if(sigma_s>M_Cohesion[cpt]-sigma_n*tan_phi)
        {
            tmp=1.0-M_Cohesion[cpt]/(sigma_s+sigma_n*tan_phi)*(1-M_damage[cpt]);

            if(tmp>damage_new[cpt])
            {
                damage_new[cpt]=tmp;
            }
        }

        /*
         * Diagnostic:
         * Recompute the internal stress
         */
        for(int i=0;i<3;i++)
        {
            if(M_damage[cpt]<1)
            {
                M_sigma[3*cpt+i] = (1.-damage_new[cpt])/(1.-M_damage[cpt])*M_sigma[3*cpt+i] ;
            }
            else
            {
                M_sigma[3*cpt+i] = 0. ;
            }
        }


        /*======================================================================
         * Update:
         * Ice damage
         * We use now a constant healing rate defined as 1/time_recovery_damage
         * so that we are now able to reset the damage to 0.
         * otherwise, it will never heal completely.
         * time_recovery_damage still depends on the temperature when themodynamics is activated.
         *======================================================================
         */
        tmp=1./(1.-damage_new[cpt]);
        tmp=tmp-1000*time_step/time_relaxation_damage;
        tmp=((tmp>1.)?(tmp):(1.));
        M_damage[cpt]=-1./tmp + 1.;

        /*======================================================================
         * Update:
         * Ice and snow thickness, and concentration using a Lagrangian or an Eulerian scheme
         *======================================================================
         */

        surface = this->measure(*it,M_mesh, UM_P);
        surface_new = this->measure(*it,M_mesh,M_UM);

        // std::cout<<"SURFACE    = "<< std::setprecision(18) << surface <<"\n";
        // std::cout<<"SURFACE_NEW= "<< std::setprecision(18) << surface_new <<"\n";

        ice_surface = M_conc[cpt]*surface;
        ice_volume = M_thick[cpt]*surface;
        snow_volume = M_snow_thick[cpt]*surface;
        thin_ice_volume = M_h_thin[cpt]*surface;
        thin_snow_volume = M_hs_thin[cpt]*surface;
        ridged_thin_ice_volume = M_h_ridged_thin_ice[cpt]*surface;
        ridged_thick_ice_volume = M_h_ridged_thick_ice[cpt]*surface;

        if(M_conc[cpt]>0.)
        {
            /* updated values */
            concentration_new[cpt]    = ice_surface/surface_new;
            thickness_new[cpt]        = ice_volume/surface_new;
            snow_thickness_new[cpt]   = snow_volume/surface_new;

            thin_thickness_new[cpt]        = thin_ice_volume/surface_new;
            thin_snow_thickness_new[cpt]   = thin_snow_volume/surface_new;

            h_ridged_thin_ice_new[cpt]    =   ridged_thin_ice_volume/surface_new;
            h_ridged_thick_ice_new[cpt]   =   ridged_thick_ice_volume/surface_new;

            /* Ridging scheme */
            if(surface_new<surface)
            {
                ridging_thin_ice=(surface-surface_new)/surface_new*M_h_thin[cpt];
                ridging_snow_thin_ice=(surface-surface_new)/surface_new*M_hs_thin[cpt];

                ridging_thin_ice=((ridging_thin_ice<M_h_thin[cpt])?(ridging_thin_ice):(M_h_thin[cpt]));
                ridging_snow_thin_ice=((ridging_snow_thin_ice<M_hs_thin[cpt])?(ridging_snow_thin_ice):(M_hs_thin[cpt])) ;

                thickness_new[cpt] += ridging_thin_ice ;
                thin_thickness_new[cpt] -= ridging_thin_ice ;

                snow_thickness_new[cpt] += ridging_snow_thin_ice ;
                thin_snow_thickness_new[cpt] -= ridging_snow_thin_ice ;

                h_ridged_thin_ice_new[cpt] += ridging_thin_ice;
                concentration_new[cpt] += ridging_thin_ice/ridge_h;

                /* upper bounds (only for the concentration) */
                ridging_thin_ice = ((concentration_new[cpt]<1.)?(0.):(thin_thickness_new[cpt])) ;
                ridging_snow_thin_ice = ((concentration_new[cpt]<1.)?(0.):(thin_snow_thickness_new[cpt])) ;

                ridging_thick_ice=((concentration_new[cpt]<1.)?(0.):(thickness_new[cpt]*(concentration_new[cpt]-1.)));
                concentration_new[cpt] = ((concentration_new[cpt]<1.)?(concentration_new[cpt]):(1.)) ;

                snow_thickness_new[cpt] += ridging_snow_thin_ice ;
                thin_snow_thickness_new[cpt] -= ridging_snow_thin_ice ;

                thickness_new[cpt] += ridging_thin_ice;
                thin_thickness_new[cpt] -= ridging_thin_ice;
            }

            /* lower bounds */
            concentration_new[cpt] = ((concentration_new[cpt]>0.)?(concentration_new[cpt] ):(0.)) ;
            thickness_new[cpt]        = ((thickness_new[cpt]>0.)?(thickness_new[cpt]     ):(0.)) ;
            snow_thickness_new[cpt]   = ((snow_thickness_new[cpt]>0.)?(snow_thickness_new[cpt]):(0.)) ;

            thin_thickness_new[cpt]        = ((thin_thickness_new[cpt]>0.)?(thin_thickness_new[cpt] ):(0.)) ;
            thin_snow_thickness_new[cpt]   = ((thin_snow_thickness_new[cpt]>0.)?(thin_snow_thickness_new[cpt]):(0.)) ;
        }
        else
        {
            concentration_new[cpt] = M_conc[cpt];
            thickness_new[cpt]     = M_thick[cpt];
            snow_thickness_new[cpt]= M_snow_thick[cpt];

            thin_thickness_new[cpt] = M_h_thin[cpt];
            thin_snow_thickness_new[cpt] = M_hs_thin[cpt];
        }

        /* For the Lagrangian scheme, we do not update the variables for the elements having one node on the open boundary. */

        if(std::find(M_neumann_flags.begin(),M_neumann_flags.end(),it->indices[0]-1) != M_neumann_flags.end() ||
           std::find(M_neumann_flags.begin(),M_neumann_flags.end(),it->indices[1]-1) != M_neumann_flags.end() ||
           std::find(M_neumann_flags.begin(),M_neumann_flags.end(),it->indices[2]-1) != M_neumann_flags.end())
        {
            //std::cout<<"--------------------------------------FOUND\n";

            concentration_new[cpt]    = M_conc[cpt] ;
            thickness_new[cpt]        = M_thick[cpt] ;
            snow_thickness_new[cpt]   = M_snow_thick[cpt] ;

            thin_thickness_new[cpt]        = M_h_thin[cpt] ;
            thin_snow_thickness_new[cpt]   = M_hs_thin[cpt] ;

            h_ridged_thin_ice_new[cpt] = M_h_ridged_thin_ice[cpt];
            h_ridged_thick_ice_new[cpt] = M_h_ridged_thick_ice[cpt];
        }

        /* Compute the redistribution of thin ice. */
        /* Returns the change in volume and concentration of thick ice as well as the
         * change in volume of thin ice. It is called after the
         * dynamics are done. */

        if(thin_thickness_new[cpt]>0.)
        {
            thin_ice_redistribute(thin_thickness_new[cpt], thin_snow_thickness_new[cpt], 0., concentration_new[cpt],
                                  tanalpha, rtanalpha, h_thin_max, &new_v_thin, &del_v, &del_c, &del_vs);

            concentration_new[cpt]       += del_c;

            thickness_new[cpt]           += del_v;
            thin_thickness_new[cpt]      -= del_v;

            snow_thickness_new[cpt]      += del_vs;
            thin_snow_thickness_new[cpt] -= del_vs;
        }
        else
        {
            snow_thickness_new[cpt] += thin_snow_thickness_new[cpt];
            thin_snow_thickness_new[cpt] = 0. ;
        }

        ++cpt;
    }

    // std::cout<<"****************************\n";
    // std::cout<<"SIGMA MIN= "<< *std::min_element(M_sigma.begin(),M_sigma.end()) <<"\n";
    // std::cout<<"SIGMA MAX= "<< *std::max_element(M_sigma.begin(),M_sigma.end()) <<"\n";

    M_thick = thickness_new;
    M_snow_thick = snow_thickness_new;
    M_conc = concentration_new;
    M_damage = damage_new;
    M_h_thin = thin_thickness_new;
    M_hs_thin = thin_snow_thickness_new;
    M_h_ridged_thin_ice = h_ridged_thin_ice_new;
    M_h_ridged_thick_ice = h_ridged_thick_ice_new;
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
    //M_solution->printMatlab("solution.m");

    Environment::logMemoryUsage("");
}

void
FiniteElement::run()
{
    this->init();

    int pcpt = 0;
    int niter = 0;
    current_time = time_init /*+ pcpt*time_step/(24*3600.0)*/;

    std::cout<<"TIMESTEP= "<< time_step <<"\n";
    std::cout<<"DURATION= "<< duration <<"\n";

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

    std::cout<<"INIT TIME= "<< time_init <<"\n";

    // main loop for nextsim program

    double displacement_factor = 1.;
    double minang = 0.;
    bool is_running = true;

    minang = this->minAngle(M_mesh);
    if (minang < vm["simul_in.regrid_angle"].as<double>())
    {
        std::cout<<"invalid regridding angle: should be smaller than the minimal angle in the intial grid\n";
        throw std::logic_error("invalid regridding angle: should be smaller than the minimal angle in the intial grid");
    }

    while (is_running)
    {
        is_running = ((pcpt+1)*time_step) < duration;

        if (pcpt > 100)
            is_running = false;

        //if(pcpt >0)
        //    this->exportResults(pcpt+3000);

        current_time = time_init + pcpt*time_step/(24*3600.0);
        //std::cout<<"TIME STEP "<< pcpt << " for "<< current_time <<"\n";
        std::cout<<"TIME STEP "<< pcpt << " for "<< current_time << " + "<< pcpt*time_step/(24*3600.0) <<"\n";

        // step 0: preparation
        // remeshing and remapping of the prognostic variables

        bool regrid_done = false;

        if (vm["simul_in.regrid"].as<std::string>() == "bamg")
        {
            minang = this->minAngle(M_mesh,M_UM,displacement_factor);
            std::cout<<"REGRID ANGLE= "<< minang <<"\n";

            if ((minang < vm["simul_in.regrid_angle"].as<double>()) || (pcpt ==0) )
            {
                regrid_done = true;
                std::cout<<"Regriding starts\n";
                this->regrid(pcpt);
                std::cout<<"Regriding done\n";
                //if(pcpt >0)
                //    this->exportResults(pcpt+1000);
            }
        }

        if (pcpt == 0)
            this->initSimulation();

        if ((pcpt==0) || (regrid_done))
        {
            //this->forcingWind(vm["simul_in.constant_u"].as<double>(),vm["simul_in.constant_v"].as<double>());
            //this->forcingOcean(0.,0.);
            this->forcingThermo(0.,0.);
            this->bathymetry();
            this->tensors();
            this->cohesion();
        }

        this->timeInterpolation(pcpt);
        this->forcingWind(regrid_done);
        this->forcingOcean(regrid_done);
        //this->timeInterpolation(pcpt);
        this->computeFactors(pcpt);

        this->assemble();
        this->solve();
        this->updateVelocity();

        //if(pcpt >0)
        //    this->exportResults(pcpt+2000);

        this->update();

        this->exportResults(pcpt+1);
        //this->asrWind();
        //this->loadTopazOcean();
        ++pcpt;
    }

    this->exportResults(1);
}

void
FiniteElement::updateVelocity()
{
    M_VTMM = M_VTM;
    M_VTM = M_VT;
    M_VT = M_solution->container();

    //std::vector<double> speed_c_scaling_test(bamgmesh->NodalElementConnectivitySize[0]);
    int elt_num;
    double c_max_nodal_neighbour;
    double speed_c_scaling;

    for (int i=0; i<bamgmesh->NodalElementConnectivitySize[0]; ++i)
    {
        std::vector<double> cloc_elts(bamgmesh->NodalElementConnectivitySize[1]);
        int j = 0;
        for (j=0; j<bamgmesh->NodalElementConnectivitySize[1]; ++j)
        {
            elt_num = bamgmesh->NodalElementConnectivity[bamgmesh->NodalElementConnectivitySize[1]*i+j]-1;

            if ((0 <= elt_num) && (elt_num < M_mesh.numTriangles()) && (elt_num != NAN))
            {
                cloc_elts[j] = M_conc[elt_num];
            }
            else
            {
                break;
            }
        }

        cloc_elts.resize(j);
        c_max_nodal_neighbour = *std::max_element(cloc_elts.begin(),cloc_elts.end());
        c_max_nodal_neighbour /= vm["simul_in.drift_limit_concentration"].as<double>();
        speed_c_scaling = std::min(1.,c_max_nodal_neighbour);
        //std::cout<<"c_max_nodal_neighbour["<< i <<"]= "<< c_max_nodal_neighbour <<"\n";
        //std::cout<<"speed_c_scaling["<< i <<"]= "<< speed_c_scaling <<"\n";
        //speed_c_scaling_test[i] = speed_c_scaling;

        // linear scaling of ice velocity
        M_VT[i] = speed_c_scaling*M_VT[i];
        M_VT[i+M_num_nodes] = speed_c_scaling*M_VT[i+M_num_nodes];
    }

    // std::cout<<"MAX SPEED= "<< *std::max_element(speed_c_scaling_test.begin(),speed_c_scaling_test.end()) <<"\n";
    // std::cout<<"MIN SPEED= "<< *std::min_element(speed_c_scaling_test.begin(),speed_c_scaling_test.end()) <<"\n";

    std::cout<<"MAX SPEED= "<< *std::max_element(M_conc.begin(),M_conc.end()) <<"\n";
    std::cout<<"MIN SPEED= "<< *std::min_element(M_conc.begin(),M_conc.end()) <<"\n";
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
FiniteElement::computeFactors(int pcpt)
{
    double welt_oce_ice = 0.;
    double welt_air_ice = 0.;
    double welt_ice = 0.;
    double welt_ssh = 0.;
    int nind;
    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        welt_oce_ice = 0.;
        welt_air_ice = 0.;
        welt_ice = 0.;
        welt_ssh = 0.;

        for (int i=0; i<3; ++i)
        {
            nind = it->indices[i]-1;
            welt_oce_ice += std::sqrt(std::pow(M_VT[nind]-M_ocean[nind],2.)+std::pow(M_VT[nind+M_num_nodes]-M_ocean[nind+M_num_nodes],2.));
            welt_air_ice += std::sqrt(std::pow(M_VT[nind]-M_wind [nind],2.)+std::pow(M_VT[nind+M_num_nodes]-M_wind [nind+M_num_nodes],2.));
            welt_ice += std::sqrt(std::pow(M_VT[nind],2.)+std::pow(M_VT[nind+M_num_nodes],2.));

            welt_ssh += M_vssh[nind];
        }

        M_norm_Voce_ice[cpt] = welt_oce_ice/3.;
        M_norm_Vair_ice[cpt] = welt_air_ice/3.;
        M_norm_Vice[cpt] = welt_ice/3.;

        M_element_ssh[cpt] = welt_ssh/3.;

        //std::cout <<"Coeff= "<< M_norm_Vice[cpt] <<"\n";
        //std::cout <<"Coeff= "<< M_element_ssh[cpt] <<"\n";

        ++cpt;
    }

    for (int i=0; i<M_Vair_factor.size(); ++i)
    {
        M_Vair_factor[i] = (vm["simul_in.lin_drag_coef_air"].as<double>()+(quad_drag_coef_air*M_norm_Vair_ice[i]));
        M_Vair_factor[i] *= Vair_coef*(vm["simul_in.rho_air"].as<double>());
        //std::cout <<"Coeff= "<< M_Vair_factor[i] <<"\n";
    }

    for (int i=0; i<M_Voce_factor.size(); ++i)
    {
        M_Voce_factor[i] = (vm["simul_in.lin_drag_coef_water"].as<double>()+(quad_drag_coef_water*M_norm_Voce_ice[i]));
        M_Voce_factor[i] *= Voce_coef*(vm["simul_in.rho_water"].as<double>());
        //std::cout <<"Coeff= "<< M_Voce_factor[i] <<"\n";
    }

    if (vm["simul_in.Lemieux_basal_k2"].as<double>() > 0 )
    {

        // for (int k = 0; k < M_num_elements; k++ )
        // {
        //     std::cout<<"DEPTH["<< k <<"]= "<< M_element_depth[k] <<"\n";
        //     std::cout<<"CONCE["<< k <<"]= "<< M_conc[k] <<"\n";
        // }


        //critical_h = std::inner_product(M_conc.begin(), M_conc.end(), M_element_depth.begin(), 0.);
        //critical_h /= (vm["simul_in.Lemieux_basal_k1"].as<double>());

        for (int i=0; i<M_basal_factor.size(); ++i)
        {
            critical_h = M_conc[i]*(M_element_depth[i]+M_element_ssh[i])/(vm["simul_in.Lemieux_basal_k1"].as<double>());
            //double _coef = ((M_thick[i]-critical_h) > 0) ? (M_thick[i]-critical_h) : 0.;
            double _coef = std::max(0., M_thick[i]-critical_h);
            M_basal_factor[i] = quad_drag_coef_air*basal_k2/(basal_drag_coef_air*(M_norm_Vice[i]+basal_u_0));
            M_basal_factor[i] *= _coef*std::exp(-basal_Cb)*(1-M_conc[i]);

            //std::cout <<"Coeff= "<< M_basal_factor[i] <<"\n";
            //std::cout <<"Coeff= "<< _coef <<"\n";
        }
    }

    // std::cout<<"critical_h   = "<< critical_h <<"\n";
    // std::cout<<"Vair_coef    = "<< Vair_coef <<"\n";
    // std::cout<<"Voce_coef    = "<< Voce_coef <<"\n";

    std::vector<double> lat = M_mesh.meanLat();
    for (int i=0; i<M_fcor.size(); ++i)
    {
        M_fcor[i] = 2*(vm["simul_in.omega"].as<double>())*std::sin(lat[i]*PI/180.);
        //std::cout <<"Coeff= "<< M_fcor[i] <<"\n";
    }

    // coriolis term
    double beta0;
    double beta1;
    double beta2;
    if (pcpt > 1)
    {
        // Adams-Bashfort 3 (AB3)
        beta0 = 23./12;
        beta1 =-16./12;
        beta2 =  5./12;
    }
    else if (pcpt == 1)
    {
        // Adams-Bashfort 2 (AB2)
        beta0 = 3/2;
        beta1 =-1/2;
        beta2 = 0  ;
    }
    else if (pcpt == 0)
    {
        // Euler explicit (Fe)
        beta0 = 1 ;
        beta1 = 0 ;
        beta2 = 0 ;
    }

    for (int i=0; i<M_Vcor.size(); ++i)
    {
        M_Vcor[i] = beta0*M_VT[i] + beta1*M_VTM[i] + beta2*M_VTMM[i];
    }

    // std::cout<<"Iter...\n";
    // std::cout<<"Max= "<< *std::max_element(M_Vcor.begin(), M_Vcor.end()) <<"\n";
    // std::cout<<"Min= "<< *std::min_element(M_Vcor.begin(), M_Vcor.end()) <<"\n";
}

void
FiniteElement::forcingWind(bool reload)//(double const& u, double const& v)
{
    switch (M_wind_type)
    {
        case forcing::WindType::CONSTANT:
            this->constantWind(vm["simul_in.constant_u"].as<double>(),vm["simul_in.constant_v"].as<double>());
            break;
        case forcing::WindType::ASR:
            this->asrWind(reload);
            break;

        default:
            std::cout << "invalid wind forcing"<<"\n";
            throw std::logic_error("invalid wind forcing");
    }
}

void
FiniteElement::constantWind(double const& u, double const& v)
{
    for (int i=0; i<M_num_nodes; ++i)
    {
        M_wind[i] = u;
        M_wind[i+M_num_nodes] = v;
    }
}

void
FiniteElement::asrWind(bool reload)
{
    if ((current_time < M_ftime_wind_range[0]) || (M_ftime_wind_range[1] < current_time) || (current_time == time_init) || reload)
    {
        if (current_time == time_init)
            std::cout<<"load forcing from ASR for initial time\n";
        else
            std::cout<<"forcing not available for the current date: load data from ASR\n";

        this->loadAsrWind();

        //std::cout<<"forcing not available for the current date\n";
        //throw std::logic_error("forcing not available for the current date");
    }

    double fdt = std::abs(M_ftime_wind_range[1]-M_ftime_wind_range[0]);
    std::vector<double> fcoeff(2);
    fcoeff[0] = std::abs(current_time-M_ftime_wind_range[1])/fdt;
    fcoeff[1] = std::abs(current_time-M_ftime_wind_range[0])/fdt;

    std::cout<<"LINEAR COEFF 1= "<< fcoeff[0] <<"\n";
    std::cout<<"LINEAR COEFF 2= "<< fcoeff[1] <<"\n";

    for (int i=0; i<2*M_num_nodes; ++i)
    {
        M_wind[i] = fcoeff[0]*M_vair[0][i] + fcoeff[1]*M_vair[1][i];

        // if (i<20)
        //     std::cout<<"data_out["<< i << "]= "<< M_wind[i] << " and "<< M_wind[i+M_num_nodes] <<"\n";
    }

}

void
FiniteElement::loadAsrWind()//(double const& u, double const& v)
{

    std::string current_timestr = to_date_string_ym(current_time);
    std::cout<<"TIMESTR= "<< current_timestr <<"\n";
    std::string asr_filename = (boost::format( "%1%/data/asr30km.comb.2d.%2%.nc" )
                                % Environment::nextsimDir().string()
                                % current_timestr ).str();

    std::cout<<"ASR FILE= "<< asr_filename <<"\n";

    double nb_timestep_day = 8;
    double asr_dt = 1./nb_timestep_day;
    double time_start = std::floor(current_time*nb_timestep_day)/nb_timestep_day;
    double time_end = std::ceil(current_time*nb_timestep_day)/nb_timestep_day;

    std::cout<<"TIME START= "<< std::setprecision(9) << time_start <<"\n";
    std::cout<<"TIME END  = "<< std::setprecision(9) << time_end <<"\n";

    // We always need at least two time steps to interpolate between
    if (time_end == time_start)
    {
        time_end = time_start + (1./nb_timestep_day);
    }

    M_ftime_wind_range.resize(0);
    for (double dt=time_start; dt<=time_end; dt+=asr_dt)
    {
        M_ftime_wind_range.push_back(dt);
    }

    for (int i=0; i<M_ftime_wind_range.size(); ++i)
    {
        std::cout<<"TIMEVEC["<< i <<"]= "<< M_ftime_wind_range[i] <<"\n";
    }

    // if ((current_time < time_start) || (time_end < current_time))
    // {
    //     std::cout<<"forcing not available for the current date\n";
    //     throw std::logic_error("forcing not available for the current date");
    // }

    int nb_forcing_step = M_ftime_wind_range.size();
    std::cout<<"NB_FORCING_STEP= "<< nb_forcing_step <<"\n";

    // read in re-analysis coordinates
    std::vector<double> XLAT(360);
    std::vector<double> XLON(360);

    std::vector<double> YLAT(360);
    std::vector<double> YLON(360);

    std::vector<size_t> index_start(2);
    std::vector<size_t> index_lat_end(2);
    std::vector<size_t> index_lon_end(2);

    std::vector<size_t> index_lat_start(2);
    std::vector<size_t> index_lon_start(2);

    index_lat_start[0] = 0;
    index_lat_start[1] = 0;

    index_lat_end[0] = 360;
    index_lat_end[1] = 1;

    index_lon_start[0] = 0;
    index_lon_start[1] = 0;

    index_lon_end[0] = 1;
    index_lon_end[1] = 360;

    std::vector<double> XTIME(248);
    std::vector<size_t> index_u10_start(3,0);
    std::vector<size_t> index_u10_end(3);

    //std::cout<<"READ NETCDF starts\n";
    netCDF::NcFile dataFile(asr_filename, netCDF::NcFile::read);
    netCDF::NcVar VXLAT = dataFile.getVar("XLAT");
    netCDF::NcVar VXLON = dataFile.getVar("XLONG");
    netCDF::NcVar VTIME = dataFile.getVar("time");
    netCDF::NcVar VU10 = dataFile.getVar("U10");
    netCDF::NcVar VV10 = dataFile.getVar("V10");
    //std::cout<<"READ NETCDF done\n";

    // VXLAT.getVar(index_start,index_lat_end,&XLAT[0]);
    // VXLON.getVar(index_start,index_lon_end,&XLON[0]);

    VXLAT.getVar(index_lon_start,index_lon_end,&XLAT[0]);
    VXLON.getVar(index_lon_start,index_lon_end,&XLON[0]);

    VXLAT.getVar(index_lat_start,index_lat_end,&YLAT[0]);
    VXLON.getVar(index_lat_start,index_lat_end,&YLON[0]);

    VTIME.getVar(&XTIME[0]);

    // VU10.getVar(index_u10_start,index_u10_end,&U10[0]);
    // VV10.getVar(index_u10_start,index_u10_end,&V10[0]);

    std::vector<double> X(360);
    std::vector<double> Y(360);

    double RE = 6378.273;
    mapx_class *map;
    std::string configfile = Environment::nextsimDir().string() + "/data/NpsASR.mpp";
    std::vector<char> str(configfile.begin(), configfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    for (int i=0; i<360; ++i)
    {
        X[i] = latLon2XY(XLAT[i], XLON[i], map, configfile)[0];
        Y[i] = latLon2XY(YLAT[i], YLON[i], map, configfile)[1];

        // if (i<10)
        // {
        //     //std::cout<<"X= "<< X[i] <<" and Y= "<< Y[i] <<"\n";
        //     std::cout<<"**********************\n";
        //     std::cout<<"XLAT= "<< XLAT[i] <<" and XLON= "<< XLON[i] <<"\n";
        //     std::cout<<"YLAT= "<< YLAT[i] <<" and YLON= "<< YLON[i] <<"\n";
        // }
    }

    // std::sort(X.begin(), X.end());
    // std::sort(Y.begin(), Y.end());

#if 1

    // rotate EB coordinates to fit the ASR coords

    double angle_stereo_mesh = -45;
    double angle_stereo_ASR = -175;
    double diff_angle = -(angle_stereo_mesh-angle_stereo_ASR)*PI/180.;

    std::cout<<"VALUE= "<< from_date_string("1901-01-01") <<"\n";
    std::for_each(XTIME.begin(), XTIME.end(), [&](double& f){ f = f/24.0+from_date_string("1901-01-01"); });
    // for (int i=0; i<248; ++i)
    // {
    //     std::cout<<"TIME["<< i <<"]= "<< XTIME[i] <<"\n";
    // }

    for (int i=0; i<M_ftime_wind_range.size(); ++i)
    {
        std::cout<<"---TIMEVEC["<< i <<"]= "<< M_ftime_wind_range[i] <<" : current_time= "<< current_time <<"\n";
    }

    std::vector<double> fvair(2*M_num_nodes);

    for (int fstep=0; fstep < nb_forcing_step; ++fstep)
    {
        double ftime = M_ftime_wind_range[fstep];

        if (to_date_string_ym(std::floor(ftime)) != to_date_string_ym(current_time))
        {
            std::string f_timestr = to_date_string_ym(std::floor(ftime));
            std::cout<<"F_TIMESTR= "<< f_timestr <<"\n";
            asr_filename = (boost::format( "%1%/data/asr30km.comb.2d.%2%.nc" )
                            % Environment::nextsimDir().string()
                            % to_date_string_ym(std::floor(ftime)) ).str();

            std::cout<<"ASR_FILENAME= "<< asr_filename <<"\n";

            netCDF::NcFile fdataFile(asr_filename, netCDF::NcFile::read);
            netCDF::NcVar FVTIME = dataFile.getVar("time");
            FVTIME.getVar(&XTIME[0]);
            std::for_each(XTIME.begin(), XTIME.end(), [&](double& f){ f = f/24.0+from_date_string("1901-01-01"); });

            // for (int i=0; i<248; ++i)
            // {
            //     std::cout<<"TIME["<< i <<"]= "<< XTIME[i] <<"\n";
            // }
        }

        auto it = std::find(XTIME.begin(), XTIME.end(), ftime);
        int index = std::distance(XTIME.begin(),it);
        std::cout<<"FIND "<< ftime <<" in index "<< index <<"\n";

        index_u10_start[0] = index;
        index_u10_start[1] = 0;
        index_u10_start[2] = 0;

        index_u10_end[0] = 1;
        index_u10_end[1] = 360;
        index_u10_end[2] = 360;

        std::vector<double> data_in_u10(360*360);
        std::vector<double> data_in_v10(360*360);

        VU10.getVar(index_u10_start,index_u10_end,&data_in_u10[0]);
        VV10.getVar(index_u10_start,index_u10_end,&data_in_v10[0]);

        // for (int i=0; i<10; ++i)
        // {
        //     for (int j=0; j<10; ++j)
        //         std::cout<<"U10["<< i << ","<< j <<"]= "<< U10[i][j] << " and "<< data_in_u10[360*i+j] <<"\n";
        // }

        // std::vector<double> data_in_u10(360*360);
        // std::vector<double> data_in_v10(360*360);

        // for (int i=0; i<360; ++i)
        // {
        //     for (int j=0; j<360; ++j)
        //     {
        //         data_in_u10[360*i+j] = U10[i][j];
        //         data_in_v10[360*i+j] = V10[i][j];
        //     }
        // }

        auto RX = M_mesh.coordX(diff_angle);
        auto RY = M_mesh.coordY(diff_angle);

        double* data_out_u10;
        double* data_out_v10;
        //int interp_type = TriangleInterpEnum;
        int interp_type = BilinearInterpEnum;
        //int interp_type = NearestInterpEnum;

        InterpFromGridToMeshx(data_out_u10, &X[0], X.size(), &Y[0], Y.size(), &data_in_u10[0], X.size(), Y.size(),
                              &RX[0], &RY[0], M_mesh.numNodes(), 1.0, interp_type);

        InterpFromGridToMeshx(data_out_v10, &X[0], X.size(), &Y[0], Y.size(), &data_in_v10[0], X.size(), Y.size(),
                              &RX[0], &RY[0], M_mesh.numNodes(), 1.0, interp_type);

        // for (int i=0; i<50; ++i)
        //     std::cout<<"data_out["<< i << "]= "<< data_out_u10[i] << " and "<< data_out_v10[i] <<"\n";

        for (int i=0; i<M_num_nodes; ++i)
        {
            fvair[i] = std::cos(-diff_angle)*data_out_u10[i] + std::sin(-diff_angle)*data_out_v10[i];
            fvair[i+M_num_nodes] = -std::sin(-diff_angle)*data_out_u10[i] + std::cos(-diff_angle)*data_out_v10[i];

            // if (i<20)
            //     std::cout<<"data_out["<< i << "]= "<< M_wind[i] << " and "<< M_wind[i+M_num_nodes] <<"\n";
        }

        M_vair[fstep] = fvair;
#if 0

        std::cout<<"MIN BOUND ASRX= "<< *std::min_element(X.begin(),X.end()) <<"\n";
        std::cout<<"MAX BOUND ASRX= "<< *std::max_element(X.begin(),X.end()) <<"\n";

        std::cout<<"MIN BOUND ASRY= "<< *std::min_element(Y.begin(),Y.end()) <<"\n";
        std::cout<<"MAX BOUND ASRY= "<< *std::max_element(Y.begin(),Y.end()) <<"\n";

        std::cout<<"DIFF ANGLE= "<< diff_angle <<"\n";

        std::cout<<"MIN BOUND MESHX= "<< *std::min_element(RX.begin(),RX.end()) <<"\n";
        std::cout<<"MAX BOUND MESHX= "<< *std::max_element(RX.begin(),RX.end()) <<"\n";

        std::cout<<"MIN BOUND MESHY= "<< *std::min_element(RY.begin(),RY.end()) <<"\n";
        std::cout<<"MAX BOUND MESHY= "<< *std::max_element(RY.begin(),RY.end()) <<"\n";
#endif

    }

    // std::cout<<"there are "<<dataFile.getVarCount()<<" variables"<<endl;
    // std::cout<<"there are "<<dataFile.getAttCount()<<" attributes"<<endl;
    // std::cout<<"there are "<<dataFile.getDimCount()<<" dimensions"<<endl;
    // std::cout<<"there are "<<dataFile.getGroupCount()<<" groups"<<endl;
    // std::cout<<"there are "<<dataFile.getTypeCount()<<" types"<<endl;

#endif

}

void
FiniteElement::forcingOcean(bool reload)//(double const& u, double const& v)
{
    switch (M_ocean_type)
    {
        case forcing::OceanType::CONSTANT:
            this->constantOcean(0.,0.);
            break;
        case forcing::OceanType::TOPAZR:
            this->topazOcean(reload);
            break;


        default:
            std::cout << "invalid ocean forcing"<<"\n";
            throw std::logic_error("invalid ocean forcing");
    }
}

void
FiniteElement::constantOcean(double const& u, double const& v)
{
    for (int i=0; i<M_num_nodes; ++i)
    {
        M_ocean[i] = u;
        M_ocean[i+M_num_nodes] = v;
    }
}

void
FiniteElement::topazOcean(bool reload)
{
    if ((current_time < M_ftime_ocean_range[0]) || (M_ftime_ocean_range[1] < current_time) || (current_time == time_init) || reload)
    {
        if (current_time == time_init)
            std::cout<<"load forcing from TOPAZ for initial time\n";
        else
            std::cout<<"forcing not available for the current date: load data from TOPAZ\n";

        this->loadTopazOcean();

        //std::cout<<"forcing not available for the current date\n";
        //throw std::logic_error("forcing not available for the current date");
    }

    double fdt = std::abs(M_ftime_ocean_range[1]-M_ftime_ocean_range[0]);
    std::vector<double> fcoeff(2);
    fcoeff[0] = std::abs(current_time-M_ftime_ocean_range[1])/fdt;
    fcoeff[1] = std::abs(current_time-M_ftime_ocean_range[0])/fdt;

    std::cout<<"TOPAZ LINEAR COEFF 1= "<< fcoeff[0] <<"\n";
    std::cout<<"TOPAZ LINEAR COEFF 2= "<< fcoeff[1] <<"\n";

    for (int i=0; i<M_num_nodes; ++i)
    {
        M_ocean[i] = fcoeff[0]*M_voce[0][i] + fcoeff[1]*M_voce[1][i];
        M_ocean[i+M_num_nodes] = fcoeff[0]*M_voce[0][i+M_num_nodes] + fcoeff[1]*M_voce[1][i+M_num_nodes];

        M_ssh[i] = ssh_coef*M_vssh[i];

        // if (i<20)
        //     std::cout<<"data_out["<< i << "]= "<< M_wind[i] << " and "<< M_wind[i+M_num_nodes] <<"\n";
    }
}

void
FiniteElement::loadTopazOcean()//(double const& u, double const& v)
{

    std::string current_timestr = to_date_string_ym(current_time);
    std::cout<<"TIMESTR= "<< current_timestr <<"\n";
    std::string topaz_filename = (boost::format( "%1%/data/TP4DAILY_%2%_30m.nc" )
                                  % Environment::nextsimDir().string()
                                  % current_timestr ).str();

    std::cout<<"TOPAZ FILE= "<< topaz_filename <<"\n";
    double time_start = std::floor(current_time);
    double time_end = std::ceil(current_time);

    std::cout<<"TOPAZ TIME START= "<< std::setprecision(9) << time_start <<"\n";
    std::cout<<"TOPAZ TIME END  = "<< std::setprecision(9) << time_end <<"\n";

    // We always need at least two time steps to interpolate between
    if (time_end == time_start)
    {
        time_end = time_start + 1.;
    }

    M_ftime_ocean_range.resize(0);
    for (double dt=time_start; dt<=time_end; dt+=1.)
    {
        M_ftime_ocean_range.push_back(dt);
    }

    for (int i=0; i<M_ftime_ocean_range.size(); ++i)
    {
        std::cout<<"TOPAZ TIMEVEC["<< i <<"]= "<< M_ftime_ocean_range[i] <<"\n";
    }

    // if ((current_time < time_start) || (time_end < current_time))
    // {
    //     std::cout<<"forcing not available for the current date\n";
    //     throw std::logic_error("forcing not available for the current date");
    // }

    int nb_forcing_step = M_ftime_ocean_range.size();
    std::cout<<"NB_FORCING_STEP= "<< nb_forcing_step <<"\n";

    // read in re-analysis coordinates
    std::vector<double> XLAT(1101);
    std::vector<double> XLON(1101);

    std::vector<double> YLAT(761);
    std::vector<double> YLON(761);

    std::vector<size_t> index_start(2);
    std::vector<size_t> index_lat_end(2);
    std::vector<size_t> index_lon_end(2);

    std::vector<size_t> index_lat_start(2);
    std::vector<size_t> index_lon_start(2);

    index_lat_start[0] = 0;
    index_lat_start[1] = 0;

    index_lat_end[0] = 1101;
    index_lat_end[1] = 1;

    index_lon_start[0] = 0;
    index_lon_start[1] = 0;

    index_lon_end[0] = 1;
    index_lon_end[1] = 761;

    std::vector<double> XTIME(31);
    std::vector<size_t> index_u_start(4,0);
    std::vector<size_t> index_u_end(4);

    std::vector<size_t> index_ssh_start(3,0);
    std::vector<size_t> index_ssh_end(3);

    std::cout<<"READ NETCDF starts\n";
    netCDF::NcFile dataFile(topaz_filename, netCDF::NcFile::read);
    netCDF::NcVar VXLAT = dataFile.getVar("latitude");
    netCDF::NcVar VXLON = dataFile.getVar("longitude");
    netCDF::NcVar VTIME = dataFile.getVar("time");
    netCDF::NcVar VU = dataFile.getVar("u");
    netCDF::NcVar VV = dataFile.getVar("v");
    netCDF::NcVar VSSH = dataFile.getVar("ssh");
    std::cout<<"READ NETCDF done\n";

    VXLAT.getVar(index_lat_start,index_lat_end,&XLAT[0]);
    VXLON.getVar(index_lat_start,index_lat_end,&XLON[0]);

    VXLAT.getVar(index_lon_start,index_lon_end,&YLAT[0]);
    VXLON.getVar(index_lon_start,index_lon_end,&YLON[0]);

    VTIME.getVar(&XTIME[0]);

    std::vector<double> X(1101);
    std::vector<double> Y(761);

    double RE = 6378.273;
    mapx_class *map;
    std::string configfile = Environment::nextsimDir().string() + "/data/NpsNextsim.mpp";
    std::vector<char> str(configfile.begin(), configfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    for (int i=0; i<1101; ++i)
    {
        X[i] = latLon2XY(XLAT[i], XLON[i], map, configfile)[1];

        // if (i<10)
        // {
        //     std::cout<<"**********************\n";
        //     std::cout<<"X= "<< X[i] <<"\n";
        //     std::cout<<"XLAT= "<< XLAT[i] <<" and XLON= "<< XLON[i] <<"\n";
        // }
    }

    for (int i=0; i<761; ++i)
    {
        Y[i] = latLon2XY(YLAT[i], YLON[i], map, configfile)[0];

        // if (i<10)
        // {
        //     std::cout<<"**********************\n";
        //     std::cout<<"Y= "<< Y[i] <<"\n";
        //     std::cout<<"YLAT= "<< YLAT[i] <<" and YLON= "<< YLON[i] <<"\n";
        // }
    }

    auto RX = M_mesh.coordX();
    auto RY = M_mesh.coordY();

    // std::cout<<"MIN BOUND TOPAZX= "<< *std::min_element(X.begin(),X.end()) <<"\n";
    // std::cout<<"MAX BOUND TOPAZX= "<< *std::max_element(X.begin(),X.end()) <<"\n";

    // std::cout<<"MIN BOUND TOPAZY= "<< *std::min_element(Y.begin(),Y.end()) <<"\n";
    // std::cout<<"MAX BOUND TOPAZY= "<< *std::max_element(Y.begin(),Y.end()) <<"\n";

    // std::cout<<"MIN BOUND MESHX= "<< *std::min_element(RX.begin(),RX.end()) <<"\n";
    // std::cout<<"MAX BOUND MESHX= "<< *std::max_element(RX.begin(),RX.end()) <<"\n";

    // std::cout<<"MIN BOUND MESHY= "<< *std::min_element(RY.begin(),RY.end()) <<"\n";
    // std::cout<<"MAX BOUND MESHY= "<< *std::max_element(RY.begin(),RY.end()) <<"\n";

    std::cout<<"VALUE= "<< from_date_string("1950-01-01") <<"\n";
    std::for_each(XTIME.begin(), XTIME.end(), [&](double& f){ f = f/24.0+from_date_string("1950-01-01"); });
    // for (int i=0; i<31; ++i)
    // {
    //     std::cout<<"TIME["<< i <<"]= "<< XTIME[i] <<"\n";
    // }

    for (int i=0; i<M_ftime_ocean_range.size(); ++i)
    {
        std::cout<<"---TIMEVEC["<< i <<"]= "<< M_ftime_ocean_range[i] <<" : current_time= "<< current_time <<"\n";
    }

    std::vector<double> fvoce(2*M_num_nodes);

    for (int fstep=0; fstep < nb_forcing_step; ++fstep)
    {
        double ftime = M_ftime_ocean_range[fstep];

        if (to_date_string_ym(std::floor(ftime)) != to_date_string_ym(current_time))
        {
            std::string f_timestr = to_date_string_ym(std::floor(ftime));
            std::cout<<"F_TIMESTR= "<< f_timestr <<"\n";
            topaz_filename = (boost::format( "%1%/data/TP4DAILY_%2%_30m.nc" )
                              % Environment::nextsimDir().string()
                              % to_date_string_ym(std::floor(ftime)) ).str();

            std::cout<<"TOPAZ_FILENAME= "<< topaz_filename <<"\n";

            netCDF::NcFile fdataFile(topaz_filename, netCDF::NcFile::read);
            netCDF::NcVar FVTIME = dataFile.getVar("time");
            FVTIME.getVar(&XTIME[0]);
            std::for_each(XTIME.begin(), XTIME.end(), [&](double& f){ f = f/24.0+from_date_string("1950-01-01"); });

            // for (int i=0; i<31; ++i)
            // {
            //     std::cout<<"TIME["<< i <<"]= "<< XTIME[i] <<"\n";
            // }
        }

        auto it = std::find(XTIME.begin(), XTIME.end(), ftime);
        int index = std::distance(XTIME.begin(),it);
        std::cout<<"FIND "<< ftime <<" in index "<< index <<"\n";

        index_u_start[0] = index;
        index_u_start[1] = 0;
        index_u_start[2] = 0;
        index_u_start[3] = 0;

        index_u_end[0] = 1;
        index_u_end[1] = 1;
        index_u_end[2] = 1101;
        index_u_end[3] = 761;

        std::vector<double> data_in_u(1101*761);
        std::vector<double> data_in_v(1101*761);

        VU.getVar(index_u_start,index_u_end,&data_in_u[0]);
        VV.getVar(index_u_start,index_u_end,&data_in_v[0]);

        index_ssh_start[0] = index;
        index_ssh_start[1] = 0;
        index_ssh_start[2] = 0;

        index_ssh_end[0] = 1;
        index_ssh_end[1] = 1101;
        index_ssh_end[2] = 761;

        std::vector<double> data_in_ssh(1101*761);
        VSSH.getVar(index_ssh_start,index_ssh_end,&data_in_ssh[0]);

        // for (int i=0; i<1101; ++i)
        // {
        //     for (int j=0; j<761; ++j)
        //     {
        //         if (i<160 && j<500)
        //             std::cout<<"U["<< i << ","<< j <<"]= "<< data_in_ssh[761*i+j]  <<"\n";
        //         //std::cout<<"U["<< i << ","<< j <<"]= "<< data_in_u[761*i+j] << "  and  " << U[i][j] <<"\n";
        //     }
        // }

        // std::cout<<"MIN DATA U= "<< *std::min_element(data_in_u.begin(),data_in_u.end()) <<"\n";
        // std::cout<<"MAX DATA U= "<< *std::max_element(data_in_u.begin(),data_in_u.end()) <<"\n";

        // std::cout<<"MIN DATA V= "<< *std::min_element(data_in_v.begin(),data_in_v.end()) <<"\n";
        // std::cout<<"MAX DATA V= "<< *std::max_element(data_in_v.begin(),data_in_v.end()) <<"\n";

        double* data_out_u;
        double* data_out_v;
        double* data_out_ssh;
        //int interp_type = TriangleInterpEnum;
        int interp_type = BilinearInterpEnum;
        //int interp_type = NearestInterpEnum;

        InterpFromGridToMeshx(data_out_u, &X[0], X.size(), &Y[0], Y.size(), &data_in_u[0], Y.size(), X.size(),
                              &RX[0], &RY[0], M_mesh.numNodes(), 1.0, interp_type);

        InterpFromGridToMeshx(data_out_v, &X[0], X.size(), &Y[0], Y.size(), &data_in_v[0], Y.size(), X.size(),
                              &RX[0], &RY[0], M_mesh.numNodes(), 1.0, interp_type);

        InterpFromGridToMeshx(data_out_ssh, &X[0], X.size(), &Y[0], Y.size(), &data_in_ssh[0], Y.size(), X.size(),
                              &RX[0], &RY[0], M_mesh.numNodes(), 1.0, interp_type);


        // for (int i=0; i<50; ++i)
        //     std::cout<<"data_out["<< i << "]= "<< data_out_u[i] << " and "<< data_out_v[i] <<"\n";

        // for (int i=0; i<50; ++i)
        //     std::cout<<"data_out["<< i << "]= "<< data_out_ssh[i] <<"\n";

        for (int i=0; i<M_num_nodes; ++i)
        {
            fvoce[i] = data_out_u[i];
            fvoce[i+M_num_nodes] = data_out_v[i];
            M_vssh[i] = data_out_ssh[i];

            // if (i<20)
            //     std::cout<<"data_out["<< i << "]= "<< M_wind[i] << " and "<< M_wind[i+M_num_nodes] <<"\n";
        }

        M_voce[fstep] = fvoce;
    }

    // std::cout<<"there are "<<dataFile.getVarCount()<<" variables"<<endl;
    // std::cout<<"there are "<<dataFile.getAttCount()<<" attributes"<<endl;
    // std::cout<<"there are "<<dataFile.getDimCount()<<" dimensions"<<endl;
    // std::cout<<"there are "<<dataFile.getGroupCount()<<" groups"<<endl;
    // std::cout<<"there are "<<dataFile.getTypeCount()<<" types"<<endl;
}

void
FiniteElement::forcingThermo(double const& u, double const& v)
{
    switch (M_thermo_type)
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
    for (int i=0; i<M_num_nodes; ++i)
    {
        M_thermo[i] = u;
        M_thermo[i+M_num_nodes] = v;
    }
}

void
FiniteElement::initConcentration()
{
    switch (M_conc_type)
    {
        case forcing::ConcentrationType::CONSTANT:
            this->constantConc();
            break;
        case forcing::ConcentrationType::TOPAZ4:
            this->topazConc();
            break;


        default:
            std::cout << "invalid initialization of concentration"<<"\n";
            throw std::logic_error("invalid initialization of concentration");
    }
}

void
FiniteElement::constantConc()
{
    std::fill(M_conc.begin(), M_conc.end(), vm["simul_in.init_concentration"].as<double>());

#if 0
    if (M_water_elements.size() == 0)
    {
        M_water_elements.resize(M_num_elements);
        double welt = 0.;
        int cpt = 0;
        for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
        {
            welt = 0.;
            for (int i=0; i<3; ++i)
            {
                welt += M_mesh.nodes()[it->indices[i]-1].coords[0];
                welt += M_mesh.nodes()[it->indices[i]-1].coords[1];
            }

            M_water_elements[cpt] = welt;

            if (welt >0.)
                M_conc[cpt] = 0.;

            ++cpt;
        }
    }
#endif
}

void
FiniteElement::topazConc()
{
}

void
FiniteElement::initThickness()
{
    switch (M_thick_type)
    {
        case forcing::ThicknessType::CONSTANT:
            this->constantThick();
            break;
        case forcing::ThicknessType::TOPAZ4:
            this->topazThick();
            break;


        default:
            std::cout << "invalid initialization of thickness"<<"\n";
            throw std::logic_error("invalid initialization of thickness");
    }
}

void
FiniteElement::constantThick()
{
    for (int i=0; i<M_num_elements; ++i)
    {
        M_thick[i] = (vm["simul_in.init_thickness"].as<double>())*M_conc[i];
    }
}

void
FiniteElement::topazThick()
{
}

void
FiniteElement::initDamage()
{
    switch (M_damage_type)
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
    std::fill(M_damage.begin(), M_damage.end(), 0.);

#if 0
    for (int i=0; i<M_num_elements; ++i)
    {
        M_damage[i] = 1.0 - M_conc[i];
    }
#endif
}

void
FiniteElement::initSnowThickness()
{
    switch (M_snow_thick_type)
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
    for (int i=0; i<M_num_elements; ++i)
    {
        M_snow_thick[i] = (vm["simul_in.init_snow_thickness"].as<double>())*M_conc[i];
    }
}


void
FiniteElement::initThermodynamics()
{
    // std::fill(damage.begin(), damage.end(), 0.);

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
    switch (M_drifter_type)
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
FiniteElement::bathymetry()
{
    // Interpolation of the bathymetry
    if (vm["simul_in.Lemieux_basal_k2"].as<double>() > 0 )
    {
        double* depth_out;

        // InterpFromMeshToMesh2dx(&depth_out,
        //                         &M_mesh_init.indexTr()[0],&M_mesh_init.coordX()[0],&M_mesh_init.coordY()[0],
        //                         M_mesh_init.numNodes(),M_mesh_init.numTriangles(),
        //                         &M_bathy_depth[0],
        //                         M_mesh_init.numNodes(),1,
        //                         &M_mesh.coordX()[0],&M_mesh.coordY()[0],M_mesh.numNodes(),
        //                         false);

        InterpFromMeshToMesh2dx(&depth_out,
                                &M_mesh_init.indexTr()[0],&M_mesh_init.coordX()[0],&M_mesh_init.coordY()[0],
                                M_mesh_init.numNodes(),M_mesh_init.numTriangles(),
                                &M_bathy_depth[0],
                                M_mesh_init.numTriangles(),1,
                                &M_mesh.bCoordX()[0],&M_mesh.bCoordY()[0],M_mesh.numTriangles(),
                                false);

        M_element_depth.resize(M_num_elements,0.);

        for (int i=0; i<M_element_depth.size(); ++i)
        {
            M_element_depth[i] = depth_out[i];
            //std::cout<<"DEPTH["<< i <<"]= "<< M_element_depth[i] <<"\n";
            //std::cout<<"DEPTH["<< i <<"]= "<< depth_out[i] <<"\n";
            //std::cout<<"DEPTH["<< i <<"]= "<< BC[i] << " and "<< M_mesh.coordX()[i] <<"\n";
        }

        // std::cout<<"MAX BATHY= "<< *std::max_element(M_element_depth.begin(),M_element_depth.end()) <<"\n";
        // std::cout<<"MIN BATHY= "<< *std::min_element(M_element_depth.begin(),M_element_depth.end()) <<"\n";

        // std::cout<<"M_element_depth= "<< M_element_depth.size() <<"\n";
        // std::cout<<"M_element      = "<< M_num_elements <<"\n";

        // std::cout<<"M_element_depth.size()= "<< M_element_depth.size() <<"\n";
        // std::cout<<"M_num_elements        = "<< M_num_elements <<"\n";
        // std::cout<<"M_num_nodes           = "<< M_num_nodes <<"\n";
        // std::cout<<"INTERP DONE\n";

        //this->nodesToElements(depth_out,M_element_depth);

        // for (int k = 0; k < M_num_elements; k++ )
        // {
        //     std::cout<<"DEPTH["<< k <<"]= "<< M_element_depth[k] <<"\n";
        // }

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
                 << "  " << setw(10) << depth_out[k] << "\n";
            //<< "  " << setw(10) << data_in[k] << "\n";
        }
#endif
    }
}

void
FiniteElement::timeInterpolation(int step)
{
    Vair_coef = 1.;
    Voce_coef = 1.;
    ssh_coef = 1.;

    if (((step+1)*time_step) < spinup_duration)
    {
        Vair_coef = ((step+1)*time_step)/spinup_duration;
        Voce_coef = ((step+1)*time_step)/spinup_duration;
        ssh_coef = ((step+1)*time_step)/spinup_duration;
    }
}

void
FiniteElement::nodesToElements(double const* depth, std::vector<double>& v)
{
    //std::vector<double> vc = v;
    //v.resize(M_num_elements);

    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        double sum = 0;
        for (int j=0; j<3; ++j)
        {
            sum += depth[it->indices[j]-1];
        }

        v[cpt] = sum/3.0;
        ++cpt;
    }
}

void
FiniteElement::equallySpacedDrifter()
{
    if (M_drifter.size() ==0)
        M_drifter.resize(M_num_elements);

    std::fill(M_drifter.begin(), M_drifter.end(), 0.);
}

void
FiniteElement::importBamg(BamgMesh const* bamg_mesh)
{
    //mesh_type mesh;
    std::vector<point_type> mesh_nodes;
    std::vector<element_type> mesh_edges;
    std::vector<element_type> mesh_triangles;
    std::vector<double> coords(3,0);

    mesh_nodes.resize(bamg_mesh->VerticesSize[0]);

    for (int id=0; id<bamg_mesh->VerticesSize[0]; ++id)
    {
        coords[0] = bamg_mesh->Vertices[3*id];
        coords[1] = bamg_mesh->Vertices[3*id+1];

        mesh_nodes[id].id = id+1;
        mesh_nodes[id].coords = coords;
    }

    int type = 2;
    int physical = 0;
    int elementary = 0;
    int numVertices = 2;
    std::vector<int> edges(numVertices);

    for (int edg=0; edg<bamg_mesh->EdgesSize[0]; ++edg)
    {
        edges[0] = bamg_mesh->Edges[3*edg];
        edges[1] = bamg_mesh->Edges[3*edg+1];
        //std::cout<< "Edges= "<< bamg_mesh->Edges[3*edg+2] <<"\n";

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

    for (int tr=0; tr<bamg_mesh->TrianglesSize[0]; ++tr)
    {
        indices[0] = bamg_mesh->Triangles[4*tr];
        indices[1] = bamg_mesh->Triangles[4*tr+1];
        indices[2] = bamg_mesh->Triangles[4*tr+2];

        element_type gmshElt( tr,
                              type,
                              physical,
                              elementary,
                              numVertices,
                              indices );

        //mesh_triangles.insert(std::make_pair(tr,gmshElt));
        mesh_triangles.push_back(gmshElt);
    }

    std::cout<<"\n";
    std::cout<<"INFO: Previous  NumNodes     = "<< M_mesh.numNodes() <<"\n";
    std::cout<<"INFO: Previous  NumTriangles = "<< M_mesh.numTriangles() <<"\n";
    std::cout<<"INFO: Previous  NumEdges     = "<< M_mesh.numEdges() <<"\n";

    M_mesh_previous = M_mesh;
    M_mesh = mesh_type(mesh_nodes,mesh_edges,mesh_triangles);
    //M_mesh.writeTofile("out.msh");

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_num_nodes = M_mesh.numNodes();

    std::cout<<"\n";
    std::cout<<"INFO: Current  NumNodes      = "<< M_mesh.numNodes() <<"\n";
    std::cout<<"INFO: Current  NumTriangles  = "<< M_mesh.numTriangles() <<"\n";
    std::cout<<"INFO: Current  NumEdges      = "<< M_mesh.numEdges() <<"\n";
    std::cout<<"\n";
}

std::vector<double>
FiniteElement::latLon2XY(double const& lat, double const& lon, mapx_class* map, std::string const& configfile)
{
    std::vector<double> xy(2);
    double x;
    double y;

    int status = forward_mapx(map,lat,lon,&x,&y);

    xy[0] = x;
    xy[1] = y;

    return xy;
}

double
FiniteElement::latLon2X(double const& lat, double const& lon, mapx_class* map, std::string const& configfile)
{
    double x;
    double y;

    int status = forward_mapx(map,lat,lon,&x,&y);

    return x;
}

double
FiniteElement::latLon2Y(double const& lat, double const& lon, mapx_class* map, std::string const& configfile)
{
    double x;
    double y;

    int status = forward_mapx(map,lat,lon,&x,&y);

    return y;
}

void
FiniteElement::exportResults(int step)
{
    vector_type mx;
    mx.init(3*M_num_elements);
    vector_type my;
    my.init(3*M_num_elements);
    vector_type mc;
    mc.init(3*M_num_elements);
    vector_type mu;
    mu.init(3*M_num_elements);
    vector_type mv;
    mv.init(3*M_num_elements);
    vector_type mfu;
    mfu.init(3*M_num_elements);
    vector_type mfv;
    mfv.init(3*M_num_elements);

    int cpt = 0;
    double sum_u = 0.;
    double sum_v = 0.;
    double sum_fu = 0.;
    double sum_fv = 0.;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        sum_u = 0.;
        sum_v = 0.;
        sum_fu = 0.;
        sum_fv = 0.;
        for (int j=0; j<3; ++j)
        {
            sum_u += M_VT[it->indices[j]-1];
            sum_v += M_VT[it->indices[j]-1+M_num_nodes];
            sum_fu += M_vector->operator()(it->indices[j]-1);
            sum_fv += M_vector->operator()(it->indices[j]-1+M_num_nodes);
        }
        sum_u /= 3.;
        sum_v /= 3.;
        sum_fu /= 3.;
        sum_fv /= 3.;

        for (int i=0; i<3; ++i)
        {
            mc(3*cpt+i) = M_conc[cpt];
            mu(3*cpt+i) = sum_u;
            mv(3*cpt+i) = sum_v;
            mfu(3*cpt+i) = sum_fu;
            mfv(3*cpt+i) = sum_fv;
            mx(3*cpt+i) = M_nodes[it->indices[i]-1].coords[0];
            my(3*cpt+i) = M_nodes[it->indices[i]-1].coords[1];
        }

        ++cpt;
    }

    std::string step_str = boost::lexical_cast<std::string>(step)+".m";
    mx.printMatlab("mx" + step_str);
    my.printMatlab("my" + step_str);
    mc.printMatlab("mc" + step_str);
    mu.printMatlab("mu" + step_str);
    mv.printMatlab("mv" + step_str);
    mfu.printMatlab("mfu" + step_str);
    mfv.printMatlab("mfv" + step_str);
}

} // Nextsim
