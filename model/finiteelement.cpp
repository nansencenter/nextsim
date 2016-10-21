/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   finiteelement.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Mon Aug 24 11:02:45 2015
 */

#include <finiteelement.hpp>
#include <constants.hpp>
#include <date.hpp>
#include <redistribute.hpp>
#include <exporter.hpp>
#include <numeric>

#define GMSH_EXECUTABLE gmsh

namespace Nextsim
{
FiniteElement::FiniteElement()
    :
    vm(Environment::vm()),
    M_mesh(),
    M_solver(),
    M_matrix(),
    M_vector()
{}

// Initialisation of the mesh and forcing
void
FiniteElement::initMesh(setup::DomainType const& domain_type, setup::MeshType const& mesh_type)
{
    switch (domain_type)
    {
        case setup::DomainType::DEFAULT:
            M_flag_fix = 10000; // free = [10001 10002];
            break;
        case setup::DomainType::KARA:
            M_flag_fix = 17; // free = [15 16];
            break;
        case setup::DomainType::BERINGKARA:
            M_flag_fix = 1; // free = [];
            break;
        case setup::DomainType::BIGKARA:
            M_flag_fix = 158; // free = 157;
            break;
        case setup::DomainType::ARCTIC:
            M_flag_fix = 174; // free = [172 173];
            break;
        case setup::DomainType::BIGARCTIC:
            M_flag_fix = 161; // free = 158:160;
            break;
        default:
            std::cout << "invalid domain type"<<"\n";
            throw std::logic_error("invalid domain type");
    }

    M_mesh.readFromFile(M_mesh_filename);

    //if (!vm["simul.wim_grid"].as<bool>())
    // setup the stereographic projection
    //M_mesh.setProjectionFile("toto.mpp");
    M_mesh.stereographicProjection();
    // M_mesh.writeTofile("copy_init_mesh.msh");

    // createGMSHMesh("hypercube.geo");
    // //M_mesh.setOrdering("gmsh");
    // M_mesh.readFromFile("hypercube.msh");

    this->initBamg();

    LOG(DEBUG) <<"Convert MESH starts\n";
    BamgConvertMeshx(
                     bamgmesh,bamggeom,
                     &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                     M_mesh.numNodes(), M_mesh.numTriangles()
                     );
    LOG(DEBUG) <<"Convert MESH done\n";


    for (auto it=M_mesh.edges().begin(), end=M_mesh.edges().end(); it!=end; ++it)
    {
        if (it->physical==M_flag_fix)
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

    // Definition of the hmin, hmax, hminVertices or hmaxVertices
    auto h = this->minMaxSide(M_mesh);

    LOG(DEBUG) <<"MESH: HMIN= "<< h[0] <<"\n";
    LOG(DEBUG) <<"MESH: HMAX= "<< h[1] <<"\n";
    LOG(DEBUG) <<"MESH: RES = "<< this->resolution(M_mesh) <<"\n";

    switch (mesh_type)
    {
        case setup::MeshType::FROM_GMSH:
            // For the other meshes, we use a constant hmin and hmax
            bamgopt->hmin = h[0];
            bamgopt->hmax = h[1];
            break;
        case setup::MeshType::FROM_SPLIT:
            bamgopt->hmin = h[0];
            bamgopt->hmax = h[1];

            M_hminVertices = this->hminVertices(M_mesh_init, bamgmesh);
            M_hmaxVertices = this->hmaxVertices(M_mesh_init, bamgmesh);

            bamgopt->hminVertices = new double[M_mesh_init.numNodes()];
            bamgopt->hmaxVertices = new double[M_mesh_init.numNodes()];
            for (int i=0; i<M_mesh_init.numNodes(); ++i)
            {
                bamgopt->hminVertices[i] = M_hminVertices[i];
                bamgopt->hmaxVertices[i] = M_hmaxVertices[i];
            }
            break;
        default:
            std::cout << "invalid mesh type"<<"\n";
            throw std::logic_error("invalid mesh type");
    }

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_num_nodes = M_mesh.numNodes();
}

// Initialise size of all physical variables with values set to zero
void
FiniteElement::initVariables()
{
    chrono_tot.restart();

    M_solver = solver_ptrtype(new solver_type());
    M_matrix = matrix_ptrtype(new matrix_type());
    M_vector = vector_ptrtype(new vector_type());
    M_solution = vector_ptrtype(new vector_type());

    M_reuse_prec = true;

    //M_matrix->init(2*M_num_nodes,2*M_num_nodes,22);
    M_matrix->init(2*M_num_nodes,2*M_num_nodes,M_graph);
    M_vector->init(2*M_num_nodes);
    M_solution->init(2*M_num_nodes);

    M_VT.resize(2*M_num_nodes,0.);
    M_VTM.resize(2*M_num_nodes,0.);
    M_VTMM.resize(2*M_num_nodes,0.);
    M_node_max_conc.resize(2*M_num_nodes,0.);

    M_sst.resize(M_num_elements);
    M_sss.resize(M_num_elements);

    M_UM.resize(2*M_num_nodes,0.);

    M_h_thin.assign(M_num_elements,0.);
    M_hs_thin.assign(M_num_elements,0.);
    M_tsurf_thin.assign(M_num_elements,0.);

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

    M_sst.resize(M_num_elements);
    M_sss.resize(M_num_elements);

    switch (M_thermo_type)
    {
        case (setup::ThermoType::ZERO_LAYER):
            M_tice.resize(1);
            break;
        case (setup::ThermoType::WINTON):
            M_tice.resize(3);
            break;
        default:
            std::cout << "thermo_type= " << (int)M_thermo_type << "\n";
            throw std::logic_error("Wrong thermo_type");
    }

    for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
        it->assign(M_num_elements,0.);

    for (int i=0; i<M_num_elements; ++i)
    {
        M_divergence_rate[i] = 0.;
        M_sigma[3*i]=0.;
        M_sigma[3*i+1]=0.;
        M_sigma[3*i+2]=0.;
        if ((M_conc[i] <= 0.) || (M_thick[i] <= 0.) )
        {
            M_conc[i] = 0.;
            M_thick[i] = 0.;
        }
    }

    M_surface.resize(M_num_elements);

    M_fcor.resize(M_num_elements);

    M_Cohesion.resize(M_num_elements);
    M_Compressive_strength.resize(M_num_elements);
    M_time_relaxation_damage.resize(M_num_elements,time_relaxation_damage);

    M_tau.resize(2*M_num_nodes,0.);

}//end initVariables

void
FiniteElement::initModelState()
{
    // Initialise the physical state of the model
    this->initIce();

    this->initSlabOcean();

    this->initDrifter();

#if defined (WAVES)
    if (M_use_wim)
        this->initNFloes();
#endif

}

void
FiniteElement::initDatasets()
{
    // Definition of the datasets
    switch(M_atmosphere_type){
        case setup::AtmosphereType::CONSTANT:
            break;

        case setup::AtmosphereType::ASR:
            M_atmosphere_nodes_dataset=DataSet("asr_nodes",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("asr_elements",M_num_elements);
            break;

        case setup::AtmosphereType::ERAi:
            M_atmosphere_nodes_dataset=DataSet("ERAi_nodes",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("ERAi_elements",M_num_elements);
            break;

        case setup::AtmosphereType::EC:
            M_atmosphere_nodes_dataset=DataSet("ec_nodes",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("ec_elements",M_num_elements);
            break;

        case setup::AtmosphereType::EC_ERAi:
            M_atmosphere_nodes_dataset=DataSet("ec_nodes",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("ec_elements",M_num_elements);
            M_atmosphere_bis_elements_dataset=DataSet("ERAi_elements",M_num_elements);
            break;

        default:        std::cout << "invalid wind forcing"<<"\n";throw std::logic_error("invalid wind forcing");
    }

    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
        break;

        case setup::OceanType::TOPAZR:
            M_ocean_nodes_dataset=DataSet("topaz_nodes",M_num_nodes);
            M_ocean_elements_dataset=DataSet("topaz_elements",M_num_elements);
            break;

        case setup::OceanType::TOPAZF:
            M_ocean_nodes_dataset=DataSet("topaz_forecast_nodes",M_num_nodes);
            M_ocean_elements_dataset=DataSet("topaz_forecast_elements",M_num_elements);
            break;

        default:        std::cout << "invalid ocean forcing"<<"\n";throw std::logic_error("invalid wind forcing");
    }

    M_ice_topaz_elements_dataset=DataSet("ice_topaz_elements",M_num_elements);

    M_ice_piomas_elements_dataset=DataSet("ice_piomas_elements",M_num_elements);

    M_ice_amsre_elements_dataset=DataSet("ice_amsre_elements",M_num_elements);

    M_ice_osisaf_elements_dataset=DataSet("ice_osisaf_elements",M_num_elements);

    M_ice_amsr2_elements_dataset=DataSet("ice_amsr2_elements",M_num_elements);

    M_ice_cs2_smos_elements_dataset=DataSet("ice_cs2_smos_elements",M_num_elements);

    M_bathymetry_elements_dataset=DataSet("etopo_elements",M_num_elements);//M_num_nodes);

#if defined (WAVES)
    M_WW3A_elements_dataset=DataSet("ww3a_elements",M_num_elements);
    M_ERAIW_1DEG_elements_dataset=DataSet("erai_waves_1deg_elements",M_num_elements);
#endif

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
    bamgopt->splitcorners      = 1; //the Devil!  Changed to 0, original 1 Phil
    bamgopt->geometricalmetric = 0;
    bamgopt->random            = true;
    bamgopt->verbose           = vm["simul.verbose"].as<int>();

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
    nu0 = vm["simul.nu0"].as<double>();
    young = vm["simul.young"].as<double>();
    rhoi = physical::rhoi;
    rhos = physical::rhos;

    days_in_sec = 24.0*3600.0;
    time_init = from_date_time_string(vm["simul.time_init"].as<std::string>());
    //std::cout<<"time_init second= "<< std::setprecision(18) << time_init <<"\n";

    output_time_step =  days_in_sec/vm["simul.output_per_day"].as<int>();
    mooring_output_time_step =  vm["simul.mooring_output_timestep"].as<double>()*days_in_sec;

    time_step = vm["simul.timestep"].as<double>();
    // output_time_step =  time_step*vm["simul.output_per_day"].as<int>(); // useful for debuging
    duration = (vm["simul.duration"].as<double>())*days_in_sec;
    restart_time_step =  vm["setup.restart_time_step"].as<double>()*days_in_sec;
    M_use_restart   = vm["setup.use_restart"].as<bool>();
    M_write_restart = vm["setup.write_restart"].as<bool>();
    if ( fmod(restart_time_step,time_step) != 0)
    {
        std::cout << restart_time_step << " " << time_step << "\n";
        throw std::runtime_error("restart_time_step not an integer multiple of time_step");
    }

    divergence_min = (1./days_in_sec)*vm["simul.divergence_min"].as<double>();
    compression_factor = vm["simul.compression_factor"].as<double>();
    exponent_compression_factor = vm["simul.exponent_compression_factor"].as<double>();
    ocean_turning_angle_rad = (PI/180.)*vm["simul.oceanic_turning_angle"].as<double>();
    ridging_exponent = vm["simul.ridging_exponent"].as<double>();

    quad_drag_coef_water = vm["simul.quad_drag_coef_water"].as<double>();

    basal_k2 = vm["simul.Lemieux_basal_k2"].as<double>();
    basal_drag_coef_air = vm["simul.Lemieux_drag_coef_air"].as<double>();
    basal_u_0 = vm["simul.Lemieux_basal_u_0"].as<double>();
    basal_Cb = vm["simul.Lemieux_basal_Cb"].as<double>();

    time_relaxation_damage = vm["simul.time_relaxation_damage"].as<double>()*days_in_sec;
    deltaT_relaxation_damage = vm["simul.deltaT_relaxation_damage"].as<double>()*days_in_sec;

    h_thin_max = vm["simul.h_thin_max"].as<double>();
    c_thin_max = vm["simul.c_thin_max"].as<double>();

    compr_strength = vm["simul.compr_strength"].as<double>();
    tract_coef = vm["simul.tract_coef"].as<double>();
    scale_coef = vm["simul.scale_coef"].as<double>();
    alea_factor = vm["simul.alea_factor"].as<double>();
    cfix = vm["simul.cfix"].as<double>();

    C_fix    = cfix*scale_coef;          // C_fix;...  : cohesion (mohr-coulomb) in MPa (40000 Pa)
    C_alea   = alea_factor*C_fix;        // C_alea;... : alea sur la cohesion (Pa)
    tan_phi = vm["simul.tan_phi"].as<double>();
    ridge_h = vm["simul.ridge_h"].as<double>();

    if ( vm["simul.newice_type"].as<int>() == 4 )
        M_ice_cat_type = setup::IceCategoryType::THIN_ICE;
    else
        M_ice_cat_type = setup::IceCategoryType::CLASSIC;

    const boost::unordered_map<const std::string, setup::ThermoType> str2thermo = boost::assign::map_list_of
        ("zero-layer", setup::ThermoType::ZERO_LAYER)
        ("winton", setup::ThermoType::WINTON);
    M_thermo_type = str2thermo.find(vm["setup.thermo-type"].as<std::string>())->second;

    const boost::unordered_map<const std::string, setup::AtmosphereType> str2atmosphere = boost::assign::map_list_of
        ("constant", setup::AtmosphereType::CONSTANT)
        ("asr", setup::AtmosphereType::ASR)
        ("erai", setup::AtmosphereType::ERAi)
        ("ec", setup::AtmosphereType::EC)
        ("ec_erai", setup::AtmosphereType::EC_ERAi);
    M_atmosphere_type = str2atmosphere.find(vm["setup.atmosphere-type"].as<std::string>())->second;

    switch(M_atmosphere_type){
        case setup::AtmosphereType::CONSTANT:   quad_drag_coef_air = vm["simul.ASR_quad_drag_coef_air"].as<double>(); break;
        case setup::AtmosphereType::ASR:        quad_drag_coef_air = vm["simul.ASR_quad_drag_coef_air"].as<double>(); break;
        case setup::AtmosphereType::ERAi:       quad_drag_coef_air = vm["simul.ERAi_quad_drag_coef_air"].as<double>(); break;
        case setup::AtmosphereType::EC:
        case setup::AtmosphereType::EC_ERAi:
                    quad_drag_coef_air = vm["simul.ECMWF_quad_drag_coef_air"].as<double>(); break;
        default:        std::cout << "invalid wind forcing"<<"\n";throw std::logic_error("invalid wind forcing");
    }

    //std::cout<<"AtmosphereType= "<< (int)M_atmosphere_type <<"\n";

    const boost::unordered_map<const std::string, setup::OceanType> str2ocean = boost::assign::map_list_of
        ("constant", setup::OceanType::CONSTANT)
        ("topaz", setup::OceanType::TOPAZR)
        ("topaz_forecast", setup::OceanType::TOPAZF);
    M_ocean_type = str2ocean.find(vm["setup.ocean-type"].as<std::string>())->second;

    //std::cout<<"OCEANTYPE= "<< (int)M_ocean_type <<"\n";

    const boost::unordered_map<const std::string, setup::IceType> str2conc = boost::assign::map_list_of
        ("constant", setup::IceType::CONSTANT)
        ("constant_partial", setup::IceType::CONSTANT_PARTIAL)
        ("target", setup::IceType::TARGET)
        ("topaz", setup::IceType::TOPAZ4)
        ("topaz_forecast", setup::IceType::TOPAZ4F)
        ("topaz_forecast_amsr2", setup::IceType::TOPAZ4FAMSR2)
        ("topaz_forecast_amsr2_osisaf", setup::IceType::TOPAZ4FAMSR2OSISAF)
        ("amsre", setup::IceType::AMSRE)
        ("amsr2", setup::IceType::AMSR2)
        ("osisaf", setup::IceType::OSISAF)
        ("piomas", setup::IceType::PIOMAS)
        ("cs2_smos", setup::IceType::CS2_SMOS);
    M_ice_type = str2conc.find(vm["setup.ice-type"].as<std::string>())->second;

#if defined (WAVES)
    const boost::unordered_map<const std::string, setup::WaveType> str2wave = boost::assign::map_list_of
        ("constant", setup::WaveType::CONSTANT)
        ("ww3a", setup::WaveType::WW3A)
        ("eraiw_1deg", setup::WaveType::ERAI_WAVES_1DEG);
    M_wave_type = str2wave.find(vm["setup.wave-type"].as<std::string>())->second;
    std::cout<<"wave forcing type "<<vm["setup.wave-type"].as<std::string>()<<"\n";
    std::cout<<"wave forcing enum "<<(int)M_wave_type<<"\n";
#endif

    const boost::unordered_map<const std::string, setup::BathymetryType> str2bathymetry = boost::assign::map_list_of
        ("constant", setup::BathymetryType::CONSTANT)
        ("etopo", setup::BathymetryType::ETOPO);
    M_bathymetry_type = str2bathymetry.find(vm["setup.bathymetry-type"].as<std::string>())->second;

    const boost::unordered_map<const std::string, setup::DrifterType> str2drifter = boost::assign::map_list_of
        ("none", setup::DrifterType::NONE)
        ("equallyspaced", setup::DrifterType::EQUALLYSPACED)
        ("iabp", setup::DrifterType::IABP);
    M_drifter_type = str2drifter.find(vm["setup.drifter-type"].as<std::string>())->second;

    //LOG(DEBUG) <<"GMSH VERSION= "<< M_mesh.version() <<"\n";

    M_mesh_filename = vm["simul.mesh_filename"].as<std::string>();

#if 0
    const boost::unordered_map<const std::string, setup::DomainType> str2domain = boost::assign::map_list_of
        ("bigarctic10km.msh", setup::DomainType::BIGARCTIC)
        ("topazreducedsplit2.msh", setup::DomainType::DEFAULT)
        ("topazreducedsplit4.msh", setup::DomainType::DEFAULT)
        ("topazreducedsplit8.msh", setup::DomainType::DEFAULT)
        ("simplesquaresplit2.msh", setup::DomainType::DEFAULT);

    M_domain_type = str2domain.find(M_mesh_filename)->second;

    const boost::unordered_map<const std::string, setup::MeshType> str2mesh = boost::assign::map_list_of
        ("bigarctic10km.msh", setup::MeshType::FROM_GMSH)
        ("topazreducedsplit2.msh", setup::MeshType::FROM_SPLIT)
        ("topazreducedsplit4.msh", setup::MeshType::FROM_SPLIT)
        ("topazreducedsplit8.msh", setup::MeshType::FROM_SPLIT)
        ("simplesquaresplit2.msh", setup::MeshType::FROM_SPLIT);

    M_mesh_type = str2mesh.find(M_mesh_filename)->second;
#endif

    if (M_mesh_filename.find("split") != std::string::npos)
    {
        M_domain_type = setup::DomainType::DEFAULT;
        M_mesh_type = setup::MeshType::FROM_SPLIT;
    }
    else if (M_mesh_filename.find("arctic") != std::string::npos)
    {
        if (M_mesh_filename.find("bigarctic") != std::string::npos)
        {
            M_domain_type = setup::DomainType::BIGARCTIC;
        }
        else
        {
            M_domain_type = setup::DomainType::ARCTIC;
        }

        M_mesh_type = setup::MeshType::FROM_GMSH;
    }

    if (M_mesh_type == setup::MeshType::FROM_SPLIT)
    {
        if (M_mesh_filename.find("wim") == std::string::npos)
        {
            //if "wim" not in name use bamg ordering
            // WIM grids are with gmsh ordering (default)
            M_mesh.setOrdering("bamg");
        }
    }
    else if (M_mesh_type == setup::MeshType::FROM_GMSH)
    {
        if (M_domain_type != setup::DomainType::ARCTIC)
        {
            M_mesh.setOrdering("bamg"); /* The .msh files bigarctic.msh,... that are on Johansen are actually using the bamg ordering*/
        }
        else
            M_mesh.setOrdering("gmsh");
    }
    else
        throw std::logic_error("Unknown setup::MeshType");


    const boost::unordered_map<const std::string, LogLevel> str2log = boost::assign::map_list_of
        ("info", INFO)
        ("warning", WARNING)
        ("debug", DEBUG)
        ("error", ERROR);

    M_log_level = str2log.find(vm["simul.log-level"].as<std::string>())->second;

    M_use_moorings =  vm["simul.use_moorings"].as<bool>();

#if 0
    M_export_path = Environment::nextsimDir().string() + "/matlab";
    // change directory for outputs if the option "output_directory" is not empty
    if ( ! (vm["simul.output_directory"].as<std::string>()).empty() )
    {
        M_export_path = vm["simul.output_directory"].as<std::string>();

        fs::path path(M_export_path);
        // add a subdirecory if needed
        // path /= "subdir";

        // create the output directory if it does not exist
        if ( !fs::exists(path) )
            fs::create_directories(path);
    }
#endif

    // this->writeLogFile(); // already called in the fonction run()
}

void
FiniteElement::createGMSHMesh(std::string const& geofilename)
{
    std::string gmshgeofile = Environment::nextsimDir().string() + "/mesh/" + geofilename;

    if (fs::exists(gmshgeofile))
    {
        //std::cout<<"NOT FOUND " << fs::absolute( gmshgeofile ).string() <<"\n";
        std::ostringstream gmshstr;
        //gmshstr << BOOST_PP_STRINGIZE( GMSH_EXECUTABLE )

        gmshstr << BOOST_PP_STRINGIZE( GMSH_EXECUTABLE )
                << " -" << 2 << " " << gmshgeofile;

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

    side[0] = std::hypot(vertex_1[0]-vertex_0[0], vertex_1[1]-vertex_0[1]);
    side[1] = std::hypot(vertex_2[0]-vertex_1[0], vertex_2[1]-vertex_1[1]);
    side[2] = std::hypot(vertex_2[0]-vertex_0[0], vertex_2[1]-vertex_0[1]);

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
    double min_angle;

#if 0
    int cpt = 0;
    for (auto it=mesh.triangles().begin(), end=mesh.triangles().end(); it!=end; ++it)
    {
        all_min_angle[cpt] = this->minAngles(*it,mesh);
        ++cpt;
    }
#endif

#if 1
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        all_min_angle[cpt] = this->minAngles(M_elements[cpt],mesh);
    }
#endif

    min_angle = *std::min_element(all_min_angle.begin(),all_min_angle.end());
    return min_angle;
}

double
FiniteElement::minAngle(mesh_type const& mesh, std::vector<double> const& um, double factor) const
{
    auto movedmesh = mesh;
    movedmesh.move(um,factor);

    std::vector<double> all_min_angle(movedmesh.numTriangles());

#if 0
    // int cpt = 0;
    // for (auto it=movedmesh.triangles().begin(), end=movedmesh.triangles().end(); it!=end; ++it)
    // {
    //     all_min_angle[cpt] = this->minAngles(*it,movedmesh);
    //     ++cpt;
    // }
#endif

#if 1
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        all_min_angle[cpt] = this->minAngles(movedmesh.triangles()[cpt],movedmesh);
    }
#endif

    return *std::min_element(all_min_angle.begin(),all_min_angle.end());;
}

bool
FiniteElement::flip(mesh_type const& mesh, std::vector<double> const& um, double factor) const
{
    auto movedmesh = mesh;
    movedmesh.move(um,factor);

    std::vector<double> area(movedmesh.numTriangles());

 #if 0
    int cpt = 0;
    for (auto it=movedmesh.triangles().begin(), end=movedmesh.triangles().end(); it!=end; ++it)
    {
        area[cpt] = this->jacobian(*it,movedmesh);
        ++cpt;
    }
#endif

#if 1
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        area[cpt] = this->jacobian(movedmesh.triangles()[cpt],movedmesh);
    }
#endif

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
        hmin[i] = std::sqrt(2.)*std::sqrt(*std::min_element(measure.begin(),measure.end()))*0.9;
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
    int substep_nb=1;
    int step_order=-1;
    bool flip = true;
    int substep = 0;

    std::vector<double> hmin_vertices_first;
    std::vector<double> hmax_vertices_first;

    if (step)
    {
        chrono.restart();
        LOG(DEBUG) <<"Flip starts\n";

        while (flip)
        {
            ++substep;
            displacement_factor /= 2.;
            step_order++;
            flip = this->flip(M_mesh,M_UM,displacement_factor);

            if (substep > 1)
                LOG(DEBUG) <<"FLIP DETECTED "<< substep-1 <<"\n";
        }

        LOG(DEBUG) <<"displacement_factor= "<< displacement_factor <<"\n";

		substep_nb=std::pow(2,step_order);

		if(substep_nb!=1)
		{
			LOG(WARNING) << substep_nb << "substeps will be needed for the remeshing!" <<"\n";
			LOG(WARNING) << "Warning: It is probably due to very high ice speed, check your fields!\n";
		}

        LOG(DEBUG) <<"Flip done in "<< chrono.elapsed() <<"s\n";
    }

	for (int substep_i = 0; substep_i < substep_nb; substep_i++ )
	{
		if(step)
		{
	        chrono.restart();
	        LOG(DEBUG) <<"Move starts\n";
			M_mesh.move(M_UM,displacement_factor);
			LOG(DEBUG) <<"Move done in "<< chrono.elapsed() <<"s\n";



	        chrono.restart();
	        LOG(DEBUG) <<"Move bamgmesh->Vertices starts\n";
		    auto RX = M_mesh.coordX();
		    auto RY = M_mesh.coordY();

			for (int id=0; id<bamgmesh->VerticesSize[0]; ++id)
	    	{
	        	bamgmesh->Vertices[3*id] = RX[id];
				bamgmesh->Vertices[3*id+1] = RY[id] ;
	    	}
			LOG(DEBUG) <<"Move bamgmesh->Vertices done in "<< chrono.elapsed() <<"s\n";
		}


		if(M_mesh_type==setup::MeshType::FROM_SPLIT)
		{
			if(step==0)
			{
		        chrono.restart();
		        LOG(DEBUG) <<"First adaptation starts\n";
				// step 1 (only for the first time step): Start by having bamg 'clean' the mesh with KeepVertices=0
				bamgopt->KeepVertices=0;
				this->adaptMesh();
				bamgopt->KeepVertices=1;
				LOG(DEBUG) <<"First adaptation done in "<< chrono.elapsed() <<"s\n";

			}

	        chrono.restart();
	        LOG(DEBUG) <<"Interpolate hminVertices starts\n";
			// Interpolate hminVertices and hmaxVertices onto the current mesh

			// NODAL INTERPOLATION
			int init_num_nodes = M_mesh_init.numNodes();

			// To avoid memory leak:
			std::vector<double> interp_Vertices_in(2*init_num_nodes);

			double* interp_Vertices_out;

			for (int i=0; i<init_num_nodes; ++i)
			{
				interp_Vertices_in[2*i]   = M_hminVertices[i];
				interp_Vertices_in[2*i+1] = M_hmaxVertices[i];
			}

			InterpFromMeshToMesh2dx(&interp_Vertices_out,
                                    &M_mesh_init.indexTr()[0],&M_mesh_init.coordX()[0],&M_mesh_init.coordY()[0],
                                    M_mesh_init.numNodes(),M_mesh_init.numTriangles(),
                                    &interp_Vertices_in[0],
                                    M_mesh_init.numNodes(),2,
                                    &M_mesh.coordX()[0],&M_mesh.coordY()[0],M_mesh.numNodes(),
                                    false);

            // No need to deallocate the memory related to hminVertices and hmaxVertices,
            // as it is done when deleting bamgopt_previous in adaptMesh
			bamgopt->hminVertices = new double[M_mesh.numNodes()];
			bamgopt->hmaxVertices = new double[M_mesh.numNodes()];

			for (int i=0; i<M_mesh.numNodes(); ++i)
			{
				bamgopt->hminVertices[i] = interp_Vertices_out[2*i];
				bamgopt->hmaxVertices[i] = interp_Vertices_out[2*i+1];
			}

			xDelete<double>(interp_Vertices_out);
			LOG(DEBUG) <<"Interpolate hmin done in "<< chrono.elapsed() <<"s\n";
		}


        if(step && (vm["simul.regrid_output_flag"].as<bool>()))
        {
            had_remeshed=true;
            mesh_adapt_step++;
            this->exportResults(200000+mesh_adapt_step);
		}

        chrono.restart();
        LOG(INFO) <<"AdaptMesh starts\n";
        this->adaptMesh();
	    LOG(INFO) <<"AdaptMesh done in "<< chrono.elapsed() <<"s\n";


		if (step)
		{
			int prv_num_elements = M_mesh_previous.numTriangles();
			int prv_num_nodes = M_mesh_previous.numNodes();

	        chrono.restart();
	        LOG(DEBUG) <<"Element Interp starts\n";
			// ELEMENT INTERPOLATION With Cavities
			int nb_var=14 + M_tice.size();

#if defined (WAVES)
            // coupling with wim
            // - only interpolate if not at a coupling time step
            // - else nfloes will just be overwritten with wimToNextsim()
            bool nfloes_interp = (M_use_wim && (!M_run_wim));

            if (nfloes_interp)
                std::cout<<"IN REGRID: "<< "interpolate nfloes\n";
            else
                std::cout<<"IN REGRID: "<< "do not interpolate nfloes\n";

            if (nfloes_interp)
                nb_var++;
#endif

			// To avoid memory leak:
			std::vector<double> interp_elt_in(nb_var*prv_num_elements);

			double* interp_elt_out;

			LOG(DEBUG) <<"ELEMENT: Interp starts\n";

			int tmp_nb_var=0;
			for (int i=0; i<prv_num_elements; ++i)
			{
				tmp_nb_var=0;

				// concentration
				interp_elt_in[nb_var*i+tmp_nb_var] = M_conc[i];
				tmp_nb_var++;

				// thickness
				interp_elt_in[nb_var*i+tmp_nb_var] = M_thick[i];
				tmp_nb_var++;

				// snow thickness
				interp_elt_in[nb_var*i+tmp_nb_var] = M_snow_thick[i];
				tmp_nb_var++;

				// integrated_stress1
				interp_elt_in[nb_var*i+tmp_nb_var] = M_sigma[3*i]*M_thick[i];
				tmp_nb_var++;

				// integrated_stress2
				interp_elt_in[nb_var*i+tmp_nb_var] = M_sigma[3*i+1]*M_thick[i];
				tmp_nb_var++;

				// integrated_stress3
				interp_elt_in[nb_var*i+tmp_nb_var] = M_sigma[3*i+2]*M_thick[i];
				tmp_nb_var++;

				// compliance
				interp_elt_in[nb_var*i+tmp_nb_var] = 1./(1.-M_damage[i]);
				tmp_nb_var++;

				// divergence_rate
				interp_elt_in[nb_var*i+tmp_nb_var] = M_divergence_rate[i];
				tmp_nb_var++;

				// h_ridged_thin_ice
				interp_elt_in[nb_var*i+tmp_nb_var] = M_h_ridged_thin_ice[i];
				tmp_nb_var++;

				// h_ridged_thick_ice
				interp_elt_in[nb_var*i+tmp_nb_var] = M_h_ridged_thick_ice[i];
				tmp_nb_var++;

				// random_number
				interp_elt_in[nb_var*i+tmp_nb_var] = M_random_number[i];
				tmp_nb_var++;

				// Ice temperature
                interp_elt_in[nb_var*i+tmp_nb_var] = M_tice[0][i];
                tmp_nb_var++;
                if ( M_thermo_type == setup::ThermoType::WINTON )
                {
                    interp_elt_in[nb_var*i+tmp_nb_var] = ( M_tice[1][i] - physical::mu*physical::si*physical::Lf/(physical::C*M_tice[1][i]) ) * M_thick[i]; // (39) times volume with f1=1
                    tmp_nb_var++;
                    interp_elt_in[nb_var*i+tmp_nb_var] = ( M_tice[2][i] ) * M_thick[i]; // (39) times volume with f1=0
                    tmp_nb_var++;
                }

				// thin ice thickness
				interp_elt_in[nb_var*i+tmp_nb_var] = M_h_thin[i];
				tmp_nb_var++;

				// snow on thin ice
				interp_elt_in[nb_var*i+tmp_nb_var] = M_hs_thin[i];
				tmp_nb_var++;

				// Ice surface temperature for thin ice
				interp_elt_in[nb_var*i+tmp_nb_var] = M_tsurf_thin[i];
				tmp_nb_var++;

#if defined (WAVES)
                // Nfloes from wim model
                if (nfloes_interp)
                {
                    interp_elt_in[nb_var*i+tmp_nb_var] = M_nfloes[i];
                    tmp_nb_var++;
                }
#endif

				if(tmp_nb_var>nb_var)
				{
					throw std::logic_error("tmp_nb_var not equal to nb_var");
				}
			}

			// To avoid memory leak:
			std::vector<double> surface_previous(prv_num_elements);
			std::vector<double> surface(M_num_elements);

			int cpt = 0;
			for (auto it=M_mesh_previous.triangles().begin(), end=M_mesh_previous.triangles().end(); it!=end; ++it)
			{
				surface_previous[cpt] = this->measure(*it,M_mesh_previous);
				++cpt;
			}

			cpt = 0;
			for (auto it=M_mesh.triangles().begin(), end=M_mesh.triangles().end(); it!=end; ++it)
			{
				surface[cpt] = this->measure(*it,M_mesh);
				++cpt;
			}

			// The interpolation with the cavities still needs to be tested on a long run.
			// By default, we then use the non-conservative MeshToMesh interpolation

			InterpFromMeshToMesh2dCavities(&interp_elt_out,&interp_elt_in[0],nb_var,
                                           &surface_previous[0], &surface[0], bamgmesh_previous, bamgmesh);

#if 0
			InterpFromMeshToMesh2dx(&interp_elt_out,
                                    &M_mesh_previous.indexTr()[0],&M_mesh_previous.coordX()[0],&M_mesh_previous.coordY()[0],
                                    M_mesh_previous.numNodes(),M_mesh_previous.numTriangles(),
                                    &interp_elt_in[0],
                                    M_mesh_previous.numTriangles(),nb_var,
                                    &M_mesh.bcoordX()[0],&M_mesh.bcoordY()[0],M_mesh.numTriangles(),
                                    false);
#endif

			M_conc.assign(M_num_elements,0.);
			M_thick.assign(M_num_elements,0.);
			M_snow_thick.assign(M_num_elements,0.);
			M_sigma.assign(3*M_num_elements,0.);
			M_damage.assign(M_num_elements,0.);

			M_divergence_rate.assign(M_num_elements,0.);
			M_h_ridged_thin_ice.assign(M_num_elements,0.);
			M_h_ridged_thick_ice.assign(M_num_elements,0.);

			M_random_number.resize(M_num_elements);

            for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
                it->assign(M_num_elements,0.);

            M_h_thin.assign(M_num_elements,0.);
            M_hs_thin.assign(M_num_elements,0.);
            M_tsurf_thin.assign(M_num_elements,0.);

#if defined (WAVES)
            if (nfloes_interp)
                M_nfloes.assign(M_num_elements,0.);
#endif

			for (int i=0; i<M_num_elements; ++i)
			{
				tmp_nb_var=0;

				// concentration
				M_conc[i] = std::max(0., std::min(1.,interp_elt_out[nb_var*i+tmp_nb_var]));
				tmp_nb_var++;

				// thickness
				M_thick[i] = std::max(0., interp_elt_out[nb_var*i+tmp_nb_var]);
				tmp_nb_var++;

				// snow thickness
				M_snow_thick[i] = std::max(0., interp_elt_out[nb_var*i+tmp_nb_var]);
				tmp_nb_var++;

				if (M_thick[i] != 0.)
				{
					// integrated_stress1
					M_sigma[3*i] = interp_elt_out[nb_var*i+tmp_nb_var]/M_thick[i];
					tmp_nb_var++;

					// integrated_stress2
					M_sigma[3*i+1] = interp_elt_out[nb_var*i+tmp_nb_var]/M_thick[i];
					tmp_nb_var++;

					// integrated_stress3
					M_sigma[3*i+2] = interp_elt_out[nb_var*i+tmp_nb_var]/M_thick[i];
					tmp_nb_var++;
				}
				else
				{
					tmp_nb_var+=3;
				}

				// compliance
				if (interp_elt_out[nb_var*i+tmp_nb_var] != 0.)
				{
					M_damage[i] = std::max(0., std::min(1.,1.-1./interp_elt_out[nb_var*i+tmp_nb_var]));
					tmp_nb_var++;
				}
				else
				{
					M_damage[i] = 0.;
					tmp_nb_var++;
				}

				// divergence_rate
				M_divergence_rate[i] = interp_elt_out[nb_var*i+tmp_nb_var];
				tmp_nb_var++;

				// // h_ridged_thin_ice
				M_h_ridged_thin_ice[i] = interp_elt_out[nb_var*i+tmp_nb_var];
				tmp_nb_var++;

				// h_ridged_thick_ice
				M_h_ridged_thick_ice[i] = interp_elt_out[nb_var*i+tmp_nb_var];
				tmp_nb_var++;

				// random_number
				M_random_number[i] = interp_elt_out[nb_var*i+tmp_nb_var];
				//M_random_number[i] = std::max(0., std::min(1.,interp_elt_in[11*i+tmp_nb_var]));
				tmp_nb_var++;

				// Ice temperature
                M_tice[0][i] = interp_elt_out[nb_var*i+tmp_nb_var];
                tmp_nb_var++;
                if ( M_thermo_type == setup::ThermoType::WINTON )
                {
                    double tmp = interp_elt_out[nb_var*i+tmp_nb_var]/M_thick[i];
                    M_tice[1][i] = 0.5*( tmp - std::sqrt(tmp*tmp + 4*physical::mu*physical::si*physical::Lf/physical::C) ); // (38) divided with volume with f1=1
                    tmp_nb_var++;
                    M_tice[2][i] = interp_elt_out[nb_var*i+tmp_nb_var]/M_thick[i]; // (40) divided with volume with f1=0
                    tmp_nb_var++;
                }

				// thin ice thickness
				M_h_thin[i] = interp_elt_out[nb_var*i+tmp_nb_var];
				tmp_nb_var++;

				// snow on thin ice
				M_hs_thin[i] = interp_elt_out[nb_var*i+tmp_nb_var];
				tmp_nb_var++;

				// Ice surface temperature for thin ice
				M_tsurf_thin[i] = interp_elt_out[nb_var*i+tmp_nb_var];
				tmp_nb_var++;

#if defined (WAVES)
                // Nfloes from wim model
                if (nfloes_interp)
                {
                    M_nfloes[i] = interp_elt_out[nb_var*i+tmp_nb_var];
                    tmp_nb_var++;
                }

#endif

				if(tmp_nb_var!=nb_var)
				{
					throw std::logic_error("tmp_nb_var not equal to nb_var");
				}
			}

			xDelete<double>(interp_elt_out);

			LOG(DEBUG) <<"ELEMENT: Interp done\n";
			LOG(DEBUG) <<"Element Interp done in "<< chrono.elapsed() <<"s\n";

	        chrono.restart();
	        LOG(DEBUG) <<"Slab Interp starts\n";
			// ELEMENT INTERPOLATION FOR SLAB OCEAN FROM OLD MESH ON ITS ORIGINAL POSITION
			nb_var=2;

			// memory leak:
			std::vector<double> interp_elt_slab_in(nb_var*prv_num_elements);

			double* interp_elt_slab_out;

			LOG(DEBUG) <<"ELEMENT SLAB: Interp starts\n";

			M_mesh_previous.move(M_UM,-displacement_factor);

			for (int i=0; i<prv_num_elements; ++i)
			{
				// Sea surface temperature
				interp_elt_slab_in[nb_var*i+0] = M_sst[i];

				// Sea surface salinity
				interp_elt_slab_in[nb_var*i+1] = M_sss[i];
			}

			InterpFromMeshToMesh2dx(&interp_elt_slab_out,
                                    &M_mesh_previous.indexTr()[0],&M_mesh_previous.coordX()[0],&M_mesh_previous.coordY()[0],
                                    M_mesh_previous.numNodes(),M_mesh_previous.numTriangles(),
                                    &interp_elt_slab_in[0],
                                    M_mesh_previous.numTriangles(),nb_var,
                                    &M_mesh.bcoordX()[0],&M_mesh.bcoordY()[0],M_mesh.numTriangles(),
                                    false);

			M_sst.resize(M_num_elements);
			M_sss.resize(M_num_elements);

			for (int i=0; i<M_mesh.numTriangles(); ++i)
			{
				// Sea surface temperature
				M_sst[i] = interp_elt_slab_out[nb_var*i+0];

				// Sea surface salinity
				M_sss[i] = interp_elt_slab_out[nb_var*i+1];
			}

			xDelete<double>(interp_elt_slab_out);

			M_mesh_previous.move(M_UM,displacement_factor);
			LOG(DEBUG) <<"ELEMENT SLAB: Interp done\n";
			LOG(DEBUG) <<"Slab Interp done in "<< chrono.elapsed() <<"s\n";

			// NODAL INTERPOLATION
			nb_var=8;

			std::vector<double> interp_in(nb_var*prv_num_nodes);

			double* interp_out;

	        chrono.restart();
	        LOG(DEBUG) <<"Nodal Interp starts\n";
			LOG(DEBUG) <<"NODAL: Interp starts\n";

			for (int i=0; i<prv_num_nodes; ++i)
			{
				// VT
				interp_in[nb_var*i] = M_VT[i];
				interp_in[nb_var*i+1] = M_VT[i+prv_num_nodes];

				// VTM
				interp_in[nb_var*i+2] = M_VTM[i];
				interp_in[nb_var*i+3] = M_VTM[i+prv_num_nodes];

				// VTMM
				interp_in[nb_var*i+4] = M_VTMM[i];
				interp_in[nb_var*i+5] = M_VTMM[i+prv_num_nodes];

				// UM
				interp_in[nb_var*i+6] = M_UM[i];
				interp_in[nb_var*i+7] = M_UM[i+prv_num_nodes];
			}

			InterpFromMeshToMesh2dx(&interp_out,
                                    &M_mesh_previous.indexTr()[0],&M_mesh_previous.coordX()[0],&M_mesh_previous.coordY()[0],
                                    M_mesh_previous.numNodes(),M_mesh_previous.numTriangles(),
                                    &interp_in[0],
                                    M_mesh_previous.numNodes(),nb_var,
                                    &M_mesh.coordX()[0],&M_mesh.coordY()[0],M_mesh.numNodes(),
                                    false);

			M_VT.assign(2*M_num_nodes,0.);
			M_VTM.assign(2*M_num_nodes,0.);
			M_VTMM.assign(2*M_num_nodes,0.);
			M_UM.assign(2*M_num_nodes,0.);

			for (int i=0; i<M_num_nodes; ++i)
			{
				// VT
				M_VT[i] = interp_out[nb_var*i];
				M_VT[i+M_num_nodes] = interp_out[nb_var*i+1];

				// VTM
				M_VTM[i] = interp_out[nb_var*i+2];
				M_VTM[i+M_num_nodes] = interp_out[nb_var*i+3];

				// VTMM
				M_VTMM[i] = interp_out[nb_var*i+4];
				M_VTMM[i+M_num_nodes] = interp_out[nb_var*i+5];

				// UM
				M_UM[i] = interp_out[nb_var*i+6];
				M_UM[i+M_num_nodes] = interp_out[nb_var*i+7];
			}

			xDelete<double>(interp_out);

			LOG(DEBUG) <<"NODAL: Interp done\n";
			LOG(DEBUG) <<"Nodal interp done in "<< chrono.elapsed() <<"s\n";

            // Drifters - if requested
            if ( M_drifter_type != setup::DrifterType::NONE )
            {
                chrono.restart();
                LOG(DEBUG) <<"Drifter starts\n";
                LOG(DEBUG) <<"DRIFTER: Interp starts\n";

                // Assemble the coordinates from the unordered_map
                std::vector<double> drifter_X(M_drifter.size());
                std::vector<double> drifter_Y(M_drifter.size());
                int j=0;
                for ( auto it = M_drifter.begin(); it != M_drifter.end(); ++it )
                {
                    drifter_X[j] = it->second[0];
                    drifter_Y[j] = it->second[1];
                    ++j;
                }

                // Interpolate the velocity and concentration onto the drifter positions
                nb_var=2;
                std::vector<double> interp_drifter_in(nb_var*M_mesh.numNodes());

                // Interpolate the velocity
                for (int i=0; i<M_mesh.numNodes(); ++i)
                {
                    interp_drifter_in[nb_var*i]   = M_UM[i];
                    interp_drifter_in[nb_var*i+1] = M_UM[i+M_mesh.numNodes()];
                }

                double* interp_drifter_out;
                InterpFromMeshToMesh2dx(&interp_drifter_out,
                                        &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                                        M_mesh.numNodes(),M_mesh.numTriangles(),
                                        &interp_drifter_in[0],
                                        M_mesh.numNodes(),nb_var,
                                        &drifter_X[0],&drifter_Y[0],M_drifter.size(),
                                        true, 0.);


                // Interpolate the concentration - re-use interp_drifter_in
                interp_drifter_in.resize(M_mesh.numTriangles());
                for (int i=0; i<M_mesh.numTriangles(); ++i)
                {
                    interp_drifter_in[i]   = M_conc[i];
                }

                double* interp_drifter_c_out;
                InterpFromMeshToMesh2dx(&interp_drifter_c_out,
                                        &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                                        M_mesh.numNodes(),M_mesh.numTriangles(),
                                        &interp_drifter_in[0],
                                        M_mesh.numTriangles(),1,
                                        &drifter_X[0],&drifter_Y[0],M_drifter.size(),
                                        true, 0.);

                // Rebuild the M_drifter map
                double clim = vm["simul.drift_limit_concentration"].as<double>();
                j=0;
                for ( auto it = M_drifter.begin(); it != M_drifter.end(); /* ++it is not allowed here, because we use 'erase' */ )
                {
                    if ( interp_drifter_c_out[j] > clim )
                    {
                        M_drifter[it->first] = std::array<double,2> {it->second[0]+interp_drifter_out[j], it->second[1]+interp_drifter_out[j+M_drifter.size()]};
                        ++it;
                    } else {
                        // Throw out drifters that drift out of the ice
                        it = M_drifter.erase(it);
                    }
                    ++j;
                }

                xDelete<double>(interp_drifter_out);
                xDelete<double>(interp_drifter_c_out);

                LOG(DEBUG) <<"DRIFTER: Interp done\n";
                LOG(DEBUG) <<"Drifter interp done in "<< chrono.elapsed() <<"s\n";
            }
		}
	}

    if (step)
    {
		for (int i=0; i<M_num_nodes; ++i)
		{
			// UM
			M_UM[i] = 0.;
			M_UM[i+M_num_nodes] = 0.;
		}

        //M_matrix->init(2*M_num_nodes,2*M_num_nodes,22);
        M_matrix->init(2*M_num_nodes,2*M_num_nodes,M_graph);
        M_vector->resize(2*M_num_nodes);
        M_solution->resize(2*M_num_nodes);
        M_reuse_prec = false;

        M_node_max_conc.resize(2*M_num_nodes,0.);

        M_fcor.resize(M_num_elements);
    }

    M_atmosphere_nodes_dataset.target_size=M_num_nodes;
    M_atmosphere_elements_dataset.target_size=M_num_elements;
    M_atmosphere_bis_elements_dataset.target_size=M_num_elements;
    M_ocean_nodes_dataset.target_size=M_num_nodes;
    M_ocean_elements_dataset.target_size=M_num_elements;

    M_ice_topaz_elements_dataset.target_size=M_num_elements;
    M_ice_piomas_elements_dataset.target_size=M_num_elements;
    M_ice_amsre_elements_dataset.target_size=M_num_elements;
    M_ice_osisaf_elements_dataset.target_size=M_num_elements;
    M_ice_amsr2_elements_dataset.target_size=M_num_elements;
    M_ice_cs2_smos_elements_dataset.target_size=M_num_elements;
    M_bathymetry_elements_dataset.target_size=M_num_elements;
#if defined (WAVES)
    M_WW3A_elements_dataset.target_size=M_num_elements;
    M_ERAIW_1DEG_elements_dataset.target_size=M_num_elements;
#endif


    M_atmosphere_nodes_dataset.reloaded=false;
    M_atmosphere_elements_dataset.reloaded=false;
    M_atmosphere_bis_elements_dataset.reloaded=false;
    M_ocean_nodes_dataset.reloaded=false;
    M_ocean_elements_dataset.reloaded=false;

    M_ice_topaz_elements_dataset.reloaded=false;
    M_ice_piomas_elements_dataset.reloaded=false;
    M_ice_amsre_elements_dataset.reloaded=false;
    M_ice_osisaf_elements_dataset.reloaded=false;
    M_ice_amsr2_elements_dataset.reloaded=false;
    M_ice_cs2_smos_elements_dataset.reloaded=false;
    M_bathymetry_elements_dataset.reloaded=false;
#if defined (WAVES)
    M_WW3A_elements_dataset.reloaded=false;
    M_ERAIW_1DEG_elements_dataset.reloaded=false;
#endif

    // for the parallel code, it will be necessary to add those lines
    // as the domain covered by the partitions changes at each remeshing/partitioning
#if 0
    M_atmosphere_nodes_dataset.grid.loaded=false;
    M_atmosphere_elements_dataset.grid.loaded=false;
    M_atmosphere_bis_elements_dataset.grid.loaded=false;
    M_ocean_nodes_dataset.grid.loaded=false;
    M_ocean_elements_dataset.grid.loaded=false;
    M_ice_topaz_elements_dataset.grid.loaded=false;
    M_ice_amsre_elements_dataset.grid.loaded=false;
    M_ice_osisaf_elements_dataset.grid.loaded=false;
    M_ice_amsr2_elements_dataset.grid.loaded=false;
    M_ice_cs2_smos_elements_dataset.grid.loaded=false;
    M_bathymetry_elements_dataset.grid.loaded=false;
#endif

    M_Cohesion.resize(M_num_elements);
    M_Compressive_strength.resize(M_num_elements);
    M_time_relaxation_damage.resize(M_num_elements,time_relaxation_damage);
}

void
FiniteElement::adaptMesh()
{
    delete bamgopt_previous;
    delete bamggeom_previous;
    delete bamgmesh_previous;

    bamgopt_previous = new BamgOpts();
    bamggeom_previous = new BamgGeom();
    bamgmesh_previous = new BamgMesh();

    *bamgmesh_previous = *bamgmesh;
    *bamggeom_previous = *bamggeom;
    *bamgopt_previous = *bamgopt;

    int fnd = 0;

    for (int edg=0; edg<bamgmesh_previous->EdgesSize[0]; ++edg)
    {
        fnd = bamgmesh_previous->Edges[3*edg]-1;

        if ((std::binary_search(M_dirichlet_flags.begin(),M_dirichlet_flags.end(),fnd)))
        {
            bamggeom_previous->Edges[3*edg+2] = M_flag_fix;
            bamgmesh_previous->Edges[3*edg+2] = M_flag_fix;
        }
    }

    //Environment::logMemoryUsage("before adaptMesh");
    Bamgx(bamgmesh,bamggeom,bamgmesh_previous,bamggeom_previous,bamgopt_previous);
    //Environment::logMemoryUsage("after adaptMesh");

    // Save the old id_node before redefining it
    std::vector<int> old_node_id=M_mesh.id();

    // Import the mesh from bamg
    this->importBamg(bamgmesh);

        // Recompute the node ids
    if(bamgopt->KeepVertices)
    {
        std::vector<int> new_nodes_id=M_mesh.id();

        int Boundary_id=0;
        int nb_new_nodes=0;

        // We mask out the boundary nodes
        std::vector<bool> mask(bamgmesh->VerticesSize[0],false) ;
        for (int vert=0; vert<bamgmesh->VerticesOnGeomVertexSize[0]; ++vert)
            mask[bamgmesh->VerticesOnGeomVertex[2*vert]-1]=true; // The factor 2 is because VerticesOnGeomVertex has 2 dimensions in bamg

        // The new id will have values higher than the previous one
        int first_new_node=*std::max_element(old_node_id.begin(),old_node_id.end())+1;

        for (int vert=0; vert<bamgmesh->VerticesSize[0]; ++vert)
        {
            if(mask[vert])
            {
                Boundary_id++;
                new_nodes_id[vert]=Boundary_id;
            }
            else
            {
                if(bamgmesh->PreviousNumbering[vert]==0)
                {
                    nb_new_nodes++;
                    new_nodes_id[vert]=first_new_node+nb_new_nodes-1;
                }
                else
                {
                    new_nodes_id[vert]=old_node_id[bamgmesh->PreviousNumbering[vert]-1];
                }
            }
        }
        M_mesh.setId(new_nodes_id);
    }


    // update dirichlet nodes
    M_boundary_flags.resize(0);
    M_dirichlet_flags.resize(0);
    for (int edg=0; edg<bamgmesh->EdgesSize[0]; ++edg)
    {
        M_boundary_flags.push_back(bamgmesh->Edges[3*edg]-1);
        //M_boundary_flags.push_back(bamgmesh->Edges[3*edg+1]-1);

        if (bamgmesh->Edges[3*edg+2] == M_flag_fix)
        {
            M_dirichlet_flags.push_back(bamgmesh->Edges[3*edg]-1);
            //M_dirichlet_flags.push_back(bamgmesh->Edges[3*edg+1]-1);

            //std::cout<<"NODES["<< edg <<"]= "<< bamgmesh->Edges[3*edg]-1 << " and "<< bamgmesh->Edges[3*edg+1]-1 <<"\n";
        }
    }

    std::sort(M_dirichlet_flags.begin(), M_dirichlet_flags.end());
    //M_dirichlet_flags.erase( std::unique(M_dirichlet_flags.begin(), M_dirichlet_flags.end() ), M_dirichlet_flags.end());

    std::sort(M_boundary_flags.begin(), M_boundary_flags.end());
    //M_boundary_flags.erase(std::unique(M_boundary_flags.begin(), M_boundary_flags.end() ), M_boundary_flags.end());

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

    M_surface.assign(M_num_elements,0.);

    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        M_surface[cpt] = this->measure(*it,M_mesh);
        ++cpt;
    }
}

void
FiniteElement::assemble(int pcpt)
{
    M_matrix->zero();
    M_vector->zero();

    std::vector<int> rhsindices(2*M_num_nodes);
    std::iota(rhsindices.begin(), rhsindices.end(), 0);
    std::vector<double> rhsdata(2*M_num_nodes, 0.);

	M_node_max_conc.assign(2*M_num_nodes,0.);

#if 1
    std::vector<double> lhsdata(M_graph.ja().size(), 0.);
#endif

    std::vector<int> extended_dirichlet_nodes = M_dirichlet_nodes;

    // ---------- Identical values for all the elements -----------
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

    double cos_ocean_turning_angle=std::cos(ocean_turning_angle_rad);
    double sin_ocean_turning_angle=std::sin(ocean_turning_angle_rad);

    // ---------- Assembling starts -----------
    LOG(DEBUG) <<"Assembling starts\n";
    chrono.restart();

    int thread_id;
    int total_threads;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    //std::cout<<"MAX THREADS= "<< max_threads <<"\n";

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        // if(thread_id == 0)
        // {
        //     std::cout<< "\nMaster thread being called\n";
        //     total_threads = omp_get_num_threads();
        //     std::cout<<"Total number of threads are "<< total_threads <<"\n";
        // }

        /* Compute the value that only depends on the element */
        double welt_oce_ice = 0.;
        double welt_air_ice = 0.;
        double welt_ice = 0.;
        double welt_ssh = 0.;
        int nind;

        for (int i=0; i<3; ++i)
        {
            nind = (M_elements[cpt]).indices[i]-1;

            welt_oce_ice += std::hypot(M_VT[nind]-M_ocean[nind],M_VT[nind+M_num_nodes]-M_ocean[nind+M_num_nodes]);
            welt_air_ice += std::hypot(M_VT[nind]-M_wind [nind],M_VT[nind+M_num_nodes]-M_wind [nind+M_num_nodes]);
            welt_ice += std::hypot(M_VT[nind],M_VT[nind+M_num_nodes]);

            welt_ssh += M_ssh[nind];
        }


        double norm_Voce_ice = welt_oce_ice/3.;
        double norm_Vair_ice = welt_air_ice/3.;
        double norm_Vice = welt_ice/3.;

        double element_ssh = welt_ssh/3.;

        double coef_Vair = (vm["simul.lin_drag_coef_air"].as<double>()+(quad_drag_coef_air*norm_Vair_ice));
        coef_Vair *= (vm["simul.rho_air"].as<double>());

        double coef_Voce = (vm["simul.lin_drag_coef_water"].as<double>()+(quad_drag_coef_water*norm_Voce_ice));
        coef_Voce *= (vm["simul.rho_water"].as<double>());


        double critical_h = M_conc[cpt]*(M_element_depth[cpt]+element_ssh)/(vm["simul.Lemieux_basal_k1"].as<double>());

        double _coef = std::max(0., M_thick[cpt]-critical_h);
        double coef_basal = quad_drag_coef_air*basal_k2/(basal_drag_coef_air*(norm_Vice+basal_u_0));
        coef_basal *= _coef*std::exp(-basal_Cb*(1.-M_conc[cpt]));

        //std::vector<double> sigma_P(3,0.); /* temporary variable for the resistance to the compression */
        //std::vector<double> B0Tj_sigma_h(2,0);
        double b0tj_sigma_hu = 0.;
        double b0tj_sigma_hv = 0.;
        double mloc = 0.;

        std::vector<double> data(36);
        std::vector<double> fvdata(6,0.);
        std::vector<int> rcindices(6);

        double duu, dvu, duv, dvv;
        int index_u, index_v;

        double tmp_thick=(0.05>M_thick[cpt]) ? 0.05 : M_thick[cpt];
        double tmp_conc=(0.01>M_conc[cpt]) ? 0.01 : M_conc[cpt];

#if 1
        //option 1 (original)
        double coef = young*(1.-M_damage[cpt])*tmp_thick*std::exp(ridging_exponent*(1.-tmp_conc));
#else
        //option 2 (we just change the value of the ridging exponent and we renamed it "damaging_exponent")
        double damaging_exponent = -80.;
        double coef = young*(1.-M_damage[cpt])*tmp_thick*std::exp(damaging_exponent*(1.-tmp_conc));
#endif
        //option 3: We change the formulation of f(A) and make it piecewise linear between limit_conc_fordamage and 1, and 0 otherwise
        //double factor = 0.;
        //double limit_conc_fordamage = 0.;
        //limit_conc_fordamage=0.95;
        //if(tmp_conc<limit_conc_fordamage)
        //{
        //factor=0.;
        //}
        //else
        //{
        //factor=(tmp_conc-limit_conc_fordamage)/(1.-limit_conc_fordamage);
        //}
        //double coef = young*(1.-M_damage[cpt])*tmp_thick*factor;

        double coef_P = 0.;
        if(M_divergence_rate[cpt] < 0.)
        {
            coef_P = compression_factor*std::pow(tmp_thick,exponent_compression_factor)*std::exp(ridging_exponent*(1.-tmp_conc));
            coef_P = coef_P/(std::abs(M_divergence_rate[cpt])+divergence_min);
        }

        double mass_e = rhoi*tmp_thick + rhos*M_snow_thick[cpt];
        mass_e = (tmp_conc > 0.) ? (mass_e/tmp_conc):0.;
        double surface_e = M_surface[cpt];

        // /* compute the x and y derivative of g*ssh */
        double g_ssh_e_x = 0.;
        double g_ssh_e_y = 0.;
        double g_ssh_e;
        for(int i=0; i<3; i++)
        {
            g_ssh_e = (vm["simul.gravity"].as<double>())*M_ssh[(M_elements[cpt]).indices[i]-1] /*g_ssh*/;   /* g*ssh at the node k of the element e */
            g_ssh_e_x += M_shape_coeff[cpt][i]*g_ssh_e; /* x derivative of g*ssh */
            g_ssh_e_y += M_shape_coeff[cpt][i+3]*g_ssh_e; /* y derivative of g*ssh */
        }

        double coef_C     = mass_e*M_fcor[cpt];              /* for the Coriolis term */
        double coef_V     = mass_e/time_step;             /* for the inertial term */
        double coef_X     = - mass_e*g_ssh_e_x;              /* for the ocean slope */
        double coef_Y     = - mass_e*g_ssh_e_y;              /* for the ocean slope */

        double Vcor_index_v, Vcor_index_u;

        if (cpt < 0)
        {
            std::cout<<"************************\n";
            std::cout<<"Coef_C    = "<< coef_C <<"\n";
            std::cout<<"Coef_V    = "<< coef_V <<"\n";
            std::cout<<"Coef_X    = "<< coef_X <<"\n";
            std::cout<<"Coef_Y    = "<< coef_Y <<"\n";
            std::cout<<"coef_Vair = "<< coef_Vair <<"\n";
            std::cout<<"coef_Voce = "<< coef_Voce <<"\n";
            std::cout<<"coef_basal= "<< coef_basal <<"\n";
        }

        /* Loop over the 6 by 6 components of the finite element integral
         * this is done smartly by looping over j=0:2 and i=0:2
         * col = (mwIndex)it[2*j]-1  , row = (mwIndex)it[2*i]-1;
         * col  , row   -> UU component
         * col  , row+1 -> VU component
         * col+1, row   -> VV component
         * col+1, row+1 -> UV component */

        for(int j=0; j<3; j++)
        {
            /* Column corresponding to indice j (we also assemble terms in col+1) */
            //col = (mwIndex)it[2*j]-1; /* -1 to use the indice convention of C */

            index_u = (M_elements[cpt]).indices[j]-1;
            index_v = (M_elements[cpt]).indices[j]-1+M_num_nodes;

            Vcor_index_v=beta0*M_VT[index_v] + beta1*M_VTM[index_v] + beta2*M_VTMM[index_v];
            Vcor_index_u=beta0*M_VT[index_u] + beta1*M_VTM[index_u] + beta2*M_VTMM[index_u];


            for(int i=0; i<3; i++)
            {
                /* Row corresponding to indice i (we also assemble terms in row+1) */
                //row = (mwIndex)it[2*i]-1; /* -1 to use the indice convention of C */

                /* Select the nodal weight values from M_loc */
                mloc = M_Mass[3*j+i];

                b0tj_sigma_hu = 0.;
                b0tj_sigma_hv = 0.;

                for(int k=0; k<3; k++)
                {
                    b0tj_sigma_hu += M_B0T[cpt][k*6+2*i]*(M_sigma[3*cpt+k]*tmp_thick/*+sigma_P[k]*/);
                    b0tj_sigma_hv += M_B0T[cpt][k*6+2*i+1]*(M_sigma[3*cpt+k]*tmp_thick/*+sigma_P[k]*/);
                }

                /* ---------- UU component */
                duu = surface_e*( mloc*(coef_Vair+coef_Voce*cos_ocean_turning_angle+coef_V+coef_basal)
                                  +M_B0T_Dunit_B0T[cpt][(2*i)*6+2*j]*coef*time_step+M_B0T_Dunit_comp_B0T[cpt][(2*i)*6+2*j]*coef_P);

                /* ---------- VU component */
                dvu = surface_e*(+M_B0T_Dunit_B0T[cpt][(2*i+1)*6+2*j]*coef*time_step+M_B0T_Dunit_comp_B0T[cpt][(2*i+1)*6+2*j]*coef_P);

                /* ---------- UV component */
                duv = surface_e*(+M_B0T_Dunit_B0T[cpt][(2*i)*6+2*j+1]*coef*time_step+M_B0T_Dunit_comp_B0T[cpt][(2*i)*6+2*j+1]*coef_P);

                /* ---------- VV component */
                dvv = surface_e*( mloc*(coef_Vair+coef_Voce*cos_ocean_turning_angle+coef_V+coef_basal)
                                  +M_B0T_Dunit_B0T[cpt][(2*i+1)*6+2*j+1]*coef*time_step+M_B0T_Dunit_comp_B0T[cpt][(2*i+1)*6+2*j+1]*coef_P);

                data[(2*i  )*6+2*j  ] = duu;
                data[(2*i+1)*6+2*j  ] = dvu;
                data[(2*i  )*6+2*j+1] = duv;
                data[(2*i+1)*6+2*j+1] = dvv;

                fvdata[2*i] += surface_e*( mloc*( +M_tau[index_u]
                                                  +coef_Vair*M_wind[index_u]
                                                  +coef_Voce*cos_ocean_turning_angle*M_ocean[index_u]
                                                  +coef_X
                                                  +coef_V*M_VT[index_u]
                                                  -coef_Voce*sin_ocean_turning_angle*(M_ocean[index_v]-M_VT[index_v])
                                                  +coef_C*Vcor_index_v)
                                           - b0tj_sigma_hu/3);

                //std::cout<<"fvdata, M_tau, M_wind (u): "<<fvdata[2*i]<<","
                //    <<M_tau[index_u]<<","<<coef_Vair*M_wind[index_u]<<"\n";


                fvdata[2*i+1] += surface_e*( mloc*( +M_tau[index_v]
                                                    +coef_Vair*M_wind[index_v]
                                                    +coef_Voce*cos_ocean_turning_angle*M_ocean[index_v]
                                                    +coef_Y
                                                    +coef_V*M_VT[index_v]
                                                    +coef_Voce*sin_ocean_turning_angle*(M_ocean[index_u]-M_VT[index_u])
                                                    -coef_C*Vcor_index_u)
                                             - b0tj_sigma_hv/3);
            }

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
// #pragma omp critical(updatematrix)
//         {
//             M_matrix->addMatrix(&rcindices[0], rcindices.size(),
//                                 &rcindices[0], rcindices.size(), &data[0]);
//         }

        for (int idf=0; idf<rcindices.size(); ++idf)
        {
#pragma omp atomic
            rhsdata[rcindices[idf]] += fvdata[idf];

#if 1
            int indexr = rcindices[idf];
            for (int idj=0; idj<rcindices.size(); ++idj)
            {
                int indexc = rcindices[idj];
                int rnnz = M_graph.nNz()[indexr];
                int start = M_graph.ia()[indexr];
                int colind = 0;

                for (int io=start; io<start+rnnz; ++io)
                {
                    if (M_graph.ja()[io] == indexc)
                    {
                        colind = io-start;
                        break;
                    }
                }

#pragma omp atomic
                lhsdata[start+colind] += data[6*idf+idj];
            }
#endif
        }

        for (int const& idn : rcindices)
            M_node_max_conc[idn] = ((M_conc[cpt]>M_node_max_conc[idn])?(M_conc[cpt] ):(M_node_max_conc[idn])) ;

#if 0
#pragma omp critical(updatematrix)
        {
            M_matrix->addMatrix(&rcindices[0], rcindices.size(),
                                &rcindices[0], rcindices.size(), &data[0]);
        }
//#pragma omp end critical
#endif
    }

    // update petsc matrix
    boost::mpi::timer petsc_chrono;
    petsc_chrono.restart();
#pragma omp parallel for num_threads(max_threads)
    for (int cptpm=0; cptpm<2*M_num_nodes; ++cptpm)
    {
        int rnnz = M_graph.nNz()[cptpm];
        std::vector<int> lrcindices(rnnz);
        std::vector<double> ldata(rnnz);
        int start = M_graph.ia()[cptpm];
        //std::cout<<"Looking for "<< indexc << " in array of size "<< rnnz << " started by "<< start <<"\n";
        for (int io=start; io<start+rnnz; ++io)
        {
            lrcindices[io-start] = M_graph.ja()[io];
            ldata[io-start] = lhsdata[io];
        }
#if 0
        if (cptpm < 10)
        {
            std::cout<<"**************Row "<< cptpm <<"*************\n";
            for (int j=0; j<rnnz; ++j)
            {
                std::cout<< std::left << std::setw(12) << ldata[j] <<" for indice "<< lrcindices[j] <<"\n";
            }
        }
#endif
        M_matrix->setMatrix(&cptpm, 1,
                            &lrcindices[0], lrcindices.size(),
                            &ldata[0]);
    }
    LOG(DEBUG) <<"SET PETSC MATRIX done in " << petsc_chrono.elapsed() <<"s\n";
    lhsdata.resize(0);

    // close petsc matrix
    M_matrix->close();

    // update petsc vector and close it
    M_vector->addVector(&rhsindices[0], rhsindices.size(), &rhsdata[0]);
    M_vector->close();

    //LOG(DEBUG) <<"Assembling done in "<< chrono.elapsed() <<"s\n";
    LOG(INFO) <<"TIMER ASSEMBLY= " << chrono.elapsed() <<"s\n";

    // extended dirichlet nodes (add nodes where M_conc <= 0)
    for (int i=0; i<2*M_num_nodes; ++i)
    {
        if(M_node_max_conc[i]<=0.)
        {
            extended_dirichlet_nodes.push_back(i);
        }
    }
    std::sort(extended_dirichlet_nodes.begin(), extended_dirichlet_nodes.end());
    extended_dirichlet_nodes.erase(std::unique( extended_dirichlet_nodes.begin(), extended_dirichlet_nodes.end() ), extended_dirichlet_nodes.end());

    chrono.restart();
    //M_matrix->on(M_dirichlet_nodes,*M_vector);
    M_matrix->on(extended_dirichlet_nodes,*M_vector);
    LOG(DEBUG) <<"TIMER DBCA= " << chrono.elapsed() <<"s\n";

    if (M_log_level == DEBUG)
    {
        LOG(DEBUG) <<"[PETSC MATRIX] CLOSED      = "<< M_matrix->closed() <<"\n";
        LOG(DEBUG) <<"[PETSC MATRIX] SIZE        = "<< M_matrix->size1() << " " << M_matrix->size2() <<"\n";
        LOG(DEBUG) <<"[PETSC MATRIX] SYMMETRIC   = "<< M_matrix->isSymmetric() <<"\n";
        LOG(DEBUG) <<"[PETSC MATRIX] NORM        = "<< M_matrix->linftyNorm() <<"\n";
        LOG(DEBUG) <<"[PETSC VECTOR] NORM        = "<< M_vector->l2Norm() <<"\n";
    }

    //M_matrix->printMatlab("stiffness.m");
    //M_vector->printMatlab("rhs.m");
}

void
FiniteElement::nodeMaxConc()
{
    int thread_id;
    int total_threads;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    //std::cout<<"MAX THREADS= "<< max_threads <<"\n";

	M_node_max_conc.assign(2*M_num_nodes,0.);

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        std::vector<int> rcindices(6);

        int index_u, index_v;

        for(int j=0; j<3; j++)
        {
            /* Column corresponding to indice j (we also assemble terms in col+1) */
            //col = (mwIndex)it[2*j]-1; /* -1 to use the indice convention of C */

            index_u = (M_elements[cpt]).indices[j]-1;
            index_v = (M_elements[cpt]).indices[j]-1+M_num_nodes;

            rcindices[2*j] = index_u;
            rcindices[2*j+1] = index_v;
        }

        for (int const& idn : rcindices)
            M_node_max_conc[idn] = ((M_conc[cpt]>M_node_max_conc[idn])?(M_conc[cpt] ):(M_node_max_conc[idn])) ;
    }
}

void
FiniteElement::tensors()
{
    M_Dunit.assign(9,0);
    M_Dunit_comp.assign(9,0);
    M_Mass.assign(9,0);

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

#if 1

    M_B0T.resize(M_num_elements);
    M_B0T_Dunit_B0T.resize(M_num_elements);
    M_B0T_Dunit_comp_B0T.resize(M_num_elements);
    M_shape_coeff.resize(M_num_elements);

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

        M_shape_coeff[cpt]        = shapecoeff;
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
    std::vector<double> UM_P = M_UM;

    for (int nd=0; nd<M_UM.size(); ++nd)
    {
        M_UM[nd] += time_step*M_VT[nd];
    }

    for (const int& nd : M_neumann_nodes)
    {
        M_UM[nd] = UM_P[nd];
    }

    LOG(DEBUG) <<"VT MIN= "<< *std::min_element(M_VT.begin(),M_VT.end()) <<"\n";
    LOG(DEBUG) <<"VT MAX= "<< *std::max_element(M_VT.begin(),M_VT.end()) <<"\n";

    int thread_id;
    int total_threads;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    //std::cout<<"MAX THREADS= "<< max_threads <<"\n";
#if defined (WAVES)
    if (M_use_wim)
        M_dfloe.assign(M_num_elements,0.);
#endif

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        double old_thick;
        double old_snow_thick;
        double old_conc;
        double old_damage;
        double old_h_ridged_thick_ice;
#if defined (WAVES)
        double old_nfloes;
#endif

        double old_h_thin;
        double old_hs_thin;
        double old_h_ridged_thin_ice;

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
        double sigma_s, sigma_n, sigma_1, sigma_2;
        double tract_max, sigma_t, sigma_c, q;
        double tmp, sigma_target;

        /* some variables used for the ice redistribution*/
        double tanalpha, rtanalpha, del_v, del_c, del_vs, new_v_thin;

        /* Indices loop*/
        int i,j;

        bool to_be_updated;

        /* set constants for the ice redistribution */
        tanalpha  = h_thin_max/c_thin_max;
        rtanalpha = 1./tanalpha;

        // beginning of the original code (without openMP)
        // Temporary memory
        old_thick = M_thick[cpt];
        old_snow_thick = M_snow_thick[cpt];
        old_conc = M_conc[cpt];
        old_damage = M_damage[cpt];
        old_h_ridged_thick_ice=M_h_ridged_thick_ice[cpt];
#if defined (WAVES)
        if (M_use_wim)
            old_nfloes = M_nfloes[cpt];
#endif

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            old_h_thin = M_h_thin[cpt];
            old_hs_thin=M_hs_thin[cpt];
            old_h_ridged_thin_ice=M_h_ridged_thin_ice[cpt];
        }

        /*======================================================================
         * Diagnostic:
         * Elastic deformation and instantaneous deformation rate
         *======================================================================
         */

        /* Compute the elastic deformation and the instantaneous deformation rate */
        for(i=0;i<3;i++)
        {
            epsilon_veloc_i = 0.0;
            for(j=0;j<3;j++)
            {
                /* deformation */
                //col = (mwIndex)it[j]-1;
                epsilon_veloc_i += M_B0T[cpt][i*6 + 2*j]*M_VT[(M_elements[cpt]).indices[j]-1]  ;
                epsilon_veloc_i += M_B0T[cpt][i*6 + 2*j + 1]*M_VT[(M_elements[cpt]).indices[j]-1+M_num_nodes]  ;
            }

            epsilon_veloc[i] = epsilon_veloc_i;
        }

        M_divergence_rate[cpt]= (epsilon_veloc[0]+epsilon_veloc[1]);

        /*======================================================================
         * Update the internal stress
         *======================================================================
         */

#if 0
        // To be uncommented if we use option 3:
        double factor = 0.;
        double limit_conc_fordamage = 0.95;

        if(limit_conc_fordamage <= old_conc)
        {
            factor = (old_conc-limit_conc_fordamage)/(1.-limit_conc_fordamage);
        }
#endif

#if 1
        //option 1 (original)
        double damaging_exponent = ridging_exponent;
#else
        //option 2
        double damaging_exponent = -80.;
#endif
        for(i=0;i<3;i++)
        {
            sigma_dot_i = 0.0;
            for(j=0;j<3;j++)
            {
            // sigma_dot_i += std::exp(ridging_exponent*(1.-old_conc))*young*(1.-old_damage)*M_Dunit[i*3 + j]*epsilon_veloc[j];
            sigma_dot_i += std::exp(damaging_exponent*(1.-old_conc))*young*(1.-old_damage)*M_Dunit[i*3 + j]*epsilon_veloc[j];
            //sigma_dot_i += factor*young*(1.-old_damage)*M_Dunit[i*3 + j]*epsilon_veloc[j];
            }

            M_sigma[3*cpt+i] += time_step*sigma_dot_i;
            sigma_pred[i]    = M_sigma[3*cpt+i] + time_step*sigma_dot_i;
        }

        /*======================================================================
         * Correct the internal stress and the damage
         *======================================================================
         */

        /* Compute the shear and normal stress, which are two invariants of the internal stress tensor */

        sigma_s=std::hypot((sigma_pred[0]-sigma_pred[1])/2.,sigma_pred[2]);
        sigma_n=-          (sigma_pred[0]+sigma_pred[1])/2.;

        sigma_1 = sigma_n+sigma_s; // max principal component following convention (positive sigma_n=pressure)
        sigma_2 = sigma_n-sigma_s; // max principal component following convention (positive sigma_n=pressure)

        q=std::pow(std::pow(std::pow(tan_phi,2.)+1,.5)+tan_phi,2.);
        sigma_c=2.*M_Cohesion[cpt]/(std::pow(std::pow(tan_phi,2.)+1,.5)-tan_phi);

        sigma_t=-sigma_c/q;

        /* minimum and maximum normal stress */
        tract_max=-tract_coef*M_Cohesion[cpt]/tan_phi;

        /* Correction of the damage */
        if(sigma_n>M_Compressive_strength[cpt])
        {
            sigma_target=M_Compressive_strength[cpt];

            tmp=1.0-sigma_target/sigma_n*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }
#if 0
        if(sigma_1<0 && sigma_2<sigma_t)
        {
            sigma_target=sigma_t;

            tmp=1.0-sigma_target/sigma_2*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }
#endif
#if 1
        if(sigma_1-q*sigma_2>sigma_c)
        {
            sigma_target=sigma_c;

            tmp=1.0-sigma_target/(sigma_1-q*sigma_2)*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }
#endif

#if 1
        if(sigma_n<tract_max)
        {
            sigma_target=tract_max;

            tmp=1.0-sigma_target/sigma_n*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }
#endif

#if 0
        if(sigma_s>M_Cohesion[cpt]-sigma_n*tan_phi)
        {
            tmp=1.0-M_Cohesion[cpt]/(sigma_s+sigma_n*tan_phi)*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }
#endif

        /*
         * Diagnostic:
         * Recompute the internal stress
         */
        for(i=0;i<3;i++)
        {
#if 1
            if(old_damage<1.0)
            {
                M_sigma[3*cpt+i] = (1.-M_damage[cpt])/(1.-old_damage)*M_sigma[3*cpt+i] ;
            }
            else
            {
                M_sigma[3*cpt+i] = 0. ;
            }
#endif
#if 0
            // test to boost the localization
            if(M_damage[cpt]!=old_damage)
            {
                M_damage[cpt]=1.;
                M_sigma[3*cpt+i] = 0. ;
            }
#endif
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
        tmp=1./(1.-M_damage[cpt]);
        tmp-= 1000*time_step/M_time_relaxation_damage[cpt];
        tmp=((tmp>1.)?(tmp):(1.));
        M_damage[cpt]=-1./tmp + 1.;

        /*======================================================================
         * Update:
         * Ice and snow thickness, and concentration using a Lagrangian or an Eulerian scheme
         *======================================================================
         */

        to_be_updated=false;
        if( M_divergence_rate[cpt]!=0.)
            to_be_updated=true;

#if 0
        /* For the Lagrangian scheme, we do not update the variables for the elements having one node on the open boundary. */
        if(std::find(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[0]-1) != M_neumann_flags.end() ||
           std::find(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[1]-1) != M_neumann_flags.end() ||
           std::find(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[2]-1) != M_neumann_flags.end())
            to_be_updated=false;
#endif
        if(std::binary_search(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[0]-1) ||
           std::binary_search(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[1]-1) ||
           std::binary_search(M_neumann_flags.begin(),M_neumann_flags.end(),(M_elements[cpt]).indices[2]-1))
            to_be_updated=false;

        if((old_conc>0.)  && to_be_updated)
        {
            surface = this->measure(M_elements[cpt],M_mesh, UM_P);
            surface_new = this->measure(M_elements[cpt],M_mesh,M_UM);

            ice_surface = old_conc*surface;
            ice_volume = old_thick*surface;
            snow_volume = old_snow_thick*surface;
            ridged_thick_ice_volume = old_h_ridged_thick_ice*surface;
#if defined (WAVES)
            if (M_use_wim)
            {
                M_nfloes[cpt]  = old_nfloes*surface/surface_new;

                // lower bounds
                M_nfloes[cpt] = ((M_nfloes[cpt]>0.)?(M_nfloes[cpt] ):(0.)) ;
            }
#endif

            M_conc[cpt]    = ice_surface/surface_new;
            M_thick[cpt]   = ice_volume/surface_new;
            M_snow_thick[cpt]   = snow_volume/surface_new;
            M_h_ridged_thick_ice[cpt]   =   ridged_thick_ice_volume/surface_new;

            if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                thin_ice_volume = old_h_thin*surface;
                thin_snow_volume = old_hs_thin*surface;
                ridged_thin_ice_volume = old_h_ridged_thin_ice*surface;

                M_h_thin[cpt]        = thin_ice_volume/surface_new;
                M_hs_thin[cpt]   = thin_snow_volume/surface_new;
                M_h_ridged_thin_ice[cpt]    =   ridged_thin_ice_volume/surface_new;
            }

            /* Ridging scheme */
            if(surface_new<surface)
            {
                if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
                {
                    ridging_thin_ice=(surface-surface_new)/surface_new*old_h_thin;
                    ridging_snow_thin_ice=(surface-surface_new)/surface_new*old_hs_thin;

                    ridging_thin_ice=((ridging_thin_ice<old_h_thin)?(ridging_thin_ice):(old_h_thin));
                    ridging_snow_thin_ice=((ridging_snow_thin_ice<old_hs_thin)?(ridging_snow_thin_ice):(old_hs_thin)) ;

                    M_thick[cpt] += ridging_thin_ice ;
                    M_h_thin[cpt] -= ridging_thin_ice ;

                    M_snow_thick[cpt] += ridging_snow_thin_ice ;
                    M_hs_thin[cpt] -= ridging_snow_thin_ice ;

                    M_h_ridged_thin_ice[cpt] += ridging_thin_ice;
                    M_conc[cpt] += ridging_thin_ice/ridge_h;

                    /* upper bounds (only for the concentration) */
                    ridging_thin_ice = ((M_conc[cpt]<1.)?(0.):(M_h_thin[cpt])) ;
                    ridging_snow_thin_ice = ((M_conc[cpt]<1.)?(0.):(M_hs_thin[cpt])) ;

                    M_snow_thick[cpt] += ridging_snow_thin_ice ;
                    M_hs_thin[cpt] -= ridging_snow_thin_ice ;

                    M_thick[cpt] += ridging_thin_ice;
                    M_h_thin[cpt] -= ridging_thin_ice;
                }
                /* upper bounds (only for the concentration) */
                ridging_thick_ice=((M_conc[cpt]<1.)?(0.):(M_thick[cpt]*(M_conc[cpt]-1.)));
                M_conc[cpt] = ((M_conc[cpt]<1.)?(M_conc[cpt]):(1.)) ;
            }

            /* lower bounds */
            M_conc[cpt] = ((M_conc[cpt]>0.)?(M_conc[cpt] ):(0.)) ;
            M_thick[cpt]        = ((M_thick[cpt]>0.)?(M_thick[cpt]     ):(0.)) ;
            M_snow_thick[cpt]   = ((M_snow_thick[cpt]>0.)?(M_snow_thick[cpt]):(0.)) ;

            if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
            {
                M_h_thin[cpt]    = ((M_h_thin[cpt]>0.)?(M_h_thin[cpt] ):(0.)) ;
                M_hs_thin[cpt]   = ((M_hs_thin[cpt]>0.)?(M_hs_thin[cpt]):(0.)) ;
            }
        }

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            /* Compute the redistribution of thin ice. */
            /* Returns the change in volume and concentration of thick ice as well as the
             * change in volume of thin ice. It is called after the
             * dynamics are done. */

            if(M_h_thin[cpt]>0.)
            {
                thin_ice_redistribute(M_h_thin[cpt], M_hs_thin[cpt], 0., M_conc[cpt],
                                      tanalpha, rtanalpha, h_thin_max, &new_v_thin, &del_v, &del_c, &del_vs);

                M_conc[cpt]       += del_c;

                M_thick[cpt]           += del_v;
                M_h_thin[cpt]      -= del_v;

                M_snow_thick[cpt]      += del_vs;
                M_hs_thin[cpt] -= del_vs;
            }
            else
            {
                M_snow_thick[cpt] += M_hs_thin[cpt];
                M_hs_thin[cpt] = 0. ;
            }
        }

#if defined (WAVES)
        if (M_use_wim)
        {
            //update Dfloe
            if (M_nfloes[cpt] > 0)
                M_dfloe[cpt] = std::sqrt(M_conc[cpt]/M_nfloes[cpt]);

            if (M_dfloe[cpt] > vm["wim.dfloepackthresh"].template as<double>())
                M_dfloe[cpt] = vm["wim.dfloepackinit"].template as<double>();

            if (M_conc[cpt] < vm["wim.cicemin"].template as<double>())
            {
                M_nfloes[cpt] = 0.;
                M_dfloe[cpt] = 0.;
            }
        }
#endif
    }
}

void
FiniteElement::solve()
{
    chrono.restart();
    M_solver->solve(_matrix=M_matrix,
                    _solution=M_solution,
                    _rhs=M_vector,
                    _ksp=vm["solver.ksp-type"].as<std::string>()/*"preonly"*/,
                    _pc=vm["solver.pc-type"].as<std::string>()/*"cholesky"*/,
                    _pcfactormatsolverpackage=vm["solver.mat-package-type"].as<std::string>()/*"cholmod"*/,
                    _reuse_prec=true,
                    _rebuild=M_regrid
                    );

    LOG(INFO) <<"TIMER SOLUTION= " << chrono.elapsed() <<"s\n";

    //M_solution->printMatlab("solution.m");
    //Environment::logMemoryUsage("");
}

// Routine for the 1D thermodynamical model
// No thin ice for now
// No stability dependent drag for now
void
FiniteElement::thermo()
{
    // There is now only one big loop for the thermodynamics so that we can use multithreading.

    // constant variables
    // Set local variable to values defined by options
    double const timeT = vm["simul.ocean_nudge_timeT"].as<double>();
    double const timeS = vm["simul.ocean_nudge_timeS"].as<double>();
    double const Qdw_const = vm["simul.constant_Qdw"].as<double>();
    double const Fdw_const = vm["simul.constant_Fdw"].as<double>();

    double const ocean_albedo=vm["simul.albedoW"].as<double>();
    double const drag_ocean_t=vm["simul.drag_ocean_t"].as<double>();
    double const drag_ocean_q=vm["simul.drag_ocean_q"].as<double>();

    double const rh0   = 1./vm["simul.hnull"].as<double>();
    double const rPhiF = 1./vm["simul.PhiF"].as<double>();

    double const tanalpha  = h_thin_max/c_thin_max;
    double const rtanalpha = 1./tanalpha;

    double const qi = physical::Lf * physical::rhoi;
    double const qs = physical::Lf * physical::rhos;

    int const newice_type = vm["simul.newice_type"].as<int>();
    int const melt_type = vm["simul.melt_type"].as<int>();
    double const PhiM = vm["simul.PhiM"].as<double>();
    double const PhiF = vm["simul.PhiF"].as<double>();

    const double aw=6.1121e2, bw=18.729, cw=257.87, dw=227.3;
    const double Aw=7.2e-4, Bw=3.20e-6, Cw=5.9e-10;

    const double alpha=0.62197, beta=0.37803;

    // initialisation of the multithreading
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int i=0; i < M_num_elements; ++i)
    {
        // -------------------------------------------------
        // 1) Initialization of the temporary variables

        double  hi=0.;     // Ice thickness (slab)
        double  hi_old=0.; // Ice thickness at the start of the time step (slab)
        double  hs=0.;     // Snow thickness (slab)

        double  hi_thin=0.;     // Thin ice thickness (slab)
        double  hi_thin_old=0.; // Thin ice thickness at the start of the time step (slab)
        double  hs_thin=0.;     // Snow thickness on thin ice (slab)

        double  del_hi=0.; // Change in ice thickness (slab only)
        double  del_hi_thin=0.; // Change in thin ice thickness (slab only)

        double  evap=0.;   // Evaporation

        double  Qdw=0.;    // Heat flux from ocean nudging
        double  Fdw=0.;    // Fresh water flux from ocean nudging

        double  Qio=0.;    // Ice-ocean heat flux
        double  Qio_thin=0.;    // Ice-ocean heat flux through thin ice
        double  Qai=0.;    // Atmosphere-ice heat flux
        double  Qow=0.;    // Open water heat flux

        // Save old _volumes_ and concentration and calculate wind speed
        double  old_vol=M_thick[i];
        double  old_snow_vol=M_snow_thick[i];
        double  old_conc=M_conc[i];

        double  old_h_thin = 0;
        double  old_hs_thin = 0;
        double  old_conc_thin=0;
        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            old_h_thin  = M_h_thin[i];
            old_hs_thin = M_hs_thin[i];
            old_conc_thin = std::min(std::min(old_h_thin/physical::hmin, std::sqrt(2.*old_h_thin*rtanalpha)), 1-old_conc);
        }

        double sum_u=0.;
        double sum_v=0.;

        for (int j=0; j<3; ++j)
        {
            // calculate wind per element
            sum_u += M_wind[M_elements[i].indices[j]-1];
            sum_v += M_wind[M_elements[i].indices[j]-1+M_num_nodes];
        }
        double  wspeed = std::hypot(sum_u, sum_v)/3.;

        // definition of the fraction of snow
        double tmp_snowfr;
        if(M_snowfr.M_initialized)
            tmp_snowfr=M_snowfr[i];
        else
        {
            if(M_tair[i]<0)
                tmp_snowfr=1.;
            else
                tmp_snowfr=0.;
        }

        double tmp_Qsw_in;
        if(M_Qsw_in.M_initialized)
            tmp_Qsw_in=M_Qsw_in[i];
        else
        {
            throw std::logic_error("The function approxSW not yet implemented, you need to initialized M_Qsw_in");
            //tmp_Qsw_in=approxSW();
        }

        double tmp_mld=( M_mld[i] > vm["simul.constant_mld"].as<double>() ) ? M_mld[i] : vm["simul.constant_mld"].as<double>();

        // -------------------------------------------------
        // 2) We calculate or set the flux due to nudging
        if ( M_ocean_type == setup::OceanType::CONSTANT )
        {
            Qdw=Qdw_const;
            Fdw=Fdw_const;
        }
        else
        {
            // nudgeFlux
            if ( M_ocean_salt[i] > physical::si )
            {
                Qdw = -(M_sst[i]-M_ocean_temp[i]) * tmp_mld * physical::rhow * physical::cpw/timeT;

                double delS = M_sss[i] - M_ocean_salt[i];
                Fdw = delS * tmp_mld * physical::rhow /(timeS*M_sss[i] - time_step*delS);
            } else {
                Qdw = Qdw_const;
                Fdw = Fdw_const;
            }
        }

        // SYL: Calculate the drag coefficients missing??

        // -------------------------------------------------
        // 3) Calculate fluxes in the open water portion (openWaterFlux in old c code, Qow_mex in matlab)
        double sphuma, sphumw;

        // Calculate atmospheric fluxes

        /* Out-going long-wave flux */
        double Qlw_out = physical::eps*physical::sigma_sb*std::pow(M_sst[i]+physical::tfrwK,4.);

        // -------------------------------------------------
        // 3.1) Specific humidity - atmosphere (calcSphumA in matlab)
        /* There are two ways to calculate this. We decide which one by
         * checking mixrat - the calling routine must set this to a negative
         * value if the dewpoint should be used. */
        if ( M_mixrat[i] < 0. )
        {
            double fa     = 1. + Aw + M_mslp[i]*1e-2*( Bw + Cw*M_dair[i]*M_dair[i] );
            double esta   = fa*aw*std::exp( (bw-M_dair[i]/dw)*M_dair[i]/(M_dair[i]+cw) );
            sphuma = alpha*fa*esta/(M_mslp[i]-beta*fa*esta) ;
        }
        else
            sphuma = M_mixrat[i]/(1.+M_mixrat[i]) ;

        // -------------------------------------------------
        // 3.2) Specific humidity - ocean surface (calcSphumW in matlab)
        double fw     = 1. + Aw + M_mslp[i]*1e-2*( Bw + Cw*M_sst[i]*M_sst[i] );
        double estw   = aw*std::exp( (bw-M_sst[i]/dw)*M_sst[i]/(M_sst[i]+cw) )*(1-5.37e-4*M_sss[i]);
        sphumw = alpha*fw*estw/(M_mslp[i]-beta*fw*estw) ;

        // -------------------------------------------------
        /* Density of air */
        double rhoair = M_mslp[i]/(physical::Ra*(M_tair[i]+tfrwK)) * (1.+sphuma)/(1.+1.609*sphuma);

        /* Sensible heat flux */
        double Qsh = drag_ocean_t*rhoair*physical::cpa*wspeed*( M_sst[i] - M_tair[i] );

        /* Latent heat flux */
        double Lv  = physical::Lv0 - 2.36418e3*M_tair[i] + 1.58927*M_tair[i]*M_tair[i] - 6.14342e-2*std::pow(M_tair[i],3.);
        double Qlh = drag_ocean_q*rhoair*Lv*wspeed*( sphumw - sphuma );

        /* Evaporation */
        evap = Qlh/(physical::rhofw*Lv);

        // Sum them up
        double tmp_Qlw_in;
        if(M_Qlw_in.M_initialized)
            tmp_Qlw_in=M_Qlw_in[i];
        else
        {
            double tsa = M_tice[0][i] + tfrwK;
            double taa = M_tair[i]  + tfrwK;
        	tmp_Qlw_in = sigma_sb*pow(taa,4) \
        			*( 1. - 0.261*exp(-7.77e-4*std::pow(taa-tfrwK,2)) ) \
        			*( 1. + 0.275*M_tcc[i] );
        }

        Qow = -tmp_Qsw_in*(1.-ocean_albedo) - tmp_Qlw_in + Qlw_out + Qsh + Qlh;

        // -------------------------------------------------
        // 4) Thickness change of the ice slab (thermoIce0 in matlab)

        switch ( M_thermo_type )
        {
            case setup::ThermoType::ZERO_LAYER:
                this->thermoIce0(i, wspeed, sphuma, M_conc[i], M_thick[i], M_snow_thick[i], tmp_Qlw_in, tmp_Qsw_in, tmp_mld, tmp_snowfr, hi, hs, hi_old, Qio, del_hi, M_tice[0][i]);
                break;
            case setup::ThermoType::WINTON:
                this->thermoWinton(i, time_step, wspeed, sphuma, M_conc[i], M_thick[i], M_snow_thick[i], tmp_Qlw_in, tmp_Qsw_in, tmp_mld, tmp_snowfr, hi, hs, hi_old, Qio, del_hi,
                        M_tice[0][i], M_tice[1][i], M_tice[2][i]);
                break;
        }

        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            this->thermoIce0(i, wspeed, sphuma, old_conc_thin, M_h_thin[i], M_hs_thin[i], tmp_Qlw_in, tmp_Qsw_in, tmp_mld, tmp_snowfr, hi_thin, hs_thin, hi_thin_old, Qio_thin, del_hi_thin, M_tsurf_thin[i]);
            M_h_thin[i]  = hi_thin * old_conc_thin;
            M_hs_thin[i] = hs_thin * old_conc_thin;
        }

        // -------------------------------------------------
        // 5) Ice growth over open water and lateral melt (thermoOW in matlab)

        /* Local variables */
        double tw_new, tfrw, newice, del_c, newsnow, h0;

        /* dT/dt due to heatflux atmos.-ocean */
        tw_new = M_sst[i] - Qow*time_step/(tmp_mld*physical::rhow*physical::cpw);
        tfrw   = -physical::mu*M_sss[i];

        /* Form new ice in case of super cooling, and reset Qow and evap */
        if ( tw_new < tfrw )
        {
            newice  = (1.-M_conc[i])*(tfrw-tw_new)*tmp_mld*physical::rhow*physical::cpw/qi;
            Qow  = -(tfrw-M_sst[i])*tmp_mld*physical::rhow*physical::cpw/time_step;
            // evap = 0.;
        } else {
            newice  = 0.;
        }

        /* Decide the change in ice fraction (del_c) */
        /* Initialise to be safe */
        del_c = 0.;
        newsnow = 0.;
        if ( newice > 0. )
        {
            /* Freezing conditions */
            switch ( newice_type )
            {
                case 1:
                    /* Hibler's (79) approach */
                    del_c = newice*rh0;
                    break;
                case 2:
                    /* Mellor and Kantha (89) */
                    if ( hi_old > 0. )
                    {
                            del_c = newice*PhiF/hi_old;
                    } else {
                            del_c = 1.;
                    }
                    break;
                case 3:
                    /* Olason and Harms (09) */
                    h0    = (1.+0.1*wspeed)/15.;
                    del_c = newice/std::max(rPhiF*hi_old,h0);
                    break;
                case 4:
                    /* Thin ice category */
                    thin_ice_redistribute(M_h_thin[i], M_hs_thin[i], newice/(1.-M_conc[i]), M_conc[i],
                                      tanalpha, rtanalpha, h_thin_max, &M_h_thin[i], &newice, &del_c, &newsnow);

                    // Change the snow _thickness_ for thick ice and _volume_ for thin ice
                    M_hs_thin[i]    -= newsnow;
                    // M_snow_thick[i] += newsnow; <- this is done properly below
                    break;
                default:
                    std::cout << "newice_type = " << newice_type << "\n";
                    throw std::logic_error("Wrong newice_type");
            }
            /* Check bounds on del_c */
            del_c = std::min( 1.-M_conc[i], del_c );
        }
        else if ( del_hi < 0. )
        {
            /* Melting conditions */
            switch ( melt_type )
            {
                case 1:
                    /* Hibler (79) using PhiM to tune. PhiM = 0.5 is
                     * equivilent Hibler's (79) approach */
                    if ( M_conc[i] < 1. )
                        del_c = del_hi*M_conc[i]*PhiM/hi_old;
                    else
                        del_c = 0.;
                    break;
                case 2:
                    /* Mellor and Kantha (89) */
                    /* Use the fraction PhiM of (1-c)*Qow to melt laterally */
                    del_c = PhiM*(1.-M_conc[i])*std::min(0.,Qow)*time_step/( hi*qi+hs*qs );
                    /* Deliver the fraction (1-PhiM) of Qow to the ocean */
                    Qow = (1.-PhiM)*Qow;
                    // This is handled below
                    // /* + Deliver excess energy to the ocean when there's no ice left */
                    //         + std::min(0., std::max(0.,M_conc[i]+del_c)*( hi*qi+hs*qs )/time_step);
                    // /* Don't suffer negative c! */
                    // del_c = std::max(del_c, -M_conc[i]);
                    break;
                default :
                    std::cout << "newice_type = " << newice_type << "\n";
                    throw std::logic_error("Wrong newice_type");
            }
        }
        else
        {
            /* There is no ice growing or melting */
            del_c = 0.;
        }

        /* New concentration */
        M_conc[i] = M_conc[i] + del_c;

        /* New thickness */
        /* We conserve volume and energy */
        if ( M_conc[i] >= physical::cmin )
        {
            hi = ( hi*old_conc + newice )/M_conc[i];
            if ( del_c < 0. )
            {
                /* We conserve the snow height, but melt away snow as the concentration decreases */
                Qow = Qow + del_c*hs*qs/time_step;
            } else {
                /* Snow volume is conserved as concentration increases */
                hs  = ( hs*old_conc + newsnow )/M_conc[i];
            }

            if ( M_thermo_type == setup::ThermoType::WINTON )
            {
                // Add newice evenly to both layers and recalculate temperature
                double f1    = M_thick[i]/(M_thick[i]+newice); // Fraction of old ice (as opposed to newice) in the upper layer
                double Tbar  = f1*( M_tice[1][i] - physical::Lf*physical::mu*physical::si/(physical::C*M_tice[1][i]) ) + (1-f1)*tfrw; // (39)
                M_tice[1][i] = ( Tbar - std::sqrt(Tbar*Tbar + 4*physical::mu*physical::si*physical::Lf/physical::C) )/2.; // (38)
                M_tice[2][i] = f1*M_tice[2][i] + (1-f1)*tfrw; // (26) slightly rewritten
            }
        }

        /* Check limits */
        if ( M_conc[i] < physical::cmin || hi < physical::hmin )
        {
            // Extract heat from the ocean corresponding to the heat in the
            // remaining ice and snow
            Qow    = Qow + M_conc[i]*hi*qi/time_step + M_conc[i]*hs*qs/time_step;
            M_conc[i]  = 0.;
            for (int j=0; j<M_tice.size(); j++)
                M_tice[j][i] = tfrw;
            M_tsurf_thin[i] = tfrw;
            hi     = 0.;
            hs     = 0.;
        }

        // -------------------------------------------------
        // 6) Calculate effective ice and snow thickness
        M_thick[i] = hi*M_conc[i];
        M_snow_thick[i] = hs*M_conc[i];

        // -------------------------------------------------
        // 7) Slab Ocean (slabOcean in matlab)

        // local variables
        double del_vi;      // Change in ice volume
        double del_vs;      // Change in snow olume
        double rain;        // Liquid precipitation
        double emp;         // Evaporation minus liquid precipitation
        double Qio_mean;    // Element mean ice-ocean heat flux
        double Qow_mean;    // Element mean open water heat flux

        // Calculate change in volume to calculate salt rejection
        del_vi = M_thick[i] - old_vol + M_h_thin[i] - old_h_thin;
        del_vs = M_snow_thick[i] - old_snow_vol + M_hs_thin[i] - old_hs_thin;

        // Rain falling on ice falls straight through. We need to calculate the
        // bulk freshwater input into the entire cell, i.e. everything in the
        // open-water part plus rain in the ice-covered part.
        rain = (1.-old_conc-old_conc_thin)*M_precip[i] + (old_conc+old_conc_thin)*(1.-tmp_snowfr)*M_precip[i];
        emp  = (evap*(1.-old_conc-old_conc_thin)-rain);

        Qio_mean = Qio*old_conc + Qio_thin*old_conc_thin;
        Qow_mean = Qow*(1.-old_conc-old_conc_thin);

        /* Heat-flux */
        M_sst[i] = M_sst[i] - time_step*( Qio_mean + Qow_mean - Qdw )/(physical::rhow*physical::cpw*tmp_mld);

        /* Change in salinity */
        M_sss[i] = M_sss[i] + ( (M_sss[i]-physical::si)*physical::rhoi*del_vi + M_sss[i]*(del_vs*physical::rhos + (emp-Fdw)*time_step) )
            / ( tmp_mld*physical::rhow - del_vi*physical::rhoi - ( del_vs*physical::rhos + (emp-Fdw)*time_step) );

        // -------------------------------------------------
        // 8) Damage manipulation (thermoDamage in matlab)

        // local variables
        double deltaT;      // Temperature difference between ice bottom and the snow-ice interface

        // Newly formed ice is undamaged - calculate damage as a weighted
        // average of the old damage and 0, weighted with volume.
        if ( M_thick[i] > old_vol )
            M_damage[i] = M_damage[i]*old_vol/M_thick[i];

        // SYL: manipulation of the ridged ice missing??

        // Set time_relaxation_damage to be inversely proportional to
        // temperature difference between bottom and snow-ice interface
        if ( M_thick[i] > 0. )
        {
            deltaT = std::max(0., -physical::mu*M_sss[i] - M_tice[0][i] )
                / ( 1. + physical::ki*M_snow_thick[i]/(physical::ks*M_thick[i]) );
            M_time_relaxation_damage[i] = std::max(time_relaxation_damage*deltaT_relaxation_damage/deltaT, time_step);
        } else {
            M_time_relaxation_damage[i] = time_relaxation_damage;
        }
        // -------------------------------------------------

    }// end for loop
}// end thermo function

// Atmospheric fluxes through bulk formula
void
FiniteElement::atmFluxBulk(int i, double Tsurf, double sphuma, double drag_ice_t, double Qsw, double Qlw_in, double wspeed,
        double &Qai, double &dQaidT, double &subl)
{
    // Constants
    const double ai=6.1115e2, bi=23.036, ci=279.82, di=333.7;
    const double Ai=2.2e-4, Bi=3.83e-6, Ci=6.4e-10;

    const double alpha=0.62197, beta=0.37803;

    // -------------------------------------------------
    // 4.1) BULK ICE
    /* Out-going long-wave flux and derivative */
    double Qlw_out =   physical::eps * physical::sigma_sb * std::pow(Tsurf+physical::tfrwK,4);
    double dQlwdT  = 4.*physical::eps * physical::sigma_sb * std::pow(Tsurf+physical::tfrwK,3);

    // -------------------------------------------------
    // calcSphumI
    /* Specific humidity - ice surface */
    double fi      = 1. + Ai + M_mslp[i]*1e-2*( Bi + Ci*Tsurf*Tsurf );
    double esti    = ai*std::exp( (bi-Tsurf/di)*Tsurf/(Tsurf+ci) );
    double sphumi = alpha*fi*esti/(M_mslp[i]-beta*fi*esti);

    /* We need the derivative of sphumi wrt. tsurf */
    double dsphumdesti = alpha/(M_mslp[i]-beta*fi*esti)*( 1. + beta*fi*esti/(M_mslp[i]-beta*fi*esti) );
    double destidT     = ( bi*ci*di-Tsurf*( 2.*ci+Tsurf) )/( di*std::pow(ci+Tsurf,2) )*esti;
    double dfidT       = 2.*Ci*Bi*Tsurf;
    double dsphumidT   = dsphumdesti*(fi*destidT+esti*dfidT);

    // -------------------------------------------------

    /* Density of air */
    double tairK  = M_tair[i] + physical::tfrwK;
    double rhoair = M_mslp[i]/(physical::Ra*tairK) * (1.+sphuma)/(1.+1.609*sphuma);

    /* Sensible heat flux and derivative */
    double Qsh    = drag_ice_t * rhoair * physical::cpa * wspeed*( Tsurf - M_tair[i] );
    double dQshdT = drag_ice_t * rhoair * physical::cpa * wspeed;

    /* Latent heat flux and derivative */
    double Qlh    = drag_ice_t*rhoair*(physical::Lf+physical::Lv0)*wspeed*( sphumi - sphuma );
    double dQlhdT = drag_ice_t*(physical::Lf+physical::Lv0)*rhoair*wspeed*dsphumidT;

    /* Sum them up */
    double Qout    = Qlw_out + Qsh + Qlh;
    dQaidT = dQlwdT + dQshdT + dQlhdT;

    /* Sublimation */
    subl    = Qlh/(physical::Lf+physical::Lv0);

    /* Sum them up */
    Qai = Qsw - Qlw_in + Qout;
}

// Ice-ocean heat flux
// We can do better than this ... but it'll wait
double
FiniteElement::iceOceanHeatflux(double sst, double sss, double mld, double dt)
{
    /* Use all excess heat to melt or grow ice. This is not
     * accurate, but will have to do for now! */
    double const Tbot = -physical::mu*sss; // Temperature at ice base (bottom), also freezing point of sea-water
    return (sst-Tbot)*physical::rhow*physical::cpw*mld/dt;
}

// Albedo
double
FiniteElement::albedo(int alb_scheme, double Tsurf, double hs, double alb_sn, double alb_ice, double I_0)
{
    double albedo;

    /* Calculate albedo - we can impliment different schemes if we want */
    switch ( alb_scheme )
    {
        case 1:
        case 2:
            /* This scheme mimics Semtner 76 and Maykut and Untersteiner 71 when
             * alb_ice = 0.64 and alb_sn = 0.85 */
            if ( hs > 0. )
            {
                /* decrease the albedo towards ice albedo for snow thinner than 0.2 m */
                if ( alb_scheme == 2 )
                    albedo = std::min(alb_sn, alb_ice + (alb_sn-alb_ice)*hs/0.2);
                else
                    albedo = alb_sn;
            } else {
                /* account for penetrating shortwave radiation */
                albedo = alb_ice + 0.4*( 1.-alb_ice )*I_0;
            }
            break;
        case 3:
            { // create a scope for abli, albs, and frac_sn
            /* Albedo scheme from ccsm 3 */
            /* The scheme is simplified using the assumption that visible solar
             * radiation accounts for 52% of the spectrum and the near-infrared for
             * 48% (the same split as used in cice when running in stand-alone
             * mode). */

            /* This is the ccsm3 scheme when alb_ice = 0.538 and alb_sn = 0.8256 */
            double albi, albs;
            if ( Tsurf > -1. )
            {
                albi = alb_ice - 0.075*(Tsurf+1.);
                albs = alb_sn  - 0.124*(Tsurf+1.);
            } else {
                albi = alb_ice;
                albs = alb_sn;
            }

            /* Snow cover fraction */
            double frac_sn = hs/(hs+0.02);

            /* Final albedo */
            albedo = frac_sn*albs + (1.-frac_sn)*albi;

            break;
            }
        default:
            std::cout << "alb_scheme = " << alb_scheme << "\n";
            throw std::logic_error("Wrong albedo_scheme");
    }

    return albedo;
}

// Winton thermo dynamics (ice temperature, growth, and melt)
void
FiniteElement::thermoWinton(int i, double dt, double wspeed, double sphuma, double conc, double voli, double vols, double Qlw_in, double Qsw_in, double mld, double snowfr,
        double &hi, double &hs, double &hi_old, double &Qio, double &del_hi, double &Tsurf, double &T1, double &T2)
{
    // Constants
    double const alb_ice = vm["simul.alb_ice"].as<double>();
    double const alb_sn  = vm["simul.alb_sn"].as<double>();
    double const I_0     = vm["simul.I_0"].as<double>();
    int const alb_scheme = vm["simul.alb_scheme"].as<int>();

    double const drag_ice_t = vm["simul.drag_ice_t"].as<double>();

    bool const flooding = vm["simul.flooding"].as<bool>();

    // Useful volumetric quantities
    double const qi   = physical::Lf * physical::rhoi;
    double const qs   = physical::Lf * physical::rhos;
    double const Crho = physical::C * physical::rhoi;

    double const Tbot     = -physical::mu*M_sss[i]; // Temperature at ice base (bottom), also freezing point of sea-water
    double const Tfr_ice  = -physical::mu*physical::si; // Freezing point of ice

    /* Local variables */
    double Qai, dQaidT, subl;

    /* Don't do anything if there's no ice */
    if ( conc <=0. )
    {
        hi       = 0.;
        hs       = 0.;
        hi_old   = 0.;
        Qio      = 0.;
        del_hi   = 0.;
        Tsurf    = Tbot;
        T1       = Tbot;
        T2       = Tbot;
    } else {
        /* Calculate the slab thickness */
        hi     = voli/conc;
        hi_old = hi;
        hs     = vols/conc;

        double const Tfr_surf = ( hs > 0 ) ? 0. : Tfr_ice; // Freezing point at surface (snow or ice)

        /*
         * Internal temperatures
         * Following Winton (2000) we start by calculating the temperature - first T1, then Tsurf, and finally T2
         * Numers in parentheses refer to equations in the paper
         */

        /* Calculate atmospheric fluxes */
        // Shortwave is modulated by the albedo
        double Qsw = -Qsw_in*(1.-FiniteElement::albedo(alb_scheme, Tsurf, hs, alb_sn, alb_ice, I_0))*(1.-I_0);
        // The rest is calculated by bulk formula
        FiniteElement::atmFluxBulk(i, Tsurf, sphuma, drag_ice_t, Qsw, Qlw_in, wspeed,
                Qai, dQaidT,subl);

        // First some coefficients based on temperatures from the previous time step
        double K12 = 4*physical::ki*physical::ks / ( physical::ks*hi + 4*physical::ki*hs ); // (5)
        double A   = Qai - Tsurf*dQaidT; // (7)
        double B   = dQaidT; // (8)
        double K32 = 2*physical::ki/hi; // (10)

        double A1 = hi*Crho/(2*dt) + K32*( 4*dt*K32 + hi*Crho ) / ( 6*dt*K32 + hi*Crho ) + K12*B/(K12+B); // (16)
        double B1 = -hi/(2*dt) * ( Crho*T1 + qi*Tfr_ice/T1 ) - I_0*Qsw
            - K32*( 4*dt*K32*Tbot + hi*Crho*T2 ) / ( 6*dt*K32 + hi*Crho ) + A*K12/(K12+B); // (17)
        double C1 = hi*qi*Tfr_ice/(2*dt); // (18)

        // Then calculate T1 and Tsurf
        T1    = - ( B1 + std::sqrt(B1*B1-4*A1*C1) ) / ( 2*A1 ); // (21)
        Tsurf = ( K12*T1 - A ) / ( K12 + B ); // (6)

        // Recalculate T1 and Tsurf if we're melting at the surface
        double Msurf = 0.;
        if ( Tsurf > Tfr_surf )
        {
            Tsurf = Tfr_surf;
            // A1 = hi*Crho/(2*dt) + K32*( 4*dt*K32 + hi*Crho ) / ( 6*dt*K32 + hi*Crho ) + K12; // (19)
            A1   += K12 - K12*B/(K12+B);
            // B1 = -hi/(2*dt) * ( Crho*T1 + qi*Tfr_ice/T1 ) - I_0*Qsw
            //     - K32 * ( 4*dt*K32*Tbot + hi*Crho*T2 ) / (6*dt*K32 + hi*Crho ) - K12*Tsurf ; // (20)
            B1   -= K12*Tsurf + A*K12/(K12+B);
            T1    = - ( B1 + std::sqrt(B1*B1-4*A1*C1) ) / ( 2*A1 ); // (21)

            // Surface melt
            Msurf = K12*(T1-Tsurf) - (A+B*Tsurf); // (22)
        }

        // Calculate T2 based on the new T1 and T2 from the previous time step
        T2 = ( 2*dt*K32*(T1+2*Tbot) + hi*Crho*T2 ) / ( 6*dt*K32 + hi*Crho ); // (15)

        /*
         * Thickness changes
         */

        double h1 = hi/2.;
        double h2 = hi/2.;
        double E1 = Crho*(T1 - Tfr_ice) - qi*( 1 - Tfr_ice/T1 ); // (1) - but I've multiplied with rhoi, because it's missing in the paper
        double E2 = Crho*(T2 - Tfr_ice) - qi; // (25) - but I've multiplied with rhoi, because it's missing in the paper

        // Snowfall
        hs += M_precip[i]*snowfr/physical::rhos*dt;

        // Bottom melt/freezing
        Qio    = FiniteElement::iceOceanHeatflux(M_sst[i], M_sss[i], mld, dt);
        double Mbot  = Qio - 4*physical::ki*(Tbot-T2)/hi; // (23)

        // Growth/melt at the ice-ocean interface
        if ( Mbot <= 0. )
        {
            // Growth
            double Ebot  = Crho*(Tbot - Tfr_ice) - qi; // (25) - but I've multiplied with rhoi, because it's missing in the paper
            double delh2 = Mbot*dt/Ebot; // (24)
            T2 = ( delh2*Tbot + h2*T2 ) / ( delh2 + h2 ); // (26)
            h2 += delh2;
        } else {
            // Melt
            double delh2 = -std::min(          -  Mbot*dt/E2,                        h2); // (31) - with added division with rhoi
            double delh1 = -std::min(std::max( -( Mbot*dt + E2*h2 )/E1,         0.), h1); // (32) - with added division with rhoi
            double delhs = -std::min(std::max(  ( Mbot*dt + E2*h2 + E1*h1 )/qs, 0.), hs); // (32) - with added division with rhoi and rhos

            // If everyting melts we need to give back to the ocean
            Qio -= std::max(Mbot*dt - qs*hs + E1*h1 + E2*h2, 0.)/dt; // (34) - with added multiplication of rhoi and rhos and division with dt

            hs += delhs;
            h1 += delh1;
            h2 += delh2;
        }

        // Sublimation at the surface
        if ( subl*dt <= hs*physical::rhos)
            hs -= subl*dt/physical::rhos;
        else if ( subl*dt - hs*physical::rhos <= h1*physical::rhoi )
        {
            hs  = 0.;
            h1 -= (subl*dt - hs*physical::rhos)/physical::rhoi;
        }
        else if ( subl*dt - h1*physical::rhoi - hs*physical::rhos <= h2*physical::rhoi )
        {
            hs  = 0.;
            h1  = 0.;
            h2 -= (subl*dt - h1*physical::rhoi - hs*physical::rhos)/physical::rhoi;
        }
        else
        {
            hs = 0.;
            h1 = 0.;
            h2 = 0.;
            double ocn_evap_err = ( subl*dt - (h1+h2)*physical::rhoi - hs*physical::rhos )/physical::rhow;
			LOG(WARNING) << "All the ice has sublimated. This shouldn't happen and will result in lack of evaporation from the ocean of "
                << ocn_evap_err*1e3 << " mm over the current time step, in element " << i << ".\n";
        }


        // Melting at the surface
        assert(Msurf >= 0); // Sanity check
        double delhs = -std::min(             Msurf*dt/qs,                          hs); // (27) - with division of rhos
        double delh1 = -std::min(std::max( -( Msurf*dt - qs*hs )/E1,           0.), h1); // (28) - with division of rhoi and rhos
        double delh2 = -std::min(std::max( -( Msurf*dt - qs*hs + E1*h1 ) / E2, 0.), h2); // (29) - with division of rhoi and rhos

        // If everyting melts we need to give back to the ocean
        Qio -= std::max(Msurf*dt - qs*hs + E1*h1 + E2*h2, 0.)/dt; // (30) - with multiplication of rhoi and rhos and division with dt

        hs += delhs;
        h1 += delh1;
        h2 += delh2;

        // Snow-to-ice conversion
        double freeboard = ( hi*(physical::rhow-physical::rhoi) - hs*physical::rhos) / physical::rhow;
        if ( flooding && freeboard < 0)
        {
            // double delhs = -std::max( ( hs - (physical::rhow-physical::rhoi)*hi/physical::rhos )*physical::rhoi/physical::rhow, 0.); // (35)
            hs += std::min( freeboard*physical::rhoi/physical::rhos, 0. );
            // double delh1 =  std::max( ( hs - (physical::rhow-physical::rhoi)*hi/physical::rhos )*physical::rhos/physical::rhow, 0.); // (36)
            double delh1 = std::max( -freeboard, 0. );

            double f1   = 1-delh1/(delh1+h1); // Fraction of layer 1 ice in the new upper layer
            double Tbar = f1*( T1 + qi*Tfr_ice/(Crho*T1) ) + (1-f1)*Tfr_ice; // (39)

            T1 = ( Tbar - std::sqrt(Tbar*Tbar - 4*Tfr_ice*qi/Crho) )/2.; // (38)
            h1 += delh1;
        }

        // Even out the layer structure and temperatures
        hi = h1 + h2;
        if ( h2 > h1 )
        {
            // Lower layer ice is added to the upper layer
            // T1 changes, but T2 not
            double f1   = h1/hi*2.; // Fraction of layer 1 ice found in the new layer 1
            double Tbar = f1*( T1 + qi*Tfr_ice/(Crho*T1) ) + (1-f1)*T2; // (39)
            T1 = ( Tbar - std::sqrt(Tbar*Tbar - 4*Tfr_ice*qi/Crho) )/2.; // (38)
        } else {
            // Upper layer ice is added to the lower layer
            // T2 changes, but T1 not
            double f1   = (2.*h1-hi)/hi; // Fraction of layer 1 ice found in new layer 2
            T2 = f1*( T1 + qi*Tfr_ice/(Crho*T1) ) + (1-f1)*T2; // (40)

            // Melt from top and bottom if T2 is too high
            if ( T2 > Tfr_ice )
            {
                // This is:
                // hi -= h2*C*(T2-Tfr_ice) / ( E1 + Ebot );
                // But h2 hasn't been updated, E1 may have changed and Ebot is not in this scope
                // so we just write it out:
                hi -= hi/2*Crho*(T2-Tfr_ice)*T1/( qi*T1 + (Crho*T1-qi)*(Tfr_ice-T1) );
                T2  = Tfr_ice;
            }
        }

        // Book keeping
        del_hi = hi-hi_old;

        /* Make sure we don't get too small hi_new */
        if ( hi < physical::hmin )
        {
            Qio   -= ( -qs*hs + (E1+E2)*hi/2. )/dt; // modified (30) - with multiplication of rhoi and rhos and division with dt

            del_hi = -hi_old;
            hi     = 0.;
            hs     = 0.;
            Tsurf  = Tbot;
            T1     = Tbot;
            T2     = Tbot;
        }
    }
}

// This is Semtner zero layer
void
FiniteElement::thermoIce0(int i, double wspeed, double sphuma, double conc, double voli, double vols, double Qlw_in, double Qsw_in, double mld, double snowfr,
        double &hi, double &hs, double &hi_old, double &Qio, double &del_hi, double &Tsurf)
{

    // Constants
    double const alb_ice = vm["simul.alb_ice"].as<double>();
    double const alb_sn  = vm["simul.alb_sn"].as<double>();
    double const I_0     = vm["simul.I_0"].as<double>();
    int const alb_scheme = vm["simul.alb_scheme"].as<int>();

    double const drag_ice_t = vm["simul.drag_ice_t"].as<double>();

    bool const flooding = vm["simul.flooding"].as<bool>();

    double const qi = physical::Lf * physical::rhoi;
    double const qs = physical::Lf * physical::rhos;

    // local variables
    double Qai = 0;

    /* Don't do anything if there's no ice */
    if ( conc <=0. )
    {
        hi      = 0.;
        hi_old  = 0.;
        hs      = 0.;
        Tsurf   = 0.;
        Qio     = 0.;
        del_hi  = 0.;
        Qai     = 0.;
    } else {
        /* Calculate the slab thickness */
        hi     = voli/conc;
        hi_old = hi;
        hs     = vols/conc;

        /* Local variables */
        double Qsw, Qout, dQaidT, subl;
        double Qic, del_hs, del_ht, del_hb, draft;

        double Qlw_out, dQlwdT;
        double tairK, sphumi;
        double rhoair, Qsh, dQshdT;
        double Qlh, dsphumidT, dQlhdT;

        double fi, esti;
        double dsphumdesti, destidT, dfidT;

        /* ---------------------------------------------------------------
        * Calculate the surface temperature within a while-loop
        * --------------------------------------------------------------- */
        double dtsurf   = 1.;
        double Tbot     = -physical::mu*M_sss[i];
        int nb_iter_while=0;
        while ( dtsurf > 1e-4 )
        {
            nb_iter_while++;

            /* Calculate atmospheric fluxes */
            // Shortwave is modulated by the albedo
            Qsw = -Qsw_in*(1.-FiniteElement::albedo(alb_scheme, Tsurf, hs, alb_sn, alb_ice, I_0))*(1.-I_0);
            // The rest is calculated by bulk formula
            FiniteElement::atmFluxBulk(i, Tsurf, sphuma, drag_ice_t, Qsw, Qlw_in, wspeed,
                    Qai, dQaidT,subl);

            // -------------------------------------------------
            /* Recalculate Tsurf */
            dtsurf = Tsurf;
            Qic    = physical::ks*( Tbot-Tsurf )/( hs + physical::ks*hi/physical::ki );
            Tsurf = Tsurf + ( Qic - Qai )/
                ( physical::ks/(hs+physical::ks*hi/physical::ki) + dQaidT );

            /* Set Tsurf to the freezing point of snow or ice */
            if ( hs > 0. )
                Tsurf = std::min(0., Tsurf);
            else
                Tsurf = std::min(-physical::mu*physical::si, Tsurf);

            /* Re-evaluate the exit condition */
            dtsurf = std::abs(dtsurf-Tsurf);
        }

        if(nb_iter_while>10)
        {
            LOG(DEBUG) << "nb_iter_while = " << nb_iter_while << "\n";
            throw std::logic_error("nb_iter_while larger than 10");
        }

        /* Conductive flux through the ice */
        Qic = physical::ks*( Tbot-Tsurf )/( hs + physical::ks*hi/physical::ki );

        /* ---------------------------------------------------------------
         * Melt and growth
         * --------------------------------------------------------------- */

        /* Top melt */
        /* Snow melt and sublimation */
        del_hs = std::min(Qai-Qic,0.)*time_step/qs - subl*time_step/physical::rhos;
        /* Use the energy left over after snow melts to melt the ice */
        del_ht = std::min(hs+del_hs,0.)*qs/qi;
        /* Can't have negative hs! */
        del_hs = std::max(del_hs,-hs);
        hs  = hs + del_hs + M_precip[i]*snowfr/physical::rhos*time_step;

        /* Heatflux from ocean */
        Qio = FiniteElement::iceOceanHeatflux(M_sst[i], M_sss[i], mld, time_step);
        /* Bottom melt/growth */
        del_hb = (Qic-Qio)*time_step/qi;

        /* Combine top and bottom */
        del_hi = del_ht+del_hb;
        hi     = hi + del_hi;

        /* Snow-to-ice conversion */
        draft = ( hi*physical::rhoi + hs*physical::rhos ) / physical::rhow;
        if ( flooding && draft > hi )
        {
            /* Subtract the mass of snow converted to ice from hs_new */
            hs = hs - ( draft - hi )*physical::rhoi/physical::rhos;
            hi = draft;
        }

        /* Make sure we don't get too small hi_new */
        if ( hi < physical::hmin )
        {
            del_hi  = del_hi-hi;
            Qio     = Qio + hi*qi/time_step + hs*qs/time_step;

            hi      = 0.;
            hs      = 0.;
            Tsurf   = Tbot;
        }
    }
}

// This is the main working function, called from main.cpp (same as perform_simul in the old code)
void
FiniteElement::run()
{
    std::string current_time_system = current_time_local();

    M_export_path = Environment::nextsimDir().string() + "/matlab";
    // change directory for outputs if the option "output_directory" is not empty
    if ( ! (vm["simul.output_directory"].as<std::string>()).empty() )
    {
        M_export_path = vm["simul.output_directory"].as<std::string>();

        fs::path path(M_export_path);
        // add a subdirecory if needed
        // path /= "subdir";

        // create the output directory if it does not exist
        if ( !fs::exists(path) )
            fs::create_directories(path);
    }

    this->writeLogFile();

    int pcpt = this->init();
    int niter = vm["simul.maxiteration"].as<int>();

    // Debug file that records the time step
    std::fstream pcpt_file;
    pcpt_file.open("Timestamp.txt", std::ios::out | std::ios::trunc);

    // main loop for nextsim program
    current_time = time_init + pcpt*time_step/(24*3600.0);
    bool is_running = true;
    while (is_running)
    {
        //std::cout<<"TIME STEP "<< pcpt << " for "<< current_time <<"\n";
        std::cout<<"---------------------- TIME STEP "<< pcpt << " : "
                 << model_time_str(vm["simul.time_init"].as<std::string>(), pcpt*time_step);

        if (!(pcpt % 20))
        {
            std::cout<<" ---------- progression: ("<< 100.0*(pcpt*time_step/duration) <<"%)"
                     <<" ---------- time spent: "<< time_spent(current_time_system);
        }

        std::cout <<"\n";

        is_running = (pcpt*time_step) < duration;

        // if (pcpt > 21)
        // if ( fmod((pcpt+1)*time_step,mooring_output_time_step) == 0 )
        //    is_running = false;

        if (pcpt == niter)
            is_running = false;

        // **********************************************************************
        // Take one time-step
        // **********************************************************************
        this->step(pcpt);

std::cout<< "pcpt= " << pcpt  <<"\n";

        pcpt_file << pcpt << "\n";
        pcpt_file << to_date_string(current_time) << "\n";
        pcpt_file.seekp(0);
    }

    pcpt_file.close();

    if ( pcpt*time_step/output_time_step < 1000 )
        this->exportResults(1000);
    LOG(INFO) <<"TIMER total = " << chrono_tot.elapsed() <<"s\n";

    this->finalise();

    LOG(INFO) << "-----------------------Simulation done on "<< current_time_local() <<"\n";
}

// Finalise everything
void
FiniteElement::finalise()
{
    // Don't forget to close the iabp file!
    if (M_drifter_type == setup::DrifterType::IABP)
    {
        M_iabp_file.close();
        M_drifters_out.close();
    }

    this->clear();
}

// Initialise everything
int
FiniteElement::init()
{
    // Initialise everything that doesn't depend on the mesh (constants, data set description, and time)
    int pcpt = 0;
    mesh_adapt_step=0;
    had_remeshed=false;
    this->initConstant();
    current_time = time_init /*+ pcpt*time_step/(24*3600.0)*/;
    this->initDatasets();

    LOG(INFO) << "-----------------------Simulation started on "<< current_time_local() <<"\n";

    LOG(INFO) <<"TIMESTEP= "<< time_step <<"\n";
    LOG(INFO) <<"DURATION= "<< duration <<"\n";

    // Initialise the mesh
    this->initMesh(M_domain_type, M_mesh_type);

#if defined (WAVES)
    // Extract the WIM grid;
    M_use_wim   = vm["simul.use_wim"].as<bool>();
    if (M_use_wim)
    {
        // initialize wim here to have access to grid information

        // instantiation of wim
        wim = wim_type(vm);

        // initialization of wim
        wim.init();

        // get wim grid
        std::cout<<"Getting WIM grid info\n";
        wim_grid    = wim.wimGrid("km");
    }
#endif

    // Check the minimum angle of the grid
    double minang = this->minAngle(M_mesh);
    if (minang < vm["simul.regrid_angle"].as<double>())
    {
        LOG(INFO) <<"invalid regridding angle: should be smaller than the minimal angle in the intial grid\n";
        throw std::logic_error("invalid regridding angle: should be smaller than the minimal angle in the intial grid");
    }
    if ( M_use_restart )
    {
        LOG(DEBUG) <<"Reading restart file\n";
        pcpt = this->readRestart(vm["setup.step_nb"].as<int>());
        current_time = time_init + pcpt*time_step/(24*3600.0);
        LOG(DEBUG) <<"Initialize forcingAtmosphere\n";
        this->forcingAtmosphere();
        LOG(DEBUG) <<"Initialize forcingOcean\n";
        this->forcingOcean();
#if defined (WAVES)
        LOG(DEBUG) <<"Initialize forcingWave\n";
        this->forcingWave();
#endif
        LOG(DEBUG) <<"Initialize bathymetry\n";
        this->bathymetry();
        chrono.restart();
        LOG(DEBUG) <<"check_and_reload starts\n";
        for ( auto it = M_external_data.begin(); it != M_external_data.end(); ++it )
            (*it)->check_and_reload(M_mesh,current_time);
        LOG(DEBUG) <<"check_and_reload in "<< chrono.elapsed() <<"s\n";
    }
    else
    {
        // Do one regrid to get the mesh right
        LOG(DEBUG) <<"Initial regrid\n";
        this->regrid(pcpt);

        // Initialise variables
        chrono.restart();
        LOG(DEBUG) <<"Initialize variables\n";
        this->initVariables();

        LOG(DEBUG) <<"Initialize forcingAtmosphere\n";
        this->forcingAtmosphere();

        LOG(DEBUG) <<"Initialize forcingOcean\n";
        this->forcingOcean();
#if defined (WAVES)
        LOG(DEBUG) <<"Initialize forcingWave\n";
        this->forcingWave();
#endif
        LOG(DEBUG) <<"Initialize bathymetry\n";
        this->bathymetry();
        chrono.restart();
        LOG(DEBUG) <<"check_and_reload starts\n";
        for ( auto it = M_external_data.begin(); it != M_external_data.end(); ++it )
            (*it)->check_and_reload(M_mesh,time_init);
        LOG(DEBUG) <<"check_and_reload in "<< chrono.elapsed() <<"s\n";
        this->initModelState();
        LOG(DEBUG) <<"initSimulation done in "<< chrono.elapsed() <<"s\n";
    }
    // Open the output file for drifters
    // TODO: Is this the right place to open the file?
    if (M_drifter_type == setup::DrifterType::IABP )
    {
        // We should tag the file name with the init time in case of a re-start.
        std::stringstream filename;
        filename << M_export_path << "/drifters_out_" << current_time << ".txt";
        M_drifters_out.open(filename.str(), std::fstream::out);
        if ( ! M_drifters_out.good() )
            throw std::runtime_error("Cannot write to file: " + filename.str());
    }

    // Initialise the moorings - if requested
    if ( M_use_moorings )
        this->initMoorings();

    return pcpt;
}

// Take one time step
void
FiniteElement::step(int &pcpt)
{



#if defined (WAVES)

    M_run_wim = !(pcpt % vm["nextwim.couplingfreq"].as<int>());

    // coupling with wim (exchange from nextsim to wim)
    if (M_use_wim)
        this->nextsimToWim(pcpt);
#endif

    // step 0: preparation
    // remeshing and remapping of the prognostic variables

    // The first time step we behave as if we just did a regrid
    M_regrid = (pcpt==0);

    if (vm["simul.regrid"].as<std::string>() == "bamg")
    {
        double displacement_factor = 1.;
        double minang = this->minAngle(M_mesh,M_UM,displacement_factor);
        LOG(DEBUG) <<"REGRID ANGLE= "<< minang <<"\n";

        if ( minang < vm["simul.regrid_angle"].as<double>() )
        {
            M_regrid = true;
            // this->writeRestart(pcpt, 0); // Write a restart before regrid - useful for debugging
            if ( M_use_moorings )
                M_moorings.updateGridMean(M_mesh);
            LOG(DEBUG) <<"Regridding starts\n";
            chrono.restart();
            this->regrid(pcpt);
            LOG(DEBUG) <<"Regridding done in "<< chrono.elapsed() <<"s\n";
            if ( M_use_moorings )
                M_moorings.resetMeshMean(M_mesh);
        }
    }

#if defined (WAVES)
    // coupling with wim (exchange from wim to nextsim)
    if (M_use_wim)
        this->wimToNextsim(pcpt);
    else if ( M_regrid || M_use_restart ) // We need to make sure M_tau is the right size
        M_tau.resize(2*M_num_nodes,0.);
#else
    // just set this vector to the right size,
    // but fill with 0
    M_tau.resize(2*M_num_nodes,0.);
#endif


    // Read in the new buoys and output
    if (M_drifter_type == setup::DrifterType::IABP && std::fmod(current_time,0.5) == 0)
    {
        this->updateIABPDrifter();
        // TODO: Do we want to output drifters at a different time interval?
        this->outputDrifter(M_drifters_out);
    }

    if ( M_regrid || M_use_restart )
    {
        chrono.restart();
        LOG(DEBUG) <<"tensors starts\n";
        this->tensors();
        LOG(DEBUG) <<"tensors done in "<< chrono.elapsed() <<"s\n";
        chrono.restart();
        LOG(DEBUG) <<"cohesion starts\n";
        this->cohesion();
        LOG(DEBUG) <<"cohesion done in "<< chrono.elapsed() <<"s\n";
        chrono.restart();
        LOG(DEBUG) <<"Coriolis starts\n";
        this->coriolis();
        LOG(DEBUG) <<"Coriolis done in "<< chrono.elapsed() <<"s\n";
    }

    chrono.restart();
    LOG(DEBUG) <<"check_and_reload starts\n";
    for ( auto it = M_external_data.begin(); it != M_external_data.end(); ++it )
        (*it)->check_and_reload(M_mesh,current_time+time_step/(24*3600.0));
    LOG(DEBUG) <<"check_and_reload in "<< chrono.elapsed() <<"s\n";

    M_use_restart = false;

#if 1
    if (pcpt == 0)
    {
        nodeMaxConc(); // needed for post-processing, recomputed in assemble() for the other steps
        chrono.restart();
        LOG(DEBUG) <<"first export starts\n";
        this->exportResults(0);
        LOG(DEBUG) <<"first export done in " << chrono.elapsed() <<"s\n";
    }
#endif

    //======================================================================
    // Do the thermodynamics
    //======================================================================
    if(vm["simul.use_thermo_forcing"].as<bool>())
    {
        chrono.restart();
        LOG(DEBUG) <<"thermo starts\n";
        this->thermo();
        LOG(DEBUG) <<"thermo done in "<< chrono.elapsed() <<"s\n";
    }

    //======================================================================
    // Assemble the matrix
    //======================================================================

    this->assemble(pcpt);

    //======================================================================
    // Solve the linear problem
    //======================================================================

    if(had_remeshed && (vm["simul.regrid_output_flag"].as<bool>()))
    {
        had_remeshed=false;
        this->exportResults(300000+mesh_adapt_step);
        had_remeshed=false;
    }

    this->solve();

    chrono.restart();
    LOG(DEBUG) <<"updateVelocity starts\n";
    this->updateVelocity();
    LOG(DEBUG) <<"updateVelocity done in "<< chrono.elapsed() <<"s\n";

    chrono.restart();
    LOG(DEBUG) <<"update starts\n";
    this->update();
    LOG(DEBUG) <<"update done in "<< chrono.elapsed() <<"s\n";

    ++pcpt;
    current_time = time_init + pcpt*time_step/(24*3600.0);

#if 1
    if(fmod(pcpt*time_step,output_time_step) == 0)
    {
        chrono.restart();
        LOG(DEBUG) <<"export starts\n";
        this->exportResults((int) pcpt*time_step/output_time_step);
        LOG(DEBUG) <<"export done in " << chrono.elapsed() <<"s\n";
    }

    if ( M_use_moorings )
    {
        this->updateMeans(M_moorings);
        if ( fmod(pcpt*time_step,mooring_output_time_step) == 0 )
        {
            M_moorings.updateGridMean(M_mesh);
            //M_moorings.exportGridMeans("_grid.dat", time_step, mooring_output_time_step);
            M_moorings.appendNetCDF(M_moorings_file, current_time, time_step, mooring_output_time_step);

            M_moorings.resetMeshMean(M_mesh);
            M_moorings.resetGridMean();
        }
    }

#endif

    if ( fmod(pcpt*time_step,restart_time_step) == 0)
    {
        std::cout << "Writing restart file after time step " <<  pcpt-1 << "\n";
        this->writeRestart(pcpt, (int) pcpt*time_step/restart_time_step);
    }
}

// Add to the mean on the mesh
void
FiniteElement::updateMeans(GridOutput &means)
{
    // Update elements
    for ( auto it=means.M_elemental_variables.begin(); it!=means.M_elemental_variables.end(); ++it )
    {
        switch (it->variableID)
        {
            case (GridOutput::variableID::conc):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_conc[i];
                break;

            case (GridOutput::variableID::thick):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_thick[i];
                break;

            case (GridOutput::variableID::damage):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_damage[i];
                break;

            case (GridOutput::variableID::snow_thick):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_snow_thick[i];
                break;

            case (GridOutput::variableID::tsurf):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_conc[i]*M_tice[0][i] + (1-M_conc[i])*M_sst[i];
                break;

            case (GridOutput::variableID::sst):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_sst[i];
                break;

            case (GridOutput::variableID::sss):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_sss[i];
                break;

            case (GridOutput::variableID::tsurf_ice):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_tice[0][i];
                break;

            case (GridOutput::variableID::t1):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_tice[1][i];
                break;

            case (GridOutput::variableID::t2):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_tice[2][i];
                break;

            default: std::logic_error("Updating of given variableID not implimented (elements)");
        }
    }

    // Update nodes
    for ( auto it=means.M_nodal_variables.begin(); it!=means.M_nodal_variables.end(); ++it )
    {
        switch (it->variableID)
        {
            case (GridOutput::variableID::VT_x):
                for (int i=0; i<M_num_nodes; i++)
                    it->data_mesh[i] += M_VT[i];
                break;

            case (GridOutput::variableID::VT_y):
                for (int i=0; i<M_num_nodes; i++)
                    it->data_mesh[i] += M_VT[i+M_num_nodes];
                break;

            default: std::logic_error("Updating of given variableID not implimented (nodes)");
        }
    }
}

// Initialise everything w.r.t. the moorings
void
FiniteElement::initMoorings()
{
    // Output dimensions
    DataSet::Dimension dimension_x={
        name:"x",
        cyclic:false
    };

    DataSet::Dimension dimension_y={
        name:"y",
        cyclic:false
    };

    DataSet::Dimension dimension_time={
        name:"time",
        cyclic:false
    };

    std::vector<DataSet::Dimension> dimensions(3);
    dimensions[0] = dimension_x;
    dimensions[1] = dimension_y;
    dimensions[2] = dimension_time;

    // Output and averaging grids
    std::vector<double> data_nodes(M_num_nodes);
    std::vector<double> data_elements(M_num_elements);
    std::vector<double> data_grid;

    // Output variables - elements
    GridOutput::Variable conc={
        name:"sic",
        longName:"Sea Ice Concentration",
        stdName:"sea_ice_area_fraction",
        dimensions: dimensions,
        Units:"1",
        data_mesh:data_elements,
        data_grid:data_grid,
        variableID: GridOutput::variableID::conc
    };

    GridOutput::Variable thick={
        name:"sit",
        longName:"Sea Ice Thickness",
        stdName:"sea_ice_thickness",
        dimensions: dimensions,
        Units:"m",
        data_mesh:data_elements,
        data_grid:data_grid,
        variableID: GridOutput::variableID::thick
    };

    GridOutput::Variable snow_thick={
        name:"snt",
        longName:"Surface Snow Thickness",
        stdName:"surface_snow_thickness",
        dimensions: dimensions,
        Units:"m",
        data_mesh:data_elements,
        data_grid:data_grid,
        variableID: GridOutput::variableID::snow_thick
    };

    std::vector<GridOutput::Variable> elemental_variables(3);
    elemental_variables[0] = conc;
    elemental_variables[1] = thick;
    elemental_variables[2] = snow_thick;

    if ( M_thermo_type == setup::ThermoType::WINTON )
    {
        // Technically, cannonical units for temperature are K, but that's just anoying here!
        GridOutput::Variable tsurf={
            name:"ts",
            longName:"Surface Temperature",
            stdName:"surface_temperature",
            dimensions: dimensions,
            Units:"degC",
            data_mesh:data_elements,
            data_grid:data_grid,
            variableID: GridOutput::variableID::tsurf
        };

        elemental_variables.push_back(tsurf);

        GridOutput::Variable t1={
            name:"t1",
            longName:"Ice Temperature 1",
            stdName:"ice_temperature_1",
            dimensions: dimensions,
            Units:"degC",
            data_mesh:data_elements,
            data_grid:data_grid,
            variableID: GridOutput::variableID::t1
        };

        elemental_variables.push_back(t1);

        GridOutput::Variable t2={
            name:"t2",
            longName:"Ice Temperature 2",
            stdName:"ice_temperature_2",
            dimensions: dimensions,
            Units:"degC",
            data_mesh:data_elements,
            data_grid:data_grid,
            variableID: GridOutput::variableID::t2
        };

        elemental_variables.push_back(t2);
    }

    // Output variables - nodes
    GridOutput::Variable siu={
        name:"siu",
        // longName:"Eastward Sea Ice Velocity",
        // stdName:"eastward_sea_ice_velocity",
        longName:"Sea Ice X Velocity",
        stdName:"sea_ice_x_velocity",
        dimensions: dimensions,
        Units:"m s-1",
        data_mesh:data_nodes,
        data_grid:data_grid,
        variableID: GridOutput::variableID::VT_x
    };

    GridOutput::Variable siv={
        name:"siv",
        // longName:"Northward Sea Ice Velocity",
        // stdName:"northward_sea_ice_velocity",
        longName:"Sea Ice Y Velocity",
        stdName:"sea_ice_y_velocity",
        dimensions: dimensions,
        Units:"m s-1",
        data_mesh:data_nodes,
        data_grid:data_grid,
        variableID: GridOutput::variableID::VT_y
    };

    std::vector<GridOutput::Variable> nodal_variables(2);
    nodal_variables[0] = siu;
    nodal_variables[1] = siv;

    // The vectorial variables are (always on the nodes) ...
    std::vector<int> siuv_id(2);
    siuv_id[0] = 0;
    siuv_id[1] = 1;

    DataSet::Vectorial_Variable siuv{
        components_Id: siuv_id,
        // east_west_oriented: true
        east_west_oriented: false
    };

    std::vector<DataSet::Vectorial_Variable> vectorial_variables(1);
    vectorial_variables[0] = siuv;

    M_moorings_file = M_export_path + "/Moorings.nc";
#if 1
    // Calculate the grid spacing (assuming a regular grid for now)
    auto RX = M_mesh.coordX();
    auto RY = M_mesh.coordY();
    auto xcoords = std::minmax_element( RX.begin(), RX.end() );
    auto ycoords = std::minmax_element( RY.begin(), RY.end() );

    double mooring_spacing = 1e3 * vm["simul.mooring_spacing"].as<double>();
    int nrows = (int) ( 0.5 + ( *xcoords.second - *xcoords.first )/mooring_spacing );
    int ncols = (int) ( 0.5 + ( *ycoords.second - *ycoords.first )/mooring_spacing );

    // Define the mooring dataset
    M_moorings = GridOutput(ncols, nrows, mooring_spacing, *xcoords.first, *ycoords.first, nodal_variables, elemental_variables, vectorial_variables);

    // Save the grid info - this is still just an ascii dump!
    std::ofstream myfile;
    myfile.open("lon_grid.dat");
    std::copy(M_moorings.M_grid.gridLON.begin(), M_moorings.M_grid.gridLON.end(), ostream_iterator<float>(myfile," "));
    myfile.close();
    myfile.open("lat_grid.dat");
    std::copy(M_moorings.M_grid.gridLAT.begin(), M_moorings.M_grid.gridLAT.end(), ostream_iterator<float>(myfile," "));
    myfile.close();

    M_moorings.initNetCDF(M_moorings_file);
#else
    // Read the grid in from file
    std::vector<DataSet::Dimension> dimensions_latlon(2);
    dimensions_latlon[0] = dimension_y;
    dimensions_latlon[1] = dimension_x;

    std::vector<std::vector<double>> data2_tmp;
    data2_tmp.resize(2);

    DataSet::Variable latitude={
        name: "latitude",
        dimensions: dimensions_latlon,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        a: 1.,
        b: 0.,
        Units: "degree_north",
        data2: data2_tmp};

    DataSet::Variable longitude={
        name: "longitude",
        dimensions: dimensions_latlon,
        land_mask_defined: false,
        land_mask_value: 0.,
        NaN_mask_defined: false,
        NaN_mask_value: 0.,
        a: 1.,
        b: 0.,
        Units: "degree_east",
        data2: data2_tmp};

    DataSet::Grid grid={
        interpolation_method: InterpolationType::None,
        interp_type: -1,
        dirname: "data",
        //filename: "TP4DAILY_200803_3m.nc",
        prefix:"TP4DAILY_200803",
        postfix: "_3m.nc",

        latitude: latitude,
        longitude: longitude,

        dimension_x: dimension_x,
        dimension_y: dimension_y,

        mpp_file: vm["simul.proj_filename"].as<std::string>(),
        interpolation_in_latlon: false,

        loaded: false,
        monthly_dataset:true,

        masking: false,
        masking_variable: latitude
    };

    // Define the mooring dataset
    M_moorings = GridOutput(grid, nodal_variables, elemental_variables, vectorial_variables);

    // Save the grid info - this is still just an ascii dump!
    std::ofstream myfile;
    myfile.open("lon_grid.dat");
    std::copy(M_moorings.M_grid.gridLON.begin(), M_moorings.M_grid.gridLON.end(), ostream_iterator<float>(myfile," "));
    myfile.close();
    myfile.open("lat_grid.dat");
    std::copy(M_moorings.M_grid.gridLAT.begin(), M_moorings.M_grid.gridLAT.end(), ostream_iterator<float>(myfile," "));
    myfile.close();
#endif

#if 0
    // Prepare the moorings grid for output
    std::vector<DataSet::Dimension> dimensions_2d(2);
    dimensions_2d[0] = dimension_x;
    dimensions_2d[1] = dimension_y;

    GridOutput::Variable lat={
        name: "lat",
        longName: "Latitude",
        stdName: "latitude",
        dimensions: dimensions_2d,
        Units: "degrees",
        data_mesh: M_mesh.lat(),
        data_grid: data_grid
    };

    GridOutput::Variable lon={
        name: "lon",
        longName: "Longitude",
        stdName: "longitude",
        dimensions: dimensions_2d,
        Units: "degrees",
        data_mesh: M_mesh.lon(),
        data_grid: data_grid
    };

    std::vector<GridOutput::Variable> grid_variables(2);
    grid_variables[0] = lon;
    grid_variables[1] = lat;

    M_moorings_grid = GridOutput(ncols, nrows, mooring_spacing, grid_variables, GridOutput::variableKind::nodal);
    M_moorings_grid.updateGridMean(M_mesh);
    M_moorings_grid.exportGridMeans("_grid.dat", 1., 1.);
#endif

} //initMoorings

void
FiniteElement::writeRestart(int pcpt, int step)
{
    Exporter exporter;
    std::string filename;

    // === Start with the mesh ===
    // First the data
    std::string directory = Environment::nextsimDir().string() + "/restart";
    // change directory for outputs if the option "output_directory" is not empty
    if ( ! (vm["simul.output_directory"].as<std::string>()).empty() )
        directory = vm["simul.output_directory"].as<std::string>() + "/restart";

    // create the output directory if it does not exist
    fs::path path(directory);
    if ( !fs::exists(path) )
        fs::create_directories(path);

    filename = (boost::format( "%1%/mesh_%2%.bin" )
            % directory
            % step ).str();

    std::fstream meshbin(filename, std::ios::binary | std::ios::out | std::ios::trunc);
    if ( ! meshbin.good() )
        throw std::runtime_error("Cannot write to file: " + filename);
    exporter.writeMesh(meshbin, M_mesh);
    meshbin.close();

    // Then the record
    filename = (boost::format( "%1%/mesh_%2%.dat" )
           % directory
           % step ).str();

    std::fstream meshrecord(filename, std::ios::out | std::ios::trunc);
    if ( ! meshrecord.good() )
        throw std::runtime_error("Cannot write to file: " + filename);
    exporter.writeRecord(meshrecord,"mesh");
    meshrecord.close();

    // === Write the prognostic variables ===
    // First the data
    filename = (boost::format( "%1%/field_%2%.bin" )
               % directory
               % step ).str();
    std::fstream outbin(filename, std::ios::binary | std::ios::out | std::ios::trunc );
    if ( ! outbin.good() )
        throw std::runtime_error("Cannot write to file: " + filename);

    std::vector<int> misc_int(2);
    misc_int[0] = pcpt;
    misc_int[1] = M_flag_fix;
    exporter.writeField(outbin, misc_int, "Misc_int");
    exporter.writeField(outbin, M_dirichlet_flags, "M_dirichlet_flags");

    exporter.writeField(outbin, M_conc, "M_conc");
    exporter.writeField(outbin, M_thick, "M_thick");
    exporter.writeField(outbin, M_snow_thick, "M_snow_thick");
    exporter.writeField(outbin, M_sigma, "M_sigma");
    exporter.writeField(outbin, M_damage, "M_damage");
    exporter.writeField(outbin, M_divergence_rate, "M_divergence_rate");
    exporter.writeField(outbin, M_h_ridged_thin_ice, "M_h_ridged_thin_ice");
    exporter.writeField(outbin, M_h_ridged_thick_ice, "M_h_ridged_thick_ice");
    exporter.writeField(outbin, M_random_number, "M_random_number");
    int i=0;
    for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
    {
        exporter.writeField(outbin, *it, "M_Tice_"+std::to_string(i));
        i++;
    }
    exporter.writeField(outbin, M_sst, "M_sst");
    exporter.writeField(outbin, M_sss, "M_sss");
    exporter.writeField(outbin, M_VT, "M_VT");
    exporter.writeField(outbin, M_VTM, "M_VTM");
    exporter.writeField(outbin, M_VTMM, "M_VTMM");
    exporter.writeField(outbin, M_UM, "M_UM");

    if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
    {
        exporter.writeField(outbin, M_h_thin, "M_h_thin");
        exporter.writeField(outbin, M_hs_thin, "M_hs_thin");
        exporter.writeField(outbin, M_tsurf_thin, "M_tsurf_thin");
    }

    if (M_drifter_type == setup::DrifterType::IABP)
    {
        std::vector<int> drifter_no(M_drifter.size());
        std::vector<double> drifter_x(M_drifter.size());
        std::vector<double> drifter_y(M_drifter.size());

        int j=0;
        for ( auto it = M_drifter.begin(); it != M_drifter.end(); ++it )
        {
            drifter_no[j] = it->first;
            drifter_x[j] = it->second[0];
            drifter_y[j] = it->second[1];
            ++j;
        }

        exporter.writeField(outbin, drifter_no, "Drifter_no");
        exporter.writeField(outbin, drifter_x, "Drifter_x");
        exporter.writeField(outbin, drifter_y, "Drifter_y");
    }

    outbin.close();

    // Then the record
    filename = (boost::format( "%1%/field_%2%.dat" )
               % directory
               % step ).str();
    std::fstream outrecord(filename, std::ios::out | std::ios::trunc);
    if ( ! outrecord.good() )
        throw std::runtime_error("Cannot write to file: " + filename);
    exporter.writeRecord(outrecord);
    outrecord.close();
}

int
FiniteElement::readRestart(int step)
{
    Exporter exp_field, exp_mesh;
    std::string filename;
    boost::unordered_map<std::string, std::vector<int>>    field_map_int;
    boost::unordered_map<std::string, std::vector<double>> field_map_dbl;

    // === Read in the mesh restart files ===
    // Start with the record
    filename = (boost::format( "%1%/restart/mesh_%2%.dat" )
               % Environment::nextsimDir().string()
               % step ).str();
    std::ifstream meshrecord(filename);
    if ( ! meshrecord.good() )
        throw std::runtime_error("File not found: " + filename);

    exp_mesh.readRecord(meshrecord);
    meshrecord.close();

    // Then onto the data itself
    filename = (boost::format( "%1%/restart/mesh_%2%.bin" )
               % Environment::nextsimDir().string()
               % step ).str();
    std::fstream meshbin(filename, std::ios::binary | std::ios::in );
    if ( ! meshbin.good() )
        throw std::runtime_error("File not found: " + filename);
    exp_mesh.loadFile(meshbin, field_map_int, field_map_dbl);
    meshbin.close();

    std::vector<int>   indexTr = field_map_int["Elements"];
    std::vector<double> coordX = field_map_dbl["Nodes_x"];
    std::vector<double> coordY = field_map_dbl["Nodes_y"];
    std::vector<int>   nodeId = field_map_int["id"];

    // === Read in the prognostic variables ===
    // Start with the record
    filename = (boost::format( "%1%/restart/field_%2%.dat" )
               % Environment::nextsimDir().string()
               % step ).str();
    std::ifstream inrecord(filename);
    if ( ! inrecord.good() )
        throw std::runtime_error("File not found: " + filename);

    exp_field.readRecord(inrecord);
    inrecord.close();

    // Then onto the data itself
    filename = (boost::format( "%1%/restart/field_%2%.bin" )
               % Environment::nextsimDir().string()
               % step ).str();
    std::fstream inbin(filename, std::ios::binary | std::ios::in );
    if ( ! inbin.good() )
        throw std::runtime_error("File not found: " + filename);
    field_map_int.clear();
    field_map_dbl.clear();
    exp_field.loadFile(inbin, field_map_int, field_map_dbl);
    inbin.close();

    // === Recreate the mesh ===
    // Create bamgmesh and bamggeom
    BamgConvertMeshx(
            bamgmesh,bamggeom,
	    &indexTr[0],&coordX[0],&coordY[0],
	    coordX.size(), indexTr.size()/3.
            );

    // Fix boundaries
    int pcpt   = field_map_int["Misc_int"].at(0);
    M_flag_fix = field_map_int["Misc_int"].at(1);
    std::vector<int> dirichlet_flags = field_map_int["M_dirichlet_flags"];
    for (int edg=0; edg<bamgmesh->EdgesSize[0]; ++edg)
    {
        int fnd = bamgmesh->Edges[3*edg]-1;
        if ((std::binary_search(dirichlet_flags.begin(),dirichlet_flags.end(),fnd)))
        {
            bamggeom->Edges[3*edg+2] = M_flag_fix;
            bamgmesh->Edges[3*edg+2] = M_flag_fix;
        }
    }

    // Import the bamg structs
    this->importBamg(bamgmesh);

    M_mesh.setId(nodeId);

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_num_nodes = M_mesh.numNodes();

    // Initialise all the variables to zero
    this->initVariables();

    // update dirichlet nodes
    M_boundary_flags.resize(0);
    M_dirichlet_flags.resize(0);
    for (int edg=0; edg<bamgmesh->EdgesSize[0]; ++edg)
    {
        M_boundary_flags.push_back(bamgmesh->Edges[3*edg]-1);
        if (bamgmesh->Edges[3*edg+2] == M_flag_fix)
            M_dirichlet_flags.push_back(bamgmesh->Edges[3*edg]-1);
    }

    std::sort(M_dirichlet_flags.begin(), M_dirichlet_flags.end());
    std::sort(M_boundary_flags.begin(), M_boundary_flags.end());

    M_neumann_flags.resize(0);
    std::set_difference(M_boundary_flags.begin(), M_boundary_flags.end(),
                        M_dirichlet_flags.begin(), M_dirichlet_flags.end(),
                        std::back_inserter(M_neumann_flags));

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

    M_surface.assign(M_num_elements,0.);
    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        M_surface[cpt] = this->measure(*it,M_mesh);
        ++cpt;
    }

    // === Set the prognostic variables ===
    M_conc       = field_map_dbl["M_conc"];
    M_thick      = field_map_dbl["M_thick"];
    M_snow_thick = field_map_dbl["M_snow_thick"];
    M_sigma      = field_map_dbl["M_sigma"];
    M_damage     = field_map_dbl["M_damage"];
    M_divergence_rate    = field_map_dbl["M_divergence_rate"];
    M_h_ridged_thin_ice  = field_map_dbl["M_h_ridged_thin_ice"];
    M_h_ridged_thick_ice = field_map_dbl["M_h_ridged_thick_ice"];
    M_random_number      = field_map_dbl["M_random_number"];
    int i=0;
    for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
    {
        *it = field_map_dbl["M_Tice_"+std::to_string(i)];
        i++;
    }
    M_sst        = field_map_dbl["M_sst"];
    M_sss        = field_map_dbl["M_sss"];
    M_VT         = field_map_dbl["M_VT"];
    M_VTM        = field_map_dbl["M_VTM"];
    M_VTMM       = field_map_dbl["M_VTMM"];
    M_UM         = field_map_dbl["M_UM"];

    //for (int i=0; i < M_thick.size(); i++)
    //{
    //  M_thick[i] *= 2.0;
    //}

    if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
    {
        M_h_thin     = field_map_dbl["M_h_thin"];
        M_hs_thin    = field_map_dbl["M_hs_thin"];
        M_tsurf_thin = field_map_dbl["M_tsurf_thin"];
    }

    if (M_drifter_type == setup::DrifterType::IABP)
    {
        std::cout << "Try reading them from the restart file\n";
        std::vector<int>    drifter_no = field_map_int["Drifter_no"];
        std::vector<double> drifter_x  = field_map_dbl["Drifter_x"];
        std::vector<double> drifter_y  = field_map_dbl["Drifter_y"];

        if (drifter_no.size() == 0)
        {
            LOG(WARNING) << "Warning: Couldn't read drifter positions from restart file. Drifter positions initialised as if there was no restart.\n";
            this->initDrifter();
        } else {
            for ( int i=0; i<drifter_no.size(); ++i )
            {
                M_drifter.emplace(drifter_no[i], std::array<double,2>{drifter_x[i], drifter_y[i]});
            }
        }
    }

    inbin.close();

    // Set the target size for the data sets
    M_atmosphere_nodes_dataset.target_size=M_num_nodes;
    M_atmosphere_elements_dataset.target_size=M_num_elements;
    M_atmosphere_bis_elements_dataset.target_size=M_num_elements;
    M_ocean_nodes_dataset.target_size=M_num_nodes;
    M_ocean_elements_dataset.target_size=M_num_elements;

    M_ice_topaz_elements_dataset.target_size=M_num_elements;
    M_ice_piomas_elements_dataset.target_size=M_num_elements;
    M_ice_amsre_elements_dataset.target_size=M_num_elements;
    M_ice_osisaf_elements_dataset.target_size=M_num_elements;
    M_ice_amsr2_elements_dataset.target_size=M_num_elements;
    M_ice_cs2_smos_elements_dataset.target_size=M_num_elements;
    M_bathymetry_elements_dataset.target_size=M_num_elements;
#if defined (WAVES)
    M_WW3A_elements_dataset.target_size=M_num_elements;
    M_ERAIW_1DEG_elements_dataset.target_size=M_num_elements;
#endif
    return pcpt;
}//readRestart

void
FiniteElement::updateVelocity()
{
    M_VTMM = M_VTM;
    M_VTM = M_VT;
    M_VT = M_solution->container();

    // TODO (updateVelocity) Sylvain: This limitation cost about 1/10 of the solver time.
    // TODO (updateVelocity) Sylvain: We could add a term in the momentum equation to avoid the need of this limitation.
    //std::vector<double> speed_c_scaling_test(bamgmesh->NodalElementConnectivitySize[0]);
    #if 1
    int elt_num, i, j;
    double c_max_nodal_neighbour;
    double speed_c_scaling;

    std::vector<double> cloc_elts(bamgmesh->NodalElementConnectivitySize[1]);

    for (i=0; i<bamgmesh->NodalElementConnectivitySize[0]; ++i)
    {
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

        c_max_nodal_neighbour = *std::max_element(cloc_elts.begin(),cloc_elts.begin()+j-1);
        c_max_nodal_neighbour /= vm["simul.drift_limit_concentration"].as<double>();
        speed_c_scaling = std::min(1.,c_max_nodal_neighbour);
        //std::cout<<"c_max_nodal_neighbour["<< i <<"]= "<< c_max_nodal_neighbour <<"\n";
        //std::cout<<"speed_c_scaling["<< i <<"]= "<< speed_c_scaling <<"\n";
        //speed_c_scaling_test[i] = speed_c_scaling;

        // linear scaling of ice velocity
        M_VT[i] *= speed_c_scaling;
        M_VT[i+M_num_nodes] *= speed_c_scaling;
    }
    #endif

    // std::cout<<"MAX SPEED= "<< *std::max_element(speed_c_scaling_test.begin(),speed_c_scaling_test.end()) <<"\n";
    // std::cout<<"MIN SPEED= "<< *std::min_element(speed_c_scaling_test.begin(),speed_c_scaling_test.end()) <<"\n";
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
FiniteElement::forcingAtmosphere()//(double const& u, double const& v)
{
    switch (M_atmosphere_type)
    {
        case setup::AtmosphereType::CONSTANT:
            M_wind=ExternalData(
                vm["simul.constant_wind_u"].as<double>(),
                vm["simul.constant_wind_v"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(vm["simul.constant_tair"].as<double>());
            M_external_data.push_back(&M_tair);

            M_mixrat=ExternalData(vm["simul.constant_mixrat"].as<double>());
            M_external_data.push_back(&M_mixrat);

            M_mslp=ExternalData(vm["simul.constant_mslp"].as<double>());
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(vm["simul.constant_Qsw_in"].as<double>());
            M_external_data.push_back(&M_Qsw_in);

            M_Qlw_in=ExternalData(vm["simul.constant_Qlw_in"].as<double>());
            M_external_data.push_back(&M_Qlw_in);

            M_snowfr=ExternalData(vm["simul.constant_snowfr"].as<double>());
            M_external_data.push_back(&M_snowfr);

            M_precip=ExternalData(vm["simul.constant_precip"].as<double>());
            M_external_data.push_back(&M_precip);

            M_dair=ExternalData(vm["simul.constant_dair"].as<double>());
            M_external_data.push_back(&M_dair);
        break;

        case setup::AtmosphereType::ASR:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_tair);

            M_mixrat=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_external_data.push_back(&M_mixrat);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_Qsw_in);

            M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            M_external_data.push_back(&M_Qlw_in);

            M_snowfr=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_external_data.push_back(&M_snowfr);

            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
            M_external_data.push_back(&M_precip);

            M_dair=ExternalData(-1.);
            M_external_data.push_back(&M_dair);
        break;

        case setup::AtmosphereType::ERAi:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            /*
            variables[0] = tair;
            variables[1] = dair;
            variables[2] = mslp;
            variables[3] = Qsw_in;
            variables[4] = tcc;
            variables[5] = precip;
            */

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_tair);

            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_external_data.push_back(&M_dair);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_Qsw_in);

            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            M_external_data.push_back(&M_tcc);

            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_external_data.push_back(&M_precip);

            M_mixrat=ExternalData(-1.);
            M_external_data.push_back(&M_mixrat);

        break;

        case setup::AtmosphereType::EC:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_tair);

            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_external_data.push_back(&M_dair);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_tcc);

            // Syl: The following two lines should be removed when approxSW will be implemented in Thermo()
            M_Qsw_in=ExternalData(vm["simul.constant_Qsw_in"].as<double>());
            M_external_data.push_back(&M_Qsw_in);

            M_precip=ExternalData(0.);
            M_external_data.push_back(&M_precip);

            M_mixrat=ExternalData(-1.);
            M_external_data.push_back(&M_mixrat);
        break;

        case setup::AtmosphereType::EC_ERAi:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_tair);

            M_dair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_external_data.push_back(&M_dair);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_tcc=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_tcc);

            M_Qsw_in=ExternalData(&M_atmosphere_bis_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_Qsw_in);

            M_precip=ExternalData(&M_atmosphere_bis_elements_dataset,M_mesh,5,false,time_init);
            M_external_data.push_back(&M_precip);

            M_mixrat=ExternalData(-1.);
            M_external_data.push_back(&M_mixrat);
        break;

        default:
            std::cout << "invalid wind forcing"<<"\n";
            throw std::logic_error("invalid wind forcing");
    }
}

void
FiniteElement::forcingOcean()//(double const& u, double const& v)
{
    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
            M_ocean=ExternalData(
                vm["simul.constant_ocean_v"].as<double>(),
                vm["simul.constant_ocean_v"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ocean);

            M_ssh=ExternalData(vm["simul.constant_ssh"].as<double>(),
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ssh);

            M_ocean_temp=ExternalData(vm["simul.constant_ocean_temp"].as<double>());
            M_external_data.push_back(&M_ocean_temp);

            M_ocean_salt=ExternalData(vm["simul.constant_ocean_salt"].as<double>());
            M_external_data.push_back(&M_ocean_salt);

            M_mld=ExternalData(vm["simul.constant_mld"].as<double>());
            M_external_data.push_back(&M_mld);
            break;
        case setup::OceanType::TOPAZR: case setup::OceanType::TOPAZF:
            M_ocean=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 0, true,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ocean);

            M_ssh=ExternalData(
                &M_ocean_nodes_dataset, M_mesh, 2, false,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_ssh);

            M_ocean_temp=ExternalData(&M_ocean_elements_dataset, M_mesh, 0,false,time_init);
            M_external_data.push_back(&M_ocean_temp);

            M_ocean_salt=ExternalData(&M_ocean_elements_dataset, M_mesh, 1,false,time_init);
            M_external_data.push_back(&M_ocean_salt);

            M_mld=ExternalData(&M_ocean_elements_dataset, M_mesh, 2,false,time_init);
            M_external_data.push_back(&M_mld);
            // SYL: there was a capping of the mld at minimum vm["simul.constant_mld"].as<double>()
            // but Einar said it is not necessary, so it is not implemented
    		break;

        default:
            std::cout << "invalid ocean forcing"<<"\n";
            throw std::logic_error("invalid ocean forcing");
    }
}

#if defined (WAVES)
void
FiniteElement::forcingWave()
{
    switch (M_wave_type)
    {
        case setup::WaveType::CONSTANT:

//            std::cout<<"simul.constant_significant_wave_height"<<"\n";
//            std::cout<<vm["simul.constant_significant_wave_height"].as<double>()<<"\n";
//            std::cout<<"simul.constant_significant_wave_height"<<"\n";
            M_SWH=ExternalData(vm["simul.constant_significant_wave_height"].as<double>());
            M_external_data.push_back(&M_SWH);

            M_MWD=ExternalData(vm["simul.constant_wave_mean_direction"].as<double>());
            M_external_data.push_back(&M_MWD);

            M_FP=ExternalData(vm["simul.constant_wave_peak_frequency"].as<double>());
            M_external_data.push_back(&M_FP);

            //wim_forcing_options = M_WW3A_elements_dataset.grid.waveOptions;
            wim_ideal_forcing   = true;
		break;

        //std::cout << age[0] << std::endl;
        case setup::WaveType::WW3A:

	        M_SWH=ExternalData(&M_WW3A_elements_dataset, M_mesh, 0,false,time_init);
            M_external_data.push_back(&M_SWH);

            M_MWD=ExternalData(&M_WW3A_elements_dataset, M_mesh, 1,false,time_init);
            M_external_data.push_back(&M_MWD);

            M_FP=ExternalData(&M_WW3A_elements_dataset, M_mesh, 2,false,time_init);
            M_external_data.push_back(&M_FP);

            wim_forcing_options = M_WW3A_elements_dataset.grid.waveOptions;
            wim_ideal_forcing   = false;

            break;

        case setup::WaveType::ERAI_WAVES_1DEG:

            M_SWH=ExternalData(&M_ERAIW_1DEG_elements_dataset, M_mesh, 0,false,time_init);
            M_external_data.push_back(&M_SWH);

            M_MWD=ExternalData(&M_ERAIW_1DEG_elements_dataset, M_mesh, 1,false,time_init);
            M_external_data.push_back(&M_MWD);

            M_FP=ExternalData(&M_ERAIW_1DEG_elements_dataset, M_mesh, 2,false,time_init);
            M_external_data.push_back(&M_FP);

            wim_forcing_options = M_ERAIW_1DEG_elements_dataset.grid.waveOptions;
            wim_ideal_forcing   = false;

            break;

        default:
            std::cout << "invalid wave forcing"<<"\n";
            throw std::logic_error("invalid wave forcing");
    }
}
#endif

void
FiniteElement::initSlabOcean()
{
    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
            std::fill(M_sst.begin(), M_sst.end(), -1.8);
            std::fill(M_sss.begin(), M_sss.end(),  1.8/physical::mu);
            break;
        case setup::OceanType::TOPAZR:
        case setup::OceanType::TOPAZF:
            for ( int i=0; i<M_num_elements; ++i)
            {
                // Make sure the erroneous salinity and temperature don't screw up the initialisation too badly
                // This can still be done much better!
                M_sss[i] = std::max(physical::si, M_ocean_salt[i]);
                M_sst[i] = std::max(-M_sss[i]*physical::mu, M_ocean_temp[i]);
            }

            break;
        default:
            std::cout << "invalid ocean initialisation"<<"\n";
            throw std::logic_error("invalid ocean forcing");
    }
#if 1
    // setting SST to the freezing point for the part cover by sea ice
    for ( int i=0; i<M_num_elements; ++i)
    {
        if(M_conc[i]>0.)
            M_sst[i] = -M_sss[i]*physical::mu;//*M_conc[i]+M_ocean_temp[i]*(1.-M_conc[i]);
    }
#endif
}

void
FiniteElement::initIce()
{
    switch (M_ice_type)
    {
        case setup::IceType::CONSTANT:
            this->constantIce();
            break;
        case setup::IceType::CONSTANT_PARTIAL:
            this->constantIce();
            break;
        case setup::IceType::TARGET:
            this->targetIce();
            break;
        case setup::IceType::TOPAZ4:
            this->topazIce();
            break;
        case setup::IceType::TOPAZ4F:
            this->topazForecastIce();
            break;
        case setup::IceType::TOPAZ4FAMSR2:
            this->topazForecastAmsr2Ice();
            break;
        case setup::IceType::TOPAZ4FAMSR2OSISAF:
            this->topazForecastAmsr2OsisafIce();
            break;
        case setup::IceType::PIOMAS:
            this->piomasIce();
            break;
        case setup::IceType::AMSRE:
            this->topazAmsreIce();
            break;
        case setup::IceType::AMSR2:
            this->topazAmsr2Ice();
            break;
        case setup::IceType::CS2_SMOS:
            this->cs2SmosIce();
            break;

        default:
            std::cout << "invalid initialization of the ice"<<"\n";
            throw std::logic_error("invalid initialization of the ice");
    }

    // It's nice to initialise the ice temperature (especially for Winton)
    for ( int i=0; i<M_num_elements; i++ )
        if ( M_snow_thick[i] > 0. )
            M_tice[0][i] = std::min(0., M_tair[i]);
        else
            M_tice[0][i] = std::min(-physical::mu*M_sss[i], M_tair[i]);

    if ( M_thermo_type == setup::ThermoType::WINTON )
    {
        for ( int i=0; i<M_num_elements; i++ )
        {
            double Tbot  = -physical::mu*M_sss[i];
            if ( M_thick[i] > 0. )
            {
                // Just a linear interpolation between bottom and snow-ice interface (i.e. a zero layer model)
                double deltaT = (Tbot - M_tice[0][i] ) / ( 1. + physical::ki*M_snow_thick[i]/(physical::ks*M_thick[i]) );
                double slope = -deltaT/M_thick[i];
                M_tice[1][i] = Tbot + 3*slope*M_thick[i]/4;
                M_tice[2][i] = Tbot +   slope*M_thick[i]/4;
            } else {
                M_tice[1][i] = Tbot;
                M_tice[2][i] = Tbot;
            }
        }
    }
}

void
FiniteElement::constantIce()
{
	LOG(DEBUG) <<"Constant Ice\n";
    std::fill(M_conc.begin(), M_conc.end(), vm["simul.init_concentration"].as<double>());
    std::fill(M_thick.begin(), M_thick.end(), vm["simul.init_thickness"].as<double>());
    std::fill(M_snow_thick.begin(), M_snow_thick.end(), vm["simul.init_snow_thickness"].as<double>());
    std::fill(M_damage.begin(), M_damage.end(), 0.);

    if (M_ice_type==setup::IceType::CONSTANT_PARTIAL)
    {
        auto Bx = M_mesh.coordX();//xmin,xmax from nodes
        double xmin = *std::min_element(Bx.begin(),Bx.end());
        double xmax = *std::max_element(Bx.begin(),Bx.end());
        double xedge = xmin + 0.3*(xmax-xmin);

        std::cout<<"In constantIce (partial cover)\n";
        std::cout<<"M_ice_type "<< (int)M_ice_type<<"\n";
        std::cout<<"Min conc = "<< *std::min_element(M_conc.begin(),M_conc.end()) <<"\n";
        std::cout<<"Max conc = "<< *std::max_element(M_conc.begin(),M_conc.end()) <<"\n";
        std::cout<<"Min thick = "<< *std::min_element(M_thick.begin(),M_thick.end()) <<"\n";
        std::cout<<"Max thick = "<< *std::max_element(M_thick.begin(),M_thick.end()) <<"\n";
        std::cout<<"xmin="<<xmin<<"\n";
        std::cout<<"xmax="<<xmax<<"\n";
        std::cout<<"xedge="<<xedge<<"\n";

        Bx = M_mesh.bcoordX();//set conc, etc on elements
        for (int i=0; i<M_conc.size(); ++i)
        {
            if (Bx[i] < xedge)
            {
                M_conc[i]       = 0.;
                M_thick[i]      = 0.;
                M_snow_thick[i] = 0.;
            }
        }
        std::cout<<"New min conc = "<< *std::min_element(M_conc.begin(),M_conc.end()) <<"\n";
        std::cout<<"New max conc = "<< *std::max_element(M_conc.begin(),M_conc.end()) <<"\n";
        std::cout<<"New min thick = "<< *std::min_element(M_thick.begin(),M_thick.end()) <<"\n";
        std::cout<<"New max thick = "<< *std::max_element(M_thick.begin(),M_thick.end()) <<"\n";
        //std::abort();
    }//partial ice cover

}

void
FiniteElement::targetIce()
{
    double y_max=300000.;
    double x_max=350000.;
    double x_min=200000.;

	double tmp_var;

    auto RX = M_mesh.bcoordX();
    auto RY = M_mesh.bcoordY();

    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_var = (RY[i]<=y_max)*(RX[i]<=x_max)*(RX[i]>=x_min);

        std::cout<<"RX: "<< RX[i] << "RY: "<< RY[i] << "tmp_var: " << tmp_var << "\n";

        M_conc[i]  = vm["simul.init_concentration"].as<double>()*tmp_var;
		M_thick[i] = vm["simul.init_thickness"].as<double>()*tmp_var;
		M_snow_thick[i] = vm["simul.init_snow_thickness"].as<double>()*tmp_var;
        M_damage[i]=0.;

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
            M_damage[i]=1.;
        }
        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
            M_damage[i]=1.;
        }
    }
}
void
FiniteElement::topazIce()
{
    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		tmp_var=std::min(1.,M_init_conc[i]);
		M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

		M_damage[i]=0.;
	}
}
void
FiniteElement::topazForecastIce()
{
    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		tmp_var=std::min(1.,M_init_conc[i]);
		M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

		M_damage[i]=0.;
	}
}

void
FiniteElement::topazForecastAmsr2Ice()
{
    double real_thickness, init_conc_tmp;

    external_data M_conc_amsr2=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
    M_conc_amsr2.check_and_reload(M_mesh,time_init);

    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
        double uncertainty;
        if(M_conc_amsr2[i]<0.1)
            uncertainty=0.1;
		else
            uncertainty=0.05;

        double diff_mod_obs=M_conc_amsr2[i]-M_init_conc[i];
        if(std::abs(diff_mod_obs)>=uncertainty)
            M_conc[i] = std::min(1.,M_conc_amsr2[i]-(diff_mod_obs)/std::abs(diff_mod_obs)*uncertainty/2.);
        else
            M_conc[i] = std::min(1.,M_init_conc[i]);

        // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_conc[i];
		init_conc_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

        // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
        if(init_conc_tmp>0.05)
        {
            real_thickness=M_thick[i]/init_conc_tmp;
            M_thick[i]=real_thickness*M_conc[i];
        }

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=vm["simul.min_c"].as<double>())
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=vm["simul.min_h"].as<double>())
        {
            M_thick[i]=0.;
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<0.1*M_conc[i])
            M_thick[i]=0.1*M_conc[i];

		M_damage[i]=0.;
	}
}
void
FiniteElement::topazForecastAmsr2OsisafIce()
{
    double real_thickness, init_conc_tmp;

    external_data M_conc_osisaf=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,0,false,time_init);
    M_conc_osisaf.check_and_reload(M_mesh,time_init);

    external_data M_confidence_osisaf=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,1,false,time_init);
    M_confidence_osisaf.check_and_reload(M_mesh,time_init);

    external_data M_conc_amsr2=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
    M_conc_amsr2.check_and_reload(M_mesh,time_init);

    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
        double weight;
        if(M_confidence_osisaf[i]>4.9)
            M_conc[i] = (std::min(1.,M_conc_osisaf[i])+M_init_conc[i])/2.;
        else
        {
            if(M_conc_amsr2[i]<0.8)
                weight=1.-std::pow(1.-M_conc_amsr2[i],2.);
            else
                weight=1.-std::pow(1.-M_conc_amsr2[i],2.);
            M_conc[i] = std::min(1.,M_conc_amsr2[i])*(weight)+M_init_conc[i]*(1.-weight);
        }
        // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_conc[i];
		init_conc_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

        // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
        if(init_conc_tmp>0.05)
        {
            real_thickness=M_thick[i]/init_conc_tmp;
            M_thick[i]=real_thickness*M_conc[i];
        }

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_thick[i]=0.;
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

		M_damage[i]=0.;
	}
}
void
FiniteElement::piomasIce()
{
    external_data M_init_conc=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,0,false,time_init);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,1,false,time_init);
    M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,2,false,time_init);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    for (int i=0; i<M_num_elements; ++i)
    {
		M_conc[i] = std::min(1.,M_init_conc[i]);
		M_thick[i] = M_init_thick[i];
        M_snow_thick[i] = M_init_snow_thick[i];

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

		M_damage[i]=0.;
	}
}
void
FiniteElement::topazAmsreIce()
{
    double real_thickness, init_conc_tmp;

    external_data M_conc_amsre=ExternalData(&M_ice_amsre_elements_dataset,M_mesh,0,false,time_init);
    M_conc_amsre.check_and_reload(M_mesh,time_init);

    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		M_conc[i] = std::min(1.,M_conc_amsre[i]);

        // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_conc[i];
		init_conc_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

        // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
        if(init_conc_tmp>0.05)
        {
            real_thickness=M_thick[i]/init_conc_tmp;
            M_thick[i]=real_thickness*M_conc[i];
        }

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_thick[i]=0.;
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

		M_damage[i]=0.;
	}
}

void
FiniteElement::topazAmsr2Ice()
{
    double real_thickness, init_conc_tmp;

    external_data M_conc_amsr2=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
    M_conc_amsr2.check_and_reload(M_mesh,time_init);

    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		M_conc[i] = std::min(1.,M_conc_amsr2[i]);

        // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_conc[i];
		init_conc_tmp = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;
		tmp_var=M_init_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.;

        // Use 0.05 to get rid of slight inconsistencies in the TOPAZ output.
        if(init_conc_tmp>0.05)
        {
            real_thickness=M_thick[i]/init_conc_tmp;
            M_thick[i]=real_thickness*M_conc[i];
        }

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_thick[i]=0.;
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

		M_damage[i]=0.;
	}
}
void
FiniteElement::cs2SmosIce()
{
    external_data M_init_conc=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,0,false,time_init);
    M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,1,false,time_init);
    M_init_thick.check_and_reload(M_mesh,time_init);

    boost::gregorian::date dt = Nextsim::parse_date(time_init);
    int month_id=dt.month().as_number(); // 1 for January, 2 for February, and so on. This will be used to compute the snow from Warren climatology

    std::cout << "month_id: " << month_id <<"\n";

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);
    M_init_snow_thick.check_and_reload(M_mesh,time_init);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		tmp_var=std::min(1.,M_init_conc[i]);
		M_conc[i] = tmp_var;
		tmp_var=M_init_thick[i];
		M_thick[i] = tmp_var ;
		tmp_var=M_init_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.

        //if either c or h equal zero, we set the others to zero as well
        if(M_conc[i]<=0.)
        {
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        if(M_thick[i]<=0.)
        {
            M_conc[i]=0.;
            M_snow_thick[i]=0.;
        }

		M_damage[i]=0.;
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
        case setup::DrifterType::NONE:
            break;

        case setup::DrifterType::EQUALLYSPACED:
            this->equallySpacedDrifter();
            break;

        case setup::DrifterType::IABP:
            this->initIABPDrifter();
            break;

        default:
            std::cout << "invalid initialization of drifter"<<"\n";
            throw std::logic_error("invalid initialization of drifter");
    }
}

void
FiniteElement::coriolis()
{
    // Interpolation of the latitude
    std::vector<double> lat = M_mesh.meanLat();

    for (int i=0; i<M_fcor.size(); ++i)
    {
        M_fcor[i] = 2*(vm["simul.omega"].as<double>())*std::sin(lat[i]*PI/180.);
    }
}

void
FiniteElement::bathymetry()//(double const& u, double const& v)
{
    switch (M_bathymetry_type)
    {
        case setup::BathymetryType::CONSTANT:
            M_element_depth=ExternalData(vm["simul.constant_bathymetry"].as<double>());
            M_external_data.push_back(&M_element_depth);
            break;
        case setup::BathymetryType::ETOPO:
            M_element_depth=ExternalData(&M_bathymetry_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_element_depth);
            break;
        default:
            std::cout << "invalid bathymetry"<<"\n";
            throw std::logic_error("invalid bathymetry");
    }
}

#if defined (WAVES)
void
FiniteElement::initNFloes()
{
    M_nfloes.resize(M_num_elements);

    for (int i=0; i<M_num_elements; ++i)
    {
        M_nfloes[i] = M_conc[i]/std::pow(vm["wim.dfloepackinit"].as<double>(),2.);
    }
    std::cout<<"initNfloes:\n";
    std::cout<<"init dfloe "<<vm["wim.dfloepackinit"].as<double>()<<"\n";
    std::cout<<"Min Nfloes = "<<*std::min_element(M_nfloes.begin(),M_nfloes.end())<<"\n";
    std::cout<<"Max Nfloes = "<<*std::max_element(M_nfloes.begin(),M_nfloes.end())<<"\n";
}
#endif


void
FiniteElement::nodesToElements(double const* depth, std::vector<double>& v)
{
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

// A simple function to output the drifters in the model, IABP or otherwise
// The output could well be prettier!
void
FiniteElement::outputDrifter(std::fstream &drifters_out)
{
    // Initialize the map
    mapx_class *map;
    std::string configfile = (boost::format( "%1%/%2%/%3%" )
                              % Environment::nextsimDir().string()
                              % "data"
                              % vm["simul.proj_filename"].as<std::string>()
                              ).str();

    std::vector<char> str(configfile.begin(), configfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    // Assemble the coordinates from the unordered_map
    std::vector<double> drifter_X(M_drifter.size());
    std::vector<double> drifter_Y(M_drifter.size());
    int j=0;
    for ( auto it = M_drifter.begin(); it != M_drifter.end(); ++it )
    {
        drifter_X[j] = it->second[0];
        drifter_Y[j] = it->second[1];
        ++j;
    }

    // Interpolate the velocity and concentration onto the drifter positions
    int nb_var=2;
    std::vector<double> interp_drifter_in(nb_var*M_mesh.numNodes());

    // Interpolate the velocity
    for (int i=0; i<M_mesh.numNodes(); ++i)
    {
        interp_drifter_in[nb_var*i]   = M_UM[i];
        interp_drifter_in[nb_var*i+1] = M_UM[i+M_mesh.numNodes()];
    }

    double* interp_drifter_out;
    InterpFromMeshToMesh2dx(&interp_drifter_out,
        &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
        M_mesh.numNodes(),M_mesh.numTriangles(),
        &interp_drifter_in[0],
        M_mesh.numNodes(),nb_var,
        &drifter_X[0],&drifter_Y[0],M_drifter.size(),
        true, 0.);

    // Loop over the map and output
    j=0;
    boost::gregorian::date           date = Nextsim::parse_date( current_time );
    boost::posix_time::time_duration time = Nextsim::parse_time( current_time );
    for ( auto it = M_drifter.begin(); it != M_drifter.end(); ++it )
    {
        double lat, lon;
        inverse_mapx(map,it->second[0]+interp_drifter_out[j],it->second[1]+interp_drifter_out[j+M_drifter.size()],&lat,&lon);
        j++;

        drifters_out << setw(4) << date.year()
            << " " << setw( 2) << date.month().as_number()
            << " " << setw( 2) << date.day().as_number()
            << " " << setw( 2) << time.hours()
            << " " << setw(16) << it->first
            << fixed << setprecision(5)
            << " " << setw( 8) << lat
            << " " << setw(10) << lon << "\n";
    }

    xDelete<double>(interp_drifter_out);
    close_mapx(map);
}

// Add the buoys that have been put into the ice and remove dead ones
void
FiniteElement::updateIABPDrifter()
{
    // Initialize the map
    mapx_class *map;
    std::string configfile = (boost::format( "%1%/%2%/%3%" )
                              % Environment::nextsimDir().string()
                              % "data"
                              % vm["simul.proj_filename"].as<std::string>()
                              ).str();

    std::vector<char> str(configfile.begin(), configfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    // Read the current buoys from file
    int pos;    // To be able to rewind one line
    double time = current_time;
    std::vector<int> keepers;
    while ( time == current_time )
    {
        // Remember where we were
        pos = M_iabp_file.tellg();

        // Read the next line
        int year, month, day, hour, number;
        double lat, lon;
        M_iabp_file >> year >> month >> day >> hour >> number >> lat >> lon;
        std::string date = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);
        time = from_date_string(date) + hour/24.;

        // Remember which buoys are in the ice according to IABP
        keepers.push_back(number);

        // Project and add the buoy to the map if it's missing
        if ( M_drifter.count(number) == 0 )
        {
            double x, y;
            forward_mapx(map,lat,lon,&x,&y);
            M_drifter.emplace(number, std::array<double,2>{x, y});

        }
    }
    close_mapx(map);

    // Go through the M_drifter map and throw out the ones which IABP doesn't
    // report as being in the ice anymore
    for ( auto model = M_drifter.begin(); model != M_drifter.end(); /* ++model is not allowed here, because we use 'erase' */ )
    {
        bool keep = false;
        // Check against all the buoys we want to keep
        for ( auto obs = keepers.begin(); obs != keepers.end(); ++obs )
        {
            if ( model->first == *obs )
            {
                keep = true;
                break;
            }
        }

        // Delete or advance the iterator
        if ( ! keep )
            model = M_drifter.erase(model);
        else
            ++model;
    }
}

// Initialise by reading all the data from '79 up to time_init
// This is too slow, but only happens once so I won't try to optimise that for now
void
FiniteElement::initIABPDrifter()
{
    std::string filename = Environment::nextsimDir().string() + "/data/IABP_buoys.txt";
    M_iabp_file.open(filename, std::fstream::in);
    if ( ! M_iabp_file.good() )
        throw std::runtime_error("File not found: " + filename);

    int pos;    // To be able to rewind one line
    double time = from_date_string("1979-01-01");
    while ( time < time_init )
    {
        // Remember where we were
        pos = M_iabp_file.tellg();

        // Read the next line
        int year, month, day, hour, number;
        double lat, lon;
        M_iabp_file >> year >> month >> day >> hour >> number >> lat >> lon;
        std::string date = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);

        time = from_date_string(date) + hour/24.;
    }

    // We must rewind one line so that updateIABPDrifter works correctly
    M_iabp_file.seekg(pos);
}

void
FiniteElement::equallySpacedDrifter()
{
    LOG(DEBUG) << "equally spaced drifters not yet implemented" << "\n";
    throw std::logic_error("equally spaced drifters not yet implemented");
    /*if (M_drifter.size() ==0)
        M_drifter.resize(M_num_elements);

    std::fill(M_drifter.begin(), M_drifter.end(), 0.);*/
}

void
FiniteElement::importBamg(BamgMesh const* bamg_mesh)
{
    std::vector<point_type> mesh_nodes;
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

    int numVertices = 3;
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

        mesh_triangles.push_back(gmshElt);
    }

    LOG(DEBUG) <<"\n";
    LOG(DEBUG) <<"INFO: Previous  NumNodes     = "<< M_mesh.numNodes() <<"\n";
    LOG(DEBUG) <<"INFO: Previous  NumTriangles = "<< M_mesh.numTriangles() <<"\n";

    M_mesh_previous = M_mesh;
    //M_mesh = mesh_type(mesh_nodes,mesh_triangles);
    M_mesh.update(mesh_nodes,mesh_triangles);
    //M_mesh.writeTofile("out.msh");

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_num_nodes = M_mesh.numNodes();

    LOG(DEBUG) <<"\n";
    LOG(DEBUG) <<"INFO: Current  NumNodes      = "<< M_mesh.numNodes() <<"\n";
    LOG(DEBUG) <<"INFO: Current  NumTriangles  = "<< M_mesh.numTriangles() <<"\n";
    LOG(DEBUG) <<"\n";

    // std::cout<<"NodalConnectivitySize[0]= "<< bamg_mesh->NodalConnectivitySize[0] <<"\n";
    // std::cout<<"NodalConnectivitySize[1]= "<< bamg_mesh->NodalConnectivitySize[1] <<"\n";

    int Nd = bamg_mesh->NodalConnectivitySize[1];
    std::vector<int> dz;
    std::vector<int> ddz_j;
    std::vector<int> ddz_i;

    for (int i=0; i<bamgmesh->NodalConnectivitySize[0]; ++i)
    {
        int Nc = bamgmesh->NodalConnectivity[Nd*(i+1)-1];
        dz.push_back(2*(Nc+1));

        std::vector<int> local_ddz;
        local_ddz.push_back(i);

        for (int j=0; j<Nc; ++j)
        {
            local_ddz.push_back(bamgmesh->NodalConnectivity[Nd*i+j]-1);
        }
        std::sort(local_ddz.begin(),local_ddz.end());

        ddz_i.push_back(ddz_j.size());

        for (int const& k : local_ddz)
        {
            ddz_j.push_back(k);
        }

        for (int const& k : local_ddz)
        {
            ddz_j.push_back(k+M_num_nodes);
        }

#if 0
        for (int j=0; j<Nc; ++j)
        {
            std::cout<<"Connectivity["<< Nd*i+j <<"]= "<< bamgmesh->NodalConnectivity[Nd*i+j]-1 <<"\n";
        }
#endif
    }

    // auto mindzit = std::min_element(dz.begin(),dz.end());
    // auto maxdzit = std::max_element(dz.begin(),dz.end());

    // std::cout<<"************MINDZ= "<< *mindzit << " at "<< std::distance(dz.begin(), mindzit) << "\n";
    // std::cout<<"************MAXDZ= "<< *maxdzit << " at "<< std::distance(dz.begin(), maxdzit) <<"\n";

    auto dzu_count = dz.size();
    dz.resize(2*dzu_count);
    std::copy_n(dz.begin(), dzu_count, dz.begin() + dzu_count);

    int ddzi_count = ddz_i.size();
    ddz_i.resize(2*ddzi_count+1,0);
    ddz_i[ddzi_count] = ddz_i[ddzi_count-1] + dz[ddzi_count-1];

    for (int ll=0; ll<ddzi_count; ++ll)
    {
        ddz_i[ddzi_count+1+ll] = ddz_i[ddzi_count+ll] + dz[ll];
    }

    auto ddzj_count = ddz_j.size();
    ddz_j.resize(2*ddzj_count);
    std::copy_n(ddz_j.begin(), ddzj_count, ddz_j.begin() + ddzj_count);

    std::vector<double> ddz_data(ddz_j.size(),0.);
    M_graph = graph_type(dz,ddz_i,ddz_j,ddz_data);

    LOG(DEBUG) <<"\n";
    LOG(DEBUG) <<"GRAPHCSR INFO: MIN NZ (per row)      = "<< *std::min_element(dz.begin(),dz.end()) <<"\n";
    LOG(DEBUG) <<"GRAPHCSR INFO: MAX NZ (per row)      = "<< *std::max_element(dz.begin(),dz.end()) <<"\n";
    LOG(DEBUG) <<"GRAPHCSR INFO: NNZ (total)           = "<< ddz_j.size() <<"\n";
    LOG(DEBUG) <<"\n";


#if 0
    auto NNZ = M_graph.nNz();
    double minNNZ = *std::min_element(NNZ.begin(),NNZ.end());
    double maxNNZ = *std::max_element(NNZ.begin(),NNZ.end());

    std::cout<<"sizeNNZ= "<< NNZ.size() <<"\n";
    std::cout<<"minNNZ= "<< minNNZ <<"\n";
    std::cout<<"maxNNZ= "<< maxNNZ <<"\n";

    auto DDZI = M_graph.ia();
    auto DDZJ = M_graph.ja();

    int minDDZI = *std::min_element(DDZI.begin(),DDZI.end());
    int maxDDZI = *std::max_element(DDZI.begin(),DDZI.end());

    int minDDZJ = *std::min_element(DDZJ.begin(),DDZJ.end());
    int maxDDZJ = *std::max_element(DDZJ.begin(),DDZJ.end());

    std::cout<<"sizeDDZI= "<< DDZI.size() <<"\n";
    std::cout<<"minDDZI= "<< minDDZI <<"\n";
    std::cout<<"maxDDZI= "<< maxDDZI <<"\n";

    std::cout<<"ACCUMULATE= "<< std::accumulate(NNZ.begin(),NNZ.end(),0) <<"\n";

    std::cout<<"sizeDDZJ= "<< DDZJ.size() <<"\n";
    std::cout<<"minDDZJ= "<< minDDZJ <<"\n";
    std::cout<<"maxDDZJ= "<< maxDDZJ <<"\n";

    auto NNZ = M_graph.nNz();
    std::cout<<"ACCUMULATE= "<< std::accumulate(NNZ.begin(),NNZ.end(),0) <<"\n";
#endif
}

void
FiniteElement::exportResults(int step, bool export_mesh, bool export_fields)
{
    Exporter exporter("float");
    std::string fileout;


    if (export_mesh)
    {
        fileout = (boost::format( "%1%/mesh_%2%.bin" )
                   % M_export_path
                   % step ).str();

        LOG(INFO) <<"MESH BINARY: Exporter Filename= "<< fileout <<"\n";

		// move the mesh for the export
		M_mesh.move(M_UM,1.);

        std::fstream meshbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
        if ( ! meshbin.good() )
            throw std::runtime_error("Cannot write to file: " + fileout);

        exporter.writeMesh(meshbin, M_mesh);
        meshbin.close();

		// move it back after the export
		M_mesh.move(M_UM,-1.);

        fileout = (boost::format( "%1%/mesh_%2%.dat" )
               % M_export_path
               % step ).str();

        LOG(INFO) <<"RECORD MESH: Exporter Filename= "<< fileout <<"\n";

        std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
        if ( ! outrecord.good() )
            throw std::runtime_error("Cannot write to file: " + fileout);
        exporter.writeRecord(outrecord,"mesh");
        outrecord.close();
    }


    if (export_fields)
    {
        fileout = (boost::format( "%1%/field_%2%.bin" )
                   % M_export_path
                   % step ).str();

        LOG(INFO) <<"BINARY: Exporter Filename= "<< fileout <<"\n";

        std::fstream outbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
        if ( ! outbin.good() )
            throw std::runtime_error("Cannot write to file: " + fileout);
        std::vector<double> timevec(1);
        timevec[0] = current_time;
        exporter.writeField(outbin, timevec, "Time");
        exporter.writeField(outbin, M_surface, "Element_area");
        exporter.writeField(outbin, M_node_max_conc, "M_node_max_conc");
        exporter.writeField(outbin, M_VT, "M_VT");
        exporter.writeField(outbin, M_conc, "Concentration");
        exporter.writeField(outbin, M_thick, "Thickness");
        exporter.writeField(outbin, M_snow_thick, "Snow");
        exporter.writeField(outbin, M_damage, "Damage");
        int i=0;
        for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
        {
            exporter.writeField(outbin, *it, "Tice_"+std::to_string(i));
            i++;
        }
        exporter.writeField(outbin, M_sst, "SST");
        exporter.writeField(outbin, M_sss, "SSS");

        if(vm["simul.save_forcing_field"].as<bool>())
        {
            // Thermodynamic and dynamic forcing
            // Atmosphere
            exporter.writeField(outbin,M_wind.getVector(), "M_wind");         // Surface wind [m/s]
            exporter.writeField(outbin,M_tair.getVector(), "M_tair");         // 2 m temperature [C]
            exporter.writeField(outbin,M_mixrat.getVector(), "M_mixrat");       // Mixing ratio
            exporter.writeField(outbin,M_mslp.getVector(), "M_mslp");         // Atmospheric pressure [Pa]
            exporter.writeField(outbin,M_Qsw_in.getVector(), "M_Qsw_in");       // Incoming short-wave radiation [W/m2]
            exporter.writeField(outbin,M_Qlw_in.getVector(), "M_Qlw_in");       // Incoming long-wave radiation [W/m2]
            exporter.writeField(outbin,M_tcc.getVector(), "M_tcc");       // Incoming long-wave radiation [W/m2]
            exporter.writeField(outbin,M_precip.getVector(), "M_precip");       // Total precipitation [m]
            exporter.writeField(outbin,M_snowfr.getVector(), "M_snowfr");       // Fraction of precipitation that is snow
            exporter.writeField(outbin,M_dair.getVector(), "M_dair");         // 2 m dew point [C]

            // Ocean
            exporter.writeField(outbin,M_ocean.getVector(), "M_ocean");        // "Geostrophic" ocean currents [m/s]
            exporter.writeField(outbin,M_ssh.getVector(), "M_ssh");          // Sea surface elevation [m]

            exporter.writeField(outbin,M_ocean_temp.getVector(), "M_ocean_temp");   // Ocean temperature in top layer [C]
            exporter.writeField(outbin,M_ocean_salt.getVector(), "M_ocean_salt");   // Ocean salinity in top layer [C]
            exporter.writeField(outbin,M_mld.getVector(), "M_mld");          // Mixed-layer depth [m]

            // Bathymetry
            exporter.writeField(outbin,M_element_depth.getVector(), "M_element_depth");
        }

#if defined (WAVES)
        if (M_use_wim)
        {
            exporter.writeField(outbin, M_tau, "Stresses");
            exporter.writeField(outbin, M_nfloes, "Nfloes");
            exporter.writeField(outbin, M_dfloe, "Dfloe");
        }
#endif

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            exporter.writeField(outbin, M_h_thin, "Thin_ice");
            exporter.writeField(outbin, M_hs_thin, "Snow_thin_ice");
            exporter.writeField(outbin, M_tsurf_thin, "Tsurf_thin_ice");

            // Re-create the diagnostic variable 'concentration of thin ice'
            std::vector<double> conc_thin(M_mesh.numTriangles());
            double const rtanalpha = c_thin_max/h_thin_max;
            for ( int i=0; i<M_mesh.numTriangles(); ++i )
            {
                conc_thin[i] = std::min(std::min(M_h_thin[i]/physical::hmin,
                            std::sqrt(2.*M_h_thin[i]*rtanalpha)), 1.-M_conc[i]);
            }
            exporter.writeField(outbin, conc_thin, "Concentration_thin_ice");
        }

        // EXPORT sigma1 sigma2
        std::vector<double> sigma1(M_mesh.numTriangles());
        std::vector<double> sigma2(M_mesh.numTriangles());
        double sigma_s, sigma_n;
        std::vector<double> sigma_pred(3);

        for ( int i=0; i<M_mesh.numTriangles(); ++i )
        {
            sigma_pred[0]=M_sigma[3*i];
            sigma_pred[1]=M_sigma[3*i+1];
            sigma_pred[2]=M_sigma[3*i+2];

            sigma_s=std::hypot((sigma_pred[0]-sigma_pred[1])/2.,sigma_pred[2]);
            sigma_n= -         (sigma_pred[0]+sigma_pred[1])/2.;

            sigma1[i] = sigma_n+sigma_s;
            sigma2[i] = sigma_n-sigma_s;
        }
        exporter.writeField(outbin, sigma1, "Sigma1");
        exporter.writeField(outbin, sigma2, "Sigma2");

        outbin.close();

        fileout = (boost::format( "%1%/field_%2%.dat" )
                   % M_export_path
                   % step ).str();

        LOG(INFO) <<"RECORD FIELD: Exporter Filename= "<< fileout <<"\n";

        std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
        if ( ! outrecord.good() )
            throw std::runtime_error("Cannot write to file: " + fileout);
        exporter.writeRecord(outrecord);
        outrecord.close();
    }
}

#if defined (WAVES)
void
FiniteElement::nextsimToWim(bool step)
{
    if (M_run_wim)
    {

        chrono.restart();
        LOG(DEBUG) <<"Element Interp starts\n";
        // ELEMENT INTERPOLATION (c, h, Nfloes)

        int nb_var;
        if ( wim_ideal_forcing )
        {
            nb_var=3;
            M_SWH_grid.resize(0);
            M_FP_grid.resize(0);
            M_MWD_grid.resize(0);
        }
        else
            nb_var=6;

        std::vector<double> interp_elt_in(nb_var*M_num_elements);

        double* interp_elt_out;

        LOG(DEBUG) <<"ELEMENT: Interp starts\n";

        std::cout <<"ELEMENT: Interp (mesh->grid) starts\n";
        std::cout<<"Min conc = "<< *std::min_element(M_conc.begin(),M_conc.end()) <<"\n";
        std::cout<<"Max conc = "<< *std::max_element(M_conc.begin(),M_conc.end()) <<"\n";
        std::cout<<"Min thick = "<< *std::min_element(M_thick.begin(),M_thick.end()) <<"\n";
        std::cout<<"Max thick = "<< *std::max_element(M_thick.begin(),M_thick.end()) <<"\n";

        int tmp_nb_var=0;
        for (int i=0; i<M_num_elements; ++i)
        {
            tmp_nb_var=0;

            // concentration
            interp_elt_in[nb_var*i+tmp_nb_var] = M_conc[i];
            tmp_nb_var++;

            // thickness
            interp_elt_in[nb_var*i+tmp_nb_var] = M_thick[i];
            tmp_nb_var++;

            // Nfloes
            interp_elt_in[nb_var*i+tmp_nb_var] = M_nfloes[i];
            tmp_nb_var++;

            if ( !wim_ideal_forcing )
            {
                // significant wave height
                interp_elt_in[nb_var*i+tmp_nb_var] = M_SWH[i];
                tmp_nb_var++;

                // wave mean direction
                interp_elt_in[nb_var*i+tmp_nb_var] = M_MWD[i];
                tmp_nb_var++;

                // wave peak period
                if ( wim_forcing_options.use_mwp )
                    interp_elt_in[nb_var*i+tmp_nb_var] = M_FP[i];
                else
                {
                    interp_elt_in[nb_var*i+tmp_nb_var] = 0.;
                    if (M_FP[i]>0.)
                        interp_elt_in[nb_var*i+tmp_nb_var] = 1/M_FP[i];
                }
                tmp_nb_var++;
            }

            if(tmp_nb_var>nb_var)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }

        // interpolation from mesh to grid
        int nx = wim_grid.nx;
        int ny = wim_grid.ny;
        double dx = wim_grid.dx;
        double dy = wim_grid.dy;

        double xmin = (wim_grid.x)[0];
        double xmax = (wim_grid.x)[nx-1];
        double ymin = (wim_grid.y)[0];
        double ymax = (wim_grid.y)[ny-1];

        auto RX = M_mesh.coordX();
        auto RY = M_mesh.coordY();

        int num_elements_grid = nx*ny;

        // move the mesh for the interpolation on to the wim grid
		M_mesh.move(M_UM,1.);

        //std::cout<<"before interp mesh2grid\n";
        InterpFromMeshToGridx(interp_elt_out,
                              &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                              M_mesh.numNodes(),M_mesh.numTriangles(),
                              &interp_elt_in[0],
                              M_mesh.numTriangles(),nb_var,
                              xmin,ymax,
                              dx,dy,
                              nx,ny,
                              0.);

        if (1)
        {
            std::cout<<"nx = "<< nx <<"\n";
            std::cout<<"ny = "<< ny <<"\n";
            std::cout<<"dx = "<< dx <<"\n";
            std::cout<<"dy = "<< dy <<"\n";
            std::cout<<"xmin = "<< xmin <<"\n";
            std::cout<<"ymax = "<< ymax <<"\n";

            std::cout<<"MIN BOUND MESHX= "<< *std::min_element(RX.begin(),RX.end()) <<"\n";
            std::cout<<"MAX BOUND MESHX= "<< *std::max_element(RX.begin(),RX.end()) <<"\n";

            std::cout<<"MIN BOUND MESHY= "<< *std::min_element(RY.begin(),RY.end()) <<"\n";
            std::cout<<"MAX BOUND MESHY= "<< *std::max_element(RY.begin(),RY.end()) <<"\n";

            std::cout<<"------------------------------------------\n";

            std::cout<<"MIN BOUND GRIDX= "<< xmin <<"\n";
            std::cout<<"MAX BOUND GRIDX= "<< xmax <<"\n";

            std::cout<<"MIN BOUND GRIDY= "<< ymin <<"\n";
            std::cout<<"MAX BOUND GRIDY= "<< ymax <<"\n";
        }

        //std::cout<<"after interp mesh2grid\n";
        std::cout<<"ideal wave forcing: "<<wim_ideal_forcing<<"\n";

        // move back the mesh after the interpolation
		M_mesh.move(M_UM,-1.);

        if (!step)
        {
            M_icec_grid.assign(num_elements_grid,0.);
            M_iceh_grid.assign(num_elements_grid,0.);
            M_nfloes_grid.assign(num_elements_grid,0.);

            if ( !wim_ideal_forcing )
            {
                M_SWH_grid.assign(num_elements_grid,0.);
                M_MWD_grid.assign(num_elements_grid,0.);
                M_FP_grid.assign(num_elements_grid,0.);
            }

            M_taux_grid.assign(num_elements_grid,0.);
            M_tauy_grid.assign(num_elements_grid,0.);
        }

        for (int i=0; i<num_elements_grid; ++i)
        {
            tmp_nb_var=0;

            // concentration
            M_icec_grid[i] = interp_elt_out[nb_var*i+tmp_nb_var];
            tmp_nb_var++;

            // thickness
            M_iceh_grid[i] = interp_elt_out[nb_var*i+tmp_nb_var];
            tmp_nb_var++;

            // Nfloes
            M_nfloes_grid[i] = interp_elt_out[nb_var*i+tmp_nb_var];
            tmp_nb_var++;

            if ( !wim_ideal_forcing )
            {
                // significant wave height
                M_SWH_grid[i] = interp_elt_out[nb_var*i+tmp_nb_var];
                tmp_nb_var++;

                // wave mean direction
                M_MWD_grid[i] = interp_elt_out[nb_var*i+tmp_nb_var];
                tmp_nb_var++;

                // wave peak frequency
                M_FP_grid[i] = interp_elt_out[nb_var*i+tmp_nb_var];
                tmp_nb_var++;
            }

            if(tmp_nb_var>nb_var)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }

        //test interp
        std::cout<<"min conc   grid= "<< *std::min_element(M_icec_grid.begin(),M_icec_grid.end() )<<"\n";
        std::cout<<"max conc   grid= "<< *std::max_element(M_icec_grid.begin(),M_icec_grid.end() )<<"\n";
        std::cout<<"min thick  grid= "<< *std::min_element(M_iceh_grid.begin(),M_iceh_grid.end() )<<"\n";
        std::cout<<"max thick  grid= "<< *std::max_element(M_iceh_grid.begin(),M_iceh_grid.end() )<<"\n";
        std::cout<<"min Nfloes grid= "<< *std::min_element(M_nfloes_grid.begin(),M_nfloes_grid.end() )<<"\n";
        std::cout<<"max Nfloes grid= "<< *std::max_element(M_nfloes_grid.begin(),M_nfloes_grid.end() )<<"\n";

        if ( !wim_ideal_forcing )
        {
            std::cout<<"min SWH_grid= "<< *std::min_element(M_SWH_grid.begin(),M_SWH_grid.end() )<<"\n";
            std::cout<<"max SWH_grid= "<< *std::max_element(M_SWH_grid.begin(),M_SWH_grid.end() )<<"\n";
            std::cout<<"min MWD_grid= "<< *std::min_element(M_MWD_grid.begin(),M_MWD_grid.end() )<<"\n";
            std::cout<<"max MWD_grid= "<< *std::max_element(M_MWD_grid.begin(),M_MWD_grid.end() )<<"\n";
            std::cout<<"min FP_grid= "<< *std::min_element(M_FP_grid.begin(),M_FP_grid.end() )<<"\n";
            std::cout<<"max FP_grid= "<< *std::max_element(M_FP_grid.begin(),M_FP_grid.end() )<<"\n";
        }

        xDelete<double>(interp_elt_out);
    }
}//nextsimToWim
#endif

#if defined (WAVES)
void
FiniteElement::wimToNextsim(bool step)
{
    bool break_on_mesh =
        ( vm["nextwim.coupling-option"].template as<std::string>() == "breaking_on_mesh");
    std::cout<<"break_on_mesh="<<break_on_mesh<<"\n";

    if (M_run_wim)
    {

        // run wim
        std::vector<double> mf1, mf2, mf3;

        if ( break_on_mesh )
        {
            wim.setMesh(M_mesh.bcoordX(),//elements of mesh (x)
                        M_mesh.bcoordY(),//elements of mesh (y)
                        M_conc,
                        M_thick,
                        M_nfloes,
                        "km");
        }

        bool TEST_INTERP_MESH = false;
        //save mesh before entering WIM:
        // mesh file can then be copied inside WIM to correct path to allow plotting
        if (TEST_INTERP_MESH)
            this->exportResults(1001,true,false);

        wim.run(M_icec_grid, M_iceh_grid, M_nfloes_grid, mf1, mf2, mf3, step);
        // test this later
        //wim.run(M_icec_grid, M_iceh_grid, M_nfloes_grid, M_SWH_grid, M_MWD_grid, M_FP_grid, step);

        M_taux_grid = wim.getTaux();
        M_tauy_grid = wim.getTauy();

        if ( !break_on_mesh )
            M_nfloes_grid = wim.getNfloes();
        else
        {
            M_nfloes        = wim.getNfloesMesh();
            auto M_broken   = wim.getBrokenMesh();//TODO check this - maybe change damage later
            wim.clearMesh();
        }
    }

    if (!M_regrid)
        M_mesh.move(M_UM,1.);

    if (M_run_wim || M_regrid)
    {
        int nx = wim_grid.nx;
        int ny = wim_grid.ny;
        double dx = wim_grid.dx;
        double dy = wim_grid.dy;

        chrono.restart();
        LOG(DEBUG) <<"Nodal Interp starts\n";
        LOG(DEBUG) <<"NODAL: Interp starts\n";

        int num_elements_grid = nx*ny;

        // NODAL INTERPOLATION
        // - taux and tauy from waves
        if ((!step) || M_regrid)
           //initialisation or after regridding (need to reset sizes)
            M_tau.assign(2*M_num_nodes,0);

        // Set type of interpolation for grid-to-mesh
        // int interptype = TriangleInterpEnum;
        int interptype = BilinearInterpEnum;
        //int interptype = NearestInterpEnum;

        if (vm["nextwim.applywavestress"].as<bool>())
        {
           // can turn off effect of wave stress for testing
           // - if this is not done, we currently interp tau_x,tau_y each time step
           // TODO rethink this? (let them be advected? - this could lead to instability perhaps)
           int nb_var=2;
           std::vector<double> interp_in(nb_var*num_elements_grid,0.);
           double* interp_out;

           for (int i=0; i<num_elements_grid; ++i)
           {
               // interp taux,tauy to nodes of mesh
               interp_in[nb_var*i] = M_taux_grid[i];
               interp_in[nb_var*i+1] = M_tauy_grid[i];
           }

           InterpFromGridToMeshx(interp_out,            //data (out)
                                 &(wim_grid.x)[0], nx,  //x vector (source), length of x vector
                                 &(wim_grid.y)[0], ny,  //y vector (source), length of y vector
                                 &interp_in[0],         //data (in)
                                 ny,nx,                 //M,N: no of grid cells in y,x directions
                                                        //(to determine if corners or centers of grid have been input)
                                 nb_var,                //no of variables
                                 &M_mesh.coordX()[0],   // x vector (target)
                                 &M_mesh.coordY()[0],   // y vector (target)
                                 M_num_nodes,0.,        //target_size,default value
                                 interptype,            //interpolation type
                                 true                   //row_major (false = fortran/matlab order)
                                 );

            std::cout<<"\nINTERP GRID TO NODES\n";


            //assign taux,tauy
            for (int i=0; i<M_num_nodes; ++i)
            {
                // tau
                M_tau[i] = interp_out[nb_var*i];//tau_x
                M_tau[i+M_num_nodes] = interp_out[nb_var*i+1];//tau_y
            }
            std::cout<<"Min tau_x on grid = "<< *std::min_element(M_taux_grid.begin(),M_taux_grid.end()) <<"\n";
            std::cout<<"Max tau_x on grid = "<< *std::max_element(M_taux_grid.begin(),M_taux_grid.end()) <<"\n";
            std::cout<<"Min tau_y on grid = "<< *std::min_element(M_tauy_grid.begin(),M_tauy_grid.end()) <<"\n";
            std::cout<<"Max tau_y on grid = "<< *std::max_element(M_tauy_grid.begin(),M_tauy_grid.end()) <<"\n";
            std::cout<<"Min tau_x on mesh = "<< *std::min_element(M_tau.begin(),M_tau.begin()+M_num_nodes) <<"\n";
            std::cout<<"Max tau_x on mesh = "<< *std::max_element(M_tau.begin(),M_tau.begin()+M_num_nodes) <<"\n";
            std::cout<<"Min tau_y on mesh = "<< *std::min_element(M_tau.begin()+M_num_nodes,M_tau.end()) <<"\n";
            std::cout<<"Max tau_y on mesh = "<< *std::max_element(M_tau.begin()+M_num_nodes,M_tau.end()) <<"\n";

            xDelete<double>(interp_out);
            //std::abort();
        }//interp taux,tauy

        if ((M_run_wim) && (!break_on_mesh))
        {
            // interpolate nfloes to elements of mesh
            double* interp_out;
            InterpFromGridToMeshx(interp_out,         //data (out)
                                &(wim_grid.x)[0], nx, //x vector (source), length of x vector
                                &(wim_grid.y)[0], ny, //y vector (source), length of y vector
                                &M_nfloes_grid[0],    //data (in)
                                ny,nx,                //M,N no of grid cells in y,x directions
                                                      //(to determine if corners or centers of grid have been input)
                                1,                    //no of variables
                                &M_mesh.bcoordX()[0], // x vector (target)
                                &M_mesh.bcoordY()[0], // y vector (target)
                                M_num_elements,0.,    //target_size,default value
                                interptype,           //interpolation type
                                true                  //row_major (false = fortran/matlab order)
                                );
            std::cout<<"\nINTERP GRID TO ELEMENTS\n";

            M_nfloes.resize(M_num_elements);
            for (int i=0; i<M_num_elements; ++i)
            {
                // nfloes
                M_nfloes[i] = interp_out[i];
            }

            xDelete<double>(interp_out);
            std::cout<<"Min Nfloes on grid = "<< *std::min_element(M_nfloes_grid.begin(),M_nfloes_grid.end()) <<"\n";
            std::cout<<"Max Nfloes on grid = "<< *std::max_element(M_nfloes_grid.begin(),M_nfloes_grid.end()) <<"\n";
            std::cout<<"Min Nfloes on mesh = "<< *std::min_element(M_nfloes.begin(),M_nfloes.end()) <<"\n";
            std::cout<<"Max Nfloes on mesh = "<< *std::max_element(M_nfloes.begin(),M_nfloes.end()) <<"\n";
            //std::abort();
        }

    }

    if (!M_regrid)
        M_mesh.move(M_UM,-1.);

#if 0
    // set dfloe each time step (can be changed due to advection of nfloes by nextsim)
    M_dfloe.assign(M_num_elements,0.);

    for (int i=0; i<M_num_elements; ++i)
    {
        if (M_nfloes[i] > 0)
            M_dfloe[i] = std::sqrt(M_conc[i]/M_nfloes[i]);

        if (M_dfloe[i] > vm["wim.dfloepackthresh"].template as<double>())
            M_dfloe[i] = vm["wim.dfloepackinit"].template as<double>();

        if (M_conc[i] < vm["wim.cicemin"].template as<double>())
            M_dfloe[i] = 0.;
    }
#endif

    std::cout<<"Finished wimToNextsim";
}//wimToNextsim
#endif

std::string
FiniteElement::gitRevision()
{
    //std::string command = "git rev-parse HEAD";
    return this->system("git rev-parse HEAD");
}

std::string
FiniteElement::system(std::string const& command)
{
    char buffer[128];
    //std::string command = "git rev-parse HEAD";
    std::string result = "";
    std::shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe)
    {
        throw std::runtime_error("popen() failed!");
    }

    while (!feof(pipe.get()))
    {
        if (fgets(buffer, 128, pipe.get()) != NULL)
        {
            // remove newline from the buffer
            int len = strlen(buffer);
            if( buffer[len-1] == '\n' )
                buffer[len-1] = 0;

            result += buffer;
        }
    }

    // return the result of the command
    return result;
}

std::string
FiniteElement::getEnv(std::string const& envname)
{
    char* senv = ::getenv(envname.c_str());
    return std::string(senv);
}

void
FiniteElement::writeLogFile()
{

    std::string logfilename = "";
    if ((vm["simul.logfile"].as<std::string>()).empty())
    {
        logfilename = "nextsim.log";
    }
    else
    {
        logfilename = vm["simul.logfile"].as<std::string>();
    }

    std::string fileout = (boost::format( "%1%/%2%" )
               % M_export_path
               % logfilename ).str();

    std::fstream logfile(fileout, std::ios::out | std::ios::trunc);
    std::cout << "Writing log file " << fileout << "...\n";

    if (logfile.is_open())
    {
        logfile << "#----------Info\n";
        logfile << std::setw(40) << std::left << "Build date "  << current_time_local() <<"\n";
        logfile << std::setw(40) << std::left << "Git revision "  << gitRevision() <<"\n";

        logfile << "#----------Compilers\n";
        logfile << std::setw(40) << std::left << "C "  << system("which gcc") << " (version "<< system("gcc -dumpversion") << ")" <<"\n";
        logfile << std::setw(40) << std::left << "C++ "  << system("which g++") << " (version "<< system("g++ -dumpversion") << ")" <<"\n";

        logfile << "#----------Environment variables\n";
        logfile << std::setw(40) << std::left << "NEXTSIMDIR "  << getEnv("NEXTSIMDIR") <<"\n";
        logfile << std::setw(40) << std::left << "SIMDATADIR "  << getEnv("SIMDATADIR") <<"\n";
        logfile << std::setw(40) << std::left << "SIMFORECASTDIR "  << getEnv("SIMFORECASTDIR") <<"\n";
        logfile << std::setw(40) << std::left << "PETSC_DIR "  << getEnv("PETSC_DIR") <<"\n";
        logfile << std::setw(40) << std::left << "BOOST_DIR "  << getEnv("BOOST_DIR") <<"\n";
        logfile << std::setw(40) << std::left << "GMSH_DIR "  << getEnv("GMSH_DIR") <<"\n";
        logfile << std::setw(40) << std::left << "NETCDF_DIR "  << getEnv("NETCDF_DIR") <<"\n";
        logfile << std::setw(40) << std::left << "OPENMPI_LIB_DIR "  << getEnv("OPENMPI_LIB_DIR") <<"\n";
        logfile << std::setw(40) << std::left << "OPENMPI_INCLUDE_DIR "  << getEnv("OPENMPI_INCLUDE_DIR") <<"\n";

        logfile << "#----------Program options\n";

        for (po::variables_map::iterator it = vm.begin(); it != vm.end(); it++)
        {
            // ignore wim options if no coupling
#if !defined (WAVES)
            if ((it->first.find("nextwim.") != std::string::npos) || (it->first.find("wim.") != std::string::npos))
            {
                continue;
            }
#endif

            logfile << std::setw(40) << std::left << it->first;

#if 0
            if (((boost::any)it->second.value()).empty())
            {
                std::cout << "(empty)";
            }
            if (vm[it->first].defaulted() || it->second.defaulted()) {
                std::cout << "(default)";
            }
#endif

            bool is_char;
            try
            {
                boost::any_cast<const char *>(it->second.value());
                is_char = true;
            }
            catch (const boost::bad_any_cast &)
            {
                is_char = false;
            }

            bool is_str;
            try
            {
                boost::any_cast<std::string>(it->second.value());
                is_str = true;
            }
            catch (const boost::bad_any_cast &)
            {
                is_str = false;
            }

            if (((boost::any)it->second.value()).type() == typeid(int))
            {
                logfile << vm[it->first].as<int>() <<"\n";
            }
            else if (((boost::any)it->second.value()).type() == typeid(bool))
            {
                logfile << vm[it->first].as<bool>() <<"\n";
            }
            else if (((boost::any)it->second.value()).type() == typeid(double))
            {
                logfile << vm[it->first].as<double>() <<"\n";
            }
            else if (is_char)
            {
                logfile << vm[it->first].as<const char * >() <<"\n";
            }
            else if (is_str)
            {
                std::string temp = vm[it->first].as<std::string>();

                logfile << temp <<"\n";
#if 0
                if (temp.size())
                {
                    logfile << temp <<"\n";
                }
                else
                {
                    logfile << "true" <<"\n";
                }
#endif
            }
            else
            { // Assumes that the only remainder is vector<string>
                try
                {
                    std::vector<std::string> vect = vm[it->first].as<std::vector<std::string> >();
                    uint i = 0;
                    for (std::vector<std::string>::iterator oit=vect.begin(); oit != vect.end(); oit++, ++i)
                    {
                        //logfile << it->first << "[" << i << "]=" << (*oit) <<"\n";
                        if (i > 0)
                            logfile << std::setw(41) << std::right;

                        logfile << "[" << i << "]=" << (*oit) <<"\n";
                    }
                }
                catch (const boost::bad_any_cast &)
                {
                    std::cout << "UnknownType(" << ((boost::any)it->second.value()).type().name() << ")" <<"\n";
                }
            }
        }
    }
}

void
FiniteElement::clear()
{
    //delete[] M_topaz_grid.pfindex;

    delete bamgmesh;
    delete bamggeom;
    delete bamgopt;
    // We need to point these to NULL because 'delete bamgopt' clears the
    // memory they were pointing to before
    bamgopt_previous->hminVertices      = NULL;
    bamgopt_previous->hmaxVertices      = NULL;

    delete bamgmesh_previous;
    delete bamggeom_previous;
    delete bamgopt_previous;

    M_matrix->clear();
    M_vector->clear();
    M_solution->clear();
    M_solver->clear();
}
} // Nextsim
