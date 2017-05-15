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

    // setup the stereographic projection
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

    // set M_flag_fix to its correct value when PhysicalNames section is present in the msh file (version 2.2)
    if (!(M_mesh.markerNames()).empty())
    {
        LOG(DEBUG) <<"M_flag_fix before being changed was: " << M_flag_fix << "\n";
        // get the id associated to the physical name "coast" and assign it to M_flag_fix
        M_flag_fix = M_mesh.markerNames().find("coast")->second[0];
        LOG(DEBUG) <<"M_flag_fix changed to: " << M_flag_fix << "\n";
    }

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

    // We mask out the boundary nodes
    M_mask.assign(bamgmesh->VerticesSize[0],false) ;
    for (int vert=0; vert<bamgmesh->VerticesOnGeomVertexSize[0]; ++vert)
        M_mask[bamgmesh->VerticesOnGeomVertex[2*vert]-1]=true; // The factor 2 is because VerticesOnGeomVertex has 2 dimensions in bamg

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
    exportInitMesh();

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

    M_nb_regrid = 0;

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

    M_sst.resize(M_num_elements);
    M_sss.resize(M_num_elements);

    M_UM.resize(2*M_num_nodes,0.);
    M_UT.resize(2*M_num_nodes,0.);

    M_h_thin.assign(M_num_elements,0.);
    M_conc_thin.assign(M_num_elements,0.);
    M_hs_thin.assign(M_num_elements,0.);
    M_tsurf_thin.assign(M_num_elements,0.);
    
    M_sigma.resize(3*M_num_elements,0.);

    M_random_number.resize(M_num_elements);
    for (int i=0; i<M_random_number.size(); ++i)
        M_random_number[i] = static_cast <double> (std::rand()) / static_cast <double> (RAND_MAX);

    M_conc.resize(M_num_elements);
    M_thick.resize(M_num_elements);
    M_damage.resize(M_num_elements);
    M_ridge_ratio.assign(M_num_elements,0.);
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
    LOG(DEBUG) << "initSlabOcean\n";
    this->initSlabOcean();

    LOG(DEBUG) << "initIce\n";
    this->initIce();

    LOG(DEBUG) << "initDrifter\n";
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

        case setup::AtmosphereType::CFSR:
            M_atmosphere_nodes_dataset=DataSet("cfsr_nodes",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("cfsr_elements",M_num_elements);
            break;

        case setup::AtmosphereType::CFSR_HI:
            M_atmosphere_nodes_dataset=DataSet("cfsr_nodes_hi",M_num_nodes);
            M_atmosphere_elements_dataset=DataSet("cfsr_elements",M_num_elements);
            break;

        default:        std::cout << "invalid wind forcing"<<"\n";throw std::logic_error("invalid wind forcing");
    }

    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
        break;

        case setup::OceanType::TOPAZR:
        case setup::OceanType::TOPAZR_atrest:
            M_ocean_nodes_dataset=DataSet("topaz_nodes",M_num_nodes);
            M_ocean_elements_dataset=DataSet("topaz_elements",M_num_elements);
            break;

        case setup::OceanType::TOPAZF:
            M_ocean_nodes_dataset=DataSet("topaz_forecast_nodes",M_num_nodes);
            M_ocean_elements_dataset=DataSet("topaz_forecast_elements",M_num_elements);
            break;

        default:        std::cout << "invalid ocean forcing"<<"\n";throw std::logic_error("invalid ocean forcing");
    }
#if defined (WAVES)
    if (M_use_wim)
    {
        switch (M_wave_type)
        {
            case setup::WaveType::SET_IN_WIM:
                break;

            case setup::WaveType::CONSTANT:
                break;

            case setup::WaveType::CONSTANT_PARTIAL:
                break;

            case setup::WaveType::WW3A:
                M_wave_elements_dataset=DataSet("ww3a_elements",M_num_elements);
                break;

            case setup::WaveType::ERAI_WAVES_1DEG:
                M_wave_elements_dataset=DataSet("erai_waves_1deg_elements",M_num_elements);
                break;

            default:
                std::cout << "invalid wave forcing"<<"\n";
                throw std::logic_error("invalid wave forcing");
        }
    }
#endif

    M_ice_topaz_elements_dataset=DataSet("ice_topaz_elements",M_num_elements);

    M_ice_icesat_elements_dataset=DataSet("ice_icesat_elements",M_num_elements);

    M_ice_piomas_elements_dataset=DataSet("ice_piomas_elements",M_num_elements);

    M_ice_amsre_elements_dataset=DataSet("ice_amsre_elements",M_num_elements);

    M_ice_osisaf_elements_dataset=DataSet("ice_osisaf_elements",M_num_elements);

    M_ice_osisaf_type_elements_dataset=DataSet("ice_osisaf_type_elements",M_num_elements);

    M_ice_amsr2_elements_dataset=DataSet("ice_amsr2_elements",M_num_elements);

    M_ice_cs2_smos_elements_dataset=DataSet("ice_cs2_smos_elements",M_num_elements);

    M_ice_smos_elements_dataset=DataSet("ice_smos_elements",M_num_elements);

    M_bathymetry_elements_dataset=DataSet("etopo_elements",M_num_elements);//M_num_nodes);


}//initDatasets

#if 0
std::vector<double>
FiniteElements::rotatedWimElementsX(double const& rotangle) const
{
    //get x coord of WIM centers (rotated)
    std::vector<double> x(num_elements_wim_grid);
    double cos_rotangle=std::cos(rotangle);
    double sin_rotangle=std::sin(rotangle);
    for (int i=0; i<num_elements_wim_grid; ++i)
    {
        x[i] = cos_rotangle*wim_grid.X[i] + sin_rotangle*wim_grid.Y[i];
    }

    return x;
}

std::vector<double>
FiniteElements::rotatedWimElementsY(double const& rotangle) const
{
    //get x coord of WIM centers (rotated)
    std::vector<double> y(num_elements_wim_grid);
    double cos_rotangle=std::cos(rotangle);
    double sin_rotangle=std::sin(rotangle);
    for (int i=0; i<num_elements_wim_grid; ++i)
    {
        y[i] = -sin_rotangle*wim_grid.X[i] + cos_rotangle*wim_grid.Y[i];
    }

    return y;
}
#endif

void
FiniteElement::checkReloadDatasets(external_data_vec const& ext_data_vec,
        double const& CRtime, std::string const& printout)
{
    std::cout<<"size of external data vector = "<<ext_data_vec.size()<<"\n";
    if ( ext_data_vec.size()==0 )
    {
        LOG(DEBUG) <<"check_and_reload ("<<printout<<"):\n";
        LOG(DEBUG) <<"nothing to do\n";
        return;
    }

    //loop over M_external_data and call check and reload for each:
    chrono.restart();
    LOG(DEBUG) <<"check_and_reload ("<<printout<<") starts\n";

    //don't rotate yet since rotation angle not always defined yet
    auto RX_nod = M_mesh.coordX ();//nodes
    auto RY_nod = M_mesh.coordY ();
    auto RX_el  = M_mesh.bcoordX();//elements
    auto RY_el  = M_mesh.bcoordY();

    for ( auto it = ext_data_vec.begin(); it != ext_data_vec.end(); ++it )
    {
        if ( (*it)->M_is_constant )
        {
            std::vector<double> RX,RY;//not needed for constant forcings
            (*it)->check_and_reload(RX,RY,CRtime);
        }
        else
        {
            //dataset & interpolation etc needed
            if ( (*it)->M_dataset->grid.target_location=="mesh_nodes" )
            {
                LOG(DEBUG)<<"in nodes: dataset = "<<(*it)->M_dataset->name<<"\n";
                (*it)->check_and_reload(RX_nod,RY_nod,CRtime);
            }
            else if ( (*it)->M_dataset->grid.target_location=="mesh_elements" )
            {
                LOG(DEBUG)<<"in elements: dataset = "<<(*it)->M_dataset->name<<"\n";
                (*it)->check_and_reload(RX_el,RY_el,CRtime);
            }
#if defined (WAVES)
            else if ( (*it)->M_dataset->grid.target_location=="wim_grid" )
            {
                //LOG(DEBUG)<<"in wim_grid: dataset = "<<(*it)->M_dataset->name<<"\n";
                std::cout<<"in wim_grid: dataset = "<<(*it)->M_dataset->name<<"\n";
                //interp to WIM grid
                (*it)->check_and_reload(wim_grid.X,wim_grid.Y,CRtime);
            }
#endif
            else
            {
                std::cout<<"Bad value for dataset->grid.target_location: "<<(*it)->M_dataset->grid.target_location<<"\n";
                std::cout<<"- set to \"mesh_nodes\", or \"mesh_elements\"\n";
#if defined (WAVES)
                std::cout<<"or \"wim_grid\"\n";
#endif
                std::abort();
            }
        }
    }

    RX_el.resize(0);
    RY_el.resize(0);
    RX_nod.resize(0);
    RY_nod.resize(0);
    LOG(DEBUG) <<"check_and_reload in "<< chrono.elapsed() <<"s\n";
}//checkReloadDatasets

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
    time_step = vm["simul.timestep"].as<double>();

    output_time_step =  (vm["simul.output_per_day"].as<int>()<0) ? time_step : days_in_sec/vm["simul.output_per_day"].as<int>();
    mooring_output_time_step =  vm["simul.mooring_output_timestep"].as<double>()*days_in_sec;
    mooring_time_factor = time_step/mooring_output_time_step;

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

    ocean_turning_angle_rad = (PI/180.)*vm["simul.oceanic_turning_angle"].as<double>();
    ridging_exponent = vm["simul.ridging_exponent"].as<double>();

    quad_drag_coef_water = vm["simul.quad_drag_coef_water"].as<double>();

    basal_k2 = vm["simul.Lemieux_basal_k2"].as<double>();
    basal_u_0 = vm["simul.Lemieux_basal_u_0"].as<double>();
    basal_Cb = vm["simul.Lemieux_basal_Cb"].as<double>();

    time_relaxation_damage = vm["simul.time_relaxation_damage"].as<double>()*days_in_sec;
    deltaT_relaxation_damage = vm["simul.deltaT_relaxation_damage"].as<double>();

    h_thin_max = vm["simul.h_thin_max"].as<double>();
    h_thin_min = vm["simul.h_thin_min"].as<double>();

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
        ("ec_erai", setup::AtmosphereType::EC_ERAi)
        ("cfsr", setup::AtmosphereType::CFSR)
        ("cfsr_hi", setup::AtmosphereType::CFSR_HI);
    M_atmosphere_type = str2atmosphere.find(vm["setup.atmosphere-type"].as<std::string>())->second;

    switch(M_atmosphere_type){
        case setup::AtmosphereType::CONSTANT:   quad_drag_coef_air = vm["simul.ASR_quad_drag_coef_air"].as<double>(); break;
        case setup::AtmosphereType::ASR:        quad_drag_coef_air = vm["simul.ASR_quad_drag_coef_air"].as<double>(); break;
        case setup::AtmosphereType::CFSR_HI:
        case setup::AtmosphereType::CFSR:       quad_drag_coef_air = vm["simul.CFSR_quad_drag_coef_air"].as<double>(); break;
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
        ("topaz_atrest", setup::OceanType::TOPAZR_atrest)
        ("topaz_forecast", setup::OceanType::TOPAZF);
    M_ocean_type = str2ocean.find(vm["setup.ocean-type"].as<std::string>())->second;

    //std::cout<<"OCEANTYPE= "<< (int)M_ocean_type <<"\n";

    const boost::unordered_map<const std::string, setup::IceType> str2conc = boost::assign::map_list_of
        ("constant", setup::IceType::CONSTANT)
        ("constant_partial", setup::IceType::CONSTANT_PARTIAL)
        ("target", setup::IceType::TARGET)
        ("binary", setup::IceType::BINARY)
        ("topaz", setup::IceType::TOPAZ4)
        ("topaz_forecast", setup::IceType::TOPAZ4F)
        ("topaz_forecast_amsr2", setup::IceType::TOPAZ4FAMSR2)
        ("topaz_forecast_amsr2_osisaf", setup::IceType::TOPAZ4FAMSR2OSISAF)
        ("amsre", setup::IceType::AMSRE)
        ("amsr2", setup::IceType::AMSR2)
        ("osisaf", setup::IceType::OSISAF)
        ("piomas", setup::IceType::PIOMAS)
        ("cs2_smos", setup::IceType::CS2_SMOS)
        ("smos", setup::IceType::SMOS)
        ("topaz_osisaf_icesat", setup::IceType::TOPAZ4OSISAFICESAT);
    M_ice_type = str2conc.find(vm["setup.ice-type"].as<std::string>())->second;

#if defined (WAVES)
    if (M_use_wim)
    {
        const boost::unordered_map<const std::string, setup::WaveType> str2wave = boost::assign::map_list_of
            ("set_in_wim", setup::WaveType::SET_IN_WIM)
            ("constant", setup::WaveType::CONSTANT)
            ("constant_partial", setup::WaveType::CONSTANT_PARTIAL)
            ("ww3a", setup::WaveType::WW3A)
            ("eraiw_1deg", setup::WaveType::ERAI_WAVES_1DEG);
        M_wave_type = str2wave.find(vm["setup.wave-type"].as<std::string>())->second;
        std::cout<<"wave forcing type "<<vm["setup.wave-type"].as<std::string>()<<"\n";
        std::cout<<"wave forcing enum "<<(int)M_wave_type<<"\n";
    }
#endif

    const boost::unordered_map<const std::string, setup::BathymetryType> str2bathymetry = boost::assign::map_list_of
        ("constant", setup::BathymetryType::CONSTANT)
        ("etopo", setup::BathymetryType::ETOPO);
    M_bathymetry_type = str2bathymetry.find(vm["setup.bathymetry-type"].as<std::string>())->second;

    const boost::unordered_map<const std::string, setup::BasalStressType> str2basal_stress= boost::assign::map_list_of
        ("lemieux", setup::BasalStressType::LEMIEUX)
        ("bouillon", setup::BasalStressType::BOUILLON);
    M_basal_stress_type = str2basal_stress.find(vm["simul.basal_stress-type"].as<std::string>())->second;

    M_use_iabp_drifters=vm["simul.use_iabp_drifters"].as<bool>();
    M_equallyspaced_drifters_output_time_step=vm["simul.equallyspaced_drifters_output_time_step"].as<double>();
    M_rgps_drifters_output_time_step=vm["simul.rgps_drifters_output_time_step"].as<double>();
    M_use_osisaf_drifters=vm["simul.use_osisaf_drifters"].as<bool>();
    
    M_use_equallyspaced_drifters=false;
    M_use_rgps_drifters=false;
    
    if(M_equallyspaced_drifters_output_time_step>0.)
        M_use_equallyspaced_drifters=true;
    if(M_rgps_drifters_output_time_step>0.)
        M_use_rgps_drifters=true;

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

    M_moorings_snapshot =  vm["simul.mooring_snapshot"].as<bool>();

    const boost::unordered_map<const std::string, GridOutput::fileLength> str2mooringsfl = boost::assign::map_list_of
        ("inf", GridOutput::fileLength::inf)
        ("daily", GridOutput::fileLength::daily)
        ("weekly", GridOutput::fileLength::weekly)
        ("monthly", GridOutput::fileLength::monthly)
        ("yearly", GridOutput::fileLength::yearly);
    M_moorings_file_length = str2mooringsfl.find(vm["simul.mooring_file_length"].as<std::string>())->second;

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
    double area_init;

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
        area_init = this->jacobian(mesh.triangles()[cpt],mesh);
        area[cpt] = this->jacobian(movedmesh.triangles()[cpt],movedmesh);

        if(area_init*area[cpt]<=0.)
        {
            LOG(DEBUG) <<"FLIP DETECTED element:"<< cpt <<"\n";
        }
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
        hmin[i] = std::sqrt(2.)*std::sqrt(*std::min_element(measure.begin(),measure.end()))*0.8;
    }

    return hmin;
}

std::vector<double>
FiniteElement::hmaxVertices(mesh_type const& mesh, BamgMesh const* bamg_mesh) const
{
    std::vector<double> hmax(bamg_mesh->NodalElementConnectivitySize[0]);

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
        hmax[i] = std::sqrt(2.)*std::sqrt(*std::max_element(measure.begin(),measure.end()))*1.2;
    }

    return hmax;
}

std::vector<double>
FiniteElement::AllMinAngle(mesh_type const& mesh, std::vector<double> const& um, double factor) const
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

//#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        all_min_angle[cpt] = this->minAngles(movedmesh.triangles()[cpt],movedmesh);
    }
#endif

    return all_min_angle;
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
    double applied_displacement_factor=0.;
    double displacement_factor = 1.;

    bool flip = true;
    double minang;
    int substep = 0;

    M_nb_regrid++;

    std::vector<double> hmin_vertices_first;
    std::vector<double> hmax_vertices_first;

    mesh_adapt_step++;

    while(applied_displacement_factor<1.)
	{
		if(step)
		{
            flip = this->flip(M_mesh,M_UM,displacement_factor);

            minang = this->minAngle(M_mesh,M_UM,displacement_factor);

            substep++;
            while (flip || minang<(vm["simul.regrid_angle"].as<double>())/10.)
            {
                displacement_factor /= 2.;
                flip = this->flip(M_mesh,M_UM,displacement_factor);
                minang = this->minAngle(M_mesh,M_UM,displacement_factor);
            }
            LOG(DEBUG) <<"displacement_factor= "<< displacement_factor <<"\n";

            if(displacement_factor<1.)
    		{
    			LOG(WARNING) << "substeps will be needed for the remeshing!" <<"\n";
    			LOG(WARNING) << "Warning: It is probably due to very high ice speed, check your fields!\n";
    		}

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
                bamgopt->splitcorners      = 1;
				this->adaptMesh();
				bamgopt->KeepVertices=1;
                bamgopt->splitcorners      = 0;
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
            this->exportResults(400000+mesh_adapt_step+substep*100000,true,true,false);
		}

        chrono.restart();
        LOG(INFO) <<"AdaptMesh starts\n";
        applied_displacement_factor=applied_displacement_factor+displacement_factor;
        this->adaptMesh();
	    LOG(INFO) <<"AdaptMesh done in "<< chrono.elapsed() <<"s\n";

        if(step && (vm["simul.regrid_output_flag"].as<bool>()))
        {
            had_remeshed=true;
            this->exportResults(4000000+mesh_adapt_step+substep*100000,true,false,false);
		}

		if (step)
		{
#if 0
	        chrono.restart();
	        LOG(DEBUG) <<"Slab Interp starts\n";

            // We need to interpolate the slab first so that we have the right
            // SSS for setting the freezing temperature in M_tice[1] &
            // M_tice[2] for ice free elements
			// ELEMENT INTERPOLATION FOR SLAB OCEAN FROM OLD MESH ON ITS ORIGINAL POSITION
            int prv_num_elements = M_mesh_previous.numTriangles();
			int nb_var=2;

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
#endif
            chrono.restart();

            LOG(DEBUG) <<"Element Interp starts\n";

            // 1) collect the variables into a single structure
            int prv_num_elements = M_mesh_previous.numTriangles();
            
            double* interp_elt_in;
            int* interp_method;
            double* diffusivity_parameters;
            
            int nb_var = this->collectVariables(&interp_elt_in, &interp_method, &diffusivity_parameters, prv_num_elements);

            // 2) Interpolate
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

            double* interp_elt_out;
			InterpFromMeshToMesh2dCavities(&interp_elt_out,&interp_elt_in[0], &interp_method[0], nb_var,
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

            // 3) Change the size of the variables
			M_conc.assign(M_num_elements,0.);
			M_thick.assign(M_num_elements,0.);
			M_snow_thick.assign(M_num_elements,0.);
			M_sigma.assign(3*M_num_elements,0.);
			M_damage.assign(M_num_elements,0.);
            M_ridge_ratio.assign(M_num_elements,0.);
            
			M_random_number.resize(M_num_elements);
            
			M_sst.assign(M_num_elements,0.);
			M_sss.assign(M_num_elements,0.);

            for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
                it->assign(M_num_elements,0.);

            M_h_thin.assign(M_num_elements,0.);
            M_conc_thin.assign(M_num_elements,0.);
            M_hs_thin.assign(M_num_elements,0.);
            M_tsurf_thin.assign(M_num_elements,0.);
            

#if defined (WAVES)
            bool nfloes_interp = M_use_wim;
            if ( !(vm["nextwim.coupling-option"].template as<std::string>() == "breaking_on_mesh"))
                bool nfloes_interp = (M_use_wim && (!M_run_wim));
            if (nfloes_interp)
                M_nfloes.assign(M_num_elements,0.);
#endif
            // 4) redistribute the interpolated values
            this->redistributeVariables(&interp_elt_out[0],nb_var);

            // 5) cleaning
			xDelete<double>(interp_elt_out);
			xDelete<double>(interp_elt_in);
            xDelete<int>(interp_method);
            xDelete<double>(diffusivity_parameters);

			LOG(DEBUG) <<"ELEMENT: Interp done\n";
			LOG(DEBUG) <<"Element Interp done in "<< chrono.elapsed() <<"s\n";

			// NODAL INTERPOLATION
			nb_var=10;

			int prv_num_nodes = M_mesh_previous.numNodes();

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

				// UT
				interp_in[nb_var*i+8] = M_UT[i];
				interp_in[nb_var*i+9] = M_UT[i+prv_num_nodes];
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
			M_UT.assign(2*M_num_nodes,0.);

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

				// UT
				M_UT[i] = interp_out[nb_var*i+8];
				M_UT[i+M_num_nodes] = interp_out[nb_var*i+9];
			}

			xDelete<double>(interp_out);

			LOG(DEBUG) <<"NODAL: Interp done\n";
			LOG(DEBUG) <<"Nodal interp done in "<< chrono.elapsed() <<"s\n";

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


        M_fcor.resize(M_num_elements);
    }

    M_atmosphere_nodes_dataset.target_size=M_num_nodes;
    M_atmosphere_elements_dataset.target_size=M_num_elements;
    M_atmosphere_bis_elements_dataset.target_size=M_num_elements;
    M_ocean_nodes_dataset.target_size=M_num_nodes;
    M_ocean_elements_dataset.target_size=M_num_elements;

    M_ice_topaz_elements_dataset.target_size=M_num_elements;
    M_ice_icesat_elements_dataset.target_size=M_num_elements;
    M_ice_piomas_elements_dataset.target_size=M_num_elements;
    M_ice_amsre_elements_dataset.target_size=M_num_elements;
    M_ice_osisaf_elements_dataset.target_size=M_num_elements;
    M_ice_osisaf_type_elements_dataset.target_size=M_num_elements;
    M_ice_amsr2_elements_dataset.target_size=M_num_elements;
    M_ice_cs2_smos_elements_dataset.target_size=M_num_elements;
    M_ice_smos_elements_dataset.target_size=M_num_elements;
    M_bathymetry_elements_dataset.target_size=M_num_elements;
#if defined (WAVES)
    if (M_use_wim)
        M_wave_elements_dataset.target_size=M_num_elements;
#endif


    M_atmosphere_nodes_dataset.interpolated=false;
    M_atmosphere_elements_dataset.interpolated=false;
    M_atmosphere_bis_elements_dataset.interpolated=false;
    M_ocean_nodes_dataset.interpolated=false;
    M_ocean_elements_dataset.interpolated=false;

    M_ice_topaz_elements_dataset.interpolated=false;
    M_ice_icesat_elements_dataset.interpolated=false;
    M_ice_piomas_elements_dataset.interpolated=false;
    M_ice_amsre_elements_dataset.interpolated=false;
    M_ice_osisaf_elements_dataset.interpolated=false;
    M_ice_osisaf_type_elements_dataset.interpolated=false;
    M_ice_amsr2_elements_dataset.interpolated=false;
    M_ice_cs2_smos_elements_dataset.interpolated=false;
    M_ice_smos_elements_dataset.interpolated=false;
    M_bathymetry_elements_dataset.interpolated=false;
#if defined (WAVES)
    if (M_use_wim)
        M_wave_elements_dataset.interpolated=false;
#endif

    // for the parallel code, it will be necessary to add these lines
    // as the domain covered by the partitions changes at each remeshing/partitioning
#if 0
    M_atmosphere_nodes_dataset.grid.interpolated=false;
    M_atmosphere_elements_dataset.grid.interpolated=false;
    M_atmosphere_bis_elements_dataset.grid.interpolated=false;
    M_ocean_nodes_dataset.grid.interpolated=false;
    M_ocean_elements_dataset.grid.interpolated=false;
    M_ice_topaz_elements_dataset.grid.interpolated=false;
    M_ice_icesat_elements_dataset.grid.interpolated=false;
    M_ice_amsre_elements_dataset.grid.interpolated=false;
    M_ice_osisaf_elements_dataset.grid.interpolated=false;
    M_ice_osisaf_type_elements_dataset.grid.interpolated=false;
    M_ice_amsr2_elements_dataset.grid.interpolated=false;
    M_ice_cs2_smos_elements_dataset.grid.interpolated=false;
    M_ice_smos_elements_dataset.grid.interpolated=false;
    M_bathymetry_elements_dataset.grid.interpolated=false;

    M_atmosphere_nodes_dataset.grid.loaded=false;
    M_atmosphere_elements_dataset.grid.loaded=false;
    M_atmosphere_bis_elements_dataset.grid.loaded=false;
    M_ocean_nodes_dataset.grid.loaded=false;
    M_ocean_elements_dataset.grid.loaded=false;
    M_ice_topaz_elements_dataset.grid.loaded=false;
    M_ice_icesat_elements_dataset.grid.loaded=false;
    M_ice_amsre_elements_dataset.grid.loaded=false;
    M_ice_osisaf_elements_dataset.grid.loaded=false;
    M_ice_osisaf_type_elements_dataset.grid.loaded=false;
    M_ice_amsr2_elements_dataset.grid.loaded=false;
    M_ice_cs2_smos_elements_dataset.grid.loaded=false;
    M_ice_smos_elements_dataset.grid.loaded=false;
    M_bathymetry_elements_dataset.grid.loaded=false;
#endif

    M_Cohesion.resize(M_num_elements);
    M_Compressive_strength.resize(M_num_elements);
    M_time_relaxation_damage.resize(M_num_elements,time_relaxation_damage);
}

void
FiniteElement::redistributeVariables(double* interp_elt_out,int nb_var)
{
	for (int i=0; i<M_num_elements; ++i)
	{
		int tmp_nb_var=0;

		// concentration
		M_conc[i] = std::max(0., interp_elt_out[nb_var*i+tmp_nb_var]);
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

            // damage
		    M_damage[i] = std::max(0., std::min(1.,interp_elt_out[nb_var*i+tmp_nb_var]));
		    tmp_nb_var++;
            
            // damage
		    M_ridge_ratio[i] = std::max(0., std::min(1.,interp_elt_out[nb_var*i+tmp_nb_var]/M_thick[i]));
		    tmp_nb_var++;
		}
		else
		{
			// integrated_stress1
			M_sigma[3*i] = 0.;
			tmp_nb_var++;

			// integrated_stress2
			M_sigma[3*i+1] = 0.;
			tmp_nb_var++;

			// integrated_stress3
			M_sigma[3*i+2] = 0.;
			tmp_nb_var++;

            // damage
		    M_damage[i] = 1.;
		    tmp_nb_var++;
            
            // damage
		    M_ridge_ratio[i] = 0.;
		    tmp_nb_var++;
		}

		// random_number
		M_random_number[i] = interp_elt_out[nb_var*i+tmp_nb_var];
		//M_random_number[i] = std::max(0., std::min(1.,interp_elt_in[11*i+tmp_nb_var]));
		tmp_nb_var++;

		// SSS
		M_sss[i] = interp_elt_out[nb_var*i+tmp_nb_var];
		tmp_nb_var++;

		// SST
		M_sst[i] = interp_elt_out[nb_var*i+tmp_nb_var];
		tmp_nb_var++;

		// Ice temperature
        M_tice[0][i] = interp_elt_out[nb_var*i+tmp_nb_var];
        tmp_nb_var++;

        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            if(M_thick[i]>0.)
            {
                double tmp = interp_elt_out[nb_var*i+tmp_nb_var]/M_thick[i];
                M_tice[1][i] = 0.5*( tmp - std::sqrt(tmp*tmp + 4*physical::mu*physical::si*physical::Lf/physical::C) ); // (38) divided with volume with f1=1
                tmp_nb_var++;

                M_tice[2][i] = interp_elt_out[nb_var*i+tmp_nb_var]/M_thick[i]; // (40) divided with volume with f1=0
                tmp_nb_var++;
            }
            else
            {
                M_tice[1][i] = -physical::mu*M_sss[i];
                tmp_nb_var++;

                M_tice[2][i] = -physical::mu*M_sss[i];
                tmp_nb_var++;
            }
        }

		// thin ice thickness
		M_h_thin[i] = std::max(0., interp_elt_out[nb_var*i+tmp_nb_var]);
		tmp_nb_var++;

		// thin ice thickness
		M_conc_thin[i] = std::max(0., interp_elt_out[nb_var*i+tmp_nb_var]);
		tmp_nb_var++;

		// snow on thin ice
		M_hs_thin[i] = std::max(0., interp_elt_out[nb_var*i+tmp_nb_var]);
		tmp_nb_var++;

		// Ice surface temperature for thin ice
		M_tsurf_thin[i] = interp_elt_out[nb_var*i+tmp_nb_var];
		tmp_nb_var++;

#if defined (WAVES)
        bool nfloes_interp = M_use_wim;
        if ( !(vm["nextwim.coupling-option"].template as<std::string>() == "breaking_on_mesh"))
            bool nfloes_interp = (M_use_wim && (!M_run_wim));

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
}

void
FiniteElement::advect(double** interp_elt_out_ptr,double* interp_elt_in, int* interp_method,int nb_var)
{

	/*Initialize output*/
	double* interp_elt_out=NULL;

    interp_elt_out=xNew<double>(nb_var*M_num_elements);

    int thread_id;
    int total_threads;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    std::vector<double> UM_P = M_UM;

    int Nd = bamgmesh->NodalConnectivitySize[1];

    int ALE_smoothing_step_nb=vm["simul.ALE_smoothing_step_nb"].as<int>();
    // ALE_smoothing_step_nb==-2 is the case with no advection M_UM is not changed and then =0 and no fluxes are computed.
    // ALE_smoothing_step_nb==-1 is the diffusive eulerian case where M_UM is not changed and then =0.
    // ALE_smoothing_step_nb=0 is the purely Lagrangian case where M_UM is updated with M_VT
    // ALE_smoothing_step_nb>0 is the ALE case where M_UM is updated with a smoothed version of M_VT

    // increment M_UT that is used for the drifters
    for (int nd=0; nd<M_UM.size(); ++nd)
    {
        M_UT[nd] += time_step*M_VT[nd]; // Total displacement (for drifters)
    }

    if(ALE_smoothing_step_nb>=0)
    {
        std::vector<double> M_VT_smoothed = M_VT;
        std::vector<double> M_VT_tmp = M_VT_smoothed;

        for (int k=0; k<ALE_smoothing_step_nb; ++k)
        {
            M_VT_tmp=M_VT_smoothed;

#pragma omp parallel for num_threads(max_threads) private(thread_id)
            for (int i=0; i<M_num_nodes; ++i)
            {
                int Nc;
                double UM_x, UM_y;

                if(M_mask_dirichlet[i]==false)
                {
                    Nc = bamgmesh->NodalConnectivity[Nd*(i+1)-1];

                    UM_x=0.;
                    UM_y=0.;
                    for (int j=0; j<Nc; ++j)
                    {
                        UM_x += M_VT_tmp[bamgmesh->NodalConnectivity[Nd*i+j]-1]  ;
                        UM_y += M_VT_tmp[bamgmesh->NodalConnectivity[Nd*i+j]-1+M_num_nodes]  ;
                    }

                    M_VT_smoothed[i             ]=UM_x/Nc;
                    M_VT_smoothed[i+M_num_nodes ]=UM_y/Nc;
                }
            }
        }
        for (int nd=0; nd<M_UM.size(); ++nd)
        {
            M_UM[nd] += time_step*M_VT_smoothed[nd];
        }

        // set back the neumann nodes (open boundaries) at their position, the fluxes will be computed thanks to the convective velocity
        for (const int& nd : M_neumann_nodes)
        {
            M_UM[nd] = UM_P[nd];
        }
    }

    LOG(DEBUG) <<"VT MIN= "<< *std::min_element(M_VT.begin(),M_VT.end()) <<"\n";
    LOG(DEBUG) <<"VT MAX= "<< *std::max_element(M_VT.begin(),M_VT.end()) <<"\n";

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        /* some variables used for the advection*/
        double surface, surface_new;
        double integrated_variable;

        /* some variables used for the advection*/
        double x[3],y[3],x_new[3],y_new[3];
        int x_ind, y_ind, neighbour_int, vertex_1, vertex_2;
        double outer_fluxes_area[3], vector_edge[2], outer_vector[2], VC_middle[2], VC_x[3], VC_y[3];
        int fluxes_source_id[3];

        double neighbour_double;
        int other_vertex[3*2]={1,2 , 2,0 , 0,1};

        /*======================================================================
         * Update:
         * Ice and snow thickness, and concentration using a Lagrangian or an Eulerian scheme
         *======================================================================
         */

        /* convective velocity */
        for(int i=0;i<3;i++)
        {
            x_ind=(M_elements[cpt]).indices[i]-1;
            y_ind=(M_elements[cpt]).indices[i]-1+M_num_nodes;

            x[i] = M_nodes[(M_elements[cpt]).indices[i]-1].coords[0];
            y[i] = M_nodes[(M_elements[cpt]).indices[i]-1].coords[1];

            /* old and new positions of the mesh */
            x_new[i]=M_nodes[(M_elements[cpt]).indices[i]-1].coords[0]+M_UM[x_ind];
            y_new[i]=M_nodes[(M_elements[cpt]).indices[i]-1].coords[1]+M_UM[y_ind];
            x[i]    =M_nodes[(M_elements[cpt]).indices[i]-1].coords[0]+UM_P[x_ind];
            y[i]    =M_nodes[(M_elements[cpt]).indices[i]-1].coords[1]+UM_P[y_ind];

            if(ALE_smoothing_step_nb==-2)
            {
                VC_x[i] =0.;
                VC_y[i] =0.;
            }
            else
            {
                VC_x[i] =M_VT[x_ind]-(M_UM[x_ind]-UM_P[x_ind])/time_step;
                VC_y[i] =M_VT[y_ind]-(M_UM[y_ind]-UM_P[y_ind])/time_step;
            }
        }
        for(int i=0;i<3;i++)
        {
            outer_fluxes_area[i]=0;

            vertex_1=other_vertex[2*i  ];
            vertex_2=other_vertex[2*i+1];

            vector_edge[0]=x[vertex_2]-x[vertex_1];
            vector_edge[1]=y[vertex_2]-y[vertex_1];

            outer_vector[0]= vector_edge[1];
            outer_vector[1]=-vector_edge[0];

            VC_middle[0] = (VC_x[vertex_2]+VC_x[vertex_1])/2.;
            VC_middle[1] = (VC_y[vertex_2]+VC_y[vertex_1])/2.;

            outer_fluxes_area[i]=outer_vector[0]*VC_middle[0]+outer_vector[1]*VC_middle[1];


            if(outer_fluxes_area[i]>0)
            {
                surface = this->measure(M_elements[cpt],M_mesh, UM_P);
                outer_fluxes_area[i]=std::min(surface/time_step/3.,outer_fluxes_area[i]);
                fluxes_source_id[i]=cpt;
            }
            else
            {
			    neighbour_double=bamgmesh->ElementConnectivity[cpt*3+i];
                neighbour_int=(int) bamgmesh->ElementConnectivity[cpt*3+i];
			    if (!std::isnan(neighbour_double) && neighbour_int>0)
                {
                    surface = this->measure(M_elements[neighbour_int-1],M_mesh, UM_P);
                    outer_fluxes_area[i]=-std::min(surface/time_step/3.,-outer_fluxes_area[i]);
                    fluxes_source_id[i]=neighbour_int-1;
                }
                else // open boundary with incoming fluxes
                    fluxes_source_id[i]=cpt;
            }
        }


        surface = this->measure(M_elements[cpt],M_mesh, UM_P);
        surface_new = this->measure(M_elements[cpt],M_mesh,M_UM);
        M_surface[cpt] = surface_new;

        for(int j=0; j<nb_var; j++)
        {
            if(interp_method[j]==1)
            {
                integrated_variable=interp_elt_in[cpt*nb_var+j]*surface -
                    (interp_elt_in[fluxes_source_id[0]*nb_var+j]*outer_fluxes_area[0]  +
                    interp_elt_in[fluxes_source_id[1]*nb_var+j]*outer_fluxes_area[1]  +
                    interp_elt_in[fluxes_source_id[2]*nb_var+j]*outer_fluxes_area[2]  )*time_step;

                interp_elt_out[cpt*nb_var+j]    = integrated_variable/surface_new;
            }
            else
            {
                interp_elt_out[cpt*nb_var+j] = interp_elt_in[cpt*nb_var+j];
            }
        }
    }
	*interp_elt_out_ptr=interp_elt_out;
}


void
    FiniteElement::diffuse(double** interp_elt_out_ptr,double* interp_elt_in, double* diffusivity_parameters,int nb_var,double dx)
{

    /*Initialize output*/
    double* interp_elt_out=NULL;

    interp_elt_out=xNew<double>(nb_var*M_num_elements);

    int thread_id;
    int total_threads;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    std::vector<double> UM_P = M_UM;

    int Nd = bamgmesh->NodalConnectivitySize[1];

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        /* some variables used for the advection*/
        int neighbour_int;
        double fluxes_source[3];
        int fluxes_source_id;

        double neighbour_double;
        int other_vertex[3*2]={1,2 , 2,0 , 0,1};
      
        for(int j=0; j<nb_var; j++)
        {
            if(diffusivity_parameters[j]==0.)
                interp_elt_out[cpt*nb_var+j] = interp_elt_in[cpt*nb_var+j];
            else
            {
                
                for(int i=0;i<3;i++)
                {
                    neighbour_double=bamgmesh->ElementConnectivity[cpt*3+i];
                    neighbour_int=(int) bamgmesh->ElementConnectivity[cpt*3+i];
                    
                    if (!std::isnan(neighbour_double) && neighbour_int>0)
                    {
                        fluxes_source_id=neighbour_int-1;
                        fluxes_source[i]=diffusivity_parameters[j]*time_step/std::pow(dx,2.)*(interp_elt_in[fluxes_source_id*nb_var+j]-interp_elt_in[cpt*nb_var+j]);
                    }
                    else // no diffusion crosses open nor closed boundaries
                        fluxes_source[i]=0.;
                }
                
                interp_elt_out[cpt*nb_var+j] = interp_elt_in[cpt*nb_var+j] + fluxes_source[0] + fluxes_source[1] + fluxes_source[2];                
            }
        }
    }
    *interp_elt_out_ptr=interp_elt_out;
}


int
FiniteElement::collectVariables(double** interp_elt_in_ptr, int** interp_method_ptr, double** diffusivity_parameters_ptr, int prv_num_elements)
{
    // ELEMENT INTERPOLATION With Cavities
	int nb_var=15 + M_tice.size();

#if defined (WAVES)
    // coupling with wim
    // - only interpolate if not at a coupling time step
    // - else nfloes will just be overwritten with wimToNextsim()
    // - EXCEPT if doing breaking on mesh
    //   - then we need to PASS IN regridded Nfloes
    bool nfloes_interp = M_use_wim;
    if ( !(vm["nextwim.coupling-option"].template as<std::string>() == "breaking_on_mesh"))
        bool nfloes_interp = (M_use_wim && (!M_run_wim));

    if (nfloes_interp)
        std::cout<<"IN REGRID: "<< "interpolate nfloes\n";

    if (nfloes_interp)
        nb_var++;
#endif

	/*Initialize output*/
	double* interp_elt_in=NULL;
	int* interp_method=NULL;
	double* diffusivity_parameters=NULL;

    interp_elt_in=xNew<double>(nb_var*prv_num_elements);
    interp_method=xNew<int>(nb_var); // 0 for non conservative method, 1 for conservative method (for variables defined in terms of blabla/per unit area)
    diffusivity_parameters=xNew<double>(nb_var); // 0 for non added diffusion, positive value for active diffusion in [m^2/s] (only non conservative implementation available)

	int tmp_nb_var=0;
	for (int i=0; i<prv_num_elements; ++i)
	{
		tmp_nb_var=0;

		// concentration
		interp_elt_in[nb_var*i+tmp_nb_var] = M_conc[i];
        interp_method[tmp_nb_var] = 1;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

		// thickness
		interp_elt_in[nb_var*i+tmp_nb_var] = M_thick[i];
        interp_method[tmp_nb_var] = 1;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

		// snow thickness
		interp_elt_in[nb_var*i+tmp_nb_var] = M_snow_thick[i];
        interp_method[tmp_nb_var] = 1;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

		// integrated_stress1
		interp_elt_in[nb_var*i+tmp_nb_var] = M_sigma[3*i]*M_thick[i];
        interp_method[tmp_nb_var] = 1;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

		// integrated_stress2
		interp_elt_in[nb_var*i+tmp_nb_var] = M_sigma[3*i+1]*M_thick[i];
        interp_method[tmp_nb_var] = 1;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

		// integrated_stress3
		interp_elt_in[nb_var*i+tmp_nb_var] = M_sigma[3*i+2]*M_thick[i];
        interp_method[tmp_nb_var] = 1;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

		// damage
		interp_elt_in[nb_var*i+tmp_nb_var] = M_damage[i];
        interp_method[tmp_nb_var] = 0;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

		// ridge_ratio
		interp_elt_in[nb_var*i+tmp_nb_var] = M_ridge_ratio[i]*M_thick[i];
        interp_method[tmp_nb_var] = 1;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

		// random_number
		interp_elt_in[nb_var*i+tmp_nb_var] = M_random_number[i];
        interp_method[tmp_nb_var] = 0;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

		// random_number
		interp_elt_in[nb_var*i+tmp_nb_var] = M_sss[i];
        interp_method[tmp_nb_var] = 0;
        diffusivity_parameters[tmp_nb_var]=vm["simul.diffusivity_sss"].as<double>();
		tmp_nb_var++;
        
		// random_number
		interp_elt_in[nb_var*i+tmp_nb_var] = M_sst[i];
        interp_method[tmp_nb_var] = 0;
        diffusivity_parameters[tmp_nb_var]=vm["simul.diffusivity_sst"].as<double>();
		tmp_nb_var++;

		// Ice temperature
        interp_elt_in[nb_var*i+tmp_nb_var] = M_tice[0][i];
        interp_method[tmp_nb_var] = 0;
        diffusivity_parameters[tmp_nb_var]=0.;
        tmp_nb_var++;

        if ( M_thermo_type == setup::ThermoType::WINTON )
        {
            interp_elt_in[nb_var*i+tmp_nb_var] = ( M_tice[1][i] - physical::mu*physical::si*physical::Lf/(physical::C*M_tice[1][i]) ) * M_thick[i]; // (39) times volume with f1=1
            interp_method[tmp_nb_var] = 1;
            diffusivity_parameters[tmp_nb_var]=0.;
            tmp_nb_var++;

            interp_elt_in[nb_var*i+tmp_nb_var] = ( M_tice[2][i] ) * M_thick[i]; // (39) times volume with f1=0
            interp_method[tmp_nb_var] = 1;
            diffusivity_parameters[tmp_nb_var]=0.;
            tmp_nb_var++;
        }

		// thin ice thickness
		interp_elt_in[nb_var*i+tmp_nb_var] = M_h_thin[i];
        interp_method[tmp_nb_var] = 1;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

		// thin ice thickness
		interp_elt_in[nb_var*i+tmp_nb_var] = M_conc_thin[i];
        interp_method[tmp_nb_var] = 1;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

		// snow on thin ice
		interp_elt_in[nb_var*i+tmp_nb_var] = M_hs_thin[i];
        interp_method[tmp_nb_var] = 1;
		tmp_nb_var++;

		// Ice surface temperature for thin ice
		interp_elt_in[nb_var*i+tmp_nb_var] = M_tsurf_thin[i];
        interp_method[tmp_nb_var] = 0;
        diffusivity_parameters[tmp_nb_var]=0.;
		tmp_nb_var++;

#if defined (WAVES)
        // Nfloes from wim model
        if (nfloes_interp)
        {
            interp_elt_in[nb_var*i+tmp_nb_var] = M_nfloes[i];
            interp_method[tmp_nb_var] = 1;
            diffusivity_parameters[tmp_nb_var]=0.;
            tmp_nb_var++;
        }
#endif

		if(tmp_nb_var>nb_var)
		{
			throw std::logic_error("tmp_nb_var not equal to nb_var");
		}
	}
	*interp_elt_in_ptr=interp_elt_in;
    *interp_method_ptr=interp_method;
    *diffusivity_parameters_ptr=diffusivity_parameters;


    return nb_var;
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
        else
        {
            bamggeom_previous->Edges[3*edg+2] = M_flag_fix+1; // we just want it to be different than M_flag_fix 
            bamgmesh_previous->Edges[3*edg+2] = M_flag_fix+1; // we just want it to be different than M_flag_fix 
        }
        
    }

    //Environment::logMemoryUsage("before adaptMesh");
    Bamgx(bamgmesh,bamggeom,bamgmesh_previous,bamggeom_previous,bamgopt_previous);
    //Environment::logMemoryUsage("after adaptMesh");

    // Save the old id_node before redefining it
    std::vector<int> old_node_id=M_mesh.id();

    // Import the mesh from bamg
    this->importBamg(bamgmesh);

    // We mask out the boundary nodes
    M_mask.assign(bamgmesh->VerticesSize[0],false) ;
    M_mask_dirichlet.assign(bamgmesh->VerticesSize[0],false) ;

    for (int vert=0; vert<bamgmesh->VerticesOnGeomVertexSize[0]; ++vert)
        M_mask[bamgmesh->VerticesOnGeomVertex[2*vert]-1]=true; // The factor 2 is because VerticesOnGeomVertex has 2 dimensions in bamg

        // Recompute the node ids
    if(bamgopt->KeepVertices)
    {
        std::vector<int> new_nodes_id=M_mesh.id();

        int Boundary_id=0;
        int nb_new_nodes=0;

        // The new id will have values higher than the previous one
        int first_new_node=*std::max_element(old_node_id.begin(),old_node_id.end())+1;

        for (int vert=0; vert<bamgmesh->VerticesSize[0]; ++vert)
        {
            if(M_mask[vert])
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
        M_mask_dirichlet[M_dirichlet_flags[i]]=true;
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

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        // if(thread_id == 0)
        // {
        //     std::cout<< "\nMaster thread being called\n";
        //     total_threads = omp_get_num_threads();
        //     std::cout<<"Total number of threads are "<< total_threads <<"\n";
        // }

        // total thickness and concentration
        double total_concentration=M_conc[cpt];
        double total_thickness=M_thick[cpt];
        double total_snow=M_snow_thick[cpt];
       
        // Add the thin ice concentration and thickness 
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {            
            total_concentration += M_conc_thin[cpt];
            total_thickness     += M_h_thin[cpt];
            total_snow          += M_hs_thin[cpt];
        }

        // Limits to avoid very small values
        total_thickness =       (vm["simul.min_h"].as<double>()>total_thickness)        ? vm["simul.min_h"].as<double>() : total_thickness;
        total_concentration =   (vm["simul.min_c"].as<double>()>total_concentration)    ? vm["simul.min_c"].as<double>() : total_concentration;

        int index_u, index_v;

        double coef_Vair    = 0.;
        double coef_Voce    = 0.;
        double coef_basal   = 0.;
        double coef         = 10.;
        double mass_e       = 0.;
        double coef_C       = 0.;
        double coef_V       = 0.;
        double coef_X       = 0.;
        double coef_Y       = 0.;

        double mloc = 0.;
        double dloc = 0.;
        
        double b0tj_sigma_hu = 0.;
        double b0tj_sigma_hv = 0.;

        double undamaged_time_relaxation_sigma=vm["simul.undamaged_time_relaxation_sigma"].as<double>();
        double exponent_relaxation_sigma=vm["simul.exponent_relaxation_sigma"].as<double>();

        double time_viscous=undamaged_time_relaxation_sigma*std::pow(1.-M_damage[cpt],exponent_relaxation_sigma-1.);
        double multiplicator=time_viscous/(time_viscous+time_step);
        
        double norm_Voce_ice=0.;
        double norm_Voce_ice_min= 0.01; // minimum value to avoid 0 water drag term.
        
        double norm_Vair_ice=0.;
        double norm_Vair_ice_min= 0.01; // minimum value to avoid 0 water drag term.
        
        double norm_Vice=0.;
        
        double element_ssh=0.;
        double critical_h = 0.;
        double max_keel_height=28; // [m] from "A comprehensive analysis of the morphology of first-year sea ice ridges"
        double ice_to_keel_factor=19.28; // from "A comprehensive analysis of the morphology of first-year sea ice ridges"
        double keel_height_estimate;
        double critical_h_mod=0.; 
        
        //if(total_concentration > vm["simul.min_c"].as<double>())
        if( (total_concentration > vm["simul.min_c"].as<double>()) && (total_thickness > vm["simul.min_h"].as<double>()) )
        {

            /* Compute the value that only depends on the element */
            double welt_ice = 0.;
            double welt_ssh = 0.;
            int nind;

            for (int i=0; i<3; ++i)
            {
                nind = (M_elements[cpt]).indices[i]-1;
                welt_ssh += M_ssh[nind];
            }
            
            element_ssh = welt_ssh/3.;
            
           
            if(M_conc[cpt]>vm["simul.min_c"].as<double>())
            {
                switch ( M_basal_stress_type )
                {
                    case setup::BasalStressType::BOUILLON:
                        // Sylvain's grounding scheme
                        keel_height_estimate = ice_to_keel_factor*std::pow(M_thick[cpt]/M_conc[cpt],0.5);
                        keel_height_estimate = ( keel_height_estimate > max_keel_height ) ? max_keel_height : keel_height_estimate;
                        critical_h      = M_conc[cpt]*std::pow((M_element_depth[cpt]+element_ssh)/ice_to_keel_factor,2.); 
                        critical_h_mod  = M_conc[cpt]*std::pow(keel_height_estimate/ice_to_keel_factor,2.); 
                        break;
                    case setup::BasalStressType::LEMIEUX:
                        // JF Lemieux's grounding
                        keel_height_estimate = vm["simul.Lemieux_basal_k1"].as<double>()*M_thick[cpt]/M_conc[cpt];
                        keel_height_estimate = ( keel_height_estimate > max_keel_height ) ? max_keel_height : keel_height_estimate;

                        critical_h      = M_conc[cpt]*(M_element_depth[cpt]+element_ssh)/(vm["simul.Lemieux_basal_k1"].as<double>());
                        critical_h_mod  = M_conc[cpt]*keel_height_estimate/(vm["simul.Lemieux_basal_k1"].as<double>()); 
                        break;
                }
            } 

    #if 1
            //option 1 (original)
            coef = multiplicator*young*(1.-M_damage[cpt])*M_thick[cpt]*std::exp(ridging_exponent*(1.-M_conc[cpt]));

    #else
            //option 2 (we just change the value of the ridging exponent and we renamed it "damaging_exponent")
            double damaging_exponent = -80.;
            double coef = young*(1.-M_damage[cpt])*M_thick[cpt]*std::exp(damaging_exponent*(1.-M_conc[cpt]));
    #endif
            //option 3: We change the formulation of f(A) and make it piecewise linear between limit_conc_fordamage and 1, and 0 otherwise
            //double factor = 0.;
            //double limit_conc_fordamage = 0.;
            //limit_conc_fordamage=0.95;
            //if(M_conc[cpt]<limit_conc_fordamage)
            //{
            //factor=0.;
            //}
            //else
            //{
            //factor=(M_conc[cpt]-limit_conc_fordamage)/(1.-limit_conc_fordamage);
            //}
            //double coef = young*(1.-M_damage[cpt])*M_thick[cpt]*factor;

            if (vm["simul.use_coriolis"].as<bool>())
                mass_e = (rhoi*total_thickness + rhos*total_snow)/total_concentration;
            else
                mass_e=0.;

            // /* compute the x and y derivative of g*ssh */
            double g_ssh_e_x = 0.;
            double g_ssh_e_y = 0.;
            double g_ssh_e;
            for(int i=0; i<3; i++)
            {
                g_ssh_e = (physical::gravity)*M_ssh[(M_elements[cpt]).indices[i]-1] /*g_ssh*/;   /* g*ssh at the node k of the element e */
                g_ssh_e_x += M_shape_coeff[cpt][i]*g_ssh_e; /* x derivative of g*ssh */
                g_ssh_e_y += M_shape_coeff[cpt][i+3]*g_ssh_e; /* y derivative of g*ssh */
            }

            coef_C     = mass_e*M_fcor[cpt];              /* for the Coriolis term */
            coef_V     = mass_e/time_step;             /* for the inertial term */
            coef_X     = - mass_e*g_ssh_e_x;              /* for the ocean slope */
            coef_Y     = - mass_e*g_ssh_e_y;              /* for the ocean slope */
        }

        /* Loop over the 6 by 6 components of the finite element integral
         * this is done smartly by looping over j=0:2 and i=0:2
         * col = (mwIndex)it[2*j]-1  , row = (mwIndex)it[2*i]-1;
         * col  , row   -> UU component
         * col  , row+1 -> VU component
         * col+1, row   -> VV component
         * col+1, row+1 -> UV component */

        double Vcor_index_v, Vcor_index_u;
        double surface_e = M_surface[cpt];
        double duu, dvu, duv, dvv;
        std::vector<double> data(36);
        std::vector<double> fvdata(6,0.);
        std::vector<int> rcindices(6);

        for(int j=0; j<3; j++)
        {
            /* Column corresponding to indice j (we also assemble terms in col+1) */
            //col = (mwIndex)it[2*j]-1; /* -1 to use the indice convention of C */

            index_u = (M_elements[cpt]).indices[j]-1;
            index_v = (M_elements[cpt]).indices[j]-1+M_num_nodes;

            Vcor_index_v=beta0*M_VT[index_v] + beta1*M_VTM[index_v] + beta2*M_VTMM[index_v];
            Vcor_index_u=beta0*M_VT[index_u] + beta1*M_VTM[index_u] + beta2*M_VTMM[index_u];

            double coef_sigma = M_thick[cpt]*multiplicator;

            for(int i=0; i<3; i++)
            {
                /* Row corresponding to indice i (we also assemble terms in row+1) */
                //row = (mwIndex)it[2*i]-1; /* -1 to use the indice convention of C */

                /* Select the nodal weight values from M_loc */
                mloc = M_Mass[3*j+i];
                dloc = M_Diag[3*j+i];

                b0tj_sigma_hu = 0.;
                b0tj_sigma_hv = 0.;

                for(int k=0; k<3; k++)
                {
                    b0tj_sigma_hu += M_B0T[cpt][k*6+2*i]*(M_sigma[3*cpt+k]*coef_sigma/*+sigma_P[k]*/);
                    b0tj_sigma_hv += M_B0T[cpt][k*6+2*i+1]*(M_sigma[3*cpt+k]*coef_sigma/*+sigma_P[k]*/);
                }

                /* ---------- UU component */
                
                norm_Voce_ice = std::hypot(M_VT[index_u]-M_ocean[index_u],M_VT[index_v]-M_ocean[index_v]);
                norm_Voce_ice = (norm_Voce_ice > norm_Voce_ice_min) ? (norm_Voce_ice):norm_Voce_ice_min;

                coef_Voce = (vm["simul.lin_drag_coef_water"].as<double>()+(quad_drag_coef_water*norm_Voce_ice));
                coef_Voce *= physical::rhow; //(vm["simul.rho_water"].as<double>());
                
                norm_Vair_ice = std::hypot(M_VT[index_u]-M_wind [index_u],M_VT[index_v]-M_wind [index_v]);
                norm_Vair_ice = (norm_Vair_ice > norm_Vair_ice_min) ? (norm_Vair_ice):norm_Vair_ice_min;

                coef_Vair = (vm["simul.lin_drag_coef_air"].as<double>()+(quad_drag_coef_air*norm_Vair_ice));
                coef_Vair *= (physical::rhoa);
                
                norm_Vice = std::hypot(M_VT[index_u],M_VT[index_v]);
                norm_Vice = (norm_Vice > basal_u_0) ? (norm_Vice):basal_u_0;
                                
                coef_basal = basal_k2/norm_Vice;
                //coef_basal *= std::max(0., M_thick[cpt]-critical_h)*std::exp(-basal_Cb*(1.-M_conc[cpt]));
                coef_basal *= std::max(0., critical_h_mod-critical_h)*std::exp(-basal_Cb*(1.-M_conc[cpt]));
                   
                duu = surface_e*( mloc*(coef_V)
                                  +dloc*(coef_Vair+coef_basal+coef_Voce*cos_ocean_turning_angle)  
                                  +M_B0T_Dunit_B0T[cpt][(2*i)*6+2*j]*coef*time_step);

                /* ---------- VU component */
                dvu = surface_e*(+M_B0T_Dunit_B0T[cpt][(2*i+1)*6+2*j]*coef*time_step);

                /* ---------- UV component */
                duv = surface_e*(+M_B0T_Dunit_B0T[cpt][(2*i)*6+2*j+1]*coef*time_step);

                /* ---------- VV component */
                dvv = surface_e*( mloc*(coef_V)
                                  +dloc*(coef_Vair+coef_basal+coef_Voce*cos_ocean_turning_angle)
                                  +M_B0T_Dunit_B0T[cpt][(2*i+1)*6+2*j+1]*coef*time_step);

                data[(2*i  )*6+2*j  ] = duu;
                data[(2*i+1)*6+2*j  ] = dvu;
                data[(2*i  )*6+2*j+1] = duv;
                data[(2*i+1)*6+2*j+1] = dvv;

                fvdata[2*i]     += surface_e*( mloc*(   
                                                    +M_tau[index_u]
                                                    +coef_X
                                                    +coef_V*M_VT[index_u]
                                                    +coef_C*Vcor_index_v
                                                        )
                                                + dloc*(   
                                                    +coef_Vair*M_wind[index_u]
                                                    +coef_Voce*cos_ocean_turning_angle*M_ocean[index_u]
                                                    -coef_Voce*sin_ocean_turning_angle*(M_ocean[index_v]-M_VT[index_v])
                                                    )
                                            - b0tj_sigma_hu/3);

                fvdata[2*i+1]   += surface_e*( mloc*( 
                                                    +M_tau[index_v]
                                                    +coef_Y
                                                    +coef_V*M_VT[index_v]
                                                    -coef_C*Vcor_index_u
                                                        )
                                                +dloc*( 
                                                    +coef_Vair*M_wind[index_v]
                                                    +coef_Voce*cos_ocean_turning_angle*M_ocean[index_v]
                                                    +coef_Voce*sin_ocean_turning_angle*(M_ocean[index_u]-M_VT[index_u])
                                                        )
                                            - b0tj_sigma_hv/3);
            }

            rcindices[2*j] = index_u;
            rcindices[2*j+1] = index_v;
        }

#if 0
            std::cout<<"************************\n";
            std::cout<<"Coef_C    = "<< coef_C <<"\n";
            std::cout<<"Coef_V    = "<< coef_V <<"\n";
            std::cout<<"Coef_X    = "<< coef_X <<"\n";
            std::cout<<"Coef_Y    = "<< coef_Y <<"\n";
            std::cout<<"coef_Vair = "<< coef_Vair <<"\n";
            std::cout<<"coef_Voce = "<< coef_Voce <<"\n";
            std::cout<<"coef_basal= "<< coef_basal <<"\n";
#endif

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

    chrono.restart();
    M_matrix->on(M_dirichlet_nodes,*M_vector);
    //M_matrix->on(extended_dirichlet_nodes,*M_vector);
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
FiniteElement::tensors()
{
    M_Dunit.assign(9,0);
    M_Mass.assign(9,0);
    M_Diag.assign(9,0);

    for (int k=0; k<6; k+=3)
    {
        for (int kk=0; kk<2; ++kk )
        {
            M_Dunit[k+kk] = (1-((k+kk)%2)*(1-nu0))/(1-std::pow(nu0,2.));
        }
    }
    M_Dunit[8] = (1-nu0)/(2.*(1-std::pow(nu0,2.)));
#if 0
    std::cout<< " M_Dunit: ";
    for (int k=0; k<9; k++)
    {
        std::cout<< M_Dunit[k] << " ";
    }
    std::cout<< " \n ";

#endif

    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            M_Mass[3*i+j] = ((i == j) ? 2.0 : 1.0)/12.0;
            M_Diag[3*i+j] = ((i == j) ? 1.0 : 0.0)/3.0;
            //std::cout<< std::left << std::setw(12) << Mass[3*i+j] <<"  ";
        }

        //std::cout<<"\n";
    }

#if 1

    M_B0T.resize(M_num_elements);
    M_B0T_Dunit_B0T.resize(M_num_elements);
    M_shape_coeff.resize(M_num_elements);

    std::vector<double> B0T(18,0);
    std::vector<double> B0Tj_Dunit(6,0);
    std::vector<double> B0T_Dunit_B0T(36,0);

    double B0Tj_Dunit_tmp0, B0Tj_Dunit_tmp1;
    double B0Tj_Dunit_B0Ti_tmp0, B0Tj_Dunit_B0Ti_tmp1, B0Tj_Dunit_B0Ti_tmp2, B0Tj_Dunit_B0Ti_tmp3;

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

                for(int kk=0; kk<3; kk++)
                {
                    B0Tj_Dunit_tmp0 += B0T[kk*6+2*j]*M_Dunit[3*i+kk];
                    B0Tj_Dunit_tmp1 += B0T[kk*6+2*j+1]*M_Dunit[3*i+kk];
                }

                B0Tj_Dunit[2*i] = B0Tj_Dunit_tmp0;
                B0Tj_Dunit[2*i+1] = B0Tj_Dunit_tmp1;
            }

            for(int i=0; i<3; i++)
            {
                /* The rigidity matrix */
                /* scalar product of B0Ti_Dunit and the first column of B0T */
                B0Tj_Dunit_B0Ti_tmp0 = 0.;
                B0Tj_Dunit_B0Ti_tmp1 = 0.;
                B0Tj_Dunit_B0Ti_tmp2 = 0.;
                B0Tj_Dunit_B0Ti_tmp3 = 0.;

                for(int kk=0; kk<3; kk++)
                {
                    B0Tj_Dunit_B0Ti_tmp0 += B0Tj_Dunit[2*kk]*B0T[kk*6+2*i];
                    B0Tj_Dunit_B0Ti_tmp1 += B0Tj_Dunit[2*kk]*B0T[kk*6+2*i+1];
                    B0Tj_Dunit_B0Ti_tmp2 += B0Tj_Dunit[2*kk+1]*B0T[kk*6+2*i];
                    B0Tj_Dunit_B0Ti_tmp3 += B0Tj_Dunit[2*kk+1]*B0T[kk*6+2*i+1];
                }

                B0T_Dunit_B0T[(2*i)*6+2*j] = B0Tj_Dunit_B0Ti_tmp0;
                B0T_Dunit_B0T[(2*i+1)*6+2*j] = B0Tj_Dunit_B0Ti_tmp1;
                B0T_Dunit_B0T[(2*i)*6+2*j+1] = B0Tj_Dunit_B0Ti_tmp2;
                B0T_Dunit_B0T[(2*i+1)*6+2*j+1] = B0Tj_Dunit_B0Ti_tmp3;
            }
        }

        /*
         * B0T_Dunit_B0T should be symmetric but is not exactly after the calcultation here above
         * because the sequence of operation is not the same for the component i,j and j,i.
         * We force the matrix to be symmetric by copying the upper part onto the lower part
         */
        for(int i=1; i<6; i++)
            for(int j=0; j<i; j++)
                B0T_Dunit_B0T[i*6+j] = B0T_Dunit_B0T[j*6+i];

        M_shape_coeff[cpt]        = shapecoeff;
        M_B0T[cpt]                = B0T;
        M_B0T_Dunit_B0T[cpt]      = B0T_Dunit_B0T;

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
    int thread_id;
    int total_threads;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    //std::cout<<"MAX THREADS= "<< max_threads <<"\n";
#if defined (WAVES)
    if (M_use_wim)
        M_dfloe.assign(M_num_elements,0.);
#endif


    // collect the variables into a single structure
    int prv_num_elements = M_mesh.numTriangles();
    double* interp_elt_in;
    int* interp_method;
    double* diffusivity_parameters;
    int nb_var = this->collectVariables(&interp_elt_in, &interp_method, &diffusivity_parameters, prv_num_elements);

    // Advection
    double* interp_elt_out;
    this->advect(&interp_elt_out,&interp_elt_in[0],&interp_method[0],nb_var);
	
    // cleaning
    xDelete<double>(interp_elt_in);
    xDelete<int>(interp_method);
            
    double* interp_elt_out_bis;
    this->diffuse(&interp_elt_out_bis,&interp_elt_out[0],&diffusivity_parameters[0],nb_var,this->resolution(M_mesh));

    // redistribute the interpolated values
    this->redistributeVariables(&interp_elt_out_bis[0],nb_var);

    // cleaning
    xDelete<double>(interp_elt_out);
    xDelete<double>(interp_elt_out_bis);
    xDelete<double>(diffusivity_parameters);
    
#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        double old_damage;

        /* deformation, deformation rate and internal stress tensor and temporary variables */
        double epsilon_veloc_i;
        std::vector<double> epsilon_veloc(3);
        double divergence_rate;
        double shear_rate;

        // ridging scheme
        double delta_ridging;
        double G_star=0.15;
        double e_factor=2.;

        std::vector<double> sigma_pred(3);
        double sigma_dot_i;

        /* invariant of the internal stress tensor and some variables used for the damaging process*/
        double sigma_s, sigma_n, sigma_1, sigma_2;
        double tract_max, sigma_t, sigma_c, q;
        double tmp, sigma_target;

        // beginning of the original code (without openMP)

        // Temporary memory
        old_damage = M_damage[cpt];
        
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
                epsilon_veloc_i += M_B0T[cpt][i*6 + 2*j]*M_VT[(M_elements[cpt]).indices[j]-1]  ;
                epsilon_veloc_i += M_B0T[cpt][i*6 + 2*j + 1]*M_VT[(M_elements[cpt]).indices[j]-1+M_num_nodes]  ;
            }

            epsilon_veloc[i] = epsilon_veloc_i;
        }

        divergence_rate= (epsilon_veloc[0]+epsilon_veloc[1]);
        shear_rate= std::hypot(epsilon_veloc[0]-epsilon_veloc[1],epsilon_veloc[2]);
        delta_ridging= std::hypot(divergence_rate,shear_rate/e_factor);

       /*======================================================================
        * Ridging scheme 
        * After the advection the concentration can be higher than 1, meaning that ridging should have occured.
        *======================================================================
        */
        //if(M_conc[cpt]>1.)
        //{
        //    M_ridge_ratio[cpt]=M_ridge_ratio[cpt]+(1.-M_ridge_ratio[cpt])*(M_conc[cpt]-1.)/M_conc[cpt];
        //    M_conc[cpt]=1.;
        //}
        //
        double open_water_concentration=1.-M_conc[cpt];
        
        /* Thin ice category */    
        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            open_water_concentration-=M_conc_thin[cpt];
        }
        // limit open_water concentration to 0.
        //open_water_concentration=(open_water_concentration<0.)?0.:open_water_concentration;

        // ridging scheme
        double opening_factor=((1.-M_conc[cpt])>G_star) ? 0. : std::pow(M_conc[cpt]/G_star,2.);
        //open_water_concentration += time_step*0.5*(delta_ridging-divergence_rate)*opening_factor;
        open_water_concentration += time_step*0.5*shear_rate/e_factor*opening_factor;
        
        /* Thin ice category */
        double new_conc_thin=0.;   
        double new_h_thin=0.;   
        double new_hs_thin=0.;   
        
        double newice = 0.;
        double del_c = 0.;
        double newsnow = 0.;
        
        double ridge_thin_ice_aspect_ratio=10.;
#if 1
        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            if(M_conc_thin[cpt]>0. && M_thick[cpt]>0.)
            {
            new_conc_thin   = std::max(1.-M_conc[cpt]-open_water_concentration,0.);
            new_h_thin      = new_conc_thin*M_h_thin[cpt]/M_conc_thin[cpt]; // so that we keep the same h0, no preferences for the ridging
            new_hs_thin     = new_conc_thin*M_hs_thin[cpt]/M_conc_thin[cpt];
 
            newice = M_h_thin[cpt]-new_h_thin;
            del_c   = (M_conc_thin[cpt]-new_conc_thin)/ridge_thin_ice_aspect_ratio;            
            newsnow = M_hs_thin[cpt]-new_hs_thin;

            M_conc_thin[cpt]= new_conc_thin;
            M_h_thin[cpt]   = new_h_thin;
            M_hs_thin[cpt]  = new_hs_thin;
            
            M_thick[cpt]        += newice;
            M_conc[cpt]         += del_c;
            M_snow_thick[cpt]   += newsnow;
            
            M_ridge_ratio[cpt]=std::max(0.,std::min(1.,(M_ridge_ratio[cpt]*(M_thick[cpt]-newice)+newice)/M_thick[cpt]));
            }
            else
            {
                M_conc_thin[cpt]=0.;
                M_h_thin[cpt]=0.;
                M_hs_thin[cpt]=0.;
                M_ridge_ratio[cpt]=0.;
            }
        }
#endif
        double new_conc=std::max(1.-M_conc_thin[cpt]-open_water_concentration+del_c,0.);
        if(new_conc<M_conc[cpt])
        {
            M_ridge_ratio[cpt]=std::max(0.,std::min(1.,(M_ridge_ratio[cpt]+(1.-M_ridge_ratio[cpt])*(M_conc[cpt]-new_conc)/M_conc[cpt])));
        }
        M_conc[cpt]=new_conc;
    
        double max_true_thickness = 50;
        if(M_conc[cpt]>0.)
        {
            double test_h_thick=M_thick[cpt]/M_conc[cpt];
            test_h_thick = (test_h_thick>max_true_thickness) ? max_true_thickness : test_h_thick ;
            M_conc[cpt]=M_thick[cpt]/test_h_thick;
        }
    else
        {
            M_ridge_ratio[cpt]=0.;
            M_thick[cpt]=0.;
            M_snow_thick[cpt]=0.;
        }
        
#if 0
        /* Initialise to be safe */
        newice = 0.;
        del_c = 0.;
        newsnow = 0.;
        /* Thin ice category */    
        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            thin_ice_redistribute(M_h_thin[cpt], M_hs_thin[cpt], 0., M_conc[cpt],
                          tanalpha, rtanalpha, h_thin_max, &M_h_thin[cpt], &newice, &del_c, &newsnow);
        
            // Change the snow _thickness_ for thick ice and _volume_ for thin ice
            M_hs_thin[cpt] -= newsnow;
            M_snow_thick[cpt] += newsnow;
            M_conc[cpt] += del_c;
            M_thick[cpt] += newice;
            
            if(M_thick[cpt]>0.)
                M_ridge_ratio[cpt]=(M_ridge_ratio[cpt]*(M_thick[cpt]-newice)+newice)/M_thick[cpt];
            else
                M_ridge_ratio[cpt]=0.;
        }
#endif   
        


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
        double undamaged_time_relaxation_sigma=vm["simul.undamaged_time_relaxation_sigma"].as<double>();
        double exponent_relaxation_sigma=vm["simul.exponent_relaxation_sigma"].as<double>();

        double time_viscous=undamaged_time_relaxation_sigma*std::pow(1.-old_damage,exponent_relaxation_sigma-1.);
        double multiplicator=time_viscous/(time_viscous+time_step);

        for(int i=0;i<3;i++)
        {
            sigma_dot_i = 0.0;
            for(int j=0;j<3;j++)
            {
            // sigma_dot_i += std::exp(ridging_exponent*(1.-old_conc))*young*(1.-old_damage)*M_Dunit[i*3 + j]*epsilon_veloc[j];
            sigma_dot_i += std::exp(damaging_exponent*(1.-M_conc[cpt]))*young*(1.-old_damage)*M_Dunit[i*3 + j]*epsilon_veloc[j];
            //sigma_dot_i += factor*young*(1.-old_damage)*M_Dunit[i*3 + j]*epsilon_veloc[j];
            }
            
            sigma_pred[i] = (M_sigma[3*cpt+i]+4*time_step*sigma_dot_i)*multiplicator;
            sigma_pred[i] = (M_conc[cpt] > vm["simul.min_c"].as<double>()) ? (sigma_pred[i]):0.;
            
            M_sigma[3*cpt+i] = (M_sigma[3*cpt+i]+time_step*sigma_dot_i)*multiplicator;
            M_sigma[3*cpt+i] = (M_conc[cpt] > vm["simul.min_c"].as<double>()) ? (M_sigma[3*cpt+i]):0.;
            
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

        double ridge_to_normal_cohesion_ratio=vm["simul.ridge_to_normal_cohesion_ratio"].as<double>();
        double effective_cohesion=M_Cohesion[cpt]*(1. + M_ridge_ratio[cpt]*(ridge_to_normal_cohesion_ratio-1.) );
        double effective_compressive_strength=M_Compressive_strength[cpt]*(1. + M_ridge_ratio[cpt]*(ridge_to_normal_cohesion_ratio-1.) );

        q=std::pow(std::pow(std::pow(tan_phi,2.)+1,.5)+tan_phi,2.);
        sigma_c=2.*effective_cohesion/(std::pow(std::pow(tan_phi,2.)+1,.5)-tan_phi);
        sigma_t=-sigma_c/q;

        /* minimum and maximum normal stress */
        tract_max=-tract_coef*effective_cohesion/tan_phi;

        /* Correction of the damage */
        if(sigma_n>effective_compressive_strength)
        {
            sigma_target=effective_compressive_strength;

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
        for(int i=0;i<3;i++)
        {
#if 0
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
         *======================================================================
         */

        /* lower bounds */
        M_conc[cpt] = ((M_conc[cpt]>0.)?(M_conc[cpt] ):(0.)) ;
        M_thick[cpt]        = ((M_thick[cpt]>0.)?(M_thick[cpt]     ):(0.)) ;
        M_snow_thick[cpt]   = ((M_snow_thick[cpt]>0.)?(M_snow_thick[cpt]):(0.)) ;

        /* Ice damage
        * We use now a constant healing rate defined as 1/time_recovery_damage
        * so that we are now able to reset the damage to 0.
        * otherwise, it will never heal completely.
        * time_recovery_damage still depends on the temperature when themodynamics is activated.
        */
        tmp=M_damage[cpt]-time_step/M_time_relaxation_damage[cpt];
        if(M_thick[cpt]==0.)
            tmp=0.;
        M_damage[cpt]=((tmp>0.)?(tmp):(0.));

#if 0
        if (M_conc[cpt]>physical::cmin)
        {
            // sanity checks on thickness and snow thickness
            if (M_thick[cpt]/M_conc[cpt]>50.)
            {
                std::cout<<"True thickness too large\n";
                std::cout<<"c,vol,h = "<<M_conc[cpt]<<","<<M_thick[cpt]<<","<<M_thick[cpt]/M_conc[cpt]<<"\n";
                throw std::runtime_error("True thickness too large");
            }
            if (M_snow_thick[cpt]/M_conc[cpt]>20.)
            {
                std::cout<<"True snow thickness too large\n";
                std::cout<<"c,vol,h = "<<M_conc[cpt]<<","<<M_snow_thick[cpt]<<M_snow_thick[cpt]/M_conc[cpt]<<"\n";
                throw std::runtime_error("True snow thickness too large");
            }
        }
#endif

    }//loop over elements
}//update

void
FiniteElement::solve()
{
    M_solver->solve(_matrix=M_matrix,
                    _solution=M_solution,
                    _rhs=M_vector,
                    _ksp=vm["solver.ksp-type"].as<std::string>()/*"preonly"*/,
                    _pc=vm["solver.pc-type"].as<std::string>()/*"cholesky"*/,
                    _pcfactormatsolverpackage=vm["solver.mat-package-type"].as<std::string>()/*"cholmod"*/,
                    _reuse_prec=true,
                    _rebuild=M_regrid
                    );

    //M_solution->printMatlab("solution.m");
    //Environment::logMemoryUsage("");
}

// Routine for the 1D thermodynamical model
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

    double const ocean_albedo = vm["simul.albedoW"].as<double>();
    double const drag_ocean_t = vm["simul.drag_ocean_t"].as<double>();
    double const drag_ocean_q = vm["simul.drag_ocean_q"].as<double>();

    double const rh0   = 1./vm["simul.hnull"].as<double>();
    double const rPhiF = 1./vm["simul.PhiF"].as<double>();

    double const qi = physical::Lf * physical::rhoi;// J m^{-3}
    double const qs = physical::Lf * physical::rhos;// J m^{-3}

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
            old_conc_thin  = M_conc_thin[i];
            old_hs_thin = M_hs_thin[i];
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

        // definition of the snow fall in kg/m^2/s
        double tmp_snowfall;
        if(M_snowfr.M_initialized)
            tmp_snowfall=M_precip[i]*M_snowfr[i];
        else if (M_snowfall.M_initialized)
            tmp_snowfall=M_snowfall[i];
        else
        {
            if(M_tair[i]<0)
                tmp_snowfall=M_precip[i];
            else
                tmp_snowfall=0.;
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
        if ( M_sphuma.M_initialized )
            sphuma = M_sphuma[i];
        else if ( M_dair.M_initialized )
        {
            double fa     = 1. + Aw + M_mslp[i]*1e-2*( Bw + Cw*M_dair[i]*M_dair[i] );
            double esta   = fa*aw*std::exp( (bw-M_dair[i]/dw)*M_dair[i]/(M_dair[i]+cw) );
            sphuma = alpha*fa*esta/(M_mslp[i]-beta*fa*esta) ;
        }
        else if ( M_mixrat.M_initialized )
            sphuma = M_mixrat[i]/(1.+M_mixrat[i]) ;
        else
            throw std::logic_error("Neither M_dair nor M_mixrat have been initialized. I cannot calculate sphuma.");

        // -------------------------------------------------
        // 3.2) Specific humidity - ocean surface (calcSphumW in matlab)
        double fw     = 1. + Aw + M_mslp[i]*1e-2*( Bw + Cw*M_sst[i]*M_sst[i] );
        double estw   = aw*std::exp( (bw-M_sst[i]/dw)*M_sst[i]/(M_sst[i]+cw) )*(1-5.37e-4*M_sss[i]);
        sphumw = alpha*fw*estw/(M_mslp[i]-beta*fw*estw) ;

        // -------------------------------------------------
        /* Density of air */
        double rhoair = M_mslp[i]/(physical::Ra*(M_tair[i]+physical::tfrwK)) * (1.+sphuma)/(1.+1.609*sphuma);

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
            double tsa = M_tice[0][i] + physical::tfrwK;
            double taa = M_tair[i]  + physical::tfrwK;
            // s.b.idso & r.d.jackson, thermal radiation from the atmosphere, j. geophys. res. 74, 5397-5403, 1969
        	tmp_Qlw_in = sigma_sb*std::pow(taa,4) \
        			*( 1. - 0.261*std::exp(-7.77e-4*std::pow(taa-physical::tfrwK,2)) ) \
        			*( 1. + 0.275*M_tcc[i] );
        }

        // Qow>0 => flux out of ocean:
        // - subtract shortwave and longwave input;
        // add heat loss from longwave radiation, sensible heat loss (temp changes)
        // and evaporation (latent heat loss - temp doesn't change, but phase changes)
        Qow = -tmp_Qsw_in*(1.-ocean_albedo) - tmp_Qlw_in + Qlw_out + Qsh + Qlh;

        // -------------------------------------------------
        // 4) Thickness change of the ice slab (thermoIce0 in matlab)

        switch ( M_thermo_type )
        {
            case setup::ThermoType::ZERO_LAYER:
                this->thermoIce0(i, wspeed, sphuma, M_conc[i], M_thick[i], M_snow_thick[i],
                        tmp_Qlw_in, tmp_Qsw_in, tmp_mld, tmp_snowfall, hi, hs, hi_old, Qio, del_hi, M_tice[0][i]);
                break;
            case setup::ThermoType::WINTON:
                this->thermoWinton(i, time_step, wspeed, sphuma, M_conc[i], M_thick[i], M_snow_thick[i],
                        tmp_Qlw_in, tmp_Qsw_in, tmp_mld, tmp_snowfall, hi, hs, hi_old, Qio, del_hi,
                        M_tice[0][i], M_tice[1][i], M_tice[2][i]);
                break;
        }

        if ( M_ice_cat_type==setup::IceCategoryType::THIN_ICE )
        {
            this->thermoIce0(i, wspeed, sphuma, old_conc_thin, M_h_thin[i], M_hs_thin[i],
                    tmp_Qlw_in, tmp_Qsw_in, tmp_mld, tmp_snowfall, hi_thin, hs_thin, hi_thin_old, Qio_thin, del_hi_thin, M_tsurf_thin[i]);
            M_h_thin[i]  = hi_thin * old_conc_thin;
            M_hs_thin[i] = hs_thin * old_conc_thin;
        }

        // -------------------------------------------------
        // 5) Ice growth over open water and lateral melt (thermoOW in matlab)

        /* Local variables */
        double tw_new, tfrw, newice, del_c, newsnow, h0;

        /* dT/dt due to heatflux ocean->atmosphere */
        tw_new = M_sst[i] - Qow*time_step/(tmp_mld*physical::rhow*physical::cpw);
        tfrw   = -physical::mu*M_sss[i];

        /* Form new ice in case of super cooling, and reset Qow and evap */
        if ( tw_new < tfrw )
        {
            newice  = (1.-M_conc[i]-M_conc_thin[i])*(tfrw-tw_new)*tmp_mld*physical::rhow*physical::cpw/qi;// m
            Qow  = -(tfrw-M_sst[i])*tmp_mld*physical::rhow*physical::cpw/time_step;
            // evap = 0.;
        } else {
            newice  = 0.;
        }

        /* Decide the change in ice fraction (del_c) */
        /* Initialise to be safe */
        del_c = 0.;
        newsnow = 0.;
        
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
                        if ( newice > 0. )
                            del_c = 1.;
                        else
                            del_c = 0.;
                    }
                    break;
                case 3:
                    /* Olason and Harms (09) */
                    h0    = (1.+0.1*wspeed)/15.;
                    del_c = newice/std::max(rPhiF*hi_old,h0);
                    break;
                case 4:
                    /* Thin ice category */

                    M_h_thin[i]+=newice;
                    M_conc_thin[i]=std::min(1.-M_conc[i],M_conc_thin[i]+newice/h_thin_min);
                    newice  = 0.;
                    newsnow = 0.;

                    if(M_conc_thin[i]>0.)
                    {    
                        /* Two cases: Thin ice fills the cell or not */
                            if ( M_h_thin[i] < h_thin_min*M_conc_thin[i] )
                                M_conc_thin[i] = M_h_thin[i]/h_thin_min;
                        else 
                        {    
                            h0 = h_thin_min + 2.*(M_h_thin[i]-h_thin_min*M_conc_thin[i])/(M_conc_thin[i]);
                            if(h0>h_thin_max)
                            {    
                                del_c = M_conc_thin[i]/(h0-h_thin_min) * (h0-h_thin_max);
                                double del_h_thin = del_c*(h0+h_thin_max)/2.;
                                double del_hs_thin = del_c*M_hs_thin[i]/M_conc_thin[i];
     
                                M_thick[i] += del_h_thin;
                                // M_conc[i]  += del_c; ; <- this is done properly below

                                newice  = del_h_thin; // Reset newice to use below
                                newsnow = del_hs_thin;
                                // M_snow_thick[i] += newsnow; <- this is done properly below

                                M_conc_thin[i] -= del_c;
                                M_h_thin[i]    -= del_h_thin;
                                M_hs_thin[i]   -= del_hs_thin;
                            }    
                        }    
                    }
                    else // we should not have thin ice, no space for it
                    {
                        M_thick[i] += M_h_thin[i];
                        
                        newice  = M_h_thin[i];
                        newsnow = M_hs_thin[i];
                        
                        M_h_thin[i] = 0.;
                        M_hs_thin[i]= 0.;
                    }    
                    break;
                default:
                    std::cout << "newice_type = " << newice_type << "\n";
                    throw std::logic_error("Wrong newice_type");
            }
            /* Check bounds on del_c */
            del_c = std::min( 1.-M_conc[i], del_c );
        
        if ( del_hi < 0. )
        {
            /* Melting conditions */
            switch ( melt_type )
            {
                case 1:
                    /* Hibler (79) using PhiM to tune. PhiM = 0.5 is
                     * equivalent to Hibler's (79) approach */
                    if ( M_conc[i] < 1. )
                        del_c += del_hi*M_conc[i]*PhiM/hi_old;
                    else
                        del_c += 0.;
                    break;
                case 2:
                    /* Mellor and Kantha (89) */
                    /* Use the fraction PhiM of (1-c)*Qow to melt laterally */
                    del_c += PhiM*(1.-M_conc[i])*std::min(0.,Qow)*time_step/( hi*qi+hs*qs );
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
#if 1
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
#endif
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
        rain = (1.-old_conc-old_conc_thin)*M_precip[i] + (old_conc+old_conc_thin)*(M_precip[i]-tmp_snowfall);
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

        // Newly formed ice is undamaged - and unridged
        // calculate damage and ridge ratio as a weighted
        // average of the old damage - ridge ratio and 0, weighted with volume.
        if ( M_thick[i] > old_vol )
        {
            M_damage[i] = M_damage[i]*old_vol/M_thick[i];
            M_ridge_ratio[i] = M_ridge_ratio[i]*old_vol/M_thick[i];
        }

        // Set time_relaxation_damage to be inversely proportional to
        // temperature difference between bottom and snow-ice interface
#if 0
        if ( M_thick[i] > 0. )
        {
            deltaT = std::max(0., -physical::mu*M_sss[i] - M_tice[0][i] )
                / ( 1. + physical::ki*M_snow_thick[i]/(physical::ks*M_thick[i]) );
            M_time_relaxation_damage[i] = std::max(time_relaxation_damage*deltaT_relaxation_damage/deltaT, time_step);
        } else {
            M_time_relaxation_damage[i] = 1e36;
        }
#endif
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
FiniteElement::iceOceanHeatflux(int cpt, double sst, double sss, double mld, double dt)
{
    /* Use all excess heat to melt or grow ice. This is not
     * accurate, but will have to do for now! */
    double const Tbot = -physical::mu*sss; // Temperature at ice base (bottom), also freezing point of sea-water
    if ( vm["simul.Qio-type"].as<std::string>() == "basic" )
    {
        return (sst-Tbot)*physical::rhow*physical::cpw*mld/dt;
    } else if ( vm["simul.Qio-type"].as<std::string>() == "exchange" ) {
        double welt_oce_ice = 0.;
        for (int i=0; i<3; ++i)
        {
            int nind = (M_elements[cpt]).indices[i]-1;
            welt_oce_ice += std::hypot(M_VT[nind]-M_ocean[nind],M_VT[nind+M_num_nodes]-M_ocean[nind+M_num_nodes]);
        }
        double norm_Voce_ice = welt_oce_ice/3.;
        double Csens_io = 1e-3;
        return (sst-Tbot)*norm_Voce_ice*Csens_io*physical::rhow*physical::cpw;
    } else {
        std::cout << "Qio-type = " << vm["simul.Qio-type"].as<std::string>() << "\n";
        throw std::logic_error("Wrong Qio-type");
    }
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
FiniteElement::thermoWinton(int i, double dt, double wspeed, double sphuma, double conc, double voli, double vols,
        double Qlw_in, double Qsw_in, double mld, double snowfall,
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
    if ( conc <=0. || voli<=0.)
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

        // Snowfall [kg/m^2/s]
        hs += snowfall/physical::rhos*dt;

        // Sublimation at the surface
        if ( subl*dt <= hs*physical::rhos)
            hs -= subl*dt/physical::rhos;
        else if ( subl*dt - hs*physical::rhos <= h1*physical::rhoi )
        {
            h1 -= (subl*dt - hs*physical::rhos)/physical::rhoi;
            hs  = 0.;
        }
        else if ( subl*dt - h1*physical::rhoi - hs*physical::rhos <= h2*physical::rhoi )
        {
            h2 -= (subl*dt - h1*physical::rhoi - hs*physical::rhos)/physical::rhoi;
            h1  = 0.;
            hs  = 0.;
        }
        else
        {
            double ocn_evap_err = ( subl*dt - (h1+h2)*physical::rhoi - hs*physical::rhos )/physical::rhow;
			LOG(WARNING) << "All the ice has sublimated. This shouldn't happen and will result in lack of evaporation from the ocean of "
                << ocn_evap_err*1e3 << " mm over the current time step, in element " << i << ".\n";
            h2 = 0.;
            h1 = 0.;
            hs = 0.;
        }

        // Bottom melt/freezing
        Qio    = FiniteElement::iceOceanHeatflux(i, M_sst[i], M_sss[i], mld, dt);
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
            if ( h2+h1+hs -delh2-delh1-delhs <= 0. )
                Qio -= std::max(Mbot*dt - qs*hs + E1*h1 + E2*h2, 0.)/dt; // (34) - with added multiplication of rhoi and rhos and division with dt

            hs += delhs;
            h1 += delh1;
            h2 += delh2;
        }

        // Melting at the surface
        assert(Msurf >= 0); // Sanity check
        double delhs = -std::min(             Msurf*dt/qs,                          hs); // (27) - with division of rhos
        double delh1 = -std::min(std::max( -( Msurf*dt - qs*hs )/E1,           0.), h1); // (28) - with division of rhoi and rhos
        double delh2 = -std::min(std::max( -( Msurf*dt - qs*hs + E1*h1 ) / E2, 0.), h2); // (29) - with division of rhoi and rhos

        // If everyting melts we need to give back to the ocean
        if ( h2+h1+hs -delh2-delh1-delhs <= 0. )
            Qio -= std::max(Msurf*dt - qs*hs + E1*h1 + E2*h2, 0.)/dt; // (30) - with multiplication of rhoi and rhos and division with dt

        hs += delhs;
        h1 += delh1;
        h2 += delh2;
        hi  = h1 + h2;

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
FiniteElement::thermoIce0(int i, double wspeed, double sphuma, double conc, double voli, double vols, double Qlw_in, double Qsw_in, double mld, double snowfall,
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
    if ( conc <=0. || voli<=0.)
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
        double albedo = this->albedo(alb_scheme, Tsurf, hs, alb_sn, alb_ice, I_0);
        double dtsurf   = 1.;
        double Tbot     = -physical::mu*M_sss[i];
        int nb_iter_while=0;
        while ( dtsurf > 1e-4 )
        {
            nb_iter_while++;

            /* Calculate atmospheric fluxes */
            // Shortwave is modulated by the albedo
            Qsw = -Qsw_in*(1.-albedo)*(1.-I_0);
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
            //throw std::logic_error("nb_iter_while larger than 10");
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
        // snowfall in kg/m^2/s
        hs  = hs + del_hs + snowfall/physical::rhos*time_step;

        /* Heatflux from ocean */
        Qio = FiniteElement::iceOceanHeatflux(i, M_sst[i], M_sss[i], mld, time_step);
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
            del_hi  = -hi_old; //del_hi-hi;
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
#ifdef WITHGPERFTOOLS
    ProfilerStart("profile.log");
#endif

    std::string current_time_system = current_time_local();

    int pcpt = this->init();
    int niter = vm["simul.maxiteration"].as<int>();

    this->writeLogFile();

    // Debug file that records the time step
    std::fstream pcpt_file;
    pcpt_file.open(M_export_path + "/Timestamp.txt", std::ios::out | std::ios::trunc);

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

        //std::cout<< "pcpt= " << pcpt  <<"\n";

        pcpt_file << pcpt << "\n";
        pcpt_file << to_date_string(current_time) << "\n";
        pcpt_file.seekp(0);
    }

    pcpt_file.close();

    if ( pcpt*time_step/output_time_step < 1000 )
        this->exportResults(1000);
    LOG(INFO) <<"TIMER total = " << chrono_tot.elapsed() <<"s\n";
    LOG(INFO) <<"nb regrid total = " << M_nb_regrid <<"\n";

#ifdef WITHGPERFTOOLS
    //ProfilerFlush();
    ProfilerStop();
#endif

    this->finalise();

    LOG(INFO) << "-----------------------Simulation done on "<< current_time_local() <<"\n";
}

// Finalise everything
void
FiniteElement::finalise()
{
    // Don't forget to close the iabp file!
    if (M_use_iabp_drifters)
    {
        M_iabp_file.close();
        M_iabp_out.close();
    }

#ifdef OASIS
    int ierror = OASIS3::terminate();
#endif

    this->clear();
}

// Initialise everything
int
FiniteElement::init()
{
    // Initialise everything that doesn't depend on the mesh (constants, data set description, and time)

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
        //wim_grid = wim.wimGrid("km");
        wim_grid = wim.wimGrid("m");

        //total number of grid cells
        num_elements_wim_grid   = wim_grid.nx*wim_grid.ny;

        //range of x,y
        xmin_wim = *std::min_element(wim_grid.X.begin(),wim_grid.X.end());
        xmax_wim = *std::max_element(wim_grid.X.begin(),wim_grid.X.end());
        ymin_wim = *std::min_element(wim_grid.Y.begin(),wim_grid.Y.end());
        ymax_wim = *std::max_element(wim_grid.Y.begin(),wim_grid.Y.end());
        std::cout<<"xmin (WIM grid) = "<<xmin_wim<<"\n";
        std::cout<<"xmax (WIM grid) = "<<xmax_wim<<"\n";
        std::cout<<"ymin (WIM grid) = "<<ymin_wim<<"\n";
        std::cout<<"ymax (WIM grid) = "<<ymax_wim<<"\n";
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
        
//        for (int i=0; i<M_num_elements; i++)
//            M_damage[i]=(M_damage[i]>0.95 ? 1. : 0.);
        
        if(fmod(pcpt*time_step,output_time_step) == 0)
        {
            LOG(DEBUG) <<"export starts\n";
            this->exportResults((int) pcpt*time_step/output_time_step);
            LOG(DEBUG) <<"export done in " << chrono.elapsed() <<"s\n";
        }
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
    }

    LOG(DEBUG) <<"Initialize forcingAtmosphere\n";
    this->forcingAtmosphere();

    LOG(DEBUG) <<"Initialize forcingOcean\n";
    this->forcingOcean();

#if defined (WAVES)
    if (M_use_wim)
    {
        LOG(DEBUG) <<"Initialize forcingWave\n";
        this->forcingWave();
    }
#endif

    LOG(DEBUG) <<"Initialize bathymetry\n";
    this->bathymetry();

    this->checkReloadDatasets(M_external_data,current_time,"init - time-dependant");
#if 0
    chrono.restart();
    LOG(DEBUG) <<"check_and_reload starts\n";
    for ( auto it = M_external_data.begin(); it != M_external_data.end(); ++it )
        (*it)->check_and_reload(M_mesh,current_time);
    LOG(DEBUG) <<"check_and_reload in "<< chrono.elapsed() <<"s\n";
#endif

    if ( ! M_use_restart )
    {
        chrono.restart();
        this->initModelState();
        LOG(DEBUG) <<"initSimulation done in "<< chrono.elapsed() <<"s\n";
    }

    // Open the output file for drifters
    // TODO: Is this the right place to open the file?
    if (M_use_iabp_drifters )
    {
        // We should tag the file name with the init time in case of a re-start.
        std::stringstream filename;
        filename << M_export_path << "/drifters_out_" << current_time << ".txt";
        M_iabp_out.open(filename.str(), std::fstream::out);
        if ( ! M_iabp_out.good() )
            throw std::runtime_error("Cannot write to file: " + filename.str());
    }

    // Initialise the moorings - if requested
    if ( M_use_moorings )
        this->initMoorings();

#ifdef OASIS
    //!!!!!!!!!!!!!!!!! OASIS_INIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    int ierror, rank;
    int  mype, npes;            // rank and  number of pe
    MPI_Comm localComm;         // local MPI communicator and Initialized
    int comp_id;                // Component identification
    const char * comp_name = "nxtsim";  // Component name (6 characters) same as in the namcouple

    ierror = OASIS3::init_comp(&comp_id, comp_name);
    if (ierror != 0) {
        std::cout << "oasis_init_comp abort by nextsim with error code " << ierror << std::endl;
        OASIS3::abort(comp_id, comp_name, "Problem calling OASIS3::init_comp");
    }

    // Unit for output messages : one file for each process
    ierror =  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    //!!!!!!!!!!!!!!!!! OASIS_GET_LOCALCOMM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ierror = OASIS3::get_localcomm(&localComm);
    if (ierror != 0) {
        std::cout << "oasis_get_localcomm abort by nextsim with error code " << ierror << std::endl;
        OASIS3::abort(comp_id, comp_name, "Problem calling OASIS3::get_localcomm");
    }


    // Get MPI size and rank
    ierror = MPI_Comm_size(localComm, &npes);
    ierror = MPI_Comm_rank (localComm, &mype);
    std::cout << "npes : " << npes << "  mype : " << mype << std::endl;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //  GRID DEFINITION
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Output and averaging grids
    std::vector<double> data_nodes(M_num_nodes);
    std::vector<double> data_elements(M_num_elements);
    std::vector<double> data_grid;

    // Output variables - elements
    GridOutput::Variable conc(GridOutput::variableID::conc, data_elements, data_grid);

    GridOutput::Variable thick(GridOutput::variableID::thick, data_elements, data_grid);

    GridOutput::Variable snow(GridOutput::variableID::snow, data_elements, data_grid);

    std::vector<GridOutput::Variable> elemental_variables(3);
    elemental_variables[0] = conc;
    elemental_variables[1] = thick;
    elemental_variables[2] = snow;

    // Calculate the grid spacing (assuming a regular grid for now)
    auto RX = M_mesh.coordX();
    auto RY = M_mesh.coordY();
    auto xcoords = std::minmax_element( RX.begin(), RX.end() );
    auto ycoords = std::minmax_element( RY.begin(), RY.end() );

    double mooring_spacing = 1e3 * vm["simul.mooring_spacing"].as<double>();
    int nrows = (int) ( 0.5 + ( *xcoords.second - *xcoords.first )/mooring_spacing );
    int ncols = (int) ( 0.5 + ( *ycoords.second - *ycoords.first )/mooring_spacing );

    // Define the mooring dataset
    M_cpl_out = GridOutput(ncols, nrows, mooring_spacing, *xcoords.first, *ycoords.first, elemental_variables, GridOutput::variableKind::elemental);
    std::vector<int> lsm = M_cpl_out.getMask(M_mesh, GridOutput::variableKind::elemental);

    /*
    // Reading global grid netcdf file

    // Reading dimensions of the global grid
    read_dimgrid_c_(&nlon,&nlat,data_filename,&data_filename_len,&w_unit);

    cout << "nlon : " << nlon << "  nlat : " << nlat << endl;

    globalgrid_lon = (double*) malloc(nlon*nlat*sizeof(double));
    globalgrid_lat = (double*) malloc(nlon*nlat*sizeof(double));
    globalgrid_clo = (double*) malloc(nlon*nlat*nc*sizeof(double));
    globalgrid_cla = (double*) malloc(nlon*nlat*nc*sizeof(double));
    globalgrid_srf = (double*) malloc(nlon*nlat*sizeof(double));
    indice_mask = (int*) malloc(nlon*nlat*sizeof(int));

    // Reading of the longitudes, latitudes, longitude and latitudes of the corners, mask of the global grid
    read_grid_c_(&nlon,&nlat,&nc, data_filename, &data_filename_len, &w_unit,
        globalgrid_lon,globalgrid_lat,
        globalgrid_clo,globalgrid_cla,
        globalgrid_srf,
        indice_mask);

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //  PARTITION DEFINITION
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Definition of the partition of the grid (calling oasis_def_partition)
    ntot=nlon*nlat;

#ifdef DECOMP_APPLE
    il_paral_size = 3;
#elif defined DECOMP_BOX
    il_paral_size = 5;
#endif

    il_paral = (int*) malloc(il_paral_size*sizeof(int));

    decomp_def_(il_paral,&il_paral_size,&nlon,&nlat,&mype,&npes,&w_unit);
    //! The data are exchanged in the global grid so you do not need to pass
    //! isize to oasis_def_partition
    ierror = OASIS3::def_partition(&part_id, il_paral, 3);
    if (ierror != 0) {
        std::cout << "oasis_def_partition abort by nextsim with error code " << ierror << std::endl;
        OASIS3::abort(comp_id, comp_name, "Problem calling OASIS3::def_partition");
    }

    */
    int part_id;                    // partition id
    int ig_paral[3];
    ig_paral[0] = 0;                // a serial partition
    ig_paral[1] = 0;
    ig_paral[2] = M_num_elements;   // the total grid size
    ierror = OASIS3::def_partition(&part_id, ig_paral, (int) sizeof(ig_paral));
    if (ierror != 0) {
        std::cout << "oasis_def_partition abort by nextsim with error code " << ierror << std::endl;
        OASIS3::abort(comp_id, comp_name, "Problem calling OASIS3::def_partition");
    }

    // (Global) grid definition for OASIS3
    // Writing of the file grids.nc and masks.nc by the processor 0 from the grid read in

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //  GRID WRITING
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    if (mype == 0) {
        OASIS3::start_grids_writing(ierror);
        OASIS3::write_grid("nxts", ncols, nrows, &M_cpl_out.M_grid.gridLON[0], &M_cpl_out.M_grid.gridLAT[0]);
        // OASIS3::write_corner("nxts", ncols, nrows, 4, globalgrid_clo, globalgrid_cla);
        // OASIS3::write_area("nxts", ncols, nrows, globalgrid_srf);
        OASIS3::write_mask("nxts", ncols, nrows, &lsm[0]);
        OASIS3::terminate_grids_writing();
    }
    std::cout << "model1_cpp ====> After grids writing" << std::endl;


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // DEFINITION OF THE LOCAL FIELDS
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //!!!!!!!!!!!!!!! !!!!!!!!! OASIS_DEF_VAR !!!!!!!!!!!!!!!!!!!!!!!!!!!

    //!  Define transient variables

    // ---------------------"12345678"
    const char var_snd1[] = "FSendCnc"; // !8 characters field sent by model1 to model2
    const char var_rcv1[] = "FRecvSST"; //! 8 characters field received by model1 from model2
    //int var_id[2];
    int var_nodims[2];
    var_nodims[0] = 2 ; //   ! Rank of the field array is 2
    var_nodims[1] = 1 ; //   ! Bundles always 1 for OASIS3
    //var_type = OASIS3::OASIS_Real;
    int var_type = OASIS3::OASIS_Double;

    int var_actual_shape[4]; // ! local dimensions of the arrays to the pe. 2 x field rank (= 4 because fields are of rank = 2)
    var_actual_shape[0] = 1;
    var_actual_shape[1] = ncols;
    var_actual_shape[2] = 1;
    var_actual_shape[3] = nrows;

    // Declaration of the field associated with the partition
    var_id.push_back(-1);
    ierror = OASIS3::def_var(&var_id[0],var_snd1, part_id,
        var_nodims, OASIS3::OASIS_Out, var_actual_shape, var_type);
    cout << "model1_cpp ====> After def_var 1 " << endl;

    var_id.push_back(-1);
    ierror = OASIS3::def_var(&var_id[1],var_rcv1, part_id,
        var_nodims, OASIS3::OASIS_In, var_actual_shape, var_type);

    cout << "model1_cpp ====> After def_var 2 " << endl;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //         TERMINATION OF DEFINITION PHASE
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //  All processes involved in the coupling must call oasis_enddef;
    //  here all processes are involved in coupling

    //!!!!!!!!!!!!!!!!!! OASIS_ENDDEF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    OASIS3::enddef();
    cout << "======>After oasis_enddef" << endl;
#endif

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

    // Update the drifters position twice a day, important to keep the same frequency as the IABP data, for the moment
    if( pcpt==0 || std::fmod(current_time,0.5)==0 )
    {   
        // Read in the new buoys and output
        if ( M_use_iabp_drifters )
        {
            this->updateIABPDrifter();
            
            chrono.restart();
            LOG(DEBUG) <<"Drifter starts\n";
            LOG(DEBUG) <<"DRIFTER: Interp starts\n";

            // Assemble the coordinates from the unordered_map
            std::vector<double> drifter_X(M_iabp_drifters.size());
            std::vector<double> drifter_Y(M_iabp_drifters.size());
            int j=0;
            for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); ++it )
            {
                drifter_X[j] = it->second[0];
                drifter_Y[j] = it->second[1];
                ++j;
            }

            // Interpolate the total displacement and concentration onto the drifter positions
            int nb_var=2;
            std::vector<double> interp_drifter_in(nb_var*M_mesh.numNodes());

            // Interpolate the velocity
            for (int i=0; i<M_mesh.numNodes(); ++i)
            {
                interp_drifter_in[nb_var*i]   = M_UT[i];
                interp_drifter_in[nb_var*i+1] = M_UT[i+M_mesh.numNodes()];
            }

            double* interp_drifter_out;
            InterpFromMeshToMesh2dx(&interp_drifter_out,
                                    &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                                    M_mesh.numNodes(),M_mesh.numTriangles(),
                                    &interp_drifter_in[0],
                                    M_mesh.numNodes(),nb_var,
                                    &drifter_X[0],&drifter_Y[0],M_iabp_drifters.size(),
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
                                    &drifter_X[0],&drifter_Y[0],M_iabp_drifters.size(),
                                    true, 0.);

            // Rebuild the M_iabp_drifters map
            double clim = vm["simul.drifter_climit"].as<double>();
            j=0;
            for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); /* ++it is not allowed here, because we use 'erase' */ )
            {
                if ( interp_drifter_c_out[j] > clim )
                {
                    M_iabp_drifters[it->first] = std::array<double,2> {it->second[0]+interp_drifter_out[nb_var*j], it->second[1]+interp_drifter_out[nb_var*j+1]};
                    ++it;
                } else {
                    // Throw out drifters that drift out of the ice
                    it = M_iabp_drifters.erase(it);
                }
                ++j;
            }

            xDelete<double>(interp_drifter_out);
            xDelete<double>(interp_drifter_c_out);

            LOG(DEBUG) <<"DRIFTER: Interp done\n";
            LOG(DEBUG) <<"Drifter interp done in "<< chrono.elapsed() <<"s\n";

            // TODO: Do we want to output drifters at a different time interval?
            this->outputDrifter(M_iabp_out);
        }
        
        if ( M_use_equallyspaced_drifters )
            M_equallyspaced_drifters.move(M_mesh, M_UT);
        if ( M_use_rgps_drifters )
            M_rgps_drifters.move(M_mesh, M_UT);
        if ( M_use_osisaf_drifters )
            for (auto it=M_osisaf_drifters.begin(); it!=M_osisaf_drifters.end(); it++)
                it->move(M_mesh, M_UT);
    
	    for (int i=0; i<M_num_nodes; ++i)
	    {
	        // UM
		    M_UT[i] = 0.;
		    M_UT[i+M_num_nodes] = 0.;
	    }
    }
    if(pcpt>0)
    {
        if ( M_use_equallyspaced_drifters && fmod(current_time,M_equallyspaced_drifters_output_time_step) == 0 )
            M_equallyspaced_drifters.appendNetCDF(current_time, M_mesh, M_UT);

        if ( M_use_rgps_drifters )
        {
            std::string time_str = vm["simul.RGPS_time_init"].as<std::string>();
            double RGPS_time_init = from_date_time_string(time_str);
        
            if( !M_rgps_drifters.isInitialised() && current_time == RGPS_time_init)
                this->updateRGPSDrifters();
            
            if( current_time != RGPS_time_init && fmod(current_time,M_rgps_drifters_output_time_step) == 0 )
                if ( M_rgps_drifters.isInitialised() )
                    M_rgps_drifters.appendNetCDF(current_time, M_mesh, M_UT);
        }
    }
     
    if ( M_use_osisaf_drifters && fmod(current_time+0.5,1.) == 0 )
    {
        // OSISAF drift is calculated as a dirfter displacement over 48 hours
        // and they have two sets of drifters in the field at all times.

        // Write out the contents of [1] if it's meaningfull
        if ( M_osisaf_drifters[1].isInitialised() )
            M_osisaf_drifters[1].appendNetCDF(current_time, M_mesh, M_UT);

        // Flip the vector so we move [0] to be [1]
        std::reverse(M_osisaf_drifters.begin(), M_osisaf_drifters.end());

        // Create a new M_drifters instance in [0], with a properly initialised netCDF file
        M_osisaf_drifters[0] = Drifters("data", "ice_drift_nh_polstere-625_multi-oi.nc", "yc", "yx", "lat", "lon", M_mesh, M_conc, vm["simul.drifter_climit"].as<double>());
        M_osisaf_drifters[0].initNetCDF(M_export_path+"/OSISAF_", current_time);
        M_osisaf_drifters[0].appendNetCDF(current_time, M_mesh, M_UT);
    }

    // remeshing and remapping of the prognostic variables

    // The first time step we behave as if we just did a regrid
    M_regrid = (pcpt==0);

    if (vm["simul.regrid"].as<std::string>() == "bamg")
    {
        double displacement_factor = 1.;
        double minang = this->minAngle(M_mesh,M_UM,displacement_factor);
        LOG(DEBUG) <<"REGRID ANGLE= "<< minang <<"\n";

        bool flip_test = this->flip(M_mesh,M_UM,displacement_factor);
        LOG(DEBUG) <<"FLIP TEST= "<< flip_test <<"\n";


        if ( minang < vm["simul.regrid_angle"].as<double>() || flip_test )
        {
            M_regrid = true;

            // this->writeRestart(pcpt, 0); // Write a restart before regrid - useful for debugging
            if ( M_use_moorings && ! M_moorings_snapshot )
                M_moorings.updateGridMean(M_mesh);
#ifdef OASIS
            M_cpl_out.updateGridMean(M_mesh);
#endif

            LOG(DEBUG) <<"Regridding starts\n";
            chrono.restart();
            this->regrid(pcpt);
            LOG(DEBUG) <<"Regridding done in "<< chrono.elapsed() <<"s\n";
            if ( M_use_moorings )
                M_moorings.resetMeshMean(M_mesh);
#ifdef OASIS
            M_cpl_out.resetMeshMean(M_mesh);
#endif
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
    this->checkReloadDatasets(M_external_data,current_time+time_step/(24*3600.0),
            "step - time-dependant");
#if 0
    LOG(DEBUG) <<"check_and_reload starts\n";
    for ( auto it = M_external_data.begin(); it != M_external_data.end(); ++it )
        (*it)->check_and_reload(M_mesh,current_time+time_step/(24*3600.0));
    LOG(DEBUG) <<"check_and_reload in "<< chrono.elapsed() <<"s\n";
#endif

    M_use_restart = false;

#if 1
    if (pcpt == 0)
    {
        chrono.restart();
        LOG(DEBUG) <<"first export starts\n";
        this->exportResults(0);
        this->writeRestart(pcpt, 0); // Write a restart before regrid - useful for debugging
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

    chrono.restart();
    this->solve();
    LOG(INFO) <<"TIMER SOLUTION= " << chrono.elapsed() <<"s\n";

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
        // If we're taking snapshots the we only call updateMeans before writing to file
        if ( ! M_moorings_snapshot )
            this->updateMeans(M_moorings, mooring_time_factor);

        if ( fmod(pcpt*time_step,mooring_output_time_step) == 0 )
        {
            if ( M_moorings_snapshot )
                // Update the snapshot
                this->updateMeans(M_moorings, 1.);

            M_moorings.updateGridMean(M_mesh);

            if ( M_moorings_snapshot )
                M_moorings.appendNetCDF(M_moorings_file, current_time);
            else
                // shift the timestamp in the file to the centre of the output interval
                M_moorings.appendNetCDF(M_moorings_file, current_time-mooring_output_time_step/86400/2);

            M_moorings.resetMeshMean(M_mesh);
            M_moorings.resetGridMean();

            double not_used;
            if ( M_moorings_file_length != GridOutput::fileLength::inf && modf(current_time, &not_used) < time_step*86400 )
            // It's a new day, so we check if we need a new file
            {
                boost::gregorian::date now = Nextsim::parse_date(current_time);
                switch (M_moorings_file_length)
                {
                    case GridOutput::fileLength::daily:
                        M_moorings_file = M_moorings.initNetCDF(M_export_path + "/Moorings", M_moorings_file_length, current_time);
                        break;
                    case GridOutput::fileLength::weekly:
                        if ( now.day_of_week().as_number() == 1 )
                            M_moorings_file = M_moorings.initNetCDF(M_export_path + "/Moorings", M_moorings_file_length, current_time);
                        break;
                    case GridOutput::fileLength::monthly:
                        if ( now.day().as_number() == 1 )
                            M_moorings_file = M_moorings.initNetCDF(M_export_path + "/Moorings", M_moorings_file_length, current_time);
                        break;
                    case GridOutput::fileLength::yearly:
                        if ( now.day_of_year() == 1 )
                            M_moorings_file = M_moorings.initNetCDF(M_export_path + "/Moorings", M_moorings_file_length, current_time);
                }
            }
        }
    }

#ifdef OASIS
    this->updateMeans(M_cpl_out, cpl_time_factor);
    if ( fmod(pcpt*time_step,cpl_time_step) == 0 )
    {
        M_cpl_out.updateGridMean(M_mesh);

        int ierror = OASIS3::put_2d(var_id[0], pcpt*time_step, &M_cpl_out.M_elemental_variables[0].data_mesh[0], M_cpl_out.M_ncols, M_cpl_out.M_nrows);

        M_cpl_out.resetMeshMean(M_mesh);
        M_cpl_out.resetGridMean();
    }
#endif


#endif

    if ( fmod(pcpt*time_step,restart_time_step) == 0)
    {
        std::cout << "Writing restart file after time step " <<  pcpt-1 << "\n";
        this->writeRestart(pcpt, (int) pcpt*time_step/restart_time_step);
    }
}

// Add to the mean on the mesh
void
FiniteElement::updateMeans(GridOutput &means, double time_factor)
{
    // Update elements and multiply with time_factor
    for ( auto it=means.M_elemental_variables.begin(); it!=means.M_elemental_variables.end(); ++it )
    {
        switch (it->varID)
        {
            case (GridOutput::variableID::conc):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_conc[i]*time_factor;
                break;

            case (GridOutput::variableID::thick):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_thick[i]*time_factor;
                break;

            case (GridOutput::variableID::damage):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_damage[i]*time_factor;
                break;

            case (GridOutput::variableID::snow):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_snow_thick[i]*time_factor;
                break;

            case (GridOutput::variableID::tsurf):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += ( M_conc[i]*M_tice[0][i] + (1-M_conc[i])*M_sst[i] )*time_factor;
                break;

            case (GridOutput::variableID::sst):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_sst[i]*time_factor;
                break;

            case (GridOutput::variableID::sss):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_sss[i]*time_factor;
                break;

            case (GridOutput::variableID::tsurf_ice):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_tice[0][i]*time_factor;
                break;

            case (GridOutput::variableID::t1):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_tice[1][i]*time_factor;
                break;

            case (GridOutput::variableID::t2):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_tice[2][i]*time_factor;
                break;

            case (GridOutput::variableID::conc_thin):
                for (int i=0; i<M_num_elements; i++)
                {
                    it->data_mesh[i] += M_conc_thin[i]*time_factor;
                }
                break;

            case (GridOutput::variableID::h_thin):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_h_thin[i]*time_factor;
                break;

            case (GridOutput::variableID::hs_thin):
                for (int i=0; i<M_num_elements; i++)
                    it->data_mesh[i] += M_hs_thin[i]*time_factor;
                break;

            default: std::logic_error("Updating of given variableID not implimented (elements)");
        }
    }

    // Update nodes
    for ( auto it=means.M_nodal_variables.begin(); it!=means.M_nodal_variables.end(); ++it )
    {
        switch (it->varID)
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
    // Output and averaging grids
    std::vector<double> data_nodes(M_num_nodes);
    std::vector<double> data_elements(M_num_elements);
    std::vector<double> data_grid;

    // Output variables - elements
    GridOutput::Variable conc(GridOutput::variableID::conc, data_elements, data_grid);
    GridOutput::Variable thick(GridOutput::variableID::thick, data_elements, data_grid);
    GridOutput::Variable snow(GridOutput::variableID::snow, data_elements, data_grid);

    std::vector<GridOutput::Variable> elemental_variables(3);
    elemental_variables[0] = conc;
    elemental_variables[1] = thick;
    elemental_variables[2] = snow;
    if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
    {
        GridOutput::Variable conc_thin(GridOutput::variableID::conc_thin, data_elements, data_grid);
        GridOutput::Variable h_thin(GridOutput::variableID::h_thin, data_elements, data_grid);
        GridOutput::Variable hs_thin(GridOutput::variableID::hs_thin, data_elements, data_grid);

        elemental_variables.push_back(conc_thin);
        elemental_variables.push_back(h_thin);
        elemental_variables.push_back(hs_thin);
    }

    // Output variables - nodes
    GridOutput::Variable siu(GridOutput::variableID::VT_x, data_nodes, data_grid);

    GridOutput::Variable siv(GridOutput::variableID::VT_y, data_nodes, data_grid);

    std::vector<GridOutput::Variable> nodal_variables(2);
    nodal_variables[0] = siu;
    nodal_variables[1] = siv;

    // The vectorial variables are (always on the nodes) ...
    std::vector<int> siuv_id(2);
    siuv_id[0] = 0;
    siuv_id[1] = 1;

    GridOutput::Vectorial_Variable siuv{
        components_Id: siuv_id,
        // east_west_oriented: true
        east_west_oriented: false
    };

    std::vector<GridOutput::Vectorial_Variable> vectorial_variables(1);
    vectorial_variables[0] = siuv;

#if 1
    // Calculate the grid spacing (assuming a regular grid for now)
    auto RX = M_mesh.coordX();
    auto RY = M_mesh.coordY();
    auto xcoords = std::minmax_element( RX.begin(), RX.end() );
    auto ycoords = std::minmax_element( RY.begin(), RY.end() );

    double mooring_spacing = 1e3 * vm["simul.mooring_spacing"].as<double>();
    int ncols = (int) ( 0.5 + ( *xcoords.second - *xcoords.first )/mooring_spacing );
    int nrows = (int) ( 0.5 + ( *ycoords.second - *ycoords.first )/mooring_spacing );

    // Define the mooring dataset
    M_moorings = GridOutput(ncols, nrows, mooring_spacing, *xcoords.first, *ycoords.first, nodal_variables, elemental_variables, vectorial_variables);
#else
    // Read the grid in from file
    // Define a grid
    GridOutput::Grid grid{
        gridFile: "ice_drift_nh_polstere-625_multi-oi.nc",
        dirname: "data",
        mpp_file: Environment::vm()["simul.proj_filename"].as<std::string>(),
        dimNameX: "yc",
        dimNameY: "xc",
        latName: "lat",
        lonName: "lon"
    };

    // Define the mooring dataset
    M_moorings = GridOutput(grid, nodal_variables, elemental_variables, vectorial_variables);

    /* Just for debuging
    // Save the grid info - this is still just an ascii dump!
    std::ofstream myfile;
    myfile.open("lon_grid.dat");
    std::copy(M_moorings.M_grid.gridLON.begin(), M_moorings.M_grid.gridLON.end(), ostream_iterator<float>(myfile," "));
    myfile.close();
    myfile.open("lat_grid.dat");
    std::copy(M_moorings.M_grid.gridLAT.begin(), M_moorings.M_grid.gridLAT.end(), ostream_iterator<float>(myfile," "));
    myfile.close();

    myfile.open("x_grid.dat");
    std::copy(M_moorings.M_grid.gridX.begin(), M_moorings.M_grid.gridX.end(), ostream_iterator<float>(myfile," "));
    myfile.close();
    myfile.open("y_grid.dat");
    std::copy(M_moorings.M_grid.gridY.begin(), M_moorings.M_grid.gridY.end(), ostream_iterator<float>(myfile," "));
    myfile.close();
    */
#endif

    M_moorings_file = M_moorings.initNetCDF(M_export_path + "/Moorings", M_moorings_file_length, time_init);

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

    std::vector<int> misc_int(4);
    misc_int[0] = pcpt;
    misc_int[1] = M_flag_fix;
    misc_int[2] = current_time;
    misc_int[3] = mesh_adapt_step;
    exporter.writeField(outbin, misc_int, "Misc_int");
    exporter.writeField(outbin, M_dirichlet_flags, "M_dirichlet_flags");

    exporter.writeField(outbin, M_conc, "M_conc");
    exporter.writeField(outbin, M_thick, "M_thick");
    exporter.writeField(outbin, M_snow_thick, "M_snow_thick");
    exporter.writeField(outbin, M_sigma, "M_sigma");
    exporter.writeField(outbin, M_damage, "M_damage");
    exporter.writeField(outbin, M_ridge_ratio, "M_ridge_ratio");
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
    exporter.writeField(outbin, M_UT, "M_UT");

    if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
    {
        exporter.writeField(outbin, M_h_thin, "M_h_thin");
        exporter.writeField(outbin, M_conc_thin, "M_conc_thin");
        exporter.writeField(outbin, M_hs_thin, "M_hs_thin");
        exporter.writeField(outbin, M_tsurf_thin, "M_tsurf_thin");
    }

    if (M_use_iabp_drifters)
    {
        std::vector<int> drifter_no(M_iabp_drifters.size());
        std::vector<double> drifter_x(M_iabp_drifters.size());
        std::vector<double> drifter_y(M_iabp_drifters.size());

        // Sort the drifters so the restart files are identical (this is just to make testing easier)
        int j=0;
        for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); ++it )
        {
            drifter_no[j] = it->first;
            ++j;
        }
        std::sort(drifter_no.begin(), drifter_no.end());

        j=0;
        for ( auto it = drifter_no.begin(); it != drifter_no.end(); it++ )
        {
            drifter_x[j] = M_iabp_drifters[drifter_no[j]][0];
            drifter_y[j] = M_iabp_drifters[drifter_no[j]][1];
            j++;
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
    mesh_adapt_step = field_map_int["Misc_int"].at(3);

    std::vector<int> dirichlet_flags = field_map_int["M_dirichlet_flags"];
    for (int edg=0; edg<bamgmesh->EdgesSize[0]; ++edg)
    {
        int fnd = bamgmesh->Edges[3*edg]-1;
        if ((std::binary_search(dirichlet_flags.begin(),dirichlet_flags.end(),fnd)))
        {
            bamggeom->Edges[3*edg+2] = M_flag_fix;
            bamgmesh->Edges[3*edg+2] = M_flag_fix;
        }
        else
        {
            bamggeom->Edges[3*edg+2] = M_flag_fix+1; // we just want it to be different than M_flag_fix
            bamgmesh->Edges[3*edg+2] = M_flag_fix+1; // we just want it to be different than M_flag_fix        
        }
    }

    // Import the bamg structs
    this->importBamg(bamgmesh);

    // We mask out the boundary nodes
    M_mask.assign(bamgmesh->VerticesSize[0],false) ;
    M_mask_dirichlet.assign(bamgmesh->VerticesSize[0],false) ;
    for (int vert=0; vert<bamgmesh->VerticesOnGeomVertexSize[0]; ++vert)
        M_mask[bamgmesh->VerticesOnGeomVertex[2*vert]-1]=true; // The factor 2 is because VerticesOnGeomVertex has 2 dimensions in bamg

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
        M_mask_dirichlet[M_dirichlet_flags[i]]=true;
    }

    M_neumann_nodes.resize(2*(M_neumann_flags.size()));
    for (int i=0; i<M_neumann_flags.size(); ++i)
    {
        M_neumann_nodes[2*i] = M_neumann_flags[i];
        M_neumann_nodes[2*i+1] = M_neumann_flags[i]+M_num_nodes;
    }

    // === Set the prognostic variables ===
    M_conc       = field_map_dbl["M_conc"];
    M_thick      = field_map_dbl["M_thick"];
    M_snow_thick = field_map_dbl["M_snow_thick"];
    M_sigma      = field_map_dbl["M_sigma"];
    M_damage     = field_map_dbl["M_damage"];
    M_ridge_ratio     = field_map_dbl["M_ridge_ratio"];
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
    M_UT         = field_map_dbl["M_UT"];

    M_surface.assign(M_num_elements,0.);
    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        M_surface[cpt] = this->measure(*it,M_mesh,M_UM);
        ++cpt;
    }

    //for (int i=0; i < M_thick.size(); i++)
    //{
    //  M_thick[i] *= 2.0;
    //}

    if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
    {
        M_h_thin     = field_map_dbl["M_h_thin"];
        M_conc_thin     = field_map_dbl["M_conc_thin"];
        M_hs_thin    = field_map_dbl["M_hs_thin"];
        M_tsurf_thin = field_map_dbl["M_tsurf_thin"];
    }

    if (M_use_iabp_drifters)
    {
        std::vector<int>    drifter_no = field_map_int["Drifter_no"];
        std::vector<double> drifter_x  = field_map_dbl["Drifter_x"];
        std::vector<double> drifter_y  = field_map_dbl["Drifter_y"];

        this->initIABPDrifter();
        if (drifter_no.size() == 0)
        {
            LOG(WARNING) << "Warning: Couldn't read drifter positions from restart file. Drifter positions initialised as if there was no restart.\n";
            this->updateIABPDrifter();
        } else {
            for ( int i=0; i<drifter_no.size(); ++i )
            {
                M_iabp_drifters.emplace(drifter_no[i], std::array<double,2>{drifter_x[i], drifter_y[i]});
            }
        }
    }

    inbin.close();

#if 0
    //TW: these lines shouldn't be necessary now - set in check_and_reload
    // Set the target size for the data sets
    M_atmosphere_nodes_dataset.target_size=M_num_nodes;
    M_atmosphere_elements_dataset.target_size=M_num_elements;
    M_atmosphere_bis_elements_dataset.target_size=M_num_elements;
    M_ocean_nodes_dataset.target_size=M_num_nodes;
    M_ocean_elements_dataset.target_size=M_num_elements;

    M_ice_topaz_elements_dataset.target_size=M_num_elements;
    M_ice_icesat_elements_dataset.target_size=M_num_elements;
    M_ice_piomas_elements_dataset.target_size=M_num_elements;
    M_ice_amsre_elements_dataset.target_size=M_num_elements;
    M_ice_osisaf_elements_dataset.target_size=M_num_elements;
    M_ice_osisaf_type_elements_dataset.target_size=M_num_elements;
    M_ice_amsr2_elements_dataset.target_size=M_num_elements;
    M_ice_cs2_smos_elements_dataset.target_size=M_num_elements;
    M_ice_smos_elements_dataset.target_size=M_num_elements;
    M_bathymetry_elements_dataset.target_size=M_num_elements;
#if defined (WAVES)
    if (M_use_wim)
        M_wave_elements_dataset.target_size=M_num_elements;
#endif
#endif
    return pcpt;
}//readRestart

void
FiniteElement::updateVelocity()
{
    M_VTMM = M_VTM;
    M_VTM  = M_VT;
    //M_VT   = M_solution->container();
    M_solution->container(&M_VT[0]); // this version is more rapid as it directly acceses the memory and does not create a new vector 0.16 ms instead of 0.23 ms for container()

    // TODO (updateVelocity) Sylvain: This limitation cost about 1/10 of the solver time.
    // TODO (updateVelocity) Sylvain: We could add a term in the momentum equation to avoid the need of this limitation.
    //std::vector<double> speed_c_scaling_test(bamgmesh->NodalElementConnectivitySize[0]);
    #if 0
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

            M_snowfall=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
            M_external_data.push_back(&M_snowfall);
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
        break;

        case setup::AtmosphereType::CFSR_HI:
        case setup::AtmosphereType::CFSR:
            M_wind=ExternalData(
                &M_atmosphere_nodes_dataset,M_mesh,0 ,true ,
                time_init, vm["simul.spinup_duration"].as<double>());
            M_external_data.push_back(&M_wind);

            M_tair=ExternalData(&M_atmosphere_elements_dataset,M_mesh,0,false,time_init);
            M_external_data.push_back(&M_tair);

            M_sphuma=ExternalData(&M_atmosphere_elements_dataset,M_mesh,1,false,time_init);
            M_external_data.push_back(&M_sphuma);

            M_mslp=ExternalData(&M_atmosphere_elements_dataset,M_mesh,2,false,time_init);
            M_external_data.push_back(&M_mslp);

            M_Qsw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,3,false,time_init);
            M_external_data.push_back(&M_Qsw_in);

            M_Qlw_in=ExternalData(&M_atmosphere_elements_dataset,M_mesh,4,false,time_init);
            M_external_data.push_back(&M_Qlw_in);

            M_precip=ExternalData(&M_atmosphere_elements_dataset,M_mesh,5,false,time_init);
            M_external_data.push_back(&M_precip);

            M_snowfr=ExternalData(&M_atmosphere_elements_dataset,M_mesh,6,false,time_init);
            M_external_data.push_back(&M_snowfr);
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

            M_ocean_temp=ExternalData(physical::ocean_freezing_temp);
            M_external_data.push_back(&M_ocean_temp);

            M_ocean_salt=ExternalData(physical::ocean_freezing_temp/physical::mu);
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
        case setup::OceanType::TOPAZR_atrest:
            M_ocean=ExternalData(
                vm["simul.constant_ocean_v"].as<double>(),
                vm["simul.constant_ocean_v"].as<double>(),
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
    wim_ideal_forcing = true;
    double xedge = xmin_wim + 0.25*(xmax_wim-xmin_wim);

    switch (M_wave_type)
    {
        case setup::WaveType::SET_IN_WIM:
            //pass in empty vectors to wim.run()
            //- then waves are set in there
            M_SWH_grid.resize(0);
            M_MWP_grid.resize(0);
            M_MWD_grid.resize(0);
            break;

        case setup::WaveType::CONSTANT:
            //set arrays to pass in to wim.run()
            M_SWH_grid.assign(num_elements_wim_grid,vm["wim.hsinc" ].as<double>());
            M_MWP_grid.assign(num_elements_wim_grid,vm["wim.tpinc" ].as<double>());
            M_MWD_grid.assign(num_elements_wim_grid,vm["wim.mwdinc"].as<double>());
            // mwd is relative to the nextsim/WIM x-y coord system in this case
            break;


        case setup::WaveType::CONSTANT_PARTIAL:
            //set arrays to pass in to wim.run()
            std::cout<<"wim elements"<<num_elements_wim_grid<<"\n";
            M_SWH_grid.assign(num_elements_wim_grid,vm["wim.hsinc" ].as<double>());
            M_MWP_grid.assign(num_elements_wim_grid,vm["wim.tpinc" ].as<double>());
            M_MWD_grid.assign(num_elements_wim_grid,vm["wim.mwdinc"].as<double>());
            std::cout<<"wim elements"<<num_elements_wim_grid<<"\n";
            for (int i=0;i<num_elements_wim_grid;i++)
            {
                if(wim_grid.X[i]>=xedge)
                {
                    M_SWH_grid[i]   = 0.;
                    M_MWP_grid[i]   = 0.;
                    M_MWD_grid[i]   = 0.;
                    // mwd is relative to the nextsim/WIM x-y coord system in this case
                }
            }
            break;
        case setup::WaveType::WW3A:

            //initialise arrays to pass in to wim.run()
            M_SWH_grid.assign(num_elements_wim_grid,0);
            M_MWP_grid.assign(num_elements_wim_grid,0);
            M_MWD_grid.assign(num_elements_wim_grid,0);

            // define external_data objects
            M_SWH        = ExternalData(&M_wave_elements_dataset, M_mesh, 0,false,time_init);
            M_MWP        = ExternalData(&M_wave_elements_dataset, M_mesh, 1,false,time_init);
            M_MWD        = ExternalData(&M_wave_elements_dataset, M_mesh, 0,true,time_init);//now a vector
            M_fice_waves = ExternalData(&M_wave_elements_dataset, M_mesh, 4,false,time_init);

            // add them to a vector for looping
            M_external_data.push_back(&M_SWH);
            M_external_data.push_back(&M_MWP);
            M_external_data.push_back(&M_MWD);
            M_external_data.push_back(&M_fice_waves);

            wim_forcing_options = M_wave_elements_dataset.grid.waveOptions;
            wim_ideal_forcing   = false;

            break;

        case setup::WaveType::ERAI_WAVES_1DEG:

            //initialise arrays to pass in to wim.run()
            M_SWH_grid.assign(num_elements_wim_grid,0);
            M_MWP_grid.assign(num_elements_wim_grid,0);
            M_MWD_grid.assign(num_elements_wim_grid,0);

            // define external_data objects
            M_SWH = ExternalData(&M_wave_elements_dataset, M_mesh, 0,false,time_init);
            M_MWP = ExternalData(&M_wave_elements_dataset, M_mesh, 1,false,time_init);
            M_MWD = ExternalData(&M_wave_elements_dataset, M_mesh, 0,true,time_init);//now a vector

            // add them to a vector for looping
            M_external_data.push_back(&M_SWH);
            M_external_data.push_back(&M_MWP);
            M_external_data.push_back(&M_MWD);

            wim_forcing_options = M_wave_elements_dataset.grid.waveOptions;
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
            //std::fill(M_sst.begin(), M_sst.end(), -1.8);
            std::fill(M_sst.begin(), M_sst.end(), 1.);
            std::fill(M_sss.begin(), M_sss.end(),  1.8/physical::mu);
            break;
        case setup::OceanType::TOPAZR:
        case setup::OceanType::TOPAZR_atrest:
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
        case setup::IceType::BINARY:
            this->binaryIce();
            break;
        case setup::IceType::TOPAZ4:
            this->topazIce();
            break;
        case setup::IceType::TOPAZ4OSISAFICESAT:
            this->topazIceOsisafIcesat();
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
        case setup::IceType::SMOS:
            this->smosIce();
            break;
        default:
            std::cout << "invalid initialization of the ice"<<"\n";
            throw std::logic_error("invalid initialization of the ice");
    }

    double hi, conc_tot, weight_conc;
    // Consistency check for the slab ocean + initialization of the ice temperature (especially for Winton)
    for ( int i=0; i<M_num_elements; i++ )
    {  
        hi=0.; 
        if(M_conc[i]>0.)
            hi = M_thick[i]/M_conc[i];
        
        if ( M_conc[i] < physical::cmin || hi < physical::hmin)
        {
            M_conc[i]=0.;
            M_thick[i]=0.;
            M_snow_thick[i]=0.;
        }
        
        hi = 0.; 
        if(M_conc_thin[i]>0.)
            hi = M_h_thin[i]/M_conc_thin[i];

        if ( M_conc_thin[i] < physical::cmin || hi < physical::hmin)
        {
            M_conc_thin[i]=0.;
            M_h_thin[i]=0.;
            M_hs_thin[i]=0.;
        }

        conc_tot=M_conc[i]+M_conc_thin[i];
        weight_conc=std::min(1.,conc_tot*2.);
        if(conc_tot>0.)
            M_sst[i] = -M_sss[i]*physical::mu*weight_conc+M_sst[i]*(1.-weight_conc);
            //M_sst[i] = -M_sss[i]*physical::mu;//*M_conc[i]+M_ocean_temp[i]*(1.-M_conc[i]);
        
        if ( M_snow_thick[i] > 0. )
            M_tice[0][i] = std::min(0., M_tair[i]);
        else
            M_tice[0][i] = std::min(-physical::mu*physical::si, M_tair[i]);
    }

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
    /*
    double y_max=300000.;
    double y_min=150000.;
    double x_max1=240000.;
    double x_min1=150000.;
    double x_max2=300000.;
    double x_min2=260000.;
    */
 
    double y_max=350000.;
    double y_min=0000.;
    double x_max=400000.;
    double x_min=200000.;

    double transition=(x_max-x_min)/10.;

	double tmp_var;

    auto RX = M_mesh.bcoordX();
    auto RY = M_mesh.bcoordY();
    double cmin= 0.;
            
    for (int i=0; i<M_num_elements; ++i)
    {
        tmp_var = (RY[i]<=y_max)*(RY[i]>=y_min)*(RX[i]<=x_max)*(RX[i]>=x_min);
        
/*        tmp_var = (RY[i]<=y_max)*(RY[i]>=y_min)*(RX[i]<=x_max1)*(RX[i]>=x_min1)
            + (RY[i]<=y_max)*(RY[i]>=y_min)*(RX[i]<=x_max2)*(RX[i]>=x_min2);
*/
            /*
            +     (RY[i]<=y_max)*(RY[i]>=y_min)*(RX[i]<x_min)*std::max(cmin,(1.-std::hypot(RX[i]-x_min,0.         )/transition))
            +     (RY[i]<=y_max)*(RY[i]>=y_min)*(RX[i]>x_max)*std::max(cmin,(1.-std::hypot(RX[i]-x_max,0.         )/transition))
            +     (RX[i]<=x_max)*(RX[i]>=x_min)*(RY[i]<y_min)*std::max(cmin,(1.-std::hypot(RY[i]-y_min,0.         )/transition))
            +     (RX[i]<=x_max)*(RX[i]>=x_min)*(RY[i]>y_max)*std::max(cmin,(1.-std::hypot(RY[i]-y_max,0.         )/transition))
            +     (RY[i]< y_min)*(RX[i]< x_min)              *std::max(cmin,(1.-std::hypot(RX[i]-x_min,RY[i]-y_min)/transition))
            +     (RY[i]< y_min)*(RX[i]> x_max)              *std::max(cmin,(1.-std::hypot(RX[i]-x_max,RY[i]-y_min)/transition))
            +     (RY[i]> y_max)*(RX[i]< x_min)              *std::max(cmin,(1.-std::hypot(RX[i]-x_min,RY[i]-y_max)/transition))
            +     (RY[i]> y_max)*(RX[i]> x_max)              *std::max(cmin,(1.-std::hypot(RX[i]-x_max,RY[i]-y_max)/transition));
            */
        std::cout<<"RX: "<< RX[i] << "RY: "<< RY[i] << "tmp_var: " << tmp_var << "\n";
	
        //M_conc[i]  = std::max(vm["simul.init_concentration"].as<double>()*tmp_var,cmin);
	//	M_thick[i] = vm["simul.init_thickness"].as<double>()*M_conc[i];
	//	M_snow_thick[i] = vm["simul.init_snow_thickness"].as<double>()*M_conc[i];

        M_conc[i]  = 1.; //vm["simul.init_concentration"].as<double>();
	
	if(i==10)
		M_conc[i]=0.;
	
	M_thick[i] = vm["simul.init_thickness"].as<double>()*M_conc[i];
	M_snow_thick[i] = vm["simul.init_snow_thickness"].as<double>()*M_conc[i];
	
	M_damage[i]=0.;
        
        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_conc_thin[i]  = vm["simul.init_thin_conc"].as<double>();
            
            M_h_thin[i]     = (vm["simul.h_thin_min"].as<double>()+(vm["simul.h_thin_max"].as<double>()-vm["simul.h_thin_min"].as<double>())/2.)*M_conc_thin[i];
            
            M_hs_thin[i]    = vm["simul.init_snow_thickness"].as<double>()*M_conc_thin[i];
        }

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
FiniteElement::binaryIce()
{
    std::fstream input;
    std::string filename;
    int length;

    // First concentration
    filename = Environment::nextsimDir().string() + "/data/initConc.dat";
    input.open(filename, std::fstream::in);

    // get length of file:
    input.seekg (0, input.end);
    length = input.tellg();
    input.seekg (0, input.beg);
    if ( length != M_num_elements*sizeof(double) )
        throw std::runtime_error("Couldn't read the correct number of elements from " + filename +
                ". File length is " + std::to_string(length) + " while M_num_elements is " + std::to_string(M_num_elements) + ".\n");

    input.read((char*) &M_conc[0], M_num_elements*sizeof(double));
    input.close();

    // Then thickness
    filename = Environment::nextsimDir().string() + "/data/initThick.dat";
    input.open(filename, std::fstream::in);

    // get length of file:
    input.seekg (0, input.end);
    length = input.tellg();
    input.seekg (0, input.beg);
    if ( length != M_num_elements*sizeof(double) )
        throw std::runtime_error("Couldn't read the correct number of elements from " + filename +
                ". File length is " + std::to_string(length) + " while M_num_elements is " + std::to_string(M_num_elements) + ".\n");

    input.read((char*) &M_thick[0], M_num_elements*sizeof(double));
    input.close();

    // Finally snow thickness
    filename = Environment::nextsimDir().string() + "/data/initSnow.dat";
    input.open(filename, std::fstream::in);

    // get length of file:
    input.seekg (0, input.end);
    length = input.tellg();
    input.seekg (0, input.beg);
    if ( length != M_num_elements*sizeof(double) )
        throw std::runtime_error("Couldn't read the correct number of elements from " + filename +
                ". File length is " + std::to_string(length) + " while M_num_elements is " + std::to_string(M_num_elements) + ".\n");

    input.read((char*) &M_snow_thick[0], M_num_elements*sizeof(double));
    input.close();
}

void
FiniteElement::topazIce()
{
    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,"init - TOPAZ ice");
    M_external_data_tmp.resize(0);

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
FiniteElement::topazIceOsisafIcesat()
{
    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);
    
    external_data M_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);
    //M_type.check_and_reload(M_mesh,time_init);
 
    external_data M_icesat_thick=ExternalData(&M_ice_icesat_elements_dataset,M_mesh,0,false,time_init);
    //M_icesat_thick.check_and_reload(M_mesh,time_init);
    
    external_data M_conc_amsre=ExternalData(&M_ice_amsre_elements_dataset,M_mesh,0,false,time_init);
    //M_conc_amsre.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    M_external_data_tmp.push_back(&M_type);
    M_external_data_tmp.push_back(&M_icesat_thick);
    M_external_data_tmp.push_back(&M_conc_amsre);
    this->checkReloadDatasets(M_external_data_tmp,time_init,"init - TOPAZ ice");
    M_external_data_tmp.resize(0);
    
    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		tmp_var=std::min(1.,M_init_conc[i]);
		M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
		tmp_var=M_init_thick[i];
		M_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
		tmp_var=M_icesat_thick[i]*M_conc[i]; // Icesat gives the actual thickness (see "Uncertainties in Arctic sea ice thickness and volume: new estimates and implications for trends")
		M_thick[i] = (M_icesat_thick[i]>physical::hmin) ? tmp_var : M_thick[i]; 
		tmp_var=M_init_snow_thick[i];
		M_snow_thick[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.

        // Correction of the value given by Warren as a function of the ice type
        //M_type[i]==1. // No ice
        //M_type[i]==2. // First-Year ice
        //M_type[i]==3. // Multi-Year ice
        //M_type[i]==4. // Mixed
        double ratio_FYI=0.3;
        double ratio_MYI=0.9;
        double ratio_Mixed=0.6;

        if((M_thick[i]>0)&&(M_conc[i])>0.2)
        {
            if(M_type[i]>1. && M_type[i]<=2)
                M_ridge_ratio[i]=(M_type[i]-1.)*ratio_FYI;
            if(M_type[i]>2. && M_type[i]<=3.)
                M_ridge_ratio[i]=(1.-(M_type[i]-2.))*ratio_FYI + (M_type[i]-2.)*ratio_MYI;
            if(M_type[i]>3. && M_type[i]<=4.)
                M_ridge_ratio[i]=(1.-(M_type[i]-3.))*ratio_MYI + (M_type[i]-3.)*ratio_Mixed;
        }
        M_ridge_ratio[i]=M_ridge_ratio[i]*std::exp(ridging_exponent*(1.-M_conc[i]));

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

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            M_conc_thin[i]=std::max(M_conc_amsre[i]-M_conc[i],0.);
            M_h_thin[i]=M_conc_thin[i]*(h_thin_min+0.5*(h_thin_max-h_thin_min));
        }

		M_damage[i]=0.;
	}
}
void
FiniteElement::topazForecastIce()
{
    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,"init - TOPAZ ice forecast");
    M_external_data_tmp.resize(0);

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
    //M_conc_amsr2.check_and_reload(M_mesh,time_init);

    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_conc_amsr2);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init - TOPAZ ice forecast + AMSR2");
    M_external_data_tmp.resize(0);

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
    //M_conc_osisaf.check_and_reload(M_mesh,time_init);

    external_data M_confidence_osisaf=ExternalData(&M_ice_osisaf_elements_dataset,M_mesh,1,false,time_init);
    //M_confidence_osisaf.check_and_reload(M_mesh,time_init);

    external_data M_conc_amsr2=ExternalData(&M_ice_amsr2_elements_dataset,M_mesh,0,false,time_init);
    //M_conc_amsr2.check_and_reload(M_mesh,time_init);

    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,4,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_conc_osisaf);
    M_external_data_tmp.push_back(&M_confidence_osisaf);
    M_external_data_tmp.push_back(&M_conc_amsr2);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init - TOPAZ ice forecast + AMSR2 + OSISAF");
    M_external_data_tmp.resize(0);

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
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,1,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_piomas_elements_dataset,M_mesh,2,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init ice - PIOMAS");
    M_external_data_tmp.resize(0);

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
    //M_conc_amsre.check_and_reload(M_mesh,time_init);

    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_conc_amsre);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init ice - TOPAZ + AMSR-E");
    M_external_data_tmp.resize(0);

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
    //M_conc_amsr2.check_and_reload(M_mesh,time_init);

    external_data M_init_conc=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,0,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,1,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_init_snow_thick=ExternalData(&M_ice_topaz_elements_dataset,M_mesh,2,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_conc_amsr2);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init ice - TOPAZ + AMSR2");
    M_external_data_tmp.resize(0);

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
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_cs2_smos_elements_dataset,M_mesh,1,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    external_data M_type=ExternalData(&M_ice_osisaf_type_elements_dataset,M_mesh,0,false,time_init);
    //M_type.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_type);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init ice - CS2 + SMOS");
    M_external_data_tmp.resize(0);

    warrenClimatology();

    double tmp_var, correction_factor_warren;
    for (int i=0; i<M_num_elements; ++i)
    {
		tmp_var=std::min(1.,M_init_conc[i]);
		M_conc[i] = tmp_var;
		tmp_var=M_init_thick[i];
		M_thick[i] = tmp_var ;

        double ratio_FYI=0.3;
        double ratio_MYI=0.9;
        double ratio_Mixed=0.6;

        if((M_thick[i]>0)&&(M_conc[i])>0.2)
        {
            if(M_type[i]>1. && M_type[i]<=2)
                M_ridge_ratio[i]=(M_type[i]-1.)*ratio_FYI;
            if(M_type[i]>2. && M_type[i]<=3.)
                M_ridge_ratio[i]=(1.-(M_type[i]-2.))*ratio_FYI + (M_type[i]-2.)*ratio_MYI;
            if(M_type[i]>3. && M_type[i]<=4.)
                M_ridge_ratio[i]=(1.-(M_type[i]-3.))*ratio_MYI + (M_type[i]-3.)*ratio_Mixed;
        }
        M_ridge_ratio[i]=M_ridge_ratio[i]*std::exp(ridging_exponent*(1.-M_conc[i]));

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

        // Correction of the value given by Warren as a function of the ice type
        //M_type[i]==1. // No ice
        //M_type[i]==2. // First-Year ice
        //M_type[i]==3. // Multi-Year ice
        //M_type[i]==4. // Mixed
        correction_factor_warren=std::max(0.,std::min(1.,(M_type[i]-1.)*0.5)); // == 1. for MY, and mixed, 0.5 for FY, 0. for No ice

        M_snow_thick[i]=correction_factor_warren*M_snow_thick[i]*M_conc[i];

		M_damage[i]=0.;

        if(M_ice_cat_type==setup::IceCategoryType::THIN_ICE)
        {
            //M_conc_thin[i]=std::max(M_conc_amsre[i]-M_conc[i],0.);
            M_conc_thin[i]=std::min(1.-M_conc[i], 0.2*M_conc[i]);
            M_h_thin[i]=M_conc_thin[i]*(h_thin_min+0.5*(h_thin_max-h_thin_min));
        }

        // Check that the snow is not so thick that the ice is flooded
        double max_snow = M_thick[i]*(physical::rhow-physical::rhoi)/physical::rhos;
        M_snow_thick[i] = std::min(max_snow, M_snow_thick[i]);

	}
}
void
FiniteElement::smosIce()
{
    external_data M_init_conc=ExternalData(&M_ocean_elements_dataset,M_mesh,3,false,time_init);
    //M_init_conc.check_and_reload(M_mesh,time_init);

    external_data M_init_thick=ExternalData(&M_ice_smos_elements_dataset,M_mesh,0,false,time_init);
    //M_init_thick.check_and_reload(M_mesh,time_init);

    boost::gregorian::date dt = Nextsim::parse_date(time_init);
    int month_id=dt.month().as_number(); // 1 for January, 2 for February, and so on. This will be used to compute the snow from Warren climatology

    std::cout << "month_id: " << month_id <<"\n";

    external_data M_init_snow_thick=ExternalData(&M_ocean_elements_dataset,M_mesh,5,false,time_init);
    //M_init_snow_thick.check_and_reload(M_mesh,time_init);

    M_external_data_tmp.resize(0);
    M_external_data_tmp.push_back(&M_init_conc);
    M_external_data_tmp.push_back(&M_init_thick);
    M_external_data_tmp.push_back(&M_init_snow_thick);
    this->checkReloadDatasets(M_external_data_tmp,time_init,
            "init ice - SMOS");
    M_external_data_tmp.resize(0);

    double tmp_var;
    for (int i=0; i<M_num_elements; ++i)
    {
		tmp_var=std::min(1.,M_init_conc[i]);
		M_conc[i] = (tmp_var>1e-14) ? tmp_var : 0.; // TOPAZ puts very small values instead of 0.
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
FiniteElement::warrenClimatology()
{
    // Coefficients for the fit from Warren et al '99
    std::vector<std::vector<double>> coeffs(12);

    coeffs[0].push_back(28.01);
    coeffs[0].push_back( 0.1270);
    coeffs[0].push_back(-1.1833);
    coeffs[0].push_back(-0.1164);
    coeffs[0].push_back(-0.0051);
    coeffs[0].push_back( 0.0243);

    coeffs[1].push_back(30.28);
    coeffs[1].push_back( 0.1056);
    coeffs[1].push_back(-0.5908);
    coeffs[1].push_back(-0.0263);
    coeffs[1].push_back(-0.0049);
    coeffs[1].push_back( 0.0044);

    coeffs[2].push_back(33.89);
    coeffs[2].push_back( 0.5486);
    coeffs[2].push_back(-0.1996);
    coeffs[2].push_back( 0.0280);
    coeffs[2].push_back( 0.0216);
    coeffs[2].push_back(-0.0176);

    coeffs[3].push_back(36.80);
    coeffs[3].push_back( 0.4046);
    coeffs[3].push_back(-0.4005);
    coeffs[3].push_back( 0.0256);
    coeffs[3].push_back( 0.0024);
    coeffs[3].push_back(-0.0641);

    coeffs[4].push_back(36.93);
    coeffs[4].push_back( 0.0214);
    coeffs[4].push_back(-1.1795);
    coeffs[4].push_back(-0.1076);
    coeffs[4].push_back(-0.0244);
    coeffs[4].push_back(-0.0142);

    coeffs[5].push_back(36.59);
    coeffs[5].push_back( 0.7021);
    coeffs[5].push_back(-1.4819);
    coeffs[5].push_back(-0.1195);
    coeffs[5].push_back(-0.0009);
    coeffs[5].push_back(-0.0603);

    coeffs[6].push_back(11.02);
    coeffs[6].push_back( 0.3008);
    coeffs[6].push_back(-1.2591);
    coeffs[6].push_back(-0.0811);
    coeffs[6].push_back(-0.0043);
    coeffs[6].push_back(-0.0959);

    coeffs[7].push_back( 4.64);
    coeffs[7].push_back( 0.3100);
    coeffs[7].push_back(-0.6350);
    coeffs[7].push_back(-0.0655);
    coeffs[7].push_back( 0.0059);
    coeffs[7].push_back(-0.0005);

    coeffs[8].push_back(15.81);
    coeffs[8].push_back( 0.2119);
    coeffs[8].push_back(-1.0292);
    coeffs[8].push_back(-0.0868);
    coeffs[8].push_back(-0.0177);
    coeffs[8].push_back(-0.0723);

    coeffs[9].push_back(22.66);
    coeffs[9].push_back( 0.3594);
    coeffs[9].push_back(-1.3483);
    coeffs[9].push_back(-0.1063);
    coeffs[9].push_back( 0.0051);
    coeffs[9].push_back(-0.0577);

    coeffs[10].push_back(25.57);
    coeffs[10].push_back( 0.1496);
    coeffs[10].push_back(-1.4643);
    coeffs[10].push_back(-0.1409);
    coeffs[10].push_back(-0.0079);
    coeffs[10].push_back(-0.0258);

    coeffs[11].push_back(26.67);
    coeffs[11].push_back(-0.1876);
    coeffs[11].push_back(-1.4229);
    coeffs[11].push_back(-0.1413);
    coeffs[11].push_back(-0.0316);
    coeffs[11].push_back(-0.0029);

    // Time interpolation
    boost::gregorian::date t = Nextsim::parse_date(time_init);
    int month  = t.month().as_number(); // 1 for January, 2 for February, and so on.
    int day    = t.day().as_number(); // Day of the month.
    int eomday = t.end_of_month().day().as_number(); // Last day of the month.

    int month2; // The other month to interpolate from (one before or after)
    double dt;
    if ( day < eomday/2. )
    { // We're in the early part of the month and interpolate to the previous month
        month2 = month-1;
        dt     = eomday/2+day;
        if(month2==0)
            month2=12;
    }
    else
    { // We're in the late part of the month and interpolate to the next month
        month2 = month+1;
        dt     = eomday/2 + eomday - day;
        if(month2==13)
            month2=1;
    }

    // Now calculate snow thickness for the current day as an inexact temporal interpolation
    // I just assume all months are as long as the current one ... but the error incurred is small
    std::vector<double> lon = M_mesh.meanLon();
    std::vector<double> lat = M_mesh.meanLat();

    for ( int i=0; i<=M_num_elements; i++ )
    {
        const double pi = std::acos(-1);

        double x = (90 - lat[i]) * std::cos(lon[i]*pi/180.);
        double y = (90 - lat[i]) * std::sin(lon[i]*pi/180.);

        M_snow_thick[i] = 1e-2*dt/eomday*std::max(0.,
                coeffs[month-1][0] + coeffs[month-1][1]*x + coeffs[month-1][2]*y + coeffs[month-1][3]*x*y + coeffs[month-1][4]*x*x + coeffs[month-1][5]*y*y)
            + 1e-2*(eomday-dt)/eomday*std::max(0.,
                coeffs[month2-1][0] + coeffs[month2-1][1]*x + coeffs[month2-1][2]*y + coeffs[month2-1][3]*x*y + coeffs[month2-1][4]*x*x + coeffs[month2-1][5]*y*y);
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
    if(M_use_equallyspaced_drifters)
        this->equallySpacedDrifter();

    if(M_use_iabp_drifters)
        this->initIABPDrifter();

    if(M_use_rgps_drifters)
        this->initRGPSDrifters();

    if(M_use_osisaf_drifters)
        this->initOSISAFDrifters();
}

void
FiniteElement::coriolis()
{
    // Interpolation of the latitude
    std::vector<double> lat = M_mesh.meanLat();

    for (int i=0; i<M_fcor.size(); ++i)
    {
        if (vm["simul.use_coriolis"].as<bool>())
        {
            M_fcor[i] = 2*(physical::omega)*std::sin(lat[i]*PI/180.);
        }
        else
        {
            M_fcor[i] = 0.;
        }
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
    M_nfloes.assign(M_num_elements,0.);
    M_dfloe.assign(M_num_elements,0.);

    for (int i=0; i<M_num_elements; ++i)
    {
        if (M_conc[i]>=vm["wim.cicemin"].as<double>())
        {
            M_nfloes[i] = M_conc[i]/std::pow(vm["wim.dfloepackinit"].as<double>(),2.);
            M_dfloe[i]  = wim.nfloesToDfloe(M_nfloes[i],M_conc[i]);
        }
    }
    std::cout<<"initNfloes:\n";
    std::cout<<"init dfloe "<<vm["wim.dfloepackinit"].as<double>()<<"\n";
    std::cout<<"Min Nfloes = "<<*std::min_element(M_nfloes.begin(),M_nfloes.end())<<"\n";
    std::cout<<"Max Nfloes = "<<*std::max_element(M_nfloes.begin(),M_nfloes.end())<<"\n";
}//initNFloes
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
    std::vector<double> drifter_X(M_iabp_drifters.size());
    std::vector<double> drifter_Y(M_iabp_drifters.size());
    int j=0;
    for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); ++it )
    {
        drifter_X[j] = it->second[0];
        drifter_Y[j] = it->second[1];
        ++j;
    }

    // Interpolate the total displacement and concentration onto the drifter positions
    int nb_var=2;
    std::vector<double> interp_drifter_in(nb_var*M_mesh.numNodes());

    // Interpolate the velocity
    for (int i=0; i<M_mesh.numNodes(); ++i)
    {
        interp_drifter_in[nb_var*i]   = M_UT[i];
        interp_drifter_in[nb_var*i+1] = M_UT[i+M_mesh.numNodes()];
    }

    double* interp_drifter_out;
    InterpFromMeshToMesh2dx(&interp_drifter_out,
        &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
        M_mesh.numNodes(),M_mesh.numTriangles(),
        &interp_drifter_in[0],
        M_mesh.numNodes(),nb_var,
        &drifter_X[0],&drifter_Y[0],M_iabp_drifters.size(),
        true, 0.);

    // Loop over the map and output
    j=0;
    boost::gregorian::date           date = Nextsim::parse_date( current_time );
    boost::posix_time::time_duration time = Nextsim::parse_time( current_time );
    for ( auto it = M_iabp_drifters.begin(); it != M_iabp_drifters.end(); ++it )
    {
        double lat, lon;
        inverse_mapx(map,it->second[0]+interp_drifter_out[nb_var*j],it->second[1]+interp_drifter_out[nb_var*j+1],&lat,&lon);
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
        if ( M_iabp_drifters.count(number) == 0 )
        {
            double x, y;
            forward_mapx(map,lat,lon,&x,&y);
            M_iabp_drifters.emplace(number, std::array<double,2>{x, y});

        }
    }
    close_mapx(map);

    // Go through the M_iabp_drifters map and throw out the ones which IABP doesn't
    // report as being in the ice anymore
    for ( auto model = M_iabp_drifters.begin(); model != M_iabp_drifters.end(); /* ++model is not allowed here, because we use 'erase' */ )
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
            model = M_iabp_drifters.erase(model);
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
    M_equallyspaced_drifters = Drifters(1e3*vm["simul.drifter_spacing"].as<double>(), M_mesh, M_conc, vm["simul.drifter_climit"].as<double>());
    M_equallyspaced_drifters.initNetCDF(M_export_path+"/Drifters_", current_time);
    M_equallyspaced_drifters.appendNetCDF(current_time, M_mesh, M_UT);
}

void
FiniteElement::initRGPSDrifters()
{
    M_rgps_drifters = Drifters();
}

void
FiniteElement::updateRGPSDrifters()
{    
    std::string time_str = vm["simul.RGPS_time_init"].as<std::string>();
    double RGPS_time_init = from_date_time_string(time_str);
    
    std::string filename = Environment::nextsimDir().string() + "/data/RGPS_" + time_str + ".txt";
    M_rgps_drifters = Drifters(filename, M_mesh, M_conc, vm["simul.drifter_climit"].as<double>(),RGPS_time_init);
    
    M_rgps_drifters.initNetCDF(M_export_path+"/RGPS_Drifters_", current_time);
    M_rgps_drifters.appendNetCDF(current_time, M_mesh, M_UT);
}

void
FiniteElement::initOSISAFDrifters()
{
    M_osisaf_drifters.resize(2);
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
FiniteElement::exportInitMesh()
{
    Exporter exporter("float");
    std::string fileout;

    fileout = (boost::format( "%1%/mesh_init.bin" )
        % M_export_path ).str();

    LOG(INFO) <<"MESH BINARY: Exporter Filename= "<< fileout <<"\n";

    std::fstream meshbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
    if ( ! meshbin.good() )
        throw std::runtime_error("Cannot write to file: " + fileout);

    exporter.writeMesh(meshbin, M_mesh_init);
    meshbin.close();

    fileout = (boost::format( "%1%/mesh_init.dat" )
        % M_export_path ).str();

    LOG(INFO) <<"RECORD MESH: Exporter Filename= "<< fileout <<"\n";

    std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
    if ( ! outrecord.good() )
        throw std::runtime_error("Cannot write to file: " + fileout);
    exporter.writeRecord(outrecord,"mesh");
    outrecord.close();


    fileout = (boost::format( "%1%/field_init.bin" )
        % M_export_path ).str();

    LOG(INFO) <<"BINARY: Exporter Filename= "<< fileout <<"\n";

    std::fstream outbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
    if ( ! outbin.good() )
        throw std::runtime_error("Cannot write to file: " + fileout);
    std::vector<double> timevec(1);
    timevec[0] = current_time;
    exporter.writeField(outbin, timevec, "Time");

    exporter.writeField(outbin, M_hminVertices, "hminVertices");
    exporter.writeField(outbin, M_hmaxVertices, "hmaxVertices");
    exporter.writeField(outbin, M_dirichlet_flags, "M_dirichlet_flags");

    outbin.close();

    fileout = (boost::format( "%1%/field_init.dat" )
        % M_export_path ).str();

    LOG(INFO) <<"RECORD FIELD: Exporter Filename= "<< fileout <<"\n";

    std::fstream outrecord2(fileout, std::ios::out | std::ios::trunc);
    if ( ! outrecord2.good() )
        throw std::runtime_error("Cannot write to file: " + fileout);
    exporter.writeRecord(outrecord2);
    outrecord2.close();
}

void
FiniteElement::exportResults(int step, bool export_mesh, bool export_fields, bool apply_displacement)
{
    //define filenames from step
    std::string meshfile    = (boost::format( "%1%/mesh_%2%" )
                               % M_export_path
                               % step ).str();

    std::string fieldfile   = (boost::format( "%1%/field_%2%" )
                               % M_export_path
                               % step ).str();

    std::vector<std::string> filenames = {meshfile,fieldfile}; 
    this->exportResults(filenames, export_mesh, export_fields, apply_displacement);
}

void
FiniteElement::exportResults(std::vector<std::string> const &filenames, bool export_mesh, bool export_fields, bool apply_displacement)
{
    Exporter exporter(vm["setup.exporter_precision"].as<std::string>());
    std::string fileout;


    if (export_mesh)
    {
        fileout = filenames[0]+".bin";
        LOG(INFO) <<"MESH BINARY: Exporter Filename= "<< fileout <<"\n";

        if(apply_displacement)
            // move the mesh for the export
            M_mesh.move(M_UM,1.);

        std::fstream meshbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
        if ( ! meshbin.good() )
            throw std::runtime_error("Cannot write to file: " + fileout);

        exporter.writeMesh(meshbin, M_mesh);
        meshbin.close();

        if(apply_displacement)
            // move it back after the export
		    M_mesh.move(M_UM,-1.);

        fileout = filenames[0]+".dat";

        LOG(INFO) <<"RECORD MESH: Exporter Filename= "<< fileout <<"\n";

        std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
        if ( ! outrecord.good() )
            throw std::runtime_error("Cannot write to file: " + fileout);
        exporter.writeRecord(outrecord,"mesh");
        outrecord.close();
    }


    if (export_fields)
    {
        fileout = filenames[1]+".bin";
        LOG(INFO) <<"BINARY: Exporter Filename= "<< fileout <<"\n";

        std::fstream outbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
        if ( ! outbin.good() )
            throw std::runtime_error("Cannot write to file: " + fileout);
        std::vector<double> timevec(1);
        timevec[0] = current_time;
        std::vector<int> regridvec(1);
        regridvec[0] = M_nb_regrid;
        exporter.writeField(outbin, timevec, "Time");
        exporter.writeField(outbin, regridvec, "M_nb_regrid");
        exporter.writeField(outbin, M_surface, "Element_area");
        exporter.writeField(outbin, M_VT, "M_VT");
        exporter.writeField(outbin, M_dirichlet_flags, "M_dirichlet_flags");
        exporter.writeField(outbin, M_conc, "Concentration");
        exporter.writeField(outbin, M_thick, "Thickness");
        exporter.writeField(outbin, M_snow_thick, "Snow");
        exporter.writeField(outbin, M_damage, "Damage");
        exporter.writeField(outbin, M_ridge_ratio, "Ridge_ratio");

        std::vector<double> AllMinAngle = this->AllMinAngle(M_mesh, M_UM, 0.);
        exporter.writeField(outbin, AllMinAngle, "AllMinAngle");

        int i=0;
        for (auto it=M_tice.begin(); it!=M_tice.end(); it++)
        {
            exporter.writeField(outbin, *it, "Tice_"+std::to_string(i));
            i++;
        }
        exporter.writeField(outbin, M_sst, "SST");
        exporter.writeField(outbin, M_sss, "SSS");

        std::vector<double> PreviousNumbering(M_mesh.numNodes());

        for ( int i=0; i<M_mesh.numNodes(); ++i )
            PreviousNumbering[i]=bamgmesh->PreviousNumbering[i];

        exporter.writeField(outbin, PreviousNumbering, "PreviousNumbering");

        if(vm["simul.save_forcing_field"].as<bool>())
        {
            // Thermodynamic and dynamic forcing
            // Atmosphere
            std::vector<std::string> ext_data_names;
#if 1
            M_external_data_tmp.push_back(&M_wind);         // Surface wind [m/s]
            ext_data_names.push_back("M_wind");
            M_external_data_tmp.push_back(&M_tair);         // 2 m temperature [C]
            ext_data_names.push_back("M_tair");
            M_external_data_tmp.push_back(&M_mixrat);       // Mixing ratio
            ext_data_names.push_back("M_mixrat");
            M_external_data_tmp.push_back(&M_mslp);         // Atmospheric pressure [Pa]
            ext_data_names.push_back("M_mslp");
            M_external_data_tmp.push_back(&M_Qsw_in);       // Incoming short-wave radiation [W/m2]
            ext_data_names.push_back("M_Qsw_in");
            M_external_data_tmp.push_back(&M_Qlw_in);       // Incoming long-wave radiation [W/m2]
            ext_data_names.push_back("M_Qlw_in");
            M_external_data_tmp.push_back(&M_tcc);          // Total cloud cover [?]
            ext_data_names.push_back("M_tcc");
            M_external_data_tmp.push_back(&M_precip);       // Total precipitation [m]
            ext_data_names.push_back("M_precip");
            M_external_data_tmp.push_back(&M_snowfr);       // Fraction of precipitation that is snow
            ext_data_names.push_back("M_snowfr");
            M_external_data_tmp.push_back(&M_snowfall);       // Fraction of precipitation that is snow
            ext_data_names.push_back("M_snowfall");
            M_external_data_tmp.push_back(&M_dair);         // 2 m dew point [C]
            ext_data_names.push_back("M_dair");
#else
            if ( M_wind.M_initialized )
                exporter.writeField(outbin,M_wind.getVector(), "M_wind");         // Surface wind [m/s]
            if ( M_tair.M_initialized )
                exporter.writeField(outbin,M_tair.getVector(), "M_tair");         // 2 m temperature [C]
            if ( M_mixrat.M_initialized )
                exporter.writeField(outbin,M_mixrat.getVector(), "M_mixrat");       // Mixing ratio
            if ( M_mslp.M_initialized )
                exporter.writeField(outbin,M_mslp.getVector(), "M_mslp");         // Atmospheric pressure [Pa]
            if ( M_Qsw_in.M_initialized )
                exporter.writeField(outbin,M_Qsw_in.getVector(), "M_Qsw_in");       // Incoming short-wave radiation [W/m2]
            if ( M_Qlw_in.M_initialized )
                exporter.writeField(outbin,M_Qlw_in.getVector(), "M_Qlw_in");       // Incoming long-wave radiation [W/m2]
            if ( M_tcc.M_initialized )
                exporter.writeField(outbin,M_tcc.getVector(), "M_tcc");       // Incoming long-wave radiation [W/m2]
            if ( M_precip.M_initialized )
                exporter.writeField(outbin,M_precip.getVector(), "M_precip");       // Total precipitation rate [kg/m^2/s]
            if ( M_snowfr.M_initialized )
                exporter.writeField(outbin,M_snowfr.getVector(), "M_snowfr");       // Fraction of precipitation that is snow
            if ( M_snowfall.M_initialized )
                exporter.writeField(outbin,M_snowfall.getVector(), "M_snowfall");       // Snow fall rate [kg/m^2/s]
            if ( M_dair.M_initialized )
                exporter.writeField(outbin,M_dair.getVector(), "M_dair");         // 2 m dew point [C]
#endif

            // Ocean
#if 1
            M_external_data_tmp.push_back(&M_ocean);        // "Geostrophic" ocean currents [m/s]
            ext_data_names.push_back("M_ocean");
            M_external_data_tmp.push_back(&M_ssh);          // Sea surface elevation [m]
            ext_data_names.push_back("M_ssh");
            M_external_data_tmp.push_back(&M_ocean_temp);   // Ocean temperature in top layer [C]
            ext_data_names.push_back("M_ocean_temp");
            M_external_data_tmp.push_back(&M_ocean_salt);   // Ocean salinity in top layer [?]
            ext_data_names.push_back("M_ocean_salt");
            M_external_data_tmp.push_back(&M_mld);           // Mixed-layer depth [m]
            ext_data_names.push_back("M_mld");
#else
            if ( M_ocean.M_initialized )
                exporter.writeField(outbin,M_ocean.getVector(), "M_ocean");        // "Geostrophic" ocean currents [m/s]
            if ( M_ssh.M_initialized )
                exporter.writeField(outbin,M_ssh.getVector(), "M_ssh");          // Sea surface elevation [m]

            if ( M_ocean_temp.M_initialized )
                exporter.writeField(outbin,M_ocean_temp.getVector(), "M_ocean_temp");   // Ocean temperature in top layer [C]
            if ( M_ocean_salt.M_initialized )
                exporter.writeField(outbin,M_ocean_salt.getVector(), "M_ocean_salt");   // Ocean salinity in top layer [C]
            if ( M_mld.M_initialized )
                exporter.writeField(outbin,M_mld.getVector(), "M_mld");          // Mixed-layer depth [m]
#endif

            // Bathymetry
#if 1
            M_external_data_tmp.push_back(&M_element_depth);           // Mixed-layer depth [m]
            ext_data_names.push_back("M_element_depth");
#else
            if ( M_element_depth.M_initialized )
                exporter.writeField(outbin,M_element_depth.getVector(), "M_element_depth");
#endif

#if 1
            //loop over external data pointers and check if they should be saved
            for (int i=0;i<M_external_data_tmp.size();i++)
            {
                if ((M_external_data_tmp[i]->M_initialized)&&
                    (!M_external_data_tmp[i]->M_is_constant))
                {
                    exporter.writeField(outbin,M_external_data_tmp[i]->getVector(), ext_data_names[i]);
                }
            }

            //clear
            M_external_data_tmp.resize(0);
            ext_data_names.resize(0);
#endif
        }//save forcing

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
            exporter.writeField(outbin, M_conc_thin, "Concentration_thin_ice");
        }

#if 1
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
#endif
        outbin.close();

        fileout = filenames[1]+".dat";
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

        int nb_var=3;
        double* interp_elt_out;//output
        std::vector<double> interp_elt_in(nb_var*M_num_elements);

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

            // effective thickness (volume)
            interp_elt_in[nb_var*i+tmp_nb_var] = M_thick[i];
            tmp_nb_var++;

            // Nfloes
            interp_elt_in[nb_var*i+tmp_nb_var] = M_nfloes[i];
            tmp_nb_var++;

            if(tmp_nb_var>nb_var)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }

        // move the mesh for the interpolation on to the wim grid
		M_mesh.move(M_UM,1.);

        //needed for interp (mesh2mesh) and assigning outputs
        //TODO: add bool "regular" to wim_grid

        if (1)//(wim_grid.regular)
        {
            std::cout<<"sim2wim: before interp mesh2grid\n";
            InterpFromMeshToGridx(interp_elt_out,
                                  &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                                  M_mesh.numNodes(),M_mesh.numTriangles(),
                                  &interp_elt_in[0],
                                  M_mesh.numTriangles(),nb_var,
                                  xmin_wim,ymax_wim,
                                  wim_grid.dx,wim_grid.dy,
                                  wim_grid.nx,wim_grid.ny,
                                  0.);
        }
        else
        {
            //std::cout<<"sim2wim: before interp mesh2mesh2d\n";
            InterpFromMeshToMesh2dx(&interp_elt_out,
                                  &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                                  M_mesh.numNodes(),M_mesh.numTriangles(),
                                  &interp_elt_in[0],
                                  M_mesh.numNodes(),nb_var,
                                  &wim_grid.X[0], &wim_grid.Y[0], num_elements_wim_grid,
                                  false);
        }

        if (1)
        {
            std::cout<<"nx = "<< wim_grid.nx <<"\n";
            std::cout<<"ny = "<< wim_grid.ny <<"\n";
            std::cout<<"dx = "<< wim_grid.dx <<"\n";
            std::cout<<"dy = "<< wim_grid.dy <<"\n";

            auto RX = M_mesh.coordX();
            auto RY = M_mesh.coordY();
            std::cout<<"MIN BOUND MESHX= "<< *std::min_element(RX.begin(),RX.end()) <<"\n";
            std::cout<<"MAX BOUND MESHX= "<< *std::max_element(RX.begin(),RX.end()) <<"\n";

            std::cout<<"MIN BOUND MESHY= "<< *std::min_element(RY.begin(),RY.end()) <<"\n";
            std::cout<<"MAX BOUND MESHY= "<< *std::max_element(RY.begin(),RY.end()) <<"\n";

            std::cout<<"------------------------------------------\n";

            std::cout<<"MIN BOUND GRIDX= "<< *std::min_element(wim_grid.X.begin(),wim_grid.X.end()) <<"\n";
            std::cout<<"MAX BOUND GRIDX= "<< *std::max_element(wim_grid.X.begin(),wim_grid.X.end()) <<"\n";

            std::cout<<"MIN BOUND GRIDY= "<< *std::min_element(wim_grid.Y.begin(),wim_grid.Y.end()) <<"\n";
            std::cout<<"MAX BOUND GRIDY= "<< *std::max_element(wim_grid.Y.begin(),wim_grid.Y.end()) <<"\n";
        }

        //std::cout<<"after interp mesh2grid\n";
        //std::cout<<"ideal wave forcing: "<<wim_ideal_forcing<<"\n";

        // move back the mesh after the interpolation
		M_mesh.move(M_UM,-1.);

        if (!step)
        {
            M_icec_grid.assign  (num_elements_wim_grid,0.);
            M_iceh_grid.assign  (num_elements_wim_grid,0.);
            M_nfloes_grid.assign(num_elements_wim_grid,0.);
            M_taux_grid.assign  (num_elements_wim_grid,0.);
            M_tauy_grid.assign  (num_elements_wim_grid,0.);
        }

        // set inputs to WIM:
        // - ice from interpolation from mesh
        // - waves from datasets if needed
        for (int i=0; i<num_elements_wim_grid; ++i)
        {
            tmp_nb_var=0;

            // concentration
            M_icec_grid[i] = interp_elt_out[nb_var*i+tmp_nb_var];
            tmp_nb_var++;

            // effective thickness (volume)
            M_iceh_grid[i] = interp_elt_out[nb_var*i+tmp_nb_var];
            tmp_nb_var++;

            // Nfloes
            M_nfloes_grid[i] = interp_elt_out[nb_var*i+tmp_nb_var];
            tmp_nb_var++;

            if(tmp_nb_var>nb_var)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }

            if ( !wim_ideal_forcing )
            {
                //get incident waves from datasets
                double cfac = 1.;
                double uwave,vwave;
                if ( wim_forcing_options.use_ice )
                {
                    //cancel waves if ice present
                    if (M_fice_waves[i]>0.)
                        cfac = 0.;
                }

                // significant wave height
                // - Hs given to the WIM should have the waves-in-ice
                // removed (so we can do our own attenuation)
                M_SWH_grid[i] = cfac*M_SWH[i];

                // mean wave direction
                uwave   = cfac*M_MWD[i];
                vwave   = cfac*M_MWD[i+num_elements_wim_grid];
                if ( std::hypot(uwave,vwave)>.5 )
                {
                    // if there are waves |(uwave,vwave)|=1,
                    // so convert to wave-from direction
                    // (degrees, clockwise from north)
                    M_MWD_grid[i] = 90.-(180./PI)*std::atan2(-uwave,-vwave);
                }

                // peak wave period
                if ( M_MWP[i]>0. )
                    if ( wim_forcing_options.use_mwp )
                        //Tp given to the WIM should have the waves-in-ice
                        //removed (so we can do our own attenuation)
                        M_MWP_grid[i] = cfac*M_MWP[i];
                    else 
                        // we are given fp, so convert to Tp,
                        // taking account of the ice
                        M_MWP_grid[i] = cfac/M_MWP[i];
                else
                {
                    //if fp or Tp are 0, set all wave inputs to 0
                    M_MWP_grid[i] = 0.;
                    M_SWH_grid[i] = 0.;
                    M_MWD_grid[i] = 0.;
                }
                    
            }//using non-ideal wave forcing
        }//loop over wim grid cells

#if 1
        //test interp
        std::cout<<"sim2wim: check ice inputs to WIM\n";
        std::cout<<"min conc   grid= "<< *std::min_element(M_icec_grid.begin(),M_icec_grid.end() )<<"\n";
        std::cout<<"max conc   grid= "<< *std::max_element(M_icec_grid.begin(),M_icec_grid.end() )<<"\n";
        std::cout<<"min thick  grid= "<< *std::min_element(M_iceh_grid.begin(),M_iceh_grid.end() )<<"\n";
        std::cout<<"max thick  grid= "<< *std::max_element(M_iceh_grid.begin(),M_iceh_grid.end() )<<"\n";
        std::cout<<"min Nfloes grid= "<< *std::min_element(M_nfloes_grid.begin(),M_nfloes_grid.end() )<<"\n";
        std::cout<<"max Nfloes grid= "<< *std::max_element(M_nfloes_grid.begin(),M_nfloes_grid.end() )<<"\n";
#endif

        LOG(DEBUG)<<"sim2wim (check wave forcing): "<<wim_ideal_forcing<<","<<M_SWH_grid.size()<<"\n";
        if (M_SWH_grid.size()>0)//( !wim_ideal_forcing )
        {
#if 0
            LOG(DEBUG)<<"min SWH_grid= "<< *std::min_element(M_SWH_grid.begin(),M_SWH_grid.end() )<<"\n";
            LOG(DEBUG)<<"max SWH_grid= "<< *std::max_element(M_SWH_grid.begin(),M_SWH_grid.end() )<<"\n";
            LOG(DEBUG)<<"min MWD_grid= "<< *std::min_element(M_MWD_grid.begin(),M_MWD_grid.end() )<<"\n";
            LOG(DEBUG)<<"max MWD_grid= "<< *std::max_element(M_MWD_grid.begin(),M_MWD_grid.end() )<<"\n";
            LOG(DEBUG)<<"min MWP_grid= "<< *std::min_element(M_MWP_grid.begin(),M_MWP_grid.end() )<<"\n";
            LOG(DEBUG)<<"max MWP_grid= "<< *std::max_element(M_MWP_grid.begin(),M_MWP_grid.end() )<<"\n";
#elif 0
            std::cout<<"nextsimToWim check wave inputs at t="<<current_time<<"\n";
            std::cout<<"min SWH_grid= "<< *std::min_element(M_SWH_grid.begin(),M_SWH_grid.end() )<<"\n";
            std::cout<<"max SWH_grid= "<< *std::max_element(M_SWH_grid.begin(),M_SWH_grid.end() )<<"\n";
            std::cout<<"min MWD_grid= "<< *std::min_element(M_MWD_grid.begin(),M_MWD_grid.end() )<<"\n";
            std::cout<<"max MWD_grid= "<< *std::max_element(M_MWD_grid.begin(),M_MWD_grid.end() )<<"\n";
            std::cout<<"min MWP_grid= "<< *std::min_element(M_MWP_grid.begin(),M_MWP_grid.end() )<<"\n";
            std::cout<<"max MWP_grid= "<< *std::max_element(M_MWP_grid.begin(),M_MWP_grid.end() )<<"\n";
#elif 0
            int icheck      = 116;
            int jcheck      = 112;
            int Icheck      = jcheck+icheck*wim_grid.ny;
            double xcheck   = wim_grid.X[Icheck];
            double ycheck   = wim_grid.Y[Icheck];
            std::cout<<"nextsimToWim check wave inputs at t="<<current_time<<"\n";
            std::cout<<"near i,j,x,y = "<<icheck<<","<<jcheck<<","<<xcheck<<","<<ycheck<<","<<"\n";
            std::cout<<"SWH_grid= "<< M_SWH_grid[Icheck]<<"\n";
            std::cout<<"MWD_grid= "<< M_MWD_grid[Icheck]<<"\n";
            std::cout<<"MWP_grid= "<< M_MWP_grid[Icheck]<<"\n";
#endif
        }

        xDelete<double>(interp_elt_out);
    }//if (M_run_wim)
}//nextsimToWim
#endif

#if defined (WAVES)
void
FiniteElement::wimToNextsim(bool step)
{
    bool break_on_mesh =
        ( vm["nextwim.coupling-option"].template as<std::string>() == "breaking_on_mesh");
    //std::cout<<"break_on_mesh="<<break_on_mesh<<"\n";

    if (!M_regrid)
        M_mesh.move(M_UM,1.);

    if (M_run_wim)
    {

        // run wim
        if ( break_on_mesh )
        {
            //give moved mesh to WIM
            wim.setMesh(M_mesh.bcoordX(),//elements of mesh (x)
                        M_mesh.bcoordY(),//elements of mesh (y)
                        M_conc,
                        M_thick,
                        M_nfloes,
                        "m");
                        //"km");
        }

        bool TEST_INTERP_MESH = false;
        //save mesh before entering WIM:
        // mesh file can then be copied inside WIM to correct path to allow plotting
        if (TEST_INTERP_MESH)
            this->exportResults(1001,true,false);

        LOG(DEBUG)<<"wim2sim (check wave forcing): "<<wim_ideal_forcing<<","<<M_SWH_grid.size()<<"\n";
        if (M_SWH_grid.size()>0)//( !wim_ideal_forcing )
        {
            LOG(DEBUG)<<"min SWH_grid= "<< *std::min_element(M_SWH_grid.begin(),M_SWH_grid.end() )<<"\n";
            LOG(DEBUG)<<"max SWH_grid= "<< *std::max_element(M_SWH_grid.begin(),M_SWH_grid.end() )<<"\n";
            LOG(DEBUG)<<"min MWD_grid= "<< *std::min_element(M_MWD_grid.begin(),M_MWD_grid.end() )<<"\n";
            LOG(DEBUG)<<"max MWD_grid= "<< *std::max_element(M_MWD_grid.begin(),M_MWD_grid.end() )<<"\n";
            LOG(DEBUG)<<"min MWP_grid= "<< *std::min_element(M_MWP_grid.begin(),M_MWP_grid.end() )<<"\n";
            LOG(DEBUG)<<"max MWP_grid= "<< *std::max_element(M_MWP_grid.begin(),M_MWP_grid.end() )<<"\n";
        }

        wim.run(M_icec_grid, M_iceh_grid, M_nfloes_grid, M_SWH_grid, M_MWP_grid, M_MWD_grid, step);
        M_taux_grid = wim.getTaux();
        M_tauy_grid = wim.getTauy();

        if ( !break_on_mesh )
            M_nfloes_grid = wim.getNfloes();
        else
        {
            M_nfloes        = wim.getNfloesMesh();
            auto M_broken   = wim.getBrokenMesh();//TODO check this - maybe change damage later
            wim.clearMesh();

            if (vm["nextwim.wim_damage_mesh"].template as<bool>())
            {
                for (int i=1;i<M_num_elements;i++)
                {
                    if (M_broken[i])
                        M_damage[i] = max(M_damage[i],
                           (vm["nextwim.wim_damage_value"].template as<double>()));
                    //std::cout<<"broken?,damage"<<M_broken[i]<<","<<M_damage[i]<<"\n";
                }
            }
        }
    }//run WIM & get outputs on grid

    //if (!M_regrid)
    //    M_mesh.move(M_UM,1.);

    bool interp_taux = vm["nextwim.applywavestress"].as<bool>();
    // can turn off effect of wave stress for testing
    // - if this is not done, we currently interp tau_x,tau_y each time step
    // TODO rethink this? (let them be advected? - this could lead to instability perhaps)

    // Set type of interpolation for grid-to-mesh
    // int interptype = TriangleInterpEnum;
    int interptype = BilinearInterpEnum;
    //int interptype = NearestInterpEnum;

    int    nx = wim_grid.nx;
    int    ny = wim_grid.ny;
    double dx = wim_grid.dx;
    double dy = wim_grid.dy;

    if (interp_taux)
    {
        chrono.restart();
        LOG(DEBUG) <<"Nodal Interp starts\n";
        LOG(DEBUG) <<"NODAL: Interp starts\n";

        // NODAL INTERPOLATION
        // - taux and tauy from waves
        if ((!step) || M_regrid)
           //initialisation or after regridding (need to reset sizes)
            M_tau.assign(2*M_num_nodes,0);


        int nb_var=2;
        std::vector<double> interp_in(nb_var*num_elements_wim_grid,0.);
        double* interp_out;

        for (int i=0; i<num_elements_wim_grid; ++i)
        {
            // interp taux,tauy to nodes of mesh
            interp_in[nb_var*i]   = M_taux_grid[i];
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
    }//interp taux,tauy

    if (M_run_wim)
    {
        //interp Nfloes
        if (!break_on_mesh)
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

            // if conc is not too low, reset Nfloes after breaking
            // (& Dfloe will also change)
            M_nfloes.assign(M_num_elements,0.);
            M_dfloe.assign(M_num_elements,0.);
            for (int i=0; i<M_num_elements; ++i)
                if (M_conc[i] >= vm["wim.cicemin"].template as<double>())
                {
                    M_nfloes[i] = interp_out[i];
                    M_dfloe[i]  = wim.nfloesToDfloe(M_nfloes[i],M_conc[i]);
                }

            xDelete<double>(interp_out);
            std::cout<<"Min Nfloes on grid = "<< *std::min_element(M_nfloes_grid.begin(),M_nfloes_grid.end()) <<"\n";
            std::cout<<"Max Nfloes on grid = "<< *std::max_element(M_nfloes_grid.begin(),M_nfloes_grid.end()) <<"\n";
            std::cout<<"Min Nfloes on mesh = "<< *std::min_element(M_nfloes.begin(),M_nfloes.end()) <<"\n";
            std::cout<<"Max Nfloes on mesh = "<< *std::max_element(M_nfloes.begin(),M_nfloes.end()) <<"\n";
            //std::abort();
        }
    }

    if (!M_regrid)//move the mesh back
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
    const char* senv = ::getenv(envname.c_str());
    if ( senv == NULL )
        senv = "NULL";
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
