/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   finiteelement.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
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
FiniteElement::init()
{
    std::cout <<"GMSH VERSION= "<< M_mesh.version() <<"\n";
    M_mesh.setOrdering("bamg");

    M_mesh_filename = vm["simul.mesh_filename"].as<std::string>();

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

    switch (M_domain_type)
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

    M_mesh.stereographicProjection();
    // M_mesh.writeTofile("copy_init_mesh.msh");

    // createGMSHMesh("hypercube.geo");
    // //M_mesh.setOrdering("gmsh");
    // M_mesh.readFromFile("hypercube.msh");
    std::cout<<"START CONSTANTS\n";
    this->initConstant();
    this->initBamg();

    std::cout<<"Convert MESH starts\n";
    BamgConvertMeshx(
                     bamgmesh,bamggeom,
                     &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                     M_mesh.numNodes(), M_mesh.numTriangles()
                     );
    std::cout<<"Convert MESH done\n";


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

    std::cout<<"MESH: HMIN= "<< h[0] <<"\n";
    std::cout<<"MESH: HMAX= "<< h[1] <<"\n";
    std::cout<<"MESH: RES = "<< this->resolution(M_mesh) <<"\n";

    switch (M_mesh_type)
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

    M_solver = solver_ptrtype(new solver_type());
    M_matrix = matrix_ptrtype(new matrix_type());
    M_vector = vector_ptrtype(new vector_type());
    M_solution = vector_ptrtype(new vector_type());

    M_reuse_prec = true;

    M_ice_type = setup::IceCategoryType::CLASSIC;

    const boost::unordered_map<const std::string, setup::AtmosphereType> str2atmosphere = boost::assign::map_list_of
        ("constant", setup::AtmosphereType::CONSTANT)
        ("asr", setup::AtmosphereType::ASR);
    M_atmosphere_type = str2atmosphere.find(vm["setup.atmosphere-type"].as<std::string>())->second;

    //std::cout<<"AtmosphereType= "<< (int)M_atmosphere_type <<"\n";

    const boost::unordered_map<const std::string, setup::OceanType> str2ocean = boost::assign::map_list_of
        ("constant", setup::OceanType::CONSTANT)
        ("topaz", setup::OceanType::TOPAZR);
    M_ocean_type = str2ocean.find(vm["setup.ocean-type"].as<std::string>())->second;

    //std::cout<<"OCEANTYPE= "<< (int)M_ocean_type <<"\n";

    const boost::unordered_map<const std::string, setup::ConcentrationType> str2conc = boost::assign::map_list_of
        ("constant", setup::ConcentrationType::CONSTANT)
        ("topaz", setup::ConcentrationType::TOPAZ4);
    M_conc_type = str2conc.find(vm["setup.concentration-type"].as<std::string>())->second;

    //std::cout<<"CONCTYPE= "<< (int)M_conc_type <<"\n";

    const boost::unordered_map<const std::string, setup::ThicknessType> str2thick = boost::assign::map_list_of
        ("constant", setup::ThicknessType::CONSTANT)
        ("topaz", setup::ThicknessType::TOPAZ4);
    M_thick_type = str2thick.find(vm["setup.thickness-type"].as<std::string>())->second;

    //std::cout<<"THICKTYPE= "<< (int)M_thick_type <<"\n";

    const boost::unordered_map<const std::string, setup::SnowThicknessType> str2snow = boost::assign::map_list_of
        ("constant", setup::SnowThicknessType::CONSTANT)
        ("topaz", setup::SnowThicknessType::TOPAZ4);
    M_snow_thick_type = str2snow.find(vm["setup.snow-thickness-type"].as<std::string>())->second;

    //std::cout<<"SNOWTHICKTYPE= "<< (int)M_snow_thick_type <<"\n";

    const boost::unordered_map<const std::string, setup::DamageType> str2damg = boost::assign::map_list_of
        ("constant", setup::DamageType::CONSTANT);
    M_damage_type = str2damg.find(vm["setup.damage-type"].as<std::string>())->second;

    //std::cout<<"DAMAGETYPE= "<< (int)M_damage_type <<"\n";


    // init options for interpolation from mesh to mesh
    // options = new Options();
}

// Initialise all physical variables to propper initial conditions
void
FiniteElement::initSimulation()
{
    chrono_tot.restart();

    //M_matrix->init(2*M_num_nodes,2*M_num_nodes,22);
    M_matrix->init(2*M_num_nodes,2*M_num_nodes,M_graph);
    M_vector->init(2*M_num_nodes);
    M_solution->init(2*M_num_nodes);

    M_VT.resize(2*M_num_nodes,0.);
    M_VTM.resize(2*M_num_nodes,0.);
    M_VTMM.resize(2*M_num_nodes,0.);
    M_vector_reduction.resize(2*M_num_nodes,0.);
    M_valid_conc.resize(2*M_num_nodes,false);

    M_wind.resize(2*M_num_nodes);
    M_ocean.resize(2*M_num_nodes);

    M_tair.resize(M_num_elements);
    M_mixrat.resize(M_num_elements);
    M_dair.resize(M_num_elements);
    M_mslp.resize(M_num_elements);
    M_Qsw_in.resize(M_num_elements);
    M_Qlw_in.resize(M_num_elements);
    M_precip.resize(M_num_elements);
    M_snowfr.resize(M_num_elements);

    M_tair2.resize(2);
    M_mixrat2.resize(2);
    M_dair2.resize(2);
    M_mslp2.resize(2);
    M_Qsw_in2.resize(2);
    M_Qlw_in2.resize(2);
    M_precip2.resize(2);
    M_snowfr2.resize(2);

    M_ocean_salt2.resize(2);
    M_ocean_temp2.resize(2);
    M_mld2.resize(2);

    M_ocean_temp.resize(M_num_elements);
    M_ocean_salt.resize(M_num_elements);
    M_mld.resize(M_num_elements);

    M_sst.resize(M_num_elements);
    M_sss.resize(M_num_elements);

    M_bathy_depth.resize(M_mesh_init.numNodes(),200.);

    M_hminVertices.resize(M_mesh_init.numNodes(),1e-100);
    M_hmaxVertices.resize(M_mesh_init.numNodes(),1e100);

    M_uair2.resize(2);
    M_vair2.resize(2);
    M_uoce2.resize(2);
    M_voce2.resize(2);
    M_ssh2.resize(2);

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
    M_tsurf.resize(M_num_elements);
    M_sst.resize(M_num_elements);
    M_sss.resize(M_num_elements);

    for (int i=0; i<M_num_elements; ++i)
    {
        M_divergence_rate[i] = 0.;
        M_sigma[i]=0.;
        M_sigma[i+M_num_elements]=0.;
        M_sigma[i+2*M_num_elements]=0.;
        if ((M_conc[i] <= 0.) || (M_thick[i] <= 0.) )
        {
            M_conc[i] = 0.;
            M_thick[i] = 0.;
        }
    }

    M_ssh.resize(M_num_nodes,0.);

    M_surface.resize(M_num_elements);

    M_fcor.resize(M_num_elements);

    Dimension asr_dimension_x={
        name:"x",
        start:0,
        end:360};

    Dimension asr_dimension_y={
        name:"y",
        start:0,
        end:360};

    Dimension asr_dimension_time={
        name:"time", // "Time"
        start:0,
        end:248};

    std::vector<Dimension> dimensions_asr_latlon(2);
    dimensions_asr_latlon[0] = asr_dimension_y;
    dimensions_asr_latlon[1] = asr_dimension_x;

    std::vector<Dimension> dimensions_asr(3);
    dimensions_asr[0] = asr_dimension_time;
    dimensions_asr[1] = asr_dimension_y;
    dimensions_asr[2] = asr_dimension_x;

    std::vector<Dimension> dimensions_asr_time(1);
    dimensions_asr_time[0] = asr_dimension_time;

    netCDF::NcVar NcVar_tmp;

    Variable asr_latitude={
        name: "XLAT",
        dimensions: dimensions_asr_latlon,
        a: 1.,
        b: 0.,
        Units: "degree_north",
        NcVar: NcVar_tmp};

    Variable asr_longitude={
        name: "XLONG",
        dimensions: dimensions_asr_latlon,
        a: 1.,
        b: 0.,
        Units: "degree_east",
        NcVar: NcVar_tmp};

    Variable asr_time={
        name: "time",
        dimensions: dimensions_asr_time,
        a: 1.,
        b: 0.,
        Units: "hours",
        NcVar: NcVar_tmp};

    // rotate EB coordinates to fit the ASR coords
    double angle_stereo_mesh = -45;
    double angle_stereo_ASR = -175;
    double diff_angle = -(angle_stereo_mesh-angle_stereo_ASR)*PI/180.;

    M_asr_grid={
        interpolation_method: setup::InterpolationType::InterpFromGridToMesh,
        dirname:"data",
        prefix:"asr30km.comb.2d.", // "asr30km.comb.2D.";
        postfix:".nc",

        latitude: asr_latitude,
        longitude: asr_longitude,

        dimension_x: asr_dimension_x,
        dimension_y: asr_dimension_y,
        
        mpp_file: "NpsASR.mpp",
        rotation_angle: diff_angle};

    // conversion factors: xnew = a*x + b
    Variable u10={
        name: "U10", // U10M
        dimensions: dimensions_asr,
        a: 1.,
        b: 0.,
        Units: "m/s",
        NcVar: NcVar_tmp};

    Variable v10={
        name: "V10", // U10M
        dimensions: dimensions_asr,
        a: 1.,
        b: 0.,
        Units: "m/s",
        NcVar: NcVar_tmp};

    std::vector<Variable> variables_tmp0(2);
        variables_tmp0[0] = u10;
        variables_tmp0[1] = v10;

    M_asr_nodes_dataset={
        case_number:0,
        dirname: "data",
        prefix: "asr30km.comb.2d.", // "asr30km.comb.2D.";
        postfix:".nc",
        reference_date: "1901-01-01",
        
        nb_timestep_day: 8,
        time: asr_time, 

        dimension_time: asr_dimension_time,
        
        variables: variables_tmp0,
        target_size: M_num_nodes,
        grid: &M_asr_grid};

    M_asr_nodes_dataset.ftime_range.resize(2,0.);

    int nb_timestep_day=8;

    Variable tair={
        name:"T2",
        dimensions: dimensions_asr,
        a:1.,
        b:-273.15,
        Units:"C"}; // T2M
    Variable mixrat={
        name:"Q2",
        dimensions: dimensions_asr,
        a:1.,
        b:0.,
        Units:""}; // Q2M
    Variable mslp={
        name:"SLP",
        dimensions: dimensions_asr,
        a:1e2,
        b:0.,
        Units:"Pa"}; //PSFC, a=1.
    Variable Qsw_in={
        name:"SWDNB",
        dimensions: dimensions_asr,
        a:1.,
        b:0.,
        Units:"W/m^2"};
    Variable Qlw_in={
        name:"LWDNB",
        dimensions: dimensions_asr,
        a:1.,
        b:0.,
        Units:"W/m^2"};
    Variable snowfr={
        name:"SR",
        dimensions: dimensions_asr,
        a:1.,
        b:0.,
        Units:""};
    Variable precip={
        name:"RAINNC",
        dimensions: dimensions_asr,
        a:nb_timestep_day/(24.*3600),
        b:0.,
        Units:"kg/m^2/s"};

    std::vector<Variable> variables_tmp1(7);
    variables_tmp1[0] = tair;
    variables_tmp1[1] = mixrat;
    variables_tmp1[2] = mslp;
    variables_tmp1[3] = Qsw_in;
    variables_tmp1[4] = Qlw_in;
    variables_tmp1[5] = snowfr;
    variables_tmp1[6] = precip;

    M_asr_elements_dataset={
        case_number:1,
        dirname:"data",
        prefix:"asr30km.comb.2d.", // "asr30km.comb.2D.";
        postfix:".nc",
        reference_date: "1901-01-01",
        
        nb_timestep_day: 8,
        time: asr_time, 
        dimension_time: asr_dimension_time,

        variables: variables_tmp1,
        target_size:M_num_elements,
        grid: &M_asr_grid};

    M_asr_elements_dataset.ftime_range.resize(2,0.);

    Dimension topaz_dimension_x={
        name:"x",
        start:0,
        end:761};

    Dimension topaz_dimension_y={
        name:"y",
        start:0,
        end:1101};

    Dimension topaz_dimension_time={
        name:"time", // "Time"
        start:0,
        end:31};

    Dimension topaz_dimension_depth={
        name:"depth", // "Time"
        start:0,
        end:1};

    std::vector<Dimension> dimensions_topaz_uv(4);
    dimensions_topaz_uv[0] = topaz_dimension_time;
    dimensions_topaz_uv[1] = topaz_dimension_depth;
    dimensions_topaz_uv[2] = topaz_dimension_y;
    dimensions_topaz_uv[3] = topaz_dimension_x;

    std::vector<Dimension> dimensions_topaz(3);
    dimensions_topaz[0] = topaz_dimension_time;
    dimensions_topaz[1] = topaz_dimension_y;
    dimensions_topaz[2] = topaz_dimension_x;

    std::vector<Dimension> dimensions_topaz_latlon(2);
    dimensions_topaz_latlon[0] = topaz_dimension_y;
    dimensions_topaz_latlon[1] = topaz_dimension_x;

    std::vector<Dimension> dimensions_topaz_time(1);
    dimensions_topaz_time[0] = topaz_dimension_time;

    Variable topaz_latitude={
        name: "latitude", 
        dimensions: dimensions_topaz_latlon,
        a: 1.,
        b: 0.,
        Units: "degree_north",
        NcVar: NcVar_tmp};

    Variable topaz_longitude={
        name: "longitude", 
        dimensions: dimensions_topaz_latlon,
        a: 1.,
        b: 0.,
        Units: "degree_east",
        NcVar: NcVar_tmp};

    Variable topaz_time={
        name: "time",
        dimensions: dimensions_topaz_time,
        a: 1.,
        b: 0.,
        Units: "hours",
        NcVar: NcVar_tmp};

    M_topaz_grid={
        interpolation_method: setup::InterpolationType::InterpFromMeshToMesh2dx,
        dirname: "data",
        prefix: "TP4DAILY_",
        postfix: "_3m.nc",

        latitude: topaz_latitude,
        longitude: topaz_longitude,

        dimension_x: topaz_dimension_x,
        dimension_y: topaz_dimension_y,
        
        mpp_file: "NpsNextsim.mpp",
        rotation_angle: 0.};

    Variable u={
        name: "u",
        dimensions: dimensions_topaz_uv,
        a: 1.,
        b: 0.,
        Units: "m/s",
        NcVar: NcVar_tmp};

    Variable v={
        name: "v",
        dimensions: dimensions_topaz_uv,
        a: 1.,
        b: 0.,
        Units: "m/s",
        NcVar: NcVar_tmp};

    Variable ssh={
        name: "ssh",
        dimensions: dimensions_topaz,
        a: 1.,
        b: 0.,
        Units: "m/s",
        NcVar: NcVar_tmp};

    std::vector<Variable> variables_tmp2(3);
    variables_tmp2[0] = u;
    variables_tmp2[1] = v;
    variables_tmp2[2] = ssh;


    M_topaz_nodes_dataset={
        case_number: 2,
        dirname: "data",
        prefix: "TP4DAILY_",
        postfix: "_30m.nc",
        reference_date: "1950-01-01",
        
        nb_timestep_day: 1,

        time: topaz_time,
        dimension_time: topaz_dimension_time,
        
        variables: variables_tmp2,
        target_size: M_num_nodes,
        grid: &M_topaz_grid};

    M_topaz_nodes_dataset.ftime_range.resize(2,0.);

    Variable sst={
        name: "temperature",
        dimensions: dimensions_topaz_uv,
        a: 1.,
        b: 0.,
        Units: "deg celsius",
        NcVar: NcVar_tmp};

    Variable sss={
        name: "salinity",
        dimensions: dimensions_topaz_uv,
        a: 1.,
        b: 0.,
        Units: "",
        NcVar: NcVar_tmp};

    Variable mld={
        name: "mlp",
        dimensions: dimensions_topaz,
        a: 1.,
        b: 0.,
        Units: "m",
        NcVar: NcVar_tmp};

    std::vector<Variable> variables_tmp3(3);
    variables_tmp3[0] = sst;
    variables_tmp3[1] = sss;
    variables_tmp3[2] = mld;

    M_topaz_elements_dataset={
        case_number: 3,
        dirname: "data",
        prefix: "TP4DAILY_",
        postfix: "_3m.nc",
        reference_date: "1950-01-01",
        
        nb_timestep_day: 1,
        time: topaz_time,
        dimension_time: topaz_dimension_time,
        
        variables: variables_tmp3,
        target_size: M_num_elements,
        grid: &M_topaz_grid};

    M_topaz_elements_dataset.ftime_range.resize(2,0.);

    loadGrid(&M_asr_grid);
    loadGrid(&M_topaz_grid);

    this->initConcentration();
    this->initThickness();
    this->initDamage();
    this->initSnowThickness();
    this->initSlabOcean();

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
    time_init = dateStr2Num(vm["simul.time_init"].as<std::string>());
    output_time_step =  days_in_sec/vm["simul.output_per_day"].as<int>();

    time_step = vm["simul.timestep"].as<double>();
    duration = (vm["simul.duration"].as<double>())*days_in_sec;
    spinup_duration = (vm["simul.spinup_duration"].as<double>())*days_in_sec;

    divergence_min = (1./days_in_sec)*vm["simul.divergence_min"].as<double>();
    compression_factor = vm["simul.compression_factor"].as<double>();
    exponent_compression_factor = vm["simul.exponent_compression_factor"].as<double>();
    ocean_turning_angle_rad = (PI/180.)*vm["simul.oceanic_turning_angle"].as<double>();
    ridging_exponent = vm["simul.ridging_exponent"].as<double>();

    quad_drag_coef_air = vm["simul.ASR_quad_drag_coef_air"].as<double>();
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
                << " -" << 2 << " -part " << 1 << " -clmax " << vm["simul.hsize"].as<double>() << " " << gmshgeofile;

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

    return *std::min_element(all_min_angle.begin(),all_min_angle.end());
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
#if 0
    BamgConvertMeshx(
                     bamgmesh_previous,bamggeom_previous,
                     &M_mesh.indexTr()[0],&M_mesh.coordX()[0],&M_mesh.coordY()[0],
                     M_mesh.numNodes(), M_mesh.numTriangles()
                     );
#endif


    if(M_mesh_type==setup::MeshType::FROM_SPLIT)
    {
        if(step==0)
        {
            // step 1 (only for the first time step): Start by having bamg 'clean' the mesh with KeepVertices=0
            bamgopt->KeepVertices=0;
            this->adaptMesh();
            bamgopt->KeepVertices=1;
        }

        // Interpolate hminVertices and hmaxVertices onto the current mesh

        // NODAL INTERPOLATION
        int init_num_nodes = M_mesh_init.numNodes();
        double* interp_Vertices_in;
        interp_Vertices_in = new double[2*init_num_nodes];

        double* interp_Vertices_out;

        for (int i=0; i<init_num_nodes; ++i)
        {
            interp_Vertices_in[2*i]   = M_hminVertices[i];
            interp_Vertices_in[2*i+1] = M_hmaxVertices[i];
        }

        InterpFromMeshToMesh2dx(&interp_Vertices_out,
                            &M_mesh_init.indexTr()[0],&M_mesh_init.coordX()[0],&M_mesh_init.coordY()[0],
                            M_mesh_init.numNodes(),M_mesh_init.numTriangles(),
                            interp_Vertices_in,
                            M_mesh_init.numNodes(),2,
                            &M_mesh.coordX()[0],&M_mesh.coordY()[0],M_mesh.numNodes(),
                            false);

        bamgopt->hminVertices = new double[M_mesh.numNodes()];
        bamgopt->hmaxVertices = new double[M_mesh.numNodes()];

        for (int i=0; i<M_mesh.numNodes(); ++i)
        {
            bamgopt->hminVertices[i] = interp_Vertices_out[2*i];
            bamgopt->hmaxVertices[i] = interp_Vertices_out[2*i+1];
        }
    }

    this->adaptMesh();

    if (step)
    {
        int prv_num_elements = M_mesh_previous.numTriangles();
        int prv_num_nodes = M_mesh_previous.numNodes();

        // ELEMENT INTERPOLATION With Cavities
        double* interp_elt_in;
        int nb_var=12;
        interp_elt_in = new double[nb_var*prv_num_elements];

        double* interp_elt_out;

        std::cout<<"ELEMENT: Interp starts\n";

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

            // Ice surface temperature
            interp_elt_in[nb_var*i+tmp_nb_var] = M_tsurf[i];
            tmp_nb_var++;

            if(tmp_nb_var>nb_var)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }

#if 1
        double* surface_previous = new double[prv_num_elements];
        double* surface = new double[M_num_elements];

        int cpt = 0;
        for (auto it=M_mesh_previous.triangles().begin(), end=M_mesh_previous.triangles().end(); it!=end; ++it)
        {
            surface_previous[cpt] = this->measure(*it,M_mesh_previous);
            ++cpt;
        }

        cpt = 0;
        for (auto it=M_mesh.triangles().begin(), end=M_mesh.triangles().end(); it!=end; ++it)
        {
            surface[cpt] = M_surface[cpt];
            ++cpt;
        }

        // The interpolation with the cavities still needs to be tested on a long run.
        // By default, we then use the non-conservative MeshToMesh interpolation

        InterpFromMeshToMesh2dCavities(&interp_elt_out,interp_elt_in,nb_var,
             surface_previous, surface, bamgmesh_previous, bamgmesh);
#endif

#if 0
        InterpFromMeshToMesh2dx(&interp_elt_out,
                                &M_mesh_previous.indexTr()[0],&M_mesh_previous.coordX()[0],&M_mesh_previous.coordY()[0],
                                M_mesh_previous.numNodes(),M_mesh_previous.numTriangles(),
                                interp_elt_in,
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

        M_tsurf.assign(M_num_elements,0.);

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

            // Ice surface temperature
            M_tsurf[i] = interp_elt_out[nb_var*i+tmp_nb_var];
            tmp_nb_var++;

            if(tmp_nb_var!=nb_var)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }

        std::cout<<"ELEMENT: Interp done\n";

        // ELEMENT INTERPOLATION FOR SLAB OCEAN FROM OLD MESH ON ITS ORIGINAL POSITION
        double* interp_elt_slab_in;
        nb_var=2;
        interp_elt_slab_in = new double[nb_var*prv_num_elements];

        double* interp_elt_slab_out;

        std::cout<<"ELEMENT SLAB: Interp starts\n";

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
                                interp_elt_slab_in,
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

        M_mesh_previous.move(M_UM,displacement_factor);
        std::cout<<"ELEMENT SLAB: Interp done\n";

        // NODAL INTERPOLATION
        double* interp_in;
        interp_in = new double[6*prv_num_nodes];

        double* interp_out;

        std::cout<<"NODAL: Interp starts\n";

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

        M_VT.assign(2*M_num_nodes,0.);
        M_VTM.assign(2*M_num_nodes,0.);
        M_VTMM.assign(2*M_num_nodes,0.);
        M_UM.assign(2*M_num_nodes,0.);

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

    }

    if (step)
    {
        //M_matrix->init(2*M_num_nodes,2*M_num_nodes,22);
        M_matrix->init(2*M_num_nodes,2*M_num_nodes,M_graph);
        M_vector->resize(2*M_num_nodes);
        M_solution->resize(2*M_num_nodes);
        M_reuse_prec = false;

        M_vector_reduction.resize(2*M_num_nodes,0.);
        M_valid_conc.resize(2*M_num_nodes,false);

        M_wind.resize(2*M_num_nodes);
        M_ocean.resize(2*M_num_nodes);

        // Atmo
        M_tair.resize(M_num_elements);
        M_dair.resize(M_num_elements);
        M_mixrat.resize(M_num_elements);
        M_mslp.resize(M_num_elements);
        M_Qsw_in.resize(M_num_elements);
        M_Qlw_in.resize(M_num_elements);
        M_precip.resize(M_num_elements);
        M_snowfr.resize(M_num_elements);

        M_ocean_temp.resize(M_num_elements);
        M_ocean_salt.resize(M_num_elements);
        M_mld.resize(M_num_elements);

        M_element_depth.assign(M_num_elements,0.);
        M_h_thin.assign(M_num_elements,0.);
        M_hs_thin.assign(M_num_elements,0.);

        M_ssh.assign(M_num_nodes,0.);

        M_fcor.resize(M_num_elements);
    }

    M_asr_nodes_dataset.target_size=M_num_nodes;
    M_asr_elements_dataset.target_size=M_num_elements;
    M_topaz_nodes_dataset.target_size=M_num_nodes;
    M_topaz_elements_dataset.target_size=M_num_elements;

    M_Cohesion.resize(M_num_elements);
    M_Compressive_strength.resize(M_num_elements);
    M_time_relaxation_damage.resize(M_num_elements,time_relaxation_damage);
}

void
FiniteElement::adaptMesh()
{
    *bamgmesh_previous = *bamgmesh;
    *bamggeom_previous = *bamggeom;
    *bamgopt_previous = *bamgopt;


    int fnd = 0;
    //int snd = 0;

#if 1
    for (int edg=0; edg<bamgmesh_previous->EdgesSize[0]; ++edg)
    {
        fnd = bamgmesh_previous->Edges[3*edg]-1;
        //snd = bamgmesh_previous->Edges[3*edg+1]-1;

        // if ((std::binary_search(M_dirichlet_flags.begin(),M_dirichlet_flags.end(),fnd))
        //     || (std::binary_search(M_dirichlet_flags.begin(),M_dirichlet_flags.end(),snd)))
        if ((std::binary_search(M_dirichlet_flags.begin(),M_dirichlet_flags.end(),fnd)))
        {
            bamggeom_previous->Edges[3*edg+2] = M_flag_fix;
            bamgmesh_previous->Edges[3*edg+2] = M_flag_fix;
        }
    }

    // for (const int& edg : M_dirichlet_flags)
    // {
    //     std::cout<<"BEFORE["<< edg << "]= ("<< M_mesh.coordX()[edg] <<","<< M_mesh.coordY()[edg] <<")\n";
    // }

    //std::sort(vec.begin(), vec.end());
    // for (int i=0; i<vec.size(); ++i)
    // {
    //     //std::cout<<"VEC["<< i <<"]= "<< vec[i] <<"\n";
    //     //std::cout<<"INIT["<< vec[i] << "]= ("<< M_mesh.coordX()[vec[i]] <<","<< M_mesh.coordY()[vec[i]] <<")\n";
    // }
#endif


    //Environment::logMemoryUsage("before adaptMesh");
    Bamgx(bamgmesh,bamggeom,bamgmesh_previous,bamggeom_previous,bamgopt_previous);
    //Environment::logMemoryUsage("after adaptMesh");
    this->importBamg(bamgmesh);

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
    //std::cout<<"Boundary size 1= "<< M_boundary_flags.size() <<"\n";
    //M_boundary_flags.erase(std::unique(M_boundary_flags.begin(), M_boundary_flags.end() ), M_boundary_flags.end());
    //std::cout<<"Boundary size 2= "<< M_boundary_flags.size() <<"\n";

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
    std::cout<<"Assembling starts\n";
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
        //double _coef = ((M_thick[i]-critical_h) > 0) ? (M_thick[i]-critical_h) : 0.;
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

        double coef = young*(1.-M_damage[cpt])*tmp_thick*std::exp(ridging_exponent*(1.-tmp_conc));

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

        /* Loop over the 6 by 6 components of the finite element intergrale
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


                fvdata[2*i] += surface_e*( mloc*( coef_Vair*M_wind[index_u]
                                                +coef_Voce*cos_ocean_turning_angle*M_ocean[index_u]
                                                +coef_X
                                                +coef_V*M_VT[index_u]
                                                -coef_Voce*sin_ocean_turning_angle*(M_ocean[index_v]-M_VT[index_v])
                                                +coef_C*Vcor_index_v)
                                - b0tj_sigma_hu/3);


                fvdata[2*i+1] += surface_e*( mloc*( +coef_Vair*M_wind[index_v]
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

        if((M_conc[cpt]>0.))
        {
            for (int const& idn : rcindices)
                M_valid_conc[idn] = true;
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
    std::cout<<"SET PETSC MATRIX done in " << petsc_chrono.elapsed() <<"s\n";
    lhsdata.resize(0);

    // close petsc matrix
    M_matrix->close();

    // update petsc vector and close it
    M_vector->addVector(&rhsindices[0], rhsindices.size(), &rhsdata[0]);
    M_vector->close();

    std::cout<<"Assembling done\n";
    std::cout<<"TIMER ASSEMBLY= " << chrono.elapsed() <<"s\n";

    // extended dirichlet nodes (add nodes where M_conc <= 0)
    for (int i=0; i<2*M_num_nodes; ++i)
    {
        if(!M_valid_conc[i])
        {
            //M_matrix->setValue(i, i, 1000000000000000.);
            //M_vector->set(i, 0.);
            extended_dirichlet_nodes.push_back(i);
        }
    }
    std::sort(extended_dirichlet_nodes.begin(), extended_dirichlet_nodes.end());
    extended_dirichlet_nodes.erase(std::unique( extended_dirichlet_nodes.begin(), extended_dirichlet_nodes.end() ), extended_dirichlet_nodes.end());

    chrono.restart();
    //M_matrix->on(M_dirichlet_nodes,*M_vector);
    M_matrix->on(extended_dirichlet_nodes,*M_vector);
    std::cout<<"TIMER DBCA= " << chrono.elapsed() <<"s\n";

    std::cout<<"[PETSC MATRIX] CLOSED      = "<< M_matrix->closed() <<"\n";
    std::cout<<"[PETSC MATRIX] SIZE        = "<< M_matrix->size1() << " " << M_matrix->size2() <<"\n";
    std::cout<<"[PETSC MATRIX] SYMMETRIC   = "<< M_matrix->isSymmetric() <<"\n";
    std::cout<<"[PETSC MATRIX] NORM        = "<< M_matrix->linftyNorm() <<"\n";

    //M_matrix->printMatlab("stiffness.m");
    //M_vector->printMatlab("rhs.m");
}

void
FiniteElement::assembleSeq(int pcpt)
{
    double coef_V, coef_Voce, coef_Vair, coef_basal, coef_X, coef_Y, coef_C;
    double coef = 0;
    double coef_P = 0.;
    double mass_e = 0.;
    double surface_e = 0.;
    double g_ssh_e_x = 0.;
    double g_ssh_e = 0.;
    double g_ssh_e_y = 0.;

    double tmp_thick, tmp_conc;

    //std::vector<double> B0T(18,0);
    // std::vector<double> B0Tj_Dunit(6,0);
    // std::vector<double> B0Tj_Dunit_B0Ti(4,0);
    // std::vector<double> B0T_Dunit_B0T(36,0);
    // std::vector<double> B0Tj_Dunit_comp(6,0);
    // std::vector<double> B0Tj_Dunit_comp_B0Ti(4,0);
    // std::vector<double> B0T_Dunit_comp_B0T(36,0);

    std::vector<double> sigma_P(3,0); /* temporary variable for the resistance to the compression */
    std::vector<double> B0Tj_sigma_h(2,0);

    // double B0Tj_Dunit_tmp0, B0Tj_Dunit_tmp1;
    // double B0Tj_Dunit_B0Ti_tmp0, B0Tj_Dunit_B0Ti_tmp1, B0Tj_Dunit_B0Ti_tmp2, B0Tj_Dunit_B0Ti_tmp3;
    // double B0Tj_Dunit_comp_tmp0, B0Tj_Dunit_comp_tmp1;
    // double B0Tj_Dunit_comp_B0Ti_tmp0, B0Tj_Dunit_comp_B0Ti_tmp1, B0Tj_Dunit_comp_B0Ti_tmp2, B0Tj_Dunit_comp_B0Ti_tmp3;
    double mloc;

    std::vector<double> B0Tj_Dunit_B0Ti(4,0);
    std::vector<double> B0Tj_Dunit_comp_B0Ti(4,0);
    // std::vector<double> B0T_Dunit_B0T(36,0);
    // std::vector<double> B0T_Dunit_comp_B0T(36,0);

    std::vector<double> data(36);
    std::vector<double> fvdata(6),fvdata_reduction(6);
    std::vector<int> rcindices(6), rcindices_i(6);

    double duu, dvu, duv, dvv, fuu, fvv;
    int index_u, index_v;

    for (int i=0; i<2*M_num_nodes; ++i)
        M_vector_reduction[i]=0.;

    M_matrix->zero();
    M_vector->zero();

    std::vector<int> rhsindices(2*M_num_nodes);
    std::iota(rhsindices.begin(), rhsindices.end(), 0);
    std::vector<double> rhsdata(2*M_num_nodes, 0.);

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

    std::cout<<"Assembling starts\n";
    chrono.restart();
    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
#if 0
        //chrono.restart();
        //std::vector<double> shapecoeff = this->shapeCoeff(*it,M_mesh);
        // if (cpt==0)
        //     std::cout<<"SHAPECOEFF elapsed in "<< chrono.elapsed() <<"s\n";

        for (int i=0; i<18; ++i)
        {
            if (i < 3)
            {
                B0T[2*i] = M_shape_coeff[cpt][i];
                B0T[12+2*i] = M_shape_coeff[cpt][i+3];
                B0T[13+2*i] = M_shape_coeff[cpt][i];
            }
            else if (i < 6)
            {
                B0T[2*i+1] = M_shape_coeff[cpt][i];
            }
        }
#endif
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

        // Factor compute now within assemble
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
        //double _coef = ((M_thick[i]-critical_h) > 0) ? (M_thick[i]-critical_h) : 0.;
        double _coef = std::max(0., M_thick[cpt]-critical_h);
        double coef_basal = quad_drag_coef_air*basal_k2/(basal_drag_coef_air*(norm_Vice+basal_u_0));
        coef_basal *= _coef*std::exp(-basal_Cb*(1.-M_conc[cpt]));

        tmp_thick=(0.05>M_thick[cpt]) ? 0.05 : M_thick[cpt];
        tmp_conc=(0.01>M_conc[cpt]) ? 0.01 : M_conc[cpt];

        coef = young*(1.-M_damage[cpt])*tmp_thick*std::exp(ridging_exponent*(1.-tmp_conc));

        coef_P = 0.;
        if(M_divergence_rate[cpt] < 0.)
        {
            coef_P = compression_factor*std::pow(tmp_thick,exponent_compression_factor)*std::exp(ridging_exponent*(1.-tmp_conc));
            coef_P = coef_P/(std::abs(M_divergence_rate[cpt])+divergence_min);
            //std::cout<<"Coeff= "<< coef_P <<"\n";
        }

        sigma_P[0]=0.;
        sigma_P[1]=0.;
        sigma_P[2]=0.;

        /* Compute the value that only depends on the element */
        mass_e = rhoi*tmp_thick + rhos*M_snow_thick[cpt];
        mass_e = (tmp_conc > 0.) ? (mass_e/tmp_conc):0.;
        surface_e = M_surface[cpt];

        // /* compute the x and y derivative of g*ssh */
        g_ssh_e_x = 0.;
        g_ssh_e_y = 0.;
        for(int i=0; i<3; i++)
        {
            g_ssh_e = (vm["simul.gravity"].as<double>())*M_ssh[it->indices[i]-1] /*g_ssh*/;   /* g*ssh at the node k of the element e */
            g_ssh_e_x += M_shape_coeff[cpt][i]*g_ssh_e; /* x derivative of g*ssh */
            g_ssh_e_y += M_shape_coeff[cpt][i+3]*g_ssh_e; /* y derivative of g*ssh */
        }

        coef_C     = mass_e*M_fcor[cpt];              /* for the Coriolis term */
        coef_V     = mass_e/time_step;             /* for the inertial term */
        coef_X     = - mass_e*g_ssh_e_x;              /* for the ocean slope */
        coef_Y     = - mass_e*g_ssh_e_y;              /* for the ocean slope */

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
#if 0
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
                    B0Tj_Dunit_tmp0 += M_B0T[cpt][kk*6+2*j]*M_Dunit[3*i+kk];
                    B0Tj_Dunit_tmp1 += M_B0T[cpt][kk*6+2*j+1]*M_Dunit[3*i+kk];

                    B0Tj_Dunit_comp_tmp0 += M_B0T[cpt][kk*6+2*j]*M_Dunit_comp[3*i+kk];
                    B0Tj_Dunit_comp_tmp1 += M_B0T[cpt][kk*6+2*j+1]*M_Dunit_comp[3*i+kk];
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
                    B0Tj_Dunit_B0Ti_tmp0 += B0Tj_Dunit[2*kk]*M_B0T[cpt][kk*6+2*i];
                    B0Tj_Dunit_B0Ti_tmp1 += B0Tj_Dunit[2*kk]*M_B0T[cpt][kk*6+2*i+1];
                    B0Tj_Dunit_B0Ti_tmp2 += B0Tj_Dunit[2*kk+1]*M_B0T[cpt][kk*6+2*i];
                    B0Tj_Dunit_B0Ti_tmp3 += B0Tj_Dunit[2*kk+1]*M_B0T[cpt][kk*6+2*i+1];

                    B0Tj_Dunit_comp_B0Ti_tmp0 += B0Tj_Dunit_comp[2*kk]*M_B0T[cpt][kk*6+2*i];
                    B0Tj_Dunit_comp_B0Ti_tmp1 += B0Tj_Dunit_comp[2*kk]*M_B0T[cpt][kk*6+2*i+1];
                    B0Tj_Dunit_comp_B0Ti_tmp2 += B0Tj_Dunit_comp[2*kk+1]*M_B0T[cpt][kk*6+2*i];
                    B0Tj_Dunit_comp_B0Ti_tmp3 += B0Tj_Dunit_comp[2*kk+1]*M_B0T[cpt][kk*6+2*i+1];
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
#endif
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

        for(int i=0; i<3; i++)
        {
            fvdata[2*i]=0.;
            fvdata[2*i+1]=0.;
        }

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
                int index_u_i = it->indices[i]-1;
                int index_v_i = it->indices[i]-1+M_num_nodes;
                /* Row corresponding to indice i (we also assemble terms in row+1) */
                //row = (mwIndex)it[2*i]-1; /* -1 to use the indice convention of C */

                B0Tj_Dunit_B0Ti[0] = M_B0T_Dunit_B0T[cpt][(2*i)*6+2*j];
                B0Tj_Dunit_B0Ti[1] = M_B0T_Dunit_B0T[cpt][(2*i+1)*6+2*j];
                B0Tj_Dunit_B0Ti[2] = M_B0T_Dunit_B0T[cpt][(2*i)*6+2*j+1];
                B0Tj_Dunit_B0Ti[3] = M_B0T_Dunit_B0T[cpt][(2*i+1)*6+2*j+1];

                B0Tj_Dunit_comp_B0Ti[0] = M_B0T_Dunit_comp_B0T[cpt][(2*i)*6+2*j];
                B0Tj_Dunit_comp_B0Ti[1] = M_B0T_Dunit_comp_B0T[cpt][(2*i+1)*6+2*j];
                B0Tj_Dunit_comp_B0Ti[2] = M_B0T_Dunit_comp_B0T[cpt][(2*i)*6+2*j+1];
                B0Tj_Dunit_comp_B0Ti[3] = M_B0T_Dunit_comp_B0T[cpt][(2*i+1)*6+2*j+1];

                /* Select the nodal weight values from M_loc */
                mloc = M_Mass[3*j+i];

                B0Tj_sigma_h[0] = 0.;
                B0Tj_sigma_h[1] = 0.;

                for(int k=0; k<3; k++)
                {
                    B0Tj_sigma_h[0] += M_B0T[cpt][k*6+2*i]*(M_sigma[3*cpt+k]*tmp_thick+sigma_P[k]);
                    B0Tj_sigma_h[1] += M_B0T[cpt][k*6+2*i+1]*(M_sigma[3*cpt+k]*tmp_thick+sigma_P[k]);
                }

                /* ---------- UU component */
                duu = surface_e*( mloc*(coef_Vair+coef_Voce*cos_ocean_turning_angle+coef_V+coef_basal)+B0Tj_Dunit_B0Ti[0]*coef*time_step+B0Tj_Dunit_comp_B0Ti[0]*coef_P);

                /* ---------- VU component */
                dvu = surface_e*(+B0Tj_Dunit_B0Ti[1]*coef*time_step+B0Tj_Dunit_comp_B0Ti[1]*coef_P);

                /* ---------- UV component */
                duv = surface_e*(+B0Tj_Dunit_B0Ti[2]*coef*time_step+B0Tj_Dunit_comp_B0Ti[2]*coef_P);

                /* ---------- VV component */
                dvv = surface_e*( mloc*(coef_Vair+coef_Voce*cos_ocean_turning_angle+coef_V+coef_basal)+B0Tj_Dunit_B0Ti[3]*coef*time_step+B0Tj_Dunit_comp_B0Ti[3]*coef_P);

                // if (cpt == 1)
                // {
                //     std::cout<<"duu= "<< duu <<"\n";
                //     std::cout<<"dvu= "<< dvu <<"\n";
                //     std::cout<<"duv= "<< duv <<"\n";
                //     std::cout<<"dvv= "<< dvv <<"\n";
                // }

                // data[12*i+2*j] = duu;
                // data[12*i+2*j+1] = duv;
                // data[(12*i+2*j)+6] = dvu;
                // data[(12*i+2*j+1)+6] = dvv;

                data[(2*i  )*6+2*j  ] = duu;
                data[(2*i+1)*6+2*j  ] = dvu;
                data[(2*i  )*6+2*j+1] = duv;
                data[(2*i+1)*6+2*j+1] = dvv;

#if 0
                fuu += surface_e*( mloc*( coef_Vair*M_wind[index_u]+coef_Voce*cos_ocean_turning_angle*M_ocean[index_u]+coef_X+coef_V*M_VT[index_u]) - B0Tj_sigma_h[0]/3);
                fuu += surface_e*( mloc*( -coef_Voce*sin_ocean_turning_angle*(M_ocean[index_u]-M_VT[index_u])-coef_C*M_Vcor[index_u]) );

                fvv += surface_e*( mloc*( coef_Vair*M_wind[index_v]+coef_Voce*cos_ocean_turning_angle*M_ocean[index_v]+coef_Y+coef_V*M_VT[index_v]) - B0Tj_sigma_h[1]/3);
                fvv += surface_e*( mloc*( -coef_Voce*sin_ocean_turning_angle*(M_ocean[index_v]-M_VT[index_v])-coef_C*M_Vcor[index_v]) );
#endif

                double Vcor_index_u=beta0*M_VT[index_u] + beta1*M_VTM[index_u] + beta2*M_VTMM[index_u];

                fvdata[2*i] += surface_e*( mloc*( coef_Vair*M_wind[index_u]+coef_Voce*cos_ocean_turning_angle*M_ocean[index_u]+coef_X+coef_V*M_VT[index_u]) - B0Tj_sigma_h[0]/3);
                fvdata[2*i+1] += surface_e*( mloc*( +coef_Voce*sin_ocean_turning_angle*(M_ocean[index_u]-M_VT[index_u])-coef_C*Vcor_index_u) );

                double Vcor_index_v=beta0*M_VT[index_v] + beta1*M_VTM[index_v] + beta2*M_VTMM[index_v];

                fvdata[2*i] += surface_e*( mloc*( -coef_Voce*sin_ocean_turning_angle*(M_ocean[index_v]-M_VT[index_v])+coef_C*Vcor_index_v) );
                fvdata[2*i+1] += surface_e*( mloc*( coef_Vair*M_wind[index_v]+coef_Voce*cos_ocean_turning_angle*M_ocean[index_v]+coef_Y+coef_V*M_VT[index_v]) - B0Tj_sigma_h[1]/3);


                rcindices_i[2*i] = index_u_i;
                rcindices_i[2*i+1] = index_v_i;
                #if 0
                M_matrix->addValue(index_u_i,index_u, duu);
                M_matrix->addValue(index_v_i,index_u, dvu);

                M_matrix->addValue(index_u_i,index_v, duv);
                M_matrix->addValue(index_v_i,index_v, dvv);


                M_vector->add(index_u_i, fvdata[2*i]);
                M_vector->add(index_v_i, fvdata[2*i+1]);
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

            //fvdata[2*j] = fuu;
            //fvdata[2*j+1] = fvv;

            fvdata_reduction[2*j] = 1.;
            fvdata_reduction[2*j+1] = 1.;

            rcindices[2*j] = index_u;
            rcindices[2*j+1] = index_v;

            //M_vector->addVector(&rcindices_i[0], rcindices_i.size(), &fvdata[0]);
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

        M_matrix->addMatrix(&rcindices_i[0], rcindices_i.size(),
                            &rcindices[0], rcindices.size(), &data[0]);


        //M_vector->addVector(&rcindices_i[0], rcindices_i.size(), &fvdata[0]);

#if 1
        for (int idf=0; idf<rcindices_i.size(); ++idf)
        {
            rhsdata[rcindices_i[idf]] += fvdata[idf];
        }
#endif
        if((M_conc[cpt]>0.))
        {
            M_vector_reduction[rcindices[0]]=1.;
            M_vector_reduction[rcindices[1]]=1.;
            M_vector_reduction[rcindices[2]]=1.;
            M_vector_reduction[rcindices[3]]=1.;
            M_vector_reduction[rcindices[4]]=1.;
            M_vector_reduction[rcindices[5]]=1.;
        }

        ++cpt;
    }



    M_vector->addVector(&rhsindices[0], rhsindices.size(), &rhsdata[0]);

    for (int i=0; i<2*M_num_nodes; ++i)
    {
        if(M_vector_reduction[i]==0.)
        {
            M_matrix->setValue(i, i, 1000000000000000.);

            M_vector->set(i, 0.);
        }
    }

    M_vector->close();
    //M_vector->printMatlab("vectorseq.m");

    M_matrix->close();
    //M_matrix->printMatlab("matrixseq.m");
    std::cout<<"Assembling done\n";
    std::cout<<"TIMER ASSEMBLY= " << chrono.elapsed() <<"s\n";

    chrono.restart();
    M_matrix->on(M_dirichlet_nodes,*M_vector);
    std::cout<<"TIMER DBCA= " << chrono.elapsed() <<"s\n";

    M_matrix->printMatlab("matrixseq.m");
    M_vector->printMatlab("vectorseq.m");

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

    std::cout<<"VT MIN= "<< *std::min_element(M_VT.begin(),M_VT.end()) <<"\n";
    std::cout<<"VT MAX= "<< *std::max_element(M_VT.begin(),M_VT.end()) <<"\n";

    int thread_id;
    int total_threads;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    //std::cout<<"MAX THREADS= "<< max_threads <<"\n";

#pragma omp parallel for num_threads(max_threads) private(thread_id)
    for (int cpt=0; cpt < M_num_elements; ++cpt)
    {
        double old_thick;
        double old_snow_thick;
        double old_conc;
        double old_damage;
        double old_h_ridged_thick_ice;

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
        double sigma_s, sigma_n;
        double tract_max;
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

        if(M_ice_type==setup::IceCategoryType::THIN_ICE)
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

        for(i=0;i<3;i++)
        {
            sigma_dot_i = 0.0;
            for(j=0;j<3;j++)
            {
                sigma_dot_i += std::exp(ridging_exponent*(1.-old_conc))*young*(1.-old_damage)*M_Dunit[i*3 + j]*epsilon_veloc[j];
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
        sigma_n=           (sigma_pred[0]+sigma_pred[1])/2.;

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

            tmp=1.0-sigma_target/sigma_n*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }

        if(sigma_s>M_Cohesion[cpt]-sigma_n*tan_phi)
        {
            tmp=1.0-M_Cohesion[cpt]/(sigma_s+sigma_n*tan_phi)*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }

        /*
         * Diagnostic:
         * Recompute the internal stress
         */
        for(i=0;i<3;i++)
        {
            if(old_damage<1.0)
            {
                M_sigma[3*cpt+i] = (1.-M_damage[cpt])/(1.-old_damage)*M_sigma[3*cpt+i] ;
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

            M_conc[cpt]    = ice_surface/surface_new;
            M_thick[cpt]   = ice_volume/surface_new; // Hold on! Isn't M_thick the effective thickness?
            M_snow_thick[cpt]   = snow_volume/surface_new;
            M_h_ridged_thick_ice[cpt]   =   ridged_thick_ice_volume/surface_new;

            if(M_ice_type==setup::IceCategoryType::THIN_ICE)
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
                if(M_ice_type==setup::IceCategoryType::THIN_ICE)
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

            if(M_ice_type==setup::IceCategoryType::THIN_ICE)
            {
                M_h_thin[cpt]    = ((M_h_thin[cpt]>0.)?(M_h_thin[cpt] ):(0.)) ;
                M_hs_thin[cpt]   = ((M_hs_thin[cpt]>0.)?(M_hs_thin[cpt]):(0.)) ;
            }
        }

        if(M_ice_type==setup::IceCategoryType::THIN_ICE)
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
    }
}

void
FiniteElement::updateSeq()
{
    double old_thick;
    double old_snow_thick;
    double old_conc;
    double old_damage;
    double old_h_ridged_thick_ice;

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
    double sigma_s, sigma_n;
    double tract_max;
    double tmp, sigma_target;

    /* some variables used for the ice redistribution*/
    double tanalpha, rtanalpha, del_v, del_c, del_vs, new_v_thin;

    /* Indices loop*/
    int i,j;

    bool to_be_updated;

    /* set constants for the ice redistribution */
	tanalpha  = h_thin_max/c_thin_max;
	rtanalpha = 1./tanalpha;

    std::vector<double> UM_P = M_UM;

    //std::cout<<"CURRENT: DIRICHLET NODES= "<< M_dirichlet_flags.size() <<"\n";
    //std::cout<<"CURRENT: NEUMANN   NODES= "<< M_neumann_flags.size() <<"\n";

    // std::cout<<"UM MIN B= "<< *std::min_element(M_UM.begin(),M_UM.end()) <<"\n";
    // std::cout<<"UM MAX B= "<< *std::max_element(M_UM.begin(),M_UM.end()) <<"\n";

    for (int nd=0; nd<M_UM.size(); ++nd)
    {
        M_UM[nd] += time_step*M_VT[nd];
    }

    for (const int& nd : M_neumann_nodes)
    {
        M_UM[nd] = UM_P[nd];
    }

    // std::cout<<"UM MIN A= "<< *std::min_element(M_UM.begin(),M_UM.end()) <<"\n";
    // std::cout<<"UM MAX A= "<< *std::max_element(M_UM.begin(),M_UM.end()) <<"\n";

    // std::cout<<"VT MIN= "<< *std::min_element(M_VT.begin(),M_VT.end()) <<"\n";
    // std::cout<<"VT MAX= "<< *std::max_element(M_VT.begin(),M_VT.end()) <<"\n";

    // std::cout<<"Jacobian init   = "<< this->measure(*M_elements.begin(),M_mesh, UM_P) <<"\n";
    // std::cout<<"Jacobian current= "<< this->measure(*M_elements.begin(),M_mesh,M_UM) <<"\n";

    int cpt = 0;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        // Temporary memory
        old_thick = M_thick[cpt];
        old_snow_thick = M_snow_thick[cpt];
        old_conc = M_conc[cpt];
        old_damage = M_damage[cpt];
        old_h_ridged_thick_ice=M_h_ridged_thick_ice[cpt];

        if(M_ice_type==setup::IceCategoryType::THIN_ICE)
        {
            old_h_thin = M_h_thin[cpt];
            old_hs_thin=M_hs_thin[cpt];
            old_h_ridged_thin_ice=M_h_ridged_thin_ice[cpt];
        }

#if 0
        std::vector<double> shapecoeff = this->shapeCoeff(*it,M_mesh);

        for (i=0; i<18; ++i)
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
#endif

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
                epsilon_veloc_i += M_B0T[cpt][i*6 + 2*j]*M_VT[it->indices[j]-1]  ;
                epsilon_veloc_i += M_B0T[cpt][i*6 + 2*j + 1]*M_VT[it->indices[j]-1+M_num_nodes]  ;
            }

            epsilon_veloc[i] = epsilon_veloc_i;
        }

        M_divergence_rate[cpt]= (epsilon_veloc[0]+epsilon_veloc[1]);

        //std::cout<<"Rate= "<< M_divergence_rate[cpt] <<"\n";

        /*======================================================================
         * Update the internal stress
         *======================================================================
         */

        for(i=0;i<3;i++)
        {
            sigma_dot_i = 0.0;
            for(j=0;j<3;j++)
            {
                sigma_dot_i += std::exp(ridging_exponent*(1.-old_conc))*young*(1.-old_damage)*M_Dunit[i*3 + j]*epsilon_veloc[j];
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

        sigma_s=std::hypot((sigma_pred[0]-sigma_pred[1])/2.,sigma_pred[2]);
        sigma_n=           (sigma_pred[0]+sigma_pred[1])/2.;

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

            tmp=1.0-sigma_target/sigma_n*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }

        if(sigma_s>M_Cohesion[cpt]-sigma_n*tan_phi)
        {
            tmp=1.0-M_Cohesion[cpt]/(sigma_s+sigma_n*tan_phi)*(1.0-old_damage);

            if(tmp>M_damage[cpt])
            {
                M_damage[cpt]=tmp;
            }
        }

        /*
         * Diagnostic:
         * Recompute the internal stress
         */
        for(i=0;i<3;i++)
        {
            if(old_damage<1.0)
            {
                M_sigma[3*cpt+i] = (1.-M_damage[cpt])/(1.-old_damage)*M_sigma[3*cpt+i] ;
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

        /* For the Lagrangian scheme, we do not update the variables for the elements having one node on the open boundary. */
        if(std::find(M_neumann_flags.begin(),M_neumann_flags.end(),it->indices[0]-1) != M_neumann_flags.end() ||
           std::find(M_neumann_flags.begin(),M_neumann_flags.end(),it->indices[1]-1) != M_neumann_flags.end() ||
           std::find(M_neumann_flags.begin(),M_neumann_flags.end(),it->indices[2]-1) != M_neumann_flags.end())
            to_be_updated=false;

        if((old_conc>0.)  && to_be_updated)
        {
            surface = this->measure(*it,M_mesh, UM_P);
            surface_new = this->measure(*it,M_mesh,M_UM);

            // std::cout<<"SURFACE    = "<< std::setprecision(18) << surface <<"\n";
            // std::cout<<"SURFACE_NEW= "<< std::setprecision(18) << surface_new <<"\n";

            ice_surface = old_conc*surface;
            ice_volume = old_thick*surface;
            snow_volume = old_snow_thick*surface;
            ridged_thick_ice_volume = old_h_ridged_thick_ice*surface;

            M_conc[cpt]    = ice_surface/surface_new;
            M_thick[cpt]   = ice_volume/surface_new;
            M_snow_thick[cpt]   = snow_volume/surface_new;
            M_h_ridged_thick_ice[cpt]   =   ridged_thick_ice_volume/surface_new;

            if(M_ice_type==setup::IceCategoryType::THIN_ICE)
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
                if(M_ice_type==setup::IceCategoryType::THIN_ICE)
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

            if(M_ice_type==setup::IceCategoryType::THIN_ICE)
            {
                M_h_thin[cpt]    = ((M_h_thin[cpt]>0.)?(M_h_thin[cpt] ):(0.)) ;
                M_hs_thin[cpt]   = ((M_hs_thin[cpt]>0.)?(M_hs_thin[cpt]):(0.)) ;
            }
        }

        if(M_ice_type==setup::IceCategoryType::THIN_ICE)
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

        ++cpt;
    }

    // std::cout<<"****************************\n";
    // std::cout<<"SIGMA MIN= "<< *std::min_element(M_sigma.begin(),M_sigma.end()) <<"\n";
    // std::cout<<"SIGMA MAX= "<< *std::max_element(M_sigma.begin(),M_sigma.end()) <<"\n";

#if 0
    std::cout<<"****************************\n";
    std::cout<<"THICK MIN= "<< *std::min_element(M_thick.begin(),M_thick.end()) <<"\n";
    std::cout<<"THICK MAX= "<< *std::max_element(M_thick.begin(),M_thick.end()) <<"\n";

    auto it = std::max_element(M_thick.begin(),M_thick.end());
    std::cout << "MAX ELEMENT AT: " << std::distance(M_thick.begin(), it) << '\n';
#endif

}

void
FiniteElement::solve()
{
    //SolverPetsc ksp;
    chrono.restart();
    //ksp.solve(M_matrix, M_solution, M_vector);
    M_solver->solve(_matrix=M_matrix,
                    _solution=M_solution,
                    _rhs=M_vector,
                    _ksp=vm["solver.ksp-type"].as<std::string>()/*"preonly"*/,
                    _pc=vm["solver.pc-type"].as<std::string>()/*"cholesky"*/,
                    _pcfactormatsolverpackage=vm["solver.mat-package-type"].as<std::string>()/*"cholmod"*/,
                    _reuse_prec=true,
                    _rebuild=M_regrid
                    );

    std::cout<<"TIMER SOLUTION= " << chrono.elapsed() <<"s\n";

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
    double const drag_ice_t=vm["simul.drag_ice_t"].as<double>();

    double const alb_ice    = vm["simul.alb_ice"].as<double>();
    double const alb_sn     = vm["simul.alb_sn"].as<double>();
    double const I_0        = vm["simul.I_0"].as<double>();
    int const alb_scheme = vm["simul.alb_scheme"].as<int>();

    double const rh0   = 1./vm["simul.hnull"].as<double>();
    double const rPhiF = 1./vm["simul.PhiF"].as<double>();

    bool const flooding = vm["simul.flooding"].as<bool>();
    
    /*double const tanalpha  = vm["simul.hi_thin_max"].as<double>()/vm["simul.c_thin_max"].as<double>();
    double const rtanalpha = 1./tanalpha;*/

    double const qi = physical::Lf * physical::rhoi;
    double const qs = physical::Lf * physical::rhos;

    int const newice_type = vm["simul.newice_type"].as<int>();
    int const melt_type = vm["simul.melt_type"].as<int>();
    double const PhiM = vm["simul.PhiM"].as<double>();
    double const PhiF = vm["simul.PhiF"].as<double>();

    const double aw=6.1121e2, bw=18.729, cw=257.87, dw=227.3;
    const double Aw=7.2e-4, Bw=3.20e-6, Cw=5.9e-10;

    const double ai=6.1115e2, bi=23.036, ci=279.82, di=333.7;
    const double Ai=2.2e-4, Bi=3.83e-6, Ci=6.4e-10;

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
        double  hs=0.;     // Snow thickness (slab)

        double  del_hi=0.; // Change in ice thickness (slab only)

        double  evap=0.;   // Evaporation

        double  Qdw=0.;    // Heat flux from ocean nudging
        double  Fdw=0.;    // Fresh water flux from ocean nudging

        double  Qio=0.;    // Ice-ocean heat flux
        double  Qai=0.;    // Atmosphere-ice heat flux
        double  Qow=0.;    // Open water heat flux

        // Save old _volumes_ and concentration and calculate wind speed
        // SYL: % thin ice missing??
        double  old_vol=M_thick[i];
        double  old_snow_vol=M_snow_thick[i];
        double  old_conc=M_conc[i];

        double sum_u=0.;
        double sum_v=0.;

        for (int j=0; j<3; ++j)
        {
            // calculate wind per element
            sum_u += M_wind[M_elements[i].indices[j]-1];
            sum_v += M_wind[M_elements[i].indices[j]-1+M_num_nodes];
        }
        double  wspeed = std::hypot(sum_u, sum_v)/3.;

        // -------------------------------------------------
        // 2) We calculate or set the flux due to nudging
        if ( M_atmosphere_type == setup::AtmosphereType::CONSTANT || M_ocean_type == setup::OceanType::CONSTANT )
        {
            Qdw=Qdw_const;
            Fdw=Fdw_const;
        }
        else
        {
            // nudgeFlux
            if ( M_ocean_salt[i] > physical::si )
            {
                Qdw = -(M_sst[i]-M_ocean_temp[i]) * M_mld[i] * physical::rhow * physical::cpw/timeT;

                double delS = M_sss[i] - M_ocean_salt[i];
                Fdw = delS * M_mld[i] * physical::rhow /(timeS*M_sss[i] - time_step*delS);
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
        Qow = -M_Qsw_in[i]*(1.-ocean_albedo) - M_Qlw_in[i] + Qlw_out + Qsh + Qlh;

        // -------------------------------------------------
        // 4) Thickness change of the ice slab (thermoIce0 in matlab)

        /* Don't do anything if there's no ice */
        if ( M_conc[i] <=0. )
        {
            hi      = 0.;
            hs      = 0.;
            M_tsurf[i] = 0.;
            Qio     = 0.;
            del_hi  = 0.;
            Qai     = 0.;
        } else {
            /* Calculate the slab thickness */
            hi = M_thick[i]/M_conc[i];
            hs = M_snow_thick[i]/M_conc[i];

            /* Local variables */
            double albedo;

            double albs, albi, frac_sn;

            double Qsw, Qout, dQaidT, subl;
            double Qic, del_hs, del_ht, del_hb, draft;

            double Qlw_out, dQlwdT;
            double tairK, sphuma, sphumi;
            double rhoair, Qsh, dQshdT;
            double Qlh, dsphumidT, dQlhdT;

            double fi, esti;
            double dsphumdesti, destidT, dfidT;

            /* ---------------------------------------------------------------
            * Calculate the surface temperature within a while-loop
            * --------------------------------------------------------------- */
            double dtsurf   = 1.;
            double tbot     = physical::mu*M_sss[i];
            int nb_iter_while=0;
            while ( dtsurf > 1e-4 )
            {
                nb_iter_while++;
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
                        /* Albedo scheme from ccsm 3 */
                        /* The scheme is simplified using the assumption that visible solar
                         * radiation accounts for 52% of the spectrum and the near-infrared for
                         * 48% (the same split as used in cice when running in stand-alone
                         * mode). */

                        /* This is the ccsm3 scheme when alb_ice = 0.538 and alb_sn = 0.8256 */
                        if ( M_tsurf[i] > -1. )
                        {
                            albi = alb_ice - 0.075*(M_tsurf[i]+1.);
                            albs = alb_sn  - 0.124*(M_tsurf[i]+1.);
                        } else {
                            albi = alb_ice;
                            albs = alb_sn;
                        }

                        /* Snow cover fraction */
                        frac_sn = hs/(hs+0.02);

                        /* Final albedo */
                        albedo = frac_sn*albs + (1.-frac_sn)*albi;

                        break;
                    default:
                        std::cout << "alb_scheme = " << alb_scheme << "\n";
                        throw std::logic_error("Wrong albedo_scheme");
                }

                /* Calculate atmospheric fluxes */
                Qsw = -M_Qsw_in[i]*(1.-albedo);

                // -------------------------------------------------
                // 4.1) BULK ICE                
                /* Out-going long-wave flux and derivative */
                Qlw_out =   physical::eps * physical::sigma_sb * std::pow(M_tsurf[i]+physical::tfrwK,4);
                dQlwdT  = 4.*physical::eps * physical::sigma_sb * std::pow(M_tsurf[i]+physical::tfrwK,3);

                // -------------------------------------------------
                // calcSphumI
                /* Specific humidity - ice surface */
                fi      = 1. + Ai + M_mslp[i]*1e-2*( Bi + Ci*M_tsurf[i]*M_tsurf[i] );
                esti    = ai*std::exp( (bi-M_tsurf[i]/di)*M_tsurf[i]/(M_tsurf[i]+ci) );
                sphumi = alpha*fi*esti/(M_mslp[i]-beta*fi*esti);

                /* We need the derivative of sphumi wrt. tsurf */
                dsphumdesti = alpha/(M_mslp[i]-beta*fi*esti)*( 1. + beta*fi*esti/(M_mslp[i]-beta*fi*esti) );
                destidT     = ( bi*ci*di-M_tsurf[i]*( 2.*ci+M_tsurf[i]) )/( di*std::pow(ci+M_tsurf[i],2) )*esti;
                dfidT       = 2.*Ci*Bi*M_tsurf[i];
                dsphumidT   = dsphumdesti*(fi*destidT+esti*dfidT);

                // -------------------------------------------------

                /* Density of air */
                tairK  = M_tair[i] + physical::tfrwK;
                rhoair = M_mslp[i]/(Ra*tairK) * (1.+sphuma)/(1.+1.609*sphuma);

                /* Sensible heat flux and derivative */
                Qsh    = drag_ice_t * rhoair * physical::cpa * wspeed*( M_tsurf[i] - M_tair[i] );
                dQshdT = drag_ice_t * rhoair * physical::cpa * wspeed;

                /* Latent heat flux and derivative */
                Qlh    = drag_ice_t*rhoair*(physical::Lf+physical::Lv0)*wspeed*( sphumi - sphuma );
                dQlhdT = drag_ice_t*(physical::Lf+physical::Lv0)*rhoair*wspeed*dsphumidT;

                /* Sum them up */
                Qout    = Qlw_out + Qsh + Qlh;
                dQaidT = dQlwdT + dQshdT + dQlhdT;

                /* Sublimation */
                subl    = Qlh/(physical::Lf+physical::Lv0);

                /* Sum them up */
                Qai = Qsw - M_Qlw_in[i] + Qout;

                // -------------------------------------------------
                /* Recalculate M_tsurf */
                dtsurf = M_tsurf[i];
                Qic    = physical::ks*( tbot-M_tsurf[i] )/( hs + physical::ks*hi/physical::ki );
                M_tsurf[i] = M_tsurf[i] + ( Qic - Qai )/
                    ( physical::ks/(hs+physical::ks*hi/physical::ki) + dQaidT );

                /* Set M_tsurf to the freezing point of snow or ice */
                if ( hs > 0. )
                    M_tsurf[i] = std::min(0., M_tsurf[i]);
                else
                    M_tsurf[i] = std::min(physical::mu*physical::si, M_tsurf[i]);

                /* Re-evaluate the exit condition */
                dtsurf = std::abs(dtsurf-M_tsurf[i]);
            }

            if(nb_iter_while>10)
            {
                std::cout << "nb_iter_while = " << nb_iter_while << "\n";
                throw std::logic_error("nb_iter_while larger than 10");
            }

            /* Conductive flux through the ice */
            Qic = physical::ks*( tbot-M_tsurf[i] )/( hs + physical::ks*hi/physical::ki );

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
            hs  = hs + del_hs + M_precip[i]*M_snowfr[i]/physical::rhos*time_step;

            /* Heatflux from ocean */
            /* Use all excess heat to melt or grow ice. This is not
             * accurate, but will have to do for now! */
            Qio = (M_sst[i]-tbot)*physical::rhow*physical::cpw*M_mld[i]/time_step;
            /* Bottom melt/growth */
            del_hb = (Qic-Qio)*time_step/qi;

            /* Combine top and bottom */
            del_hi = del_ht+del_hb;
            hi     = hi + del_hi;

            /* Make sure we don't get too small hi_new */
            if ( hi < hmin )
            {
                del_hi  = del_hi-hi;
                Qio     = Qio + hi*qi/time_step + hs*qs/time_step;

                hi      = 0.;
                hs      = 0.;
                M_tsurf[i] = 0.;
            }

            /* Snow-to-ice conversion */
            draft = ( hi*physical::rhoi + hs*physical::rhos ) / physical::rhow;
            if ( flooding && draft > hi )
            {
                /* Subtract the mass of snow converted to ice from hs_new */
                hs = hs - ( draft - hi )*physical::rhoi/physical::rhos;
                hi = draft;
            }
        }

        // -------------------------------------------------
        // 5) Ice growth over open water and lateral melt (thermoOW in matlab)

        /* Local variables */
        double hi_old, tw_new, tfrw, newice, del_c, newsnow, h0;
        
        /* dT/dt due to heatflux atmos.-ocean */
        tw_new = M_sst[i] - Qow*time_step/(M_mld[i]*physical::rhow*physical::cpw);
        tfrw   = physical::mu*M_sss[i];

        /* Form new ice in case of super cooling, and reset Qow and evap */
        if ( tw_new < tfrw )
        {
            newice  = (1.-M_conc[i])*(tfrw-tw_new)*M_mld[i]*physical::rhow*physical::cpw/qi;
            Qow  = -(tfrw-M_sst[i])*M_mld[i]*physical::rhow*physical::cpw/time_step;
            evap = 0.;
        } else {
            newice  = 0.;
        }

        /* Decide the change in ice fraction (del_c) */
        hi_old = hi-del_hi;
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
                    /* Thin ice category -- not for now!
                    thinIceRedistribute(v_thin[i], vs_thin[i], newice/(1-c[i]), c[i], tanalpha, rtanalpha, hi_thin_max, v_thin_new[i], newice, del_c, newsnow);
                    // Change the snow _thickness_ for thick ice and _volume_ for thin ice
                    vs_thin_new[i] = vs_thin[i] - newsnow;
                    break; */
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
                    /* + Deliver excess energy to the ocean when there's no ice left */
                    Qow = (1.-PhiM)*Qow
                            + std::min(0., std::max(0.,M_conc[i]+del_c)*( hi*qi+hs*qs )/time_step);
                    /* Don't suffer negative c! */
                    del_c = std::max(del_c, -M_conc[i]);
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
            hi = ( hi*(M_conc[i]-del_c) + newice )/M_conc[i];
            if ( del_c < 0. )
            {
                /* We conserve the snow height, but melt away snow as the concentration decreases */
                Qow = Qow + del_c*hs*qs/time_step;
            } else {
                /* Snow volume is conserved as concentration increases */
                hs  = ( hs*(M_conc[i]-del_c) + newsnow )/M_conc[i];
            }
        }

        /* Check limits */
        if ( M_conc[i] < cmin || hi < hmin )
        {
            //Qow    = Qow + (M_conc[i]-del_c)*hi*qi/time_step + (M_conc[i]-del_c)*hs*qs/time_step; 
            //SYL: is this commented line right??
            // I think irt is wrong as we already modify Qow, I would do:
            Qow    = Qow + M_conc[i]*hi*qi/time_step + M_conc[i]*hs*qs/time_step;
            M_conc[i] = 0.;
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
        del_vi = M_thick[i] - old_vol;
        del_vs = M_snow_thick[i] - old_snow_vol;

        // Rain falling on ice falls straight through. We need to calculate the
        // bulk freshwater input into the entire cell, i.e. everything in the
        // open-water part plus rain in the ice-covered part.
        rain = (1.-old_conc)*M_precip[i] + old_conc*(1.-M_snowfr[i])*M_precip[i];
        emp  = (evap*(1.-old_conc)-rain);

        Qio_mean = Qio*old_conc;
        Qow_mean = Qow*(1.-old_conc);

        /* Heat-flux */
        M_sst[i] = M_sst[i] - time_step*( Qio_mean + Qow_mean - Qdw )/(physical::rhow*physical::cpw*M_mld[i]);

        /* Change in salinity */
        M_sss[i] = M_sss[i] + ( (M_sss[i]-physical::si)*physical::rhoi*del_vi + M_sss[i]*(del_vs*physical::rhos + (emp-Fdw)*time_step) )
            / ( M_mld[i]*physical::rhow - del_vi*physical::rhoi - ( del_vs*physical::rhos + (emp-Fdw)*time_step) );

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
            deltaT = std::max(0., physical::mu*M_sss[i] - M_tsurf[i] )
                / ( 1. + physical::ki*M_snow_thick[i]/(physical::ks*M_thick[i]) );
            M_time_relaxation_damage[i] = std::max(time_relaxation_damage*deltaT_relaxation_damage/deltaT, time_step);
        } else {
            M_time_relaxation_damage[i] = time_relaxation_damage;
        }
        // -------------------------------------------------

    }// end for loop
}// end thermo function

// This is the main working function, called from main.cpp (same as perform_simul in the old code)
void
FiniteElement::run()
{
    // Initialise grid and forcing
    this->init();

    // Initialise time
    int ind;
    int pcpt = 0;
    int niter = 0;
    current_time = time_init /*+ pcpt*time_step/(24*3600.0)*/;

    std::cout<<"TIMESTEP= "<< time_step <<"\n";
    std::cout<<"DURATION= "<< duration <<"\n";
#if 0
    gregorian::date epoch = date_time::parse_date<gregorian::date>(
                                                                   vm["simul.time_init"].as<std::string>(),
                                                                   //date_time::ymd_order_dmy
                                                                   date_time::ymd_order_iso
                                                                   );

    std::string time_init_ym = to_iso_string(epoch).substr(0,6);
    std::string init_topaz_file = (boost::format( "TP4DAILY_%1%_3m.nc" ) % time_init_ym ).str();
    std::cout<<"INIT_TOPAZ_FILE "<< init_topaz_file <<"\n";

    std::string init_mit_file = (boost::format( "MITgcm_%1%_3m.nc" ) % time_init_ym ).str();
    std::cout<<"INIT_MIT_FILE "<< init_mit_file <<"\n";

    std::cout<<"INIT TIME= "<< to_iso_string(epoch) <<"\n";
#endif

    double displacement_factor = 1.;
    double minang = 0.;
    bool is_running = true;

    // Check the minimum angle of the grid
    minang = this->minAngle(M_mesh);
    if (minang < vm["simul.regrid_angle"].as<double>())
    {
        std::cout<<"invalid regridding angle: should be smaller than the minimal angle in the intial grid\n";
        throw std::logic_error("invalid regridding angle: should be smaller than the minimal angle in the intial grid");
    }

    // main loop for nextsim program
    while (is_running)
    {
        is_running = ((pcpt+1)*time_step) < duration;

        // if (pcpt > 35)
        //     is_running = false;

        current_time = time_init + pcpt*time_step/(24*3600.0);
        //std::cout<<"TIME STEP "<< pcpt << " for "<< current_time <<"\n";
        std::cout<<"TIME STEP "<< pcpt << " for "<< current_time << " + "<< pcpt*time_step/(24*3600.0) <<"\n";

        // step 0: preparation
        // remeshing and remapping of the prognostic variables

        M_regrid = false;

        if (vm["simul.regrid"].as<std::string>() == "bamg")
        {
            minang = this->minAngle(M_mesh,M_UM,displacement_factor);
            //std::cout<<"REGRID ANGLE= "<< minang <<"\n";

            if ((minang < vm["simul.regrid_angle"].as<double>()) || (pcpt ==0) )
            {
                M_regrid = true;
                chrono.restart();
                std::cout<<"Regriding starts\n";
                this->regrid(pcpt);
                std::cout<<"Regriding done in "<< chrono.elapsed() <<"s\n";
            }
        }

        if (pcpt == 0)
            this->initSimulation();

        if ((pcpt==0) || (M_regrid))
        {
            std::cout<<"bathymetry starts\n";
            this->bathymetry();
            std::cout<<"bathymetry done in "<< chrono.elapsed() <<"s\n";
            chrono.restart();
            std::cout<<"tensors starts\n";
            this->tensors();
            std::cout<<"tensors done in "<< chrono.elapsed() <<"s\n";
            chrono.restart();
            std::cout<<"cohesion starts\n";
            this->cohesion();
            std::cout<<"cohesion done in "<< chrono.elapsed() <<"s\n";
            chrono.restart();
            std::cout<<"Coriolis starts\n";
            this->coriolis();
            std::cout<<"Coriolis done in "<< chrono.elapsed() <<"s\n";
        }

        this->timeInterpolation(pcpt);

        chrono.restart();
        std::cout<<"forcingAtmosphere starts\n";
        this->forcingAtmosphere(M_regrid);
        std::cout<<"forcingAtmosphere done in "<< chrono.elapsed() <<"s\n";

        chrono.restart();
        //std::cout<<"forcingOcean starts\n";
        this->forcingOcean(M_regrid);
        //std::cout<<"forcingOcean done in "<< chrono.elapsed() <<"s\n";

#if 1
        if (pcpt == 0)
        {
            chrono.restart();
            std::cout<<"first export starts\n";
            this->exportResults(0);
            std::cout<<"first export done in " << chrono.elapsed() <<"s\n";
            ind=1;
        }
#endif

        //======================================================================
        // Do the thermodynamics
        //======================================================================
        chrono.restart();
        std::cout<<"thermo starts\n";
        this->thermo();
        std::cout<<"thermo done in "<< chrono.elapsed() <<"s\n";

        //======================================================================
        // Assemble the matrix
        //======================================================================

        this->assemble(pcpt);

        //======================================================================
        // Solve the linear problem
        //======================================================================

        this->solve();

        chrono.restart();
        std::cout<<"updateVelocity starts\n";
        this->updateVelocity();
        std::cout<<"updateVelocity done in "<< chrono.elapsed() <<"s\n";

        chrono.restart();
        std::cout<<"update starts\n";
        this->update();
        std::cout<<"update done in "<< chrono.elapsed() <<"s\n";

#if 1

    if(fmod((pcpt+1)*time_step,output_time_step) == 0)
    {
        chrono.restart();
        std::cout<<"export starts\n";
        this->exportResults(ind);
        std::cout<<"export done in " << chrono.elapsed() <<"s\n";
        ind+=1;
    }

#endif
        ++pcpt;
    }

    this->exportResults(1000);
    std::cout<<"TIMER total = " << chrono_tot.elapsed() <<"s\n";
}

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
FiniteElement::forcingAtmosphere(bool reload)//(double const& u, double const& v)
{
    switch (M_atmosphere_type)
    {
        case setup::AtmosphereType::CONSTANT:
            this->constantAtmosphere();
            break;
        case setup::AtmosphereType::ASR:
            this->asrAtmosphere(reload);
            break;

        default:
            std::cout << "invalid wind forcing"<<"\n";
            throw std::logic_error("invalid wind forcing");
    }
}

void
FiniteElement::constantAtmosphere()
{
    for (int i=0; i<M_num_nodes; ++i)
    {
        M_wind[i]             = Vair_coef*vm["simul.constant_wind_u"].as<double>();
        M_wind[i+M_num_nodes] = Vair_coef*vm["simul.constant_wind_v"].as<double>();
    }

    M_tair.assign(M_num_elements,vm["simul.constant_tair"].as<double>());
    std::cout << "simul.constant_tair:   " << vm["simul.constant_tair"].as<double>() << "\n";

    M_dair.assign(M_num_elements,vm["simul.constant_dair"].as<double>());
    std::cout << "simul.constant_dair:   " << vm["simul.constant_dair"].as<double>() << "\n";

    M_mixrat.assign(M_num_elements,vm["simul.constant_mixrat"].as<double>());
    std::cout << "simul.constant_mixrat: " << vm["simul.constant_mixrat"].as<double>() << "\n";

    M_mslp.assign(M_num_elements,vm["simul.constant_mslp"].as<double>());
    std::cout << "simul.constant_mslp:   " << vm["simul.constant_mslp"].as<double>() << "\n";

    M_Qsw_in.assign(M_num_elements,vm["simul.constant_Qsw_in"].as<double>());
    std::cout << "simul.constant_Qsw_in: " << vm["simul.constant_Qsw_in"].as<double>() << "\n";

    M_Qlw_in.assign(M_num_elements,vm["simul.constant_Qlw_in"].as<double>());
    std::cout << "simul.constant_Qlw_in: " << vm["simul.constant_Qlw_in"].as<double>() << "\n";

    M_precip.assign(M_num_elements,vm["simul.constant_precip"].as<double>());
    std::cout << "simul.constant_precip: " << vm["simul.constant_precip"].as<double>() << "\n";

    M_snowfr.assign(M_num_elements,vm["simul.constant_snowfr"].as<double>());
    std::cout << "simul.constant_snowfr: " << vm["simul.constant_snowfr"].as<double>() << "\n";
}

void
FiniteElement::asrAtmosphere(bool reload)
{

    if ((current_time < M_asr_elements_dataset.ftime_range[0]) || (M_asr_elements_dataset.ftime_range[1] < current_time) || (current_time == time_init) || reload)
    {
        if (current_time == time_init)
            std::cout<<"load forcing from ASR for initial time\n";
        else
            std::cout<<"forcing not available for the current date: load data from ASR\n";

        this->loadDataset(&M_asr_elements_dataset);
        this->loadDataset(&M_asr_nodes_dataset);

        //std::cout<<"forcing not available for the current date\n";
        //throw std::logic_error("forcing not available for the current date");
    }

    double fdt = std::abs(M_asr_elements_dataset.ftime_range[1]-M_asr_elements_dataset.ftime_range[0]);
    std::vector<double> fcoeff(2);
    fcoeff[0] = std::abs(current_time-M_asr_elements_dataset.ftime_range[1])/fdt;
    fcoeff[1] = std::abs(current_time-M_asr_elements_dataset.ftime_range[0])/fdt;

    std::cout<<"LINEAR COEFF 1= "<< fcoeff[0] <<"\n";
    std::cout<<"LINEAR COEFF 2= "<< fcoeff[1] <<"\n";

    double angle_stereo_mesh = -45;
    double angle_stereo_ASR = -175;
    double diff_angle = -(angle_stereo_mesh-angle_stereo_ASR)*PI/180.;

    double cos_m_diff_angle=std::cos(-diff_angle);
    double sin_m_diff_angle=std::sin(-diff_angle);

    double u10_tmp[2];
    double v10_tmp[2];

    std::cout<<"lbefore uv\n";
    for (int i=0; i<M_num_nodes; ++i)
    {
        for(int j=0; j<2; ++j)
        {
            u10_tmp[j]= cos_m_diff_angle*M_uair2[j][i] + sin_m_diff_angle*M_vair2[j][i];
            v10_tmp[j]=-sin_m_diff_angle*M_uair2[j][i] + cos_m_diff_angle*M_vair2[j][i];
        }

        M_wind[i            ] = Vair_coef*(fcoeff[0]*u10_tmp[0] + fcoeff[1]*u10_tmp[1]);
        M_wind[i+M_num_nodes] = Vair_coef*(fcoeff[0]*v10_tmp[0] + fcoeff[1]*v10_tmp[1]);

        // if (i<20)
        //     std::cout<<"data_out["<< i << "]= "<< M_wind[i] << " and "<< M_wind[i+M_num_nodes] <<"\n";
    }

    std::cout<<"lbefore thermo\n";

    for (int i=0; i<M_num_elements; ++i)
    {
        M_tair[i] = fcoeff[0]*M_tair2[0][i] + fcoeff[1]*M_tair2[1][i];
        M_mixrat[i] = fcoeff[0]*M_mixrat2[0][i] + fcoeff[1]*M_mixrat2[1][i];
        M_mslp[i] = fcoeff[0]*M_mslp2[0][i] + fcoeff[1]*M_mslp2[1][i];
        M_Qsw_in[i] = fcoeff[0]*M_Qsw_in2[0][i] + fcoeff[1]*M_Qsw_in2[1][i];
        M_Qlw_in[i] = fcoeff[0]*M_Qlw_in2[0][i] + fcoeff[1]*M_Qlw_in2[1][i];
        M_precip[i] = fcoeff[0]*M_precip2[0][i] + fcoeff[1]*M_precip2[1][i];
        M_snowfr[i] = fcoeff[0]*M_snowfr2[0][i] + fcoeff[1]*M_snowfr2[1][i];
    }

    M_dair.assign(M_num_elements,vm["simul.constant_dair"].as<double>());
    std::cout << "simul.constant_dair:   " << vm["simul.constant_dair"].as<double>() << "\n";
}

void
FiniteElement::forcingOcean(bool reload)//(double const& u, double const& v)
{
    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
            this->constantOcean();
            break;
        case setup::OceanType::TOPAZR:
            this->topazOcean(reload);
            break;


        default:
            std::cout << "invalid ocean forcing"<<"\n";
            throw std::logic_error("invalid ocean forcing");
    }
}

void
FiniteElement::constantOcean()
{
    for (int i=0; i<M_num_nodes; ++i)
    {
        M_ocean[i] = Voce_coef*vm["simul.constant_ocean_v"].as<double>();
        M_ocean[i+M_num_nodes] = Voce_coef*vm["simul.constant_ocean_v"].as<double>();
    }

    M_mld.assign(M_num_elements,vm["simul.constant_mld"].as<double>());
}

void
FiniteElement::topazOcean(bool reload)
{
    if ((current_time < M_topaz_nodes_dataset.ftime_range[0]) || (M_topaz_nodes_dataset.ftime_range[1] < current_time) || (current_time == time_init) || reload)
    {
        if (current_time == time_init)
            std::cout<<"load forcing from TOPAZ for initial time\n";
        else
            std::cout<<"forcing not available for the current date: load data from TOPAZ\n";

        this->loadDataset(&M_topaz_nodes_dataset);
        this->loadDataset(&M_topaz_elements_dataset);

        //std::cout<<"forcing not available for the current date\n";
        //throw std::logic_error("forcing not available for the current date");
    }

    double fdt = std::abs(M_topaz_nodes_dataset.ftime_range[1]-M_topaz_nodes_dataset.ftime_range[0]);
    std::vector<double> fcoeff(2);
    fcoeff[0] = std::abs(current_time-M_topaz_nodes_dataset.ftime_range[1])/fdt;
    fcoeff[1] = std::abs(current_time-M_topaz_nodes_dataset.ftime_range[0])/fdt;

    // std::cout<<"TOPAZ LINEAR COEFF 1= "<< fcoeff[0] <<"\n";
    // std::cout<<"TOPAZ LINEAR COEFF 2= "<< fcoeff[1] <<"\n";

    for (int i=0; i<M_num_nodes; ++i)
    {
        M_ocean[i] = Voce_coef*(fcoeff[0]*M_uoce2[0][i] + fcoeff[1]*M_uoce2[1][i]);
        M_ocean[i+M_num_nodes] = Voce_coef*(fcoeff[0]*M_voce2[0][i] + fcoeff[1]*M_voce2[1][i]);

        M_ssh[i] = ssh_coef*(fcoeff[0]*M_ssh2[0][i] + fcoeff[1]*M_ssh2[1][i]);

        // if (i<20)
        //     std::cout<<"data_out["<< i << "]= "<< M_wind[i] << " and "<< M_wind[i+M_num_nodes] <<"\n";
    }

    for (int i=0; i<M_num_elements; ++i)
    {
        M_ocean_salt[i] = fcoeff[0]*M_ocean_salt2[0][i] + fcoeff[1]*M_ocean_salt2[1][i];
        M_ocean_temp[i] = fcoeff[0]*M_ocean_temp2[0][i] + fcoeff[1]*M_ocean_temp2[1][i];
        M_mld[i] = fcoeff[0]*M_mld2[0][i] + fcoeff[1]*M_mld2[1][i];
    }
    M_mld.assign(M_num_elements,vm["simul.constant_mld"].as<double>());

}

void
FiniteElement::loadDataset(Dataset *dataset)//(double const& u, double const& v)
{

    std::string current_timestr = to_date_string_ym(current_time);
    std::cout<<"TIMESTR= "<< current_timestr <<"\n";

    std::string filename = (boost::format( "%1%/%2%/%3%%4%%5%" )
                                % Environment::nextsimDir().string()
                                % dataset->dirname
                                % dataset->prefix
                                % current_timestr
                                % dataset->postfix
                                ).str();

    std::cout<<"FILE= "<< filename <<"\n";

    double file_dt = 1./dataset->nb_timestep_day;
    double time_start = std::floor(current_time*dataset->nb_timestep_day)/dataset->nb_timestep_day;
    double time_end = std::ceil(current_time*dataset->nb_timestep_day)/dataset->nb_timestep_day;

    // We always need at least two time steps to interpolate between
    if (time_end == time_start)
    {
        time_end = time_start + (1./dataset->nb_timestep_day);
    }

    dataset->ftime_range.resize(0);
    for (double dt=time_start; dt<=time_end; dt+=file_dt)
    {
        dataset->ftime_range.push_back(dt);
    }

    for (int i=0; i<dataset->ftime_range.size(); ++i)
    {
        std::cout<<"TIMEVEC["<< i <<"]= "<< dataset->ftime_range[i] <<"\n";
    }

    // if ((current_time < time_start) || (time_end < current_time))
    // {
    //     std::cout<<"forcing not available for the current date\n";
    //     throw std::logic_error("forcing not available for the current date");
    // }

    int nb_forcing_step = dataset->ftime_range.size();
    std::cout<<"NB_FORCING_STEP= "<< nb_forcing_step <<"\n";

    std::vector<size_t> index_start(4);
    std::vector<size_t> index_end(4);

    std::vector<double> XTIME(dataset->time.dimensions[0].end);

    netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

    netCDF::NcVar VTIME = dataFile.getVar(dataset->time.name);

    for(int j=0; j<dataset->variables.size(); ++j)
        dataset->variables[j].NcVar = dataFile.getVar(dataset->variables[j].name);

    VTIME.getVar(&XTIME[0]);


    std::cout<<"VALUE= "<< from_date_string(dataset->reference_date) <<"\n";
    std::for_each(XTIME.begin(), XTIME.end(), [&](double& f){ f = f/24.0+from_date_string(dataset->reference_date); });

    for (int i=0; i<dataset->ftime_range.size(); ++i)
    {
        std::cout<<"---TIMEVEC["<< i <<"]= "<< dataset->ftime_range[i] <<" : current_time= "<< current_time <<"\n";
    }

    std::vector<double> tmp_interpolated_field(dataset->target_size);

    int N_data =dataset->variables.size();
    int M  =dataset->grid->dimension_y.end;
    int N  = dataset->grid->dimension_x.end;
    int MN = M*N;

    double* data_in = new double[N_data*nb_forcing_step*MN];

    std::vector<double> data_in_tmp(MN);

    netCDF::NcVarAtt att;
    double scale_factor;
    double add_offset;

    for (int fstep=0; fstep < nb_forcing_step; ++fstep)
    {
        double ftime = dataset->ftime_range[fstep];

        if (to_date_string_ym(std::floor(ftime)) != to_date_string_ym(current_time))
        {
            std::string f_timestr = to_date_string_ym(std::floor(ftime));
            std::cout<<"F_TIMESTR= "<< f_timestr <<"\n";

            filename = (boost::format( "%1%/%2%/%3%%4%%5%" )
                                % Environment::nextsimDir().string()
                                % dataset->dirname
                                % dataset->prefix
                                % f_timestr
                                % dataset->postfix
                                ).str();

            std::cout<<"FILENAME= "<< filename <<"\n";

            netCDF::NcFile fdataFile(filename, netCDF::NcFile::read);
            netCDF::NcVar FVTIME = dataFile.getVar(dataset->time.name);
            FVTIME.getVar(&XTIME[0]);
            std::for_each(XTIME.begin(), XTIME.end(), [&](double& f){ f = f/24.0+from_date_string(dataset->reference_date); });

            // for (int i=0; i<31; ++i)
            // {
            //     std::cout<<"TIME["<< i <<"]= "<< XTIME[i] <<"\n";
            // }
        }

        auto it = std::find(XTIME.begin(), XTIME.end(), ftime);
        int index = std::distance(XTIME.begin(),it);
        std::cout<<"FIND "<< ftime <<" in index "<< index <<"\n";

        for(int j=0; j<dataset->variables.size(); ++j)
        {
            index_start.resize(dataset->variables[j].dimensions.size());
            index_end.resize(dataset->variables[j].dimensions.size());

            for(int k=0; k<dataset->variables[j].dimensions.size(); ++k)
            {
                index_start[k] = dataset->variables[j].dimensions[k].start;
                index_end[k] = dataset->variables[j].dimensions[k].end;
            }
            index_start[0] = index;
            index_end[0] = 1;

            dataset->variables[j].NcVar.getVar(index_start,index_end,&data_in_tmp[0]);

            // Need to multiply with scale factor and add offset - these are stored as variable attributes
            scale_factor=1.;
            try
            {
                att = dataset->variables[j].NcVar.getAtt("scale_factor");
                att.getValues(&scale_factor);
            }
            catch(netCDF::exceptions::NcException& e)
            {}

            add_offset=0.;
            try
            {
                att = dataset->variables[j].NcVar.getAtt("add_offset");
                att.getValues(&add_offset);
            }
            catch(netCDF::exceptions::NcException& e)
            {}

            for (int i=0; i<(MN); ++i)
                data_in[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j]=data_in_tmp[i]*scale_factor + add_offset;
        }
    }

    double* data_out;
    double tmp_data;

    auto RX = M_mesh.coordX(dataset->grid->rotation_angle);
    auto RY = M_mesh.coordY(dataset->grid->rotation_angle);

    switch(dataset->case_number)
    {
        case 1: case 3:
            RX = M_mesh.bcoordX(dataset->grid->rotation_angle);
            RY = M_mesh.bcoordY(dataset->grid->rotation_angle);
            break;
    }

    //int interp_type = TriangleInterpEnum;
    int interp_type = BilinearInterpEnum;
    //int interp_type = NearestInterpEnum;

    std::cout<<"before interp " <<"\n";

    switch(dataset->grid->interpolation_method)
    {
        case setup::InterpolationType::InterpFromGridToMesh:
            InterpFromGridToMeshx(  data_out, &dataset->grid->gridX[0], dataset->grid->gridX.size(), &dataset->grid->gridY[0], dataset->grid->gridY.size(), 
                                  &data_in[0], dataset->grid->gridX.size(), dataset->grid->gridY.size(), 
                                  dataset->variables.size()*nb_forcing_step,
                                 &RX[0], &RY[0], dataset->target_size, 1.0, interp_type);
        break;
        case setup::InterpolationType::InterpFromMeshToMesh2dx:
            InterpFromMeshToMesh2dx(&data_out,
                                        dataset->grid->pfindex,&dataset->grid->gridX[0],&dataset->grid->gridY[0],
                                        dataset->grid->gridX.size(),dataset->grid->pfnels,
                                        &data_in[0],
                                        dataset->grid->gridX.size(),N_data*nb_forcing_step,
                                        &RX[0], &RY[0], dataset->target_size,
                                        false /*options*/);
        break;
        default:
            std::cout << "invalid interpolation type:" <<"\n";
            throw std::logic_error("invalid interpolation type");
    }

std::cout<<"after interp " <<"\n";

    for (int fstep=0; fstep < nb_forcing_step; ++fstep)
    {
        for(int j=0; j<dataset->variables.size(); ++j)
        {
            for (int i=0; i<dataset->target_size; ++i)
            {
                tmp_data=data_out[(dataset->variables.size()*nb_forcing_step)*i+fstep*dataset->variables.size()+j];
                tmp_interpolated_field[i]=dataset->variables[j].a*tmp_data+dataset->variables[j].b;
            }

            switch(dataset->case_number)
            {
                case 0:
                    switch(j)
                    {
                    case 0:
                        M_uair2[fstep]=tmp_interpolated_field;
                        break;
                    case 1:
                        M_vair2[fstep]=tmp_interpolated_field;
                        break;
                    }
                break;

                case 1:
                    switch(j)
                    {
                    case 0:
                        M_tair2[fstep]=tmp_interpolated_field;
                        break;
                    case 1:
                        M_mixrat2[fstep]=tmp_interpolated_field;
                        break;
                    case 2:
                        M_mslp2[fstep]=tmp_interpolated_field;
                        break;
                    case 3:
                        M_Qsw_in2[fstep]=tmp_interpolated_field;
                        break;
                    case 4:
                        M_Qlw_in2[fstep]=tmp_interpolated_field;
                        break;
                    case 5:
                        M_snowfr2[fstep]=tmp_interpolated_field;
                        break;
                    case 6:
                        M_precip2[fstep]=tmp_interpolated_field;
                        break;
                    default:
                        break;
                    }
                break;

                case 2:
                    switch(j)
                    {
                    case 0:
                        M_uoce2[fstep]=tmp_interpolated_field;
                        break;
                    case 1:
                        M_voce2[fstep]=tmp_interpolated_field;
                        break;
                    case 2:
                        M_ssh2[fstep]=tmp_interpolated_field;
                        break;
                    }
                break;

                case 3:
                    switch(j)
                    {
                    case 0:
                        M_ocean_temp2[fstep]=tmp_interpolated_field;
                        break;
                    case 1:
                        M_ocean_salt2[fstep]=tmp_interpolated_field;
                        break;
                    case 2:
                        M_mld2[fstep]=tmp_interpolated_field;
                        break;
                    }
                break;
            }
        }
    }

    std::cout<<"end load" <<"\n";
}

void
FiniteElement::loadGrid(Grid *grid)
{
    std::string current_timestr = to_date_string_ym(current_time);
    std::cout<<"TIMESTR= "<< current_timestr <<"\n";
    std::string filename = (boost::format( "%1%/%2%/%3%%4%%5%" )
                                % Environment::nextsimDir().string()
                                % grid->dirname
                                % grid->prefix
                                % current_timestr 
                                % grid->postfix
                                ).str();

    // read in re-analysis coordinates
    std::vector<size_t> index_x_end(2);
    std::vector<size_t> index_y_end(2);

    std::vector<size_t> index_x_start(2);
    std::vector<size_t> index_y_start(2);

    index_y_start[0] = 0;
    index_y_start[1] = 0;

    index_y_end[0] = grid->dimension_y.end;
    index_y_end[1] = grid->dimension_x.end;

    index_x_start[0] = 0;
    index_x_start[1] = 0;

    index_x_end[0] = grid->dimension_y.end;
    index_x_end[1] = grid->dimension_x.end;

    if(grid->interpolation_method==setup::InterpolationType::InterpFromGridToMesh)
    {
        index_y_end[1] = 1;
        index_x_end[0] = 1;
    }

    std::vector<double> XLAT(index_x_end[0]*index_x_end[1]);
    std::vector<double> XLON(index_x_end[0]*index_x_end[1]);
    std::vector<double> YLAT(index_y_end[0]*index_y_end[1]);
    std::vector<double> YLON(index_y_end[0]*index_y_end[1]);

    std::cout<<"GRID : READ NETCDF starts\n";
    netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
    netCDF::NcVar VLAT = dataFile.getVar(grid->latitude.name);
    netCDF::NcVar VLON = dataFile.getVar(grid->longitude.name);
    std::cout<<"GRID : READ NETCDF done\n";

    VLAT.getVar(index_x_start,index_x_end,&XLAT[0]);
    VLON.getVar(index_x_start,index_x_end,&XLON[0]);

    VLAT.getVar(index_y_start,index_y_end,&YLAT[0]);
    VLON.getVar(index_y_start,index_y_end,&YLON[0]);

    grid->gridX.resize(index_x_end[0]*index_x_end[1]);
    grid->gridY.resize(index_y_end[0]*index_y_end[1]);

    double RE = 6378.273;
    mapx_class *map;
    std::string configfile = (boost::format( "%1%/%2%/%3%" )
                                % Environment::nextsimDir().string()
                                % grid->dirname
                                % grid->mpp_file
                                ).str();

    std::vector<char> str(configfile.begin(), configfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    for (int i=0; i<index_x_end[0]; ++i)
    {
        for (int j=0; j<index_x_end[1]; ++j)
        {
            grid->gridX[index_x_end[1]*i+j]=latLon2XY(XLAT[index_x_end[1]*i+j], XLON[index_x_end[1]*i+j], map, configfile)[0];
        }
    }

    for (int i=0; i<index_y_end[0]; ++i)
    {
        for (int j=0; j<index_y_end[1]; ++j)
        {
            grid->gridY[index_y_end[1]*i+j]=latLon2XY(YLAT[index_y_end[1]*i+j], YLON[index_y_end[1]*i+j], map, configfile)[1];
        }
    }

    close_mapx(map);

    if(grid->interpolation_method==setup::InterpolationType::InterpFromMeshToMesh2dx)
    {
        std::cout<<"GRID : Triangulate starts\n";
        BamgTriangulatex(&grid->pfindex,&grid->pfnels,&grid->gridX[0],&grid->gridY[0],grid->gridX.size());
        std::cout<<"GRID : NUMTRIANGLES= "<< grid->pfnels <<"\n";
        std::cout<<"GRID : Triangulate done\n";
    }
}

void
FiniteElement::initSlabOcean()
{
    switch (M_ocean_type)
    {
        case setup::OceanType::CONSTANT:
            std::fill(M_sst.begin(), M_sst.end(), -1.8);
            std::fill(M_sss.begin(), M_sss.end(), -1.8/physical::mu);
            break;
        case setup::OceanType::TOPAZR:
            this->topazOcean(1); // This is lazy re-use of code
            for ( int i=0; i<M_num_elements; ++i)
            {
                // Make sure the erroneous salinity and temperature don't screw up the initialisation too badly
                // This can still be done much better!
                M_sss[i] = std::max(physical::si, M_ocean_salt[i]);
                M_sst[i] = std::max(M_sss[i]*physical::mu, M_ocean_temp[i]);
            }
            break;



        default:
            std::cout << "invalid ocean initialisation"<<"\n";
            throw std::logic_error("invalid ocean forcing");
    }
}

void
FiniteElement::initConcentration()
{
    switch (M_conc_type)
    {
        case setup::ConcentrationType::CONSTANT:
            this->constantConc();
            break;
        case setup::ConcentrationType::TOPAZ4:
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
    std::fill(M_conc.begin(), M_conc.end(), vm["simul.init_concentration"].as<double>());

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
    std::string init_timestr = to_date_string_ym(time_init);
    std::cout<<"TIMEINITSTR= "<< init_timestr <<"\n";
    std::string init_filename = (boost::format( "%1%/data/TP4DAILY_%2%_3m.nc" )
                                       % Environment::nextsimDir().string()
                                       % init_timestr ).str();

    std::cout<<"TOPAZ INIT FILE= "<< init_filename <<"\n";

    // read in latitude and longitude
    std::vector<size_t> index_start(2);
    std::vector<size_t> index_end(2);

    index_start[0] = 0;
    index_start[1] = 0;

    index_end[0] = 1101;
    index_end[1] = 761;

    std::vector<double> LAT(index_end[0]*index_end[1]);
    std::vector<double> LON(index_end[0]*index_end[1]);

    std::vector<double> XTIME(31);
    std::vector<size_t> index_fhice_start(3,0);
    std::vector<size_t> index_fhice_end(3);

    std::cout<<"READING NETCDF file "<< init_filename <<" starts\n";
    //netCDF::NcFile::FileFormat format = netCDF::NcFile::classic;
    netCDF::NcFile dataFile(init_filename, netCDF::NcFile::read);
    netCDF::NcVar VLAT = dataFile.getVar("latitude");
    netCDF::NcVar VLON = dataFile.getVar("longitude");
    netCDF::NcVar VFICE;
    netCDF::NcVar VHICE;
    netCDF::NcVar VSNOW;

    if (M_conc_type == setup::ConcentrationType::TOPAZ4)
    {
        VFICE = dataFile.getVar("fice");
    }

    if (M_thick_type == setup::ThicknessType::TOPAZ4)
    {
        VHICE = dataFile.getVar("hice");
    }

    if (M_snow_thick_type == setup::SnowThicknessType::TOPAZ4)
    {
        VSNOW = dataFile.getVar("hsnow");
    }

    netCDF::NcVar VTIME = dataFile.getVar("time");
    std::cout<<"READING NETCDF "<< init_filename << " done\n";

    VLAT.getVar(index_start,index_end,&LAT[0]);
    VLON.getVar(index_start,index_end,&LON[0]);

    VTIME.getVar(&XTIME[0]);

    std::vector<double> X(index_end[0]*index_end[1]);
    std::vector<double> Y(index_end[0]*index_end[1]);

    mapx_class *map;
    std::string configfile = Environment::nextsimDir().string() + "/data/NpsNextsim.mpp";
    std::vector<char> str(configfile.begin(), configfile.end());
    str.push_back('\0');
    map = init_mapx(&str[0]);

    std::vector<double> xy(2);

    for (int i=0; i<index_end[0]; ++i)
    {
        for (int j=0; j<index_end[1]; ++j)
        {
            xy=latLon2XY(LAT[index_end[1]*i+j], LON[index_end[1]*i+j], map, configfile);
            X[index_end[1]*i+j] = xy[0];
            Y[index_end[1]*i+j] = xy[1];
        }
    }

    close_mapx(map);

    auto RX = M_mesh.bcoordX();
    auto RY = M_mesh.bcoordY();

#if 0
    std::cout<<"MIN BOUND TOPAZX= "<< *std::min_element(X.begin(),X.end()) <<"\n";
    std::cout<<"MAX BOUND TOPAZX= "<< *std::max_element(X.begin(),X.end()) <<"\n";

    std::cout<<"MIN BOUND TOPAZY= "<< *std::min_element(Y.begin(),Y.end()) <<"\n";
    std::cout<<"MAX BOUND TOPAZY= "<< *std::max_element(Y.begin(),Y.end()) <<"\n";

    std::cout<<"MIN BOUND MESHX= "<< *std::min_element(RX.begin(),RX.end()) <<"\n";
    std::cout<<"MAX BOUND MESHX= "<< *std::max_element(RX.begin(),RX.end()) <<"\n";

    std::cout<<"MIN BOUND MESHY= "<< *std::min_element(RY.begin(),RY.end()) <<"\n";
    std::cout<<"MAX BOUND MESHY= "<< *std::max_element(RY.begin(),RY.end()) <<"\n";
#endif

    std::cout<<"VALUE= "<< from_date_string("1950-01-01") <<"\n";
    double target = (time_init - from_date_string("1950-01-01"))*24.0;
    //std::for_each(XTIME.begin(), XTIME.end(), [&](double& f){ f = f/24.0+from_date_string("1950-01-01"); });
    std::cout<<"TARGET= "<< target <<"\n";

    // for (int i=0; i<31; ++i)
    // {
    //     std::cout<<"TIME["<< i <<"]= "<< XTIME[i] <<"\n";
    // }

    if (std::find(XTIME.begin(), XTIME.end(), target) == XTIME.end())
    {
        std::cout<<"forcing not available for this initial time\n";
        std::cout<<"take the largest integer value not greater than initial time\n";

        target = (std::floor(time_init) - from_date_string("1950-01-01"))*24.0;
    }

    auto it = std::find(XTIME.begin(), XTIME.end(), target);
    int index = std::distance(XTIME.begin(),it);
    std::cout<<"INIT TIME TOPAZ FOUND "<< target <<" in index "<< index <<"\n";

    index_fhice_start[0] = index;
    index_fhice_start[1] = index_start[0];
    index_fhice_start[2] = index_start[1];

    index_fhice_end[0] = 1;
    index_fhice_end[1] = index_end[0];
    index_fhice_end[2] = index_end[1];

 #if 0
    std::cout<<"NETCDF INFO: "<<dataFile.getVarCount()<<" variables\n";
    std::cout<<"NETCDF INFO: "<<dataFile.getAttCount()<<" attributes\n";
    std::cout<<"NETCDF INFO: "<<dataFile.getDimCount()<<" dimensions\n";
    std::cout<<"NETCDF INFO: "<<dataFile.getGroupCount()<<" groups\n";
    std::cout<<"NETCDF INFO: "<<dataFile.getTypeCount()<<" types\n";
#endif

    // void* data_void[1101*index_end[1]];
    // VFICE.getVar(index_fhice_start,index_fhice_end,data_void);

    //float* data_values = (float*)data_void;
    //float* data_values = reinterpret_cast<float*>(data_void);
    // double* data_values = (double*)data_void;

    std::vector<double> data_in_fice;
    std::vector<double> data_in_hice;
    std::vector<double> data_in_snow;

    std::vector<double> reduced_data_in_fice;
    std::vector<double> reduced_FX;
    std::vector<double> reduced_FY;

    std::vector<double> reduced_data_in_hice;
    std::vector<double> reduced_HX;
    std::vector<double> reduced_HY;

    std::vector<double> reduced_data_in_snow;
    std::vector<double> reduced_SX;
    std::vector<double> reduced_SY;

    netCDF::NcVarAtt att;
    double scale_factor_fice, add_offset_fice;
    double scale_factor_hice, add_offset_hice;
    double scale_factor_snow, add_offset_snow;
    int FillValue_fice, FillValue_hice, FillValue_snow;

    // Need to multiply with scale factor and add offset - these are stored as variable attributes
    if (M_conc_type == setup::ConcentrationType::TOPAZ4)
    {
        data_in_fice.resize(index_end[0]*index_end[1]);
        VFICE.getVar(index_fhice_start,index_fhice_end,&data_in_fice[0]);

	att = VFICE.getAtt("scale_factor");
	att.getValues(&scale_factor_fice);
	att = VFICE.getAtt("add_offset");
	att.getValues(&add_offset_fice);
	att = VFICE.getAtt("_FillValue");
	att.getValues(&FillValue_fice);
    }

    if (M_thick_type == setup::ThicknessType::TOPAZ4)
    {
        data_in_hice.resize(index_end[0]*index_end[1]);
        VHICE.getVar(index_fhice_start,index_fhice_end,&data_in_hice[0]);

	att = VHICE.getAtt("scale_factor");
	att.getValues(&scale_factor_hice);
	att = VHICE.getAtt("add_offset");
	att.getValues(&add_offset_hice);
	att = VHICE.getAtt("_FillValue");
	att.getValues(&FillValue_hice);
    }

    if (M_snow_thick_type == setup::SnowThicknessType::TOPAZ4)
    {
        data_in_snow.resize(index_end[0]*index_end[1]);
        VSNOW.getVar(index_fhice_start,index_fhice_end,&data_in_snow[0]);

	att = VSNOW.getAtt("scale_factor");
	att.getValues(&scale_factor_snow);
	att = VSNOW.getAtt("add_offset");
	att.getValues(&add_offset_snow);
	att = VSNOW.getAtt("_FillValue");
	att.getValues(&FillValue_snow);
    }

    for (int i=0; i<index_end[0]; ++i)
    {
        for (int j=0; j<index_end[1]; ++j)
        {
            if (M_conc_type == setup::ConcentrationType::TOPAZ4)
            {
                // maskvfh = data_in_fice[index_end[1]*i+j]*scale_factor_fice+add_offset_fice;
                // maskvfh = std::abs(maskvfh);

                if (data_in_fice[index_end[1]*i+j] != FillValue_fice)
                {
                    reduced_data_in_fice.push_back(data_in_fice[index_end[1]*i+j]*scale_factor_fice+add_offset_fice);
                    reduced_FX.push_back(X[index_end[1]*i+j]);
                    reduced_FY.push_back(Y[index_end[1]*i+j]);
                }
            }

            if (M_thick_type == setup::ThicknessType::TOPAZ4)
            {
                // maskvfh = data_in_hice[index_end[1]*i+j]*scale_factor_hice+add_offset_hice;
                // maskvfh = std::abs(maskvfh);

                if (data_in_fice[index_end[1]*i+j] != FillValue_hice)
                {
                    reduced_data_in_hice.push_back(data_in_hice[index_end[1]*i+j]*scale_factor_hice+add_offset_hice);
                    reduced_HX.push_back(X[index_end[1]*i+j]);
                    reduced_HY.push_back(Y[index_end[1]*i+j]);
                }
            }

            if (M_snow_thick_type == setup::SnowThicknessType::TOPAZ4)
            {
                // maskvfh = data_in_snow[index_end[1]*i+j]*scale_factor_snow+add_offset_snow;
                // maskvfh = std::abs(maskvfh);

                if (data_in_fice[index_end[1]*i+j] != FillValue_snow)
                {
                    reduced_data_in_snow.push_back(data_in_snow[index_end[1]*i+j]*scale_factor_snow+add_offset_snow);
                    reduced_SX.push_back(X[index_end[1]*i+j]);
                    reduced_SY.push_back(Y[index_end[1]*i+j]);
                }
            }
        }
    }

#if 0
    if (M_conc_type == setup::ConcentrationType::TOPAZ4)
    {
        std::cout<<"MIN DATA_IN FICE= "<< *std::min_element(data_in_fice.begin(),data_in_fice.end()) <<"\n";
        std::cout<<"MAX DATA_IN FICE= "<< *std::max_element(data_in_fice.begin(),data_in_fice.end()) <<"\n";
    }

    if (M_thick_type == setup::ThicknessType::TOPAZ4)
    {
        std::cout<<"MIN DATA_IN HICE= "<< *std::min_element(data_in_hice.begin(),data_in_hice.end()) <<"\n";
        std::cout<<"MAX DATA_IN HICE= "<< *std::max_element(data_in_hice.begin(),data_in_hice.end()) <<"\n";
    }
#endif

#if 1
    // bamg triangulation

    if (M_conc_type == setup::ConcentrationType::TOPAZ4)
    {
        std::cout<<"SIZE REDUCED_FX= "<< reduced_FX.size() <<"\n";
        std::cout<<"SIZE REDUCED_FY= "<< reduced_FY.size() <<"\n";

        std::cout<<"MIN DATA_RF_IN = "<< *std::min_element(reduced_data_in_fice.begin(),reduced_data_in_fice.end()) <<"\n";
        std::cout<<"MAX DATA_RF_IN = "<< *std::max_element(reduced_data_in_fice.begin(),reduced_data_in_fice.end()) <<"\n";
    }

    if (M_thick_type == setup::ThicknessType::TOPAZ4)
    {
        std::cout<<"SIZE REDUCED_HX= "<< reduced_HX.size() <<"\n";
        std::cout<<"SIZE REDUCED_HY= "<< reduced_HY.size() <<"\n";

        std::cout<<"MIN DATA_RH_IN = "<< *std::min_element(reduced_data_in_hice.begin(),reduced_data_in_hice.end()) <<"\n";
        std::cout<<"MAX DATA_RH_IN = "<< *std::max_element(reduced_data_in_hice.begin(),reduced_data_in_hice.end()) <<"\n";
    }

    if (M_snow_thick_type == setup::SnowThicknessType::TOPAZ4)
    {
        std::cout<<"SIZE REDUCED_SX= "<< reduced_SX.size() <<"\n";
        std::cout<<"SIZE REDUCED_SY= "<< reduced_SY.size() <<"\n";

        std::cout<<"MIN DATA_RS_IN = "<< *std::min_element(reduced_data_in_snow.begin(),reduced_data_in_snow.end()) <<"\n";
        std::cout<<"MAX DATA_RS_IN = "<< *std::max_element(reduced_data_in_snow.begin(),reduced_data_in_snow.end()) <<"\n";
    }

    int* pfindex;
    int pfnels;

    if (M_conc_type == setup::ConcentrationType::TOPAZ4)
    {
        std::cout<<"FICE: Triangulate starts\n";
        BamgTriangulatex(&pfindex,&pfnels,&reduced_FX[0],&reduced_FY[0],reduced_FX.size());
        std::cout<<"FICE: NUMTRIANGLES= "<< pfnels <<"\n";
        std::cout<<"FICE: Triangulate done\n";
        // for (int i=0; i<Y.size(); ++i)
        // {
        //     std::cout<<"Point["<< i <<"]= ("<< RX[i] << " , "<< RY[i] <<")\n";
        // }
    }

    int* phindex;
    int phnels;

    if (M_thick_type == setup::ThicknessType::TOPAZ4)
    {
        std::cout<<"HICE: Triangulate starts\n";
        BamgTriangulatex(&phindex,&phnels,&reduced_HX[0],&reduced_HY[0],reduced_HX.size());
        std::cout<<"HICE: NUMTRIANGLES= "<< phnels <<"\n";
        std::cout<<"HICE: Triangulate done\n";
        // for (int i=0; i<Y.size(); ++i)
        // {
        //     std::cout<<"Point["<< i <<"]= ("<< RX[i] << " , "<< RY[i] <<")\n";
        // }
    }

    int* psindex;
    int psnels;

    if (M_snow_thick_type == setup::SnowThicknessType::TOPAZ4)
    {
        std::cout<<"HSNOW: Triangulate starts\n";
        BamgTriangulatex(&psindex,&psnels,&reduced_SX[0],&reduced_SY[0],reduced_SX.size());
        std::cout<<"HSNOW: NUMTRIANGLES= "<< psnels <<"\n";
        std::cout<<"HSNOW: Triangulate done\n";
        // for (int i=0; i<Y.size(); ++i)
        // {
        //     std::cout<<"Point["<< i <<"]= ("<< RX[i] << " , "<< RY[i] <<")\n";
        // }
    }
#endif

    //int interp_type = TriangleInterpEnum;
    //int interp_type = BilinearInterpEnum;
    int interp_type = NearestInterpEnum;

    // std::vector<double> data_out_fice_tmp;
    // std::vector<double> data_out_hice_tmp;
    // std::vector<double> data_out_snow_tmp;

    if (M_conc_type == setup::ConcentrationType::TOPAZ4)
    {
        double* data_out_fice;
        //data_out_fice_tmp.resize(M_num_elements);

        // InterpFromGridToMeshx(data_out_fice, &X[0], X.size(), &Y[0], Y.size(), &data_in_fice[0], Y.size(), X.size(),
        //                       &RX[0], &RY[0], M_mesh.numTriangles(), 1.0, interp_type);

        InterpFromMeshToMesh2dx(&data_out_fice,
                                pfindex,&reduced_FX[0],&reduced_FY[0],
                                reduced_FX.size(),pfnels,
                                &reduced_data_in_fice[0],
                                reduced_FX.size(),1,
                                &RX[0], &RY[0], M_num_elements,
                                false /*options*/);


        for (int i=0; i<M_num_elements; ++i)
        {
            //data_out_fice_tmp[i] = data_out_fice[i];
            // M_conc[i] = data_out_fice[i];
            M_conc[i] = (data_out_fice[i]>1e-14) ? data_out_fice[i] : 0.;
            //std::cout<<"MCONC["<< i <<"]= "<< M_conc[i] <<"\n";
        }
    }

    if (M_thick_type == setup::ThicknessType::TOPAZ4)
    {
        double* data_out_hice;
        //data_out_hice_tmp.resize(M_num_elements);

        // InterpFromGridToMeshx(data_out_hice, &X[0], X.size(), &Y[0], Y.size(), &data_in_hice[0], Y.size(), X.size(),
        //                       &RX[0], &RY[0], M_mesh.numTriangles(), 1.0, interp_type);

        InterpFromMeshToMesh2dx(&data_out_hice,
                                phindex,&reduced_HX[0],&reduced_HY[0],
                                reduced_HX.size(),phnels,
                                &reduced_data_in_hice[0],
                                reduced_HX.size(),1,
                                &RX[0], &RY[0], M_num_elements,
                                false /*options*/);


        for (int i=0; i<M_num_elements; ++i)
        {
            //data_out_hice_tmp[i] = data_out_hice[i];
            //M_thick[i] = data_out_hice[i];
            M_thick[i] = (data_out_hice[i]>1e-14) ? data_out_hice[i] : 0.;

            //std::cout<<"MTHICKC["<< i <<"]= "<< M_thick[i] <<"\n";
        }
    }

    if (M_snow_thick_type == setup::SnowThicknessType::TOPAZ4)
    {
        double* data_out_snow;
        //data_out_snow_tmp.resize(M_num_elements);

        // InterpFromGridToMeshx(data_out_hice, &X[0], X.size(), &Y[0], Y.size(), &data_in_hice[0], Y.size(), X.size(),
        //                       &RX[0], &RY[0], M_num_elements, 1.0, interp_type);

        InterpFromMeshToMesh2dx(&data_out_snow,
                                psindex,&reduced_SX[0],&reduced_SY[0],
                                reduced_SX.size(),psnels,
                                &reduced_data_in_snow[0],
                                reduced_SX.size(),1,
                                &RX[0], &RY[0], M_num_elements,
                                false /*options*/);


        for (int i=0; i<M_num_elements; ++i)
        {
            //data_out_snow_tmp[i] = data_out_snow[i];
            // M_snow_thick[i] = data_out_snow[i];
            M_snow_thick[i] = (data_out_snow[i]>1e-14) ? data_out_snow[i] : 0.;

            //std::cout<<"MTHICKC["<< i <<"]= "<< M_snow_thick[i] <<"\n";
        }
    }

    if (M_thick_type == setup::ThicknessType::TOPAZ4)
    {
        for (int i=0; i<M_num_elements; ++i)
        {
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
            //std::cout<<"MTHICKC["<< i <<"]= "<< M_thick[i] <<"\n";
        }
    }


#if 0
    if (M_conc_type == forcing::ConcentrationType::TOPAZ4)
    {
        std::cout<<"MIN DATA_OUT FICE= "<< *std::min_element(data_out_fice_tmp.begin(),data_out_fice_tmp.end()) <<"\n";
        std::cout<<"MAX DATA_OUT FICE= "<< *std::max_element(data_out_fice_tmp.begin(),data_out_fice_tmp.end()) <<"\n";
    }

    if (M_thick_type == setup::ThicknessType::TOPAZ4)
    {
        std::cout<<"MIN DATA_OUT HICE= "<< *std::min_element(data_out_hice_tmp.begin(),data_out_hice_tmp.end()) <<"\n";
        std::cout<<"MAX DATA_OUT HICE= "<< *std::max_element(data_out_hice_tmp.begin(),data_out_hice_tmp.end()) <<"\n";
    }

    if (M_snow_thick_type == setup::SnowThicknessType::TOPAZ4)
    {
        std::cout<<"MIN DATA_OUT SNOW= "<< *std::min_element(data_out_snow_tmp.begin(),data_out_snow_tmp.end()) <<"\n";
        std::cout<<"MAX DATA_OUT SNOW= "<< *std::max_element(data_out_snow_tmp.begin(),data_out_snow_tmp.end()) <<"\n";
    }
#endif
}

void
FiniteElement::initThickness()
{
    switch (M_thick_type)
    {
        case setup::ThicknessType::CONSTANT:
            this->constantThick();
            break;
        case setup::ThicknessType::TOPAZ4:
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
        M_thick[i] = (vm["simul.init_thickness"].as<double>())*M_conc[i];
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
        case setup::DamageType::CONSTANT:
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
        case setup::SnowThicknessType::CONSTANT:
            this->constantSnowThick();
            break;
        case setup::SnowThicknessType::TOPAZ4:
            this->topazSnowThick();
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
        M_snow_thick[i] = (vm["simul.init_snow_thickness"].as<double>())*M_conc[i];
    }
}

void
FiniteElement::topazSnowThick()
{
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
        case setup::DrifterType::EQUALLYSPACED:
            this->equallySpacedDrifter();
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
        //std::cout <<"Coeff= "<< M_fcor[i] <<"\n";
    }
}

void
FiniteElement::bathymetry()
{
    // Interpolation of the bathymetry
    if (vm["simul.Lemieux_basal_k2"].as<double>() > 0 )
    {
        double* depth_out;

        InterpFromMeshToMesh2dx(&depth_out,
                                &M_mesh_init.indexTr()[0],&M_mesh_init.coordX()[0],&M_mesh_init.coordY()[0],
                                M_mesh_init.numNodes(),M_mesh_init.numTriangles(),
                                &M_bathy_depth[0],
                                M_mesh_init.numNodes(),1,
                                &M_mesh.bcoordX()[0],&M_mesh.bcoordY()[0],M_num_elements,
                                false);

        M_element_depth.assign(M_num_elements,0.);

        for (int i=0; i<M_num_elements; ++i)
        {
            M_element_depth[i] = depth_out[i];
        }
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
    //std::vector<element_type> mesh_edges;
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
 #if 0
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
#endif

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

        //mesh_triangles.insert(std::make_pair(tr,gmshElt));
        mesh_triangles.push_back(gmshElt);
    }

    std::cout<<"\n";
    std::cout<<"INFO: Previous  NumNodes     = "<< M_mesh.numNodes() <<"\n";
    std::cout<<"INFO: Previous  NumTriangles = "<< M_mesh.numTriangles() <<"\n";
    //std::cout<<"INFO: Previous  NumEdges     = "<< M_mesh.numEdges() <<"\n";

    M_mesh_previous = M_mesh;
    //M_mesh = mesh_type(mesh_nodes,mesh_edges,mesh_triangles);
    M_mesh = mesh_type(mesh_nodes,mesh_triangles);
    //M_mesh.writeTofile("out.msh");

    M_elements = M_mesh.triangles();
    M_nodes = M_mesh.nodes();

    M_num_elements = M_mesh.numTriangles();
    M_num_nodes = M_mesh.numNodes();

    std::cout<<"\n";
    std::cout<<"INFO: Current  NumNodes      = "<< M_mesh.numNodes() <<"\n";
    std::cout<<"INFO: Current  NumTriangles  = "<< M_mesh.numTriangles() <<"\n";
    //std::cout<<"INFO: Current  NumEdges      = "<< M_mesh.numEdges() <<"\n";
    std::cout<<"\n";

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

    std::cout<<"\n";
    std::cout<<"GRAPHCSR INFO: MIN NZ (per row)      = "<< *std::min_element(dz.begin(),dz.end()) <<"\n";
    std::cout<<"GRAPHCSR INFO: MAX NZ (per row)      = "<< *std::max_element(dz.begin(),dz.end()) <<"\n";
    std::cout<<"GRAPHCSR INFO: NNZ (total)           = "<< ddz_j.size() <<"\n";
    std::cout<<"\n";


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
FiniteElement::exportResults(int step, bool export_mesh)
{
#if 0
    vector_type mx;
    mx.init(3*M_num_elements);
    vector_type my;
    my.init(3*M_num_elements);
    vector_type mc;
    mc.init(3*M_num_elements);
    vector_type mh;
    mh.init(3*M_num_elements);
    vector_type mhs;
    mhs.init(3*M_num_elements);
    vector_type mssh;
    mssh.init(3*M_num_elements);
    vector_type mu;
    mu.init(3*M_num_elements);
    vector_type mv;
    mv.init(3*M_num_elements);

    int cpt = 0;
    double sum_u = 0.;
    double sum_v = 0.;
    double sum_ssh = 0.;
    for (auto it=M_elements.begin(), end=M_elements.end(); it!=end; ++it)
    {
        sum_u = 0.;
        sum_v = 0.;
        //sum_ssh = 0.;
        for (int j=0; j<3; ++j)
        {
            // export ocean field
            // sum_u += M_ocean[it->indices[j]-1];
            // sum_v += M_ocean[it->indices[j]-1+M_num_nodes];

            // export wind field
            // sum_u += M_wind[it->indices[j]-1];
            // sum_v += M_wind[it->indices[j]-1+M_num_nodes];

            // export velocity
            sum_u += M_VT[it->indices[j]-1];
            sum_v += M_VT[it->indices[j]-1+M_num_nodes];
            sum_ssh += M_ssh[it->indices[j]-1];
        }
        sum_u /= 3.;
        sum_v /= 3.;
        sum_ssh /= 3.;

        for (int i=0; i<3; ++i)
        {
            mc(3*cpt+i) = M_conc[cpt];
            mh(3*cpt+i) = M_thick[cpt];
            mhs(3*cpt+i) = M_snow_thick[cpt];
            mssh(3*cpt+i) = sum_ssh;
            mu(3*cpt+i) = sum_u;
            mv(3*cpt+i) = sum_v;
            mx(3*cpt+i) = M_nodes[it->indices[i]-1].coords[0];
            my(3*cpt+i) = M_nodes[it->indices[i]-1].coords[1];
        }

        ++cpt;
    }

    std::string step_str = boost::lexical_cast<std::string>(step)+".m";
    mx.printMatlab("mx" + step_str);
    my.printMatlab("my" + step_str);
    mc.printMatlab("mc" + step_str);
    mh.printMatlab("mh" + step_str);
    mhs.printMatlab("mhs" + step_str);
    mssh.printMatlab("mssh" + step_str);
    mu.printMatlab("mu" + step_str);
    mv.printMatlab("mv" + step_str);
#endif

#if 1
    Exporter exporter;
    std::string fileout;

    if (export_mesh)
    {
        fileout = (boost::format( "%1%/matlab/mesh_%2%.bin" )
                   % Environment::nextsimDir().string()
                   % step ).str();

        std::cout<<"MESH BINARY: Exporter Filename= "<< fileout <<"\n";

        std::fstream meshbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
        exporter.writeMesh(meshbin, M_mesh);
        meshbin.close();

        fileout = (boost::format( "%1%/matlab/mesh_%2%.dat" )
               % Environment::nextsimDir().string()
               % step ).str();

        std::cout<<"RECORD MESH: Exporter Filename= "<< fileout <<"\n";

        std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
        exporter.writeRecord(outrecord,"mesh");
        outrecord.close();
    }


    fileout = (boost::format( "%1%/matlab/field_%2%.bin" )
               % Environment::nextsimDir().string()
               % step ).str();

    std::cout<<"BINARY: Exporter Filename= "<< fileout <<"\n";

    std::fstream outbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
    exporter.writeField(outbin, M_mld, "mld");
    exporter.writeField(outbin, M_VT, "Velocity");
    exporter.writeField(outbin, M_conc, "Concentration");
    exporter.writeField(outbin, M_thick, "Thickness");
    exporter.writeField(outbin, M_wind, "Wind");
    exporter.writeField(outbin, M_ocean, "Ocean");
    exporter.writeField(outbin, M_damage, "Damage");
    exporter.writeField(outbin, M_sst, "SST");
    exporter.writeField(outbin, M_sss, "SSS");
    outbin.close();

    fileout = (boost::format( "%1%/matlab/field_%2%.dat" )
               % Environment::nextsimDir().string()
               % step ).str();

    std::cout<<"RECORD FIELD: Exporter Filename= "<< fileout <<"\n";

    std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
    exporter.writeRecord(outrecord);
    outrecord.close();
#endif
}
} // Nextsim
