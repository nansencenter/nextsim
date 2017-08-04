/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   wimdiscr.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug  3 11:52:35 2015
 */

#include <wimdiscr.hpp>
#include <date_wim.hpp>
#include <tools.hpp>
#ifdef __cplusplus
extern "C"
{
#endif
#include <RTparam_outer.h>
#include <mapx.h>
#ifdef __cplusplus
}
#endif

namespace Wim
{

template<typename T>
void WimDiscr<T>::gridProcessing()
{
    // * sets the following arrays:
    // X_array.resize(boost::extents[nx][ny]);
    // Y_array.resize(boost::extents[nx][ny]);
    // SCUY_array.resize(boost::extents[nx][ny]);
    // SCVX_array.resize(boost::extents[nx][ny]);
    // SCP2_array.resize(boost::extents[nx][ny]);
    // SCP2I_array.resize(boost::extents[nx][ny]);
    // LANDMASK_array.resize(boost::extents[nx][ny]);
    //
    // * if wim.gridfilename is given in the config file,
    // read it from file and apply the stereographic projection
    // to get it in (x,y) coords
    // * else set it manually
    wim_gridfile    = vm["wim.gridfilename"].template as<std::string>();
    if ( wim_gridfile != "" )
    {
        std::cout<<"Getting WIM grid from file: "<<wim_gridfile<<"\n";
        this->readGridFromFile();
    }
    else
    {
        std::cout<<"Generating WIM grid manually...\n";
        chrono.restart();

        nx = vm["wim.nx"].template as<int>();
        ny = vm["wim.ny"].template as<int>();
        dx = vm["wim.dx"].template as<double>();
        dy = vm["wim.dy"].template as<double>();
        x0 = vm["wim.xmin"].template as<double>();
        y0 = vm["wim.ymin"].template as<double>();

        this->gridFromParameters();
        std::cout<<"Grid generation done in "<< chrono.elapsed() <<"s\n";
    }

    this->gridPostProcessing();//save grid, set some other variables

}//gridProcessing

template<typename T>
void WimDiscr<T>::gridProcessing(mesh_type const &mesh_in)
{
    // * sets the following arrays:
    // X_array.resize(boost::extents[nx][ny]);
    // Y_array.resize(boost::extents[nx][ny]);
    // SCUY_array.resize(boost::extents[nx][ny]);
    // SCVX_array.resize(boost::extents[nx][ny]);
    // SCP2_array.resize(boost::extents[nx][ny]);
    // SCP2I_array.resize(boost::extents[nx][ny]);
    // LANDMASK_array.resize(boost::extents[nx][ny]);
    //

    std::cout<<"Generating WIM grid from mesh...\n";
    chrono.restart();

    auto xnod = mesh_in.coordX();
    auto ynod = mesh_in.coordY();
    auto res = NextsimTools::resolution(mesh_in);
    auto x1 = *std::max_element(xnod.begin(),xnod.end());
    auto y1 = *std::max_element(ynod.begin(),ynod.end());
    std::cout<<"Resolution (km) = "<<res/1.e3<<"\n";

    //parameters required for regular grid
    x0  = *std::min_element(xnod.begin(),xnod.end());
    y0  = *std::min_element(ynod.begin(),ynod.end());
    nx  = std::ceil((x1-x0)/res);//round up to lower resolution
    ny  = std::ceil((y1-y0)/res);
    dx  = (x1-x0)/nx;
    dy  = (y1-y0)/ny;

    this->gridFromParameters();
    std::cout<<"Grid generation done in "<< chrono.elapsed() <<"s\n";

    this->gridPostProcessing();//save grid, set some other variables

    //TODO use InterpFromMeshToGrid (or InterpFromMeshToMesh) to set LANDMASK_array

}//gridProcessing

template<typename T>
void WimDiscr<T>::gridPostProcessing()
{
    bool DoSaveGrid = (vm["wim.checkprog"].template as<bool>())
                   || (vm["wim.checkinit"].template as<bool>())
                   || (vm["wim.checkfinal"].template as<bool>())
                   || (vm["nextwim.exportresults"].template as<bool>());

    //std::cout<<" ---before saving\n";
    if (DoSaveGrid)
       this->saveGrid(); //save grid to binary
    //std::cout<<" ---after saving\n";

    //for use in wim_grid
    x_col.resize(nx);
    y_row.resize(ny);
#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < nx; i++)
        x_col[i] = X_array[ny*i];

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int j = 0; j < ny; j++)
        y_row[j] = Y_array[j];

    //Do triangulation
    if (!M_regular)
    {
        //wet cells
        for (int i=0;i<num_p_wim;i++)
            if (LANDMASK_array[i]<.5)
            {
                M_wim_triangulation.nodes_x.push_back(X_array[i]);
                M_wim_triangulation.nodes_y.push_back(Y_array[i]);
                wet_indices.push_back(i);
            }

        int Nnod = wet_indices.size();
        std::cout<<"#nodes for triangulation = "<<Nnod<<"\n";

        //do triangulation
        int* index;
        int Nel  = 0;
		BamgTriangulatex(&index,    //pointer to index              (output)
                         &Nel,      //pointer to num elements       (output)
                         &(M_wim_triangulation.nodes_x)[0],  //pointer to x-coord of nodes   (input)
                         &(M_wim_triangulation.nodes_y)[0],  //pointer to y-coord of nodes   (input)
                         Nnod);     //num nodes                     (input)

        std::cout<<"#elements from triangulation = "<<Nel<<"\n";
        M_wim_triangulation.index.resize(3*Nel);
        for (int i=0; i<3*Nel; i++)
        {
            //std::cout<<"index["<<i<<"] = "<<index[i]<<"\n";
            M_wim_triangulation.index[i]   = index[i];
        }
        xDelete<int>(index);

        //finish M_wim_triangulation (NB don't need coords of elements)
        M_wim_triangulation.initialised     = true;
        M_wim_triangulation.num_nodes       = Nnod;
        M_wim_triangulation.num_elements    = Nel;
    }

}//gridPostProcessing

template<typename T>
void WimDiscr<T>::gridFromParameters()
{
    // * sets the following arrays:
    // X_array.resize(boost::extents[nx][ny]);
    // Y_array.resize(boost::extents[nx][ny]);
    // SCUY_array.resize(boost::extents[nx][ny]);
    // SCVX_array.resize(boost::extents[nx][ny]);
    // SCP2_array.resize(boost::extents[nx][ny]);
    // SCP2I_array.resize(boost::extents[nx][ny]);
    // LANDMASK_array.resize(boost::extents[nx][ny]);

    std::cout<<"nx,ny = "<<nx<<","<<ny<<"\n";
    std::cout<<"dx,dy = "<<dx<<","<<dy<<"\n";
    std::cout<<"xmin,ymin = "<<x0<<","<<y0<<"\n";

    
    num_p_wim    = nx*ny;//number of p points
    num_q_wim    = (nx+1)*(ny+1);//number of q points
    num_u_wim    = (nx+1)*ny;//number of u points
    num_v_wim    = nx*(ny+1);//number of v points
    X_array.resize(num_p_wim);
    Y_array.resize(num_p_wim);
    SCUY_array.resize(num_u_wim);
    SCVX_array.resize(num_v_wim);
    SCP2_array.resize(num_p_wim);
    SCP2I_array.resize(num_p_wim);
    LANDMASK_array.resize(num_p_wim);


    // int thread_id;
    // int total_threads;
    // std::cout<<"MAX THREADS= "<< max_threads <<"\n";

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            X_array[ny*i+j] = x0 + i*dx+.5*dx;
            Y_array[ny*i+j] = y0 + j*dy+.5*dy;
            SCP2_array[ny*i+j] = dx*dy;
            SCP2I_array[ny*i+j] = 1./(dx*dy);

            //add land on 3 edges (upper,lower,RH)
            if (vm["wim.landon3edges"].template as<bool>())
            {
               if (i==nx-1)
               {
                   LANDMASK_array[ny*i+j] = 1.;
               }

               if ((j==0) || (j==ny-1))
               {
                   LANDMASK_array[ny*i+j] = 1.;
               }
            }

        }
    }

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < num_u_wim; i++)
            SCUY_array[i] = dy;

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < num_v_wim; i++)
            SCVX_array[i] = dx;

}//gridFromParameters

template<typename T>
void WimDiscr<T>::saveGrid()
{
    //save grid to binary
    std::string str = vm["wim.outparentdir"].template as<std::string>();
    fs::path path(str);
    path /= "binaries";

    if ( !fs::exists(path) )
        fs::create_directories(path);


    std::string fileout = (boost::format( "%1%/wim_grid.a" ) % path.string()).str();
    std::fstream out(fileout, std::ios::binary | std::ios::out | std::ios::trunc);

    if (out.is_open())
    {
        for (int i = 0; i < X_array.size(); i++)
            out.write((char *)&X_array[i], sizeof(value_type));

        for (int i = 0; i < Y_array.size(); i++)
            out.write((char *)&Y_array[i], sizeof(value_type));

        for (int i = 0; i < SCUY_array.size(); i++)
            out.write((char *)&SCUY_array[i], sizeof(value_type));

        for (int i = 0; i < SCVX_array.size(); i++)
            out.write((char *)&SCVX_array[i], sizeof(value_type));

        for (int i = 0; i < SCP2_array.size(); i++)
            out.write((char *)&SCP2_array[i], sizeof(value_type));

        for (int i = 0; i < SCP2I_array.size(); i++)
            out.write((char *)&SCP2I_array[i], sizeof(value_type));

        for (int i = 0; i < LANDMASK_array.size(); i++)
            out.write((char *)&LANDMASK_array[i], sizeof(value_type));

        out.close();
    }
    else
    {
        std::cout << "Cannot open " << fileout  << "\n";
        std::cerr << "error: open file " << fileout << " for output failed!" <<"\n";
        std::abort();
    }



    // export the txt file for grid field information
    std::string fileoutb = (boost::format( "%1%/wim_grid.b" ) % path.string()).str();
    std::fstream outb(fileoutb, std::ios::out | std::ios::trunc);

    std::string nxstr = std::string(4-std::to_string(nx).length(),'0') + std::to_string(nx);
    std::string nystr = std::string(4-std::to_string(ny).length(),'0') + std::to_string(ny);

    // std::cout<<"-----------nx= "<< nxstr <<"\n";
    // std::cout<<"-----------ny= "<< nystr <<"\n";

    if (outb.is_open())
    {
        outb << std::setw(15) << std::left << "07"  << "    Nrecs    # "<< "Number of records" <<"\n";
        outb << std::setw(15) << std::left << "0"   << "    Norder   # "<< "Storage order [column-major (F/matlab) = 1; row-major (C) = 0]" <<"\n";
        outb << std::setw(15) << std::left << nxstr << "    nx       # "<< "Record length in x direction (elements)" <<"\n";
        outb << std::setw(15) << std::left << nystr << "    ny       # "<< "Record length in y direction (elements)" <<"\n";

        outb <<"\n";

        outb << "Record number and name:" <<"\n";
        outb << std::setw(9) << std::left << "01" << "X" <<"\n";
        outb << std::setw(9) << std::left << "02" << "Y" <<"\n";
        outb << std::setw(9) << std::left << "03" << "scuy" <<"\n";
        outb << std::setw(9) << std::left << "04" << "scvx" <<"\n";
        outb << std::setw(9) << std::left << "05" << "scp2" <<"\n";
        outb << std::setw(9) << std::left << "06" << "scp2i" <<"\n";
        outb << std::setw(9) << std::left << "07" << "LANDMASK" <<"\n";
    }
    else
    {
        std::cout << "Cannot open " << fileoutb  << "\n";
        std::cerr << "error: open file " << fileoutb << " for output failed!" <<"\n";
        std::abort();
    }
}//saveGrid

template<typename T>
void WimDiscr<T>::readGridFromFile()
{
    std::cout<<"Reading grid starts...\n";

    char * senv = ::getenv( "WIMGRIDPATH" );
    if ( (senv == NULL) || (senv[0] == '\0') )
    {
        std::cout << "you must define 'WIMGRIDPATH' environment variable for grid file directory needed for WIM"<<"\n";
        throw std::logic_error("invalid environment variable");
    }

    // start reading first the record file (.b)
    std::string str = std::string(senv);
    fs::path path(str);

    //str = vm["wim.gridfilename"].template as<std::string>();
    str = wim_gridfile;
    std::size_t found = str.find(".");
    if (found != std::string::npos)
    {
        str.erase(std::next( str.begin(), found), str.end());
        std::cout<<"binary filename for wim grid= "<< str <<"\n";
    }
    //str += ".b";
    std::string afilename = (boost::format("%1%/%2%.a") % path.string() % str).str();
    std::string bfilename = (boost::format("%1%/%2%.b") % path.string() % str).str();

    std::ifstream brecord ( bfilename.c_str(), std::ios::in );
    std::vector<int> record(5);
    if (brecord.is_open())
    {
        for (int i=0; i<record.size(); ++i)
        {
            brecord >> record[i];
            brecord.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
    else
    {
        std::cout << "Cannot open " << bfilename  << "\n";
        std::cerr << "error: open file " << bfilename << " for input failed!" <<"\n";
        std::abort();
    }

    bool column_major = record[1];
    nx = record[2];
    ny = record[3];
    std::cout<<"nx= "<< nx <<"\n";
    std::cout<<"ny= "<< ny <<"\n";
    int nbytes = record[4];
    std::cout<<"nbytes= "<< nbytes <<"\n";

    num_p_wim    = nx*ny;//number of p points
    num_q_wim    = (nx+1)*(ny+1);//number of q points
    num_u_wim    = (nx+1)*ny;//number of u points
    num_v_wim    = nx*(ny+1);//number of v points
    value_type_vec PLat_array, PLon_array;

    std::fstream in( afilename, std::ios::binary | std::ios::in);

    int off = num_q_wim*nbytes*2/* skip qlon and qlat*/;

    this->readFromBinary(in, PLon_array, off, std::ios::beg);
    this->readFromBinary(in, PLat_array);

    off = num_u_wim*nbytes*2/* skip ulon and ulat*/;
    off += num_v_wim*nbytes*2/* skip vlon and vlat*/;

    this->readFromBinary(in, SCUY_array, off, std::ios::cur, 1, 0);
    this->readFromBinary(in, SCVX_array, 0, std::ios::cur, 0, 1);
    this->readFromBinary(in, SCP2_array);
    this->readFromBinary(in, LANDMASK_array);

    X_array.resize(num_p_wim);
    Y_array.resize(num_p_wim);
    SCP2I_array.resize(num_p_wim);


    // polar stereographic projection
    mapx_class *map;
    std::string filename = (boost::format("%1%/%2%") % path.string() % "NpsNextsim.mpp").str();
    std::cout<<"stereographic description file= "<< filename <<"\n";

    std::vector<char> _str(filename.begin(), filename.end());
    _str.push_back('\0');

    map = init_mapx(&_str[0]);

    value_type dx_min = 1.e30;
    value_type dy_min = 1.e30;
    value_type dx_max = -1.e30;
    value_type dy_max = -1.e30;
    for (int k = 0; k < num_p_wim; k++)
    {
        double x, y;
        int status = forward_mapx(map,PLat_array[k],PLon_array[k],&x,&y);

        X_array[k]      = x;
        Y_array[k]      = y;
        SCP2I_array[k]  = 1./SCP2_array[k];

        //k=i*ny+j
        int i = k/ny;
        int j = k%ny;
        if (i>1)
        {
            value_type dx_ = X_array[ny*(i)+j]-X_array[ny*(i-1)+j];
            dx_min  = std::min(dx_min,dx_);
            dx_max  = std::max(dx_max,dx_);
        }
        if (j>1)
        {
            value_type dy_ = Y_array[ny*(i)+j]-Y_array[ny*(i)+j-1];
            dy_min  = std::min(dy_min,dy_);
            dy_max  = std::max(dy_max,dy_);
        }
    }
    close_mapx(map);

    dx = .5*(dx_min+dx_max);
    dy = .5*(dy_min+dy_max);
    //dx = X_array[ny]-X_array[0];
    //dy = Y_array[1]-Y_array[0];

    // is the grid regular?
    value_type tol = 1.;//threshold for range in dx,dy (m)
    if ((dx_max-dx_min>tol)||(dy_max-dy_min>tol))
        M_regular = false;

    std::cout<<"dx= "<< dx <<"\n";
    std::cout<<"dy= "<< dy <<"\n";

#if 0
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            //std::cout<<"X_array["<< i <<"]["<< j <<"]= "<< X_array[ny*i+j] <<"\n";
            //std::cout<<"Y_array["<< i <<"]["<< j <<"]= "<< Y_array[ny*i+j] <<"\n";

            //std::cout<<"PLon_array["<< i <<"]["<< j <<"]= "<< PLon_array[ny*i+j] <<"\n";
            //std::cout<<"PLat_array["<< i <<"]["<< j <<"]= "<< PLat_array[ny*i+j] <<"\n";

            //std::cout<<"SCUY_array["<< i <<"]["<< j <<"]= "<< SCUY_array[ny*i+j] <<"\n";
            //std::cout<<"SCVX_array["<< i <<"]["<< j <<"]= "<< SCVX_array[ny*i+j] <<"\n";
            //std::cout<<"SCP2_array["<< i <<"]["<< j <<"]= "<< SCP2_array[ny*i+j] <<"\n";
            //std::cout<<"SCP2I_array["<< i <<"]["<< j <<"]= "<< SCP2I_array[ny*i+j] <<"\n";
            //std::cout<<"LANDMASK_array["<< i <<"]["<< j <<"]= "<< LANDMASK_array[ny*i+j] <<"\n";
        }
    }
#endif

    std::cout<<"Reading grid done...\n";
}//readGridFromFile

template<typename T>
void WimDiscr<T>::readFromBinary(std::fstream &in, value_type_vec& in_array, int off, std::ios_base::seekdir direction, int addx, int addy)
{
    if (off && (in.is_open()))
    {
        in.seekg(off, direction); // skip from the direction (beginning/current/end) position of the file
    }

    int nx_in = nx+addx;
    int ny_in = ny+addy;
    in_array.resize(nx_in*ny_in);

    if (in.is_open())
    {
        // NB assumes files are saved with fortran ordering
        for (int j = 0; j < ny_in; j++)
            for (int i = 0; i < nx_in; i++)
                in.read((char *)&in_array[ny_in*i+j], sizeof(value_type));
    }
    else
    {
        std::cout << "Cannot open " << in << "\n";
        std::cerr << "error: open file " << in << " for input failed!" <<"\n";
        std::abort();
    }
}//readFromBinary

template<typename T>
void WimDiscr<T>::init(int nextsim_cpt)
{
    max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    //set initialised to false for all MeshInfo objects
    nextsim_mesh        = mesh_info_tmp;
    nextsim_mesh_old    = mesh_info_tmp;
    M_wim_triangulation = mesh_info_tmp;

    // wim grid generation/reading
    this->gridProcessing();

    //parameters
    this->initConstant(nextsim_cpt);

    // call assign to set sizes of arrays
    // and initialise some arrays that are stationary in time
    this->assign();
    this->assignSpatial();

    //set global counter to 0
    M_cpt = 0;

    std::cout<<"wim.init() finished\n";

}//end ::init()

template<typename T>
void WimDiscr<T>::init(mesh_type const &mesh_in,int nextsim_cpt)
{
    max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    //set initialised to false for all MeshInfo objects
    nextsim_mesh        = mesh_info_tmp;
    nextsim_mesh_old    = mesh_info_tmp;
    M_wim_triangulation = mesh_info_tmp;

    //set global counter to 0
    M_cpt = 0;

    // can either init from MESH or ON mesh
    if (!M_wim_on_mesh)
        this->gridProcessing(mesh_in);
    else
        this->setMesh(mesh_in);

    //parameters
    this->initConstant(nextsim_cpt);

    // call assign to set sizes of arrays
    // and initialise some arrays that are stationary in time
    this->assign();
    this->assignSpatial();

    std::cout<<"wim.init() finished\n";

}//end ::init()

template<typename T>
void WimDiscr<T>::initConstant(int nextsim_cpt)
{

    // needed to set sizes of variables
    // - grid should have been set before this function has been called
    M_num_elements  = X_array.size();

    // wave parameters
    nwavedirn   = vm["wim.nwavedirn"].template as<int>();
    nwavefreq   = vm["wim.nwavefreq"].template as<int>();
    Tmin        = vm["wim.tmin"].template as<double>(); /* 2.5 */
    Tmax        = vm["wim.tmax"].template as<double>(); /* 25. */
    ref_Hs_ice  = vm["wim.refhsice"].template as<bool>();
    atten       = vm["wim.atten"].template as<bool>();
    useicevel   = vm["wim.useicevel"].template as<bool>();
    M_steady    = vm["wim.steady"].template as<bool>();
    scatmod     = vm["wim.scatmod"].template as<std::string>();

    // ice parameters
    breaking            = vm["wim.breaking"].template as<bool>();
    fsdopt              = vm["wim.fsdopt"].template as<std::string>();
    dfloe_pack_init     = vm["wim.dfloepackinit"].template as<double>(); /* 300.0 */
    dfloe_pack_thresh   = vm["wim.dfloepackthresh"].template as<double>(); /* 400.0 */
    young               = vm["wim.young"].template as<double>();
    drag_rp             = vm["wim.dragrp"].template as<double>();

    // ==============================================================================
    //local diagnostics
    int wim_itest = vm["wim.itest"].template as<int>();
    int wim_jtest = vm["wim.jtest"].template as<int>();

    M_itest  = -1;
    if(M_wim_on_mesh)
        M_itest  = wim_itest;//could still be -1
    else if((wim_itest>=0)&&(wim_jtest>=0))
    {
        //(!M_wim_on_mesh) && positive itest,jtest input
        std::cout<<"\nitest,jtest: "<<wim_itest<<","<<wim_jtest<<"\n";
        M_itest  = wim_itest*ny+wim_jtest;
        if ((wim_itest>=nx)||(wim_jtest>=ny))
        {
            std::cout<<"\nitest,jtest: "<<wim_itest<<","<<wim_jtest<<"\n";
            std::cout<<"\nnx,ny: "<<nx<<","<<ny<<"\n";
            throw std::runtime_error("wim.itest/jtest out of range");
        }
    }

    //global test index
    std::cout<<"\nWIM diagnostics to be done at global index: "
        <<M_itest<<"(of "<<M_num_elements<<")\n";
    if (M_itest>=M_num_elements)
        throw std::runtime_error("M_itest out of range");
    // ==============================================================================


    //numerical parameters
    M_advdim    = vm["wim.advdim"].template as<int>();
    M_advopt    = vm["wim.advopt"].template as<std::string>();
    M_cfl       = vm["wim.cfl"].template as<double>();

    //need nghost>=4 for WENO advection
    nghost  = 4;
    nbdx    = nghost;
    nbdy    = nghost;

    if (useicevel)
        throw std::runtime_error("useicevel=true not implemented\n");

    if(M_steady && (M_advopt=="xy-periodic"))
    {
        std::string tmps = "advopt = xy-periodic and steady options incompatible - ";
        tmps += "use y-periodic with steady or turn off steady\n";
        throw std::runtime_error(tmps);
    }

    if (M_advdim == 1)
        nbdy = 0;

    nxext = nx+2*nbdx;
    nyext = ny+2*nbdy;//ny if M_advdim==1

    gravity             = 9.81;
    rhowtr              = 1025.;
    rhoice              = 922.5;
    poisson             = 0.3;
    dmin                = 20.;
    xi                  = 2.;
    fragility           = 0.9;
    cice_min            = vm["wim.cicemin"].template as<double>();
    dfloe_miz_thresh    = 200.;

    vbf = 0.1;//brine volume fraction
    vb = vbf;
    sigma_c  = (1.76e+6)*std::exp(-5.88*std::sqrt(vbf));//flexural strength (Pa)
    epsc = sigma_c/young;//breaking strain
    flex_rig_coeff = young/(12.0*(1-std::pow(poisson,2.)));

    //some options need to be disabled if being called from nextsim
    docoupling = !( vm.count("simul.use_wim")==0 );
    M_restart_time  = 0.;
    if (!docoupling)
    {
        //set duration of call to wim.run() from wim.duration
        duration = vm["wim.duration"].template as<double>();

        //get initial time from wim.initialtime
        init_time_str = vm["wim.initialtime"].as<std::string>();

        break_on_mesh   =  false;

        //save options from wimoptions.cpp
        this->saveOptionsLog();
    }
    else
    {
        //set duration of call to wim.run() from nextwim.couplingfreq
        value_type nextsim_time_step = vm["simul.timestep"].template as<double>();
        duration   = vm["nextwim.couplingfreq"].template as<int>()*nextsim_time_step;
                     
        //get initial time from simul.time_init
        init_time_str  = vm["simul.time_init"].template as<std::string>();

        //if using restart, need to calculate shift from initial
        //time
        M_restart_time = nextsim_cpt*nextsim_time_step;//model time of current call to wim
 
        break_on_mesh   =  ( vm["nextwim.coupling-option"].template as<std::string>() == "breaking_on_mesh");
    }

    M_update_time   = M_restart_time;//reset at last update

#if 1
    //print initial & restart times
    std::cout<<"initial time = "<<init_time_str<<"\n";
    auto time_str1     = ptime(init_time_str,M_restart_time);
    std::cout<<"restart time = "<<time_str1<<"\n";
#endif

    //no of cosines/sines to use - for isotropic scattering code
    ncs = std::round(nwavedirn/2);

}//end ::initConstant()


template<typename T>
void WimDiscr<T>::assign()
{
    // this doesn't need to be called each time wim.run is called
    // - sets arrays with dimensions depending on frequency, wave dirn (never change)
    wt_simp.resize(nwavefreq);
    wt_om.resize(nwavefreq);

    freq_vec.resize(nwavefreq);
    vec_period.resize(nwavefreq);
    wlng.resize(nwavefreq);
    ag.resize(nwavefreq);
    ap.resize(nwavefreq);

    //4d vec
    sdf_dir.resize(nwavefreq);

    // =============================================
    // set frequencies to use (freq_vec)
    value_type tp_in = vm["wim.tpinc"].template as<double>();
    if (nwavefreq == 1)
        freq_vec[0] = 1./tp_in;
    else
    {
        // multiple frequencies
        fmin = 1./Tmax;
        fmax = 1./Tmin;
        df = (fmax-fmin)/(nwavefreq-1);

        for (int fq = 0; fq < nwavefreq; fq++)
            freq_vec[fq] = fmin+fq*df;
    }
    // =============================================


    // =============================================
    // set directions to use (wavedir)
    value_type avgdir = vm["wim.mwdinc"].template as<double>();
    if (nwavedirn == 1)
        wavedir.push_back(avgdir);
    else
    {
        value_type theta_max = 90.;
        value_type theta_min = -270.;
        value_type dtheta = (theta_min-theta_max)/nwavedirn;

        for (int nth = 0; nth < nwavedirn; nth++)
            wavedir.push_back(theta_max+nth*dtheta);
    }
    // =============================================


    // =============================================
    // weights for integration with respect to frequency
    if (nwavefreq==1)
        wt_om[0] = 1.;
    else
    {
        // (wt_simp  = weights for Simpson's rule)
        std::fill(wt_simp.begin(), wt_simp.end(), 2.);
        wt_simp[0] = 1.;
        wt_simp[nwavefreq-1] = 1.;

        int w = 1;
        while (w < nwavefreq-1)
        {
            wt_simp[w] = 4.;
            w +=2;
        }

        dom   = 2*PI*(freq_vec[nwavefreq-1]-freq_vec[0])/(nwavefreq-1);
        wt_om = wt_simp;
        std::for_each(wt_om.begin(), wt_om.end(), [&](value_type& f){ f = dom*f/3.0; });
    }
    // =============================================


    // =============================================
    // vector of periods, open water wavelength, phase and group velocities

    // periods
    vec_period = freq_vec;
    std::for_each(vec_period.begin(), vec_period.end(), [&](value_type& f){ f = 1./f; });

    // open water wavelengths
    wlng = freq_vec;
    std::for_each(wlng.begin(), wlng.end(), [&](value_type& f){ f = gravity/(2*PI*std::pow(f,2.)); });

    // open water phase velocities
    ap = wlng;
    std::for_each(ap.begin(), ap.end(), [&](value_type& f){ f = std::sqrt(gravity*f/(2*PI)); });

    // open water group velocities
    ag = ap;
    std::for_each(ag.begin(), ag.end(), [&](value_type& f){ f = f/2. ; });
    // =============================================

}//end: assign()

template<typename T>
void WimDiscr<T>::assignSpatial()
{
    // this needs to be called each time grid or mesh changes
    // ie initially, and if(M_wim_on_mesh), after regridding
    // set sizes of arrays, initialises some others that are constant in time

    value_type_vec vtmp(M_num_elements,0.);//temp vector of correct size, but filled with zeros;

    //3D var's
    // - space and dirn
    //sdf3d_dir_temp.resize(boost::extents[num_p_wim][nwavedirn]);

    // - space and freq
#if 0
    ag_eff.assign(nwavefreq,vtmp);
    ap_eff.assign(nwavefreq,vtmp);
    wlng_ice.assign(nwavefreq,vtmp);
    atten_nond.assign(nwavefreq,vtmp);
    damping.assign(nwavefreq,vtmp);
    disp_ratio.assign(nwavefreq,vtmp);
#else
    ag_eff.resize(boost::extents[M_num_elements][nwavefreq]);
    ap_eff.resize(boost::extents[M_num_elements][nwavefreq]);
    wlng_ice.resize(boost::extents[M_num_elements][nwavefreq]);
    atten_nond.resize(boost::extents[M_num_elements][nwavefreq]);
    damping.resize(boost::extents[M_num_elements][nwavefreq]);
    disp_ratio.resize(boost::extents[M_num_elements][nwavefreq]);
#endif

    //2D var's
    ag2d_eff_temp   = vtmp;
    wave_mask       = vtmp;

    wim_ice.mask    = vtmp;
    wim_ice.conc    = vtmp;
    wim_ice.thick   = vtmp;
    wim_ice.dfloe   = vtmp;
    wim_ice.nfloes  = vtmp;

    dave        = vtmp;
    atten_dim   = vtmp;
    damp_dim    = vtmp;

    tau_x           = vtmp;
    tau_y           = vtmp;
    stokes_drift_x  = vtmp;
    stokes_drift_y  = vtmp;
    Hs              = vtmp;
    Tp              = vtmp;
    mwd             = vtmp;

    S_freq              = vtmp;
    taux_om             = vtmp;
    tauy_om             = vtmp;
    stokes_drift_x_om   = vtmp;
    stokes_drift_y_om   = vtmp;

    if (sdf_dir[0].size()==0)
        //set dir spec to 0. on 1st call to init/assign
        //std::cout<<"Init sdf_dir in wim.assign()\n";
        for (auto it=sdf_dir.begin();it!=sdf_dir.end();it++)
            it->assign(nwavedirn,vtmp);
    // =============================================

}//end: assignSpatial()


template<typename T>
void WimDiscr<T>::update()
{

    //====================================================
    // set incident wave spec where wave_mask==1;
    // this also sets sdf_dir in this region
    this->setIncWaveSpec();
    //====================================================


    //====================================================
    // update attenuation coefficients, wavelengths and phase/group velocities
    this->updateWaveMedium();
    //====================================================


    // ====================================================================================
    // set time step
    // - this can change with time if using group velocity for ice
    // - NB needs to be done after updateWaveMedium
    amax = *std::max_element(ag_eff.data(),ag_eff.data() + ag_eff.num_elements());
    //std::cout<<"dx,amax,M_cfl= "<< dx<<","<<amax<<","<<M_cfl <<"\n";
    M_timestep = M_cfl*dx/amax;

    //reduce time step slightly (if necessary) to make duration an integer multiple of M_timestep
    nt = std::ceil(duration/M_timestep);
    M_timestep = duration/nt;
    //std::cout<<"M_timestep,nt= "<< M_timestep<<","<<nt<<"\n";
    // ====================================================================================


#if 0
    std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    std::cout<<"icec_max= "<< *std::max_element(icec.data(), icec.data()+icec.num_elements()) <<"\n";
    std::cout<<"icec_min= "<< *std::min_element(icec.data(), icec.data()+icec.num_elements()) <<"\n";
    std::cout<<"icec_acc= "<< std::accumulate(icec.data(), icec.data()+icec.num_elements(),0.) <<"\n";

    std::cout<<"--------------------------------------------------------\n";

    std::cout<<"iceh_max= "<< *std::max_element(iceh.data(), iceh.data()+iceh.num_elements()) <<"\n";
    std::cout<<"iceh_min= "<< *std::min_element(iceh.data(), iceh.data()+iceh.num_elements()) <<"\n";
    std::cout<<"iceh_acc= "<< std::accumulate(iceh.data(), iceh.data()+iceh.num_elements(),0.) <<"\n";

    std::cout<<"--------------------------------------------------------\n";

    std::cout<<"dfloe_max= "<< *std::max_element(dfloe.begin(), dfloe.end()) <<"\n";
    std::cout<<"dfloe_min= "<< *std::min_element(dfloe.begin(), dfloe.end()) <<"\n";
    std::cout<<"dfloe_acc= "<< std::accumulate(dfloe.begin(), dfloe.end(),0.) <<"\n";

    std::cout<<"--------------------------------------------------------\n";

    std::cout<<"nfloes_max= "<< *std::max_element(nfloes_in.begin(), nfloes_in.end()) <<"\n";
    std::cout<<"nfloes_min= "<< *std::min_element(nfloes_in.begin(), nfloes_in.end()) <<"\n";
    std::cout<<"nfloes_acc= "<< std::accumulate(nfloes_in.begin(), nfloes_in.end(),0.) <<"\n";
#endif


}//end: update()


template<typename T>
void WimDiscr<T>::updateWaveMedium()
{
    //updates attenuation coefficients, wavelengths and phase/group velocities

    // =============================================================================================
    //std::cout<<"attenuation loop starts (big loop)\n";
    std::fill( disp_ratio.data(), disp_ratio.data() + disp_ratio.num_elements(), 1. );
    for (int fq = 0; fq < nwavefreq; fq++)
    {
#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
        {
            double params[5];
            params[0] = young;
            params[1] = gravity;
            params[2] = rhowtr;
            params[3] = rhoice;
            params[4] = poisson;

            double outputs[8];

            if (wim_ice.mask[i] >.5)
            {
                om = 2*PI*freq_vec[fq];
                guess = std::pow(om,2.)/gravity;

                //if (i==M_itest)
                //{
                //    std::cout<<"fq,om,freq,period = "<<fq<<","<<om<<","<<freq_vec[fq]<<","<<1./freq_vec[fq]<<"\n";
                //    std::cout<<"guess (open water) = "<<guess<<"\n";
                //}
                if (fq > 0)
                {
                    //improve guess by using output of last call;
                    guess = 2*PI/wlng_ice[i][fq-1];
                }

                RTparam_outer(outputs,wim_ice.thick[i],double(om),double(drag_rp),double(guess),params);

                value_type kice, kwtr, int_adm, modT, argR, argT;

                damping[i][fq] = outputs[0];
                kice = outputs[1];
                kwtr = outputs[2];
                int_adm = outputs[3];
                atten_nond[i][fq] = outputs[4];
                modT = outputs[5];
                argR = outputs[6];
                argT = outputs[7];

                double tmp = 1.;
                for (int io = 0; io<8; io++)
                    tmp *= outputs[io];

                if (std::isnan(tmp))
                {
                    std::cout<<"found NaN in some outputs of RTparam_outer at i,fq="<<i<<","<<fq<<"\n";
                    std::cout<<"\ninputs to RTparam_outer:\n";
                    std::cout<<"h = "<<wim_ice.thick[i]<<"\n";
                    std::cout<<"om = "<<om<<"\n";
                    std::cout<<"drag_rp = "<<drag_rp<<"\n";
                    std::cout<<"guess = "<<guess<<"\n";
                    //
                    std::cout<<"\noutputs from RTparam_outer:\n";
                    std::cout<<"damping = "<<damping[i][fq]<<"\n";
                    std::cout<<"kice = "<<kice<<"\n";
                    std::cout<<"kwtr = "<<kwtr<<"\n";
                    std::cout<<"int_adm = "<<int_adm<<"\n";
                    std::cout<<"atten_nond = "<<atten_nond[i][fq]<<"\n";
                    std::cout<<"modT = "<<modT<<"\n";
                    std::cout<<"argR = "<<argR<<"\n";
                    std::cout<<"argT = "<<argT<<"\n";
                    throw std::runtime_error("some outputs of RTparam_outer have NaN");
                }

                //convert amplitude in water to amplitude in ice
                disp_ratio[i][fq] = (kice*modT)/kwtr;

                //wavelength to use in ice
                if (1)
                {
                   //use ice wavelength TODO make an option?
                   wlng_ice[i][fq] = 2*PI/kice;
                }
                else
                {
                   //use water wavelength instead of ice wavelength
                   wlng_ice[i][fq] = wlng[fq];
                }

                //group and phase velocities to use in ice
                if (!useicevel)
                {
                   //water group and phase velocities
                   //(ice ones not implemented)
                   ag_eff[i][fq] = ag[fq];
                   ap_eff[i][fq] = ap[fq];
                }

            }//ice
            else
            {
                ag_eff[i][fq] = ag[fq];
                ap_eff[i][fq] = ap[fq];
                wlng_ice[i][fq] = wlng[fq];
            }//water
        }//end i loop
    }//end freq loop
    //std::cout<<"big loop done\n";
    // =============================================================================================

    if (!atten)
    {
        std::fill( atten_nond.data(), atten_nond.data() + atten_nond.num_elements(), 0. );
        std::fill( damping.data(), damping.data() + damping.num_elements(), 0. );
    }

}//end: update()


template<typename T>
void WimDiscr<T>::idealWaveFields(value_type const xfac)
{
    M_initialised_waves = true;

    value_type x0   = *std::min_element(X_array.begin(),X_array.end());
    value_type xmax = *std::max_element(X_array.begin(),X_array.end());

    //waves initialised for x<x_edge
    value_type x_edge = 0.5*(x0+xmax)-xfac*(0.5*(xmax-x0));

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {
        if ((X_array[i] < x_edge) && (LANDMASK_array[i]<1.))
        {
           wave_mask[i] = 1.;
           Hs [i] = vm["wim.hsinc"].template as<double>();
           Tp [i] = vm["wim.tpinc"].template as<double>();
           mwd[i] = vm["wim.mwdinc"].template as<double>();
           //std::cout<<Hs[i]<<" "<<Tp[i]<<" "<<mwd[i];
        }
        else
           wave_mask[i] = 0.;
    }
}//idealWaveFields


template<typename T>
void WimDiscr<T>::inputWaveFields(value_type_vec const& swh_in,
                                  value_type_vec const& mwp_in,
                                  value_type_vec const& mwd_in)
{
    M_initialised_waves = true;

    bool checkincwaves = vm["wim.checkincwaves"].template as<bool>();
    if (checkincwaves)
    {
        //these arrays are only needed for diagnostics
        swh_in_array.resize(num_p_wim);
        mwp_in_array.resize(num_p_wim);
        mwd_in_array.resize(num_p_wim);
    }

#if 0
    // test routine against ideal version
    // - also need to uncomment printout loop after Hs, Tp, mwd & wave_mask are set
    std::cout<<"assign wave fields 2\n";
    auto Hs_old=Hs;
    auto Tp_old=Tp;
    auto mwd_old=mwd;
    auto WM_old=wave_mask;
    this->idealWaveFields(.8);//.8 determines ice edge location
    auto Hs_ideal = Hs;
    auto Tp_ideal = Tp;
    auto mwd_ideal = mwd;
    auto WM_ideal=wave_mask;
    Hs  = Hs_old;
    Tp  = Tp_old;
    mwd = mwd_old;
    wave_mask = WM_old;
#endif


    double Hs_min=100.;
    double Hs_max=0.;
    double Hs_min_i=100.;
    double Hs_max_i=0.;
    double Hs_min_ice=100.;
    double Hs_max_ice=0.;

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < num_p_wim; i++)
    {
        if (checkincwaves)
        {
            //only needed for diagnostics
            swh_in_array[i] = swh_in[i];
            mwp_in_array[i] = mwp_in[i];
            mwd_in_array[i] = mwd_in[i];
        }

        if ((wim_ice.mask[i]<.5)                    //not ice
            &&(LANDMASK_array[i]<.5)                //not land
            &&(swh_in[i]>1.e-3)&&(mwp_in[i]>1.e-8)  //some waves
            &&(mwp_in[i]<1.5*Tmax))                 //wave period not too big
        {
           wave_mask[i] = 1.;
           Hs [i] = swh_in[i];
           Tp [i] = mwp_in[i];
           mwd[i] = mwd_in[i];
           //std::cout<<Hs[i]<<" "<<Tp[i]<<" "<<mwd[i];
        }
        else
           wave_mask[i] = 0.;
        
        if (wim_ice.mask[i]>0.)
        {
            Hs_min_ice    = std::min(Hs_min_ice,Hs[i]);
            Hs_max_ice    = std::max(Hs_max_ice,Hs[i]);
        }
        Hs_min = std::min(Hs_min,Hs[i]);
        Hs_max = std::max(Hs_max,Hs[i]);
        Hs_min_i = std::min(Hs_min_i,swh_in[i]);
        Hs_max_i = std::max(Hs_max_i,swh_in[i]);
    }

    std::cout<<"Hs range from outside:\n";
    std::cout<<"Hs_min = "<<Hs_min_i<<"\n";
    std::cout<<"Hs_max = "<<Hs_max_i<<"\n";
    std::cout<<"Hs_min (processed) = "<<Hs_min<<"\n";
    std::cout<<"Hs_max (processed) = "<<Hs_max<<"\n";
    std::cout<<"Hs_min (in ice) = "<<Hs_min_ice<<"\n";
    std::cout<<"Hs_max (in ice) = "<<Hs_max_ice<<"\n";

#if 0
    // print out for test vs ideal case
#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < num_p_wim; i++)
    {
        std::cout<<"Hs["<<i<<"]="<<Hs[i]<<","<<Hs_ideal[i]<<","<<swh_in[i]<<"\n";
         std::cout<<"Tp["<<i<<"]="<<Tp[i]<<","<<Tp_ideal[i]<<","<<mwp_in[i]<<"\n";
         std::cout<<"mwd["<<i<<"]="<<mwd[i]<<","<<mwd_ideal[i]<<","<<mwd_in[i]<<"\n";
         std::cout<<"wave_mask["<<i<<"]="<<wave_mask[i]<<","<<WM_ideal[i]<<"\n\n";
    }
    std::abort();
#endif
}//inputWaveFields()

template<typename T>
void WimDiscr<T>::setIncWaveSpec()
{

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < num_p_wim; i++)
    {
        if (wave_mask[i] == 1.)
        {
            std::vector<value_type> Sfreq(nwavefreq);
            std::vector<value_type> theta_fac(nwavedirn);
            value_type f1, f2, f3, t_m, om_m, chi, om;

            // ============================================================
            // frequency spectrum
            if (nwavefreq == 1)
                Sfreq[0] = std::pow(Hs[i]/4.0,2.);
            else
            {
                for (int fq = 0; fq < nwavefreq; fq++)
                {
                    om        = 2*PI*freq_vec[fq];
                    t_m       = 2*PI/om;
                    om_m      = 2*PI/Tp[i];
                    f1        = (5.0/16.0)*std::pow(Hs[i],2.)*std::pow(om_m,4.);
                    f2        = 1.0/std::pow(om,5.);
                    f3        = std::exp(-1.25*std::pow(t_m/Tp[i],4.));
                    Sfreq[fq] = f1*f2*f3;

#if 0
                    if (i==M_itest)
                    {
                        std::cout<<"i,fq="<<i<<","<<fq<<"\n";
                        std::cout<<"om,f1,f2,f3="<<om<<","<<f1<<","<<f2<<","<<f3<<"\n";
                        std::cout<<"Tp,Hs,Sfreq="
                            <<Tp[i]<<","<<Hs[i]<<","<<Sfreq[fq]<<"\n";
                    }
#endif
                }
            }//multiple frequencies
            // ============================================================


            // ============================================================
            // directional spreading
            if (nwavedirn == 1)
                theta_fac[0]    = 1.;
            else
            {
                value_type dtheta = std::abs(wavedir[1]-wavedir[0]);

                //if (mwd[i]!=0.)
                //   std::cout<<"dir-frac ("<<i<<")"<<std::endl;
                for (int nth = 0; nth < nwavedirn; nth++)
                {
#if 0
                    //less accurate way of calculating spreading
                    //(sample cos^2 at mid-point of interval)
                    chi = PI*(wavedir[nth]-mwd[i])/180.0;
                    if (std::cos(chi) > 0.)
                        theta_fac[nth] = 2.0*std::pow(std::cos(chi),2.)/PI;
                    else
                        theta_fac[nth] = 0.;
#else
                    //more accurate way of calculating spreading
                    //(integrate cos^2 over interval)
                    theta_fac[nth] = 180./(PI*dtheta)*thetaDirFrac(
                            wavedir[nth]-dtheta/2., dtheta, mwd[i] );
#endif

                    //if (Hs[i*ny+j]!=0.)
                    //   std::cout<<wavedir[nth]<<" "<<mwd[i]<<" "
                    //            <<theta_fac[nth]<<std::endl;
                }
            }//multiple directions
            // ============================================================


            // ============================================================
            // combine freq and dir
            for (int fq = 0; fq < nwavefreq; fq++)
                for (int nth = 0; nth < nwavedirn; nth++)
                {
                    // set sdf_dir to inc waves each time new waves are input called
                    // NB but only inside the wave mask
                    sdf_dir[fq][nth][i] = Sfreq[fq]*theta_fac[nth];

#if 0
                    if (i==M_itest)
                    {
                        std::cout<<"fq,nth="<<fq<<","<<nth<<"\n";
                        std::cout<<"wave_mask,Sfreq,theta_fac="
                            <<wave_mask[i]<<","<<Sfreq[fq]<<","<<theta_fac[nth]<<"\n";
                    }
#endif

                }
            // ============================================================
        }//inside wave_mask

#if 0
        //test thetaDirFrac function
        int Ntst            = 110;
        value_type dth_tst  = 360./Ntst;
        value_type dir_tst  = 0.;

        value_type integral1   = thetaDirFrac(dir_tst,360.,270.);
        value_type integral2   = thetaDirFrac(dir_tst,360.,90.);
        value_type integral3   = thetaDirFrac(dir_tst,360.,330.);
        value_type integral4   = thetaDirFrac(dir_tst,360.,0.);
        value_type integral5   = thetaDirFrac(dir_tst,360.,360.);
        std::cout<<"test total thetaDirFrac integral (=1?) = "<<integral1<<"\n";
        std::cout<<"test total thetaDirFrac integral (=1?) = "<<integral2<<"\n";
        std::cout<<"test total thetaDirFrac integral (=1?) = "<<integral3<<"\n";
        std::cout<<"test total thetaDirFrac integral (=1?) = "<<integral4<<"\n";
        std::cout<<"test total thetaDirFrac integral (=1?) = "<<integral5<<"\n";
        std::cout<<"\n";

        integral1   = 0.;
        integral2   = 0.;
        integral3   = 0.;
        integral4   = 0.;
        integral5   = 0.;
        for (int itst=0; itst<Ntst; itst++)
        {
            integral1  += thetaDirFrac(dir_tst,dth_tst,270.);
            integral2  += thetaDirFrac(dir_tst,dth_tst,90.);
            integral3  += thetaDirFrac(dir_tst,dth_tst,330.);
            integral4  += thetaDirFrac(dir_tst,dth_tst,0.);
            integral5  += thetaDirFrac(dir_tst,dth_tst,360.);
            dir_tst    += dth_tst;
        }
        std::cout<<"test thetaDirFrac integral (=1?) = "<<integral1<<"\n";
        std::cout<<"test thetaDirFrac integral (=1?) = "<<integral2<<"\n";
        std::cout<<"test thetaDirFrac integral (=1?) = "<<integral3<<"\n";
        std::cout<<"test thetaDirFrac integral (=1?) = "<<integral4<<"\n";
        std::cout<<"test thetaDirFrac integral (=1?) = "<<integral5<<"\n";
        std::abort();
#endif

    }//end i loop

#if 0
    if(M_itest>=0)
    {
        std::cout<<"i="<<M_itest<<"\n";
        for (int fq = 0; fq < nwavefreq; fq++)
            for (int nth = 0; nth < nwavedirn; nth++)
            {
                std::cout<<"fq,nth="<<fq<<","<<nth<<"\n";
                std::cout<<"sdf_dir (setIncWaveSpec) ="<<sdf_dir[M_itest][nth][fq]<<"\n";
            }
    }
#endif

}//setIncWaveSpec


template<typename T>
void WimDiscr<T>::idealIceFields(value_type const xfac)
{
   value_type x0   = *std::min_element(X_array.begin(),X_array.end());
   value_type xmax = *std::max_element(X_array.begin(),X_array.end());

   //ice initialised for x>=x_edge
   value_type x_edge = 0.5*(x0+xmax)-xfac*(0.5*(xmax-x0));

   (wim_ice.conc  ).assign(M_num_elements,0.);
   (wim_ice.vol   ).assign(M_num_elements,0.);
   (wim_ice.nfloes).assign(M_num_elements,0.);
    value_type unifc = vm["wim.unifc"].template as<double>(); /* 0.7 */
    value_type unifh = vm["wim.unifh"].template as<double>(); /* 2.0 */
#pragma omp parallel for num_threads(max_threads) collapse(1)
   for (int i = 0; i < M_num_elements; i++)
   {
       if ((X_array[i] >= x_edge) && (LANDMASK_array[i]<.5))
       {
           wim_ice.conc[i]      = unifc;
           wim_ice.vol[i]       = unifc*unifh;
           wim_ice.nfloes[i]    = unifc/std::pow(dfloe_pack_init,2);
           //std::cout<<Hs[i]<<" "<<Tp[i]<<" "<<mwd[i];
       }
   }
   this->transformIce(wim_ice);
}//idealIceFields()


template<typename T>
void WimDiscr<T>::inputIceFields(value_type_vec const& icec_in,   // conc
                                 value_type_vec const& iceh_in,   // effective thickness
                                 value_type_vec const& nfloes_in) // c/Dmax^2
{

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {
        if (icec_in[i] < vm["wim.cicemin"].template as<double>())
        {
            wim_ice.mask[i]     = 0.;
            wim_ice.conc[i]     = 0.;
            wim_ice.thick[i]    = 0.;
            wim_ice.nfloes[i]   = 0.;
        }//water
        else
        {
            wim_ice.mask[i]     = 1.;
            wim_ice.conc[i]     = icec_in[i];
            wim_ice.thick[i]    = iceh_in[i]/icec_in[i];//convert to true thickness
            wim_ice.nfloes[i]   = nfloes_in[i];
        }//ice

        //nfloe->dfloe
        wim_ice.dfloe[i] = this->nfloesToDfloe(wim_ice.nfloes[i],wim_ice.conc[i]);

#if 1
        //check ranges of inputs
        if ((wim_ice.thick[i]<0.)||(wim_ice.thick[i]>50.))
        {
            std::cout<<"thickness on WIM grid out of range for i="<<i<<"\n";
            std::cout<<"thick,h,c="<<iceh_in[i]<<","<<wim_ice.thick[i]<<","<<wim_ice.conc[i]<<"\n";
            throw std::runtime_error("thickness on WIM grid out of range");
        }
        if ((wim_ice.conc[i]<0.)||(wim_ice.conc[i]>1.))
        {
            std::cout<<"conc on WIM grid out of range for i="<<i<<"\n";
            std::cout<<"c="<<wim_ice.conc[i]<<"\n";
            throw std::runtime_error("conc on WIM grid out of range");
        }
#endif

    }//end i loop
}//inputIceFields


#if 0
template<typename T>
void getWimCenters(value_type &x,value_type &y,value_type const& rotangle)
{
    //get x coord of nodes (rotated)
    std::vector<value_type> x(nx*ny);
    std::vector<value_type> y(nx*ny);
    int kpt = 0;
    double cos_rotangle=std::cos(rotangle);
    double sin_rotangle=std::sin(rotangle);
    for (int i = 0; i<nx; i++)
    {
        for (int j = 0; j<nx; j++)
        {
        x[kpt] = cos_rotangle*X_array[ny*i+j] + sin_rotangle*Y_array[ny*i+j];
        y[kpt] = -sin_rotangle*X_array[ny*i+j] + cos_rotangle*Y_array[ny*i+j];
        ++kpt;
    }
}

template<typename T>
void getWimShape() const
{
    //get x coord of nodes (rotated)
    std::vector<int> shape={nx,ny};
    return shape;
}
#endif

template<typename T>
void WimDiscr<T>::timeStep()
{
    std::fill( tau_x.begin(), tau_x.end(), 0. );
    std::fill( tau_y.begin(), tau_y.end(), 0. );
    std::fill( stokes_drift_x.begin(), stokes_drift_x.end(), 0. );
    std::fill( stokes_drift_y.begin(), stokes_drift_y.end(), 0. );

    std::vector<value_type> mom0, mom2, var_strain, mom0w, mom2w;
    mom0      .resize(M_num_elements);
    mom2      .resize(M_num_elements);
    var_strain.resize(M_num_elements);
    mom0w     .resize(M_num_elements);
    mom2w     .resize(M_num_elements);

    value_type E_tot, wlng_crest, Dc;
    value_type F, kicel, om, ommin, ommax, om1, lam1, lam2, dom, c1d, tmp;
    int jcrest;
    bool break_criterion,test_ij;

    value_type t_step = this->getModelTime();//seconds from initial time
    std::string timestpstr = ptime(init_time_str, t_step);

    // dump local diagnostic file
    // - directory to put it
    std::string outdir = vm["wim.outparentdir"].template as<std::string>();
    fs::path path(outdir);
    path /= "diagnostics/local";
    if ( !fs::exists(path) )
        fs::create_directories(path);

    //set file name and open
    std::string diagfile   = (boost::format( "%1%/WIMdiagnostics_local%2%.txt" )
            % path.string() % timestpstr).str();

    if ( M_dump_diag )
    {

       std::fstream diagID(diagfile, std::ios::out | std::ios::trunc);
       if (diagID.is_open())
       {
            //
            //diagID << "20150101" << " # date";
            //diagID << "04:06:35" << " # time";
            //diagID << "042003" << " # model day";
            //diagID << "14794.52055" << " # model second";
            diagID << timestpstr << " # model time\n";
            if (M_wim_on_mesh)
                diagID << std::setw(16) << std::left
                   << M_itest << " # itest\n";
            else
            {
                int wim_itest = vm["wim.itest"].template as<int>();
                int wim_jtest = vm["wim.jtest"].template as<int>();
                diagID << std::setw(16) << std::left
                   << wim_itest << " # itest\n";
                diagID << std::setw(16) << std::left
                   << wim_jtest << " # jtest\n";
            }
            diagID << std::setw(16) << std::left
               << wim_ice.mask[M_itest] << " # ICE_MASK\n";
       }
       else
       {
         std::cout << "Cannot open " << diagfile  << "\n";
         std::cerr << "error: open file " << diagfile << " for output failed!" <<"\n";
         std::abort();
       }
       diagID.close();
    }//end M_dump_diag

    dom = 2*PI*(freq_vec[nwavefreq-1]-freq_vec[0])/(nwavefreq-1);

    if (M_steady)
    {
        for (int fq = 0; fq < nwavefreq; fq++)
            for (int nth = 0; nth < nwavedirn; nth++)
            {
                value_type adv_dir = (-PI/180.)*(wavedir[nth]+90.);

                if (std::cos(adv_dir) >= 0.)
#pragma omp parallel for num_threads(max_threads) collapse(1)
                    for (int i = 0; i < M_num_elements; i++)
                        if (M_steady_mask[i] > 0.5)
                            sdf_dir[fq][nth][i] = sdf_inc[fq][nth][i];
            }
    }//M_steady

    //calc mean floe size outside of frequency loop;
    dave.assign(M_num_elements,0.);

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {
        if (wim_ice.mask[i] > .5)
        {
            //std::cout<<"setting dave, Dmax="<<wim_ice.dfloe[i]<<"\n";
            if (wim_ice.dfloe[i] <dfloe_miz_thresh)
            {
                if ( fsdopt == "RG" )
                    this->floeScaling(wim_ice.dfloe[i],1,dave[i]);
                else if ( fsdopt == "PowerLawSmooth" )
                    this->floeScalingSmooth(wim_ice.dfloe[i],1,dave[i]);
                //std::cout<<"dave ("<<fsdopt<<") = "<<dave[i]<<"\n";
            }
            else
            {
                //just use uniform dist
                dave[i] = wim_ice.dfloe[i];
                //std::cout<<"dave (uniform) = "<<dave[i]<<"\n";
            }
        }


        if ( M_dump_diag && (i==M_itest) )
        {
            std::fstream diagID(diagfile, std::ios::out | std::ios::app);
            diagID << "\n# Ice info: pre-breaking\n";
            diagID << std::setw(16) << std::left
                << wim_ice.mask[i] << " # ice mask\n";
            diagID << std::setw(16) << std::left
                << wim_ice.conc[i] << " # conc\n";
            diagID << std::setw(16) << std::left
                << wim_ice.thick[i] << " # h, m\n";
            diagID << std::setw(16) << std::left
                << dave[i] << " # D_av, m\n";
            diagID << std::setw(16) << std::left
                << wim_ice.dfloe[i] << " # D_max, m\n";

            if (atten && (wim_ice.mask[i]>.5))
                //no attenuation outside ice
                diagID << "\n# period, s | atten_dim, m^{-1}| damp_dim, m^{-1}\n";

            diagID.close();
        }
    }//end spatial loop i - have dave

    for (int fq = 0; fq < nwavefreq; fq++)
    {
        std::fill( atten_dim.begin(), atten_dim.end(), 0. );
        std::fill( damp_dim.begin(), damp_dim.end(), 0. );

#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
        {
            value_type c1d;
            if ((wim_ice.mask[i] >.5) && (atten))
            {
                // floes per unit length
                c1d = wim_ice.conc[i]/dave[i];

                // scattering
                atten_dim[i] = atten_nond[i][fq]*c1d;

                // damping
                damp_dim[i] = 2*damping[i][fq]*wim_ice.conc[i];

                if ( M_dump_diag && (i==M_itest) )
                {
                   std::fstream diagID(diagfile, std::ios::out | std::ios::app);
                   diagID << vec_period[fq] << "   "
                          << atten_dim[i] << "   "
                          << damp_dim[i] << "\n";
                   diagID.close();
                }
            }//end of ice check
        }//end of spatial loop i


        // set group velocity and directional spectrum to pass into advAtten*()
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < M_num_elements; i++)
        {
            for (int nth = 0; nth < nwavedirn; nth++)
            {
                ag2d_eff_temp[i] = ag_eff[i][fq];//move inside loop to use OMP

#if 0
                sdf3d_dir_temp[i][nth] = sdf_dir[i][nth][fq];
                if (std::isnan(sdf_dir[i][nth][fq]))
                {
                    std::cout<<"found NaN in sdf_dir at i,nth,fq="<<i<<","<<nth<<","<<fq<<"\n";
                    std::cout<<"ag_eff="<<ag_eff[i][fq]<<"\n";
                    throw std::runtime_error("sdf_dir has NaN (before advection & attenuation)");
                }
#endif
            }
        }

        //advect all directions
        if (!M_wim_on_mesh)
            this->advectDirections(sdf_dir[fq],ag2d_eff_temp);
        else
            this->advectDirectionsMesh(sdf_dir[fq],ag2d_eff_temp);

        //do attenuation &/or scattering, and integrate over directions 
        if (scatmod == "dissipated")
            this->attenSimple(sdf_dir[fq], S_freq, taux_om, tauy_om,
                    stokes_drift_x_om, stokes_drift_y_om,ag2d_eff_temp);
        else if (scatmod == "isotropic")
            this->attenIsotropic(sdf_dir[fq], S_freq, taux_om, tauy_om,
                    stokes_drift_x_om, stokes_drift_y_om,ag2d_eff_temp);

        // update after application of advAtten*()
#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
        {
            ag_eff[i][fq] = ag2d_eff_temp[i];//is this necessary?

#if 0
            for (int nth = 0; nth < nwavedirn; nth++)
            {
                sdf_dir[i][nth][fq] = sdf3d_dir_temp[i][nth];
                if (std::isnan(sdf_dir[i][nth][fq]))
                {
                    std::cout<<"found NaN in sdf_dir at i,nth,fq="<<i<<","<<nth<<","<<fq<<"\n";
                    std::cout<<"ag_eff="<<ag_eff[i][fq]<<"\n";
                    throw std::runtime_error("sdf_dir has NaN (after advection & attenuation)");
                }
                // std::cout<<"AFTER: SDIR["<< i << "]= "<< sdf_dir[i][nth][fq] <<"\n";
            }
#endif

            //std::cout<<"AFTER: taux_om["<< i << "]= "<< taux_om[i] <<"\n";
#if 1
            if (std::isnan(S_freq[i]))
            {
                std::cout<<"found NaN in S_freq at i="<<i<<"\n";
                throw std::runtime_error("S_freq has NaN (after advection & attenuation)");
            }
            if (std::isnan(taux_om[i]))
            {
                std::cout<<"found NaN in taux_om at i="<<i<<"\n";
                throw std::runtime_error("taux_om has NaN (after advection & attenuation)");
            }
            if (std::isnan(tauy_om[i]))
            {
                std::cout<<"found NaN in tauy_om at i="<<i<<"\n";
                throw std::runtime_error("tauy_om has NaN (after advection & attenuation)");
            }
#endif
        }//update sdf_dir after advAtten*()

        // integrate stress and stokes drift densities over frequency
#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
        {
            value_type tmp1 = rhowtr*gravity*taux_om[i]/ap_eff[i][fq];
            tau_x[i] += wt_om[fq]*tmp1;

            tmp1 = rhowtr*gravity*tauy_om[i]/ap_eff[i][fq];
            tau_y[i] += wt_om[fq]*tmp1;

            //2*\omega*k*\int_0^\pi S(\omega,\theta)\cos(\theta)d\theta
            tmp1 = 2*(2*PI*freq_vec[fq])*(2*PI/wlng_ice[i][fq])*stokes_drift_x_om[i];
            stokes_drift_x[i] += wt_om[fq]*tmp1;

            //2*\omega*k*\int_0^\pi S(\omega,\theta)\sin(\theta)d\theta
            tmp1 = 2*(2*PI*freq_vec[fq])*(2*PI/wlng_ice[i][fq])*stokes_drift_y_om[i];
            stokes_drift_y[i] += wt_om[fq]*tmp1;
        }

        // integrals for breaking program_options
#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
        {
            // convert from water amp's to ice amp's
            value_type F     = disp_ratio[i][fq];
            value_type kicel = 2*PI/wlng_ice[i][fq];
            value_type om    = 2*PI*freq_vec[fq];

            if (std::isnan(S_freq[i]))
            {
                std::cout<<"found NaN in S_freq at i,fq="<<i<<","<<fq<<"\n";
                std::cout<<"F,kicel,om,wt_om="<<F<<","<<kicel<<","<<om<<","<<wt_om[fq]<<"\n";
                throw std::runtime_error("S_freq has NaN");
            }
            if (std::isnan(F))
            {
                std::cout<<"found NaN in disp_ratio at i,fq="<<i<<","<<fq<<"\n";
                throw std::runtime_error("disp_ratio has NaN");
            }


            // ======================================================
            // 0-th spectral moments
            // - take abs as small errors can make S_freq negative
            value_type tmp = wt_om[fq]*S_freq[i];

            // variance of displacement (water)
            mom0w[i] += std::abs(tmp);

            // variance of displacement (ice)
            mom0[i] += std::abs(tmp*std::pow(F,2.));
            // ======================================================


            // ======================================================
            // 2-nd spectral moments
            tmp = wt_om[fq]*std::pow(om,2.)*S_freq[i];

            // variance of speed (water)
            mom2w[i] += std::abs(tmp);

            // variance of speed (ice)
            mom2[i] += std::abs(tmp*std::pow(F,2.));
            // ======================================================


            // ======================================================
            // variance of strain
            if (wim_ice.mask[i] == 1.)
            {
                // strain conversion factor
                // = k^2*h/2*F
                tmp = F*std::pow(kicel,2.)*wim_ice.thick[i]/2.0;

                // strain density
                tmp = wt_om[fq]*S_freq[i]*std::pow(tmp,2.);
                var_strain[i] += std::abs(tmp);
            }
            // ======================================================
        }//end i loop
    }//end freq loop

    // value_type _min = *std::min_element(mom0w.begin(),mom0w.end());
    // value_type _max = *std::max_element(mom0w.begin(),mom0w.end());
    // std::cout<<"Min f= " << _min <<"\n";
    // std::cout<<"Max f= " << _max <<"\n";


    // for (int i = 0; i < M_num_elements; i++)
    //     std::cout << "VRT[" << i < "]= " << var_strain[i] <<"\n";

    std::fill( Tp.begin(), Tp.end(), 0. );

    if (ref_Hs_ice)
    {
#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
        {
            Hs[i] = 4*std::sqrt(mom0[i]);
            if (mom2[i] > 0.)
                Tp[i] = 2*PI*std::sqrt(mom0[i]/mom2[i]);
        }
    }
    else
    {
#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
        {
            Hs[i] = 4*std::sqrt(mom0w[i]);
            if (mom2w[i] > 0.)
                Tp[i] = 2*PI*std::sqrt(mom0w[i]/mom2w[i]);
        }
    }

    //update mwd
    this->calcMWD();

    if ( M_dump_diag )
    {
       int i =  M_itest;
       std::fstream diagID(diagfile, std::ios::out | std::ios::app);

       diagID << std::setw(16) << std::left
          << mom0w[i] << " # mom0w, m^2\n";
       diagID << std::setw(16) << std::left
          << mom2w[i] <<" # mom2w, m^2/s^2\n";
       diagID << std::setw(16) << std::left
          << mom0[i] <<" # mom0, m^2\n";
       diagID << std::setw(16) << std::left
          << mom2[i] <<" # mom2, m^2/s^2\n";
       diagID << std::setw(16) << std::left
          << Hs[i] <<" # Hs, m\n";
       diagID << std::setw(16) << std::left
          << Tp[i] <<" # Tp, s\n";
       diagID << std::setw(16) << std::left
          << mwd[i] <<" # mwd, deg\n";
       diagID << std::setw(16) << std::left
          << tau_x[i] <<" # tau_x, Pa\n";
       diagID << std::setw(16) << std::left
          << tau_y[i] <<" # tau_y, Pa\n";
       diagID.close();
    }


    if (!(M_steady) && !(breaking))
    {
       //check energy conservation
       auto temparray = Hs;
       std::for_each(
             temparray.begin(), temparray.end(),
             [&](value_type& f){ f *= f; });
       E_tot = std::accumulate(
             temparray.begin(), temparray.end(),0.);

        // std::fill( var_strain.data(), var_strain.data()+var_strain.num_elements(), 1. );
        // E_tot = std::accumulate(var_strain.data(), var_strain.data()+var_strain.num_elements(),0.);
        // std::cout<<"Sum= "<< E_tot <<"\n";
    }

    // finally do floe breaking

    //std::cout<<"max_threads= "<< max_threads <<"\n";

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {
        value_type E_tot, wlng_crest, Dc;
        value_type adv_dir, F, kicel, om, ommin, ommax, om1, lam1, lam2, c1d, tmp;
        int jcrest;
        bool break_criterion;

        //std::cout << "MASK[" << i << "," << j << "]= " << wim_ice.mask[i] << " and "<< mom0[i]  <<"\n";

        bool broken_this_timestep = false;
        if ((wim_ice.mask[i] == 1.) && (mom0[i] >= 0.))
        {
            BreakInfo breakinfo =
            {
                conc:       wim_ice.conc[i],
                thick:      wim_ice.thick[i],
                mom0:       mom0[i],
                mom2:       mom2[i],
                var_strain: var_strain[i],
                dfloe:      wim_ice.dfloe[i],
                broken:     false
            };

            this->doBreaking(breakinfo);

            broken_this_timestep = breakinfo.broken;
            if (broken_this_timestep)
            {
                wim_ice.dfloe[i]     = breakinfo.dfloe;
                wim_ice.broken[i]    = 1.;
            }

        }//end test for ice and waves

        //wim_ice.nfloes[i] = 0.;

        if (wim_ice.dfloe[i] > 0.)
            wim_ice.nfloes[i] = wim_ice.conc[i]/std::pow(wim_ice.dfloe[i],2.);

        if (M_dump_diag && (wim_ice.mask[i] > .5) && (i==M_itest))
        {
           //dump diagnostic even if no waves (but only if ice)
           std::fstream diagID(diagfile, std::ios::out | std::ios::app);
           diagID << "\n# Ice info: post-breaking\n";
           diagID << std::setw(16) << std::left
              << wlng_crest << " # peak wavelength, m\n";
           diagID << std::setw(16) << std::left
              << wim_ice.dfloe[i] << " # D_max, m\n";
           diagID << std::setw(16) << std::left
              << broken_this_timestep << " # broken_this_timestep, bool\n";
           diagID.close();
        }

        if (std::isnan(mom0[i]))
        {
            std::cout<<"found NaN in mom0 at i="<<i<<"\n";
            throw std::runtime_error("mom0 has NaN");
        }
    }//end spatial loop i


    if (break_on_mesh) // breaking on nextsim mesh
    {

        // =================================================================
        bool TEST_INTERP_MESH   = false;
            //set to true if want to save mesh quantities inside WIM for testing
            //TODO get this option working  

        std::vector<value_type> mesh_dfloe_old;
        if (TEST_INTERP_MESH)
        {
            mesh_dfloe_old = nextsim_ice.dfloe;//copy before breaking

            std::vector<std::vector<value_type>> vectors(3);
            std::vector<std::string> names(3);
            vectors[0]  = mom0;
            vectors[1]  = mom2;
            vectors[2]  = var_strain;
            names[0]    = "mom0";
            names[1]    = "mom2";
            names[2]    = "var_strain";
            this->testInterp("grid",t_step,vectors,names);
            std::cout<<"TEST_INTERP_MESH: after export of grid results\n";
        }
        // =================================================================

        // =================================================================
        std::cout<<"break_on_mesh: before interp grid to mesh\n";
        // do interpolation
        // - set input data
        value_type_vec_ptrs input_data = {&mom0,&mom2,&var_strain};

        // - set output data
        int Ne = nextsim_mesh.num_elements;
        value_type_vec mom0_mesh      (Ne,0.);
        value_type_vec mom2_mesh      (Ne,0.);
        value_type_vec var_strain_mesh(Ne,0.);
        value_type_vec_ptrs output_data = {&mom0_mesh,&mom2_mesh,&var_strain_mesh};

        // - call routine
        this->gridToPoints(output_data,input_data,
                nextsim_mesh.elements_x, nextsim_mesh.elements_y);
        std::cout<<"break_on_mesh: after interp grid to mesh\n";
        // =================================================================


        // =================================================================
        //do breaking
        for (int i=0; i<Ne; ++i)
        {
            // set inputs to doBreaking
            BreakInfo breakinfo =
            {
                conc:       nextsim_ice.conc[i],
                thick:      nextsim_ice.thick[i],
                mom0:       mom0_mesh[i],
                mom2:       mom2_mesh[i],
                var_strain: var_strain_mesh[i],
                dfloe:      nextsim_ice.dfloe[i],
                broken:     false
            };

            //do breaking
            this->doBreaking(breakinfo);

            //update mesh vectors
            if (breakinfo.broken)
            {
                nextsim_ice.dfloe[i]   = breakinfo.dfloe;
                nextsim_ice.broken[i]  = 1.;
            }

        }//finish loop over elements
        std::cout<<"break_on_mesh: after breaking\n";
        // =================================================================

        // =================================================================
        if (TEST_INTERP_MESH)
        {   
            std::vector<std::vector<value_type>> vectors(7);
            std::vector<std::string> names(7);
            vectors[0]  = nextsim_ice.conc;
            vectors[1]  = nextsim_ice.thick;
            vectors[2]  = nextsim_ice.dfloe;
            vectors[3]  = mom0_mesh;
            vectors[4]  = mom2_mesh;
            vectors[5]  = var_strain_mesh;
            vectors[6]  = mesh_dfloe_old;
            names[0]    = "Concentration";
            names[1]    = "Thickness";
            names[2]    = "Dfloe";
            names[3]    = "Mom0";
            names[4]    = "Mom2";
            names[5]    = "Var_strain";
            names[6]    = "Dfloe_old";
            this->testInterp("mesh",t_step,vectors,names);
            std::cout<<"TEST_INTERP_MESH: after export of mesh results\n";
        }
        // =================================================================
    }

#if 0
    double Hs_max=0;
    for (int i = 0; i < nx; i++)
    {
        int j = ny-1;
        Hs_max  = max(Hs_max,Hs[i*ny+j]);
    }
    std::cout<<"Hs_max (j=ny) = "<< Hs_max <<"\n";

    Hs_max=0;
    for (int i = 0; i < nx; i++)
        {
            int j = 0;
            Hs_max  = max(Hs_max,Hs[i*ny+j]);
        }
    std::cout<<"Hs_max (j=0) = "<< Hs_max <<"\n";
#endif

    std::cout<<"Hs_max= "<< *std::max_element(Hs.begin(), Hs.end()) <<"\n";

    double taux_min  = *std::min_element(tau_x.begin(), tau_x.end());
    std::cout<<"taux_min= "
             <<std::setprecision(10)<< taux_min <<"\n";
    double taux_max  = *std::max_element(tau_x.begin(), tau_x.end());
    std::cout<<"taux_max= "
             <<std::setprecision(10)<< taux_max <<"\n";

    // std::cout<<"------------------------------------------------------\n";
    // std::cout<<"dfloe_max= "<< *std::max_element(dfloe.data(), dfloe.data()+dfloe.num_elements()) <<"\n";
    // std::cout<<"dfloe_min= "<< *std::min_element(dfloe.data(), dfloe.data()+dfloe.num_elements()) <<"\n";
}//timeStep

template<typename T>
void WimDiscr<T>::setMesh(mesh_type const &movedmesh)
{
    M_time_mesh_set     = this->getModelTime();//used in check when ice fields are set on mesh
    nextsim_mesh_old    = nextsim_mesh;

    //update nextsim_mesh with moved mesh
    this->resetMesh(movedmesh);
}

template<typename T>
void WimDiscr<T>::setMesh(mesh_type const &mesh_in,value_type_vec const &um_in)
{
    auto movedmesh = mesh_in;
    movedmesh.move(um_in,1.);
    this->setMesh(movedmesh);
}//setMesh

template<typename T>
void WimDiscr<T>::setMesh(mesh_type const &mesh_in,value_type_vec const &um_in,BamgMesh *bamgmesh)
{
    //update nextsim_mesh with moved mesh
    auto movedmesh = mesh_in;
    movedmesh.move(um_in,1.);
    this->setMesh(movedmesh);

    // ================================================================================
    // this interface should be called after regridding (if M_wim_on_mesh)
    //
    // regridding procedure INSIDE NEXTSIM:
    // *BEFORE REGRID:
    // 1) um0 = wim.RelativeMeshDisplacement(M_mesh_old,M_UM_old);
    //    - relative displacement between moved mesh just before displacement
    //      & mesh at last call to WIM
    // 2) interp um0 -> new mesh (after regrid)
    //    - this gives um1
    // 3) wim.resetMesh(M_mesh,um1)
    // 4) interp sdf_dir,taux,tauy -> new mesh
    // 5) integrate sdf_dir? (puts Hs etc on new mesh)
    // ================================================================================

    // get relative displacement of nodes since last call
    // - M_UM may already be nonzero if regridding has happened
    // - it is reset to zero at end of wim.run() and at initialisation
    // - used to correct group velocity when waves are advected
    int Nn = nextsim_mesh.num_nodes;
    if (M_UM.size()==0)
        M_UM.assign(2*Nn,0.);

    for (int i=0;i<Nn;i++)
    {
        M_UM[i]    += nextsim_mesh.nodes_x[i]-nextsim_mesh_old.nodes_x[i];
        M_UM[i+Nn] += nextsim_mesh.nodes_y[i]-nextsim_mesh_old.nodes_y[i];
    }

    //caluculate the surface areas,
    //get the element connectivity from bamgmesh
    int Nels = nextsim_mesh.num_elements;
    nextsim_mesh.surface.assign(Nels,0);
    nextsim_mesh.element_connectivity.assign(3*Nels,0);
    for (int i=0;i<Nels;i++)
    {
        value_type_vec xnods(3);
        value_type_vec ynods(3);
        for (int k=0;k<3;k++)
        {
            int ind  = nextsim_mesh.index[3*i+k];
            xnods[k] = nextsim_mesh.nodes_x[ind];
            ynods[k] = nextsim_mesh.nodes_y[ind];

            nextsim_mesh.element_connectivity[3*i+k]
                = bamgmesh->ElementConnectivity[3*i+k]-1;//NB bamg indices go from 1 to Nels
        }
        nextsim_mesh.surface[i] = NextsimTools::measure(
                xnods[0],ynods[0],xnods[1],ynods[1],xnods[2],ynods[2]);
    }
    // ================================================================================

}//setMesh

template<typename T>
void WimDiscr<T>::resetMesh(mesh_type const &mesh,value_type_vec const &um_in)
{
    //move the mesh then set nextsim_mesh
    auto movedmesh = mesh;
    movedmesh.move(um_in,1.);
    this->resetMesh(movedmesh);
}//resetMesh()

template<typename T>
void WimDiscr<T>::resetMesh(mesh_type const &mesh_in)
{
    //sets the variable "nextsim_mesh"
    nextsim_mesh.initialised    = true;
    nextsim_mesh.num_nodes      = mesh_in.numNodes();
    nextsim_mesh.num_elements   = mesh_in.numTriangles();
    nextsim_mesh.index          = mesh_in.indexTr();
    nextsim_mesh.nodes_x        = mesh_in.coordX();
    nextsim_mesh.nodes_y        = mesh_in.coordY();
    nextsim_mesh.elements_x     = mesh_in.bcoordX();
    nextsim_mesh.elements_y     = mesh_in.bcoordY();

}//resetMesh

template<typename T>
typename WimDiscr<T>::value_type_vec
WimDiscr<T>::relativeMeshDisplacement(mesh_type const &mesh_in,value_type_vec const &um_in) const
{
    auto movedmesh = mesh_in;
    movedmesh.move(um_in,1);
    auto nodes_x = movedmesh.coordX();
    auto nodes_y = movedmesh.coordY();
    int Nn = nodes_x.size();

    if (!nextsim_mesh.initialised)
        throw std::runtime_error("relativeMeshDisplacement: nextsim_mesh not initialised yet");
    if (nextsim_mesh.num_nodes!=Nn)
        throw std::runtime_error("relativeMeshDisplacement: mesh_in and nextsim_mesh have different sizes");

    value_type_vec um(2*Nn);
    for (int i=0;i<Nn;i++)
    {
        um[i]    = nodes_x[i]-nextsim_mesh.nodes_x[i];
        um[i+Nn] = nodes_y[i]-nextsim_mesh.nodes_y[i];
    }

    return um;
}//relativeMeshDisplacement

template<typename T>
void WimDiscr<T>::setIceFields(
                          std::vector<value_type> const& m_conc,  // conc
                          std::vector<value_type> const& m_vol,   // ice vol or effective thickness (conc*thickness)
                          std::vector<value_type> const& m_nfloes,// Nfloes=conc/Dmax^2
                          bool pre_regrid)
{
    if(M_time_mesh_set != this->getModelTime())
        throw std::runtime_error("setIceFields: setting ice without setting mesh first");

    // pre-regrid options:
    if(pre_regrid&&M_wim_on_mesh)
        //do nothing
        return;

    //if here, not (M_wim_on_mesh && pre_regrid)
    M_initialised_ice = true;
    if (pre_regrid)//not M_wim_on_mesh
    {
        //interp from mesh to grid
        nextsim_ice.conc    = m_conc;
        nextsim_ice.vol     = m_vol;
        nextsim_ice.nfloes  = m_nfloes;
        this->transformIce(nextsim_ice);
#if 0
        value_type_vec v1 = m_conc;
        value_type_vec v2 = m_vol;
        value_type_vec v3 = m_nfloes;
        value_type_vec_ptrs input_data = {&v1,&v2,&v3};
        icec  .assign(M_num_elements,0);
        nfloes.assign(M_num_elements,0);
        value_type_vec icevol(M_num_elements,0);
        value_type_vec_ptrs output_data = {&icec,&icevol,&nfloes};
#endif
        value_type_vec_ptrs input_data = {&(nextsim_ice.conc),&(nextsim_ice.vol),&(nextsim_ice.nfloes)};
        value_type_vec_ptrs output_data = {&(wim_ice.conc),&(wim_ice.vol),&(wim_ice.nfloes)};
        this->meshToGrid(output_data,input_data);
        this->transformIce(wim_ice);
#if 0
        ice_mask    = this->transformIce(iceh,icec,icevol,nfloes);
        dfloe   = this->nfloesToDfloe(nfloes,icec);
#endif
#if 1
        //test interp
        std::cout<<"setIceFields (pre_regrid): check ice inputs to WIM\n";
        std::cout<<"min conc   grid = "<< *std::min_element((wim_ice.conc).begin()  ,(wim_ice.conc).end()   )<<"\n";
        std::cout<<"max conc   grid = "<< *std::max_element((wim_ice.conc).begin()  ,(wim_ice.conc).end()   )<<"\n";
        std::cout<<"min thick  grid = "<< *std::min_element((wim_ice.thick).begin() ,(wim_ice.thick).end()  )<<"\n";
        std::cout<<"max thick  grid = "<< *std::max_element((wim_ice.thick).begin() ,(wim_ice.thick).end()  )<<"\n";
        std::cout<<"min dfloe  grid = "<< *std::min_element((wim_ice.dfloe).begin() ,(wim_ice.dfloe).end()  )<<"\n";
        std::cout<<"max dfloe  grid = "<< *std::max_element((wim_ice.dfloe).begin() ,(wim_ice.dfloe).end()  )<<"\n";
        std::cout<<"min Nfloes grid = "<< *std::min_element((wim_ice.nfloes).begin(),(wim_ice.nfloes).end() )<<"\n";
        std::cout<<"max Nfloes grid = "<< *std::max_element((wim_ice.nfloes).begin(),(wim_ice.nfloes).end() )<<"\n";
#endif
    }
    // end of pre-regrid options
    // ================================================================

    // ================================================================
    // post-regrid options:
    else if (M_wim_on_mesh)
    {
        //ice fields already where we need them
        wim_ice.conc    = m_conc;
        wim_ice.vol     = m_vol;
        wim_ice.nfloes  = m_nfloes;
        this->transformIce(wim_ice);
#if 0
        icec        = m_conc;
        nfloes      = m_nfloes;
        auto icevol = m_vol;
        ice_mask    = this->transformIce(iceh,icec,icevol,nfloes);
        dfloe       = this->nfloesToDfloe(nfloes,icec);
#endif
    }
    else if (break_on_mesh)
    {
        nextsim_ice.conc    = m_conc;
        nextsim_ice.vol     = m_vol;
        nextsim_ice.nfloes  = m_nfloes;
        this->transformIce(nextsim_ice);
#if 0
        auto tmp_nfloes = m_nfloes;
        auto tmp_vol    = m_vol;
        mesh_conc       = m_conc;
        auto tmp_mask   = this->transformIce(mesh_thick,mesh_conc,tmp_vol,tmp_nfloes);
        mesh_dfloe      = this->nfloesToDfloe(tmp_nfloes,mesh_conc);

        //need extra field to return to nextsim
        mesh_broken.assign(nextsim_mesh.num_elements,0);
        std::fill(mesh_broken.begin(),mesh_broken.end(),false);
#endif
    }
    // end of post-regrid options
    // ================================================================
}//setIceFields()

template<typename T>
void WimDiscr<T>::transformIce(IceInfo &ice_info)
{
    //checks for too low concentration,
    //calculates true thickness from ice volume
    int Nel = (ice_info.conc).size();
    (ice_info.thick).assign(Nel,0.);
    (ice_info.dfloe).assign(Nel,0.);
    (ice_info.mask).assign(Nel,1.);
    (ice_info.broken).assign(Nel,0.);
    for (int i=0;i<Nel;i++)
    {
        if (ice_info.conc[i]>=vm["wim.cicemin"].template as<double>())
            ice_info.thick[i]  = ice_info.vol[i]/ice_info.conc[i];//convert to actual thickness
        else
        {
            ice_info.conc[i]   = 0;
            ice_info.vol[i]    = 0;
            ice_info.nfloes[i] = 0;
            ice_info.mask[i]   = 0.;
        }
        ice_info.dfloe[i]  = this->nfloesToDfloe(ice_info.nfloes[i],ice_info.conc[i]);

#if 1
        //check ranges of inputs
        if ((ice_info.thick[i]<0.)||(ice_info.thick[i]>50.))
        {
            std::cout<<"thickness out of range for i="<<i<<"\n";
            std::cout<<"vol,h,c="<<ice_info.vol[i]<<","<<ice_info.thick[i]<<","<<ice_info.conc[i]<<"\n";
            throw std::runtime_error("thickness out of range");
        }
        if ((ice_info.conc[i]<0.)||(ice_info.conc[i]>1.))
        {
            std::cout<<"conc out of range for i="<<i<<"\n";
            std::cout<<"c="<<ice_info.conc[i]<<"\n";
            throw std::runtime_error("conc out of range");
        }
#endif
    }
}//transformIce()


template<typename T>
void WimDiscr<T>::clearMeshFields()
{
    (nextsim_ice.conc  ).resize(0);
    (nextsim_ice.thick ).resize(0);
    (nextsim_ice.dfloe ).resize(0);
    (nextsim_ice.broken).resize(0);
}

template<typename T>
void WimDiscr<T>::gridToPoints(
        value_type_vec_ptrs &output_data,       //output data
        value_type_vec_ptrs const &input_data,  //input data
        value_type_vec &Rx, value_type_vec &Ry) //output locations
{
    int nb_var      = input_data.size();
    int Ninterp     = (*(input_data[0])).size();//get pointer, then get size
    int target_size = Rx.size();

    value_type_vec interp_in(Ninterp*nb_var);   //input to interp routine
    for (int i=0;i<Ninterp;i++)
        for (int p=0;p<nb_var;p++)
            interp_in[nb_var*i+p]   = (*(input_data[p]))[i];

    value_type* interp_out;

    //bool regular_source = (source_grid.grid_type == "RegularGrid")&&M_regular;
    if (M_regular)
    {
        std::cout<<"gridToPoints: InterpFromGridToMeshx\n";
        int interptype = BilinearInterpEnum;
        InterpFromGridToMeshx(interp_out,       //data (out)
                              &x_col[0], nx,    //x vector (source), length of x vector
                              &y_row[0], ny,    //y vector (source), length of y vector
                              &interp_in[0],    //data (in)
                              ny, nx,           //M,N: no of grid cells in y,x directions
                                                //(to determine if corners or centers of grid have been input)
                              nb_var,           //no of variables
                              &Rx[0],           // x vector (target) (NB already a reference)
                              &Ry[0],           // x vector (target) (NB already a reference)
                              target_size,0.,   //target_size,default value
                              interptype,       //interpolation type
                              true              //row_major (false = fortran/matlab order)         
                              );
    }
    else
    {
        std::cout<<"gridToPoints: InterpFromMeshToMesh2dx\n";
        InterpFromMeshToMesh2dx(&interp_out,                        //output data
                              &(M_wim_triangulation.index)[0],      //index 
                              &(M_wim_triangulation.nodes_x)[0],    //input location (x-coord)
                              &(M_wim_triangulation.nodes_y)[0],    //input location (y-coord)
                              M_wim_triangulation.num_nodes,        //num nodes
                              M_wim_triangulation.num_elements,     //num elements
                              &interp_in[0],                        //input data
                              M_wim_triangulation.num_nodes,        //num input locations
                              nb_var,                               //num input variables
                              &Rx[0],                               //output location (x-coord)
                              &Ry[0],                               //output location (y-coord)
                              Rx.size(),                            //num output locations
                              false);                               //use default value if outside mesh (use nearest)
    }
    

    //output
    for (int p=0;p<nb_var;p++)
    {
        output_data[p]->resize(target_size,0.);
        for (int i=0;i<target_size;i++)
            (*(output_data[p]))[i]  = interp_out[nb_var*i+p];
    }

    xDelete<value_type>(interp_out);

}//gridToPoints

template<typename T>
void WimDiscr<T>::meshToGrid(
        value_type_vec_ptrs &output_data,       //output data
        value_type_vec_ptrs const &input_data)  //input data
{
    int nb_var      = input_data.size();
    int Ninterp     = (*(input_data[0])).size();//get pointer, then get size

    value_type_vec interp_in(Ninterp*nb_var);   //input to interp routine
    for (int i=0;i<Ninterp;i++)
        for (int p=0;p<nb_var;p++)
            interp_in[nb_var*i+p]   = (*(input_data[p]))[i];

    value_type* interp_out;
    int target_size = X_array.size();

    //bool regular_source = (source_grid.grid_type == "RegularGrid")&&M_regular;
    if (M_regular)
    {
        value_type ymax = *std::max_element(Y_array.begin(),Y_array.end());
        std::cout<<"sim2wim: before interp mesh2grid\n";
        InterpFromMeshToGridx(interp_out,                   //output data (pointer)
                              &(nextsim_mesh.index)[0],     //mesh index (element->node map)
                              &(nextsim_mesh.nodes_x)[0],   //node positions (x-coords)
                              &(nextsim_mesh.nodes_y)[0],   //node positions (y-coords)
                              nextsim_mesh.num_nodes,       //num nodes
                              nextsim_mesh.num_elements,    //num elements
                              &interp_in[0],                //input data
                              Ninterp,                      //data length
                              nb_var,                       //no of variables
                              x0,                           //xmin (of grid elements' positions)
                              ymax,                         //ymax (of grid elements' positions)
                              dx,                           //grid dx
                              dy,                           //grid dy
                              nx,                           //grid nx
                              ny,                           //grid ny
                              0.);                          //default value (given to points outside mesh)
    }
    else
    {
        value_type_vec xdata = nextsim_mesh.nodes_x;
        value_type_vec ydata = nextsim_mesh.nodes_y;
        if(Ninterp==nextsim_mesh.num_nodes)
        {
            xdata = nextsim_mesh.elements_x;
            ydata = nextsim_mesh.elements_y;
        }

        std::cout<<"meshToGrid: InterpFromMeshToMesh2dx\n";
        InterpFromMeshToMesh2dx(&interp_out,             // output data
                              &(nextsim_mesh.index)[0],  // index 
                              &xdata[0],                 // input location (x-coord)
                              &ydata[0],                 // input location (y-coord)
                              nextsim_mesh.num_nodes,    // num nodes
                              nextsim_mesh.num_elements, // num elements
                              &interp_in[0],             // input data
                              Ninterp,                   // num input locations
                              nb_var,                    // num input variables
                              &X_array[0],               // output location (x-coord)
                              &Y_array[0],               // output location (y-coord)
                              target_size,               // num output locations
                              false);                    // use default value if outside mesh (use nearest)
    }
    

    //output
    for (int p=0;p<nb_var;p++)
    {
        output_data[p]->resize(target_size);
        for (int i=0;i<target_size;i++)
            (*(output_data[p]))[i]  = interp_out[nb_var*i+p];
    }

    xDelete<value_type>(interp_out);

}//meshToGrid

template<typename T>
void WimDiscr<T>::meshToPoints(
        value_type_vec_ptrs &output_data,       //output data
        value_type_vec_ptrs const &input_data,  //input data
        value_type_vec &Rx,                     //location of output data (x-coord)
        value_type_vec &Ry)                     //location of output data (y-coord)
{
    int nb_var      = input_data.size();
    int Ninterp     = (*(input_data[0])).size();//get pointer, then get size

    value_type_vec interp_in(Ninterp*nb_var);   //input to interp routine
    for (int i=0;i<Ninterp;i++)
        for (int p=0;p<nb_var;p++)
            interp_in[nb_var*i+p]   = (*(input_data[p]))[i];

    value_type* interp_out;
    int target_size = Rx.size();


    value_type_vec xdata = nextsim_mesh.nodes_x;
    value_type_vec ydata = nextsim_mesh.nodes_y;
    if(Ninterp==nextsim_mesh.num_nodes)
    {
        xdata = nextsim_mesh.elements_x;
        ydata = nextsim_mesh.elements_y;
    }

    std::cout<<"meshToPoints: InterpFromMeshToMesh2dx\n";
    InterpFromMeshToMesh2dx(&interp_out,             // output data
                          &(nextsim_mesh.index)[0],  // index 
                          &xdata[0],                 // input location (x-coord)
                          &ydata[0],                 // input location (y-coord)
                          nextsim_mesh.num_nodes,    // num nodes
                          nextsim_mesh.num_elements, // num elements
                          &interp_in[0],             // input data
                          Ninterp,                   // num input locations
                          nb_var,                    // num input variables
                          &Rx[0],                    // output location (x-coord)
                          &Ry[0],                    // output location (y-coord)
                          target_size,               // num output locations
                          false);                    // use default value if outside mesh (use nearest)
    

    //output
    for (int p=0;p<nb_var;p++)
    {
        output_data[p]->resize(target_size);
        for (int i=0;i<target_size;i++)
            (*(output_data[p]))[i]  = interp_out[nb_var*i+p];
    }

    xDelete<value_type>(interp_out);

}//meshToPoints


template<typename T>
typename WimDiscr<T>::unord_map_vecs_type
WimDiscr<T>::returnFieldsNodes(std::vector<std::string> const & fields,
        mesh_type const &mesh_in,value_type_vec const &um_in)
{
    auto movedmesh  = mesh_in;
    movedmesh.move(um_in,1.);
    auto xnod = movedmesh.coordX();
    auto ynod = movedmesh.coordY();
    return this->returnFieldsNodes(fields,xnod,ynod);
}

template<typename T>
typename WimDiscr<T>::unord_map_vecs_type
WimDiscr<T>::returnFieldsElements(std::vector<std::string> const &fields,
        mesh_type const &mesh_in,value_type_vec const &um_in)
{
    auto movedmesh  = mesh_in;
    movedmesh.move(um_in,1.);
    auto xel  = movedmesh.bcoordX();
    auto yel  = movedmesh.bcoordY();
    return this->returnFieldsElements(fields,xel,yel);
}

template<typename T>
typename WimDiscr<T>::unord_map_vecs_type
WimDiscr<T>::returnFieldsNodes(std::vector<std::string> const &fields,
        value_type_vec &xnod, value_type_vec &ynod)
{
    // return fields to nodes of nextsim mesh
    // - usually to export diagnostic fields on nextsim mesh
    unord_map_vecs_type output_nodes;
    int Nnod = xnod.size();

    if (fields.size()==0)
        //nothing to do
        return output_nodes;
    else
        //initialise outputs
        for (auto it=fields.begin();it!=fields.end();it++)
        {
            value_type_vec tmp(2*Nnod,0.);
            output_nodes.emplace(*it,tmp);
        }

    // ==========================================================================================
    //nodes - vectors

    //input to and output from interpolation routines
    value_type_vec_ptrs input_nodes, out_nodes;

    //need to make some dummy variables since vector components handled individually
    value_type_vec tx_out,ty_out,sdfx_out,sdfy_out;
    for (typename unord_map_vecs_type::iterator it = output_nodes.begin(); it != output_nodes.end(); it++)
    {
        //get inputs
        if(it->first=="Stress_waves_ice")
        {
            input_nodes.push_back(&tau_x);
            input_nodes.push_back(&tau_y);
            out_nodes.push_back(&tx_out);
            out_nodes.push_back(&ty_out);
        }
        else if(it->first=="Stokes_drift")
        {
            input_nodes.push_back(&stokes_drift_x);
            input_nodes.push_back(&stokes_drift_y);
            out_nodes.push_back(&sdfx_out);
            out_nodes.push_back(&sdfy_out);
        }
        else
            throw std::runtime_error("returnFieldsNodes: unknown variable name - "+it->first+"\n");
    }
    // ==========================================================================================

    // ==========================================================================================
    // Do the interpolation
    if(M_wim_on_mesh)
        //from elements of last mesh to nodes of the input mesh
        this->meshToPoints(out_nodes,input_nodes,xnod,ynod);
    else
        //from grid elements to mesh nodes
        this->gridToPoints(out_nodes,input_nodes,xnod,ynod);
    // ==========================================================================================


    // ==========================================================================================
    // assign the outputs
    for (typename unord_map_vecs_type::iterator it = output_nodes.begin(); it != output_nodes.end(); it++)
    {
        if(it->first=="Stress_waves_ice")
            for (int i=0;i<Nnod;i++)
            {
                (it->second)[i]         = tx_out[i];
                (it->second)[i+Nnod]    = ty_out[i];
            }
        else if(it->first=="Stokes_drift")
            for (int i=0;i<Nnod;i++)
            {
                (it->second)[i]         = sdfx_out[i];
                (it->second)[i+Nnod]    = sdfy_out[i];
            }
    }
    // ==========================================================================================

    return output_nodes;
}//returnFieldsNodes

template<typename T>
typename WimDiscr<T>::unord_map_vecs_type
WimDiscr<T>::returnFieldsElements(std::vector<std::string> const &fields,
        value_type_vec &xel, value_type_vec &yel)
{
    // return fields to elements of nextsim_mesh
    // - usually to export diagnostic fields on nextsim mesh
    unord_map_vecs_type output_els;
    int Nels = xel.size();

    if (fields.size()==0)
        //do nothing
        return output_els;
    else
        //initialise outputs
        for (auto it=fields.begin();it!=fields.end();it++)
        {
            value_type_vec tmp(Nels,0.);
            output_els.emplace(*it,tmp);
        }


    // ==========================================================================================
    //elements - scalars

    //input to and output from interpolation routines
    value_type_vec_ptrs input_els, out_els;
    for (typename unord_map_vecs_type::iterator it = output_els.begin(); it != output_els.end(); it++)
    {
        //get inputs
        if(it->first=="Hs")
        {
            if (M_wim_on_mesh)
                it->second = Hs;
            else
            {
                input_els.push_back(&Hs);
                out_els.push_back(&(it->second));
            }
        }
        else if(it->first=="Tp")
        {
            if (M_wim_on_mesh)
                it->second = Tp;
            else
            {
                input_els.push_back(&Tp);
                out_els.push_back(&(it->second));
            }
        }
        else if(it->first=="MWD")
        {
            if (M_wim_on_mesh)
                it->second = mwd;
            else
            {
                input_els.push_back(&mwd);
                out_els.push_back(&(it->second));
            }
        }
        else
            throw std::runtime_error("returnFieldsElements): unknown variable name - "+it->first+"\n");
    }
    // ==========================================================================================

    //interp to elements if necessary
    if (!M_wim_on_mesh)
        this->gridToPoints(out_els,input_els,xel,yel);//out_els already points to output_els

    return output_els;
}//returnFieldsElements


template<typename T>
void WimDiscr<T>::returnWaveStress(value_type_vec &M_tau,mesh_type const &mesh_in,value_type_vec const &um_in)
{
    auto movedmesh  = mesh_in;
    movedmesh.move(um_in,1.);
    auto xnod = movedmesh.coordX();
    auto ynod = movedmesh.coordY();
    this->returnWaveStress(M_tau,xnod,ynod);
}

template<typename T>
void WimDiscr<T>::returnWaveStress(value_type_vec &M_tau,value_type_vec &xnod,value_type_vec &ynod)
{
    //return wave stress on nodes of nextsim mesh

    int Nnod = xnod.size();

    //nodes - vectors
    value_type_vec tx_out,ty_out;
    value_type_vec_ptrs input_nodes = {&tau_x,&tau_y};
    value_type_vec_ptrs out_nodes   = {&tx_out,&ty_out};

    //interp to nodes
    if(M_wim_on_mesh)
        this->meshToPoints(out_nodes,input_nodes,xnod,ynod);
    else
        this->gridToPoints(out_nodes,input_nodes,xnod,ynod);

    M_tau.resize(2*Nnod,0.);
    for (int i=0;i<Nnod;i++)
    {
        M_tau[i]        = tx_out[i];
        M_tau[i+Nnod]   = ty_out[i];
    }

}//returnWaveStress


template<typename T>
void WimDiscr<T>::doBreaking(BreakInfo const& breakinfo)
{

    bool break_criterion = (breakinfo.conc > 0) /*ice present*/ && ((2*breakinfo.var_strain) > std::pow(epsc, 2.)) /*big enough waves*/;

    if (break_criterion)
    {
        double params[5];
        params[0] = young;
        params[1] = gravity;
        params[2] = rhowtr;
        params[3] = rhoice;
        params[4] = poisson;

        double outputs[8];

        value_type om = std::sqrt(breakinfo.mom2/breakinfo.mom0);
        value_type guess = std::pow(om,2.)/gravity;

        RTparam_outer(outputs,breakinfo.thick,double(om),double(drag_rp),double(guess),params);

        value_type kice = outputs[1];
        value_type lam = 2*PI/kice;

        if (lam < (2*breakinfo.dfloe))
        {
            breakinfo.dfloe     = std::max<value_type>(dmin,lam/2.);
            breakinfo.broken    = true;
        }
    }
}//doBreaking


template<typename T>
typename WimDiscr<T>::value_type
WimDiscr<T>::dfloeToNfloes(value_type const& dfloe_in,
                           value_type const& conc_in)
{
    value_type nfloes_out   = 0.;

    if ( (dfloe_in>0) &&(conc_in >= cice_min) )
    {
        //conc high enough & dfloe OK
        nfloes_out = conc_in/std::pow(dfloe_in,2.);
    }

    return nfloes_out;
}//dfloesToNfloes


template<typename T>
typename WimDiscr<T>::value_type_vec
WimDiscr<T>::dfloeToNfloes(value_type_vec const& dfloe_in,
                           value_type_vec const& conc_in)
{
    int N   = conc_in.size();
    value_type_vec nfloes_out(N);

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i=0;i<N;i++)
        nfloes_out[i]   = this->dfloeToNfloes(dfloe_in[i],conc_in[i]);

    return nfloes_out;
}//dfloesToNfloes


template<typename T>
typename WimDiscr<T>::value_type
WimDiscr<T>::nfloesToDfloe(value_type const& nfloes_in,
                           value_type const& conc_in)
{
        value_type dfloe_out  = 0.;
        if ( (nfloes_in>0)
                &&(conc_in >= vm["wim.cicemin"].template as<double>()) )
        {
            //conc high enough & Nfloes OK
            dfloe_out   = std::sqrt(conc_in/nfloes_in);
        }

        //dfloe shouldn't get too big
        if ( dfloe_out>=vm["wim.dfloepackthresh"].template as<double>() )
            dfloe_out = vm["wim.dfloepackinit"].template as<double>();

    return dfloe_out;
}//nfloesToDfloe


template<typename T>
typename WimDiscr<T>::value_type_vec
WimDiscr<T>::nfloesToDfloe(value_type_vec const& nfloes_in,
                           value_type_vec const& conc_in)
{
    int N   = conc_in.size();
    value_type_vec dfloe_out(N);

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i=0;i<N;i++)
        dfloe_out[i] = this->nfloesToDfloe(nfloes_in[i],conc_in[i]);

    return dfloe_out;
}//nfloesToDfloe


template<typename T>
void WimDiscr<T>::getFsdMesh(value_type_vec &nfloes_out,value_type_vec &dfloe_out,value_type_vec &broken)
{
    if (M_wim_on_mesh)
    {
        nfloes_out  = wim_ice.nfloes;
        dfloe_out   = wim_ice.dfloe;
        broken      = wim_ice.broken;
    }
    else if (break_on_mesh)
    {
        broken      = nextsim_ice.broken;
        dfloe_out   = nextsim_ice.dfloe;
        nfloes_out  = dfloeToNfloes(nextsim_ice.dfloe,nextsim_ice.conc);
    }
    else
        throw std::runtime_error("getFsdMesh: using wrong interface");
}//getFsdMesh

template<typename T>
void WimDiscr<T>::getFsdMesh(value_type_vec &nfloes_out,value_type_vec &dfloe_out,value_type_vec &broken,
        value_type_vec const & conc_tot, mesh_type const &mesh_in,value_type_vec const &um_in)
{
    if((M_wim_on_mesh)||(break_on_mesh))
        throw std::runtime_error("getFsdMesh: using wrong interface");

    // set nextsim_mesh (need to know where to interpolate to)
    // - NB set in FiniteElement::wimPreRegrid(),
    // but this function is called from FiniteElement::wimPostRegrid(),
    // and mesh could have changed due to regridding
    this->resetMesh(mesh_in,um_in);

    //do interpolation
    value_type_vec cinterp;//interp conc as well, to correct for interpolation error
    value_type_vec_ptrs input_data = {&(wim_ice.nfloes),&(wim_ice.conc),&(wim_ice.broken)};
    value_type_vec_ptrs output_data = {&nfloes_out,&cinterp,&broken};
    this->gridToPoints(output_data,input_data,nextsim_mesh.elements_x,nextsim_mesh.elements_y);

    int Nel = nextsim_mesh.num_elements;
    dfloe_out.assign(Nel,0.);
    for (int i=0;i<Nel;i++)
    {
        if(conc_tot[i]>=vm["wim.cicemin"].template as<double>())
        {
            if(cinterp[i]>0.)
            {
                // if conc_tot is high enough:
                // nfloes  has      been interpolated: old mesh->grid->new mesh
                // cinterp has      been interpolated: old mesh->grid->new mesh
                // ctot    may have been interpolated: old mesh->new mesh (if no regridding, old = new & nextsim_ice.conc = ctot)
                // cinterp and ctot should be the same, apart from errors in interpolation,
                // so hopefully the errors in nfloes can be estimated from this process
                value_type cfac = conc_tot[i]/cinterp[i];
                nfloes_out[i]   = cfac*nfloes_out[i];
            }

            //else keep nfloes_out the same
            dfloe_out[i]    = nfloesToDfloe(nfloes_out[i],conc_tot[i]);
            broken[i]       = std::round(broken[i]);
        }
        else
        {
            nfloes_out[i]   = 0.;
            broken[i]       = 0.;
        }
    }
}//getFsdMesh


template<typename T>
void WimDiscr<T>::run()
{

    //====================================================
    //set ice & wave conditions
    if (!M_initialised_ice)
    {
        M_initialised_ice   = true;
        this->idealIceFields(0.7);
    }
    if (!M_initialised_waves)
    {
        this->idealWaveFields(0.8);
    }
    // ===================================================


    // set incident wave spectrum,
    // attenuation coefficients and wave speeds/lengths
    this->update();
    if(M_steady&&(M_cpt==0))
    {
        //at first time step, set these, and top up initial area
        M_steady_mask = wave_mask;
        sdf_inc       = sdf_dir;
    }


    int lcpt  = 0;//local counter

    std::cout << "-----------------------Simulation started at "<< current_time_local() <<"\n";

    //init_time_str is human readable time eg "2015-01-01 00:00:00"
    //init_time is "formal" time format eg "20150101T000000Z"

    std::string init_time = ptime(init_time_str);

    value_type t_in  = this->getModelTime(lcpt);//model time of current call to wim (relative to when it was first called)
    std::string call_time = ptime(init_time_str,t_in);
    std::cout<<"---------------INITIAL TIME: "<< init_time <<"\n";
    std::cout<<"---------------CALLING TIME: "<< call_time <<"\n";


    std::cout<<"Running starts\n";
    chrono.restart();

    std::cout<<"duration= "<< duration <<"\n";
    std::cout<<"amax= "<< amax <<"\n";
    std::cout<<"M_timestep= "<< M_timestep <<"\n";
    std::cout<<"nt= "<< nt <<"\n";

    if (vm["wim.checkinit"].template as<bool>())
        this->exportResults("init",t_in);

    if (vm["wim.checkincwaves"].template as<bool>()
        &&(swh_in_array.size()>0))
        this->exportResults("incwaves",t_in);

#if 1
    if (swh_in_array.size()>0)
    {
        //test inc waves
        value_type _min = *std::min_element(swh_in_array.begin(),swh_in_array.end());
        value_type _max = *std::max_element(swh_in_array.begin(),swh_in_array.end());
        std::cout<<"Min swh in = " << _min <<"\n";
        std::cout<<"Max swh in = " << _max <<"\n";
        //
        _min = *std::min_element(mwp_in_array.begin(),mwp_in_array.end());
        _max = *std::max_element(mwp_in_array.begin(),mwp_in_array.end());
        std::cout<<"Min mwp in = " << _min <<"\n";
        std::cout<<"Max mwp in = " << _max <<"\n";
        //
        _min = *std::min_element(mwd_in_array.begin(),mwd_in_array.end());
        _max = *std::max_element(mwd_in_array.begin(),mwd_in_array.end());
        std::cout<<"Min mwd in = " << _min <<"\n";
        std::cout<<"Max mwd in = " << _max <<"\n";
    }
#endif
    std::cout<<"checked incwaves...\n";

    while (lcpt < nt)
    {
        std::cout <<  ":[WIM2D TIME STEP]^"<< lcpt+1
           <<" (out of "<<nt<<")"<<"\n";

        //export progress results
        bool export_now = false;
        int dump_freq   = vm["wim.dumpfreq"].template as<int>();//output "prog" and diagnostic text files after this no of time steps
        if (dump_freq>0)
            export_now = !(M_cpt % dump_freq);
        bool exportProg = export_now && (vm["wim.checkprog"].template as<bool>());
        if ( exportProg )
            this->exportResults("prog",this->getModelTime(lcpt));

        //integrate model
        M_dump_diag = export_now && (M_itest>0); 
        this->timeStep();

        ++lcpt;//local counter incremented in wim.run()
        ++M_cpt;//global counter incremented in wim.run()
    }

    M_update_time   = this->getModelTime(lcpt);//next time run is called, lcpt will be relative to this

    if (vm["wim.checkfinal"].template as<bool>())
       this->exportResults("final",M_update_time);

    // save diagnostic file
    if (vm["wim.savelog"].template as<bool>())
       this->saveLog(t_in);

    if (M_wim_on_mesh)
        //set M_UM to zero again
        std::fill(M_UM.begin(),M_UM.end(),0.);

    std::cout<<"Running done in "<< chrono.elapsed() <<"s\n";

    std::cout << "-----------------------Simulation completed at "<< current_time_local() <<"\n";
}//run

template<typename T>
void WimDiscr<T>::floeScaling(
      value_type const& dmax, int const& moment, value_type& dave_)
{
    value_type nm,nm1,dm,nsum,ndsum,r;

    int mm;
    value_type ffac     = fragility*std::pow(xi,2);

    dave_ = std::max(std::pow(dmin,moment),std::pow(dmax,moment));
    if (dmax>=xi*dmin)
    {
       //determine number of breaks
       r    = dmax/dmin;
       mm   = 0;
       while (r >= xi)
       {
          //r<2,mm=0 => doesn't break in 2
          //dave stays at dmax^moment;
          r  = r/xi;
          ++mm;
       }

       if (mm > 0)
       {
          nm1   = 1.;
          dm    = dmax; //floe length
          nsum  = 0.;   //eventually \sum_m=0^mm.n_m
          ndsum = 0.;   //eventually \sum_m=0^mm.n_m.d_m^moment

          for (int m=0; m<mm; ++m)
          {
              //no of floes of length dm;
              nm     = nm1*(1-fragility);
              nsum  += nm;
              ndsum += nm*std::pow(dm,moment);
              //std::cout<<"nsum,dm: "<<nsum<<" , "<<dm<<"\n";

              nm1   *= ffac;
              dm    /= xi;
          }

          //m=mm:
          nsum   += nm1;
          ndsum  += nm1*std::pow(dm,moment);
          dave_   = ndsum/nsum;
          //std::cout<<"nsum,dm: "<<nsum<<" , "<<dm<<"\n";
       }
    }
}//floeScaling

template<typename T>
void WimDiscr<T>::floeScalingSmooth(
      value_type const& dmax,int const& moment, value_type& dave_)
{

    value_type fsd_exp,b,A;

    fsd_exp = 2+log(fragility)/log(xi);//power law exponent: P(d>D)=(D_min/D)^fsd_exp;
    b       = moment-fsd_exp;

    // calculate <D^moment> from Dmax
    // - uniform dist for larger floes
    dave_ = std::pow(dmax,moment);
    //std::cout<<"dave (1) = "<<dave_<<"\n";

    if (dmax<=dmin)
    {
        // small floes
        dave_ = std::pow(dmin,moment);
        //std::cout<<"dave (2) = "<<dave_<<"\n";
    }
    else
    {
        // bigger floes
        A     = (fsd_exp*std::exp(fsd_exp*(std::log(dmin)+std::log(dmax))));
        A     = A/(std::exp(fsd_exp*std::log(dmax))-std::exp(fsd_exp*std::log(dmin)));
        dave_ = -(A/b)*(std::exp(b*std::log(dmin))-exp(b*std::log(dmax)));
        //std::cout<<"dave (3) = "<<dave_<<"\n";
    }
}//floeScalingSmooth

template<typename T>
void WimDiscr<T>::advectDirections(value_type_vec2d& Sdir,value_type_vec const& ag2d_eff)
//void WimDiscr<T>::advectDirections(array2_type& Sdir,value_type_vec const& ag2d_eff)
{

	for (int nth = 0; nth < nwavedirn; nth++)
    {
        value_type adv_dir      = -PI*(90.0+wavedir[nth])/180.0;
        value_type_vec uwave   = ag2d_eff;
        value_type_vec vwave   = ag2d_eff;

        //set wave speeds
        //TODO if M_wim_on_mesh, subtract average mesh velocity
        std::for_each(uwave.begin(), uwave.end(), [&](value_type& f){ f *= std::cos(adv_dir); });
        if (M_advdim == 2)
            std::for_each(vwave.begin(), vwave.end(), [&](value_type& f){ f *= std::sin(adv_dir); });

#if 0
        // copy from 3D input array to 2D temporary array
        value_type_vec temp(M_num_elements,0.);
#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
            temp[i] = Sdir[nth][i];
#endif

        //do advection
        //TODO if M_wim_on_mesh call NextsimTools::advect here
        this->waveAdvWeno(Sdir[nth],uwave,vwave);

#if 0
        // copy from 2D temporary array back to 3D input array
#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
            Sdir[i][nth] = temp[i];
#endif
    }//advection of each direction done
}//advectDirections()

template<typename T>
void WimDiscr<T>::advectDirectionsMesh(value_type_vec2d& Sdir,value_type_vec & ag2d_eff)
//void WimDiscr<T>::advectDirectionsMesh(array2_type& Sdir,value_type_vec& ag2d_eff)
{

    //interpolate ag2d_eff from elements to nodes
    int Nnod = nextsim_mesh.num_nodes;
    value_type_vec agnod(Nnod,0.);
    value_type_vec_ptrs input  = {&ag2d_eff};
    value_type_vec_ptrs output = {&agnod};
    this->meshToPoints(output,input,nextsim_mesh.nodes_x,nextsim_mesh.nodes_y);
    value_type* advect_out;
    int nb_var  = 1;                    //have to advect 1 vbl at a time
    std::vector<int> adv_method = {1};  //alternative (0) is do nothing

    //advect the directions
	for (int nth = 0; nth < nwavedirn; nth++)
    {
        value_type adv_dir = -PI*(90.0+wavedir[nth])/180.0;
        value_type_vec VC(2*Nnod,0.);

        // set wave speeds
        // - subtract average mesh velocity
        // (average over length of call to wim = "duration")
#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i=0;i<Nnod;i++)
        {
            VC[i]      = agnod[i]*std::cos(adv_dir)-M_UM[i]/duration;
            VC[i+Nnod] = agnod[i]*std::sin(adv_dir)-M_UM[i+Nnod]/duration;
        }

#if 0
        // copy from 3D input array to 2D temporary array
        value_type_vec advect_in(M_num_elements,0.);
#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
            advect_in[i] = Sdir[i][nth];
#endif

        //do advection
        //auto advect_in  = Sdir[nth];
        //NextsimTools::advect(&advect_out,&advect_in[0],&nextsim_mesh,
        NextsimTools::advect(&advect_out,&(Sdir[nth])[0],&nextsim_mesh,
            &VC[0],&adv_method[0],nb_var,M_timestep);

        // copy from 2D temporary array back to 3D input array
#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
            Sdir[nth][i] = advect_out[i];
    }//advection of each direction done

    xDelete<value_type>(advect_out);
}//advectDirectionsMesh()

template<typename T>
//void WimDiscr<T>::attenSimple(array2_type& Sdir, value_type_vec& Sfreq,
void WimDiscr<T>::attenSimple(value_type_vec2d& Sdir, value_type_vec& Sfreq,
        value_type_vec& taux_omega, value_type_vec& tauy_omega,
        value_type_vec& sdx_omega, value_type_vec& sdy_omega,
        value_type_vec const& ag2d_eff)
{

	std::vector<value_type> wt_theta(nwavedirn);
	value_type S_th, tmp, alp_dim, source;

    if (nwavedirn == 1)
        std::fill( wt_theta.begin(), wt_theta.end(), 1. );
    else
        std::fill( wt_theta.begin(), wt_theta.end(), 2*PI/(1.0*nwavedirn) );

    std::fill( Sfreq.begin()      ,Sfreq.end()      ,0. );
    std::fill( taux_omega.begin() ,taux_omega.end() ,0. );
    std::fill( tauy_omega.begin() ,tauy_omega.end() ,0. );
    std::fill( sdx_omega.begin()  ,sdx_omega.end()  ,0. );
    std::fill( sdy_omega.begin()  ,sdy_omega.end()  ,0. );

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {
        if (wim_ice.mask[i] > .5)
        {
            for (int nth = 0; nth < nwavedirn; nth++)
            {
                value_type adv_dir = -PI*(90.0+wavedir[nth])/180.0;
                value_type S_th = Sdir[nth][i];
                value_type alp_dim = atten_dim[i]+damp_dim[i];

                // stress calculation
                value_type source = -alp_dim*ag2d_eff[i]*S_th;
                value_type tmp = -std::cos(adv_dir)*wt_theta[nth]*source;
                taux_omega[i] += tmp;
                tmp = -std::sin(adv_dir)*wt_theta[nth]*source;
                tauy_omega[i] += tmp;
#if 0
                if (i==M_itest)
                {
                    std::cout<<"i,nth,adv_dir = "<<i<<","<<nth<<","<<adv_dir<<"\n";
                    std::cout<<"atten_dim,damp_dim = "<<atten_dim[i]<<","<<damp_dim[i]<<"\n";
                    std::cout<<"alp_dim,source,wt_theta = "<<alp_dim<<","<<source<<","<<wt_theta[nth]<<"\n";
                    std::cout<<"taux_omega = "<<taux_omega[i]<<"\n";
                    std::cout<<"tauy_omega = "<<tauy_omega[i]<<"\n";
                }
#endif

                // do attenuation
                Sdir[nth][i] = S_th*std::exp(-alp_dim*ag2d_eff[i]*M_timestep);

                //std::cout<<"tau_x["<< i << "]= "<< atten_dim[i] <<"\n";
            }//loop over directions
        }//end "if ice"

        // integrals to be done everywhere (not just in ice)
        for (int nth = 0; nth < nwavedirn; nth++)
        {
            //frequency spectrum
            Sfreq[i] += wt_theta[nth]*Sdir[nth][i];

            //Stoke's drift at surface (x)
            value_type adv_dir = -PI*(90.0+wavedir[nth])/180.0;
            value_type tmp     = std::cos(adv_dir)*wt_theta[nth]*Sdir[nth][i];
            sdx_omega[i] += tmp;

            //Stoke's drift at surface (y)
            tmp = std::sin(adv_dir)*wt_theta[nth]*Sdir[nth][i];
            sdy_omega[i] += tmp;
        }

#if 0
        if (i==M_itest)
        {
            std::cout<<"i = "<<i<<"\n";
            std::cout<<"Sfreq = "<<Sfreq[i]<<"\n";
            std::cout<<"taux_omega = "<<taux_omega[i]<<"\n";
            std::cout<<"tauy_omega = "<<tauy_omega[i]<<"\n";
        }
#endif
    }
}//attenSimple


template<typename T>
//void WimDiscr<T>::attenIsotropic(array2_type& Sdir, value_type_vec& Sfreq,
void WimDiscr<T>::attenIsotropic(value_type_vec2d& Sdir, value_type_vec& Sfreq,
        value_type_vec& taux_omega, value_type_vec& tauy_omega,
        value_type_vec& sdx_omega, value_type_vec& sdy_omega,
        value_type_vec const& ag2d_eff)
{
    std::vector<value_type> nvec(nwavedirn);
	std::vector<value_type> K_fou(nwavedirn), S_th(nwavedirn), theta_vec(nwavedirn), wt_theta(nwavedirn);
	std::vector<value_type> tmp1(nwavedirn), evals_x(nwavedirn);
	value_type tmp, alp_dim, source;

	std::vector<std::complex<value_type> > S_fou(nwavedirn);
	std::complex<value_type> zi, src_fou_p1, src_fou_m1;

	std::vector<value_type> S_cos(ncs), S_sin(ncs);
	value_type cg, q_scat, q_abs, q_tot, src_cos_1, src_sin_1;
    int n, jp1, jm1;

	zi = std::complex<value_type>(0.,1.);

	std::fill( wt_theta.begin(), wt_theta.end(), 2*PI/(1.0*nwavedirn) );//TODO ok for ndir==1?

    std::fill( Sfreq.begin()      ,Sfreq.end()      ,0. );
    std::fill( taux_omega.begin() ,taux_omega.end() ,0. );
    std::fill( tauy_omega.begin() ,tauy_omega.end() ,0. );
    std::fill( sdx_omega.begin()  ,sdx_omega.end()  ,0. );
    std::fill( sdy_omega.begin()  ,sdy_omega.end()  ,0. );


	for (int i = 0; i < M_num_elements; i++)
    {
        for (int nth = 0; nth < nwavedirn; nth++)
            S_th[nth] = Sdir[nth][i];

        std::fill( S_fou.begin(), S_fou.end(), zi );

        //S_fou[0] = std::complex<value_type>( sum(wt_theta*S_th) );
        S_fou[0] = std::complex<value_type>( std::inner_product(wt_theta.begin(), wt_theta.end(), S_th.begin(), 0.) );

        if (wim_ice.mask[i] > 0.5 )
        {
            if (wim_ice.dfloe[i] < dfloe_pack_init)
            {
                q_scat = atten_dim[i];
                q_abs = damp_dim[i];
            }
            else
            {
                q_scat = 0;
                q_abs = atten_dim[i]+damp_dim[i];
            }

            q_tot = q_scat+q_abs;
            cg = ag2d_eff[i];

            std::fill(K_fou.begin(), K_fou.end(), 0.);
            K_fou[0] = q_scat;

            std::fill(evals_x.begin(), evals_x.end(), -q_tot);
            evals_x[0] = -q_abs;


            for (int nth = 0; nth < ncs; nth++)
            {
                std::vector<value_type> prodtmp = theta_vec;

                std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::cos((nth+1)*f); });

                std::transform(prodtmp.begin(), prodtmp.end(), S_th.begin(), prodtmp.begin(),
                              std::multiplies<value_type>());

                S_cos[nth] = std::inner_product(prodtmp.begin(), prodtmp.end(), wt_theta.begin(), 0.);

                prodtmp.clear();
                prodtmp.resize(theta_vec.size());
                prodtmp = theta_vec;
                std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::sin((nth+1)*f); });

                std::transform(prodtmp.begin(), prodtmp.end(), S_th.begin(), prodtmp.begin(),
                               std::multiplies<value_type>());

                S_sin[nth] = std::inner_product(prodtmp.begin(), prodtmp.end(), wt_theta.begin(), 0.);

                S_fou[nth+1] = std::complex<value_type>(S_cos[nth],S_sin[nth]);
                //S_fou[nwavedirn-nth] = std::complex<value_type>(S_cos[nth],-S_sin[nth]);

                if (nth != ncs-1)
                {
                    S_fou[nwavedirn-(nth+1)] = std::complex<value_type>(S_cos[nth],-S_sin[nth]);
                }

                // std::cout<<"nth= " << nth << ": and ncs+nth= "<< ncs+nth <<"\n";
                // std::cout<<"taux= "<< S_sin[nth] <<"\n";
            }

            // stresses
            jp1 = 1;
            jm1 = nwavedirn-1;

            src_fou_p1 = cg*(-q_tot*S_fou[jp1]+K_fou[jp1]*S_fou[jp1]);
            src_fou_m1 = cg*(-q_tot*S_fou[jm1]+K_fou[jm1]*S_fou[jm1]);
            src_cos_1 = std::real(std::complex<value_type>(0.5)*(src_fou_p1+src_fou_m1));
            src_sin_1 = std::real(-zi*std::complex<value_type>(0.5)*(src_fou_p1-src_fou_m1));

            taux_omega[i] = -src_cos_1;
            tauy_omega[i] = -src_sin_1;

            sdx_omega[i] = std::real(.5*(S_fou[jp1]+S_fou[jm1]));
            sdy_omega[i] = std::real(-.5*zi*(S_fou[jp1]-S_fou[jm1]));

            // if (i==1)
            // {
            //     std::cout<<"taux_omega= "<< taux_omega[i] <<"\n";
            //     std::cout<<"tauy_omega= "<< tauy_omega[i] <<"\n";
            // }

            std::vector<value_type> prodtmp = evals_x;
            std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::exp(cg*M_timestep*f); });
            std::transform(S_fou.begin(), S_fou.end(), prodtmp.begin(), S_fou.begin(),
                           std::multiplies<std::complex<value_type> >());

            std::vector<value_type> Sfoutempcos(S_fou.size());// = S_fou;
            std::vector<value_type> Sfoutempsin(S_fou.size());// = S_fou;

            for (int nth = 0; nth < nwavedirn; nth++)
            {
                // prodtmp = std::vector<value_type>(nvec.begin(), nvec.end()); // also works

                for (int ss=0; ss<S_fou.size(); ++ss)
                {
                    Sfoutempcos[ss] = std::real(S_fou[ss]);
                    Sfoutempsin[ss] = std::imag(S_fou[ss]);
                }

                prodtmp.clear();
                prodtmp.resize(nwavedirn);
                prodtmp = nvec;

                std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::cos(theta_vec[nth]*f); });

                std::transform(Sfoutempcos.begin(), Sfoutempcos.end(), prodtmp.begin(), Sfoutempcos.begin(),
                               std::multiplies<value_type>());


                // prodtmp = std::vector<value_type>(nvec.begin(), nvec.end()); // also works
                prodtmp.clear();
                prodtmp.resize(nwavedirn);
                prodtmp = nvec;

                std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::sin(theta_vec[nth]*f); });

                std::transform(Sfoutempsin.begin(), Sfoutempsin.end(), prodtmp.begin(), Sfoutempsin.begin(),
                               std::multiplies<value_type>());

                for (int kth = 0; kth < nwavedirn; kth++)
                {
                    tmp1[kth] = Sfoutempcos[kth]-Sfoutempsin[kth];
                }

                Sdir[nth][i] = std::accumulate(tmp1.begin(), tmp1.end(), 0.0)/(2*PI);
                S_th[nth]    = Sdir[nth][i];
            }
        }

        Sfreq[i] = std::real(S_fou[0]);
    }

#if 0
    // value_type _min = *std::min_element(Sdir.data(),Sdir.data() + Sdir.num_elements());
    // value_type _max = *std::max_element(Sdir.data(),Sdir.data() + Sdir.num_elements());
    // std::cout<<"Min OUT= " << _min <<"\n";
    // std::cout<<"Max OUT= " << _max <<"\n";


    value_type _min = *std::min_element(Sfreq.data(),Sfreq.data() + Sfreq.num_elements());
    value_type _max = *std::max_element(Sfreq.data(),Sfreq.data() + Sfreq.num_elements());
    std::cout<<"Min OUT= " << _min <<"\n";
    std::cout<<"Max OUT= " << _max <<"\n";
#endif
}//attenIsotropic

template<typename T>
void WimDiscr<T>::waveAdvWeno(value_type_vec& h, value_type_vec const& u, value_type_vec const& v)
{
    value_type_vec sao;
    value_type_vec u_pad, v_pad, scp2_pad, scp2i_pad, scuy_pad, scvx_pad, h_pad;
    int num_p_ext = nxext*nyext;
    hp.resize(num_p_ext);
    sao.resize(num_p_ext);

    value_type_vec hp_temp;
    hp_temp.resize(num_p_wim);

    padVar(u, u_pad, "xy-periodic");
    padVar(v, v_pad,"xy-periodic");
    padVar(SCP2_array, scp2_pad,"xy-periodic");
    padVar(SCP2I_array, scp2i_pad,"xy-periodic");
    padVar(SCUY_array, scuy_pad,"xy-periodic");
    padVar(SCVX_array, scvx_pad,"xy-periodic");

#if 1
    padVar(h, h_pad,M_advopt);
#else
    //apply "steady" forcing here by adjusting the far-left ghost cells
    //doesn't work for some reason if nwavefreq>1
    padVar(h, h_pad,M_advopt,M_steady);
#endif

    // prediction step
    weno3pdV2(h_pad, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad, sao);

    if (nghost<4)
        throw std::runtime_error("Advection (WENO): 'nghost' should be >=4\n");


#pragma omp parallel for num_threads(max_threads) collapse(1)
   for (int i = 0; i < h_pad.size(); i++)
       hp[i] = h_pad[i]+M_timestep*sao[i];

    // correction step
    weno3pdV2(hp, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad, sao);


    //final output
    //std::cout<<"in WENO\n";
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
        {
            h[ny*i+j] = 0.5*(h_pad[(i+nbdx)*nyext+j+nbdy]+hp[(i+nbdx)*nyext+j+nbdy]+M_timestep*sao[(i+nbdx)*nyext+j+nbdy]);

            //mask land cells
            h[ny*i+j] *= 1-LANDMASK_array[ny*i+j];
        }

    //std::cout<<"min LANDMASK"
    //   <<*std::min_element(LANDMASK_array.data(),LANDMASK_array.data() + LANDMASK_array.num_elements())
    //   <<"\n";
    //std::cout<<"max LANDMASK"
    //   <<*std::max_element(LANDMASK_array.data(),LANDMASK_array.data() + LANDMASK_array.num_elements())
    //   <<"\n";
    //std::cout<<"advected thing at [nx-1,ny-1]: "<<h[nx-1][ny-1]<<"\n";

}//waveAdvWeno

template<typename T>
void WimDiscr<T>::weno3pdV2(value_type_vec const& gin, value_type_vec const& u, value_type_vec const& v, value_type_vec const& scuy,
                       value_type_vec const& scvx, value_type_vec const& scp2i, value_type_vec const& scp2, value_type_vec& saoout)
{

	value_type cq00=-1./2 ,cq01=3./2, cq10=1./2, cq11=1./2, ca0=1./3, ca1=2./3, eps=1e-12;
	value_type q0, q1, a0, a1, q;
	int im1, im2, ip1, jm1, jm2, jp1, ymargin;

    value_type_vec ful, fuh, fvl, fvh, gt;

    int num_p_ext = nxext*nyext;
    ful.resize(num_p_ext);
    fuh.resize(num_p_ext);
    fvl.resize(num_p_ext);
    fvh.resize(num_p_ext);
    gt.resize (num_p_ext);

    if (M_advdim == 2)
        ymargin = 1;
    else
        ymargin = 0;

    // fluxes in x direction
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 2; i < nxext-1; i++)
    {
        for (int j = 0; j < nyext; j++)
        {
            value_type q0, q1, a0, a1, q;
            int im1, im2, ip1, jm1, jm2, jp1, ymargin;

            im1 = i-1;

            if (u[nyext*i+j] > 0.)
            {
                // coefficents to calc higher-order fluxes
                im2 = im1-1;
                q0 = cq00*gin[im2*nyext+j]+cq01*gin[im1*nyext+j];
                q1 = cq10*gin[im1*nyext+j]+cq11*gin[nyext*i+j];
                a0 = ca0;
                a1 = ca1*(std::abs(gin[im2*nyext+j]-gin[im1*nyext+j])+eps)/(std::abs(gin[im1*nyext+j]-gin[nyext*i+j])+eps);

                // lower-order fluxes
                ful[nyext*i+j] = u[nyext*i+j]*gin[im1*nyext+j]*scuy[nyext*i+j];

            }
            else
            {
                // coefficents to calc higher-order fluxes
                ip1 = i+1;
                q0 = cq11*gin[im1*nyext+j]+cq10*gin[i*nyext+j];
                q1 = cq01*gin[nyext*i+j]+cq00*gin[ip1*nyext+j];
                a0 = ca1;
                a1 = ca0*(abs(gin[im1*nyext+j]-gin[nyext*i+j])+eps)/(abs(gin[nyext*i+j]-gin[ip1*nyext+j])+eps);

                // lower-order fluxes
                ful[nyext*i+j] = u[nyext*i+j]*gin[nyext*i+j]*scuy[nyext*i+j];
            }

            // higher-order fluxes
            fuh[nyext*i+j] = (u[nyext*i+j]*(a0*q0+a1*q1)*scuy[nyext*i+j]/(a0+a1))-ful[nyext*i+j];
        }
    }

    // fluxes in y direction
    if (M_advdim == 2)
    {
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nxext; i++)
        {
            for (int j = 2; j < nyext-1; j++)
            {
                value_type q0, q1, a0, a1, q;
                int im1, im2, ip1, jm1, jm2, jp1, ymargin;

                jm1 = j-1;

                if (v[nyext*i+j] > 0.)
                {
                    jm2 = jm1-1;
                    q0 = cq00*gin[nyext*i+jm2]+cq01*gin[nyext*i+jm1];
                    q1 = cq10*gin[nyext*i+jm1]+cq11*gin[nyext*i+j];
                    a0 = ca0;
                    a1 = ca1*(std::abs(gin[nyext*i+jm2]-gin[nyext*i+jm1])+eps)/(std::abs(gin[nyext*i+jm1]-gin[nyext*i+j])+eps);
                    fvl[nyext*i+j] = v[nyext*i+j]*gin[nyext*i+jm1]*scvx[nyext*i+j];
                }
                else
                {
                    jp1 = j+1;
                    q0 = cq11*gin[nyext*i+jm1]+cq10*gin[nyext*i+j];
                    q1 = cq01*gin[nyext*i+j]+cq00*gin[nyext*i+jp1];
                    a0 = ca1;
                    a1 = ca0*(abs(gin[nyext*i+jm1]-gin[nyext*i+j])+eps)/(abs(gin[nyext*i+j]-gin[nyext*i+jp1])+eps);
                    fvl[nyext*i+j] = v[nyext*i+j]*gin[nyext*i+j]*scvx[nyext*i+j];
                }

                fvh[nyext*i+j] = (v[nyext*i+j]*(a0*q0+a1*q1)*scvx[nyext*i+j]/(a0+a1))-fvl[nyext*i+j];
            }
        }
    }

    // update field with low order fluxes
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nxext-1; i++)
    {
        for (int j = 0; j < nyext-ymargin; j++)//nyext-1 if 2d advection; else nyext
        {
            if (M_advdim == 2)
            {
                gt[nyext*i+j] = gin[nyext*i+j]-M_timestep*(ful[(i+1)*nyext+j]-ful[nyext*i+j]+fvl[nyext*i+j+1]-fvl[nyext*i+j])*scp2i[nyext*i+j];
            }
            else if (M_advdim == 1)
            {
                gt[nyext*i+j] = gin[nyext*i+j]-M_timestep*(ful[(i+1)*nyext+j]-ful[nyext*i+j])*scp2i[nyext*i+j];
            }
        }
    }

    q = 0.25/M_timestep;

    // obtain fluxes with limited high order correction fluxes
    // - x dirn
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 1; i < nxext; i++)
    {
        for (int j = 0; j < nyext; j++)
        {
            fuh[nyext*i+j] = ful[nyext*i+j]+
               std::max(-q*gt[nyext*i+j]*scp2[nyext*i+j],
                        std::min(q*gt[(i-1)*nyext+j]*scp2[(i-1)*nyext+j],fuh[nyext*i+j]));
        }
    }

    // obtain fluxes with limited high order correction fluxes
    // - y dirn
    if (M_advdim == 2)
    {
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nxext; i++)
        {
            for (int j = 1; j < nyext; j++)
            {
                fvh[nyext*i+j]=fvl[nyext*i+j]+
                   std::max(-q*gt[nyext*i+j]*scp2[nyext*i+j],
                            std::min(q*gt[nyext*i+j-1]*scp2[nyext*i+j-1],fvh[nyext*i+j]));
            }
        }
    }

#if 1
    // compute the spatial advective operator
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nxext-1; i++)
    {
        for (int j = 0; j < nyext-ymargin; j++)
        {
            if (M_advdim == 2)
            {
                saoout[nyext*i+j] = -(fuh[(i+1)*nyext+j]-fuh[nyext*i+j]+fvh[nyext*i+j+1]-fvh[nyext*i+j])*scp2i[nyext*i+j];
            }
            else if (M_advdim == 1)
            {
                saoout[nyext*i+j] = -(fuh[(i+1)*nyext+j]-fuh[nyext*i+j])*scp2i[nyext*i+j];
            }
        }
    }
#endif

}//weno3pdV2

template<typename T>
void WimDiscr<T>::padVar(value_type_vec const& u, value_type_vec& upad,
        std::string const & advopt, bool const & steady)
{
    int num_p_ext   = nxext*nyext;
    upad.resize(num_p_ext);

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nxext; i++)
    {
        for (int j = 0; j < nyext; j++)
        {

            bool i_inner = ((nbdx-1 < i) && (i < nx+nbdx));
            bool j_inner = ((nbdy-1 < j) && (j < ny+nbdy));//also works for adv_dim==1 (-1<j<ny: ie all j)

            // interior cells
            if ( i_inner && j_inner )
                upad[nyext*i+j] = u[(i-nbdx)*ny+j-nbdy];

            // apply steady conditions here by setting the far-left ghost cells
            // to be the same as the far-left "real" cells
            if( steady && i<nbdx )
            {
                int ju = std::max(0,std::min(ny,j-nbdy));
                upad[nyext*i+j] = u[nyext*nbdx+ju];
            }

            if (M_advdim == 1)
            {
                if (advopt == "xy-periodic")
                {
                    // make periodic in i
                    if ((i < nbdx) && j_inner && (!steady) )
                        //far-left cells
                        upad[nyext*i+j] = u[(nx-nbdx+i)*ny+j-nbdy];

                    if ((nx+nbdx-1 < i) && j_inner )
                        //far-right cells
                        upad[nyext*i+j] = u[(i-nx-nbdx)*ny+j-nbdy];
                }
            }
            else if (M_advdim == 2)
            {
                if (advopt != "notperiodic")
                {
                    // ie either y-periodic or xy-periodic

                    // make periodic in j
                    // - lower cells
                    if ((j < nbdy) && i_inner)
                        upad[nyext*i+j] = u[(i-nbdx)*ny+ny-nbdy+j];

                    // - upper cells
                    if ((ny+nbdy-1 < j) && i_inner)
                        upad[nyext*i+j] = u[(i-nbdx)*ny+j-ny-nbdy];
                }

                if (advopt == "xy-periodic")
                {
                    // make periodic in i
                    // - far-left cells
                    if ((i < nbdx) && j_inner )
                        upad[nyext*i+j] = u[(nx-nbdx+i)*ny+j-nbdy];

                    // - far-right cells
                    if ((nx+nbdx-1 < i) && j_inner )
                        upad[nyext*i+j] = u[(i-nx-nbdx)*ny+j-nbdy];

                    // TL
                    if ((i < nbdx) && (ny+nbdy-1 < j))
                        upad[nyext*i+j] = u[(i+nx-nbdx)*ny+j-ny-nbdy];

                    // BL
                    if ((i < nbdx) && (j < nbdy))
                        upad[nyext*i+j] = u[(i+nx-nbdx)*ny+j];

                    // TR
                    if ((nx+nbdx-1 < i) && (ny+nbdy-1 < j))
                        upad[nyext*i+j] = u[(i-nx-nbdx)*ny+j-ny-nbdy];

                    // BR
                    if ((nx+nbdx-1 < i) && (j < nbdy))
                        upad[nyext*i+j] = u[(i-nx-nbdx)*ny+ny-nbdy+j];

                }//advopt=="xy-periodic"
            }//M_advdim==2
        }//j
    }//i
}//padVar

template<typename T>
void WimDiscr<T>::calcMWD()
{
    value_type adv_dir, wt_theta, om, CF;
    value_type_vec cmom0,cmom_dir,CSfreq, cmom_dir0;
    cmom0.resize(M_num_elements);
    cmom_dir.resize(M_num_elements);
    CSfreq.resize(M_num_elements);
    cmom_dir0.resize(M_num_elements);

    if (nwavedirn == 1)
        wt_theta = 1.;
    else
        wt_theta = 2.0*PI/(nwavedirn);

    // spectral moments
    std::fill( cmom0.begin(), cmom0.end(), 0. );
    std::fill( cmom_dir.begin(), cmom_dir.end(), 0. );

    for (int fq = 0; fq < nwavefreq; fq++)
    {
        om = 2*PI*freq_vec[fq];

        std::fill( CSfreq.begin(), CSfreq.end(), 0. );
        std::fill( cmom_dir0.begin(), cmom_dir0.end(), 0. );

        for (int nth = 0; nth < nwavedirn; nth++)
        {
            adv_dir = -PI*(90.0+wavedir[nth])/180.0;

#pragma omp parallel for num_threads(max_threads) collapse(1)
            for (int i = 0; i < M_num_elements; i++)
            {
                CSfreq[i] += wt_theta*sdf_dir[fq][nth][i];
                cmom_dir0[i] += wt_theta*sdf_dir[fq][nth][i]*adv_dir;
            }
        }//end loop over directions (still in freq loop)

#pragma omp parallel for num_threads(max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
        {
            if (ref_Hs_ice)
                CF = std::pow(disp_ratio[i][fq],2);
            else
                CF = 1.;

            cmom0[i]    += std::abs(wt_om[fq]*CF*CSfreq[i]);
            cmom_dir[i] += std::abs(wt_om[fq]*CF*cmom_dir0[i]);
        }//end spatial loop i
    }//end freq loop

    //assign mwd
    std::fill( mwd.begin(), mwd.end(), 0. );

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {
        if (cmom0[i] > 0.)
            mwd[i] = -90.-180*(cmom_dir[i]/cmom0[i])/PI;
    }//end spatial loop i

}//end calcMWD

template<typename T>
typename WimDiscr<T>::value_type
WimDiscr<T>::thetaDirFrac(value_type const& th1_, value_type const& dtheta_, value_type const& mwd_)
{
    //get mwd\pm90 inside [th1_,th1_+360)
    value_type phi1  = thetaInRange(mwd_-90.,th1_);//>=th1_
    value_type phi2  = thetaInRange(mwd_+90.,th1_);//>=th1_
    value_type th2_  = th1_+dtheta_;

    value_type integral = 0.;
    if (phi2>phi1)
    {
       // th1_,phi1,phi2, and th2_
       value_type L1    = std::max(th1_,phi1);
       value_type L2    = std::min(th2_,phi2);
       L2               = std::max(L1,L2); //make L2>=L1
       value_type chi1  = PI*(L1-mwd_)/180.;
       value_type chi2  = PI*(L2-mwd_)/180.;
       integral        += 2.*(chi2-chi1)+std::sin(2*chi2)-std::sin(2*chi1);
    }
    else
    {
       // th1_,phi2,phi1, and th2_
       // 1st consider (th1_,phi2) interval
       value_type L1    = th1_;
       value_type L2    = min(th2_,phi2);
       value_type chi1  = PI*(L1-mwd_)/180.;
       value_type chi2  = PI*(L2-mwd_)/180.;
       integral        += 2.*(chi2-chi1)+std::sin(2*chi2)-std::sin(2*chi1);
      
       // 2nd consider (phi1,th2_) interval
       L1           = phi1;
       L2           = max(L1,th2_);      //make L2>=L1
       chi1         = PI*(L1-mwd_)/180.;
       chi2         = PI*(L2-mwd_)/180.;
       integral    += 2.*(chi2-chi1)+std::sin(2*chi2)-std::sin(2*chi1);
    }

    value_type theta_dirfrac  = integral/(2.*PI);
    return theta_dirfrac;

}//thetaDirFrac

template<typename T>
typename WimDiscr<T>::value_type
WimDiscr<T>::thetaInRange(value_type const& th_, value_type const& th1, bool const& close_on_right)
{
    value_type th2, dth, th;
    int njump;

    th2   = th1 + 360.;
    if (th_ < th1)
    {
        dth   = th1 - th_;
        njump = std::ceil(dth/360.);
        th    = th_ + njump*360.;
    }
    else if (th_ > th2)
    {
        dth   = th_ - th2;
        njump = std::ceil(dth/360.);
        th = th_ - njump*360.;
    }
    else if (th_ == th2)
    {
        th = th1;
    }
    else
    {
        th = th_;
    }

    if (close_on_right && th == th1)
        th = th2;

    return th;
}//thetaInRange


template<typename T> 
typename WimDiscr<T>::WimGrid WimDiscr<T>::wimGrid(std::string const& units)
{
    value_type fac = 1.;
    if (units == "m")
        fac = 1.;
    else if (units == "km")
        fac = 1.e-3;
    else
    {
        std::cout<<"Units <<"<<units<<">> not implemented\n";
    }

    std::vector<value_type> X(nx*ny);
    std::vector<value_type> Y(nx*ny);
    std::vector<value_type> x(nx);
    std::vector<value_type> y(ny);

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i=0;i<nx;i++)
    {
        for (int j=0;j<ny;j++)
        {
            X[i*ny+j]  = fac*X_array[ny*i+j];//row major (C)
            Y[i*ny+j]  = fac*Y_array[ny*i+j];//row major (C)
            x[i]       = fac*x_col[i];
            y[j]       = fac*y_row[j];
        }
    }

    WimGrid wim_grid =
    {
        nx : nx,
        ny : ny,
        dx : fac*dx,
        dy : fac*dy,
        X  : X,
        Y  : Y,
        x  : x,
        y  : y
    };

    return wim_grid;
}//WimGrid

template<typename T>
void WimDiscr<T>::readDataFromFile(std::string const& filein)
{
    Fdmax.resize(num_p_wim);
    Ftaux.resize(num_p_wim);
    Ftauy.resize(num_p_wim);
    Fhs.resize(num_p_wim);
    Ftp.resize(num_p_wim);

    char * senv = ::getenv( "WIM2D_PATH" );
    std::string str = std::string( senv ) + "/fortran/run/out/binaries/prog";
    fs::path path(str);

    std::string _filein = (boost::format("%1%/%2%") % path.string() % filein).str();

    // s::path path(str);
    // path /= "outputs/binaries/prog";

    std::fstream in(_filein, std::ios::binary | std::ios::in);

    // NB data is in fortran order (column major)
    // - convert to row major
    if (in.is_open())
    {
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Fdmax[ny*i+j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Ftaux[ny*i+j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Ftauy[ny*i+j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Fhs[ny*i+j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Ftp[ny*i+j], sizeof(int));

        in.close();
    }
    else
    {
        std::cout << "Cannot open " << _filein  << "\n";
        std::cerr << "error: open file " << _filein << " for input failed!" <<"\n";
        std::abort();
    }
}//readDataFromFile

template<typename T>
void WimDiscr<T>::exportResults(std::string const& output_type,
                                 value_type const& t_out) const
{

    typedef struct extractFields {
        bool icec;
        bool iceh;
        bool Dmax;
        bool taux;
        bool tauy;
        bool sdx;
        bool sdy;
        bool swh;
        bool mwp;
        bool mwd;
        bool swh_in;
        bool mwp_in;
        bool mwd_in;
    } extractFields;

    extractFields extract_fields;
    int Nrecs,recno;

    std::string str = vm["wim.outparentdir"].template as<std::string>();
    fs::path path(str);
    path   /= "binaries/"+output_type;

    std::string init_time  = ptime(init_time_str);
    std::string timestpstr = ptime(init_time_str, t_out);
    std::string fileout,fileoutb;
    if ( output_type == "prog" )
    {
        Nrecs   = 8;
        fileout = (boost::format( "%1%/wim_prog%2%" ) % path.string() % timestpstr).str();

        //fields to extract
        extract_fields  =
        {
            icec:   false,
            iceh:   false,
            Dmax:   true,
            taux:   true,
            tauy:   true,
            sdx:    true,
            sdy:    true,
            swh:    true,
            mwp:    true,
            mwd:    true,
            swh_in: false,
            mwp_in: false,
            mwd_in: false
        };
    }
    else if ( output_type == "final" )
    {
        Nrecs   = 8;
        fileout = (boost::format( "%1%/wim_out%2%" ) % path.string() % timestpstr).str();

        //fields to extract
        extract_fields  =
        {
            icec:   false,
            iceh:   false,
            Dmax:   true,
            taux:   true,
            tauy:   true,
            sdx:    true,
            sdy:    true,
            swh:    true,
            mwp:    true,
            mwd:    true,
            swh_in: false,
            mwp_in: false,
            mwd_in: false
        };
    }
    else if ( output_type == "init" )
    {
        Nrecs   = 6;
        fileout = (boost::format( "%1%/wim_init%2%" ) % path.string() % timestpstr).str();

        //fields to extract
        extract_fields  =
        {
            icec:   true,
            iceh:   true,
            Dmax:   true,
            taux:   false,
            tauy:   false,
            sdx:    false,
            sdy:    false,
            swh:    true,
            mwp:    true,
            mwd:    true,
            swh_in: false,
            mwp_in: false,
            mwd_in: false
        };
    }
    else if ( output_type == "incwaves" )
    {
        Nrecs   = 3;
        fileout = (boost::format( "%1%/wim_inc%2%" ) % path.string() % timestpstr).str();

        //fields to extract
        extract_fields  =
        {
            icec:   false,
            iceh:   false,
            Dmax:   false,
            taux:   false,
            tauy:   false,
            sdx:    false,
            sdy:    false,
            swh:    false,
            mwp:    false,
            mwd:    false,
            swh_in: true,
            mwp_in: true,
            mwd_in: true
        };
    }
    else if ( output_type == "nextwim" )
    {
        Nrecs   = 8;
        fileout = (boost::format( "%1%/nextwim%2%" ) % path.string() % timestpstr).str();

        //fields to extract
        extract_fields  =
        {
            icec:   false,
            iceh:   false,
            Dmax:   true,
            taux:   true,
            tauy:   true,
            sdx:    true,
            sdy:    true,
            swh:    true,
            mwp:    true,
            mwd:    true,
            swh_in: false,
            mwp_in: false,
            mwd_in: false
        };
    }

    fileoutb   = fileout+".b";
    fileout    = fileout+".a";
    if ( !fs::exists(path) )
        fs::create_directories(path);

    std::fstream out(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!out.is_open())
    {
       std::cout << "Cannot open " << fileout  << "\n";
       std::cerr << "error: open file " << fileout << " for output failed!" <<"\n";
       std::abort();
    }

    // export the txt file for grid field information
    std::fstream outb(fileoutb, std::ios::out | std::ios::trunc);
    if (!outb.is_open())
    {
        std::cout << "Cannot open " << fileoutb  << "\n";
        std::cerr << "error: open file " << fileoutb << " for output failed!" <<"\n";
        std::abort();
    }

    // ==================================================================================================
    //.b file header
    std::string rstr   = std::string(2-std::to_string(Nrecs).length(),'0')
                          + std::to_string(Nrecs);
    outb << std::setw(15) << std::left << rstr  << "    Nrecs    # "<< "Number of records" <<"\n";
    outb << std::setw(15) << std::left << "0"   << "    Norder   # "
         << "Storage order [column-major (F/matlab) = 1; row-major (C) = 0]" <<"\n";
    outb << std::setw(15) << std::left << nx    << "    nx       # "<< "Record length in x direction (elements)" <<"\n";
    outb << std::setw(15) << std::left << ny    << "    ny       # "<< "Record length in y direction (elements)" <<"\n";

    outb <<"\n";
    outb << std::left << init_time << "    t_start    # "<< "Model time of WIM call" <<"\n";
    outb << std::left << timestpstr << "    t_out    # "<< "Model time of output" <<"\n";

    outb <<"\n";
    outb << "Record number and name:" <<"\n";
    // ==================================================================================================

    recno = 1;
    if ( extract_fields.icec )
    {
        //icec,iceh,dfloe,swh,mwp,mwd
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&wim_ice.conc[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "icec" <<"\n";
        ++recno;
    }

    if ( extract_fields.iceh )
    {
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&wim_ice.thick[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "iceh" <<"\n";
        ++recno;
    }

    if ( extract_fields.Dmax )
    {
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&wim_ice.dfloe[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "Dmax" <<"\n";
        ++recno;
    }

    if ( extract_fields.taux )
    {
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&tau_x[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "tau_x" <<"\n";
        ++recno;
    }

    if ( extract_fields.tauy )
    {
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&tau_y[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "tau_y" <<"\n";
        ++recno;
    }

    if ( extract_fields.sdx )
    {
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&stokes_drift_x[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "stokes_drift_x" <<"\n";
        ++recno;
    }

    if ( extract_fields.sdy )
    {
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&stokes_drift_y[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "stokes_drift_y" <<"\n";
        ++recno;
    }

    if ( extract_fields.swh )
    {
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&Hs[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "Hs" <<"\n";
        ++recno;
    }

    if ( extract_fields.mwp )
    {
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&Tp[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "Tp" <<"\n";
        ++recno;
    }

    if ( extract_fields.mwd )
    {
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&mwd[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "mwd" <<"\n";
        ++recno;
    }

    if ( extract_fields.swh_in )
    {
        //swh_in,mwp_in,mwd_in: test incident waves
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&swh_in_array[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "Hs" <<"\n";
        ++recno;
    }

    if ( extract_fields.mwp_in )
    {
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&mwp_in_array[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "Tp" <<"\n";
        ++recno;
    }

    if ( extract_fields.mwd_in )
    {
        for (int i = 0; i < M_num_elements; i++)
            out.write((char *)&mwd_in_array[i], sizeof(value_type));

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "mwd" <<"\n";
        ++recno;
    }

    out.close();
    outb.close();
}//exportResults

template<typename T>
void WimDiscr<T>::testInterp(std::string const& output_type,
                             value_type const& t_out,
                             std::vector<std::vector<value_type>> const& vectors,
                             std::vector<std::string> const& names
                             ) const
{

    int Nrecs     = vectors.size();
    int Nelements = vectors[0].size();
    int recno;

    std::string str = vm["wim.outparentdir"].template as<std::string>();

    std::string init_time  = ptime(init_time_str);
    std::string timestpstr = ptime(init_time_str, t_out);

    if (output_type == "grid")
    {
        fs::path path(str);
        path   /= "binaries/test_interp_grid";
        std::string fileout,fileoutb;

        fileout     = (boost::format( "%1%/wim_test_interp_grid%2%" ) % path.string() % timestpstr).str();
        fileoutb   = fileout+".b";
        fileout    = fileout+".a";
        if ( !fs::exists(path) )
            fs::create_directories(path);

        std::fstream out(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
        if (!out.is_open())
        {
           std::cout << "Cannot open " << fileout  << "\n";
           std::cerr << "error: open file " << fileout << " for output failed!" <<"\n";
           std::abort();
        }

        // export the txt file for grid field information
        std::fstream outb(fileoutb, std::ios::out | std::ios::trunc);
        if (!outb.is_open())
        {
            std::cout << "Cannot open " << fileoutb  << "\n";
            std::cerr << "error: open file " << fileoutb << " for output failed!" <<"\n";
            std::abort();
        }

        // ==================================================================================================
        //.b file header
        std::string rstr   = std::string(2-std::to_string(Nrecs).length(),'0')
                              + std::to_string(Nrecs);
        outb << std::setw(15) << std::left << rstr  << "    Nrecs    # "<< "Number of records" <<"\n";
        outb << std::setw(15) << std::left << "0"   << "    Norder   # "
             << "Storage order [column-major (F/matlab) = 1; row-major (C) = 0]" <<"\n";
        outb << std::setw(15) << std::left << nx    << "    nx       # "<< "Record length in x direction (elements)" <<"\n";
        outb << std::setw(15) << std::left << ny    << "    ny       # "<< "Record length in y direction (elements)" <<"\n";

        outb <<"\n";
        outb << std::left << init_time << "    t_start    # "<< "Model time of WIM call" <<"\n";
        outb << std::left << timestpstr << "    t_out    # "<< "Model time of output" <<"\n";

        outb <<"\n";
        outb << "Record number and name:" <<"\n";
        // ==================================================================================================

        // =================================================================
        //loop over vectors & names
        for (int recno=0;recno<Nrecs;recno++)
        {
            for (int i = 0; i < Nelements; i++)
            {
                out.write((char *)&vectors[recno][i], sizeof(value_type));
            }

            int r   = recno+1;
            rstr   = std::string(2-std::to_string(r).length(),'0')
               + std::to_string(r);
            outb << std::setw(9) << rstr << names[recno] <<"\n";
        }//finish writing to files
        // =================================================================

        out.close();
        outb.close();
    }//export on grid
    else
    {
        //export on mesh
        fs::path path(str);
        path   /= "binaries/test_interp_mesh";
        std::string fileout,fileoutb;

        fileout     = (boost::format( "%1%/field_%2%" ) % path.string() % timestpstr).str();
        fileoutb   = fileout+".dat";
        fileout    = fileout+".bin";
        if ( !fs::exists(path) )
            fs::create_directories(path);

        std::fstream out(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
        if (!out.is_open())
        {
           std::cout << "Cannot open " << fileout  << "\n";
           std::cerr << "error: open file " << fileout << " for output failed!" <<"\n";
           std::abort();
        }

        // export the txt file for grid field information
        std::fstream outb(fileoutb, std::ios::out | std::ios::trunc);
        if (!outb.is_open())
        {
            std::cout << "Cannot open " << fileoutb  << "\n";
            std::cerr << "error: open file " << fileoutb << " for output failed!" <<"\n";
            std::abort();
        }

        // =================================================================
        //loop over vectors & names
        for (int recno=0;recno<Nrecs;recno++)
        {
            std::string prec;
            if (sizeof(value_type)==8)
                prec    = "double";
            else
                prec    = "float";

            out.write((char*) &Nelements, sizeof(Nelements));
            for (int i = 0; i < Nelements; i++)
            {
                out.write((char *)&vectors[recno][i], sizeof(value_type));
            }

            outb << names[recno] << " " << prec <<"\n";
        }//finish writing to files
        // =================================================================

        out.close();
        outb.close();

        //get mesh file & copy to correct location:
        std::string meshdir   = vm["simul.output_directory"].template as<std::string>();
        std::string meshfile  = meshdir+"/mesh_1001.bin";
        std::string meshfile2 = meshdir+"/mesh_1001.dat";
        fileout     = (boost::format( "%1%/mesh_%2%" ) % path.string() % timestpstr).str();
        fileoutb    = fileout+".dat";
        fileout     = fileout+".bin";

        //cp meshfile fileout;
        //cp meshfile2 fileoutb;
        fs::copy_file(fs::path(meshfile),fs::path(fileout),fs::copy_option::overwrite_if_exists);
        fs::copy_file(fs::path(meshfile2),fs::path(fileoutb),fs::copy_option::overwrite_if_exists);
    }//export on grid
}//testInterp

template<typename T>
void WimDiscr<T>::saveLog(value_type const& t_out) const
{
    std::string str = vm["wim.outparentdir"].template as<std::string>();
    fs::path path(str);
    path /= "diagnostics/global";
    if ( !fs::exists(path) )
       fs::create_directories(path);

    std::string init_time  = ptime(init_time_str);
    std::string timestpstr = ptime(init_time_str, t_out);
    std::string fileout    = (boost::format( "%1%/WIMdiagnostics%2%.txt" )
                                % path.string()
                                % timestpstr).str();

    std::fstream out(fileout, std::ios::out | std::ios::trunc);
    if ( !out.is_open() )
    {
        std::cout << "Cannot open " << fileout  << "\n";
        std::cerr << "error: open file " << fileout << " for output failed!" <<"\n";
        std::abort();
    }

    int log_width = 34;

    out << "***********************************************\n";
    out << "Outer subroutine:\n";
    out << ">> " << "wimdiscr.cpp\n\n";
    out << std::left << std::setw(log_width) << "Start time"<<" : " << init_time << "\n";
    out << std::left << std::setw(log_width) << "Call time" <<" : " << timestpstr << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "Main parameters:" << "\n";
    out << std::left << std::setw(log_width) << "SCATMOD"<<" : " << scatmod << "\n";
    out << std::left << std::setw(log_width) << "ADV_DIM"<<" : " << M_advdim << "\n";
    out << std::left << std::setw(log_width) << "ADV_OPT"<<" : " << M_advopt << "\n";
#if 0
    //TODO implement brkopt
    out << std::left << std::setw(log_width) << "BRK_OPT:" << brkopt << "\n";
    if (BRK_OPT.eq.0) then
       write(fid,'(a)'),'(No breaking)'
    elseif (BRK_OPT.eq.1) then
       write(fid,'(a)'),'(Williams et al, 2013, Oc Mod)'
    elseif (BRK_OPT.eq.2) then
       write(fid,'(a)'),'(Marchenko)'
    elseif (BRK_OPT.eq.3) then
       write(fid,'(a)'),'(Mohr-Coulomb)'
    end if
#endif
    out << std::left << std::setw(log_width) << "STEADY"<<" : " << M_steady << "\n";
    out << std::left << std::setw(log_width) << "DO_ATTEN"<<" : " << atten << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "Other integer parameters:" << "\n";
    out << std::left << std::setw(log_width) << "FSD_OPT"<<" : " << fsdopt << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "WIM parameters:" << "\n";
    out << std::left << std::setw(log_width) << "Brine volume fraction"  <<" : " << vbf << "\n";
    out << std::left << std::setw(log_width) << "Youngs modulus (Pa)"    <<" : " << young << "\n";
    out << std::left << std::setw(log_width) << "Flexural strength (Pa)" <<" : " << sigma_c << "\n";
//    out << std::left << std::setw(log_width) << "Breaking stress (Pa)"<<" : " << stress_c << "\n";
    out << std::left << std::setw(log_width) << "Breaking strain"        <<" : " << epsc << "\n";
    out << std::left << std::setw(log_width) << "Damping (Pa.s/m)"       <<" : " << drag_rp << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "FSD parameters:" << "\n";
    out << std::left << std::setw(log_width) << "Dmin (m)"  <<" : " << dmin << "\n";
    out << std::left << std::setw(log_width) << "xi"        <<" : " << xi << "\n";
    out << std::left << std::setw(log_width) << "fragility" <<" : " << fragility << "\n";
    out << std::left << std::setw(log_width) << "Dthresh"   <<" : " << dfloe_miz_thresh << "\n";
    out << std::left << std::setw(log_width) << "cice_min"  <<" : " << cice_min << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "Other parameters:" << "\n";
    out << std::left << std::setw(log_width) << "Time step (s)" <<" : "     << M_timestep << "\n";
    out << std::left << std::setw(log_width) << "CFL number"                <<" : " << M_cfl << "\n";
    out << std::left << std::setw(log_width) << "Max wave group vel (m/s)"  <<" : " << amax << "\n";
    out << std::left << std::setw(log_width) << "Number of time steps"      <<" : " << nt << "\n";
    out << std::left << std::setw(log_width) << "Time interval (h)"         <<" : " << duration/60.0/60.0 << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << std::left << std::setw(log_width) << "Grid dimensions"        <<" : " << nx << ", " << ny << "\n";
    out << std::left << std::setw(log_width) << "Spatial resolution (km)"<<" : " << dx/1.0e3 << ", " << dy/1.0e3 << "\n";
    out << std::left << std::setw(log_width) << "Extent of domain (km)"  <<" : "
        << nx*dx/1.e3 << ", " << ny*dy/1.e3 << "\n";
    out << std::left << std::setw(log_width) << "Minimum period (s)"          <<" : " << 1.0/freq_vec[nwavefreq-1] << "\n";
    out << std::left << std::setw(log_width) << "Maximum period (s)"          <<" : " << 1.0/freq_vec[0] << "\n";
    out << std::left << std::setw(log_width) << "Number of wave frequencies"  <<" : " << nwavefreq << "\n";
    out << std::left << std::setw(log_width) << "Number of wave directions"   <<" : "  << nwavedirn << "\n";
    out << std::left << std::setw(log_width) << "Directional resolution (deg)"<<" : " << 360.0/nwavedirn << "\n";
    out << "***********************************************\n";

    value_type taux_min  = *std::min_element(tau_x.begin(), tau_x.end());
    value_type taux_max  = *std::max_element(tau_x.begin(), tau_x.end());
    value_type tauy_min  = *std::min_element(tau_y.begin(), tau_y.end());
    value_type tauy_max  = *std::max_element(tau_y.begin(), tau_y.end());
    value_type Hs_max    = *std::max_element(Hs.begin(), Hs.end());

    //MIZ diagnostics
    value_type Dmax_min = 10.e3;
    value_type Dmax_max = 0.e3;
    int Nmiz   = 0;

#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int j = 0; j < wim_ice.dfloe.size(); j++)
    {
       if ((wim_ice.dfloe[j]<dfloe_pack_init)&&(wim_ice.dfloe[j]>0))
       {
          ++Nmiz;
          Dmax_min   = std::min(wim_ice.dfloe[j],Dmax_min);
          Dmax_max   = std::max(wim_ice.dfloe[j],Dmax_max);
       }
    }

#define WIMDIAG1D
#if defined (WIMDIAG1D)
    //this definition of MIZ only works in 1d geometries
    value_type W_miz;
    if ( ny == 1 )
    {
       W_miz = (Nmiz*dx);
    }
    else if ( vm["wim.landon3edges"].template as<bool>() )
    {
       W_miz = (Nmiz*dx)/(ny-2);
    }
    else
    {
       W_miz = (Nmiz*dx)/ny;
    }
#endif

    out << "\n***********************************************\n";
    out << "Diagnostics:\n";
#if defined (WIMDIAG1D)
    out << std::left << std::setw(log_width) << "MIZ width (km)"<<" : "         << W_miz/1.e3 << "\n";
#endif
    out << std::left << std::setw(log_width) << "Dmax range in MIZ (m)" <<" : " << Dmax_min << ", " << Dmax_max << "\n";
    out << std::left << std::setw(log_width) << "tau_x range (Pa)"      <<" : " << taux_min << ", " << taux_max << "\n";
    out << std::left << std::setw(log_width) << "tau_y range (Pa)"      <<" : " << tauy_min << ", " << tauy_max << "\n";
    out << std::left << std::setw(log_width) << "Hs max (m)"            <<" : " << Hs_max << "\n";
    out << "***********************************************\n";

    out.close();
}//saveLog

template<typename T>
void WimDiscr<T>::saveOptionsLog()
{
    std::string fileout = vm["wim.outparentdir"].template as<std::string>();
    fileout += "/diagnostics/global";
    fs::path path(fileout);
    if ( !fs::exists(path) )
       fs::create_directories(path);

    fileout += "/wimoptions.log";

    std::fstream logfile(fileout, std::ios::out | std::ios::trunc);
    std::cout << "Writing log file " << fileout << "...\n";

    int log_width = 55;
    if (logfile.is_open())
    {
        for (po::variables_map::iterator it = vm.begin(); it != vm.end(); it++)
        {
            logfile << std::setw(log_width) << std::left << it->first;

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
            }
            else
            { // Assumes that the only remainder is vector<string>
                try
                {
                    std::vector<std::string> vect = vm[it->first].as<std::vector<std::string> >();
                    int i = 0;
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
        }//loop over options in variable map
    }//check if file opens
}//saveOptionsLog

// instantiate wim class for type float
//template class WimDiscr<float>;

// instantiate wim class for type double
template class WimDiscr<double>;

} // namespace WIM2D
