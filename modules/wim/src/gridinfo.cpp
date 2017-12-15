/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   wimdiscr.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug  3 11:52:35 2015
 */

#include <gridinfo.hpp>
#include <date_wim.hpp>
#include <exporter.hpp>
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
GridInfo<T>::GridInfo(T_vmap const& vmIn)
{
    M_max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
    vm = vmIn;
    M_initialised = true;
    M_use_regular = (vm["wim.useregulargridtools"].template as<bool>());

    M_gridfile    = vm["wim.gridfilename"].template as<std::string>();
    if ( M_gridfile != "" )
    {
        std::cout<<"Getting WIM grid from file: "<<M_gridfile<<"\n";
        this->readGridFromFile();
    }
    else
    {
        std::cout<<"Generating WIM grid manually...\n";
        chrono.restart();

        M_num_px = vm["wim.M_num_px"].template as<int>();
        M_num_py = vm["wim.M_num_py"].template as<int>();
        M_dx     = vm["wim.M_dx"].template as<double>();
        M_dy     = vm["wim.M_dy"].template as<double>();
        M_xmin   = vm["wim.M_xmin"].template as<double>();
        M_ymin   = vm["wim.M_ymin"].template as<double>();
        M_xmax   = M_xmin+(M_num_px-1)*M_dx;
        M_ymax   = M_ymin+(M_num_py-1)*M_dy;

        this->gridFromParameters();

        // land mask
        // - add land on 3 edges (upper,lower,RH)
        if (vm["wim.landon3edges"].template as<bool>())
        {
#pragma omp parallel for num_threads(M_max_threads) collapse(1)
            for (int j = 0; j < M_num_py; j++)
            {
                int i = M_num_px-1;
                M_land_mask[M_num_py*i+j] = 1.;
            }
#pragma omp parallel for num_threads(M_max_threads) collapse(1)
            for (int i = 0; i < M_num_px; i++)
            {
                int j = 0;
                M_land_mask[M_num_py*i+j] = 1.;
                j = M_num_py-1;
                M_land_mask[M_num_py*i+j] = 1.;
            }
        }

        std::cout<<"Grid generation done in "<< chrono.elapsed() <<"s\n";
    }//manual grid

    this->gridPostProcessing();//save grid, set some other variables

}//GridInfo


template<typename T>
GridInfo<T>::GridInfo(T_vmap const &vmIn,T_mesh &mesh_in)
{

    M_max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
    vm = vmIn;
    M_initialised = true;
    M_use_regular = (vm["wim.useregulargridtools"].template as<bool>());

    std::cout<<"Generating WIM grid from mesh...\n";
    chrono.restart();

    M_xmin = *std::min_element(mesh_in.M_nodes_x.begin(),mesh_in.M_nodes_x.end());
    M_ymin = *std::min_element(mesh_in.M_nodes_y.begin(),mesh_in.M_nodes_y.end());
    M_xmax = *std::max_element(mesh_in.M_nodes_x.begin(),mesh_in.M_nodes_x.end());
    M_ymax = *std::max_element(mesh_in.M_nodes_y.begin(),mesh_in.M_nodes_y.end());

    //parameters required for regular grid
    std::cout<<"Resolution (km) = "<<mesh_in.M_resolution/1.e3<<"\n";
    M_num_px  = std::ceil((M_xmax-M_xmin)/mesh_in.M_resolution);//round up to lower resolution
    M_num_py  = std::ceil((M_ymax-M_ymin)/mesh_in.M_resolution);
    M_dx  = (M_xmax-M_xmin)/M_num_px;
    M_dy  = (M_ymax-M_ymin)/M_num_py;

    this->gridFromParameters();

    //set land mask
    T_val_vec mask_in(mesh_in.M_num_nodes,0.);
    T_val_vec_ptrs interp_out = {&M_land_mask};
    T_val_vec_ptrs interp_in  = {&mask_in};
    //if(M_use_regular&&M_regular)
    if(false)
        //default value doesn't work well for InterpMeshToGrid
        mesh_in.interpToGrid(interp_out,interp_in,M_xmin,M_ymax,M_num_px,M_num_py,M_dx,M_dy,1.);
    else
    {
        //seems to work
        std::vector<int> tmp = {};
        mesh_in.interpToPoints(interp_out,interp_in,M_px,M_py,tmp,true,1.);
    }

    std::cout<<"Grid generation done in "<< chrono.elapsed() <<"s\n";

    this->gridPostProcessing();//save grid, set some other variables

}//GridInfo


template<typename T>
void GridInfo<T>::gridFromParameters()
{
    // * sets the following arrays:
    // M_px.resize(boost::extents[M_num_px][M_num_py]);
    // M_py.resize(boost::extents[M_num_px][M_num_py]);
    // M_scuy.resize(boost::extents[M_num_px][M_num_py]);
    // M_scvx.resize(boost::extents[M_num_px][M_num_py]);
    // M_scp2.resize(boost::extents[M_num_px][M_num_py]);
    // M_scp2i.resize(boost::extents[M_num_px][M_num_py]);
    // M_land_mask.resize(boost::extents[M_num_px][M_num_py]);

    M_regular   = true;

    std::cout<<"M_num_px,M_num_py = "<<M_num_px<<","<<M_num_py<<"\n";
    std::cout<<"M_dx,M_dy = "<<M_dx<<","<<M_dy<<"\n";
    std::cout<<"M_xmin,M_ymin = "<<M_xmin<<","<<M_ymin<<"\n";

    
    M_num_p    = M_num_px*M_num_py;//number of p points
    M_num_q    = (M_num_px+1)*(M_num_py+1);//number of q points
    M_num_u    = (M_num_px+1)*M_num_py;//number of u points
    M_num_v    = M_num_px*(M_num_py+1);//number of v points
    M_px       .assign(M_num_p,0.);
    M_py       .assign(M_num_p,0.);
    M_scuy     .assign(M_num_u,0.);
    M_scvx     .assign(M_num_v,0.);
    M_scp2     .assign(M_num_p,0.);
    M_scp2i    .assign(M_num_p,0.);
    M_land_mask.assign(M_num_p,0.);


#pragma omp parallel for num_threads(M_max_threads) collapse(2)
    for (int i = 0; i < M_num_px; i++)
        for (int j = 0; j < M_num_py; j++)
        {
            M_px[M_num_py*i+j] = M_xmin + i*M_dx+.5*M_dx;
            M_py[M_num_py*i+j] = M_ymin + j*M_dy+.5*M_dy;
            M_scp2[M_num_py*i+j] = M_dx*M_dy;
            M_scp2i[M_num_py*i+j] = 1./(M_dx*M_dy);
        }

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < M_num_u; i++)
            M_scuy[i] = M_dy;

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < M_num_v; i++)
            M_scvx[i] = M_dx;

}//gridFromParameters


template<typename T>
void GridInfo<T>::gridPostProcessing()
{
    M_use_regular   = (vm["wim.useregulargridtools"].template as<bool>());
    bool DoSaveGrid = (vm["wim.checkprog"].template as<bool>())
                   || (vm["wim.checkinit"].template as<bool>())
                   || (vm["wim.checkfinal"].template as<bool>())
                   || (vm["nextwim.exportresults"].template as<bool>());

    //std::cout<<" ---before saving\n";
    if (DoSaveGrid)
       this->saveGrid(); //save grid to binary
    //std::cout<<" ---after saving\n";

    //for use by interpToPoints
    M_px_vec.resize(M_num_px);
    M_py_vec.resize(M_num_py);
#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < M_num_px; i++)
        M_px_vec[i] = M_px[M_num_py*i];

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int j = 0; j < M_num_py; j++)
        M_py_vec[j] = M_py[j];

    //Do triangulation
    if (!(M_regular&&M_use_regular))
    {
        // wet cells
        // - centres become nodes for a new mesh that can be used for interpolation
        T_val_vec nodes_x,nodes_y;
        for (int i=0;i<M_num_p;i++)
            if (M_land_mask[i]<.5)
            {
                nodes_x.push_back(M_px[i]);
                nodes_y.push_back(M_py[i]);
                M_wet_indices.push_back(i);
            }

        M_triangulation = T_mesh(nodes_x,nodes_y);
    }

    //global variable needed by assign(), all loops
    M_num_elements  = M_px.size();
    std::cout<<"on grid, M_num_elements = "<<M_num_elements<<"\n";

    //length scale to determine the time step from (CFL criterion)
    M_resolution = std::min(M_dx,M_dy);

    this->setupAdvection();
}//gridPostProcessing


template<typename T>
void GridInfo<T>::setupAdvection()
{
    //advection stuff
    M_advdim = vm["wim.advdim"].template as<int>();
    M_advopt = vm["wim.advopt"].template as<std::string>();

    //need nghost>=4 for WENO advection
    T_val nghost  = 4;
    M_nbdy_x    = nghost;
    M_nbdy_y    = nghost;

    if (M_advdim == 1)
        M_nbdy_y = 0;

    M_num_px_ext = M_num_px+2*M_nbdy_x;
    M_num_py_ext = M_num_py+2*M_nbdy_y;//M_num_py if M_advdim==1

    //temp variables, but don't want to assign memory each time advect is called
    int num_p_ext = M_num_px_ext*M_num_py_ext;
    Mtmp_sao      .resize(num_p_ext);
    Mtmp_hp       .resize(num_p_ext);
    Mtmp_hp_tmp   .resize(num_p_ext);
    Mtmp_u_pad    .resize(num_p_ext);
    Mtmp_v_pad    .resize(num_p_ext);
    Mtmp_scp2_pad .resize(num_p_ext);
    Mtmp_scp2i_pad.resize(num_p_ext);
    Mtmp_scuy_pad .resize(num_p_ext);
    Mtmp_scvx_pad .resize(num_p_ext);
    Mtmp_h_pad    .resize(num_p_ext);

    //weno3pdV2
    Mtmp_ful.resize(num_p_ext);
    Mtmp_fuh.resize(num_p_ext);
    Mtmp_fvl.resize(num_p_ext);
    Mtmp_fvh.resize(num_p_ext);
    Mtmp_gt .resize(num_p_ext);
}//setupAdvection()


template<typename T>
void GridInfo<T>::saveGrid()
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
        for (int i = 0; i < M_px.size(); i++)
            out.write((char *)&M_px[i], sizeof(T_val));

        for (int i = 0; i < M_py.size(); i++)
            out.write((char *)&M_py[i], sizeof(T_val));

        for (int i = 0; i < M_scuy.size(); i++)
            out.write((char *)&M_scuy[i], sizeof(T_val));

        for (int i = 0; i < M_scvx.size(); i++)
            out.write((char *)&M_scvx[i], sizeof(T_val));

        for (int i = 0; i < M_scp2.size(); i++)
            out.write((char *)&M_scp2[i], sizeof(T_val));

        for (int i = 0; i < M_scp2i.size(); i++)
            out.write((char *)&M_scp2i[i], sizeof(T_val));

        for (int i = 0; i < M_land_mask.size(); i++)
            out.write((char *)&M_land_mask[i], sizeof(T_val));

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

    std::string nxstr = std::string(4-std::to_string(M_num_px).length(),'0') + std::to_string(M_num_px);
    std::string nystr = std::string(4-std::to_string(M_num_py).length(),'0') + std::to_string(M_num_py);

    // std::cout<<"-----------M_num_px= "<< nxstr <<"\n";
    // std::cout<<"-----------M_num_py= "<< nystr <<"\n";

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
void GridInfo<T>::readGridFromFile()
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
    str = M_gridfile;
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
    M_num_px = record[2];
    M_num_py = record[3];
    std::cout<<"M_num_px= "<< M_num_px <<"\n";
    std::cout<<"M_num_py= "<< M_num_py <<"\n";
    int nbytes = record[4];
    std::cout<<"nbytes= "<< nbytes <<"\n";

    M_num_p    = M_num_px*M_num_py;//number of p points
    M_num_q    = (M_num_px+1)*(M_num_py+1);//number of q points
    M_num_u    = (M_num_px+1)*M_num_py;//number of u points
    M_num_v    = M_num_px*(M_num_py+1);//number of v points
    T_val_vec PLat_array, PLon_array;

    std::fstream in( afilename, std::ios::binary | std::ios::in);

    int off = M_num_q*nbytes*2/* skip qlon and qlat*/;

    this->readFromBinary(in, PLon_array, off, std::ios::beg);
    this->readFromBinary(in, PLat_array);

    off = M_num_u*nbytes*2/* skip ulon and ulat*/;
    off += M_num_v*nbytes*2/* skip vlon and vlat*/;

    this->readFromBinary(in, M_scuy, off, std::ios::cur, 1, 0);
    this->readFromBinary(in, M_scvx, 0, std::ios::cur, 0, 1);
    this->readFromBinary(in, M_scp2);
    this->readFromBinary(in, M_land_mask);

    M_px.resize(M_num_p);
    M_py.resize(M_num_p);
    M_scp2i.resize(M_num_p);


    // polar stereographic projection
    mapx_class *map;
    std::string filename = (boost::format("%1%/%2%") % path.string() % "NpsNextsim.mpp").str();
    std::cout<<"stereographic description file= "<< filename <<"\n";

    std::vector<char> _str(filename.begin(), filename.end());
    _str.push_back('\0');

    map = init_mapx(&_str[0]);

    T_val dx_min = 1.e30;
    T_val dy_min = 1.e30;
    T_val dx_max = -1.e30;
    T_val dy_max = -1.e30;
    for (int k = 0; k < M_num_p; k++)
    {
        double x, y;
        int status = forward_mapx(map,PLat_array[k],PLon_array[k],&x,&y);

        M_px[k]      = x;
        M_py[k]      = y;
        M_scp2i[k]  = 1./M_scp2[k];

        //k=i*M_num_py+j
        int i = k/M_num_py;
        int j = k%M_num_py;
        if (i>1)
        {
            T_val dx_ = M_px[M_num_py*(i)+j]-M_px[M_num_py*(i-1)+j];
            dx_min  = std::min(dx_min,dx_);
            dx_max  = std::max(dx_max,dx_);
        }
        if (j>1)
        {
            T_val dy_ = M_py[M_num_py*(i)+j]-M_py[M_num_py*(i)+j-1];
            dy_min  = std::min(dy_min,dy_);
            dy_max  = std::max(dy_max,dy_);
        }
    }
    close_mapx(map);

    M_dx = dx_min;
    M_dy = dy_min;

    // is the grid regular?
    T_val tol = 1.;//threshold for range in M_dx,M_dy (m)
    if ((dx_max-dx_min>tol)||(dy_max-dy_min>tol))
        M_regular = false;

    std::cout<<"M_dx= "<< M_dx <<"\n";
    std::cout<<"M_dy= "<< M_dy <<"\n";

    std::cout<<"Reading grid done...\n";
}//readGridFromFile


template<typename T>
void GridInfo<T>::readFromBinary(std::fstream &in, T_val_vec& in_array, int off, std::ios_base::seekdir direction, int addx, int addy)
{
    if (off && (in.is_open()))
    {
        in.seekg(off, direction); // skip from the direction (beginning/current/end) position of the file
    }

    int nx_in = M_num_px+addx;
    int ny_in = M_num_py+addy;
    in_array.resize(nx_in*ny_in);

    if (in.is_open())
    {
        // NB assumes files are saved with fortran ordering
        for (int j = 0; j < ny_in; j++)
            for (int i = 0; i < nx_in; i++)
                in.read((char *)&in_array[ny_in*i+j], sizeof(T_val));
    }
    else
    {
        std::cout << "readFromBinary: Cannot open file\n";
        std::cerr << "error: opening file for input failed!" <<"\n";
        std::abort();
    }
}//readFromBinary


template<typename T>
void GridInfo<T>::interpToPoints(
        T_val_vec_ptrs &output_data,       //output data
        T_val_vec_ptrs const &input_data,  //input data
        T_val_vec &Rx, T_val_vec &Ry) //output locations
{

    if (!(M_regular&&M_use_regular))
    {
        M_triangulation.interpToPoints(output_data,input_data,Rx,Ry,M_wet_indices);
        return;
    }
        
    int nb_var      = input_data.size();
    int Ninterp     = (*(input_data[0])).size();//get pointer, then get size
    int target_size = Rx.size();

    T_val_vec interp_in(Ninterp*nb_var);   //input to interp routine
    for (int i=0;i<Ninterp;i++)
        for (int p=0;p<nb_var;p++)
            interp_in[nb_var*i+p]   = (*(input_data[p]))[i];

    T_val* interp_out;

    std::cout<<"gridToPoints: InterpFromGridToMeshx\n";
    int interptype = BilinearInterpEnum;
    InterpFromGridToMeshx(interp_out,       //data (out)
                          &M_px_vec[0], M_num_px,    //x vector (source), length of x vector
                          &M_py_vec[0], M_num_py,    //y vector (source), length of y vector
                          &interp_in[0],    //data (in)
                          M_num_py, M_num_px,           //M,N: no of grid cells in y,x directions
                                            //(to determine if corners or centers of grid have been input)
                          nb_var,           //no of variables
                          &Rx[0],           // x vector (target) (NB already a reference)
                          &Ry[0],           // x vector (target) (NB already a reference)
                          target_size,0.,   //target_size,default value
                          interptype,       //interpolation type
                          true              //row_major (false = fortran/matlab order)         
                          );

    //output
    for (int p=0;p<nb_var;p++)
    {
        output_data[p]->assign(target_size,0.);
        for (int i=0;i<target_size;i++)
            (*(output_data[p]))[i]  = interp_out[nb_var*i+p];
    }

    xDelete<T_val>(interp_out);

}//interpToPoints


template<typename T>
void GridInfo<T>::interpFromMesh(T_mesh &mesh,
        T_val_vec_ptrs &output_data,       //output data
        T_val_vec_ptrs const &input_data)  //input data
{
    if(this->isRegular())
        mesh.interpToGrid(output_data,input_data,
                M_xmin,M_ymax,M_num_px,M_num_py,M_dx,M_dy);
    else
        mesh.interpToPoints(output_data,input_data,M_px,M_py);
}//meshToGrid


template<typename T>
void GridInfo<T>::waveAdvWeno(T_val_vec& h, T_val_vec const& u, T_val_vec const& v,T_val const& timestep)
{
    std::fill( Mtmp_hp    .begin(),Mtmp_hp    .end(),0.);
    std::fill( Mtmp_sao   .begin(),Mtmp_sao   .end(),0.);
    std::fill( Mtmp_hp_tmp.begin(),Mtmp_hp_tmp.end(),0.);

    //pad variables
    padVar(u, Mtmp_u_pad, "xy-periodic");
    padVar(v, Mtmp_v_pad,"xy-periodic");
    padVar(M_scp2, Mtmp_scp2_pad,"xy-periodic");
    padVar(M_scp2i, Mtmp_scp2i_pad,"xy-periodic");
    padVar(M_scuy, Mtmp_scuy_pad,"xy-periodic");
    padVar(M_scvx, Mtmp_scvx_pad,"xy-periodic");

#if 1
    padVar(h, Mtmp_h_pad,M_advopt);
#else
    //apply "steady" forcing here by adjusting the far-left ghost cells
    //doesn't work for some reason if nwavefreq>1
    padVar(h, Mtmp_h_pad,M_advopt,M_steady);
#endif

    // prediction step
    weno3pdV2(Mtmp_h_pad, Mtmp_u_pad, Mtmp_v_pad, Mtmp_scuy_pad, Mtmp_scvx_pad, Mtmp_scp2i_pad, Mtmp_scp2_pad, Mtmp_sao,timestep);

    if(M_advdim==1)
        if (M_nbdy_x<4)
            throw std::runtime_error("Advection (WENO): 'nghost' should be >=4\n");
    else
        if ((M_nbdy_x<4)||(M_nbdy_y<4))
            throw std::runtime_error("Advection (WENO): 'nghost' should be >=4\n");


#pragma omp parallel for num_threads(M_max_threads) collapse(1)
   for (int i = 0; i < Mtmp_h_pad.size(); i++)
       Mtmp_hp[i] = Mtmp_h_pad[i]+timestep*Mtmp_sao[i];

    // correction step
    weno3pdV2(Mtmp_hp, Mtmp_u_pad, Mtmp_v_pad, Mtmp_scuy_pad, Mtmp_scvx_pad, Mtmp_scp2i_pad, Mtmp_scp2_pad, Mtmp_sao,timestep);


    //final output
    //std::cout<<"in WENO\n";
#pragma omp parallel for num_threads(M_max_threads) collapse(2)
    for (int i = 0; i < M_num_px; i++)
        for (int j = 0; j < M_num_py; j++)
        {
            h[M_num_py*i+j] = 0.5*(Mtmp_h_pad[(i+M_nbdy_x)*M_num_py_ext+j+M_nbdy_y]+Mtmp_hp[(i+M_nbdy_x)*M_num_py_ext+j+M_nbdy_y]+timestep*Mtmp_sao[(i+M_nbdy_x)*M_num_py_ext+j+M_nbdy_y]);

            //mask land cells
            h[M_num_py*i+j] *= 1-M_land_mask[M_num_py*i+j];
        }

    //std::cout<<"min LANDMASK"
    //   <<*std::min_element(M_land_mask.data(),M_land_mask.data() + M_land_mask.num_elements())
    //   <<"\n";
    //std::cout<<"max LANDMASK"
    //   <<*std::max_element(M_land_mask.data(),M_land_mask.data() + M_land_mask.num_elements())
    //   <<"\n";
    //std::cout<<"advected thing at [M_num_px-1,M_num_py-1]: "<<h[M_num_px-1][M_num_py-1]<<"\n";

}//waveAdvWeno


template<typename T>
void GridInfo<T>::weno3pdV2(T_val_vec const& gin, T_val_vec const& u, T_val_vec const& v, T_val_vec const& scuy,
                       T_val_vec const& scvx, T_val_vec const& scp2i, T_val_vec const& scp2, T_val_vec& saoout, T_val const& timestep)
{

	T_val cq00=-1./2 ,cq01=3./2, cq10=1./2, cq11=1./2, ca0=1./3, ca1=2./3, eps=1e-12;
	T_val q0, q1, a0, a1, q;
	int im1, im2, ip1, jm1, jm2, jp1, ymargin;

    std::fill( Mtmp_ful.begin(),Mtmp_ful.end(),0.);
    std::fill( Mtmp_fuh.begin(),Mtmp_fuh.end(),0.);
    std::fill( Mtmp_fvl.begin(),Mtmp_fvl.end(),0.);
    std::fill( Mtmp_fvh.begin(),Mtmp_fvh.end(),0.);
    std::fill( Mtmp_gt .begin(),Mtmp_gt .end(),0.);

    if (M_advdim == 2)
        ymargin = 1;
    else
        ymargin = 0;

    // fluxes in x direction
#pragma omp parallel for num_threads(M_max_threads) collapse(2)
    for (int i = 2; i < M_num_px_ext-1; i++)
    {
        for (int j = 0; j < M_num_py_ext; j++)
        {
            T_val q0, q1, a0, a1, q;
            int im1, im2, ip1, jm1, jm2, jp1, ymargin;

            im1 = i-1;

            if (u[M_num_py_ext*i+j] > 0.)
            {
                // coefficents to calc higher-order fluxes
                im2 = im1-1;
                q0 = cq00*gin[im2*M_num_py_ext+j]+cq01*gin[im1*M_num_py_ext+j];
                q1 = cq10*gin[im1*M_num_py_ext+j]+cq11*gin[M_num_py_ext*i+j];
                a0 = ca0;
                a1 = ca1*(std::abs(gin[im2*M_num_py_ext+j]-gin[im1*M_num_py_ext+j])+eps)/(std::abs(gin[im1*M_num_py_ext+j]-gin[M_num_py_ext*i+j])+eps);

                // lower-order fluxes
                Mtmp_ful[M_num_py_ext*i+j] = u[M_num_py_ext*i+j]*gin[im1*M_num_py_ext+j]*scuy[M_num_py_ext*i+j];

            }
            else
            {
                // coefficents to calc higher-order fluxes
                ip1 = i+1;
                q0 = cq11*gin[im1*M_num_py_ext+j]+cq10*gin[i*M_num_py_ext+j];
                q1 = cq01*gin[M_num_py_ext*i+j]+cq00*gin[ip1*M_num_py_ext+j];
                a0 = ca1;
                a1 = ca0*(abs(gin[im1*M_num_py_ext+j]-gin[M_num_py_ext*i+j])+eps)/(abs(gin[M_num_py_ext*i+j]-gin[ip1*M_num_py_ext+j])+eps);

                // lower-order fluxes
                Mtmp_ful[M_num_py_ext*i+j] = u[M_num_py_ext*i+j]*gin[M_num_py_ext*i+j]*scuy[M_num_py_ext*i+j];
            }

            // higher-order fluxes
            Mtmp_fuh[M_num_py_ext*i+j] = (u[M_num_py_ext*i+j]*(a0*q0+a1*q1)*scuy[M_num_py_ext*i+j]/(a0+a1))-Mtmp_ful[M_num_py_ext*i+j];
        }
    }

    // fluxes in y direction
    if (M_advdim == 2)
    {
#pragma omp parallel for num_threads(M_max_threads) collapse(2)
        for (int i = 0; i < M_num_px_ext; i++)
        {
            for (int j = 2; j < M_num_py_ext-1; j++)
            {
                T_val q0, q1, a0, a1, q;
                int im1, im2, ip1, jm1, jm2, jp1, ymargin;

                jm1 = j-1;

                if (v[M_num_py_ext*i+j] > 0.)
                {
                    jm2 = jm1-1;
                    q0 = cq00*gin[M_num_py_ext*i+jm2]+cq01*gin[M_num_py_ext*i+jm1];
                    q1 = cq10*gin[M_num_py_ext*i+jm1]+cq11*gin[M_num_py_ext*i+j];
                    a0 = ca0;
                    a1 = ca1*(std::abs(gin[M_num_py_ext*i+jm2]-gin[M_num_py_ext*i+jm1])+eps)/(std::abs(gin[M_num_py_ext*i+jm1]-gin[M_num_py_ext*i+j])+eps);
                    Mtmp_fvl[M_num_py_ext*i+j] = v[M_num_py_ext*i+j]*gin[M_num_py_ext*i+jm1]*scvx[M_num_py_ext*i+j];
                }
                else
                {
                    jp1 = j+1;
                    q0 = cq11*gin[M_num_py_ext*i+jm1]+cq10*gin[M_num_py_ext*i+j];
                    q1 = cq01*gin[M_num_py_ext*i+j]+cq00*gin[M_num_py_ext*i+jp1];
                    a0 = ca1;
                    a1 = ca0*(abs(gin[M_num_py_ext*i+jm1]-gin[M_num_py_ext*i+j])+eps)/(abs(gin[M_num_py_ext*i+j]-gin[M_num_py_ext*i+jp1])+eps);
                    Mtmp_fvl[M_num_py_ext*i+j] = v[M_num_py_ext*i+j]*gin[M_num_py_ext*i+j]*scvx[M_num_py_ext*i+j];
                }

                Mtmp_fvh[M_num_py_ext*i+j] = (v[M_num_py_ext*i+j]*(a0*q0+a1*q1)*scvx[M_num_py_ext*i+j]/(a0+a1))-Mtmp_fvl[M_num_py_ext*i+j];
            }
        }
    }

    // update field with low order fluxes
#pragma omp parallel for num_threads(M_max_threads) collapse(2)
    for (int i = 0; i < M_num_px_ext-1; i++)
    {
        for (int j = 0; j < M_num_py_ext-ymargin; j++)//M_num_py_ext-1 if 2d advection; else M_num_py_ext
        {
            if (M_advdim == 2)
            {
                Mtmp_gt[M_num_py_ext*i+j] = gin[M_num_py_ext*i+j]-timestep*(Mtmp_ful[(i+1)*M_num_py_ext+j]-Mtmp_ful[M_num_py_ext*i+j]+Mtmp_fvl[M_num_py_ext*i+j+1]-Mtmp_fvl[M_num_py_ext*i+j])*scp2i[M_num_py_ext*i+j];
            }
            else if (M_advdim == 1)
            {
                Mtmp_gt[M_num_py_ext*i+j] = gin[M_num_py_ext*i+j]-timestep*(Mtmp_ful[(i+1)*M_num_py_ext+j]-Mtmp_ful[M_num_py_ext*i+j])*scp2i[M_num_py_ext*i+j];
            }
        }
    }

    q = 0.25/timestep;

    // obtain fluxes with limited high order correction fluxes
    // - x dirn
#pragma omp parallel for num_threads(M_max_threads) collapse(2)
    for (int i = 1; i < M_num_px_ext; i++)
    {
        for (int j = 0; j < M_num_py_ext; j++)
        {
            Mtmp_fuh[M_num_py_ext*i+j] = Mtmp_ful[M_num_py_ext*i+j]+
               std::max(-q*Mtmp_gt[M_num_py_ext*i+j]*scp2[M_num_py_ext*i+j],
                        std::min(q*Mtmp_gt[(i-1)*M_num_py_ext+j]*scp2[(i-1)*M_num_py_ext+j],Mtmp_fuh[M_num_py_ext*i+j]));
        }
    }

    // obtain fluxes with limited high order correction fluxes
    // - y dirn
    if (M_advdim == 2)
    {
#pragma omp parallel for num_threads(M_max_threads) collapse(2)
        for (int i = 0; i < M_num_px_ext; i++)
        {
            for (int j = 1; j < M_num_py_ext; j++)
            {
                Mtmp_fvh[M_num_py_ext*i+j]=Mtmp_fvl[M_num_py_ext*i+j]+
                   std::max(-q*Mtmp_gt[M_num_py_ext*i+j]*scp2[M_num_py_ext*i+j],
                            std::min(q*Mtmp_gt[M_num_py_ext*i+j-1]*scp2[M_num_py_ext*i+j-1],Mtmp_fvh[M_num_py_ext*i+j]));
            }
        }
    }

#if 1
    // compute the spatial advective operator
#pragma omp parallel for num_threads(M_max_threads) collapse(2)
    for (int i = 0; i < M_num_px_ext-1; i++)
    {
        for (int j = 0; j < M_num_py_ext-ymargin; j++)
        {
            if (M_advdim == 2)
            {
                saoout[M_num_py_ext*i+j] = -(Mtmp_fuh[(i+1)*M_num_py_ext+j]-Mtmp_fuh[M_num_py_ext*i+j]+Mtmp_fvh[M_num_py_ext*i+j+1]-Mtmp_fvh[M_num_py_ext*i+j])*scp2i[M_num_py_ext*i+j];
            }
            else if (M_advdim == 1)
            {
                saoout[M_num_py_ext*i+j] = -(Mtmp_fuh[(i+1)*M_num_py_ext+j]-Mtmp_fuh[M_num_py_ext*i+j])*scp2i[M_num_py_ext*i+j];
            }
        }
    }
#endif

}//weno3pdV2


template<typename T>
void GridInfo<T>::padVar(T_val_vec const& u, T_val_vec& upad,
        std::string const & advopt, bool const & steady)
{

#pragma omp parallel for num_threads(M_max_threads) collapse(2)
    for (int i = 0; i < M_num_px_ext; i++)
    {
        for (int j = 0; j < M_num_py_ext; j++)
        {

            bool i_inner = ((M_nbdy_x-1 < i) && (i < M_num_px+M_nbdy_x));
            bool j_inner = ((M_nbdy_y-1 < j) && (j < M_num_py+M_nbdy_y));//also works for adv_dim==1 (-1<j<M_num_py: ie all j)

            // interior cells
            if ( i_inner && j_inner )
                upad[M_num_py_ext*i+j] = u[(i-M_nbdy_x)*M_num_py+j-M_nbdy_y];

            // apply steady conditions here by setting the far-left ghost cells
            // to be the same as the far-left "real" cells
            if( steady && i<M_nbdy_x )
            {
                int ju = std::max(0,std::min(M_num_py,j-M_nbdy_y));
                upad[M_num_py_ext*i+j] = u[M_num_py_ext*M_nbdy_x+ju];
            }

            if (M_advdim == 1)
            {
                if (advopt == "xy-periodic")
                {
                    // make periodic in i
                    if ((i < M_nbdy_x) && j_inner && (!steady) )
                        //far-left cells
                        upad[M_num_py_ext*i+j] = u[(M_num_px-M_nbdy_x+i)*M_num_py+j-M_nbdy_y];

                    if ((M_num_px+M_nbdy_x-1 < i) && j_inner )
                        //far-right cells
                        upad[M_num_py_ext*i+j] = u[(i-M_num_px-M_nbdy_x)*M_num_py+j-M_nbdy_y];
                }
            }
            else if (M_advdim == 2)
            {
                if (advopt != "notperiodic")
                {
                    // ie either y-periodic or xy-periodic

                    // make periodic in j
                    // - lower cells
                    if ((j < M_nbdy_y) && i_inner)
                        upad[M_num_py_ext*i+j] = u[(i-M_nbdy_x)*M_num_py+M_num_py-M_nbdy_y+j];

                    // - upper cells
                    if ((M_num_py+M_nbdy_y-1 < j) && i_inner)
                        upad[M_num_py_ext*i+j] = u[(i-M_nbdy_x)*M_num_py+j-M_num_py-M_nbdy_y];
                }

                if (advopt == "xy-periodic")
                {
                    // make periodic in i
                    // - far-left cells
                    if ((i < M_nbdy_x) && j_inner )
                        upad[M_num_py_ext*i+j] = u[(M_num_px-M_nbdy_x+i)*M_num_py+j-M_nbdy_y];

                    // - far-right cells
                    if ((M_num_px+M_nbdy_x-1 < i) && j_inner )
                        upad[M_num_py_ext*i+j] = u[(i-M_num_px-M_nbdy_x)*M_num_py+j-M_nbdy_y];

                    // TL
                    if ((i < M_nbdy_x) && (M_num_py+M_nbdy_y-1 < j))
                        upad[M_num_py_ext*i+j] = u[(i+M_num_px-M_nbdy_x)*M_num_py+j-M_num_py-M_nbdy_y];

                    // BL
                    if ((i < M_nbdy_x) && (j < M_nbdy_y))
                        upad[M_num_py_ext*i+j] = u[(i+M_num_px-M_nbdy_x)*M_num_py+j];

                    // TR
                    if ((M_num_px+M_nbdy_x-1 < i) && (M_num_py+M_nbdy_y-1 < j))
                        upad[M_num_py_ext*i+j] = u[(i-M_num_px-M_nbdy_x)*M_num_py+j-M_num_py-M_nbdy_y];

                    // BR
                    if ((M_num_px+M_nbdy_x-1 < i) && (j < M_nbdy_y))
                        upad[M_num_py_ext*i+j] = u[(i-M_num_px-M_nbdy_x)*M_num_py+M_num_py-M_nbdy_y+j];

                }//advopt=="xy-periodic"
            }//M_advdim==2
        }//j
    }//i
}//padVar


// instantiate wim class for type float
//template class GridInfo<float>;

// instantiate wim class for type double
template class GridInfo<double>;

} // namespace Wim
