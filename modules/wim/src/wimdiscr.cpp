/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   wimdiscr.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug  3 11:52:35 2015
 */

#include <wimdiscr.hpp>

namespace Wim
{

template<typename T>
void WimDiscr<T>::gridProcessing()
{
    X_array.resize(boost::extents[nx][ny]);
    Y_array.resize(boost::extents[nx][ny]);
    SCUY_array.resize(boost::extents[nx][ny]);
    SCVX_array.resize(boost::extents[nx][ny]);
    SCP2_array.resize(boost::extents[nx][ny]);
    SCP2I_array.resize(boost::extents[nx][ny]);
    LANDMASK_array.resize(boost::extents[nx][ny]);

    dx = vm["wim.dx"].template as<double>();
    dy = vm["wim.dy"].template as<double>();

    x0 = vm["wim.xmin"].template as<double>();
    y0 = vm["wim.ymin"].template as<double>();

    // int thread_id;
    // int total_threads;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
    // std::cout<<"MAX THREADS= "<< max_threads <<"\n";

    std::cout<<"grid generation starts\n";
    chrono.restart();

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            X_array[i][j] = x0 + i*dx+.5*dx;
            Y_array[i][j] = y0 + j*dy+.5*dy;
            SCUY_array[i][j] = dy;
            SCVX_array[i][j] = dx;
            SCP2_array[i][j] = dx*dy;
            SCP2I_array[i][j] = 1./(dx*dy);

            //add land on 3 edges (upper,lower,RH)
            if (vm["wim.landon3edges"].template as<bool>())
            {
               if (i==nx-1)
               {
                   LANDMASK_array[i][j] = 1.;
               }

               if ((j==0) || (j==ny-1))
               {
                   LANDMASK_array[i][j] = 1.;
               }
            }

        }
    }

    bool critter = (vm["wim.checkprog"].template as<bool>())
               || (vm["nextwim.exportresults"].template as<bool>());
    if (critter)
    {
       //save grid to binary
        std::string str = vm["wim.outparentdir"].template as<std::string>();

        //char * senv = ::getenv( "WIM2D_PATH" );
        //if ( (str == ".") && (senv != NULL) && (senv[0] != '\0') )
        //{
        //    str = std::string( senv ) + "/CXX";
        //}

        fs::path path(str);
        path /= "binaries";

        if ( !fs::exists(path) )
            fs::create_directories(path);


        std::string fileout = (boost::format( "%1%/wim_grid.a" ) % path.string()).str();
        std::fstream out(fileout, std::ios::binary | std::ios::out | std::ios::trunc);

        if (out.is_open())
        {
            for (int i = 0; i < X_array.shape()[0]; i++)
                for (int j = 0; j < X_array.shape()[1]; j++)
                    out.write((char *)&X_array[i][j], sizeof(value_type));

            for (int i = 0; i < Y_array.shape()[0]; i++)
                for (int j = 0; j < Y_array.shape()[1]; j++)
                    out.write((char *)&Y_array[i][j], sizeof(value_type));

            for (int i = 0; i < SCUY_array.shape()[0]; i++)
                for (int j = 0; j < SCUY_array.shape()[1]; j++)
                    out.write((char *)&SCUY_array[i][j], sizeof(value_type));

            for (int i = 0; i < SCVX_array.shape()[0]; i++)
                for (int j = 0; j < SCVX_array.shape()[1]; j++)
                    out.write((char *)&SCVX_array[i][j], sizeof(value_type));

            for (int i = 0; i < SCP2_array.shape()[0]; i++)
                for (int j = 0; j < SCP2_array.shape()[1]; j++)
                    out.write((char *)&SCP2_array[i][j], sizeof(value_type));

            for (int i = 0; i < SCP2I_array.shape()[0]; i++)
                for (int j = 0; j < SCP2I_array.shape()[1]; j++)
                    out.write((char *)&SCP2I_array[i][j], sizeof(value_type));

            for (int i = 0; i < LANDMASK_array.shape()[0]; i++)
                for (int j = 0; j < LANDMASK_array.shape()[1]; j++)
                    out.write((char *)&LANDMASK_array[i][j], sizeof(value_type));

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
    }

    std::cout<<"grid generation done in "<< chrono.elapsed() <<"s\n";
}

template<typename T>
void WimDiscr<T>::readGridFromFile(std::string const& filein)
{
    //array2_type X_array(boost::extents[nx][ny],boost::fortran_storage_order());
    X_array.resize(boost::extents[nx][ny]);
    Y_array.resize(boost::extents[nx][ny]);
    SCUY_array.resize(boost::extents[nx][ny]);
    SCVX_array.resize(boost::extents[nx][ny]);
    SCP2_array.resize(boost::extents[nx][ny]);
    SCP2I_array.resize(boost::extents[nx][ny]);
    LANDMASK_array.resize(boost::extents[nx][ny]);

    dx = vm["wim.dx"].template as<double>();
    dy = vm["wim.dy"].template as<double>();


    std::fstream in(filein, std::ios::binary | std::ios::in);

    if (in.is_open())
    {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&X_array[i][j], sizeof(value_type));

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&Y_array[i][j], sizeof(value_type));


        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&SCUY_array[i][j], sizeof(value_type));


        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&SCVX_array[i][j], sizeof(value_type));


        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&SCP2_array[i][j], sizeof(value_type));


        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&SCP2I_array[i][j], sizeof(value_type));


        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&LANDMASK_array[i][j], sizeof(value_type));

        in.close();
    }
    else
    {
        std::cout << "Cannot open " << filein  << "\n";
        std::cerr << "error: open file " << filein << " for input failed!" <<"\n";
        std::abort();
    }
}

template<typename T>
void WimDiscr<T>::init()
{
    // wim grid generation
    this->gridProcessing();
    //std::cout<<"\nHI\n";
    //this->readGridFromFile("wim_grid.a");

    // parameters
    nwavedirn = vm["wim.nwavedirn"].template as<int>();
    nwavefreq = vm["wim.nwavefreq"].template as<int>();
    advdim = vm["wim.advdim"].template as<int>();
    ref_Hs_ice = vm["wim.refhsice"].template as<bool>();
    atten = vm["wim.atten"].template as<bool>();
    useicevel = vm["wim.useicevel"].template as<bool>();
    steady = vm["wim.steady"].template as<bool>();
    breaking = vm["wim.breaking"].template as<bool>();
    scatmod = vm["wim.scatmod"].template as<std::string>();
    advopt = vm["wim.advopt"].template as<std::string>();
    fsdopt = vm["wim.fsdopt"].template as<std::string>();
    cfl = vm["wim.cfl"].template as<double>();
    Hs_inc = vm["wim.hsinc"].template as<double>(); /* 2.0 */
    Tp_inc = vm["wim.tpinc"].template as<double>(); /* 12.0 */
    mwd_inc = vm["wim.mwdinc"].template as<double>(); /* -90. */ //-135.;//-90.;
    Tmin = vm["wim.tmin"].template as<double>(); /* 2.5 */
    Tmax = vm["wim.tmax"].template as<double>(); /* 25. */
    unifc = vm["wim.unifc"].template as<double>(); /* 0.7 */
    unifh = vm["wim.unifh"].template as<double>(); /* 2.0 */
    dfloe_pack_init = vm["wim.dfloepackinit"].template as<double>(); /* 300.0 */
    dfloe_pack_thresh = vm["wim.dfloepackthresh"].template as<double>(); /* 400.0 */
    young = vm["wim.young"].template as<double>();
    visc_rp = vm["wim.viscrp"].template as<double>();
    wim_itest = vm["wim.itest"].template as<int>();
    wim_jtest = vm["wim.jtest"].template as<int>();
    std::cout<<"\nitest,jtest: "<<wim_itest<<","<<wim_jtest<<"\n";

    //nghost==3 needs padVar (boundary conditions) between prediction/advection step
    //nghost>3 doesn't
    //nghost<3 is not possible
    nghost  = 4;
    nbdx    = nghost;
    nbdy    = nghost;

    if (useicevel)
    {
      std::cerr << std::endl
      << "useicevel=true not implemented"
      << std::endl;
      std::abort();
    }

    if (advdim == 1)
        nbdy = 0;

    nxext = nx+2*nbdx;
    nyext = ny+2*nbdy;//ny if advdim==1

    gravity = 9.81;
    rhowtr = 1025.;
    rhoice = 922.5;
    poisson = 0.3;
    dmin = 20.;
    xi = 2.;
    fragility = 0.9;

    vbf = 0.1;//brine volume fraction
    vb = vbf;
    sigma_c  = (1.76e+6)*std::exp(-5.88*std::sqrt(vbf));//flexural strength (Pa)
    epsc = sigma_c/young;//breaking strain
    flex_rig_coeff = young/(12.0*(1-std::pow(poisson,2.)));
}

template<typename T>
void WimDiscr<T>::assign(std::vector<value_type> const& ice_c, std::vector<value_type> const& ice_h, std::vector<value_type> const& n_floes, bool step) // reset
{
    wt_simp.resize(nwavefreq);
    wt_om.resize(nwavefreq);

    freq_vec.resize(nwavefreq);
    vec_period.resize(nwavefreq);
    wlng.resize(nwavefreq);
    ag.resize(nwavefreq);
    ap.resize(nwavefreq);

    steady_mask.resize(boost::extents[nx][ny]);
    wave_mask.resize(boost::extents[nx][ny]);

    if (!step)
    {
        sdf_dir.resize(boost::extents[nx][ny][nwavedirn][nwavefreq]);
    }

    sdf_inc.resize(boost::extents[nx][ny][nwavedirn][nwavefreq]);

    ag_eff.resize(boost::extents[nx][ny][nwavefreq]);
    ap_eff.resize(boost::extents[nx][ny][nwavefreq]);
    wlng_ice.resize(boost::extents[nx][ny][nwavefreq]);

    ag2d_eff_temp.resize(boost::extents[nx][ny]);
    sdf3d_dir_temp.resize(boost::extents[nx][ny][nwavedirn]);

    atten_nond.resize(boost::extents[nx][ny][nwavefreq]);
    damping.resize(boost::extents[nx][ny][nwavefreq]);
    disp_ratio.resize(boost::extents[nx][ny][nwavefreq]);

    ice_mask.resize(boost::extents[nx][ny]);
    wtr_mask.resize(boost::extents[nx][ny]);

    icec.resize(boost::extents[nx][ny]);
    iceh.resize(boost::extents[nx][ny]);
    dfloe.resize(nx*ny);
    nfloes.resize(nx*ny);

    dave.resize(boost::extents[nx][ny]);
    atten_dim.resize(boost::extents[nx][ny]);
    damp_dim.resize(boost::extents[nx][ny]);

    tau_x.resize(nx*ny);
    tau_y.resize(nx*ny);

    S_freq.resize(boost::extents[nx][ny]);
    taux_om.resize(boost::extents[nx][ny]);
    tauy_om.resize(boost::extents[nx][ny]);
    hp.resize(boost::extents[nxext][nyext]);

    Hs.resize(boost::extents[nx][ny]);
    Tp.resize(boost::extents[nx][ny]);
    mwd.resize(boost::extents[nx][ny]);

    //std::fill( steady_mask.data(), steady_mask.data() + steady_mask.num_elements(), 10 );
    // steady_mask
    if (steady)
        std::fill( &steady_mask[0][0], &steady_mask[3][0], 1. );

    std::fill( disp_ratio.data(), disp_ratio.data() + disp_ratio.num_elements(), 1. );

    // freq_vec
    if (nwavefreq == 1)
    {
        freq_vec[0] = 1./Tp_inc;
    }
    else
    {
        // multiple frequencies
        fmin = 1./Tmax;
        fmax = 1./Tmin;
        df = (fmax-fmin)/(nwavefreq-1);

        for (int i = 0; i < nwavefreq; i++)
        {
            freq_vec[i] = fmin+i*df;
            //freq_vec[i] = fmax-i*df;
        }
    }

    // set directional grid
    if (nwavedirn == 1)
    {
        wavedir.push_back(mwd_inc);
    }
    else
    {
        value_type theta_max = 90.;
        value_type theta_min = -270.;
        value_type dtheta = (theta_min-theta_max)/nwavedirn;

        for (int i = 0; i < nwavedirn; i++)
        {
            wavedir.push_back(theta_max+i*dtheta);
        }
    }


    // wt_simp (weights for Simpson's rule)
    if (nwavefreq >1)
    {
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
        //wt_om = dom*wt_simp/3.0;
        wt_om = wt_simp;
        std::for_each(wt_om.begin(), wt_om.end(), [&](value_type& f){ f = dom*f/3.0; });
    }
    else
    {
        wt_om[0] = 1.;
    }


    // vec_period
    vec_period = freq_vec;
    std::for_each(vec_period.begin(), vec_period.end(), [&](value_type& f){ f = 1./f; });

    // wlng
    wlng = freq_vec;
    std::for_each(wlng.begin(), wlng.end(), [&](value_type& f){ f = gravity/(2*PI*std::pow(f,2.)); });

    // ap
    ap = wlng;
    std::for_each(ap.begin(), ap.end(), [&](value_type& f){ f = std::sqrt(gravity*f/(2*PI)); });

    // ag
    ag = ap;
    std::for_each(ag.begin(), ag.end(), [&](value_type& f){ f = f/2. ; });


    //====================================================
    //ideal ice/waves TODO sep function
    x0 = X_array[0][0];
    xmax = X_array[nx-1][ny-1]; //x0+(nx-1)*dx;

    y0 = Y_array[0][0];
    ym = Y_array[0][0]+(ny-1)*dy;
    //value_type ymax = Y_array[nx-1][ny-1];

    x_edge = 0.5*(x0+xmax)-0.8*(0.5*(xmax-x0));

    //std::cout<<Hs_inc<<" "<<Tp_inc<<" "<<mwd_inc<<std::endl;
    //std::exit(1);
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (X_array[i][j] < x_edge)
            {
                wave_mask[i][j] = 1.;
                Hs[i][j] = Hs_inc;
                Tp[i][j] = Tp_inc;
                mwd[i][j] = mwd_inc;
                //std::cout<<Hs[i][j]<<" "<<Tp[i][j]<<" "<<mwd[i][j];
            }
        }
    }
    //end: ideal ice/waves TODO sep function
    //====================================================


    //====================================================
    //set incident wave spec TODO sep function
    om = 2*PI*freq_vec[0];

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            std::vector<value_type> Sfreq(nwavefreq);
            std::vector<value_type> theta_fac(nwavedirn,0.);
            value_type f1, f2, f3, t_m, om_m, chi, om;

            if (wave_mask[i][j] == 1.)
            {
                if (nwavefreq > 1)
                {
                    for (int fq = 0; fq < nwavefreq; fq++)
                    {
                        om = 2*PI*freq_vec[fq];
                        t_m = 2*PI/om;
                        om_m  = 2*PI/Tp[i][j];
                        f1 = (5.0/16.0)*std::pow(Hs[i][j],2.)*std::pow(om_m,4.);
                        f2 = 1.0/std::pow(om,5.);
                        f3 = std::exp(-1.25*std::pow(t_m/Tp[i][j],4.));
                        Sfreq[fq] = f1*f2*f3;
                    }
                }
                else
                {
                    Sfreq[0] = std::pow(Hs[i][j]/4.0,2.);
                }
            }

            // direcional spreading
            if (nwavedirn == 1)
            {
                std::fill(theta_fac.begin(), theta_fac.end(), 1.);
            }
            else
            {
                value_type dtheta = std::abs(wavedir[1]-wavedir[0]);

                //if (mwd[i][j]!=0.)
                //   std::cout<<"dir-frac ("<<i<<","<<j<<")"<<std::endl;
                for (int dn = 0; dn < nwavedirn; dn++)
                {
#if 0
                    //less accurate way of calculating spreading
                    //(sample cos^2 at mid-point of interval)
                    chi = PI*(wavedir[dn]-mwd[i][j])/180.0;
                    if (std::cos(chi) > 0.)
                        theta_fac[dn] = 2.0*std::pow(std::cos(chi),2.)/PI;
                    else
                        theta_fac[dn] = 0.;
#else
                    //more accurate way of calculating spreading
                    //(integrate cos^2 over interval)
                    theta_fac[dn] = thetaDirFrac(wavedir[dn]-dtheta/2.,dtheta,mwd[i][j]);
                    theta_fac[dn] = theta_fac[dn]*180./(PI*dtheta);
#endif

                    //if (Hs[i][j]!=0.)
                    //   std::cout<<wavedir[dn]<<" "<<mwd[i][j]<<" "
                    //            <<theta_fac[dn]<<std::endl;
                }
            }

            // combine freq and dir
            for (int fq = 0; fq < nwavefreq; fq++)
            {
                for (int dn = 0; dn < nwavedirn; dn++)
                {
                    sdf_inc[i][j][dn][fq] = Sfreq[fq]*theta_fac[dn];

                    // if ((i==1) && (j==1))
                    //     adv_dir = (-PI/180.0)*(wavedir[dn]+90.0);
                }
            }
        }
    }
    //end: set incident wave spec TODO sep function
    //====================================================


    //====================================================
    //ideal ice conditions TODO sep function
    x_edge = 0.5*(x0+xmax)-0.7*(0.5*(xmax-x0));

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (X_array[i][j] < x_edge)
                wtr_mask[i][j] = 1.;

            ice_mask[i][j] = (1-wtr_mask[i][j])*(1-LANDMASK_array[i][j]);

            if (ice_c.size() == 0)
            {
                icec[i][j] = unifc*ice_mask[i][j];
                iceh[i][j] = unifh*ice_mask[i][j];
                dfloe[ny*i+j] = dfloe_pack_init*ice_mask[i][j];
            }
            else
            {
                icec[i][j] = ice_c[ny*i+j];
                nfloes[ny*i+j] = n_floes[ny*i+j];

                dfloe[ny*i+j] = 0.;

                if (icec[i][j] < vm["wim.cicemin"].template as<double>())
                {
                    ice_mask[i][j] = 0.;
                    icec[i][j] = 0.;
                    iceh[i][j] = 0.;
                    nfloes[ny*i+j] = 0.;
                }
                else
                {
                    ice_mask[i][j] = 1.;

                    iceh[i][j] = ice_h[ny*i+j]/icec[i][j];

                    dfloe[ny*i+j] = std::sqrt(icec[i][j]/nfloes[ny*i+j]);
                }

                if (dfloe[ny*i+j] > dfloe_pack_thresh)
                    dfloe[ny*i+j] = dfloe_pack_init;

                //std::cout<<"----------------------complex case\n";
            }

            //std::cout<<"ICE_MASK["<< i  << "," << j << "]= " << ice_mask[i][j] <<"\n";
            //ice_mask[i][j] = (1-wtr_mask[i][j])*(1-LANDMASK_array[i][j]);

            if ((LANDMASK_array[i][j] == 1.) /*&& (!step)*/)
                ice_mask[i][j] = 0.;
        }
    }
    //end: ideal ice conditions TODO sep function
    //====================================================

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

    std::cout<<"nfloes_max= "<< *std::max_element(n_floes.begin(), n_floes.end()) <<"\n";
    std::cout<<"nfloes_min= "<< *std::min_element(n_floes.begin(), n_floes.end()) <<"\n";
    std::cout<<"nfloes_acc= "<< std::accumulate(n_floes.begin(), n_floes.end(),0.) <<"\n";
#endif

    //std::cout<<"attenuation loop starts (big loop)\n";
    for (int fq = 0; fq < nwavefreq; fq++)
    {
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {

                double params[5];
                params[0] = young;
                params[1] = gravity;
                params[2] = rhowtr;
                params[3] = rhoice;
                params[4] = poisson;

                double outputs[8];

                if (ice_mask[i][j] == 1.)
                {
                    om = 2*PI*freq_vec[fq];

                    if (fq == 0)
                        guess = std::pow(om,2.)/gravity;
                    else
                        guess = 2*PI/wlng_ice[i][j][fq-1];

                    RTparam_outer(outputs,iceh[i][j],double(om),double(visc_rp),double(guess),params);

                    value_type kice, kwtr, int_adm, modT, argR, argT;

                    damping[i][j][fq] = outputs[0];
                    kice = outputs[1];
                    kwtr = outputs[2];
                    int_adm = outputs[3];
                    atten_nond[i][j][fq] = outputs[4];
                    modT = outputs[5];
                    argR = outputs[6];
                    argT = outputs[7];

                    disp_ratio[i][j][fq] = (kice*modT)/kwtr;
                    //wavelength to use in ice
                    if (1)
                    {
                       //use ice wavelength TODO make an option?
                       wlng_ice[i][j][fq] = 2*PI/kice;
                    }
                    else
                    {
                       //use water wavelength instead of ice wavelength
                       wlng_ice[i][j][fq] = wlng[fq];
                    }

                    //group and phase velocities to use in ice
                    if (!useicevel)
                    {
                       //water group and phase velocities
                       //(ice ones not implemented)
                       ag_eff[i][j][fq] = ag[fq];
                       ap_eff[i][j][fq] = ap[fq];
                    }

                }
                else
                {
                    ag_eff[i][j][fq] = ag[fq];
                    ap_eff[i][j][fq] = ap[fq];
                    wlng_ice[i][j][fq] = wlng[fq];
                }
            }
        }

        if (fq == 0)
        {
            amax = *std::max_element(ag_eff.data(),ag_eff.data() + ag_eff.num_elements());
        }

        if (fq == nwavefreq-1)
        {
            amin = *std::min_element(ag_eff.data(),ag_eff.data() + ag_eff.num_elements());
        }
    }
    //std::cout<<"big loop done\n";

    if (!atten)
    {
        std::fill( atten_nond.data(), atten_nond.data() + atten_nond.num_elements(), 0. );
        std::fill( damping.data(), damping.data() + damping.num_elements(), 0. );
    }

    //int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/


    // ========================================================================
    // initialise sdf_dir from sdf_inc
    // - only do this if step==0
    if (!step)
    {
        std::fill( sdf_dir.data(), sdf_dir.data() + sdf_dir.num_elements(), 0. );

        //#pragma omp parallel for num_threads(max_threads) collapse(4)
        for (int i = 0; i < nwavefreq; i++)
        {
            for (int j = 0; j < nwavedirn; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    for (int l = 0; l < ny; l++)
                    {
                        if (wave_mask[k][l] == 1.)
                        {
                            sdf_dir[k][l][j][i] = sdf_inc[k][l][j][i];
                            //std::cout<<"sdf_dir[k][l][j][i]= "<< sdf_dir[k][l][j][i] <<"\n";
                        }
                    }
                }
            }
        }
    }
    // end: initialise sdf_dir from sdf_inc
    // ========================================================================


    //dt = cfl*std::min(dx,dy)/amax;
    amax = *std::max_element(ag_eff.data(),ag_eff.data() + ag_eff.num_elements());
    std::cout<<"dx= "<< dx <<"\n";
    std::cout<<"amax= "<< amax <<"\n";
    std::cout<<"cfl= "<< cfl <<"\n";
    dt = cfl*dx/amax;

    //reduce time step slightly to make duration an integer multiple of dt
    duration = vm["wim.duration"].template as<double>();
    nt = std::ceil(duration/dt);
    dt = duration/nt;


    //ncs = std::ceil(nwavedirn/2);
    ncs = std::round(nwavedirn/2);
}//end: init()


template<typename T>
void WimDiscr<T>::timeStep(bool step)
{
    std::fill( tau_x.begin(), tau_x.end(), 0. );
    std::fill( tau_y.begin(), tau_y.end(), 0. );

    array2_type tmp1, mom0, mom2, var_strain, mom0w, mom2w;
    tmp1.resize(boost::extents[nx][ny]);
    mom0.resize(boost::extents[nx][ny]);
    mom2.resize(boost::extents[nx][ny]);
    var_strain.resize(boost::extents[nx][ny]);
    mom0w.resize(boost::extents[nx][ny]);
    mom2w.resize(boost::extents[nx][ny]);

    value_type E_tot, sig_strain, Pstrain, P_crit, wlng_crest, Dc;
    value_type adv_dir, F, kicel, om, ommin, ommax, om1, lam1, lam2, dom, c1d, tmp;
    int jcrest;
    bool break_criterion,test_ij;

    value_type t_step = cpt*dt;//seconds from start time
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

    if ( dumpDiag )
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
         diagID << std::setw(16) << std::left
            << wim_itest << " # itest\n";
         diagID << std::setw(16) << std::left
            << wim_jtest << " # jtest\n";
         diagID << std::setw(16) << std::left
            << ice_mask[wim_itest][wim_jtest] << " # ICE_MASK\n";
       }
       else
       {
         std::cout << "Cannot open " << diagfile  << "\n";
         std::cerr << "error: open file " << diagfile << " for output failed!" <<"\n";
         std::abort();
       }
       diagID.close();
    }//end dumpDiag

    dom = 2*PI*(freq_vec[nwavefreq-1]-freq_vec[0])/(nwavefreq-1);

    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/


    if (vm["wim.steady"].template as<bool>())
    {
        for (int i = 0; i < nwavefreq; i++)
        {
            for (int j = 0; j < nwavedirn; j++)
            {
                adv_dir = (-PI/180.)*(wavedir[j]+90.);

                if (std::cos(adv_dir) >= 0.)
                {
#pragma omp parallel for num_threads(max_threads) collapse(2)
                    for (int k = 0; k < nx; k++)
                    {
                        for (int l = 0; l < ny; l++)
                        {
                            if (steady_mask[k][l] > 0.)
                            {
                                sdf_dir[k][l][j][i] = sdf_inc[k][l][j][i];
                                //std::cout<<"sdf_dir[k][l][j][i]= "<< sdf_dir[k][l][j][i] <<"\n";
                            }
                        }
                    }
                }

            }
        }
    }

    //calc mean floe size outside of frequency loop;
    std::fill( dave.data(), dave.data() + dave.num_elements(), 0. );
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
       for (int j = 0; j < ny; j++)
       {
          if (ice_mask[i][j] == 1.)
          {
             if (dfloe[ny*i+j] <200.)
             {
                if ( fsdopt == "RG" )
                {
                    floeScaling(dfloe[ny*i+j],1,dave[i][j]);
                }
                else if ( fsdopt == "PowerLawSmooth" )
                {
                    floeScalingSmooth(dfloe[ny*i+j],1,dave[i][j]);
                }
             }
             else
             {
                //just use uniform dist
                dave[i][j] = dfloe[ny*i+j];
             }
          }


          test_ij = (i==wim_itest) && (j==wim_jtest);
          if ( dumpDiag && test_ij )
          {
             std::fstream diagID(diagfile, std::ios::out | std::ios::app);
             diagID << "\nIce info: pre-breaking\n";
             diagID << std::setw(16) << std::left
                << icec[i][j] << " # conc\n";
             diagID << std::setw(16) << std::left
                << iceh[i][j] << " # h, m\n";
             diagID << std::setw(16) << std::left
                << dave[i][j] << " # D_av, m\n";
             diagID << std::setw(16) << std::left
                << dfloe[ny*i+j] << " # D_max, m\n";

             if (atten)
               diagID << "\n# period, s | atten_dim, m^{-1}| damp_dim, m^{-1}\n";

             diagID.close();
          }
       }
    }//end spatial loop i - have dave

    for (int fq = 0; fq < nwavefreq; fq++)
    {
        std::fill( atten_dim.data(), atten_dim.data() + atten_dim.num_elements(), 0. );
        std::fill( damp_dim.data(), damp_dim.data() + damp_dim.num_elements(), 0. );

#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                value_type c1d;
                if ((ice_mask[i][j] == 1.) && (atten))
                {
                    // floes per unit length
                    c1d = icec[i][j]/dave[i][j];

                    // scattering
                    atten_dim[i][j] = atten_nond[i][j][fq]*c1d;

                    // damping
                    damp_dim[i][j] = 2*damping[i][j][fq]*icec[i][j];

                    test_ij = (i==wim_itest) && (j==wim_jtest);
                    if ( dumpDiag && test_ij )
                    {
                       std::fstream diagID(diagfile, std::ios::out | std::ios::app);
                       diagID << vec_period[fq] << "   "
                              << atten_dim[i][j] << "   "
                              << damp_dim[i][j] << "\n";
                       diagID.close();
                    }
                }//end of ice check
            }
        }//end of spatial loop i


        // copy for application of advAttenSimple
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                ag2d_eff_temp[i][j] = ag_eff[i][j][fq];

                for (int dn = 0; dn < nwavedirn; dn++)
                {
                    sdf3d_dir_temp[i][j][dn] = sdf_dir[i][j][dn][fq];
                }
            }
        }

        // std::cout<<"applied advection starts\n";
        if (scatmod == "dissipated")
        {
            advAttenSimple(sdf3d_dir_temp, S_freq, taux_om, tauy_om, ag2d_eff_temp);
        }
        else if (scatmod == "isotropic")
        {
            advAttenIsotropic(sdf3d_dir_temp, S_freq, taux_om, tauy_om, ag2d_eff_temp);
        }
        // std::cout<<"applied advection done\n";

        // update after application of advAttenSimple
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                ag_eff[i][j][fq] = ag2d_eff_temp[i][j];

                for (int dn = 0; dn < nwavedirn; dn++)
                {
                    sdf_dir[i][j][dn][fq] = sdf3d_dir_temp[i][j][dn];
                    // std::cout<<"AFTER: SDIR["<< i << "," << j << "]= "<< sdf_dir[i][j][dn][fq] <<"\n";
                }

                //std::cout<<"AFTER: taux_om["<< i << "," << j << "]= "<< taux_om[i][j] <<"\n";
            }
        }

        // integrate stress densities over frequency
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                tmp1[i][j] = rhowtr*gravity*taux_om[i][j]/ap_eff[i][j][fq];
                tau_x[ny*i+j] += wt_om[fq]*tmp1[i][j];

                tmp1[i][j] = rhowtr*gravity*tauy_om[i][j]/ap_eff[i][j][fq];
                tau_y[ny*i+j] += wt_om[fq]*tmp1[i][j];

                //std::cout<<"tau_x["<< i << "," << j << "]= "<< taux_om[i][j] <<"\n";
            }
        }

        // integrals for breaking program_options
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                value_type adv_dir, F, kicel, om, tmp;
                // convert from water amp's to ice amp's
                F = disp_ratio[i][j][fq];
                kicel = 2*PI/wlng_ice[i][j][fq];

                // spectral moments: take abs as small errors can make S_freq negative
                om = 2*PI*freq_vec[fq];
                tmp = wt_om[fq]*S_freq[i][j];

                // variance of displacement (water)
                mom0w[i][j] += std::abs(tmp);

                // variance of displacement (ice)
                mom0[i][j] += std::abs(tmp*std::pow(F,2.));
                tmp = wt_om[fq]*std::pow(om,2.)*S_freq[i][j];

                // variance of speed (water)
                mom2w[i][j] += std::abs(tmp);

                // variance of speed (ice)
                mom2[i][j] += std::abs(tmp*std::pow(F,2.));

                // variance of strain
                if (ice_mask[i][j] == 1.)
                {
                    // strain conversion factor
                    tmp = F*std::pow(kicel,2.)*iceh[i][j]/2.0;

                    // strain density
                    tmp = wt_om[fq]*S_freq[i][j]*std::pow(tmp,2.);
                    var_strain[i][j] += std::abs(tmp);
                }
            }
        }
    }

    // value_type _min = *std::min_element(mom0w.data(),mom0w.data() + mom0w.num_elements());
    // value_type _max = *std::max_element(mom0w.data(),mom0w.data() + mom0w.num_elements());
    // std::cout<<"Min f= " << _min <<"\n";
    // std::cout<<"Max f= " << _max <<"\n";


    // for (int i = 0; i < nx; i++)
    //     for (int j = 0; j < ny; j++)
    //         std::cout << "VRT[" << i << "," << j << "]= " << var_strain[i][j] <<"\n";

    std::fill( Tp.data(), Tp.data() + Tp.num_elements(), 0. );

    if (ref_Hs_ice)
    {
        Hs = mom0;
        std::for_each(Hs.data(), Hs.data()+Hs.num_elements(), [&](value_type& f){ f = 4*std::sqrt(f); });
        // std::fill( Tp.data(), Tp.data() + Tp.num_elements(), 0. );

#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if (mom2[i][j] > 0.)
                    Tp[i][j] = 2*PI*std::sqrt(mom0[i][j]/mom2[i][j]);
            }
        }
    }
    else
    {
        Hs = mom0w;
        std::for_each(Hs.data(), Hs.data()+Hs.num_elements(), [&](value_type& f){ f = 4*std::sqrt(f); });

        // std::fill( Tp.data(), Tp.data() + Tp.num_elements(), 0. );
        //std::fill( Hs.data(), Hs.data() + Hs.num_elements(), 0. );

#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                //Hs[i][j] = 4*std::sqrt(mom0w[i][j]);
                if (mom2w[i][j] > 0.)
                    Tp[i][j] = 2*PI*std::sqrt(mom0w[i][j]/mom2w[i][j]);
            }
        }
    }

    //update mwd
    calcMWD();

    if ( dumpDiag )
    {
       int i =  wim_itest;
       int j =  wim_jtest;
       std::fstream diagID(diagfile, std::ios::out | std::ios::app);

       diagID << std::setw(16) << std::left
          << mom0w[i][j] << " # mom0w, m^2\n";
       diagID << std::setw(16) << std::left
          << mom2w[i][j] <<" # mom2w, m^2/s^2\n";
       diagID << std::setw(16) << std::left
          << mom0[i][j] <<" # mom0, m^2\n";
       diagID << std::setw(16) << std::left
          << mom2[i][j] <<" # mom2, m^2/s^2\n";
       diagID << std::setw(16) << std::left
          << Hs[i][j] <<" # Hs, m\n";
       diagID << std::setw(16) << std::left
          << Tp[i][j] <<" # Tp, s\n";
       diagID << std::setw(16) << std::left
          << mwd[i][j] <<" # mwd, deg\n";
       diagID << std::setw(16) << std::left
          << tau_x[ny*i+j] <<" # tau_x, Pa\n";
       diagID << std::setw(16) << std::left
          << tau_y[ny*i+j] <<" # tau_y, Pa\n";
       diagID.close();
    }


    if (!(steady) && !(breaking))
    {
       //check energy conservation
       auto temparray = Hs;
       std::for_each(
             temparray.data(), temparray.data()+temparray.num_elements(),
             [&](value_type& f){ f *= f; });
       E_tot = std::accumulate(
             temparray.data(), temparray.data()+temparray.num_elements(),0.);

        // std::fill( var_strain.data(), var_strain.data()+var_strain.num_elements(), 1. );
        // E_tot = std::accumulate(var_strain.data(), var_strain.data()+var_strain.num_elements(),0.);
        // std::cout<<"Sum= "<< E_tot <<"\n";
    }

    // finally do floe breaking

    //std::cout<<"max_threads= "<< max_threads <<"\n";

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            value_type E_tot, sig_strain, Pstrain, P_crit, wlng_crest, Dc;
            value_type adv_dir, F, kicel, om, ommin, ommax, om1, lam1, lam2, c1d, tmp;
            int jcrest;
            bool break_criterion;

            //std::cout << "MASK[" << i << "," << j << "]= " << ice_mask[i][j] << " and "<< mom0[i][j]  <<"\n";

            Pstrain      = 0.;
            P_crit       = std::exp(-1.0);
            if ((ice_mask[i][j] == 1.) && (mom0[i][j] >= 0.))
            {
                // significant strain amp
                sig_strain = 2*std::sqrt(var_strain[i][j]);

                // probability of critical strain
                // being exceeded from Rayleigh distribution
                Pstrain = std::exp( -std::pow(epsc,2.)/(2*var_strain[i][j]) );

                break_criterion = (Pstrain >= P_crit) && breaking;

                if (break_criterion)
                {
                    //std::cout << "TP[" << i << "," << j << "]= " << Tp[i][j] <<"\n";
                    om    = 2*PI/Tp[i][j];
                    ommin = 2*PI*freq_vec[0];
                    ommax = 2*PI*freq_vec[nwavefreq-1];

                    if (om <= ommin)
                        wlng_crest = wlng_ice[i][j][0];
                    else if (om >= ommax)
                        wlng_crest = wlng_ice[i][j][nwavefreq-1];
                    else
                    {
                        jcrest = std::floor((om-ommin+dom)/dom);
                        // std::cout<<"jrest= "<< jcrest <<"\n";
                        om1 = 2*PI*freq_vec[jcrest-1];
                        lam1 = wlng_ice[i][j][jcrest-1];
                        lam2 = wlng_ice[i][j][jcrest];
                        wlng_crest = lam1+(om-om1)*(lam2-lam1)/dom;
                    }

                    Dc = std::max<value_type>(dmin,wlng_crest/2.0);
                    dfloe[ny*i+j] = std::min<value_type>(Dc,dfloe[ny*i+j]);
                    //std::cout<<"DMAX= std::MAX("<< Dc << "," << dfloe[ny*i+j] <<")\n";
                }
            }//end test for ice and waves
            else if (wtr_mask[i][j] == 1.)
            {
                dfloe[ny*i+j] = 0;
            }

            //nfloes[ny*i+j] = 0.;

            if (dfloe[ny*i+j] > 0.)
                nfloes[ny*i+j] = icec[i][j]/std::pow(dfloe[ny*i+j],2.);

            test_ij  = (i==wim_itest)&&(j==wim_jtest);
            if (dumpDiag && (ice_mask[i][j] == 1.) && test_ij)
            {
               //dump diagnostic even if no waves (but only if ice)
               std::fstream diagID(diagfile, std::ios::out | std::ios::app);
               diagID << "\nIce info: post-breaking\n";
               diagID << std::setw(16) << std::left
                  << Pstrain << " # P_strain\n";
               diagID << std::setw(16) << std::left
                  << P_crit << " # P_crit\n";
               diagID << std::setw(16) << std::left
                  << wlng_crest << " # peak wavelength, m\n";
               diagID << std::setw(16) << std::left
                  << dfloe[ny*i+j] << " # D_max, m\n";
               diagID.close();
            }
        }
    }//end spatial loop i


    // for (int i = 0; i < nx; i++)
    //     for (int j = 0; j < ny; j++)
    //         std::cout << "Dmax[" << i << "," << j << "]= " << dfloe[i][j] <<"\n";

    // for (int i = 0; i < nx; i++)
    //     for (int j = 0; j < ny; j++)
    //         std::cout<<"Hs["<< i << "," << j << "]= "<< Hs[i][j] <<"\n";

    // std::cout<<"Hs_max= "<< *std::max_element(Hs.data(), Hs.data()+Hs.num_elements()) <<"\n";
    // std::cout<<"Hs_min= "<< *std::min_element(Hs.data(), Hs.data()+Hs.num_elements()) <<"\n";

    double taux_min  = *std::min_element(tau_x.begin(), tau_x.end());
    std::cout<<"taux_min= "
             <<std::setprecision(10)<< taux_min <<"\n";
    double taux_max  = *std::max_element(tau_x.begin(), tau_x.end());
    std::cout<<"taux_max= "
             <<std::setprecision(10)<< taux_max <<"\n";

    // std::cout<<"------------------------------------------------------\n";
    // std::cout<<"dfloe_max= "<< *std::max_element(dfloe.data(), dfloe.data()+dfloe.num_elements()) <<"\n";
    // std::cout<<"dfloe_min= "<< *std::min_element(dfloe.data(), dfloe.data()+dfloe.num_elements()) <<"\n";
}

template<typename T>
void WimDiscr<T>::run(std::vector<value_type> const& ice_c, std::vector<value_type> const& ice_h, std::vector<value_type> const& n_floes, bool step)
{
    std::cout << "-----------------------Simulation started on "<< current_time_local() <<"\n";

    //init_time_str is human readable time eg "2015-01-01 00:00:00"
    //init_time is "formal" time format eg "20150101T000000Z"
    init_time_str = vm["wim.initialtime"].as<std::string>();
    std::string init_time = ptime(init_time_str);
    std::cout<<"---------------INITIAL TIME= "<< init_time <<"\n";

    this->assign(ice_c,ice_h,n_floes,step);

    bool critter;

    std::cout<<"Running starts\n";
    chrono.restart();

    std::cout<<"duration= "<< duration <<"\n";
    std::cout<<"amax= "<< amax <<"\n";
    std::cout<<"dt= "<< dt <<"\n";
    std::cout<<"nt= "<< nt <<"\n";


    //if (vm["wim.checkinit"].template as<bool>())
        exportResults("init",0);

    if ((!vm["nextwim.docoupling"].template as<bool>()) || (!step))
        cpt = 0;

    while (cpt < nt)
    {
        std::cout <<  ":[WIM2D TIME STEP]^"<< cpt+1 <<"\n";
        value_type t_out = dt*cpt;

        critter = !(cpt % vm["wim.dumpfreq"].template as<int>())
           && (vm["wim.checkprog"].template as<bool>());
        dumpDiag  = !(cpt % vm["wim.dumpfreq"].template as<int>())
           && (wim_itest>0) && (wim_jtest>0);

        if ( critter )
        {
            if (vm["nextwim.exportresults"].template as <bool>())
                exportResults("prog",t_out);//global counter from nextsim
           //else
           //exportResults(cpt,t_out);//local counter from wim
        }

        // if (cpt == 30)
        // {
        //     for (int i = 0; i < nx; i++)
        //         for (int j = 0; j < ny; j++)
        //             std::cout<<"Hs["<< i << "," << j << "]= "<< Tp[i][j] /*- Tp.data()[ny*i+j]*/ <<"\n";

        // }

        timeStep(step);

        ++cpt;//local counter incremented in wim.run()
    }

    //if (vm["wim.checkfinal"].template as<bool>())
       exportResults("out",nt*dt);

    // save diagnostic file
    save_log(nt*dt);

    std::cout<<"Running done in "<< chrono.elapsed() <<"s\n";

    std::cout << "-----------------------Simulation completed on "<< current_time_local() <<"\n";
}

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
}

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

    if (dmax<=dmin)
    {
       // small floes
       dave_ = std::pow(dmin,moment);
    }
    else
    {
       // bigger floes
       A     = (fsd_exp*std::exp(fsd_exp*(std::log(dmin)+std::log(dmax))));
       A     = A/(std::exp(fsd_exp*std::log(dmax))-std::exp(fsd_exp*std::log(dmin)));
       dave_ = -(A/b)*(std::exp(b*std::log(dmin))-exp(b*std::log(dmax)));
    }
}

template<typename T>
void WimDiscr<T>::advAttenSimple(array3_type& Sdir, array2_type& Sfreq, array2_type& taux_omega, array2_type& tauy_omega, array2_type const& ag2d_eff)
{
    array2_type uwave, vwave, temp;
    uwave.resize(boost::extents[nx][ny]);
    vwave.resize(boost::extents[nx][ny]);
    temp.resize(boost::extents[nx][ny]);

	std::vector<value_type> wt_theta(nwavedirn);
	value_type adv_dir, S_th, tmp, alp_dim, source;

    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

	for (int k = 0; k < nwavedirn; k++)
    {
        adv_dir = -PI*(90.0+wavedir[k])/180.0;
        uwave = ag2d_eff;
        std::for_each(uwave.data(), uwave.data()+uwave.num_elements(), [&](value_type& f){ f *= std::cos(adv_dir); });

        if (advdim == 2)
        {
            vwave = ag2d_eff;
            std::for_each(vwave.data(), vwave.data()+vwave.num_elements(), [&](value_type& f){ f *= std::sin(adv_dir); });
        }

        // copy from 3D input array to 2D temporary array
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                temp[i][j] = Sdir[i][j][k];
            }
        }

        waveAdvWeno(temp,uwave,vwave);

        // copy from 2D temporary array to 3D input array
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                Sdir[i][j][k] = temp[i][j];
            }
        }
    }

    if (nwavedirn == 1)
        std::fill( wt_theta.begin(), wt_theta.end(), 1. );
    else
        std::fill( wt_theta.begin(), wt_theta.end(), 2*PI/(1.0*nwavedirn) );


    std::fill( Sfreq.data(), Sfreq.data() + Sfreq.num_elements(), 0. );
    std::fill( taux_omega.data(), taux_omega.data() + taux_omega.num_elements(), 0. );
    std::fill( tauy_omega.data(), tauy_omega.data() + tauy_omega.num_elements(), 0. );

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            value_type adv_dir, S_th, tmp, alp_dim, source;

            if (ice_mask[i][j] > 0.)
            {
                for (int wnd = 0; wnd < nwavedirn; wnd++)
                {
                    adv_dir = -PI*(90.0+wavedir[wnd])/180.0;
                    S_th = Sdir[i][j][wnd];
                    alp_dim = atten_dim[i][j]+damp_dim[i][j];

                    // stress calculation
                    source = -alp_dim*ag2d_eff[i][j]*S_th;
                    tmp = -std::cos(adv_dir)*wt_theta[wnd]*source;
                    taux_omega[i][j] += tmp;
                    tmp = -std::sin(adv_dir)*wt_theta[wnd]*source;
                    tauy_omega[i][j] += tmp;

                    // do attenuation
                    Sdir[i][j][wnd] = S_th*std::exp(-alp_dim*ag2d_eff[i][j]*dt);

                    //std::cout<<"tau_x["<< i << "," << j << "]= "<< atten_dim[i][j] <<"\n";
                }
            }
#if 0
            // integrate spectrum over direction
            for (int wnd = 0; wnd < nwavedirn; wnd++)
            {
                Sfreq[i][j] += wt_theta[wnd]*Sdir[i][j][wnd];
            }
            //std::cout<<"taux_om["<< i << "," << j << "]= "<< taux_om[i][j] <<"\n";
#endif
        }
    }

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            // integrate spectrum over direction
            for (int wnd = 0; wnd < nwavedirn; wnd++)
            {
                Sfreq[i][j] += wt_theta[wnd]*Sdir[i][j][wnd];
            }
            //std::cout<<"taux_om["<< i << "," << j << "]= "<< taux_om[i][j] <<"\n";
        }
    }
}


template<typename T>
void WimDiscr<T>::advAttenIsotropic(array3_type& Sdir, array2_type& Sfreq, array2_type& taux_omega, array2_type& tauy_omega, array2_type const& ag2d_eff)
{
    array2_type uwave, vwave, temp;
    uwave.resize(boost::extents[nx][ny]);
    vwave.resize(boost::extents[nx][ny]);
    temp.resize(boost::extents[nx][ny]);

    std::vector<value_type> nvec(nwavedirn);
	std::vector<value_type> K_fou(nwavedirn), S_th(nwavedirn), theta_vec(nwavedirn), wt_theta(nwavedirn);
	std::vector<value_type> tmp1(nwavedirn), evals_x(nwavedirn);
	value_type adv_dir, tmp, alp_dim, source;

	std::vector<std::complex<value_type> > S_fou(nwavedirn);
	std::complex<value_type> zi, src_fou_p1, src_fou_m1;

	std::vector<value_type> S_cos(ncs), S_sin(ncs);
	value_type cg, q_scat, q_abs, q_tot, src_cos_1, src_sin_1;
    int n, jp1, jm1;

	zi = std::complex<value_type>(0.,1.);

	std::fill( wt_theta.begin(), wt_theta.end(), 2*PI/(1.0*nwavedirn) );

	for (int k = 0; k < nwavedirn; k++)
    {
        adv_dir = -PI*(90.0+wavedir[k])/180.0;

        uwave = ag2d_eff;

        std::for_each(uwave.data(), uwave.data()+uwave.num_elements(), [&](value_type& f){ f *= std::cos(adv_dir); });

        if (advdim == 2)
        {
	        vwave = ag2d_eff;
	        std::for_each(vwave.data(), vwave.data()+vwave.num_elements(), [&](value_type& f){ f *= std::sin(adv_dir); });
        }

        // copy from 3D input array to 2D temporary array
        for (int i = 0; i < nx; i++)
        {
	        for (int j = 0; j < ny; j++)
	        {
		        temp[i][j] = Sdir[i][j][k];
	        }
        }

        // advection
        waveAdvWeno(temp,uwave,vwave);

        // copy from 2D temporary array to 3D input array
        for (int i = 0; i < nx; i++)
        {
	        for (int j = 0; j < ny; j++)
	        {
		        Sdir[i][j][k] = temp[i][j];
	        }
        }

        theta_vec[k] = adv_dir;

        //std::cout<<"theta_vec["<< k <<"]= "<< std::setprecision(9) << std::sin((k+1)*theta_vec[k]) <<"\n";

        nvec[k] = (value_type)k;
    }

    std::fill( Sfreq.data(), Sfreq.data() + Sfreq.num_elements(), 0. );
    std::fill( taux_omega.data(), taux_omega.data() + taux_omega.num_elements(), 0. );
    std::fill( tauy_omega.data(), tauy_omega.data() + tauy_omega.num_elements(), 0. );

    int cpt=0;

	for (int i = 0; i < nx; i++)
    {
	    for (int j = 0; j < ny; j++)
	    {
		    for (int k = 0; k < nwavedirn; k++)
			    S_th[k] = Sdir[i][j][k];

		    std::fill( S_fou.begin(), S_fou.end(), zi );

		    //S_fou[0] = std::complex<value_type>( sum(wt_theta*S_th) );
            S_fou[0] = std::complex<value_type>( std::inner_product(wt_theta.begin(), wt_theta.end(), S_th.begin(), 0.) );

            if (ice_mask[i][j] > 0.)
            {
                ++cpt;

                if (dfloe[ny*i+j] < dfloe_pack_init)
                {
                    q_scat = atten_dim[i][j];
                    q_abs = damp_dim[i][j];
                }
                else
                {
                    q_scat = 0;
                    q_abs = atten_dim[i][j]+damp_dim[i][j];
                }

                q_tot = q_scat+q_abs;
                cg = ag2d_eff[i][j];

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

                taux_omega[i][j] = -src_cos_1;
                tauy_omega[i][j] = -src_sin_1;

                // if (cpt==1)
                // {
                //     std::cout<<"taux_omega= "<< taux_omega[i][j] <<"\n";
                //     std::cout<<"tauy_omega= "<< tauy_omega[i][j] <<"\n";
                // }

                std::vector<value_type> prodtmp = evals_x;
                std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::exp(cg*dt*f); });
                std::transform(S_fou.begin(), S_fou.end(), prodtmp.begin(), S_fou.begin(),
                               std::multiplies<std::complex<value_type> >());

                std::vector<value_type> Sfoutempcos(S_fou.size());// = S_fou;
                std::vector<value_type> Sfoutempsin(S_fou.size());// = S_fou;

                for (int k = 0; k < nwavedirn; k++)
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

                    std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::cos(theta_vec[k]*f); });

                    std::transform(Sfoutempcos.begin(), Sfoutempcos.end(), prodtmp.begin(), Sfoutempcos.begin(),
                                   std::multiplies<value_type>());


                    // prodtmp = std::vector<value_type>(nvec.begin(), nvec.end()); // also works
                    prodtmp.clear();
                    prodtmp.resize(nwavedirn);
                    prodtmp = nvec;

                    std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::sin(theta_vec[k]*f); });

                    std::transform(Sfoutempsin.begin(), Sfoutempsin.end(), prodtmp.begin(), Sfoutempsin.begin(),
                                   std::multiplies<value_type>());

                    for (int kin = 0; kin < nwavedirn; kin++)
                    {
                        tmp1[kin] = Sfoutempcos[kin]-Sfoutempsin[kin];
                    }

                    Sdir[i][j][k] = std::accumulate(tmp1.begin(), tmp1.end(), 0.0)/(2*PI);

                    // std::cout<<"std::numeric_limits<short>::max()=  "<< std::numeric_limits<float>::min() <<"\n";
                    S_th[k] = Sdir[i][j][k];
                }
            }

            Sfreq[i][j] = std::real(S_fou[0]);
	    }
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
}

template<typename T>
void WimDiscr<T>::waveAdvWeno(array2_type& h, array2_type const& u, array2_type const& v)
{
    array2_type sao;
    array2_type u_pad, v_pad, scp2_pad, scp2i_pad, scuy_pad, scvx_pad, h_pad;
    sao.resize(boost::extents[nxext][nyext]);

    array2_type hp_temp;
    hp_temp.resize(boost::extents[nx][ny]);

    padVar(u, u_pad);
    padVar(v, v_pad);
    padVar(SCP2_array, scp2_pad);
    padVar(SCP2I_array, scp2i_pad);
    padVar(SCUY_array, scuy_pad);
    padVar(SCVX_array, scvx_pad);
    padVar(h, h_pad);
      // TODO in fortran/matlab code,
      // h is treated differently from u,v etc
      // - they are always doubly periodic BC's
      // - may need to check with other BS's

    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    // prediction step
    weno3pdV2(h_pad, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad, sao);

    if (nghost==3)
    {
       // if only using nghost==3, need to apply boundary conditions
       // before correction step
#pragma omp parallel for num_threads(max_threads) collapse(2)
       for (int i = 0; i < nx; i++)
       {
           for (int j = 0; j < ny; j++)
           {
               hp[i+nbdx][j+nbdy] = h_pad[i+nbdx][j+nbdy]+dt*sao[i+nbdx][j+nbdy];
               hp_temp[i][j]      = hp[i+nbdx][j+nbdy];
           }
       }

       padVar(hp_temp,hp);//apply boundary conditions
    }
    else if (nghost>3)
    {
       // if using nghost>3, don't need to apply boundary conditions
       //  before correction step,
       // but need to loop over full padded domain

       //std::cout<<"not using padVar\n";
#pragma omp parallel for num_threads(max_threads) collapse(2)
       for (int i = 0; i < nxext; i++)
       {
           for (int j = 0; j < nyext; j++)
           {
               hp[i][j] = h_pad[i][j]+dt*sao[i][j];
           }
       }
    }
    else
    {
       std::cerr<<std::endl
                <<"Advection (WENO): 'nghost' should be >=3"
                <<std::endl;
       std::abort();
    }

    // correction step
    weno3pdV2(hp, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad, sao);


    //final output
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            h[i][j] = 0.5*(h_pad[i+nbdx][j+nbdy]+hp[i+nbdx][j+nbdy]+dt*sao[i+nbdx][j+nbdy]);

            //mask land cells
            h[i][j] *= 1-LANDMASK_array[i][j];
        }
    }
}

template<typename T>
void WimDiscr<T>::weno3pdV2(array2_type const& gin, array2_type const& u, array2_type const& v, array2_type const& scuy,
                       array2_type const& scvx, array2_type const& scp2i, array2_type const& scp2, array2_type& saoout)
{

	value_type cq00=-1./2 ,cq01=3./2, cq10=1./2, cq11=1./2, ca0=1./3, ca1=2./3, eps=1e-12;
	value_type q0, q1, a0, a1, q;
	int im1, im2, ip1, jm1, jm2, jp1, ymargin;

    array2_type ful, fuh, fvl, fvh, gt;
    ful.resize(boost::extents[nxext][nyext]);
    fuh.resize(boost::extents[nxext][nyext]);
    fvl.resize(boost::extents[nxext][nyext]);
    fvh.resize(boost::extents[nxext][nyext]);
    gt.resize(boost::extents[nxext][nyext]);

    if (advdim == 2)
        ymargin = 1;
    else
        ymargin = 0;

    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    // fluxes in x direction
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 2; i < nxext-1; i++)
    {
        for (int j = 0; j < nyext; j++)
        {
            value_type q0, q1, a0, a1, q;
            int im1, im2, ip1, jm1, jm2, jp1, ymargin;

            im1 = i-1;

            if (u[i][j] > 0.)
            {
                // coefficents to calc higher-order fluxes
                im2 = im1-1;
                q0 = cq00*gin[im2][j]+cq01*gin[im1][j];
                q1 = cq10*gin[im1][j]+cq11*gin[i][j];
                a0 = ca0;
                a1 = ca1*(std::abs(gin[im2][j]-gin[im1][j])+eps)/(std::abs(gin[im1][j]-gin[i][j])+eps);

                // lower-order fluxes
                ful[i][j] = u[i][j]*gin[im1][j]*scuy[i][j];

            }
            else
            {
                // coefficents to calc higher-order fluxes
                ip1 = i+1;
                q0 = cq11*gin[im1][j]+cq10*gin[i][j];
                q1 = cq01*gin[i][j]+cq00*gin[ip1][j];
                a0 = ca1;
                a1 = ca0*(abs(gin[im1][j]-gin[i][j])+eps)/(abs(gin[i][j]-gin[ip1][j])+eps);

                // lower-order fluxes
                ful[i][j] = u[i][j]*gin[i][j]*scuy[i][j];
            }

            // higher-order fluxes
            fuh[i][j] = (u[i][j]*(a0*q0+a1*q1)*scuy[i][j]/(a0+a1))-ful[i][j];
        }
    }

    // fluxes in y direction
    if (advdim == 2)
    {
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nxext; i++)
        {
            for (int j = 2; j < nyext-1; j++)
            {
                value_type q0, q1, a0, a1, q;
                int im1, im2, ip1, jm1, jm2, jp1, ymargin;

                jm1 = j-1;

                if (v[i][j] > 0.)
                {
                    jm2 = jm1-1;
                    q0 = cq00*gin[i][jm2]+cq01*gin[i][jm1];
                    q1 = cq10*gin[i][jm1]+cq11*gin[i][j];
                    a0 = ca0;
                    a1 = ca1*(std::abs(gin[i][jm2]-gin[i][jm1])+eps)/(std::abs(gin[i][jm1]-gin[i][j])+eps);
                    fvl[i][j] = v[i][j]*gin[i][jm1]*scvx[i][j];
                }
                else
                {
                    jp1 = j+1;
                    q0 = cq11*gin[i][jm1]+cq10*gin[i][j];
                    q1 = cq01*gin[i][j]+cq00*gin[i][jp1];
                    a0 = ca1;
                    a1 = ca0*(abs(gin[i][jm1]-gin[i][j])+eps)/(abs(gin[i][j]-gin[i][jp1])+eps);
                    fvl[i][j] = v[i][j]*gin[i][j]*scvx[i][j];
                }

                fvh[i][j] = (v[i][j]*(a0*q0+a1*q1)*scvx[i][j]/(a0+a1))-fvl[i][j];
            }
        }
    }

    // update field with low order fluxes
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nxext-1; i++)
    {
        for (int j = 0; j < nyext-ymargin; j++)//nyext-1 if 2d advection; else nyext
        {
            if (advdim == 2)
            {
                gt[i][j] = gin[i][j]-dt*(ful[i+1][j]-ful[i][j]+fvl[i][j+1]-fvl[i][j])*scp2i[i][j];
            }
            else if (advdim == 1)
            {
                gt[i][j] = gin[i][j]-dt*(ful[i+1][j]-ful[i][j])*scp2i[i][j];
            }
        }
    }

    q = 0.25/dt;

    // // obtain fluxes with limited high order correction fluxes
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 1; i < nxext; i++)
    {
        for (int j = 0; j < nyext; j++)
        {
            fuh[i][j] = ful[i][j]+
               std::max(-q*gt[i][j]*scp2[i][j],
                        std::min(q*gt[i-1][j]*scp2[i-1][j],fuh[i][j]));
        }
    }

    // obtain fluxes with limited high order correction fluxes
    if (advdim == 2)
    {
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nxext; i++)
        {
            for (int j = 1; j < nyext; j++)
            {
                fvh[i][j]=fvl[i][j]+
                   std::max(-q*gt[i][j]*scp2[i][j],
                            std::min(q*gt[i][j-1]*scp2[i][j-1],fvh[i][j]));
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
            if (advdim == 2)
            {
                saoout[i][j] = -(fuh[i+1][j]-fuh[i][j]+fvh[i][j+1]-fvh[i][j])*scp2i[i][j];
            }
            else if (advdim == 1)
            {
                saoout[i][j] = -(fuh[i+1][j]-fuh[i][j])*scp2i[i][j];
            }
        }
    }
#endif

}

template<typename T>
void WimDiscr<T>::padVar(array2_type const& u, array2_type& upad)
{
    upad.resize(boost::extents[nxext][nyext]);

    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nxext; i++)
    {
        for (int j = 0; j < nyext; j++)
        {

            if ((nbdx-1 < i) && (i < nx+nbdx) && (nbdy-1 < j) && (j < ny+nbdy))
            {
                upad[i][j] = u[i-nbdx][j-nbdy];
            }

            if (advdim == 1)
            {
                if (advopt != "notperiodic")
                {
                   // make periodic in i
                   if ((i < nbdx) && (nbdy-1 < j) && (j < ny+nbdy))
                   {
                       upad[i][j] = u[nx-nbdx+i][j-nbdy];
                   }

                   if ((nx+nbdx-1 < i) && (nbdy-1 < j) && (j < ny+nbdy))
                   {
                       upad[i][j] = u[i-nx-nbdx][j-nbdy];
                   }
                }
            }
            else if (advdim == 2)
            {
                if (advopt != "notperiodic")
                {
                    // make periodic in j
                    // if ((j < nbdy) && (nbdx-1 < i) && (i < nx+nbdx))
                    //     upad[i][j] = u[i-nbdx][ny-nbdy+j];

                    if ((j < nbdy) && (nbdx-1 < i) && (i < nx+nbdx))
                        upad[i][j] = u[i-nbdx][ny-nbdy+j-1];

                    if ((ny+nbdy-1 < j) && (nbdx-1 < i) && (i < nx+nbdx))
                        upad[i][j] = u[i-nbdx][j-ny-nbdy];
                }

                if (advopt == "xy-periodic")
                {
                    // make periodic in i
                    if ((i < nbdx) && (nbdy-1 < j) && (j < ny+nbdy))
                        upad[i][j] = u[nx-nbdx+i][j-nbdy];


                    if ((nx+nbdx-1 < i) && (nbdy-1 < j) && (j < ny+nbdy))
                        upad[i][j] = u[i-nx-nbdx][j-nbdy];


                    // // make periodic in j
                    // if ((j < nbdy) && (nbdy-1 < i) && (i < nx+nbdy))
                    //     upad[i][j] = u[i-nbdy][ny-nbdy+j];


                    if ((ny+nbdy-1 < j) && (nbdx-1 < i) && (i < nx+nbdx))
                        upad[i][j] = u[i-nbdx][j-ny-nbdy];


                    // BR, TL
                    if ((nx+nbdx-1 < i) && (ny+nbdy-1 < j))
                        upad[i][j] = u[i-nx-nbdx][j-ny-nbdy];


                    if ((i < nbdx) && (j < nbdy))
                        upad[i][j] = u[i+nx-nbdx][j];


                    // BL, TR
                    if ((nx+nbdx-1 < i) && (j < nbdy))
                        upad[i][j] = u[i-nx-nbdx][ny-nbdy+j];


                    if ((i < nbdx) && (ny+nbdy-1 < j))
                        upad[i][j] = u[i+nx-nbdx][j-ny-nbdy];

                }
            }
        }
    }
}

template<typename T>
void WimDiscr<T>::calcMWD()
{
    value_type adv_dir, wt_theta, om, CF;
    array2_type cmom0,cmom_dir,CSfreq, cmom_dir0;
    cmom0.resize(boost::extents[nx][ny]);
    cmom_dir.resize(boost::extents[nx][ny]);
    CSfreq.resize(boost::extents[nx][ny]);
    cmom_dir0.resize(boost::extents[nx][ny]);

    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    if (nwavedirn == 1)
        wt_theta = 1.;
    else
        wt_theta = 2.0*PI/(nwavedirn);

    // spectral moments
    std::fill( cmom0.data(), cmom0.data() + cmom0.num_elements(), 0. );
    std::fill( cmom_dir.data(), cmom_dir.data() + cmom_dir.num_elements(), 0. );

    for (int fq = 0; fq < nwavefreq; fq++)
    {
        om = 2*PI*freq_vec[fq];

        std::fill( CSfreq.data(), CSfreq.data() + CSfreq.num_elements(), 0. );
        std::fill( cmom_dir0.data(), cmom_dir0.data() + cmom_dir0.num_elements(), 0. );

        for (int dn = 0; dn < nwavedirn; dn++)
        {
            adv_dir = -PI*(90.0+wavedir[dn])/180.0;

#pragma omp parallel for num_threads(max_threads) collapse(2)
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    CSfreq[i][j] += wt_theta*sdf_dir[i][j][dn][fq];
                    cmom_dir0[i][j] += wt_theta*sdf_dir[i][j][dn][fq]*adv_dir;
                }
            }
        }//end loop over directions (still in freq loop)

#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if (ref_Hs_ice)
                    CF = std::pow(disp_ratio[i][j][fq],2);
                else
                    CF = 1.;

                cmom0[i][j]    += std::abs(wt_om[fq]*CF*CSfreq[i][j]);
                cmom_dir[i][j] += std::abs(wt_om[fq]*CF*cmom_dir0[i][j]);
            }
        }//end spatial loop i
    }//end freq loop

    //assign mwd
    std::fill( mwd.data(), mwd.data() + mwd.num_elements(), 0. );

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (cmom0[i][j] > 0.)
                mwd[i][j] = -90.-180*(cmom_dir[i][j]/cmom0[i][j])/PI;
        }
    }//end spatial loop i

}//end calcMWD

template<typename T>
typename WimDiscr<T>::value_type
WimDiscr<T>::thetaDirFrac(value_type const& th1_, value_type const& dtheta_, value_type const& mwd_)
{
    value_type mwd = thetaInRange(mwd_,th1_);   // th1_<=mwd<th1_+360
    value_type th2   = th1_ + dtheta_;

    if ((mwd > th2) && ((mwd - th2) > std::abs(mwd-360-th1_)))
    {
        mwd   = mwd - 360.;
    }

    value_type th1 = std::max<value_type>(mwd-90.,th1_);
    th2 = std::min<value_type>(mwd+90.,th2);
    th2 = std::max<value_type>(th1,th2); //!make th2>=th1

    value_type chi1 = PI*(th1-mwd)/180.;
    value_type chi2 = PI*(th2-mwd)/180.;

    value_type theta_dirfrac  = 2*(chi2-chi1)+std::sin(2*chi2)-std::sin(2*chi1);
    theta_dirfrac  = theta_dirfrac/(2*PI);

    return theta_dirfrac;
}

template<typename T>
typename WimDiscr<T>::value_type
WimDiscr<T>::thetaInRange(value_type const& th_, value_type const& th1)
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

    return th;
}

template<typename T>
void WimDiscr<T>::readDataFromFile(std::string const& filein)
{
    Fdmax.resize(boost::extents[nx][ny]);
    Ftaux.resize(boost::extents[nx][ny]);
    Ftauy.resize(boost::extents[nx][ny]);
    Fhs.resize(boost::extents[nx][ny]);
    Ftp.resize(boost::extents[nx][ny]);

    char * senv = ::getenv( "WIM2D_PATH" );
    std::string str = std::string( senv ) + "/fortran/run/out/binaries/prog";
    fs::path path(str);

    std::string _filein = (boost::format("%1%/%2%") % path.string() % filein).str();

    // s::path path(str);
    // path /= "outputs/binaries/prog";

    std::fstream in(_filein, std::ios::binary | std::ios::in);

    if (in.is_open())
    {
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Fdmax[i][j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Ftaux[i][j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Ftauy[i][j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Fhs[i][j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Ftp[i][j], sizeof(int));

        in.close();
    }
    else
    {
        std::cout << "Cannot open " << _filein  << "\n";
        std::cerr << "error: open file " << _filein << " for input failed!" <<"\n";
        std::abort();
    }
}

template<typename T>
void WimDiscr<T>::exportResults(std::string const& output_type,
                                 value_type const& t_out) const
{

    std::string str = vm["wim.outparentdir"].template as<std::string>();

    //char * senv = ::getenv( "WIM2D_PATH" );
    //if ( (str == ".") && (senv != NULL) && (senv[0] != '\0') )
    //{
    //    str = std::string( senv ) + "/CXX";
    //}

    fs::path path(str);
    std::string init_time  = ptime(init_time_str);
    std::string timestpstr = ptime(init_time_str, t_out);
    std::string fileout,fileoutb;
    if ( output_type == "prog" )
    {
       path   /= "binaries/prog";
       fileout = (boost::format( "%1%/wim_prog%2%" ) % path.string() % timestpstr).str();
    }
    else if ( output_type == "out" )
    {
       path   /= "binaries";
       fileout = (boost::format( "%1%/wim_out%2%" ) % path.string() % timestpstr).str();
    }
    else if ( output_type == "init" )
    {
       path   /= "binaries";
       fileout = (boost::format( "%1%/wim_init%2%" ) % path.string() % init_time).str();
    }

    fileoutb   = fileout+".b";
    fileout    = fileout+".a";
    if ( !fs::exists(path) )
        fs::create_directories(path);

    std::fstream out(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
    if (out.is_open())
    {
       if ( output_type == "init" )
       {
           for (int i = 0; i < icec.shape()[0]; i++)
               for (int j = 0; j < icec.shape()[1]; j++)
                   out.write((char *)&icec[i][j], sizeof(value_type));

           for (int i = 0; i < iceh.shape()[0]; i++)
               for (int j = 0; j < iceh.shape()[1]; j++)
                   out.write((char *)&iceh[i][j], sizeof(value_type));
       }

       for (int i = 0; i < nx; i++)
           for (int j = 0; j < ny; j++)
               out.write((char *)&dfloe[ny*i+j], sizeof(value_type));

       if ( output_type == "init" )
       {
          for (int i = 0; i < nx; i++)
              for (int j = 0; j < ny; j++)
                  out.write((char *)&tau_x[ny*i+j], sizeof(value_type));

          for (int i = 0; i < nx; i++)
              for (int j = 0; j < ny; j++)
                  out.write((char *)&tau_y[ny*i+j], sizeof(value_type));
       }


       for (int i = 0; i < Hs.shape()[0]; i++)
           for (int j = 0; j < Hs.shape()[1]; j++)
               out.write((char *)&Hs[i][j], sizeof(value_type));

       for (int i = 0; i < Tp.shape()[0]; i++)
           for (int j = 0; j < Tp.shape()[1]; j++)
               out.write((char *)&Tp[i][j], sizeof(value_type));

       for (int i = 0; i < mwd.shape()[0]; i++)
           for (int j = 0; j < mwd.shape()[1]; j++)
               out.write((char *)&mwd[i][j], sizeof(value_type));

       out.close();
    }
    else
    {
       std::cout << "Cannot open " << fileout  << "\n";
       std::cerr << "error: open file " << fileout << " for output failed!" <<"\n";
       std::abort();
    }

    // export the txt file for grid field information
    std::fstream outb(fileoutb, std::ios::out | std::ios::trunc);

    if (outb.is_open())
    {
        outb << std::setw(15) << std::left << "06"  << "    Nrecs    # "<< "Number of records" <<"\n";
        outb << std::setw(15) << std::left << "0"   << "    Norder   # "<< "Storage order [column-major (F/matlab) = 1; row-major (C) = 0]" <<"\n";
        outb << std::setw(15) << std::left << nx    << "    nx       # "<< "Record length in x direction (elements)" <<"\n";
        outb << std::setw(15) << std::left << ny    << "    ny       # "<< "Record length in y direction (elements)" <<"\n";

        outb <<"\n";
        outb << std::left << init_time << "    t_start    # "<< "Model time of WIM call" <<"\n";
        outb << std::left << timestpstr << "    t_out    # "<< "Model time of output" <<"\n";

        outb <<"\n";
        outb << "Record number and name:" <<"\n";

        int recno = 1;
        std::string rstr;

        if ( output_type == "init" )
        {
           rstr   = std::string(2-std::to_string(recno).length(),'0')
              + std::to_string(recno);
           outb << std::setw(9) << rstr << "icec" <<"\n";
           ++recno;
           //
           rstr   = std::string(2-std::to_string(recno).length(),'0')
              + std::to_string(recno);
           outb << std::setw(9) << rstr << "iceh" <<"\n";
           ++recno;
        }

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "Dmax" <<"\n";
        ++recno;

        if ( output_type != "init" )
        {
           rstr   = std::string(2-std::to_string(recno).length(),'0')
              + std::to_string(recno);
           outb << std::setw(9) << rstr << "tau_x" <<"\n";
           ++recno;
           //
           rstr   = std::string(2-std::to_string(recno).length(),'0')
              + std::to_string(recno);
           outb << std::setw(9) << rstr << "tau_y" <<"\n";
           ++recno;
        }

        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "Hs" <<"\n";
        ++recno;
        //
        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "Tp" <<"\n";
        ++recno;
        //
        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << "mwd" <<"\n";
        ++recno;
        outb.close();
    }
    else
    {
        std::cout << "Cannot open " << fileoutb  << "\n";
        std::cerr << "error: open file " << fileoutb << " for output failed!" <<"\n";
        std::abort();
    }
}

template<typename T>
void WimDiscr<T>::save_log(value_type const& t_out) const
{
    std::string str = vm["wim.outparentdir"].template as<std::string>();
    fs::path path(str);
    path /= "diagnostics/global";
    if ( !fs::exists(path) )
       fs::create_directories(path);

    std::string init_time  = ptime(init_time_str);
    std::string timestpstr = ptime(init_time_str, t_out);
    std::string fileout    = (boost::format( "%1%/WIMdiagnostics%2%.txt" ) % path.string() % init_time).str();

    std::fstream out(fileout, std::ios::out | std::ios::trunc);
    if ( !out.is_open() )
    {
        std::cout << "Cannot open " << fileout  << "\n";
        std::cerr << "error: open file " << fileout << " for output failed!" <<"\n";
        std::abort();
    }

    out << "***********************************************\n";
    out << "Outer subroutine:\n";
    out << ">> " << "wimdiscr.cpp\n\n";
    out << std::left << std::setw(32) << "Start time:  " << init_time << "\n";
    out << std::left << std::setw(32) << "Output time: " << timestpstr << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "Main parameters:" << "\n";
    out << std::left << std::setw(32) << "SCATMOD:" << scatmod << "\n";
    out << std::left << std::setw(32) << "ADV_DIM:" << advdim << "\n";
    out << std::left << std::setw(32) << "ADV_OPT:" << advopt << "\n";
#if 0
    //TODO implement brkopt
    out << std::left << std::setw(32) << "BRK_OPT:" << brkopt << "\n";
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
    out << std::left << std::setw(32) << "STEADY:" << steady << "\n";
    out << std::left << std::setw(32) << "DO_ATTEN:" << atten << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "Other integer parameters:" << "\n";
    out << std::left << std::setw(32) << "FSD_OPT:" << fsdopt << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "WIM parameters:" << "\n";
    out << std::left << std::setw(32) << "Brine volume fraction:" << vbf << "\n";
    out << std::left << std::setw(32) << "Youngs modulus (Pa):" << young << "\n";
    out << std::left << std::setw(32) << "Flexural strength (Pa):" << sigma_c << "\n";
//    out << std::left << std::setw(32) << "Breaking stress (Pa):" << stress_c << "\n";
    out << std::left << std::setw(32) << "Breaking strain:" << epsc << "\n";
    out << std::left << std::setw(32) << "Damping (Pa.s/m):" << visc_rp << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "Other parameters:" << "\n";
    out << std::left << std::setw(32) << "Time step (s):" << dt << "\n";
    out << std::left << std::setw(32) << "CFL number:" << cfl << "\n";
    out << std::left << std::setw(32) << "Max wave group vel (m/s):" << amax << "\n";
    out << std::left << std::setw(32) << "Number of time steps:" << nt << "\n";
    out << std::left << std::setw(32) << "Time interval (h):" << duration/60.0/60.0 << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << std::left << std::setw(32) << "Grid dimensions:" << nx << ", " << ny << "\n";
    out << std::left << std::setw(32) << "Spatial resolution (km):" << dx/1.0e3 << ", " << dy/1.0e3 << "\n";
    out << std::left << std::setw(32) << "Extent of domain (km):" << nx*dx/1.e3 << ", " << ny*dy/1.e3 << "\n";
    out << std::left << std::setw(32) << "Minimum period (s):" << 1.0/freq_vec[nwavefreq-1] << "\n";
    out << std::left << std::setw(32) << "Maximum period (s):" << 1.0/freq_vec[0] << "\n";
    out << std::left << std::setw(32) << "Number of wave frequencies:" << nwavefreq << "\n";
    out << std::left << std::setw(32) << "Number of wave directions:"  << nwavedirn << "\n";
    out << std::left << std::setw(32) << "Directional resolution (deg):" << 360.0/nwavedirn << "\n";
    out << "***********************************************\n";

    value_type taux_min  = *std::min_element(tau_x.begin(), tau_x.end());
    value_type taux_max  = *std::max_element(tau_x.begin(), tau_x.end());
    value_type tauy_min  = *std::min_element(tau_y.begin(), tau_y.end());
    value_type tauy_max  = *std::max_element(tau_y.begin(), tau_y.end());

    //MIZ diagnostics
    value_type Dmax_min = 10.e3;
    value_type Dmax_max = 0.e3;
    int Nmiz   = 0;

    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
#pragma omp parallel for num_threads(max_threads) collapse(1)
    for (int j = 0; j < dfloe.size(); j++)
    {
       if ((dfloe[j]<dfloe_pack_init)&&(dfloe[j]>0))
       {
          ++Nmiz;
          Dmax_min   = std::min(dfloe[j],Dmax_min);
          Dmax_max   = std::max(dfloe[j],Dmax_max);
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
    out << std::left << std::setw(32) << "MIZ width (km):"         << W_miz/1.e3 << "\n";
#endif
    out << std::left << std::setw(32) << "Dmax range in MIZ (m):"  << Dmax_min << ", " << Dmax_max << "\n";
    out << std::left << std::setw(32) << "tau_x range (Pa):"       << taux_min << ", " << taux_max << "\n";
    out << std::left << std::setw(32) << "tau_y range (Pa):"       << tauy_min << ", " << tauy_max << "\n";
    out << "***********************************************\n";

    out.close();
}

// instantiate wim class for type float
template class WimDiscr<float>;

// instantiate wim class for type double
template class WimDiscr<double>;

} // namespace WIM2D
