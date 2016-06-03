/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   wimdiscr.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug  3 11:52:35 2015
 */

#include <wimdiscr.hpp>

namespace Wim
{

template<typename T>
void WimDiscr<T>::gridProssessing()
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
        }
    }

    if (vm["wim.exportresults"].template as<bool>())
    {

        std::string str = vm["wim.outparentdir"].template as<std::string>();

        char * senv = ::getenv( "WIM2D_PATH" );
        if ( (str == ".") && (senv != NULL) && (senv[0] != '\0') )
        {
            str = std::string( senv ) + "/CXX";
        }

        fs::path path(str);
        path /= "outputs/binaries";

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

        std::string nxstr = std::string(3-std::to_string(nx).length(),'0') + std::to_string(nx);
        std::string nystr = std::string(3-std::to_string(ny).length(),'0') + std::to_string(ny);

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
	this->gridProssessing();
    //this->readGridFromFile("wim_grid.a");

    // parameters
	nwavedirn = vm["wim.nwavedirn"].template as<int>();
    nwavefreq = vm["wim.nwavefreq"].template as<int>();
    advdim = vm["wim.advdim"].template as<int>();
    ref_Hs_ice = vm["wim.refhsice"].template as<bool>();
    atten = vm["wim.atten"].template as<bool>();
    icevel = vm["wim.icevel"].template as<bool>();
    steady = vm["wim.steady"].template as<bool>();
    breaking = vm["wim.breaking"].template as<bool>();
    scatmod = vm["wim.scatmod"].template as<std::string>();
    advopt = vm["wim.advopt"].template as<std::string>();
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
    nbdx = vm["wim.nbdx"].template as<int>();
    nbdy = vm["wim.nbdy"].template as<int>();

    if (advdim == 1)
        nbdy = 0;

    nxext = nx+2*nbdx;
    nyext = ny+2*nbdy;

    gravity = 9.81;
    rhowtr = 1025.;
    rhoice = 922.5;
    poisson = 0.3;
    dmin = 20.;
    xi = 2.;
    fragility = 0.9;

    vbf = 0.1;
    vb = vbf;
    sigma_c  = (1.76e+6)*std::exp(-5.88*std::sqrt(vbf));
    epsc = sigma_c/young;
    flex_rig_coeff = young/(12.0*(1-std::pow(poisson,2.)));
}

template<typename T>
void WimDiscr<T>::assign(std::vector<value_type> const& ice_c, std::vector<value_type> const& ice_h, std::vector<value_type> const& n_floes) // reset
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

    sdf_dir.resize(boost::extents[nx][ny][nwavedirn][nwavefreq]);
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


    // wt_simp
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
    std::for_each(ag.begin(), ag.end(), [&](value_type& f){ f = f/2 ; });


    x0 = X_array[0][0];
    xmax = X_array[nx-1][ny-1]; //x0+(nx-1)*dx;

    y0 = Y_array[0][0];
    ym = Y_array[0][0]+(ny-1)*dy;
    //value_type ymax = Y_array[nx-1][ny-1];

    x_edge = 0.5*(x0+xmax)-0.8*(0.5*(xmax-x0));

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
            }
        }
    }

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
                for (int dn = 0; dn < nwavedirn; dn++)
                {
                    chi = PI*(wavedir[dn]-mwd[i][j])/180.0;
                    if (std::cos(chi) > 0.)
                        theta_fac[dn] = 2.0*std::pow(cos(chi),2.)/PI;
                    else
                        theta_fac[dn] = 0.;
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

                iceh[i][j] = 0.;

                if (ice_c[ny*i+j] > 0)
                    iceh[i][j] = ice_h[ny*i+j]/ice_c[ny*i+j];

                dfloe[ny*i+j] = 0.;

                if (n_floes[ny*i+j] > 0)
                    dfloe[ny*i+j] = std::sqrt(ice_c[ny*i+j]/n_floes[ny*i+j]);

                if (dfloe[ny*i+j] > dfloe_pack_thresh)
                    dfloe[ny*i+j] = dfloe_pack_init;


                if (icec[i][j] < vm["wim.cicemin"].template as<double>())
                    dfloe[ny*i+j] = 0.;

                if (icec[i][j] < vm["wim.cicemin"].template as<double>())
                    ice_mask[i][j] = 0.;
                else
                    ice_mask[i][j] = 1.;

                //std::cout<<"----------------------complex case\n";
            }

            //std::cout<<"ICE_MASK["<< i  << "," << j << "]= " << ice_mask[i][j] <<"\n";
        }
    }

    //std::cout<<"big loop starts\n";
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
                    wlng_ice[i][j][fq] = 2*PI/kice;

                    ag_eff[i][j][fq] = ag[fq];
                    ap_eff[i][j][fq] = ap[fq];

                    if (icevel)
                    {
                        wlng_ice[i][j][fq] = wlng[fq];
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

    dt = cfl*dx/amax;
    //ncs = std::ceil(nwavedirn/2);
    ncs = std::round(nwavedirn/2);
}

template<typename T>
void WimDiscr<T>::timeStep()
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
    value_type adv_dir, F, kicel, om, ommin, ommax, om1, lam1, lam2, dom, dave, c1d, tmp;
    int jcrest;
    bool break_criterion;

    dom = 2*PI*(freq_vec[nwavefreq-1]-freq_vec[0])/(nwavefreq-1);

    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    if (vm["wim.steady"].template as<bool>())
    {
        for (int i = 0; i < nwavefreq; i++)
        {
            for (int j = 0; j < nwavedirn; j++)
            {
                adv_dir = (-PI/180)*(wavedir[j]+90.);

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

    dave = 0;

    for (int fq = 0; fq < nwavefreq; fq++)
    {
        std::fill( atten_dim.data(), atten_dim.data() + atten_dim.num_elements(), 0. );
        std::fill( damp_dim.data(), damp_dim.data() + damp_dim.num_elements(), 0. );

#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                value_type dave, c1d;
                if ((ice_mask[i][j] == 1.) && (atten))
                {
                    if (dfloe[ny*i+j] <200.)
                        floeScaling(dfloe[ny*i+j],dave);
                    else
                        dave = dfloe[ny*i+j];

                    // floes per unit length
                    c1d = icec[i][j]/dave;

                    // scattering
                    atten_dim[i][j] = atten_nond[i][j][fq]*c1d;

                    // damping
                    damp_dim[i][j] = 2*damping[i][j][fq]*icec[i][j];
                }
            }
        }


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

    calcMWD();

    if (!(steady) && !(breaking))
    {
        auto temparray = Hs;
        std::for_each(temparray.data(), temparray.data()+temparray.num_elements(), [&](value_type& f){ f *= f; });
        E_tot = std::accumulate(temparray.data(), temparray.data()+temparray.num_elements(),0.);

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
            value_type adv_dir, F, kicel, om, ommin, ommax, om1, lam1, lam2, /*dom,*/ dave, c1d, tmp;
            int jcrest;
            bool break_criterion;

            //std::cout << "MASK[" << i << "," << j << "]= " << ice_mask[i][j] << " and "<< mom0[i][j]  <<"\n";

            if ((ice_mask[i][j] == 1.) && (mom0[i][j] >= 0.))
            {
                // significant strain amp
                sig_strain = 2*std::sqrt(var_strain[i][j]);

                // probability of critical strain
                // being exceeded from Rayleigh distribution
                Pstrain = std::exp( -std::pow(epsc,2.)/(2*var_strain[i][j]) );
                P_crit = std::exp(-1.0);

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
            }
            else if (wtr_mask[i][j] == 1.)
            {
                dfloe[ny*i+j] = 0;
            }

            nfloes[ny*i+j] = 0.;

            if (dfloe[ny*i+j] > 0.)
                nfloes[ny*i+j] = icec[i][j]/std::pow(dfloe[ny*i+j],2.);
        }
    }


    // for (int i = 0; i < nx; i++)
    //     for (int j = 0; j < ny; j++)
    //         std::cout << "Dmax[" << i << "," << j << "]= " << dfloe[i][j] <<"\n";

    // for (int i = 0; i < nx; i++)
    //     for (int j = 0; j < ny; j++)
    //         std::cout<<"Hs["<< i << "," << j << "]= "<< Hs[i][j] <<"\n";

    std::cout<<"Hs_max= "<< *std::max_element(Hs.data(), Hs.data()+Hs.num_elements()) <<"\n";
    std::cout<<"Hs_min= "<< *std::min_element(Hs.data(), Hs.data()+Hs.num_elements()) <<"\n";

    // std::cout<<"------------------------------------------------------\n";
    // std::cout<<"dfloe_max= "<< *std::max_element(dfloe.data(), dfloe.data()+dfloe.num_elements()) <<"\n";
    // std::cout<<"dfloe_min= "<< *std::min_element(dfloe.data(), dfloe.data()+dfloe.num_elements()) <<"\n";
}

template<typename T>
void WimDiscr<T>::run(std::vector<value_type> const& ice_c, std::vector<value_type> const& ice_h, std::vector<value_type> const& n_floes)
{
    this->assign(ice_c,ice_h,n_floes);

    value_type duration;
    int nt;
    bool critter;

    std::cout<<"Running starts\n";
    chrono.restart();

    //duration = 1.0e3*x_ext/u_ref;
    duration = vm["wim.duration"].template as<double>();
    std::cout<<"duration= "<< duration <<"\n";
    std::cout<<"dt= "<< dt <<"\n";

    //nt = std::floor(duration/dt);
    nt = std::round(duration/dt);
    std::cout<<"nt= "<< nt <<"\n";
    int cpt = 0;

    while (cpt < nt)
    {
        std::cout <<  ":[WIM2D TIME STEP]^"<< cpt+1 <<"\n";
        value_type t_out = dt*cpt;

        critter = !(cpt % vm["wim.reps"].template as<int>()) && (vm["wim.checkprog"].template as<bool>());

        if ((vm["wim.exportresults"].template as<bool>()) && (critter))
            exportResults(cpt,t_out);

        // if (cpt == 30)
        // {
        //     for (int i = 0; i < nx; i++)
        //         for (int j = 0; j < ny; j++)
        //             std::cout<<"Hs["<< i << "," << j << "]= "<< Tp[i][j] /*- Tp.data()[ny*i+j]*/ <<"\n";

        // }

        timeStep();

        ++cpt;
    }

    std::cout<<"Running done in "<< chrono.elapsed() <<"s\n";
}

template<typename T>
void WimDiscr<T>::floeScaling(value_type const& dmax, value_type& dave)
{
    int mm = 0;
    value_type n, nsum, nd, ndsum, r, dfac;

    r = dmax/dmin;
    n = nsum = nd = ndsum = dfac = 0;

    while (r > xi)
    {
        r  = r/xi;
        ++mm;
    }

    if (mm > 0)
    {
        for (int m=0; m<mm+1; ++m)
        {
            n = (1.0-fragility)*std::pow((fragility*std::pow(xi,2.)),m);
            nd = n/(std::pow(xi,m));
            nsum += n;
            ndsum += nd;
            dfac = ndsum/nsum;
        }

        dave = dfac*dmax;
    }
    else
    {
        dave = dmin;
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

    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    // prediction step
    weno3pdV2(h_pad, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad, sao);

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            hp[i+nbdx][j+nbdy] = h_pad[i+nbdx][j+nbdy]+dt*sao[i+nbdx][j+nbdy];
        }
    }

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            hp_temp[i][j] = hp[i+nbdx][j+nbdy];
        }
    }

    padVar(hp_temp,hp);

    // correction step
    weno3pdV2(hp, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad, sao);

#if 0
    for (int i = 0; i < nxext; i++)
    {
        for (int j = 0; j < nyext; j++)
        {
            h_pad[i][j] = 0.5*(h_pad[i][j]+hp[i][j]+dt*sao[i][j]);
        }
    }

    for (int i = nbdx; i < nx+nbdx; i++)
    {
        for (int j = nbdy; j < ny+nbdy; j++)
        {
            h[i-nbdx][j-nbdy] = h_pad[i][j];

            // mask land (no waves on land)
            h[i-nbdx][j-nbdy] = h[i-nbdx][j-nbdy]*(1-LANDMASK_array[i-nbdx][j-nbdy]);
        }
    }
#endif

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            h[i][j] = 0.5*(h_pad[i+nbdx][j+nbdy]+hp[i+nbdx][j+nbdy]+dt*sao[i+nbdx][j+nbdy]);
        }
    }

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            h[i][j] = h[i][j]*(1-LANDMASK_array[i][j]);
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

    // fluxes in x directional
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = nbdx-1; i < nx+nbdx+2; i++)
    {
        for (int j = nbdy-ymargin; j < ny+nbdy+ymargin; j++)
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
        for (int i = nbdx-1; i < nx+nbdx+1; i++)
        {
            for (int j = nbdy-1; j < ny+nbdy+2; j++)
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
    for (int i = nbdx-1; i < nx+nbdx+1; i++)
    {
        for (int j = nbdy-ymargin; j < ny+nbdy+ymargin; j++)
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
    for (int i = nbdx; i < nx+nbdx+1; i++)
    {
        for (int j = nbdy; j < ny+nbdy; j++)
        {
            fuh[i][j] = ful[i][j]+std::max(-q*gt[i][j]*scp2[i][j],std::min(q*gt[i-1][j]*scp2[i-1][j],fuh[i][j]));
        }
    }

    // obtain fluxes with limited high order correction fluxes
    if (advdim == 2)
    {
#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = nbdx; i < nx+nbdx; i++)
        {
            for (int j = nbdy; j < ny+nbdy+1; j++)
            {
                fvh[i][j]=fvl[i][j]+std::max(-q*gt[i][j]*scp2[i][j],std::min(q*gt[i][j-1]*scp2[i][j-1],fvh[i][j]));
            }
        }
    }

#if 1
    // compute the spatial advective operator
#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = nbdx; i < nx+nbdx; i++)
    {
        for (int j = nbdy; j < ny+nbdy; j++)
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
    value_type adv_dir, wt_theta, om;
    array2_type cmom0,cmom_dir,CSfreq, cmom_dir0, CF;
    cmom0.resize(boost::extents[nx][ny]);
    cmom_dir.resize(boost::extents[nx][ny]);
    CSfreq.resize(boost::extents[nx][ny]);
    cmom_dir0.resize(boost::extents[nx][ny]);
    CF.resize(boost::extents[nx][ny]);

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
        }

#pragma omp parallel for num_threads(max_threads) collapse(2)
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if (ref_Hs_ice)
                    CF[i][j] = disp_ratio[i][j][fq];
                else
                    CF[i][j] = 1;

                cmom0[i][j] += std::abs(wt_om[fq]*std::pow(CF[i][j],2.)*CSfreq[i][j]);
                cmom_dir[i][j] += std::abs(wt_om[fq]*std::pow(CF[i][j],2.)*cmom_dir0[i][j]);
            }
        }
    }

    std::fill( mwd.data(), mwd.data() + mwd.num_elements(), 0. );

#pragma omp parallel for num_threads(max_threads) collapse(2)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (cmom0[i][j] > 0.)
                mwd[i][j] = -90.-180*(cmom_dir[i][j]/cmom0[i][j])/PI;
        }
    }

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
void WimDiscr<T>::exportResults(size_type const& timestp, value_type const& t_out) const
{

    std::string str = vm["wim.outparentdir"].template as<std::string>();

    char * senv = ::getenv( "WIM2D_PATH" );
    if ( (str == ".") && (senv != NULL) && (senv[0] != '\0') )
    {
        str = std::string( senv ) + "/CXX";
    }

    fs::path path(str);
    path /= "outputs/binaries/prog";

    if ( !fs::exists(path) )
        fs::create_directories(path);

    std::string timestpstr = std::string(3-std::to_string(timestp).length(),'0') + std::to_string(timestp);

    std::string fileout = (boost::format( "%1%/wim_prog%2%.a" ) % path.string() % timestpstr).str();
    std::fstream out(fileout, std::ios::binary | std::ios::out | std::ios::trunc);

    if (out.is_open())
    {
        for (int i = 0; i < icec.shape()[0]; i++)
            for (int j = 0; j < icec.shape()[1]; j++)
                out.write((char *)&icec[i][j], sizeof(value_type));

        for (int i = 0; i < iceh.shape()[0]; i++)
            for (int j = 0; j < iceh.shape()[1]; j++)
                out.write((char *)&iceh[i][j], sizeof(value_type));
#if 0
        for (int i = 0; i < dfloe.shape()[0]; i++)
            for (int j = 0; j < dfloe.shape()[1]; j++)
                out.write((char *)&dfloe[i][j], sizeof(value_type));

        for (int i = 0; i < tau_x.shape()[0]; i++)
            for (int j = 0; j < tau_x.shape()[1]; j++)
                out.write((char *)&tau_x[i][j], sizeof(value_type));

        for (int i = 0; i < tau_y.shape()[0]; i++)
            for (int j = 0; j < tau_y.shape()[1]; j++)
                out.write((char *)&tau_y[i][j], sizeof(value_type));
#endif

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                out.write((char *)&dfloe[ny*i+j], sizeof(value_type));

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                out.write((char *)&tau_x[ny*i+j], sizeof(value_type));

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                out.write((char *)&tau_y[ny*i+j], sizeof(value_type));


        for (int i = 0; i < Hs.shape()[0]; i++)
            for (int j = 0; j < Hs.shape()[1]; j++)
                out.write((char *)&Hs[i][j], sizeof(value_type));

        for (int i = 0; i < Tp.shape()[0]; i++)
            for (int j = 0; j < Tp.shape()[1]; j++)
                out.write((char *)&Tp[i][j], sizeof(value_type));

        out.close();
    }
    else
    {
        std::cout << "Cannot open " << fileout  << "\n";
        std::cerr << "error: open file " << fileout << " for output failed!" <<"\n";
        std::abort();
    }

    // export the txt file for grid field information
    std::string fileoutb = (boost::format( "%1%/wim_prog%2%.b" ) % path.string() % timestpstr).str();
    std::fstream outb(fileoutb, std::ios::out | std::ios::trunc);

    if (outb.is_open())
    {
        outb << std::setw(15) << std::left << "07"  << "    Nrecs    # "<< "Number of records" <<"\n";
        outb << std::setw(15) << std::left << "0"   << "    Norder   # "<< "Storage order [column-major (F/matlab) = 1; row-major (C) = 0]" <<"\n";
        outb << std::setw(15) << std::left << nx    << "    nx       # "<< "Record length in x direction (elements)" <<"\n";
        outb << std::setw(15) << std::left << ny    << "    ny       # "<< "Record length in y direction (elements)" <<"\n";
        outb << std::setw(15) << std::left << t_out << "    t_out    # "<< "Model time of output (s)" <<"\n";
        //outb << std::setw(15) << std::left << nwavefreq << "          "<< "Number of wave frequencies" <<"\n";
        //outb << std::setw(15) << std::left << nwavedirn << "          "<< "Number of wave directions" <<"\n";

        outb <<"\n";

        outb << "Record number and name:" <<"\n";
        outb << std::setw(9) << std::left << "01" << "icec" <<"\n";
        outb << std::setw(9) << std::left << "02" << "iceh" <<"\n";
        outb << std::setw(9) << std::left << "03" << "Dmax" <<"\n";
        outb << std::setw(9) << std::left << "04" << "tau_x" <<"\n";
        outb << std::setw(9) << std::left << "05" << "tau_y" <<"\n";
        outb << std::setw(9) << std::left << "06" << "Hs" <<"\n";
        outb << std::setw(9) << std::left << "07" << "Tp" <<"\n";
    }
    else
    {
        std::cout << "Cannot open " << fileoutb  << "\n";
        std::cerr << "error: open file " << fileoutb << " for output failed!" <<"\n";
        std::abort();
    }
}

// instantiate wim class for type float
template class WimDiscr<float>;

// instantiate wim class for type double
template class WimDiscr<double>;

} // namespace WIM2D
