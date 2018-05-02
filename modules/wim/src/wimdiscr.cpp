/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   wimdiscr.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug  3 11:52:35 2015
 */

#include <wimdiscr.hpp>
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
WimDiscr<T>::WimDiscr(po::variables_map const& vmIn,int const& nextsim_cpt)
{
    vm  = vmIn;
    this->initConstant(nextsim_cpt);

    if(!M_wim_on_mesh)
    {
        // wim grid generation/reading
        // NB if M_wim_on_mesh, setMesh2 before wim.run() and at regridding
        // time
        M_grid = T_grid(vm);
    }

    this->initRemaining();
}//WimDiscr()


template<typename T>
WimDiscr<T>::WimDiscr(po::variables_map const& vmIn,T_gmsh const &mesh_in,int const& nextsim_cpt)
{
    vm  = vmIn;
    this->initConstant(nextsim_cpt);

    // init grid FROM mesh
    T_mesh mesh(mesh_in);//tmp mesh object
    M_grid  = T_grid(vm,mesh);

    this->initRemaining();
}//WimDiscr()


template<typename T>
void WimDiscr<T>::initConstant(int const& nextsim_cpt)
{

    //OMP max threads
    M_max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/

    //set global counter to 0
    M_cpt = 0;

    // wave parameters
    nwavedirn = vm["wimsetup.nwavedirn"].template as<int>();
    nwavefreq = vm["wimsetup.nwavefreq"].template as<int>();
    Tmin      = vm["wimsetup.tmin"].template as<double>(); /* 2.5 */
    Tmax      = vm["wimsetup.tmax"].template as<double>(); /* 25. */

    M_ref_Hs_ice = vm["wim.refhsice"].template as<bool>();
    M_atten      = vm["wim.atten"].template as<bool>();
    M_useicevel  = vm["wim.useicevel"].template as<bool>();
    M_steady     = vm["wim.steady"].template as<bool>();
    M_scatmod    = vm["wim.scatmod"].template as<std::string>();

    // ice parameters
    M_breaking        = vm["wim.breaking"].template as<bool>();
    M_dfloe_pack_init = vm["wim.dfloepackinit"].template as<double>(); /* 300.0 */
    M_ice_params      = T_icep(vm);

    //numerical parameters
    M_advdim = vm["wim.advdim"].template as<int>();
    M_advopt = vm["wim.advopt"].template as<std::string>();
    M_cfl    = vm["wim.cfl"].template as<double>();

    if (M_useicevel)
        throw std::runtime_error("M_useicevel=true not implemented\n");

    if(M_steady && (M_advopt=="xy-periodic"))
    {
        std::string tmps = "advopt = xy-periodic and steady options incompatible - ";
        tmps += "use y-periodic with steady or turn off steady\n";
        throw std::runtime_error(tmps);
    }

    //some options need to be disabled if being called from nextsim
    docoupling = !( vm.count("simul.duration")==0 );
    M_restart_time = 0.;
    if (!docoupling)
    {
        //set duration of call to wim.run() from wim.duration
        M_duration = vm["wimsetup.duration"].template as<double>();

        //get initial time from wimsetup.initialtime
        M_init_time_str = vm["wimsetup.initialtime"].as<std::string>();

        //save options from wimoptions.cpp
        this->saveOptionsLog();
    }
    else
    {
        //set duration of call to wim.run() from nextwim.couplingfreq
        T_val nextsim_time_step = vm["simul.timestep"].template as<double>();
        M_duration = vm["nextwim.couplingfreq"].template as<int>()*nextsim_time_step;
                     
        //get initial time from simul.time_init
        M_init_time_str = vm["simul.time_init"].template as<std::string>();

        //if using restart, may need to calculate shift from initial time
        M_restart_time = nextsim_cpt*nextsim_time_step;//model time of current call to wim
 
        if ( vm["nextwim.coupling-option"].template as<std::string>() == "run_on_mesh")
            //run on mesh takes priority
            M_wim_on_mesh = true;
        else if ( vm["nextwim.coupling-option"].template as<std::string>() == "break_on_mesh")
            M_break_on_mesh =  true;
    }

    M_current_time = M_restart_time;//reset at last update
    M_update_time  = M_restart_time;//reset at last update

#if 1
    //print initial & restart times
    auto time_str1 = ptime(M_init_time_str,M_restart_time);
    std::cout<<"initial time = "<<M_init_time_str<<"\n";
    std::cout<<"restart time = "<<time_str1<<"\n";
#endif
}//end ::initConstant()


template<typename T>
void WimDiscr<T>::initRemaining()
{
    //after grid/mesh are set
    if(!M_wim_on_mesh)
    {
        M_num_elements  = M_grid.M_num_elements;
        M_length_cfl    = M_grid.M_resolution;
    }

    // ==============================================================================
    //local diagnostics
    int wim_itest = vm["wimdiag.itest"].template as<int>();
    int wim_jtest = vm["wimdiag.jtest"].template as<int>();

    M_itest  = -1;
    if(!M_wim_on_mesh)
        if((wim_itest>=0)&&(wim_jtest>=0))
        {
            //(!M_wim_on_mesh) && positive itest,jtest input
            std::cout<<"\nitest,jtest: "<<wim_itest<<","<<wim_jtest<<"\n";
            M_itest  = wim_itest*M_grid.M_num_py+wim_jtest;
            if ((wim_itest>=M_grid.M_num_px)||(wim_jtest>=M_grid.M_num_py))
            {
                std::cout<<"\nitest,jtest: "<<wim_itest<<","<<wim_jtest<<"\n";
                std::cout<<"\nnx,ny: "<<M_grid.M_num_px<<","<<M_grid.M_num_py<<"\n";
                throw std::runtime_error("wimdiag.itest/jtest out of range");
            }
        }

    //global test index
    if (M_itest>=0)
        std::cout<<"\nWIM diagnostics to be done at global index: "
            <<M_itest<<"(of "<<M_num_elements<<")\n";
    if (M_itest>=M_num_elements)
        throw std::runtime_error("M_itest out of range");
    // ==============================================================================

    // call assign to set sizes of some arrays (not depending on space)
    this->assign();

    if(!M_wim_on_mesh)
    {
        // if(M_wim_on_mesh), assignSpatial() called in run()
        // - since mesh is changing each time
        M_num_elements  = M_grid.M_num_elements;
        M_land_mask     = M_grid.M_land_mask;
        this->assignSpatial();
    }

    std::cout<<"wim instantiation finished\n";
}//initRemaining()


template<typename T>
void WimDiscr<T>::assign()
{
    // this doesn't need to be called each time wim.run is called
    // - sets arrays with dimensions depending on frequency, wave dirn (never change)

    // initialise IceInfo types 1st
    // - NB don't set vectors yet
    std::cout<<"init M_ice["<<M_ice.size()<<"]: wim\n";
    M_ice.push_back(T_ice(M_ice_params,"wim"));
    if(!M_wim_on_mesh)
    {
        //if not on mesh we need a separate T_ice object for the
        //sea ice model
        std::cout<<"init M_ice["<<M_ice.size()<<"]: sim\n";
        M_ice.push_back(T_ice(M_ice_params,"sim"));
    }

    M_quadrature_wt_freq.resize(nwavefreq);
    M_freq_vec.resize(nwavefreq);

    M_wlng_wtr.resize(nwavefreq);
    M_ag_wtr.resize(nwavefreq);
    M_ap_wtr.resize(nwavefreq);

    //4d vec
    M_sdf_dir.resize(nwavefreq);

    // =============================================
    // set frequencies to use (M_freq_vec)
    // and quadrature weights to integrate wrt freq
    T_val tp_in = vm["wim.tpinc"].template as<double>();
    if (nwavefreq == 1)
    {
        M_freq_vec[0] = 1./tp_in;
        M_quadrature_wt_freq[0] = 1.;
    }
    else
    {
        // multiple frequencies
        T_val fmin = 1./Tmax;
        T_val fmax = 1./Tmin;
        T_val df = (fmax-fmin)/(nwavefreq-1);

        for (int fq = 0; fq < nwavefreq; fq++)
            M_freq_vec[fq] = fmin+fq*df;

        // weights for Simpson's rule
        // {1,4,2,4...,2,4,1}*dom/3.
        std::fill(M_quadrature_wt_freq.begin(), M_quadrature_wt_freq.end(), 2.);
        M_quadrature_wt_freq[0] = 1.;
        M_quadrature_wt_freq[nwavefreq-1] = 1.;

        int w = 1;
        while (w < nwavefreq-1)
        {
            M_quadrature_wt_freq[w] = 4.;
            w +=2;
        }

        T_val dom = 2*PI*(M_freq_vec[nwavefreq-1]-M_freq_vec[0])/(nwavefreq-1);
        std::for_each(M_quadrature_wt_freq.begin(), M_quadrature_wt_freq.end(), [&](T_val& f){ f = dom*f/3.0; });
    }
    // =============================================


    // =============================================
    // set directions to use (M_wavedir)
    T_val avgdir = vm["wim.mwdinc"].template as<double>();
    M_wavedir.assign(nwavedirn,avgdir);
    M_quadrature_wt_dir.assign(nwavedirn,1.);
    if (nwavedirn > 1)
    {
        T_val theta_max = 90.;
        T_val theta_min = -270.;
        T_val dtheta = (theta_min-theta_max)/nwavedirn;

        for (int nth = 0; nth < nwavedirn; nth++)
            M_wavedir[nth]    = theta_max+nth*dtheta;

        std::fill( M_quadrature_wt_dir.begin(), M_quadrature_wt_dir.end(), (2*PI)/nwavedirn );
    }
    // =============================================


    // =============================================
    // vector of periods, open water wavelength, phase and group velocities

    // open water wavelengths
    T_val gravity = M_ice_params.gravity();
    M_wlng_wtr = M_freq_vec;
    std::for_each(M_wlng_wtr.begin(), M_wlng_wtr.end(), [&](T_val& f){ f = gravity/(2*PI*std::pow(f,2.)); });

    // open water phase velocities
    M_ap_wtr = M_wlng_wtr;
    std::for_each(M_ap_wtr.begin(), M_ap_wtr.end(), [&](T_val& f){ f = std::sqrt(gravity*f/(2*PI)); });

    // open water group velocities
    M_ag_wtr = M_ap_wtr;
    std::for_each(M_ag_wtr.begin(), M_ag_wtr.end(), [&](T_val& f){ f = f/2. ; });
    // =============================================

    //if running on mesh, take these values as input from open boundaries
    //default is zero, unless using ideal waves and "steady" option
    T_val_vec ztmp(nwavedirn,0.);
    M_open_boundary_vals.assign(nwavefreq,ztmp);

}//end: assign()


template<typename T>
void WimDiscr<T>::assignSpatial()
{
    // this needs to be called each time grid or mesh changes
    // ie initially, and if(M_wim_on_mesh), after regridding
    // set sizes of arrays, initialises some others that are constant in time

    M_assigned      = true;
    M_num_elements  = this->getNumElements();

    //2D var's
    M_dave.assign(M_num_elements,0.);

    // NB this clears wave diagnostics
    // so take care to reset them after regrid
    // if running WIM on mesh
    M_Hs.assign(M_num_elements,0.);
    M_Tp.assign(M_num_elements,0.);
    M_mwd.assign(M_num_elements,0.);
    M_stokes_drift_x.assign(M_num_elements,0.);
    M_stokes_drift_y.assign(M_num_elements,0.);
    M_tau_x.assign(M_num_elements,0.);
    M_tau_y.assign(M_num_elements,0.);
    M_mwd_x.assign(M_num_elements,0.);
    M_mwd_y.assign(M_num_elements,0.);

    //these are only temporary vectors, but they are global in order to
    //save creating and destroying them extremely often
    Mtmp_sdf_freq.assign(M_num_elements,0.);
    Mtmp_stokes_drift_x_om.assign(M_num_elements,0.);
    Mtmp_stokes_drift_y_om.assign(M_num_elements,0.);
    Mtmp_mom0      .assign( M_num_elements, 0. );
    Mtmp_mom2      .assign( M_num_elements, 0. );
    Mtmp_var_strain.assign( M_num_elements, 0. );
    Mtmp_mom0w     .assign( M_num_elements, 0. );
    Mtmp_mom2w     .assign( M_num_elements, 0. );
    if(M_atten)
    {
        Mtmp_atten_dim.assign(M_num_elements,0.);
        Mtmp_damp_dim.assign(M_num_elements,0.);
        Mtmp_taux_om.assign(M_num_elements,0.);
        Mtmp_tauy_om.assign(M_num_elements,0.);
    }
    
    //3D var's
    // - space and freq
    T_val_vec ztmp(M_num_elements,0.);
    M_ag_eff.assign(nwavefreq,ztmp);
    M_agnod_eff.assign(nwavefreq,{});//if(M_wim_on_mesh), interp group vel from elements to nodes
    M_ap_eff.assign(nwavefreq,ztmp);
    M_wlng_ice.assign(nwavefreq,ztmp);
    M_disp_ratio.assign(nwavefreq,ztmp);
    if(M_atten)
    {
        M_atten_nond.assign(nwavefreq,ztmp);
        M_damping.assign(nwavefreq,ztmp);
    }

    if (M_sdf_dir[0].size()==0)
        //set dir spec to 0. on 1st call to init/assign
        //std::cout<<"Init M_sdf_dir in wim.assign()\n";
        for (auto it=M_sdf_dir.begin();it!=M_sdf_dir.end();it++)
            it->assign(nwavedirn,ztmp);
    // =============================================

}//end: assignSpatial()


template<typename T>
void WimDiscr<T>::update()
{

    //====================================================
    // update attenuation coefficients, wavelengths and phase/group velocities
    this->updateWaveMedium();
    //====================================================


    // ====================================================================================
    // set time step
    // - this can change with time if using group velocity for ice
    // - NB needs to be done after updateWaveMedium
    M_timestep = M_cfl*M_length_cfl/M_max_cg;

    //reduce time step slightly (if necessary) to make M_duration an integer multiple of M_timestep
    M_num_timesteps = std::ceil(M_duration/M_timestep);
    M_timestep      = M_duration/M_num_timesteps;
    //std::cout<<"M_timestep,M_num_timesteps= "<< M_timestep<<","<<M_num_timesteps<<"\n";
    // ====================================================================================

}//end: update()


template<typename T>
void WimDiscr<T>::updateWaveMedium()
{
    //updates attenuation coefficients, wavelengths and phase/group velocities

    // =============================================================================================
    //std::cout<<"attenuation loop starts (big loop)\n";
    M_max_cg    = -1;
    T_val_vec_ptrs ag_ptrs = {};   //input to interp
    T_val_vec_ptrs agnod_ptrs = {};//output from interp

    std::vector<double> atten_params = M_ice_params.getAttenParams();
    T_val drag_rp = M_ice_params.dragRp();
    T_val gravity = M_ice_params.gravity();
    for (int fq = 0; fq < nwavefreq; fq++)
    {
        if(M_wim_on_mesh)
        {
            ag_ptrs.push_back(&(M_ag_eff[fq]));//input to interp;
            agnod_ptrs.push_back(&(M_agnod_eff[fq]));//output of interp;
        }

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
        {
            double outputs[8];

            if (M_ice[IceType::wim].M_mask[i] >.5)
            {
                T_val om = 2*PI*M_freq_vec[fq];
                T_val guess = std::pow(om,2.)/gravity;

                //if (i==M_itest)
                //{
                //    std::cout<<"fq,om,freq,period = "<<fq<<","<<om<<","<<M_freq_vec[fq]<<","<<1./M_freq_vec[fq]<<"\n";
                //    std::cout<<"guess (open water) = "<<guess<<"\n";
                //}
                if (fq > 0)
                {
                    //improve guess by using output of last call;
                    guess = 2*PI/M_wlng_ice[fq-1][i];
                }

                RTparam_outer(outputs,M_ice[IceType::wim].M_thick[i],double(om),double(drag_rp),double(guess),&atten_params[0]);
                T_val kice    = outputs[1];
                T_val kwtr    = outputs[2];
                T_val int_adm = outputs[3];
                T_val modT    = outputs[5];
                T_val argR    = outputs[6];
                T_val argT    = outputs[7];

                if(M_atten)
                {
                    M_damping[fq][i]      = outputs[0];
                    M_atten_nond[fq][i]   = outputs[4];
                }

                double tmp = 1.;
                for (int io = 0; io<8; io++)
                    tmp *= outputs[io];

                if (std::isnan(tmp))
                {
                    std::cout<<"found NaN in some outputs of RTparam_outer at i,fq="<<i<<","<<fq<<"\n";
                    std::cout<<"\ninputs to RTparam_outer:\n";
                    std::cout<<"h = "<<M_ice[IceType::wim].M_thick[i]<<"\n";
                    std::cout<<"om = "<<om<<"\n";
                    std::cout<<"drag_rp = "<<M_ice_params.dragRp()<<"\n";
                    std::cout<<"guess = "<<guess<<"\n";
                    //
                    std::cout<<"\noutputs from RTparam_outer:\n";
                    std::cout<<"M_damping = "<<M_damping[fq][i]<<"\n";
                    std::cout<<"kice = "<<kice<<"\n";
                    std::cout<<"kwtr = "<<kwtr<<"\n";
                    std::cout<<"int_adm = "<<int_adm<<"\n";
                    std::cout<<"M_atten_nond = "<<M_atten_nond[fq][i]<<"\n";
                    std::cout<<"modT = "<<modT<<"\n";
                    std::cout<<"argR = "<<argR<<"\n";
                    std::cout<<"argT = "<<argT<<"\n";
                    throw std::runtime_error("some outputs of RTparam_outer have NaN");
                }

                //convert amplitude in water to amplitude in ice
                M_disp_ratio[fq][i] = (kice*modT)/kwtr;

                //wavelength to use in ice
                if (1)
                   //use ice wavelength TODO make an option?
                   M_wlng_ice[fq][i] = 2*PI/kice;
                else
                   //use water wavelength instead of ice wavelength
                   M_wlng_ice[fq][i] = M_wlng_wtr[fq];

                //group and phase velocities to use in ice
                if (!M_useicevel)
                {
                   //water group and phase velocities
                   //(ice ones not implemented)
                   M_ag_eff[fq][i] = M_ag_wtr[fq];
                   M_ap_eff[fq][i] = M_ap_wtr[fq];
                }

            }//ice
            else
            {
                M_ag_eff  [fq][i]   = M_ag_wtr  [fq];
                M_ap_eff  [fq][i]   = M_ap_wtr  [fq];
                M_wlng_ice[fq][i]   = M_wlng_wtr[fq];
                M_disp_ratio[fq][i] = 1.;
            }//water

            M_max_cg    = std::max(M_ag_eff[fq][i],M_max_cg);
        }//end i loop
    }//end freq loop
    // =============================================================================================

    if(M_wim_on_mesh)
        //get group velocity on nodes of mesh
        M_mesh.elementsToNodes(agnod_ptrs,ag_ptrs);

}//end: update()


template<typename T>
void WimDiscr<T>::idealWaveFields(T_val const xfac)
{
    M_initialised_waves = true;

    //waves initialised for x<x_edge
    auto xvec = this->getX();
    T_val xmin = *std::min_element(xvec.begin(),xvec.end());
    T_val xmax = *std::max_element(xvec.begin(),xvec.end());
    T_val x_edge = 0.5*(xmin+xmax)-xfac*(0.5*(xmax-xmin));
    T_val_vec wave_mask(M_num_elements,0.); 
    int i_wave = -1;

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {
        if ((xvec[i] < x_edge) && (M_land_mask[i]<.5))
        {
            i_wave  = i;// just need index of one element inside wave mask for testing
            wave_mask[i] = 1.;
            M_Hs [i] = vm["wim.hsinc"].template as<double>();
            M_Tp [i] = vm["wim.tpinc"].template as<double>();
            M_mwd[i] = vm["wim.mwdinc"].template as<double>();
            //std::cout<<M_Hs[i]<<" "<<M_Tp[i]<<" "<<M_mwd[i];
        }
    }

    this->setIncWaveSpec(wave_mask);

    if(M_steady)
    {
        if(i_wave<0)
            throw std::runtime_error("idealWaveFields(): i_wave<0\n");

        for(int fq=0;fq<nwavefreq;fq++)
            for(int nth=0;nth<nwavedirn;nth++)
                // can use any element from inside wave mask since ideal fields are constant in space
                M_open_boundary_vals[fq][nth] = M_sdf_dir[fq][nth][i_wave];
    }

}//idealWaveFields


template<typename T>
void WimDiscr<T>::setWaveFields(T_val_vec const& swh_in,
                                  T_val_vec const& mwp_in,
                                  T_val_vec const& mwd_in)
{
    M_initialised_waves = true;

    bool checkincwaves = vm["wimdiag.checkincwaves"].template as<bool>();
    if (checkincwaves)
    {
        //these arrays are only needed for diagnostics
        M_swh_in = swh_in;
        M_mwp_in = mwp_in;
        M_mwd_in = mwd_in;
    }

#if 0
    // test routine against ideal version
    // - also need to uncomment printout loop after Hs, Tp, mwd & wave_mask are set
    std::cout<<"assign wave fields 2\n";
    auto Hs_old=M_Hs;
    auto Tp_old=M_Tp;
    auto mwd_old=M_mwd;
    auto WM_old=wave_mask;
    this->idealWaveFields(.8);//.8 determines ice edge location
    auto Hs_ideal = M_Hs;
    auto Tp_ideal = M_Tp;
    auto mwd_ideal = M_mwd;
    auto WM_ideal=wave_mask;
    M_Hs  = Hs_old;
    M_Tp  = Tp_old;
    M_mwd = mwd_old;
    wave_mask = WM_old;
#endif


    double Hs_min=100.;
    double Hs_max=0.;
    double Hs_min_i=100.;
    double Hs_max_i=0.;
    double Hs_min_ice=100.;
    double Hs_max_ice=0.;
    T_val_vec wave_mask(M_num_elements,0.);

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {

        if ((M_ice[IceType::wim].M_mask[i]<.5)      //not ice
            &&(M_land_mask[i]<.5)                   //not land
            &&(swh_in[i]>1.e-3)&&(mwp_in[i]>1.e-8)  //some waves
            &&(mwp_in[i]<1.5*Tmax))                 //wave period not too big
        {
           wave_mask[i] = 1.;
           M_Hs [i] = swh_in[i];
           M_Tp [i] = mwp_in[i];
           M_mwd[i] = mwd_in[i];
           //std::cout<<M_Hs[i]<<" "<<M_Tp[i]<<" "<<M_mwd[i];
        }
        
        if (M_ice[IceType::wim].M_mask[i]>0.5)
        {
            Hs_min_ice    = std::min(Hs_min_ice,M_Hs[i]);
            Hs_max_ice    = std::max(Hs_max_ice,M_Hs[i]);
        }
        Hs_min = std::min(Hs_min,M_Hs[i]);
        Hs_max = std::max(Hs_max,M_Hs[i]);
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
#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < num_p_wim; i++)
    {
         std::cout<<"M_Hs["<<i<<"]="<<M_Hs[i]<<","<<Hs_ideal[i]<<","<<swh_in[i]<<"\n";
         std::cout<<"M_Tp["<<i<<"]="<<M_Tp[i]<<","<<Tp_ideal[i]<<","<<mwp_in[i]<<"\n";
         std::cout<<"M_mwd["<<i<<"]="<<M_mwd[i]<<","<<mwd_ideal[i]<<","<<mwd_in[i]<<"\n";
         std::cout<<"wave_mask["<<i<<"]="<<wave_mask[i]<<","<<WM_ideal[i]<<"\n\n";
    }
    std::abort();
#endif

    this->setIncWaveSpec(wave_mask);
}//setWaveFields()


template<typename T>
void WimDiscr<T>::setIncWaveSpec(T_val_vec const& wave_mask)
{

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {
        if (wave_mask[i] == 1.)
        {
            T_val_vec Sfreq(nwavefreq);
            T_val_vec theta_fac(nwavedirn);

            // ============================================================
            // frequency spectrum
            if (nwavefreq == 1)
                Sfreq[0] = std::pow(M_Hs[i]/4.0,2.);
            else
            {
                for (int fq = 0; fq < nwavefreq; fq++)
                {
                    T_val om   = 2*PI*M_freq_vec[fq];
                    T_val t_m  = 2*PI/om;
                    T_val om_m = 2*PI/M_Tp[i];
                    T_val f1   = (5.0/16.0)*std::pow(M_Hs[i],2.)*std::pow(om_m,4.);
                    T_val f2   = 1.0/std::pow(om,5.);
                    T_val f3   = std::exp(-1.25*std::pow(t_m/M_Tp[i],4.));
                    Sfreq[fq]       = f1*f2*f3;

#if 0
                    if (i==M_itest)
                    {
                        std::cout<<"i,fq="<<i<<","<<fq<<"\n";
                        std::cout<<"om,f1,f2,f3="<<om<<","<<f1<<","<<f2<<","<<f3<<"\n";
                        std::cout<<"M_Tp,M_Hs,Sfreq="
                            <<Tp[i]<<","<<M_Hs[i]<<","<<Sfreq[fq]<<"\n";
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
                T_val dtheta = std::abs(M_wavedir[1]-M_wavedir[0]);

                //if (mwd[i]!=0.)
                //   std::cout<<"dir-frac ("<<i<<")"<<std::endl;
                for (int nth = 0; nth < nwavedirn; nth++)
                {
#if 0
                    //less accurate way of calculating spreading
                    //(sample cos^2 at mid-point of interval)
                    T_val chi = PI*(M_wavedir[nth]-M_mwd[i])/180.0;
                    if (std::cos(chi) > 0.)
                        theta_fac[nth] = 2.0*std::pow(std::cos(chi),2.)/PI;
                    else
                        theta_fac[nth] = 0.;
#else
                    //more accurate way of calculating spreading
                    //(integrate cos^2 over interval)
                    theta_fac[nth] = 180./(PI*dtheta)*thetaDirFrac(
                            M_wavedir[nth]-dtheta/2., dtheta, M_mwd[i] );
#endif

                    //if (M_Hs[i*M_grid.M_num_py+j]!=0.)
                    //   std::cout<<M_wavedir[nth]<<" "<<M_mwd[i]<<" "
                    //            <<theta_fac[nth]<<std::endl;
                }
            }//multiple directions
            // ============================================================


            // ============================================================
            // combine freq and dir
            for (int fq = 0; fq < nwavefreq; fq++)
                for (int nth = 0; nth < nwavedirn; nth++)
                {
                    // set M_sdf_dir to inc waves each time new waves are input
                    // NB but only inside the wave mask
                    M_sdf_dir[fq][nth][i] = Sfreq[fq]*theta_fac[nth];

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
    }//end i loop

#if 0
    if(M_itest>=0)
    {
        std::cout<<"i="<<M_itest<<"\n";
        for (int fq = 0; fq < nwavefreq; fq++)
            for (int nth = 0; nth < nwavedirn; nth++)
            {
                std::cout<<"fq,nth="<<fq<<","<<nth<<"\n";
                std::cout<<"M_sdf_dir (setIncWaveSpec) ="<<M_sdf_dir[M_itest][nth][fq]<<"\n";
            }
    }
#endif

    if(M_steady&&(M_cpt==0)&&(!M_wim_on_mesh))
    {
        M_steady_mask = wave_mask;
        M_sdf_dir_inc = M_sdf_dir;
    }

}//setIncWaveSpec


template<typename T>
void WimDiscr<T>::idealIceFields(T_val const xfac)
{

    M_initialised_ice    = true;

    auto xvec = this->getX();
    T_val xmin = *std::min_element(xvec.begin(),xvec.end());
    T_val xmax = *std::max_element(xvec.begin(),xvec.end());

    //ice initialised for x>=x_edge
    T_val x_edge = 0.5*(xmin+xmax)-xfac*(0.5*(xmax-xmin));

    T_val_vec conc  (M_num_elements,0.);
    T_val_vec vol   (M_num_elements,0.);
    T_val_vec nfloes(M_num_elements,0.);
    T_val unifc = vm["wim.unifc"].template as<double>(); /* 0.7 */
    T_val unifh = vm["wim.unifh"].template as<double>(); /* 2.0 */
#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {
        if ((xvec[i] >= x_edge) && (M_land_mask[i]<.5))
        {
            conc[i]      = unifc;
            vol[i]       = unifc*unifh;
            nfloes[i]    = unifc/std::pow(M_dfloe_pack_init,2);
            //std::cout<<M_Hs[i]<<" "<<M_Tp[i]<<" "<<M_mwd[i];
        }
    }
    M_ice[IceType::wim].setFields(conc,vol,nfloes);

}//idealIceFields()


template<typename T>
void WimDiscr<T>::timeStep()
{
    std::fill( M_tau_x.begin()         , M_tau_x.end(), 0. );
    std::fill( M_tau_y.begin()         , M_tau_y.end(), 0. );
    std::fill( M_mwd_x.begin()         , M_mwd_x.end(), 0. );
    std::fill( M_mwd_y.begin()         , M_mwd_y.end(), 0. );
    std::fill( M_stokes_drift_x.begin(), M_stokes_drift_x.end(), 0. );
    std::fill( M_stokes_drift_y.begin(), M_stokes_drift_y.end(), 0. );

    std::fill( Mtmp_mom0      .begin(), Mtmp_mom0      .end(), 0. );
    std::fill( Mtmp_mom2      .begin(), Mtmp_mom2      .end(), 0. );
    std::fill( Mtmp_var_strain.begin(), Mtmp_var_strain.end(), 0. );
    std::fill( Mtmp_mom0w     .begin(), Mtmp_mom0w     .end(), 0. );
    std::fill( Mtmp_mom2w     .begin(), Mtmp_mom2w     .end(), 0. );

    T_val wlng_crest, Dc;
    T_val F, kicel, om, ommin, ommax, om1, lam1, lam2, tmp;
    int jcrest;
    bool break_criterion,test_ij;

    std::string timestpstr = ptime(M_init_time_str, M_current_time);

    // dump local diagnostic file
    // - directory to put it
    std::string outdir = vm["wimdiag.outparentdir"].template as<std::string>();
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
                int wim_itest = vm["wimdiag.itest"].template as<int>();
                int wim_jtest = vm["wimdiag.jtest"].template as<int>();
                diagID << std::setw(16) << std::left
                    << wim_itest << " # itest\n";
                diagID << std::setw(16) << std::left
                    << wim_jtest << " # jtest\n";
            }
            diagID << std::setw(16) << std::left
                << M_ice[IceType::wim].M_mask[M_itest] << " # ICE_MASK\n";
        }
        else
        {
            std::cout << "Cannot open " << diagfile  << "\n";
            std::cerr << "error: open file " << diagfile << " for output failed!" <<"\n";
            std::abort();
        }
        diagID.close();
    }//end M_dump_diag


    if (M_steady&&(!M_wim_on_mesh))//if wim on mesh, use M_open_boundary_vals instead
    {
        for (int fq = 0; fq < nwavefreq; fq++)
            for (int nth = 0; nth < nwavedirn; nth++)
            {
                T_val adv_dir = (-PI/180.)*(M_wavedir[nth]+90.);

                if (std::cos(adv_dir) >= 0.)
#pragma omp parallel for num_threads(M_max_threads) collapse(1)
                    for (int i = 0; i < M_num_elements; i++)
                        if (M_steady_mask[i] > 0.5)
                            M_sdf_dir[fq][nth][i] = M_sdf_dir_inc[fq][nth][i];
            }
    }//M_steady


    //calc mean floe size outside of frequency loop;
    //std::cout<<"calculating <D>\n";
    M_ice[IceType::wim].getDave(M_dave,1.);
    if ( M_dump_diag && (M_itest>-1) )
    {
        int i = M_itest;
        std::fstream diagID(diagfile, std::ios::out | std::ios::app);
        diagID << "\n# Ice info: pre-breaking\n";
        diagID << std::setw(16) << std::left
            << M_ice[IceType::wim].M_mask[i] << " # ice mask\n";
        diagID << std::setw(16) << std::left
            << M_ice[IceType::wim].M_conc[i] << " # conc\n";
        diagID << std::setw(16) << std::left
            << M_ice[IceType::wim].M_thick[i] << " # h, m\n";
        diagID << std::setw(16) << std::left
            << M_dave[i] << " # D_av, m\n";
        diagID << std::setw(16) << std::left
            << M_ice[IceType::wim].M_dfloe[i] << " # D_max, m\n";

        if (M_atten && (M_ice[IceType::wim].M_mask[i]>.5))
            //no attenuation outside ice
            diagID << "\n# period, s | Mtmp_atten_dim, m^{-1}| damp_dim, m^{-1}\n";

        diagID.close();
    }

    for (int fq = 0; fq < nwavefreq; fq++)
    {
        //std::cout<<"calculating dimensional M_atten\n";
        if(M_atten)
        {
            std::fill( Mtmp_atten_dim.begin(), Mtmp_atten_dim.end(), 0. );
            std::fill( Mtmp_damp_dim .begin(), Mtmp_damp_dim .end(), 0. );
        }

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
        {
            if ((M_ice[IceType::wim].M_mask[i] >.5) && (M_atten))
            {
                // floes per unit length
                T_val c1d = M_ice[IceType::wim].M_conc[i]/M_dave[i];

                // scattering
                Mtmp_atten_dim[i] = M_atten_nond[fq][i]*c1d;

                // damping
                Mtmp_damp_dim[i] = 2*M_damping[fq][i]*M_ice[IceType::wim].M_conc[i];

                if ( M_dump_diag && (i==M_itest) )
                {
                   std::fstream diagID(diagfile, std::ios::out | std::ios::app);
                   diagID << 1/M_freq_vec[fq] << "   "
                          << Mtmp_atten_dim[i] << "   "
                          << Mtmp_damp_dim[i] << "\n";
                   diagID.close();
                }
            }//end of ice check
        }//end of spatial loop i


        //advect all directions
        if (!M_wim_on_mesh)
            this->advectDirections(M_sdf_dir[fq],M_ag_eff[fq]);
        else
            this->advectDirectionsMesh(M_sdf_dir[fq],M_agnod_eff[fq],M_open_boundary_vals[fq]);

        //do attenuation &/or scattering, and integrate over directions 
        //std::cout<<"attenuating\n";
        if(!M_atten)
            this->intDirns(M_sdf_dir[fq], Mtmp_sdf_freq,
                    Mtmp_stokes_drift_x_om, Mtmp_stokes_drift_y_om);
        else if (M_scatmod == "dissipated")
            this->attenSimple(M_sdf_dir[fq], Mtmp_sdf_freq, Mtmp_taux_om, Mtmp_tauy_om,
                    Mtmp_stokes_drift_x_om, Mtmp_stokes_drift_y_om,M_ag_eff[fq]);
        else if (M_scatmod == "isotropic")
            this->attenIsotropic(M_sdf_dir[fq], Mtmp_sdf_freq, Mtmp_taux_om, Mtmp_tauy_om,
                    Mtmp_stokes_drift_x_om, Mtmp_stokes_drift_y_om,M_ag_eff[fq]);


        // integrate stress and stokes drift densities over frequency
        // integrals for breaking
        T_val tau_fac = M_ice_params.rhowtr()*M_ice_params.gravity();
#pragma omp parallel for num_threads(M_max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
        {
            T_val om    = 2*PI*M_freq_vec[fq];    // radial freq
            T_val cp    = M_ap_eff[fq][i];        // phase velocity
            T_val kicel = 2*PI/M_wlng_ice[fq][i]; // ice wave number (just water wavelength if no ice)
            T_val F     = M_disp_ratio   [fq][i]; // convert from water amp's to ice amp's
            T_val tmp1  = 0.;
            T_val F2    = 1.;
            if (M_ref_Hs_ice)
                F2  = std::pow(F,2);//for outputs only

            // ================================================================================
            // integrate stress, MWD and stokes drift densities over frequency
            if(M_atten)
            {
                tmp1 = tau_fac*Mtmp_taux_om[i]/cp;
                M_tau_x[i] += M_quadrature_wt_freq[fq]*tmp1;

                tmp1 = tau_fac*Mtmp_tauy_om[i]/cp;
                M_tau_y[i] += M_quadrature_wt_freq[fq]*tmp1;
            }

            //integrals for MWD
            M_mwd_x[i] += M_quadrature_wt_freq[fq]*F2*Mtmp_stokes_drift_x_om[i];
            M_mwd_y[i] += M_quadrature_wt_freq[fq]*F2*Mtmp_stokes_drift_y_om[i];

            //2*\omega*k*\int_0^\pi S(\omega,\theta)\cos(\theta)d\theta
            tmp1 = 2*om*kicel*F2*Mtmp_stokes_drift_x_om[i];
            M_stokes_drift_x[i] += M_quadrature_wt_freq[fq]*tmp1;

            //2*\omega*k*\int_0^\pi S(\omega,\theta)\sin(\theta)d\theta
            tmp1 = 2*om*kicel*F2*Mtmp_stokes_drift_y_om[i];
            M_stokes_drift_y[i] += M_quadrature_wt_freq[fq]*tmp1;
            // ================================================================================

            // ================================================================================
            // integrals for breaking
            
            // -----------------------------------------------------------------
            // 0-th spectral moments
            // - take abs as small errors can make Mtmp_sdf_freq negative
            tmp1    = M_quadrature_wt_freq[fq]*Mtmp_sdf_freq[i];

            // variance of displacement (water)
            Mtmp_mom0w[i] += std::abs(tmp1);

            // variance of displacement (ice)
            Mtmp_mom0[i] += std::abs(tmp1*std::pow(F,2.));
            // -----------------------------------------------------------------

            // -----------------------------------------------------------------
            // 2-nd spectral moments
            tmp1    = M_quadrature_wt_freq[fq]*std::pow(om,2.)*Mtmp_sdf_freq[i];

            // variance of speed (water)
            Mtmp_mom2w[i] += std::abs(tmp1);

            // variance of speed (ice)
            Mtmp_mom2[i] += std::abs(tmp1*std::pow(F,2.));
            // -----------------------------------------------------------------

            // -----------------------------------------------------------------
            // variance of strain
            if (M_ice[IceType::wim].M_mask[i] == 1.)
            {
                // strain conversion factor
                // = k^2*h/2*F
                tmp1 = F*std::pow(kicel,2.)*M_ice[IceType::wim].M_thick[i]/2.0;

                // strain density
                tmp1           = M_quadrature_wt_freq[fq]*Mtmp_sdf_freq[i]*std::pow(tmp1,2.);
                Mtmp_var_strain[i] += std::abs(tmp1);
            }
            // -----------------------------------------------------------------

            // ================================================================================


            // ================================================================================
            //do some checks
            if (std::isnan(Mtmp_sdf_freq[i]))
            {
                std::cout<<"fq = "<<fq<<"\n";
                std::cout<<"found NaN in Mtmp_sdf_freq at i = "<<i<<"\n";
                std::cout<<"F,kicel,om,M_quadrature_wt_freq="<<F<<","<<kicel<<","<<om<<","<<M_quadrature_wt_freq[fq]<<"\n";
                throw std::runtime_error("Mtmp_sdf_freq has NaN (after advection & attenuation)");
            }
            if (std::isnan(F))
            {
                std::cout<<"fq = "<<fq<<"\n";
                std::cout<<"found NaN in M_disp_ratio at i = "<<i<<"\n";
                throw std::runtime_error("M_disp_ratio has NaN");
            }

            if(M_atten)
            {
                if (std::isnan(Mtmp_taux_om[i]))
                {
                    std::cout<<"fq = "<<fq<<"\n";
                    std::cout<<"found NaN in Mtmp_taux_om at i = "<<i<<"\n";
                    throw std::runtime_error("taux_om has NaN (after advection & attenuation)");
                }
                if (std::isnan(Mtmp_tauy_om[i]))
                {
                    std::cout<<"fq = "<<fq<<"\n";
                    std::cout<<"found NaN in Mtmp_tauy_om at i = "<<i<<"\n";
                    throw std::runtime_error("Mtmp_tauy_om has NaN (after advection & attenuation)");
                }
            }
            // ================================================================================
        }//end i loop
    }//end freq loop

    // T_val _min = *std::min_element(mom0w.begin(),mom0w.end());
    // T_val _max = *std::max_element(mom0w.begin(),mom0w.end());
    // std::cout<<"Min f= " << _min <<"\n";
    // std::cout<<"Max f= " << _max <<"\n";


    // for (int i = 0; i < M_num_elements; i++)
    //     std::cout << "VRT[" << i < "]= " << var_strain[i] <<"\n";

    //update integrated variables
    std::fill( M_Tp .begin(), M_Tp .end(), 0. );
    std::fill( M_mwd.begin(), M_mwd.end(), 0. );

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {
        if (M_ref_Hs_ice)
        {
            M_Hs[i] = 4*std::sqrt(Mtmp_mom0[i]);
            if (Mtmp_mom2[i] > 0.)
            {
                M_Tp[i]   = 2*PI*std::sqrt(Mtmp_mom0[i]/Mtmp_mom2[i]);

                //mwd: waves-from dirn and degrees
                M_mwd[i]  = -90.-(180./PI)*std::atan2(M_mwd_y[i],M_mwd_x[i]);
            }
        }
        else
        {
            M_Hs[i] = 4*std::sqrt(Mtmp_mom0w[i]);
            if (Mtmp_mom2w[i] > 0.)
            {
                M_Tp[i]   = 2*PI*std::sqrt(Mtmp_mom0w[i]/Mtmp_mom2w[i]);

                //mwd: waves-from dirn and degrees
                M_mwd[i]  = -90.-(180./PI)*std::atan2(M_mwd_y[i],M_mwd_x[i]);
            }
        }
    }


    if ( M_dump_diag )
    {
       int i =  M_itest;
       std::fstream diagID(diagfile, std::ios::out | std::ios::app);

       diagID << std::setw(16) << std::left
          << Mtmp_mom0w[i] << " # mom0w, m^2\n";
       diagID << std::setw(16) << std::left
          << Mtmp_mom2w[i] <<" # mom2w, m^2/s^2\n";
       diagID << std::setw(16) << std::left
          << Mtmp_mom0[i] <<" # mom0, m^2\n";
       diagID << std::setw(16) << std::left
          << Mtmp_mom2[i] <<" # mom2, m^2/s^2\n";
       diagID << std::setw(16) << std::left
          << M_Hs[i] <<" # Hs, m\n";
       diagID << std::setw(16) << std::left
          << M_Tp[i] <<" # Tp, s\n";
       diagID << std::setw(16) << std::left
          << M_mwd[i] <<" # mwd, deg\n";
       diagID << std::setw(16) << std::left
          << M_tau_x[i] <<" # tau_x, Pa\n";
       diagID << std::setw(16) << std::left
          << M_tau_y[i] <<" # tau_y, Pa\n";
       diagID.close();
    }


    // finally do floe breaking
    if(M_breaking)
        M_ice[IceType::wim].doBreaking(Mtmp_mom0,Mtmp_mom2,Mtmp_var_strain);

    if (M_break_on_mesh&&M_breaking) // breaking on nextsim mesh
    {

        // =================================================================
        //std::cout<<"break_on_mesh: before interp grid to mesh\n";
        // do interpolation
        // - set input data
        T_val_vec_ptrs input_data = {&Mtmp_mom0,&Mtmp_mom2,&Mtmp_var_strain};

        // - set output data
        // - these are automatically resized in M_grid.interpToPoints()
        T_val_vec mom0_mesh;
        T_val_vec mom2_mesh;
        T_val_vec var_strain_mesh;
        T_val_vec_ptrs output_data = {&mom0_mesh,&mom2_mesh,&var_strain_mesh};

        // - call routine
        M_grid.interpToPoints(output_data,input_data,
                M_mesh.M_elements_x, M_mesh.M_elements_y);
        //std::cout<<"break_on_mesh: after interp grid to mesh\n";
        // =================================================================


        // =================================================================
        //do breaking
        M_ice[IceType::sim].doBreaking(mom0_mesh,mom2_mesh,var_strain_mesh);
        //std::cout<<"break_on_mesh: after breaking\n";
    }//break_on_mesh


    this->printRange("Hs",M_Hs);
    this->printRange("tau_x",M_tau_x,10);
}//timeStep


template<typename T>
void WimDiscr<T>::setMesh(T_gmsh const &movedmesh)
{
    M_time_mesh_set = M_update_time;//used in check when ice fields are set on mesh
    M_mesh_old      = M_mesh;

    //update M_mesh with moved mesh
    M_mesh  = T_mesh(movedmesh);
}


template<typename T>
void WimDiscr<T>::setMesh(T_gmsh const &mesh_in,T_val_vec const &um_in)
{
    auto movedmesh = mesh_in;
    movedmesh.move(um_in,1.);
    this->setMesh(movedmesh);
}//setMesh


template<typename T>
void WimDiscr<T>::setMesh2(T_gmsh const &mesh_in,
        T_val_vec const &um_in,BamgMesh* bamgmesh,int const& flag_fix,bool const& regridding)
{
    //interface for M_wim_on_mesh
    auto movedmesh = mesh_in;
    movedmesh.move(um_in,1.);
    this->setMesh2(movedmesh,bamgmesh,flag_fix,regridding);
}


template<typename T>
void WimDiscr<T>::setMesh2(T_gmsh const &movedmesh,BamgMesh* bamgmesh,int const& flag_fix,bool const& regridding)
{
    //interface for M_wim_on_mesh

    if(regridding)
        //need to call assignSpatial() to resize vectors
        M_assigned = false;
    else
        M_mesh_old = M_mesh;

    M_time_mesh_set = M_update_time;//used in check when ice fields are set on mesh
    M_mesh          = T_mesh(movedmesh,bamgmesh,flag_fix);

    // get relative displacement of nodes since last call
    // - M_UM may already be nonzero if regridding has happened
    // - it is reset to zero at end of wim.run() and at initialisation
    // - used to correct group velocity when waves are advected
    int Nn = M_mesh.M_num_nodes;
    int sz = M_UM.size(); 
    bool update_mesh_disp = (!regridding);
    if(sz==0)
    {
        update_mesh_disp = false;
        M_UM.assign(2*Nn,0.);
    }
    else if(sz!=2*Nn)
    {
        std::string msg  = (boost::format(  "M_UM is wrong size (%1%) - should be %2%\n" )
                % sz % (2*Nn)).str();
        throw std::runtime_error(msg);
    }

    if(update_mesh_disp)
        //NB if regridding don't need to add any correction to M_UM
        for (int i=0;i<Nn;i++)
        {
            //nextsim_mesh_old is either from last WIM call or last regrid
            M_UM[i]    += M_mesh.M_nodes_x[i]-M_mesh_old.M_nodes_x[i];
            M_UM[i+Nn] += M_mesh.M_nodes_y[i]-M_mesh_old.M_nodes_y[i];
        }

    // ================================================================================

    M_num_elements  = M_mesh.M_num_elements;
    std::cout<<"on mesh, M_num_elements = "<<M_num_elements<<"\n";

    if(!M_assigned)
    {
        //need to set sizes each time mesh changes: init,regrid
        std::cout<<"calling assignSpatial() inside setMesh2()\n";
        this->assignSpatial();
        M_land_mask.assign(M_num_elements,0.);
    }

    M_length_cfl = M_mesh.lengthCfl();
}//setMesh2


template<typename T>
typename WimDiscr<T>::T_val_vec
WimDiscr<T>::getSurfaceFactor(T_gmsh const &movedmesh)
{
    // wave spectrum needs to be updated if mesh changes due to divergence of mesh velocity
    // ie element surface area changes need to be taken into account;
    // call this before setMesh2() at regrid time or before call to WIM
    auto nodes_x = movedmesh.coordX();
    auto nodes_y = movedmesh.coordY();
    auto index   = movedmesh.indexTr();

    int Nels = movedmesh.numTriangles();
    T_val_vec surface_fac(Nels,0.);
    for (int i=0;i<Nels;i++)
    {
        T_val_vec xnods(3);
        T_val_vec ynods(3);
        for (int k=0;k<3;k++)
        {
            int ind  = index[3*i+k]-1;//NB bamg index starts at 1
            xnods[k] = nodes_x[ind];
            ynods[k] = nodes_y[ind];
        }

        T_val area = .5*MeshTools::jacobian(
                xnods[0],ynods[0],xnods[1],ynods[1],xnods[2],ynods[2]);
        if (area>0)
            surface_fac[i] = area/M_mesh.M_surface[i];
        else
        {
            std::cout<<"Area of triangle "<<i<<" <0: "<<area<<"\n";
            throw std::runtime_error("getSurfaceFactor: found negative area\n");
        }
    }

    return surface_fac;
}//getSurfaceFactor()


template<typename T>
void WimDiscr<T>::updateWaveSpec(T_gmsh const &movedmesh)
{
    // wave spectrum needs to be updated if mesh changes due to divergence of mesh velocity
    // ie element surface area changes need to be taken into account;
    // call this before setMesh2() at regrid time or before call to WIM
    auto nodes_x = movedmesh.coordX();
    auto nodes_y = movedmesh.coordY();
    auto index   = movedmesh.indexTr();

    std::cout<<"updateWaveSpec: calling getSurfaceFactor";
    auto surface_fac = this->getSurfaceFactor(movedmesh);

    int Nels = movedmesh.numTriangles();
    std::fill(M_Tp.begin() ,M_Tp.end() ,0.);
    std::fill(M_mwd.begin(),M_mwd.end(),0.);
    for (int i=0;i<Nels;i++)
    {
        //integrate wave spectrum here
        T_val mom0 = 0.;
        T_val mom2 = 0.;
        T_val momc = 0.;
        T_val moms = 0.;
        T_val sdfx = 0.;
        T_val sdfy = 0.;

        for(int fq=0;fq<nwavefreq;fq++)
        {
            T_val kice = 2*PI/M_wlng_ice[fq][i];
            T_val om   = 2*PI*M_freq_vec[fq];
            T_val om2  = std::pow(om,2);
            T_val F2   = 1.;
            if(M_ref_Hs_ice)
                F2   = std::pow(M_disp_ratio[fq][i],2);//TODO should stokes drift be a relative thing? maybe should take conc-weighted average?

            for(int nth=0;nth<nwavedirn;nth++)
            {
                T_val adv_dir          = -PI*(90.0+M_wavedir[nth])/180.0;
                M_sdf_dir[fq][nth][i] *= surface_fac[i];
                T_val sdf              = M_sdf_dir[fq][nth][i];

                mom0 += M_quadrature_wt_freq[fq]*M_quadrature_wt_dir[nth]*sdf*F2;
                mom2 += M_quadrature_wt_freq[fq]*M_quadrature_wt_dir[nth]*sdf*F2*om2;
                //
                T_val tmp = M_quadrature_wt_freq[fq]*M_quadrature_wt_dir[nth]*sdf*F2*std::cos(adv_dir);
                momc     += tmp;
                sdfx     += 2*om*kice*tmp;
                //
                tmp   = M_quadrature_wt_freq[fq]*M_quadrature_wt_dir[nth]*sdf*F2*std::sin(adv_dir);
                moms += tmp;
                sdfy += 2*om*kice*tmp;
            }//nth
        }//fq

        M_Hs[i]   = 4*std::sqrt(mom0);
        if(mom2>0.);
        {
            M_Tp[i]   = 2*PI*std::sqrt(mom0/mom2);
            M_mwd[i]  = std::atan2(moms,momc);
        }
        M_stokes_drift_x[i]   = sdfx;
        M_stokes_drift_y[i]   = sdfy;
    }//loop over elements
}//updateWaveSpec


template<typename T>
void WimDiscr<T>::updateWaveSpec(T_gmsh const &mesh_in,T_val_vec const &um_in)
{
    auto movedmesh = mesh_in;
    movedmesh.move(um_in,1);
    this->updateWaveSpec(movedmesh);
}//updateWaveSpec


template<typename T>
typename WimDiscr<T>::T_val_vec
WimDiscr<T>::getRelativeMeshDisplacement(T_gmsh const &mesh_in,T_val_vec const &um_in) const
{
    auto movedmesh = mesh_in;
    movedmesh.move(um_in,1);
    this->getRelativeMeshDisplacement(movedmesh);
}//getRelativeMeshDisplacement


template<typename T>
typename WimDiscr<T>::T_val_vec
WimDiscr<T>::getRelativeMeshDisplacement(T_gmsh const &movedmesh) const
{
    auto nodes_x = movedmesh.coordX();
    auto nodes_y = movedmesh.coordY();
    int Nn = nodes_x.size();

    if (M_mesh.M_mesh_type==T_mesh::E_mesh_type::uninitialised)
        throw std::runtime_error("relativeMeshDisplacement: M_mesh not initialised yet");
    if (M_mesh.M_num_nodes!=Nn)
        throw std::runtime_error("relativeMeshDisplacement: mesh_in and M_mesh have different sizes");

    auto um_out = M_UM;//in case there has been another regrid already

    for (int i=0;i<Nn;i++)
    {
        um_out[i]    += nodes_x[i]-M_mesh.M_nodes_x[i];
        um_out[i+Nn] += nodes_y[i]-M_mesh.M_nodes_y[i];
    }

    return um_out;
}//relativeMeshDisplacement


template<typename T>
void WimDiscr<T>::setIceFields(
                          std::vector<T_val> const& conc,  // conc
                          std::vector<T_val> const& vol,   // ice vol or effective thickness (conc*thickness)
                          std::vector<T_val> const& nfloes,// Nfloes=conc/Dmax^2
                          bool pre_regrid)
{
    if(M_time_mesh_set != M_current_time)
    {
        std::cout<<"M_time_mesh_set,M_current_time = "<<M_time_mesh_set<<","<<M_current_time<<"\n";
        throw std::runtime_error("setIceFields: setting ice without setting mesh first");
    }

    // pre-regrid options:
    if(pre_regrid&&M_wim_on_mesh)
        //do nothing
        return;

    //if here, not (M_wim_on_mesh && pre_regrid)
    M_initialised_ice = true;
    if (pre_regrid)//not M_wim_on_mesh
    {
        //interp from mesh to grid
        M_ice[IceType::sim].setFields(conc,vol,nfloes);
        this->interpIceMeshToGrid();
        return;
    }
    // end of pre-regrid options
    // ================================================================

    // ================================================================
    // post-regrid options:
    if (M_wim_on_mesh)
        // ice fields already where we need them,
        M_ice[IceType::wim].setFields(conc,vol,nfloes);
    else if (M_break_on_mesh)
        // if(M_break_on_mesh), ice fields already where we need them,
        M_ice[IceType::sim].setFields(conc,vol,nfloes);

    // end of post-regrid options
    // ================================================================

}//setIceFields()


template<typename T>
void WimDiscr<T>::interpIceMeshToGrid()
{
    T_val_vec_ptrs input_data  = {&(M_ice[IceType::sim].M_conc),&(M_ice[IceType::sim].M_vol),&(M_ice[IceType::sim].M_nfloes)};
    T_val_vec_ptrs output_data = {&(M_ice[IceType::wim].M_conc),&(M_ice[IceType::wim].M_vol),&(M_ice[IceType::wim].M_nfloes)};
    M_grid.interpFromMesh(M_mesh,output_data,input_data);
    M_ice[IceType::wim].updateFields();
}//iceMeshToGrid


template<typename T>
typename WimDiscr<T>::T_map_vec
WimDiscr<T>::returnFieldsNodes(std::vector<std::string> const & fields,
        T_gmsh const &movedmesh)
{
    auto xnod = movedmesh.coordX();
    auto ynod = movedmesh.coordY();
    return this->returnFieldsNodes(fields,xnod,ynod);
}


template<typename T>
typename WimDiscr<T>::T_map_vec
WimDiscr<T>::returnFieldsNodes(std::vector<std::string> const & fields,
        T_gmsh const &mesh_in,T_val_vec const &um_in)
{
    auto movedmesh  = mesh_in;
    movedmesh.move(um_in,1.);
    return this->returnFieldsNodes(fields,movedmesh);
}


template<typename T>
typename WimDiscr<T>::T_map_vec
WimDiscr<T>::returnFieldsElements(std::vector<std::string> const &fields,
        T_gmsh const &movedmesh)
{
    auto xel  = movedmesh.bcoordX();
    auto yel  = movedmesh.bcoordY();

    T_val_vec surface_fac(xel.size(),1.);
    if(M_wim_on_mesh)
    {
        std::cout<<"returnFieldsElements: calling getSurfaceFactor";
        surface_fac = this->getSurfaceFactor(movedmesh);
    }

    return this->returnFieldsElements(fields,xel,yel,surface_fac);
}


template<typename T>
typename WimDiscr<T>::T_map_vec
WimDiscr<T>::returnFieldsElements(std::vector<std::string> const &fields,
        T_gmsh const &mesh_in,T_val_vec const &um_in)
{
    auto movedmesh  = mesh_in;
    movedmesh.move(um_in,1.);
    return this->returnFieldsElements(fields,movedmesh);
}


template<typename T>
typename WimDiscr<T>::T_map_vec
WimDiscr<T>::returnFieldsNodes(std::vector<std::string> const &fields,
        T_val_vec &xnod, T_val_vec &ynod)
{
    // return fields to nodes of nextsim mesh
    // - usually to export diagnostic fields on nextsim mesh
    T_map_vec output_nodes;
    int Nnod = xnod.size();

    if (fields.size()==0)
        //nothing to do
        return output_nodes;
    else
        //initialise outputs
        for (auto it=fields.begin();it!=fields.end();it++)
        {
            T_val_vec tmp(2*Nnod,0.);
            output_nodes.emplace(*it,tmp);
#if 0
            //TODO add check for M_wim_on_mesh after M_sdf_dir is reset at regrid time
            //check if need to integrate spectrum before the export
            if ( (!M_wavespec_integrated) && (*it=="Stokes_drift") )
                this->intWaveSpec(); TODO define this function
#endif
        }

    // ==========================================================================================
    //nodes - vectors

    //input to and output from interpolation routines
    T_val_vec_ptrs input_nodes, out_nodes;

    //need to make some dummy variables since vector components handled individually
    T_val_vec tx_out,ty_out,sdfx_out,sdfy_out;
    for (auto it = output_nodes.begin(); it != output_nodes.end(); it++)
    {
        //get inputs
        if(it->first=="Stress_waves_ice")
        {
            input_nodes.push_back(&M_tau_x);
            input_nodes.push_back(&M_tau_y);
            out_nodes.push_back(&tx_out);
            out_nodes.push_back(&ty_out);
        }
        else if(it->first=="Stokes_drift")
        {
            input_nodes.push_back(&M_stokes_drift_x);
            input_nodes.push_back(&M_stokes_drift_y);
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
        M_mesh.elementsToNodes(out_nodes,input_nodes);
    else
        //from grid elements to mesh nodes
        M_grid.interpToPoints(out_nodes,input_nodes,xnod,ynod);
    // ==========================================================================================


    // ==========================================================================================
    // assign the outputs
    for (auto it = output_nodes.begin(); it != output_nodes.end(); it++)
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
typename WimDiscr<T>::T_map_vec
WimDiscr<T>::returnFieldsElements(std::vector<std::string> const &fields,
        T_val_vec &xel, T_val_vec &yel, T_val_vec const& surface_fac)
{
    // return fields on elements of M_mesh
    // - usually to export diagnostic fields on nextsim mesh
    T_map_vec output_els;
    int Nels = xel.size();

    if (fields.size()==0)
        //do nothing
        return output_els;
    else
        //initialise outputs
        for (auto it=fields.begin();it!=fields.end();it++)
        {
            T_val_vec tmp(Nels,0.);
            output_els.emplace(*it,tmp);
        }

    // ==========================================================================================
    //elements - scalars

    //input to and output from interpolation routines
    T_val_vec_ptrs input_els, out_els;
    for (auto it = output_els.begin(); it != output_els.end(); it++)
    {
        //get inputs
        if(it->first=="Hs")
        {
            if (M_wim_on_mesh)
                for(int i=0;i<M_num_elements;i++)
                    it->second[i] = std::sqrt(surface_fac[i])*M_Hs[i];//NB SDF scales with surface area, so Hs scales by sqrt(SDF)
            else
            {
                input_els.push_back(&M_Hs);
                out_els.push_back(&(it->second));
            }
        }
        else if(it->first=="Tp")
        {
            if (M_wim_on_mesh)
                it->second = M_Tp;//NB indep of surface area
            else
            {
                input_els.push_back(&M_Tp);
                out_els.push_back(&(it->second));
            }
        }
        else if(it->first=="MWD")
        {
            if (M_wim_on_mesh)
                it->second = M_mwd;//NB indep of surface area
            else
            {
                input_els.push_back(&M_mwd);
                out_els.push_back(&(it->second));
            }
        }
        else
            throw std::runtime_error("returnFieldsElements): unknown variable name - "+it->first+"\n");
    }
    // ==========================================================================================

    //interp to elements if necessary
    if (!M_wim_on_mesh)
        M_grid.interpToPoints(out_els,input_els,xel,yel);//out_els already points to output_els

    return output_els;
}//returnFieldsElements


template<typename T>
void WimDiscr<T>::returnWaveStress(T_val_vec &M_tau,T_gmsh const &mesh_in,T_val_vec const &um_in)
{
    auto movedmesh  = mesh_in;
    movedmesh.move(um_in,1.);
    this->returnWaveStress(M_tau,movedmesh);
}


template<typename T>
void WimDiscr<T>::returnWaveStress(T_val_vec &M_tau,T_gmsh const &movedmesh)
{
    auto xnod = movedmesh.coordX();
    auto ynod = movedmesh.coordY();
    this->returnWaveStress(M_tau,xnod,ynod);
}


template<typename T>
void WimDiscr<T>::returnWaveStress(T_val_vec &M_tau,T_val_vec &xnod,T_val_vec &ynod)
{
    //return wave stress on nodes of nextsim mesh

    int Nnod = xnod.size();

    //nodes - vectors
    T_val_vec tx_out,ty_out;
    T_val_vec_ptrs input_nodes = {&M_tau_x,&M_tau_y};
    T_val_vec_ptrs out_nodes   = {&tx_out,&ty_out};

    //interp to nodes
    if(!M_wim_on_mesh)
        M_grid.interpToPoints(out_nodes,input_nodes,xnod,ynod);
    else
        throw std::runtime_error("returnWaveStress: using wrong interface for M_wim_on_mesh");

    M_tau.resize(2*Nnod,0.);
    for (int i=0;i<Nnod;i++)
    {
        M_tau[i]        = tx_out[i];
        M_tau[i+Nnod]   = ty_out[i];
    }

}//returnWaveStress


template<typename T>
void WimDiscr<T>::returnWaveStress(T_val_vec &M_tau)
{
    //nodes - vectors
    T_val_vec tx_out,ty_out;
    T_val_vec_ptrs input_nodes = {&M_tau_x,&M_tau_y};
    T_val_vec_ptrs out_nodes   = {&tx_out,&ty_out};
    M_mesh.elementsToNodes(out_nodes,input_nodes);
}


template<typename T>
typename WimDiscr<T>::T_val
WimDiscr<T>::dfloeToNfloes(T_val const& dfloe_in,
                           T_val const& conc_in)
{
    T_val nfloes_out   = 0.;

    if ( (dfloe_in>0) &&(conc_in >= M_ice_params.ciceMin()) )
    {
        //conc high enough & dfloe OK
        nfloes_out = conc_in/std::pow(dfloe_in,2.);
    }

    return nfloes_out;
}//dfloesToNfloes


template<typename T>
typename WimDiscr<T>::T_val_vec
WimDiscr<T>::dfloeToNfloes(T_val_vec const& dfloe_in,
                           T_val_vec const& conc_in)
{
    int N = conc_in.size();
    T_val_vec nfloes_out(N);

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i=0;i<N;i++)
        nfloes_out[i] = this->dfloeToNfloes(dfloe_in[i],conc_in[i]);

    return nfloes_out;
}//dfloesToNfloes


template<typename T>
typename WimDiscr<T>::T_val
WimDiscr<T>::nfloesToDfloe(T_val const& nfloes_in,
                           T_val const& conc_in)
{
        T_val dfloe_out  = 0.;
        if ( (nfloes_in>0)
                &&(conc_in >= vm["wim.cicemin"].template as<double>()) )
        {
            //conc high enough & Nfloes OK
            dfloe_out = std::sqrt(conc_in/nfloes_in);
        }

        //dfloe shouldn't get too big
        if ( dfloe_out>=vm["wim.dfloepackthresh"].template as<double>() )
            dfloe_out = vm["wim.dfloepackinit"].template as<double>();

    return dfloe_out;
}//nfloesToDfloe


template<typename T>
typename WimDiscr<T>::T_val_vec
WimDiscr<T>::nfloesToDfloe(T_val_vec const& nfloes_in,
                           T_val_vec const& conc_in)
{
    int N   = conc_in.size();
    T_val_vec dfloe_out(N);

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i=0;i<N;i++)
        dfloe_out[i] = this->nfloesToDfloe(nfloes_in[i],conc_in[i]);

    return dfloe_out;
}//nfloesToDfloe


template<typename T>
void WimDiscr<T>::getFsdMesh(T_val_vec &nfloes_out,T_val_vec &dfloe_out,T_val_vec &broken)
{
    if (M_wim_on_mesh)
    {
        nfloes_out  = M_ice[IceType::wim].M_nfloes;
        dfloe_out   = M_ice[IceType::wim].M_dfloe;
        broken      = M_ice[IceType::wim].M_broken;
    }
    else if (M_break_on_mesh)
    {
        nfloes_out  = M_ice[IceType::sim].M_nfloes;
        dfloe_out   = M_ice[IceType::sim].M_dfloe;
        broken      = M_ice[IceType::sim].M_broken;
    }
    else
        throw std::runtime_error("getFsdMesh: using wrong interface");
}//getFsdMesh


template<typename T>
void WimDiscr<T>::getFsdMesh(T_val_vec &nfloes_out,T_val_vec &dfloe_out,T_val_vec &broken,
        T_val_vec const & conc_tot, T_gmsh const &mesh_in,T_val_vec const &um_in)
{
    if((M_wim_on_mesh)||(M_break_on_mesh))
        throw std::runtime_error("getFsdMesh: using wrong interface");

    // set M_mesh (need to know where to interpolate to)
    // - NB set in FiniteElement::wimPreRegrid(),
    // but this function is called from FiniteElement::wimPostRegrid(),
    // and mesh could have changed due to regridding
    auto movedmesh  = mesh_in;
    movedmesh.move(um_in,1.);
    this->getFsdMesh(nfloes_out,dfloe_out,broken,conc_tot,movedmesh);
}


template<typename T>
void WimDiscr<T>::getFsdMesh(T_val_vec &nfloes_out,T_val_vec &dfloe_out,T_val_vec &broken,
        T_val_vec const & conc_tot, T_gmsh const &movedmesh)
{
    if((M_wim_on_mesh)||(M_break_on_mesh))
        throw std::runtime_error("getFsdMesh: using wrong interface");

    // set M_mesh (need to know where to interpolate to)
    // - NB set in FiniteElement::wimPreRegrid(),
    // but this function is called from FiniteElement::wimPostRegrid(),
    // and mesh could have changed due to regridding
    M_mesh    = T_mesh(movedmesh);

    //do interpolation
    T_val_vec cinterp;//interp conc as well, to correct for interpolation error
    T_val_vec_ptrs input_data = {&(M_ice[IceType::wim].M_nfloes),&(M_ice[IceType::wim].M_conc),&(M_ice[IceType::wim].M_broken)};
    T_val_vec_ptrs output_data = {&nfloes_out,&cinterp,&broken};
    M_grid.interpToPoints(output_data,input_data,M_mesh.M_elements_x,M_mesh.M_elements_y);

    int Nel = M_mesh.M_num_elements;
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
                // ctot    may have been interpolated: old mesh->new mesh (if no regridding, old = new & M_ice[IceType::sim].conc = ctot)
                // cinterp and ctot should be the same, apart from errors in interpolation,
                // so hopefully the errors in nfloes can be estimated from this process
                T_val cfac = conc_tot[i]/cinterp[i];
                nfloes_out[i] = cfac*nfloes_out[i];
            }

            //else keep nfloes_out the same
            dfloe_out[i] = nfloesToDfloe(nfloes_out[i],conc_tot[i]);
            broken[i]    = std::round(broken[i]);
        }
        else
        {
            nfloes_out[i] = 0.;
            broken[i]     = 0.;
        }
    }
}//getFsdMesh


template<typename T>
void WimDiscr<T>::run()
{

    //====================================================
    // set ice & wave conditions
    // - incident wave spectrum set in here now
    // (or in setWaveFields)
    if (!M_initialised_ice)
    {
        std::cout<<"WIM: Calling idealIceFields()";
        this->idealIceFields(0.7);
    }

    if (!M_initialised_waves)
    {
        std::cout<<"WIM: Calling idealWaveFields()";
        this->idealWaveFields(0.8);
    }
    // ===================================================


    // set attenuation coefficients and wave speeds/lengths
    // - also set M_timestep
    this->update();

    std::cout << "-----------------------Simulation started at "<< Wim::current_time_local() <<"\n";

    //M_init_time_str is human readable time eg "2015-01-01 00:00:00"
    //init_time is "formal" time format eg "20150101T000000Z"
    std::string init_time = ptime(M_init_time_str);

    int lcpt  = 0;//local counter
    M_current_time  = this->getModelTime(lcpt);//model time of current call to wim (relative to init_time)
    std::string call_time = ptime(M_init_time_str,M_current_time);
    std::cout<<"---------------INITIAL TIME: "<< init_time <<"\n";
    std::cout<<"---------------CALLING TIME: "<< call_time <<"\n";


    std::cout<<"Running starts\n";
    chrono.restart();

    std::cout<<"M_duration = "<< M_duration <<"\n";
    std::cout<<"M_max_cg = "<< M_max_cg <<"\n";
    std::cout<<"M_length_cfl = "<< M_length_cfl <<"\n";
    std::cout<<"M_timestep = "<< M_timestep <<"\n";
    std::cout<<"M_num_timesteps = "<< M_num_timesteps <<"\n";

    if (vm["wimdiag.checkinit"].template as<bool>())
        this->exportResults("init");

    if (vm["wimdiag.checkincwaves"].template as<bool>()
        &&(M_swh_in.size()>0))
        this->exportResults("incwaves");

#if 1
    if (M_swh_in.size()>0)
    {
        //test inc waves
        this->printRange("swh (in)",M_swh_in);
        this->printRange("mwp (in)",M_mwp_in);
        this->printRange("mwd (in)",M_mwd_in);
    }
#endif

    // check sizes of wave arrays
    if (M_Hs.size()!=M_Tp.size()||M_Hs.size()!=M_num_elements)
    {
        std::cout<<"M_cpt = "<<M_cpt<<"\n";
        std::cout<<"size Hs = "<<M_Hs.size()<<"\n";
        std::cout<<"size Tp = "<<M_Tp.size()<<"\n";
        std::cout<<"M_num_elements = "<<M_num_elements<<"\n";
        throw std::runtime_error("Wave arrays wrong size\n");
    }

    while (lcpt < M_num_timesteps)
    {
        std::cout <<  ":[WIM2D TIME STEP]^"<< lcpt+1
           <<" (out of "<<M_num_timesteps<<")"<<"\n";

        //export progress results
        bool export_now = false;
        int dump_freq   = vm["wimdiag.dumpfreq"].template as<int>();//output "prog" and diagnostic text files after this no of time steps
        if (dump_freq>0)
            export_now = !(M_cpt % dump_freq);
        bool exportProg = export_now && (vm["wimdiag.checkprog"].template as<bool>());
        if ( exportProg )
            this->exportResults("prog");

        //integrate model
        M_dump_diag = export_now && (M_itest>0); 
        this->timeStep();

        ++lcpt;//local counter incremented here now
        ++M_cpt;//global counter incremented here now
        //std::cout<<"lcpt,M_cpt = "<<lcpt<<","<<M_cpt<<"\n";

        M_current_time  = this->getModelTime(lcpt);//model time of current call to wim (relative to init_time)
    }


    if (vm["wimdiag.checkfinal"].template as<bool>())
       this->exportResults("final");

    // save diagnostic file
    if (vm["wimdiag.savelog"].template as<bool>())
       this->saveLog(M_update_time);

    if (M_wim_on_mesh)
        //set M_UM to zero again
        std::fill(M_UM.begin(),M_UM.end(),0.);

    M_update_time   = M_current_time;//next time run is called, lcpt will be relative to this

    std::cout<<"Running done in "<< chrono.elapsed() <<"s\n";

    std::cout << "-----------------------Simulation completed at "<< current_time_local() <<"\n";
}//run


template<typename T>
void WimDiscr<T>::advectDirections(T_val_vec2d& Sdir,T_val_vec const& ag2d_eff)
{

	for (int nth = 0; nth < nwavedirn; nth++)
    {
        T_val adv_dir      = -PI*(90.0+M_wavedir[nth])/180.0;
        T_val_vec uwave   = ag2d_eff;
        T_val_vec vwave   = ag2d_eff;

        //set wave speeds
        //TODO if M_wim_on_mesh, subtract average mesh velocity
        std::for_each(uwave.begin(), uwave.end(), [&](T_val& f){ f *= std::cos(adv_dir); });
        if (M_advdim == 2)
            std::for_each(vwave.begin(), vwave.end(), [&](T_val& f){ f *= std::sin(adv_dir); });

        //do advection
        //TODO if M_wim_on_mesh call MeshTools::advect here
        M_grid.waveAdvWeno(Sdir[nth],uwave,vwave,M_timestep);
    }//advection of each direction done

}//advectDirections()


template<typename T>
void WimDiscr<T>::advectDirectionsMesh(T_val_vec2d& Sdir,T_val_vec & agnod,
        T_val_vec const &boundary_vals)
{

    int Nnod = M_mesh.M_num_nodes;
#if 0
    std::cout<<"advectDirectionsMesh: calling testMesh\n";
    this->testMesh();
#endif
    T_val* advect_out;
    int nb_var  = 1;                    //have to advect 1 vbl at a time
    std::vector<int> adv_method = {1};  //alternative (0) is do nothing

    auto test_vec = Sdir[0];//energy pre-advection

    //advect the directions
	for (int nth = 0; nth < nwavedirn; nth++)
    {
        T_val adv_dir = -PI*(90.0+M_wavedir[nth])/180.0;
        T_val_vec VC(2*Nnod,0.);
        T_val_vec bvals    = {boundary_vals[nth]};

        // set wave speeds
        // - subtract average mesh velocity
        // (average over length of call to wim = "M_duration")
#pragma omp parallel for num_threads(M_max_threads) collapse(1)
        for (int i=0;i<Nnod;i++)
        {
            VC[i]      = agnod[i]*std::cos(adv_dir)-M_UM[i]/M_duration;
            VC[i+Nnod] = agnod[i]*std::sin(adv_dir)-M_UM[i+Nnod]/M_duration;
        }

        //do advection
        //std::cout<<"advectDirectionsMesh: calling MeshTools::advect()\n";
        M_mesh.advect(&advect_out,&(Sdir[nth])[0],&VC[0],
                &adv_method[0],nb_var,M_timestep,&bvals[0]);

        // copy from 2D temporary array back to 3D input array
#pragma omp parallel for num_threads(M_max_threads) collapse(1)
        for (int i = 0; i < M_num_elements; i++)
            Sdir[nth][i] = advect_out[i];
    }//advection of each direction done

    xDelete<T_val>(advect_out);

#if 0
    std::cout<<"export: test advection\n";

    //choose the variables
    T_map_vec_ptrs extract_fields;
    extract_fields.emplace("agnod",&agnod);
    extract_fields.emplace("M_UM",&M_UM);
    extract_fields.emplace("E_pre",&(test_vec));
    extract_fields.emplace("E_post",&(Sdir[0]));

    //filenames
    std::string pathstr = vm["wimdiag.outparentdir"].template as<std::string>();
    pathstr += "/binaries/test_advection";
    fs::path path(pathstr);
    if ( !fs::exists(path) )
        fs::create_directories(path);
    std::string mfile  = (boost::format(  "%1%/mesh_%2%" ) % pathstr % M_cpt).str();
    std::string ffile  = (boost::format( "%1%/field_%2%" ) % pathstr % M_cpt).str();
    std::vector<std::string> filenames = {mfile,ffile};

    //export
    this->exportResultsMesh(extract_fields,filenames);//TODO fix time to nextsim standard
#endif

}//advectDirectionsMesh()


template<typename T>
void WimDiscr<T>::intWaveSpec()
{
    std::fill( M_mwd_x.begin()          , M_mwd_x.end()         , 0. );
    std::fill( M_mwd_y.begin()          , M_mwd_y.end()         , 0. );
    std::fill( M_stokes_drift_x.begin() , M_stokes_drift_x.end(), 0. );
    std::fill( M_stokes_drift_y.begin() , M_stokes_drift_y.end(), 0. );

    // NB mom0,mom2 are local to this routine
    // - not so important here since intWaveSpec() is only called at regridding time
    //   (ie not so often)
    T_val_vec mom0(M_num_elements,0.);
    T_val_vec mom2(M_num_elements,0.);
    for (int fq=0;fq<nwavefreq;fq++)
    {
        T_val om    = 2*PI*M_freq_vec[fq];    // radial freq
        this->intDirns(M_sdf_dir[fq], Mtmp_sdf_freq,
                Mtmp_stokes_drift_x_om, Mtmp_stokes_drift_y_om);

        for (int i=0;i<M_num_elements;i++)
        {
            T_val kicel = 2*PI/M_wlng_ice[fq][i]; // ice wave number (just water wavelenght if no ice)
            T_val F     = M_disp_ratio   [fq][i]; // convert from water amp's to ice amp's
            T_val F2    = 1.;
            if (M_ref_Hs_ice)
                F2  = std::pow(F,2);//for outputs only

            //integrals for Hs,Tp
            T_val tmp1  = M_quadrature_wt_freq[fq]*F2*Mtmp_sdf_freq[i];
            mom0[i] += tmp1;
            mom2[i] += tmp1*std::pow(om,2);

            //integrals for MWD
            M_mwd_x[i] += M_quadrature_wt_freq[fq]*F2*Mtmp_stokes_drift_x_om[i];
            M_mwd_y[i] += M_quadrature_wt_freq[fq]*F2*Mtmp_stokes_drift_y_om[i];

            //2*\omega*k*\int_0^\pi S(\omega,\theta)\cos(\theta)d\theta
            tmp1 = 2*om*kicel*F2*Mtmp_stokes_drift_x_om[i];
            M_stokes_drift_x[i] += M_quadrature_wt_freq[fq]*tmp1;

            //2*\omega*k*\int_0^\pi S(\omega,\theta)\sin(\theta)d\theta
            tmp1 = 2*om*kicel*F2*Mtmp_stokes_drift_y_om[i];
            M_stokes_drift_y[i] += M_quadrature_wt_freq[fq]*tmp1;
        }//i loop
    }//freq loop

    std::fill( M_Tp.begin() , M_Tp.end() , 0. );
    std::fill( M_mwd.begin(), M_mwd.end(), 0. );
    for (int i=0;i<M_num_elements;i++)
    {
        M_Hs[i]   = 4*std::sqrt(mom0[i]);
        if(mom2[i]>0)
        {
            M_Tp[i]   = 2*PI*std::sqrt(mom0[i]/mom2[i]);
            M_mwd[i]  = -90.-(180./PI)*std::atan2(M_mwd_y[i],M_mwd_x[i]);
        }
    }
}//intWaveSpec()


template<typename T>
void WimDiscr<T>::intDirns(T_val_vec2d const& Sdir, T_val_vec& Sfreq,
        T_val_vec& sdx_omega, T_val_vec& sdy_omega)
{

    std::fill( Sfreq.begin()      ,Sfreq.end()      ,0. );
    std::fill( sdx_omega.begin()  ,sdx_omega.end()  ,0. );
    std::fill( sdy_omega.begin()  ,sdy_omega.end()  ,0. );

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {

        // integrals to be done everywhere (not just in ice)
        for (int nth = 0; nth < nwavedirn; nth++)
        {
            //frequency spectrum
            Sfreq[i] += M_quadrature_wt_dir[nth]*Sdir[nth][i];

            //Stoke's drift at surface (x)
            T_val adv_dir = -PI*(90.0+M_wavedir[nth])/180.0;
            T_val tmp     = std::cos(adv_dir)*M_quadrature_wt_dir[nth]*Sdir[nth][i];
            sdx_omega[i] += tmp;

            //Stoke's drift at surface (y)
            tmp = std::sin(adv_dir)*M_quadrature_wt_dir[nth]*Sdir[nth][i];
            sdy_omega[i] += tmp;
        }

#if 0
        if (i==M_itest)
        {
            std::cout<<"i = "<<i<<"\n";
            std::cout<<"Mtmp_sdf_freq = "<<Mtmp_sdf_freq[i]<<"\n";
        }
#endif
    }
}//intDirns


template<typename T>
void WimDiscr<T>::attenSimple(T_val_vec2d& Sdir, T_val_vec& Sfreq,
        T_val_vec& taux_omega, T_val_vec& tauy_omega,
        T_val_vec& sdx_omega, T_val_vec& sdy_omega,
        T_val_vec const& ag2d_eff)
{

	T_val S_th, tmp, alp_dim, source;

    std::fill( Sfreq.begin()      ,Sfreq.end()      ,0. );
    std::fill( taux_omega.begin() ,taux_omega.end() ,0. );
    std::fill( tauy_omega.begin() ,tauy_omega.end() ,0. );
    std::fill( sdx_omega.begin()  ,sdx_omega.end()  ,0. );
    std::fill( sdy_omega.begin()  ,sdy_omega.end()  ,0. );

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < M_num_elements; i++)
    {
        if (M_ice[IceType::wim].M_mask[i] > .5)
        {
            for (int nth = 0; nth < nwavedirn; nth++)
            {
                T_val adv_dir = -PI*(90.0+M_wavedir[nth])/180.0;
                T_val S_th = Sdir[nth][i];
                T_val alp_dim = Mtmp_atten_dim[i]+Mtmp_damp_dim[i];

                // stress calculation
                T_val source = -alp_dim*ag2d_eff[i]*S_th;
                T_val tmp = -std::cos(adv_dir)*M_quadrature_wt_dir[nth]*source;
                taux_omega[i] += tmp;
                tmp = -std::sin(adv_dir)*M_quadrature_wt_dir[nth]*source;
                tauy_omega[i] += tmp;
#if 0
                if (i==M_itest)
                {
                    std::cout<<"i,nth,adv_dir = "<<i<<","<<nth<<","<<adv_dir<<"\n";
                    std::cout<<"Mtmp_atten_dim,damp_dim = "<<Mtmp_atten_dim[i]<<","<<Mtmp_damp_dim[i]<<"\n";
                    std::cout<<"alp_dim,source,M_quadrature_wt_dir = "<<alp_dim<<","<<source<<","<<M_quadrature_wt_dir[nth]<<"\n";
                    std::cout<<"taux_omega = "<<taux_omega[i]<<"\n";
                    std::cout<<"tauy_omega = "<<tauy_omega[i]<<"\n";
                }
#endif

                // do attenuation
                Sdir[nth][i] = S_th*std::exp(-alp_dim*ag2d_eff[i]*M_timestep);

                //std::cout<<"tau_x["<< i << "]= "<< Mtmp_atten_dim[i] <<"\n";
            }//loop over directions
        }//end "if ice"

        // integrals to be done everywhere (not just in ice)
        for (int nth = 0; nth < nwavedirn; nth++)
        {
            //frequency spectrum
            Sfreq[i] += M_quadrature_wt_dir[nth]*Sdir[nth][i];

            //Stoke's drift at surface (x)
            T_val adv_dir = -PI*(90.0+M_wavedir[nth])/180.0;
            T_val tmp     = std::cos(adv_dir)*M_quadrature_wt_dir[nth]*Sdir[nth][i];
            sdx_omega[i] += tmp;

            //Stoke's drift at surface (y)
            tmp = std::sin(adv_dir)*M_quadrature_wt_dir[nth]*Sdir[nth][i];
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
void WimDiscr<T>::attenIsotropic(T_val_vec2d& Sdir, T_val_vec& Sfreq,
        T_val_vec& taux_omega, T_val_vec& tauy_omega,
        T_val_vec& sdx_omega, T_val_vec& sdy_omega,
        T_val_vec const& ag2d_eff)
{
    std::vector<T_val> nvec(nwavedirn);
	std::vector<T_val> K_fou(nwavedirn), S_th(nwavedirn), theta_vec(nwavedirn);
	std::vector<T_val> tmp1(nwavedirn), evals_x(nwavedirn);
	T_val tmp, alp_dim, source;

	std::vector<std::complex<T_val> > S_fou(nwavedirn);
	std::complex<T_val> zi, src_fou_p1, src_fou_m1;

    //no of cosines/sines to use
    int ncs = std::round(nwavedirn/2);
	std::vector<T_val> S_cos(ncs), S_sin(ncs);
	T_val cg, q_scat, q_abs, q_tot, src_cos_1, src_sin_1;
    int n, jp1, jm1;

	zi = std::complex<T_val>(0.,1.);

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

        //S_fou[0] = std::complex<T_val>( sum(M_quadrature_wt_dir*S_th) );
        S_fou[0] = std::complex<T_val>( std::inner_product(M_quadrature_wt_dir.begin(), M_quadrature_wt_dir.end(), S_th.begin(), 0.) );

        if (M_ice[IceType::wim].M_mask[i] > 0.5 )
        {
            if (M_ice[IceType::wim].M_dfloe[i] < M_dfloe_pack_init)
            {
                q_scat = Mtmp_atten_dim[i];
                q_abs = Mtmp_damp_dim[i];
            }
            else
            {
                q_scat = 0;
                q_abs = Mtmp_atten_dim[i]+Mtmp_damp_dim[i];
            }

            q_tot = q_scat+q_abs;
            cg = ag2d_eff[i];

            std::fill(K_fou.begin(), K_fou.end(), 0.);
            K_fou[0] = q_scat;

            std::fill(evals_x.begin(), evals_x.end(), -q_tot);
            evals_x[0] = -q_abs;


            for (int nth = 0; nth < ncs; nth++)
            {
                std::vector<T_val> prodtmp = theta_vec;

                std::for_each(prodtmp.begin(), prodtmp.end(), [&](T_val& f){ f = std::cos((nth+1)*f); });

                std::transform(prodtmp.begin(), prodtmp.end(), S_th.begin(), prodtmp.begin(),
                              std::multiplies<T_val>());

                S_cos[nth] = std::inner_product(prodtmp.begin(), prodtmp.end(), M_quadrature_wt_dir.begin(), 0.);

                prodtmp.clear();
                prodtmp.resize(theta_vec.size());
                prodtmp = theta_vec;
                std::for_each(prodtmp.begin(), prodtmp.end(), [&](T_val& f){ f = std::sin((nth+1)*f); });

                std::transform(prodtmp.begin(), prodtmp.end(), S_th.begin(), prodtmp.begin(),
                               std::multiplies<T_val>());

                S_sin[nth] = std::inner_product(prodtmp.begin(), prodtmp.end(), M_quadrature_wt_dir.begin(), 0.);

                S_fou[nth+1] = std::complex<T_val>(S_cos[nth],S_sin[nth]);
                //S_fou[nwavedirn-nth] = std::complex<T_val>(S_cos[nth],-S_sin[nth]);

                if (nth != ncs-1)
                {
                    S_fou[nwavedirn-(nth+1)] = std::complex<T_val>(S_cos[nth],-S_sin[nth]);
                }

                // std::cout<<"nth= " << nth << ": and ncs+nth= "<< ncs+nth <<"\n";
                // std::cout<<"taux= "<< S_sin[nth] <<"\n";
            }

            // stresses
            jp1 = 1;
            jm1 = nwavedirn-1;

            src_fou_p1 = cg*(-q_tot*S_fou[jp1]+K_fou[jp1]*S_fou[jp1]);
            src_fou_m1 = cg*(-q_tot*S_fou[jm1]+K_fou[jm1]*S_fou[jm1]);
            src_cos_1 = std::real(std::complex<T_val>(0.5)*(src_fou_p1+src_fou_m1));
            src_sin_1 = std::real(-zi*std::complex<T_val>(0.5)*(src_fou_p1-src_fou_m1));

            taux_omega[i] = -src_cos_1;
            tauy_omega[i] = -src_sin_1;

            sdx_omega[i] = std::real(.5*(S_fou[jp1]+S_fou[jm1]));
            sdy_omega[i] = std::real(-.5*zi*(S_fou[jp1]-S_fou[jm1]));

            // if (i==1)
            // {
            //     std::cout<<"taux_omega= "<< taux_omega[i] <<"\n";
            //     std::cout<<"tauy_omega= "<< tauy_omega[i] <<"\n";
            // }

            std::vector<T_val> prodtmp = evals_x;
            std::for_each(prodtmp.begin(), prodtmp.end(), [&](T_val& f){ f = std::exp(cg*M_timestep*f); });
            std::transform(S_fou.begin(), S_fou.end(), prodtmp.begin(), S_fou.begin(),
                           std::multiplies<std::complex<T_val> >());

            std::vector<T_val> Sfoutempcos(S_fou.size());// = S_fou;
            std::vector<T_val> Sfoutempsin(S_fou.size());// = S_fou;

            for (int nth = 0; nth < nwavedirn; nth++)
            {
                // prodtmp = std::vector<T_val>(nvec.begin(), nvec.end()); // also works

                for (int ss=0; ss<S_fou.size(); ++ss)
                {
                    Sfoutempcos[ss] = std::real(S_fou[ss]);
                    Sfoutempsin[ss] = std::imag(S_fou[ss]);
                }

                prodtmp.clear();
                prodtmp.resize(nwavedirn);
                prodtmp = nvec;

                std::for_each(prodtmp.begin(), prodtmp.end(), [&](T_val& f){ f = std::cos(theta_vec[nth]*f); });

                std::transform(Sfoutempcos.begin(), Sfoutempcos.end(), prodtmp.begin(), Sfoutempcos.begin(),
                               std::multiplies<T_val>());


                // prodtmp = std::vector<T_val>(nvec.begin(), nvec.end()); // also works
                prodtmp.clear();
                prodtmp.resize(nwavedirn);
                prodtmp = nvec;

                std::for_each(prodtmp.begin(), prodtmp.end(), [&](T_val& f){ f = std::sin(theta_vec[nth]*f); });

                std::transform(Sfoutempsin.begin(), Sfoutempsin.end(), prodtmp.begin(), Sfoutempsin.begin(),
                               std::multiplies<T_val>());

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

}//attenIsotropic


template<typename T>
typename WimDiscr<T>::T_val
WimDiscr<T>::thetaDirFrac(T_val const& th1_, T_val const& dtheta_, T_val const& mwd_) const
{
    //get mwd\pm90 inside [th1_,th1_+360)
    T_val phi1  = thetaInRange(mwd_-90.,th1_);//>=th1_
    T_val phi2  = thetaInRange(mwd_+90.,th1_);//>=th1_
    T_val th2_  = th1_+dtheta_;

    T_val integral = 0.;
    if (phi2>phi1)
    {
       // th1_,phi1,phi2, and th2_
       T_val L1    = std::max(th1_,phi1);
       T_val L2    = std::min(th2_,phi2);
       L2               = std::max(L1,L2); //make L2>=L1
       T_val chi1  = PI*(L1-mwd_)/180.;
       T_val chi2  = PI*(L2-mwd_)/180.;
       integral        += 2.*(chi2-chi1)+std::sin(2*chi2)-std::sin(2*chi1);
    }
    else
    {
       // th1_,phi2,phi1, and th2_
       // 1st consider (th1_,phi2) interval
       T_val L1    = th1_;
       T_val L2    = min(th2_,phi2);
       T_val chi1  = PI*(L1-mwd_)/180.;
       T_val chi2  = PI*(L2-mwd_)/180.;
       integral        += 2.*(chi2-chi1)+std::sin(2*chi2)-std::sin(2*chi1);
      
       // 2nd consider (phi1,th2_) interval
       L1           = phi1;
       L2           = max(L1,th2_);      //make L2>=L1
       chi1         = PI*(L1-mwd_)/180.;
       chi2         = PI*(L2-mwd_)/180.;
       integral    += 2.*(chi2-chi1)+std::sin(2*chi2)-std::sin(2*chi1);
    }

    T_val theta_dirfrac  = integral/(2.*PI);
    return theta_dirfrac;

}//thetaDirFrac


template<typename T>
typename WimDiscr<T>::T_val
WimDiscr<T>::thetaInRange(T_val const& th_, T_val const& th1, bool const& close_on_right) const
{
    T_val th2, dth, th;
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
void WimDiscr<T>::exportResults(std::string const& output_type)
{

    T_map_vec_ptrs extract_fields;

    std::string pathstr = vm["wimdiag.outparentdir"].template as<std::string>();
    pathstr += "/binaries/"+output_type;

    std::string prefix = output_type;
    int step = 0;

    //select output_type
    if ( output_type == "prog" )
    {
        prefix  = "wim_prog";

        //fields to extract
        extract_fields.emplace("MWD"           ,&(M_mwd));
        extract_fields.emplace("Tp"            ,&(M_Tp));
        extract_fields.emplace("Hs"            ,&(M_Hs));
        extract_fields.emplace("stokes_drift_y",&(M_stokes_drift_y));
        extract_fields.emplace("stokes_drift_x",&(M_stokes_drift_x));
        extract_fields.emplace("tau_y"         ,&(M_tau_y));
        extract_fields.emplace("tau_x"         ,&(M_tau_x));
        extract_fields.emplace("Dmax"          ,&(M_ice[IceType::wim].M_dfloe));

        step = M_nb_export_prog;
        M_nb_export_prog++;
    }
    else if ( output_type == "final" )
    {
        prefix  = "wim_out";

        //fields to extract
        extract_fields.emplace("MWD"           ,&(M_mwd));
        extract_fields.emplace("Tp"            ,&(M_Tp));
        extract_fields.emplace("Hs"            ,&(M_Hs));
        extract_fields.emplace("stokes_drift_y",&(M_stokes_drift_y));
        extract_fields.emplace("stokes_drift_x",&(M_stokes_drift_x));
        extract_fields.emplace("tau_y"         ,&(M_tau_y));
        extract_fields.emplace("tau_x"         ,&(M_tau_x));
        extract_fields.emplace("Dmax"          ,&(M_ice[IceType::wim].M_dfloe));

        step = M_nb_export_final;
        M_nb_export_final++;
    }
    else if ( output_type == "init" )
    {
        prefix  = "wim_init";

        //fields to extract
        extract_fields.emplace("MWD" ,&(M_mwd));
        extract_fields.emplace("Tp"  ,&(M_Tp));
        extract_fields.emplace("Hs"  ,&(M_Hs));
        extract_fields.emplace("Dmax",&(M_ice[IceType::wim].M_dfloe));
        extract_fields.emplace("iceh",&(M_ice[IceType::wim].M_thick));
        extract_fields.emplace("icec",&(M_ice[IceType::wim].M_conc));

        step = M_nb_export_init;
        M_nb_export_init++;
    }
    else if ( output_type == "incwaves" )
    {
        prefix  = "wim_inc";

        //fields to extract
        extract_fields.emplace("MWD",&(M_mwd_in));
        extract_fields.emplace("Tp" ,&(M_mwp_in));
        extract_fields.emplace("Hs" ,&(M_swh_in));

        step = M_nb_export_inc;
        M_nb_export_inc++;
    }
    else if ( output_type == "nextwim" )
    {
        //fields to extract
        extract_fields.emplace("MWD"           ,&(M_mwd));
        extract_fields.emplace("Tp"            ,&(M_Tp));
        extract_fields.emplace("Hs"            ,&(M_Hs));
        extract_fields.emplace("stokes_drift_y",&(M_stokes_drift_y));
        extract_fields.emplace("stokes_drift_x",&(M_stokes_drift_x));
        extract_fields.emplace("tau_y"         ,&(M_tau_y));
        extract_fields.emplace("tau_x"         ,&(M_tau_x));
        extract_fields.emplace("Dmax"          ,&(M_ice[IceType::wim].M_dfloe));

        step = M_nb_export_nextwim;
        M_nb_export_nextwim++;
    }//select output_type

    fs::path path(pathstr);
    if ( !fs::exists(path) )
        fs::create_directories(path);

    if((!M_wim_on_mesh)&&(extract_fields.size()>0))
    {
        std::string init_time  = Wim::ptime(M_init_time_str);
        std::string timestpstr = Wim::ptime(M_init_time_str, M_current_time);
        std::string fileout    = (boost::format( "%1%/%2%%3%" ) % pathstr % prefix % timestpstr).str();
        std::vector<std::string> export_strings = {fileout,init_time,timestpstr};
        this->exportResultsGrid(extract_fields,export_strings);
    }
    else
    {
        std::string mfile  = (boost::format(  "%1%/mesh_%2%" ) % pathstr % step).str();
        std::string ffile  = (boost::format( "%1%/field_%2%" ) % pathstr % step).str();
        std::vector<std::string> filenames = {mfile,ffile};
        this->exportResultsMesh(extract_fields,filenames);
    }
}//exportResults


template<typename T>
void WimDiscr<T>::exportResultsGrid(T_map_vec_ptrs & extract_fields,
        std::vector<std::string> const& strings)
{

    // ==========================================================================================
    //filenames
    std::string fileout     = strings[0]+".a";
    std::string fileoutb    = strings[0]+".b";
    std::string init_time   = strings[1];
    std::string timestpstr  = strings[2];

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
    // ==========================================================================================

    // ==================================================================================================
    //.b file header
    int Nrecs = extract_fields.size();
    std::string rstr   = std::string(2-std::to_string(Nrecs).length(),'0')
                          + std::to_string(Nrecs);
    outb << std::setw(15) << std::left << rstr  << "    Nrecs    # "<< "Number of records" <<"\n";
    outb << std::setw(15) << std::left << "0"   << "    Norder   # "
         << "Storage order [column-major (F/matlab) = 1; row-major (C) = 0]" <<"\n";
    outb << std::setw(15) << std::left << M_grid.M_num_px    << "    nx       # "<< "Record length in x direction (elements)" <<"\n";
    outb << std::setw(15) << std::left << M_grid.M_num_py    << "    ny       # "<< "Record length in y direction (elements)" <<"\n";

    outb <<"\n";
    outb << std::left << init_time << "    t_start    # "<< "Model time of WIM call" <<"\n";
    outb << std::left << timestpstr << "    t_out    # "<< "Model time of output" <<"\n";

    outb <<"\n";
    outb << "Record number and name:" <<"\n";
    // ==================================================================================================

    int recno = 0;
    for(auto it=extract_fields.begin();it!=extract_fields.end();it++)
    {
        auto vname = it->first;
        //std::cout<<vname<<"\n";
        //auto vtmp = *(it->second);
        //std::cout<<"vtmp size = "<<vtmp.size()<<"\n";
        //std::cout<<"M_num_elements = "<<M_num_elements<<"\n";
        for (int i = 0; i < M_num_elements; i++)
        {
            //data
            T_val tmp = (*(it->second))[i];//it->second points to a vector
            out.write((char *)&tmp, sizeof(T_val));
            //out.write((char *)&vtmp[i], sizeof(T_val));
        }

        //record
        ++recno;
        rstr   = std::string(2-std::to_string(recno).length(),'0')
           + std::to_string(recno);
        outb << std::setw(9) << rstr << vname <<"\n";
    }

    out.close();
    outb.close();
}//exportResultsGrid


template<typename T>
void WimDiscr<T>::exportResultsMesh(T_map_vec_ptrs & extract_fields,
        std::vector<std::string> const &filenames,
        bool export_mesh, bool export_fields)
{

    //export the mesh
    if (export_mesh)
        this->exportMesh(filenames[0]);

    //can stop here if not exporting the fields
    if (!export_fields)
        return;

    //need to convert some names from WIM to neXtSIM conventions
    T_map dict;
    dict.emplace("icec","Concentration");
    dict.emplace("iceh","Thickness");
    dict.emplace("Dmax","Dfloe");

    std::string fileout;
    fileout = filenames[1]+".bin";
    std::fstream outbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
    if ( ! outbin.good() )
        throw std::runtime_error("Cannot write to file: " + fileout);

    //construct exporter
    Nextsim::Exporter exporter(vm["output.exporter_precision"].as<std::string>());
    std::vector<double> timevec = {this->getNextsimTime()};
    exporter.writeField(outbin, timevec, "Time");
    exporter.writeField(outbin, M_mesh.M_surface, "Element_area");

    //loop over the fields and write them
    for (auto it=extract_fields.begin();it!=extract_fields.end();it++)
    {
        auto vname = it->first;//"first" is a string
        if(dict.count(vname))
            vname=dict[vname];
        exporter.writeField(outbin, *(it->second), vname);//"second" points to a vector
    }
    outbin.close();

    //write the .dat file
    fileout = filenames[1]+".dat";
    std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
    if ( ! outrecord.good() )
        throw std::runtime_error("Cannot write to file: " + fileout);

    exporter.writeRecord(outrecord);
    outrecord.close();

}//exportResultsMesh()


template<typename T>
void WimDiscr<T>::testMesh()
{
    std::string pathstr = vm["wimdiag.outparentdir"].template as<std::string>();
    pathstr += "/binaries/test_mesh";
    fs::path path(pathstr);
    if ( !fs::exists(path) )
        fs::create_directories(path);

    std::string filename = (boost::format( "%1%/mesh_%2%" ) % pathstr  % M_nb_mesh_test).str();
    this->exportMesh(filename);
    M_nb_mesh_test++;
}


template<typename T>
void WimDiscr<T>::exportMesh(std::string const &filename)
{
    Nextsim::Exporter exporter(vm["output.exporter_precision"].as<std::string>());
    std::string fileout;

    fileout = filename+".bin";
    std::fstream meshbin(fileout, std::ios::binary | std::ios::out | std::ios::trunc);
    if ( ! meshbin.good() )
        throw std::runtime_error("Cannot write to file: " + fileout);

    exporter.writeMesh(meshbin, M_mesh.M_nodes_x, M_mesh.M_nodes_y,
            M_mesh.M_id,M_mesh.M_index);
    meshbin.close();

    fileout = filename+".dat";
    std::fstream outrecord(fileout, std::ios::out | std::ios::trunc);
    if ( ! outrecord.good() )
        throw std::runtime_error("Cannot write to file: " + fileout);

    exporter.writeRecord(outrecord,"mesh");
    outrecord.close();
}//exportMesh()


template<typename T>
void WimDiscr<T>::saveLog(T_val const& t_out) const
{
    std::string str = vm["wimdiag.outparentdir"].template as<std::string>();
    fs::path path(str);
    path /= "diagnostics/global";
    if ( !fs::exists(path) )
       fs::create_directories(path);

    std::string init_time  = ptime(M_init_time_str);
    std::string timestpstr = ptime(M_init_time_str, t_out);
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
    out << std::left << std::setw(log_width) << "SCATMOD"<<" : " << M_scatmod << "\n";
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
    out << std::left << std::setw(log_width) << "DO_ATTEN"<<" : " << M_atten << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "Other integer parameters:" << "\n";
    out << std::left << std::setw(log_width) << "FSD_OPT"<<" : " << M_ice_params.fsdopt() << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "WIM parameters:" << "\n";
    out << std::left << std::setw(log_width) << "Brine volume fraction"  <<" : " << M_ice_params.vbf() << "\n";
    out << std::left << std::setw(log_width) << "Youngs modulus (Pa)"    <<" : " << M_ice_params.young() << "\n";
    out << std::left << std::setw(log_width) << "Flexural strength (Pa)" <<" : " << M_ice_params.sigmac() << "\n";
//    out << std::left << std::setw(log_width) << "Breaking stress (Pa)"<<" : " << M_ice_params.stressc() << "\n";
    out << std::left << std::setw(log_width) << "Breaking strain"        <<" : " << M_ice_params.epsc() << "\n";
    out << std::left << std::setw(log_width) << "Damping (Pa.s/m)"       <<" : " << M_ice_params.dragRp() << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "FSD parameters:" << "\n";
    out << std::left << std::setw(log_width) << "Dmin (m)"  <<" : " << M_ice_params.dmin() << "\n";
    out << std::left << std::setw(log_width) << "xi"        <<" : " << M_ice_params.xi() << "\n";
    out << std::left << std::setw(log_width) << "fragility" <<" : " << M_ice_params.fragility() << "\n";
    out << std::left << std::setw(log_width) << "Dthresh"   <<" : " << M_ice_params.dfloeMizThresh() << "\n";
    out << std::left << std::setw(log_width) << "cice_min"  <<" : " << M_ice_params.ciceMin() << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    out << "Other parameters:" << "\n";
    out << std::left << std::setw(log_width) << "Time step (s)"             <<" : " << M_timestep << "\n";
    out << std::left << std::setw(log_width) << "CFL length (km)"           <<" : " << M_length_cfl/1.e3 << "\n";
    out << std::left << std::setw(log_width) << "CFL number"                <<" : " << M_cfl << "\n";
    out << std::left << std::setw(log_width) << "Max wave group vel (m/s)"  <<" : " << M_max_cg << "\n";
    out << std::left << std::setw(log_width) << "Number of time steps"      <<" : " << M_num_timesteps << "\n";
    out << std::left << std::setw(log_width) << "Time interval (h)"         <<" : " << M_duration/60.0/60.0 << "\n";
    out << "***********************************************\n";

    out << "\n***********************************************\n";
    if(!M_wim_on_mesh)
    {
        out << std::left << std::setw(log_width) << "Grid dimensions"        <<" : " << M_grid.M_num_px << ", " << M_grid.M_num_py << "\n";
        out << std::left << std::setw(log_width) << "Spatial resolution (km)"<<" : " << M_grid.M_dx/1.0e3 << ", " << M_grid.M_dy/1.0e3 << "\n";
        out << std::left << std::setw(log_width) << "Extent of domain (km)"  <<" : "
            << (M_grid.M_xmax-M_grid.M_xmin)/1.e3 << ", " << (M_grid.M_ymax-M_grid.M_ymin)/1.e3 << "\n";
    }
    out << std::left << std::setw(log_width) << "Minimum period (s)"          <<" : " << 1.0/M_freq_vec[nwavefreq-1] << "\n";
    out << std::left << std::setw(log_width) << "Maximum period (s)"          <<" : " << 1.0/M_freq_vec[0] << "\n";
    out << std::left << std::setw(log_width) << "Number of wave frequencies"  <<" : " << nwavefreq << "\n";
    out << std::left << std::setw(log_width) << "Number of wave directions"   <<" : "  << nwavedirn << "\n";
    out << std::left << std::setw(log_width) << "Directional resolution (deg)"<<" : " << 360.0/nwavedirn << "\n";
    out << "***********************************************\n";

    T_val taux_min  = *std::min_element(M_tau_x.begin(), M_tau_x.end());
    T_val taux_max  = *std::max_element(M_tau_x.begin(), M_tau_x.end());
    T_val tauy_min  = *std::min_element(M_tau_y.begin(), M_tau_y.end());
    T_val tauy_max  = *std::max_element(M_tau_y.begin(), M_tau_y.end());
    T_val Hs_max    = *std::max_element(M_Hs.begin()   , M_Hs.end());

    //MIZ diagnostics
    T_val Dmax_min = 10.e3;
    T_val Dmax_max = 0.e3;
    int Nmiz   = 0;

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i = 0; i < (M_ice[IceType::wim].M_dfloe).size(); i++)
       if ((M_ice[IceType::wim].M_dfloe[i]<M_dfloe_pack_init)&&(M_ice[IceType::wim].M_dfloe[i]>0))
       {
          ++Nmiz;
          Dmax_min   = std::min(M_ice[IceType::wim].M_dfloe[i],Dmax_min);
          Dmax_max   = std::max(M_ice[IceType::wim].M_dfloe[i],Dmax_max);
       }

//#define WIMDIAG1D
#if defined (WIMDIAG1D)
    //this definition of MIZ only works in 1d geometries
    //and calculation only works for wim on grid
    T_val W_miz;
    if(!M_wim_on_mesh)
    {
        //this definition of MIZ only works in 1d geometries
        if ( M_grid.M_num_py == 1 )
            W_miz = (Nmiz*M_grid.M_dx);
        else if ( vm["wimgrid.landon3edges"].template as<bool>() )
            W_miz = (Nmiz*M_grid.M_dx)/(M_grid.M_num_py-2);
        else
            W_miz = (Nmiz*M_grid.M_dx)/M_grid.M_num_py;
    }
#endif

    out << "\n***********************************************\n";
    out << "Diagnostics:\n";
#if defined (WIMDIAG1D)
    if(!M_wim_on_mesh)
        out << std::left << std::setw(log_width) << "MIZ width (km)"<<" : "
            << W_miz/1.e3 << "\n";
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
    std::string fileout = vm["wimdiag.outparentdir"].template as<std::string>();
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


template<typename T>
typename WimDiscr<T>::T_val
WimDiscr<T>::getNextsimTime() const
{
    T_val t0 = Wim::dateStr2Num(M_init_time_str); //days from ref time (1901-1-1) to init_time
    return t0+M_current_time/(24*3600.);          //days from ref time (1901-1-1) to model time
}


template<typename T>
void WimDiscr<T>::getRangeXY(T_val &xmin,T_val &xmax,T_val &ymin,T_val &ymax) const
{
    if(!M_wim_on_mesh)
    {
        this->getRange(M_grid.M_px,xmin,xmax);
        this->getRange(M_grid.M_py,ymin,ymax);
    }
    else
    {
        this->getRange(M_mesh.M_elements_x,xmin,xmax);
        this->getRange(M_mesh.M_elements_y,ymin,ymax);
    }
}


template<typename T>
void WimDiscr<T>::getRange(T_val_vec const &vec, T_val &vmin, T_val &vmax) const
{

    vmin = *std::min_element(vec.begin(),vec.end());
    vmax = *std::max_element(vec.begin(),vec.end());
}

template<typename T>
void WimDiscr<T>::printRange(std::string const &name,T_val_vec const &vec, int const & prec) const
{

    T_val vmin = *std::min_element(vec.begin(),vec.end());
    T_val vmax = *std::max_element(vec.begin(),vec.end());
    if(prec==0)
    {
        //automatic precision
        std::cout<<"Min "<<name<<" = "<<vmin<<"\n";
        std::cout<<"Max "<<name<<" = "<<vmax<<"\n";
    }
    else
    {
        //specified precision
        std::cout<<"Min "<<name<<" = "<<std::setprecision(prec)<<vmin<<"\n";
        std::cout<<"Max "<<name<<" = "<<std::setprecision(prec)<<vmax<<"\n";
    }
}


// instantiate wim class for type float
//template class WimDiscr<float>;

// instantiate wim class for type double
template class WimDiscr<double>;

} // namespace WIM2D
