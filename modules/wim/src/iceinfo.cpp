/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   wimdiscr.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug  3 11:52:35 2015
 */

#include <iceinfo.hpp>
#ifdef __cplusplus
extern "C"
{
#endif
#include <RTparam_outer.h>
#ifdef __cplusplus
}
#endif

namespace Wim
{


template<typename T>
IceParams<T>::IceParams(po::variables_map const& vm)
{
    M_cice_min          = vm["wim.cicemin"].template as<double>();
    M_fsdopt            = vm["wim.fsdopt"].template as<std::string>();
    M_dfloe_pack_init   = vm["wim.dfloepackinit"].template as<double>();
    M_dfloe_pack_thresh = vm["wim.dfloepackthresh"].template as<double>();
    M_young             = vm["wim.young"].template as<double>();
    M_drag_rp           = vm["wim.dragrp"].template as<double>();

    this->init();
}


template<typename T>
void IceParams<T>::init()
{
    M_sigma_c  = (1.76e+6)*std::exp(-5.88*std::sqrt(M_vbf));//flexural strength (Pa)
    M_epsc = M_sigma_c/M_young;//breaking strain
}//end ::init()


template<typename T>
std::vector<double> IceParams<T>::getAttenParams()
{
    std::vector<double> params(5);
    params[0] = M_young;
    params[1] = M_gravity;
    params[2] = M_rhowtr;
    params[3] = M_rhoice;
    params[4] = M_poisson;
    return params;
}


template<typename T>
IceInfo<T>::IceInfo(T_icep ice_params, std::string const &name)
{
    M_initialized = true;
    M_max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
    M_ice_params  = ice_params;
    M_name        = name;
    std::cout<<"Init IceInfo(), name = "<<M_name<<"\n";
    std::cout<<"cicemin (IceInfo) = "<<M_ice_params.M_cice_min<<"\n";
}


template<typename T>
void IceInfo<T>::clearFields()
{
    // this needs to be called each time grid or mesh changes
    // ie initially, and if(M_wim_on_mesh), after regridding
    // set sizes of arrays, initialises some others that are constant in time

    M_conc.resize(0);
    M_vol.resize(0);
    M_nfloes.resize(0);
    M_dfloe.resize(0);
    M_thick.resize(0);
    M_mask.resize(0);
    M_broken.resize(0);
}//clearFields()

template<typename T>
void IceInfo<T>::setFields(T_val_vec const& conc,
            T_val_vec const& vol,
            T_val_vec const& nfloes)
{
    // this needs to be called each time grid or mesh changes
    // ie initially, and if(M_wim_on_mesh), after regridding
    // set sizes of arrays, initialises some others that are constant in time

    M_conc          = conc;
    M_vol           = vol;
    M_nfloes        = nfloes;
    this->updateFields();
}//setFields()


template<typename T>
void IceInfo<T>::updateFields()
{
    // this needs to be called each time grid or mesh changes
    // ie initially, and if(M_wim_on_mesh), after regridding
    // set sizes of arrays, initialises some others that are constant in time

    if(!M_initialized)
        throw std::runtime_error("IceInfo objected (name="+M_name+") not initialised yet - can't call updateFields()\n");
    M_num_elements  = M_conc.size();

    //2D var's
    M_dfloe.assign (M_num_elements ,0.);
    M_thick.assign (M_num_elements ,0.);
    M_broken.assign(M_num_elements ,0.);
    M_mask.assign  (M_num_elements ,0.);
    //std::cout<<"cicemin (updateFields) = "<<M_ice_params.M_cice_min<<"\n";

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i=0;i<M_num_elements;i++)
    {
        if (M_conc[i]>=M_ice_params.M_cice_min)
        {
            //set dependant variables
            M_thick[i]  = M_vol[i]/M_conc[i];//convert to actual thickness
            M_dfloe[i]  = this->nfloesToDfloe(M_nfloes[i],M_conc[i]);
            M_mask[i]   = 1.;
        }
        else
        {
            //set input variables to 0.
            M_conc[i]   = 0;
            M_vol[i]    = 0;
            M_nfloes[i] = 0;
        }

#if 1
        //check ranges of inputs
        if ((M_thick[i]<0.)||(M_thick[i]>50.))
        {
            std::cout<<"thickness out of range for i="<<i<<"\n";
            std::cout<<"vol,h,c="<<M_vol[i]<<","<<M_thick[i]<<","<<M_conc[i]<<"\n";
            throw std::runtime_error("thickness out of range");
        }
        if ((M_conc[i]<0.)||(M_conc[i]>1.))
        {
            std::cout<<"conc out of range for i="<<i<<"\n";
            std::cout<<"c="<<M_conc[i]<<"\n";
            throw std::runtime_error("conc out of range");
        }
#endif

        //std::cout<<"c,Nfloes,dfloe,mask = "<<M_conc[i]<<","<<M_nfloes[i]<<","<<M_dfloe[i]<<","<<M_mask[i]<<"\n";

    }//loop over elements

#if 1
    //test fields
    std::cout<<"IceInfo::updateFields ("<<M_name<<"): check ice fields\n";
    this->printRange( "conc      (iceinfo)" , M_conc   );
    this->printRange( "thickness (iceinfo)" , M_thick  );
    this->printRange( "Nfloes    (iceinfo)" , M_nfloes );
    this->printRange( "dfloe     (iceinfo)" , M_dfloe  );
    this->printRange( "M_mask    (iceinfo)" , M_mask   );
#endif

}//updateFields()

template<typename T>
void IceInfo<T>::doBreaking(T_val_vec const& mom0,T_val_vec const& mom2,T_val_vec const& var_strain)
{

    std::vector<double> params = M_ice_params.getAttenParams();//currently independant of space

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i=0;i<M_num_elements;i++)
    {
        bool break_criterion = (M_mask[i] > 0.5)          //ice present
            && ((2*var_strain[i]) > std::pow(M_ice_params.M_epsc, 2.));//big enough waves

        if (break_criterion)
        {
            double outputs[8];

            T_val om    = std::sqrt(mom2[i]/mom0[i]);
            T_val guess = std::pow(om,2.)/M_ice_params.M_gravity;

            RTparam_outer(outputs,M_thick[i],double(om),double(M_ice_params.M_drag_rp),double(guess),&params[0]);

            T_val kice = outputs[1];
            T_val lam = 2*PI/kice;

            if (lam < (2*M_dfloe[i]))
            {
                M_dfloe[i]  = std::max<T_val>(M_ice_params.M_dmin,lam/2.);
                M_broken[i] = 1.;
                M_nfloes[i] = this->dfloeToNfloes(M_dfloe[i],M_conc[i]);
            }
        }
    }
}//doBreaking




template<typename T>
typename IceInfo<T>::T_val
IceInfo<T>::dfloeToNfloes(T_val const& dfloe_in,
                           T_val const& conc_in)
{
    T_val nfloes_out   = 0.;

    if ( (dfloe_in>0) &&(conc_in >= M_ice_params.M_cice_min) )
    {
        //conc high enough & dfloe OK
        nfloes_out = conc_in/std::pow(dfloe_in,2.);
    }

    return nfloes_out;
}//dfloesToNfloes


template<typename T>
typename IceInfo<T>::T_val_vec
IceInfo<T>::dfloeToNfloes(T_val_vec const& dfloe_in,
                           T_val_vec const& conc_in)
{
    int N   = conc_in.size();
    T_val_vec nfloes_out(N);

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i=0;i<N;i++)
        nfloes_out[i]   = this->dfloeToNfloes(dfloe_in[i],conc_in[i]);

    return nfloes_out;
}//dfloesToNfloes


template<typename T>
typename IceInfo<T>::T_val
IceInfo<T>::nfloesToDfloe(T_val const& nfloes_in,
                           T_val const& conc_in)
{
        T_val dfloe_out  = 0.;
        if ( (nfloes_in>0) && (conc_in >= M_ice_params.M_cice_min) )
        {
            //conc high enough & Nfloes OK
            dfloe_out   = std::sqrt(conc_in/nfloes_in);
        }

        //dfloe shouldn't get too big
        dfloe_out   = std::min(M_ice_params.M_dfloe_pack_thresh,dfloe_out);

    return dfloe_out;
}//nfloesToDfloe


template<typename T>
typename IceInfo<T>::T_val_vec
IceInfo<T>::nfloesToDfloe(T_val_vec const& nfloes_in,
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
void IceInfo<T>::floeScaling(
      T_val const& dmax, int const& moment, T_val& dave)
{
    T_val nm,nm1,dm,nsum,ndsum,r;

    int mm;
    T_val ffac     = M_ice_params.M_fragility*std::pow(M_ice_params.M_xi,2);

    dave = std::max(std::pow(M_ice_params.M_dmin,moment),std::pow(dmax,moment));
    if (dmax>=M_ice_params.M_xi*M_ice_params.M_dmin)
    {
       //determine number of breaks
       r    = dmax/M_ice_params.M_dmin;
       mm   = 0;
       while (r >= M_ice_params.M_xi)
       {
          //r<2,mm=0 => doesn't break in 2
          //dave stays at dmax^moment;
          r  = r/M_ice_params.M_xi;
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
              nm     = nm1*(1-M_ice_params.M_fragility);
              nsum  += nm;
              ndsum += nm*std::pow(dm,moment);
              //std::cout<<"nsum,dm: "<<nsum<<" , "<<dm<<"\n";

              nm1   *= ffac;
              dm    /= M_ice_params.M_xi;
          }

          //m=mm:
          nsum   += nm1;
          ndsum  += nm1*std::pow(dm,moment);
          dave    = ndsum/nsum;
          //std::cout<<"nsum,dm: "<<nsum<<" , "<<dm<<"\n";
       }
    }
}//floeScaling

template<typename T>
void IceInfo<T>::floeScalingSmooth(
      T_val const& dmax,int const& moment, T_val& dave)
{

    T_val fsd_exp,b,A;

    fsd_exp = 2+log(M_ice_params.M_fragility)/log(M_ice_params.M_xi);//power law exponent: P(d>D)=(D_min/D)^fsd_exp;
    b       = moment-fsd_exp;

    // calculate <D^moment> from Dmax
    // - uniform dist for larger floes
    dave = std::pow(dmax,moment);
    //std::cout<<"dave (1) = "<<dave<<"\n";

    if (dmax<=M_ice_params.M_dmin)
    {
        // small floes
        dave = std::pow(M_ice_params.M_dmin,moment);
        //std::cout<<"dave (2) = "<<dave<<"\n";
    }
    else
    {
        // bigger floes
        A    = (fsd_exp*std::exp(fsd_exp*(std::log(M_ice_params.M_dmin)+std::log(dmax))));
        A    = A/(std::exp(fsd_exp*std::log(dmax))-std::exp(fsd_exp*std::log(M_ice_params.M_dmin)));
        dave = -(A/b)*(std::exp(b*std::log(M_ice_params.M_dmin))-exp(b*std::log(dmax)));
        //std::cout<<"dave (3) = "<<dave_<<"\n";
    }
}//floeScalingSmooth

template<typename T>
void IceInfo<T>::getDave(T_val_vec &dave,int const& moment)
{

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i=0;i<M_num_elements;i++)
    {
        if(M_mask[i]<.5)//not ice
            dave[i] = 0.;
        else//ice
            if (M_dfloe[i] <M_ice_params.M_dfloe_miz_thresh)
            {
                if ( M_ice_params.M_fsdopt == "RG" )
                    this->floeScaling(M_dfloe[i],moment,dave[i]);
                else if ( M_ice_params.M_fsdopt == "PowerLawSmooth" )
                    this->floeScalingSmooth(M_dfloe[i],moment,dave[i]);
                //std::cout<<"dave ("<<fsdopt<<") = "<<dave[i]<<"\n";
            }
            else
            {
                //just use uniform dist
                dave[i] = M_dfloe[i];
                //std::cout<<"dave (uniform) = "<<dave[i]<<"\n";
            }
    }
}

template<typename T>
void IceInfo<T>::printRange(std::string const &name,T_val_vec const &vec, int const & prec) const
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

// instantiate class for type float
//template class IceInfo<float>;

// instantiate class for type double
template class IceInfo<double>;

// instantiate class for type double
template class IceParams<double>;

} // namespace Wim
