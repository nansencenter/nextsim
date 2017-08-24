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


template<typename T> IceInfo<T>::~IceInfo() {}
template<typename T> IceInfo<T>::IceInfo() {}

template<typename T> IceInfo<T>::IceInfo(po::variables_map const& vm)
{
    this->update(vm);
}

template<typename T>
void IceInfo<T>::update(po::variables_map const& vm)
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
void IceInfo<T>::init()
{
    M_max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
    M_sigma_c  = (1.76e+6)*std::exp(-5.88*std::sqrt(M_vbf));//flexural strength (Pa)
    M_epsc = M_sigma_c/M_young;//breaking strain
}//end ::init()

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
}//setFields()

template<typename T>
void IceInfo<T>::setFields(value_type_vec const& conc,
            value_type_vec const& vol,
            value_type_vec const& nfloes)
{
    // this needs to be called each time grid or mesh changes
    // ie initially, and if(M_wim_on_mesh), after regridding
    // set sizes of arrays, initialises some others that are constant in time

    M_conc          = conc;
    M_vol           = vol;
    M_nfloes        = nfloes;
    this->setFields();
}//setFields()


template<typename T>
void IceInfo<T>::setFields()
{
    // this needs to be called each time grid or mesh changes
    // ie initially, and if(M_wim_on_mesh), after regridding
    // set sizes of arrays, initialises some others that are constant in time

    M_num_elements  = M_conc.size();

    //2D var's
    M_dfloe.assign(M_num_elements,0.);
    M_thick.assign(M_num_elements,0.);
    M_broken.assign(M_num_elements,0.);
    M_mask.assign(M_num_elements,1.);

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i=0;i<M_num_elements;i++)
    {
        if (M_conc[i]>=M_cice_min)
            M_thick[i]  = M_vol[i]/M_conc[i];//convert to actual thickness
        else
        {
            M_conc[i]   = 0;
            M_vol[i]    = 0;
            M_nfloes[i] = 0;
            M_mask[i]   = 0.;
        }
        M_dfloe[i]  = this->nfloesToDfloe(M_nfloes[i],M_conc[i]);

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
    }
}//setFields()

template<typename T>
void IceInfo<T>::doBreaking(value_type_vec const& mom0,value_type_vec const& mom2,value_type_vec const& var_strain)
{

    std::vector<double> params = this->getAttenParams();//currently independant of space

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i=0;i<M_num_elements;i++)
    {
        bool break_criterion = (M_mask[i] > 0.5)          //ice present
            && ((2*var_strain[i]) > std::pow(M_epsc, 2.));//big enough waves

        if (break_criterion)
        {
            double outputs[8];

            value_type om    = std::sqrt(mom2[i]/mom0[i]);
            value_type guess = std::pow(om,2.)/M_gravity;

            RTparam_outer(outputs,M_thick[i],double(om),double(M_drag_rp),double(guess),&params[0]);

            value_type kice = outputs[1];
            value_type lam = 2*PI/kice;

            if (lam < (2*M_dfloe[i]))
            {
                M_dfloe[i]  = std::max<value_type>(M_dmin,lam/2.);
                M_broken[i] = 1.;
                M_nfloes[i] = this->dfloeToNfloes(M_dfloe[i],M_conc[i]);
            }
        }
    }
}//doBreaking


template<typename T>
std::vector<double> IceInfo<T>::getAttenParams()
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
typename IceInfo<T>::value_type
IceInfo<T>::dfloeToNfloes(value_type const& dfloe_in,
                           value_type const& conc_in)
{
    value_type nfloes_out   = 0.;

    if ( (dfloe_in>0) &&(conc_in >= M_cice_min) )
    {
        //conc high enough & dfloe OK
        nfloes_out = conc_in/std::pow(dfloe_in,2.);
    }

    return nfloes_out;
}//dfloesToNfloes


template<typename T>
typename IceInfo<T>::value_type_vec
IceInfo<T>::dfloeToNfloes(value_type_vec const& dfloe_in,
                           value_type_vec const& conc_in)
{
    int N   = conc_in.size();
    value_type_vec nfloes_out(N);

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i=0;i<N;i++)
        nfloes_out[i]   = this->dfloeToNfloes(dfloe_in[i],conc_in[i]);

    return nfloes_out;
}//dfloesToNfloes


template<typename T>
typename IceInfo<T>::value_type
IceInfo<T>::nfloesToDfloe(value_type const& nfloes_in,
                           value_type const& conc_in)
{
        value_type dfloe_out  = 0.;
        if ( (nfloes_in>0) && (conc_in >= M_cice_min) )
        {
            //conc high enough & Nfloes OK
            dfloe_out   = std::sqrt(conc_in/nfloes_in);
        }

        //dfloe shouldn't get too big
        dfloe_out   = std::min(M_dfloe_pack_thresh,dfloe_out);

    return dfloe_out;
}//nfloesToDfloe


template<typename T>
typename IceInfo<T>::value_type_vec
IceInfo<T>::nfloesToDfloe(value_type_vec const& nfloes_in,
                           value_type_vec const& conc_in)
{
    int N   = conc_in.size();
    value_type_vec dfloe_out(N);

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i=0;i<N;i++)
        dfloe_out[i] = this->nfloesToDfloe(nfloes_in[i],conc_in[i]);

    return dfloe_out;
}//nfloesToDfloe


template<typename T>
void IceInfo<T>::floeScaling(
      value_type const& dmax, int const& moment, value_type& dave)
{
    value_type nm,nm1,dm,nsum,ndsum,r;

    int mm;
    value_type ffac     = M_fragility*std::pow(M_xi,2);

    dave = std::max(std::pow(M_dmin,moment),std::pow(dmax,moment));
    if (dmax>=M_xi*M_dmin)
    {
       //determine number of breaks
       r    = dmax/M_dmin;
       mm   = 0;
       while (r >= M_xi)
       {
          //r<2,mm=0 => doesn't break in 2
          //dave stays at dmax^moment;
          r  = r/M_xi;
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
              nm     = nm1*(1-M_fragility);
              nsum  += nm;
              ndsum += nm*std::pow(dm,moment);
              //std::cout<<"nsum,dm: "<<nsum<<" , "<<dm<<"\n";

              nm1   *= ffac;
              dm    /= M_xi;
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
      value_type const& dmax,int const& moment, value_type& dave)
{

    value_type fsd_exp,b,A;

    fsd_exp = 2+log(M_fragility)/log(M_xi);//power law exponent: P(d>D)=(D_min/D)^fsd_exp;
    b       = moment-fsd_exp;

    // calculate <D^moment> from Dmax
    // - uniform dist for larger floes
    dave = std::pow(dmax,moment);
    //std::cout<<"dave (1) = "<<dave<<"\n";

    if (dmax<=M_dmin)
    {
        // small floes
        dave = std::pow(M_dmin,moment);
        //std::cout<<"dave (2) = "<<dave<<"\n";
    }
    else
    {
        // bigger floes
        A    = (fsd_exp*std::exp(fsd_exp*(std::log(M_dmin)+std::log(dmax))));
        A    = A/(std::exp(fsd_exp*std::log(dmax))-std::exp(fsd_exp*std::log(M_dmin)));
        dave = -(A/b)*(std::exp(b*std::log(M_dmin))-exp(b*std::log(dmax)));
        //std::cout<<"dave (3) = "<<dave_<<"\n";
    }
}//floeScalingSmooth

template<typename T>
void IceInfo<T>::getDave(value_type_vec &dave,int const& moment)
{

#pragma omp parallel for num_threads(M_max_threads) collapse(1)
    for (int i=0;i<M_num_elements;i++)
    {
        if(M_mask[i]<.5)//not ice
            dave[i] = 0.;
        else//ice
            if (M_dfloe[i] <M_dfloe_miz_thresh)
            {
                if ( M_fsdopt == "RG" )
                    this->floeScaling(M_dfloe[i],moment,dave[i]);
                else if ( M_fsdopt == "PowerLawSmooth" )
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

// instantiate class for type float
//template class IceInfo<float>;

// instantiate class for type double
template class IceInfo<double>;


} // namespace Wim
