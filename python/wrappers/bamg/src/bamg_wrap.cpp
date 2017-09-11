/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   bamgwrap.cpp
 * @author Timothy Williams <Timothy.Williams@nersc.no>
 * @date   Mon Sep 11 2017
 */

#include <bamg_wrap.hpp>
#include <InterpFromMeshToMesh2dx.h>
#include <MemOps.h>
#include <iostream>

namespace PyWrap
{

std::vector<std::vector<double>> interpMeshToPointsCpp(
      std::vector<int> index,
      std::vector<double> xnods,
      std::vector<double> ynods,
      std::vector<std::vector<double>> data,
      std::vector<double> xout,
      std::vector<double> yout,
      bool isdefault, double defaultvalue)
{

    int Ne       = index.size()/3; // num of mesh elements
    int Nn       = xnods.size();   // num of mesh nodes
    int Nvars    = data.size();    // num of vars to interp
    int Npts_in  = data[0].size(); // number of spatial locations that input data are defined on (to determine if they are on nodes or elements)
    int Npts_out = xout.size();    // number of spatial locations that output data are defined on
 
    std::vector<double> interp_in(Npts_out*Nvars,0.);
    for (int i=0;i<Npts_in;i++)
        for (int p=0;p<Nvars;p++)
            interp_in[Nvars*i+p] = data[p][i];
 
    double *interp_out;
    InterpFromMeshToMesh2dx(&interp_out,&index[0],&xnods[0],&ynods[0],Nn,Ne,
                             &interp_in[0],Npts_in,Nvars,&xout[0],&yout[0],Npts_out,
                             isdefault,defaultvalue);
 
    std::vector<double> tmp(Npts_out,0.);
    std::vector<std::vector<double>> output(Nvars,tmp);
    for (int i=0;i<Npts_out;i++)
        for (int p=0;p<Nvars;p++)
            output[p][i] = interp_out[Nvars*i+p];
 
    xDelete<double>(interp_out);
 
    return output;
         
}//interpMeshToPoints()

}//namespace PyWrappers
