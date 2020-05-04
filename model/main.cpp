/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   main.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug 24 11:05:03 2015
 */

#include <finiteelement.hpp>
using namespace rheolef;
//using namespace std;
#include "torus.icc"
// namespace fs = boost::filesystem;
// namespace po = boost::program_options;

// namespace Nextsim
// {
//     po::options_description descrOptions();
// }

// ----------------------------------------------------------------------------------------
//! Main.cpp: launches the model run by calling the FiniteElement::run() function
int main(int argc, char** argv )
{
    // using namespace Nextsim;

    // Environment env(argc, argv, descrOptions());

    // if ( Environment::vm().count("help")  )
    // {
    //     std::cout<< "BOOST VERSION= "<< BOOST_LIB_VERSION <<"\n";
    //     std::cout << "Basic Command Line Parameter Application" <<"\n"
    //               << descrOptions() << "\n";
    //     return 0;
    // }


    environment rheolef (argc, argv);

    // geo lambda (argv[1]);
    geo lambda ("hypercube.geo");
    // size_t d = lambda.dimension();
    // space Xh (lambda, "P1");
    // field phi_h = interpolate(Xh, phi);
    // band gamma_h (phi_h);
    // field phi_h_band = phi_h [gamma_h.band()];
    // space Bh (gamma_h.band(), "P1");
    // Bh.block ("isolated");
    // Bh.unblock ("zero");
    // trial u (Bh); test v (Bh);
    // form  a  = integrate (gamma_h, dot(grad_s(u),grad_s(v)));
    // field c  = integrate (gamma_h, v);
    // field lh = integrate (gamma_h, f(d)*v);

    // vector<field> b (gamma_h.n_connected_component());
    // vector<Float> z (gamma_h.n_connected_component(), 0);
    // for (size_t i = 0; i < b.size(); i++) {
    //     const domain& cci = gamma_h.band() ["cc"+itos(i)];
    //     field phi_h_cci (Bh, 0);
    //     phi_h_cci [cci] = phi_h_band [cci];
    //     b[i] = phi_h_cci;
    // }
    // form  A = {{    a,  trans(b), c },
    //            {    b,     0,     0 },
    //            { trans(c), 0,     0 }};
    // field Fh = {    lh,    z,     0 };
    // field Uh (Fh.get_space(), 0);

    // A.set_symmetry(true);
    // problem pa (A);
    // pa.solve (Fh, Uh);
    // dout << catchmark("phi") << phi_h
    //      << catchmark("u")   << Uh[0];

    //FiniteElement FE;
    //FE.run();
}//main
