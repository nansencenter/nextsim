/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

#include <finiteelement.hpp>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

namespace Nextsim
{
    po::options_description descrOptions();
}

int main(int argc, char** argv )
{
    using namespace Nextsim;

    Environment env(argc, argv, descrOptions());

    if ( Environment::vm().count("help")  )
    {
        std::cout<< "BOOST VERSION= "<< BOOST_LIB_VERSION <<"\n";
        std::cout << "Basic Command Line Parameter Application" <<"\n"
                  << descrOptions() << "\n";
        return 0;
    }

    FiniteElement FE;
    FE.run();
}
