/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

#include <wimdiscr.hpp>
#include <iomanip>
#include <options_wim.hpp>

namespace
{
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
} // namespace

int main(int argc, char** argv )
{
    using namespace Wim;

    po::options_description desc = descrWimOptions();
    po::variables_map vm;

    try
    {
        po::store(po::parse_command_line(argc, argv, desc),vm);

        if ( vm.count("help")  )
        {
            std::cout << "Basic Command Line Parameter Application" <<"\n"
                      << desc << "\n";
            return SUCCESS;
        }

        if ( vm.count( "config-file" ) )
        {
            if ( fs::exists( vm["config-file"].as<std::string>() ) )
            {
                std::ifstream ifs( vm["config-file"].as<std::string>().c_str() );
                po::store( parse_config_file( ifs, desc, true ), vm );
                po::notify( vm );
            }
            else
            {
                std::cout << "Cannot find " << "config-file" << "\n";
            }
        }

        po::notify(vm);

    }
    catch(po::error& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl <<"\n";
        std::cerr << desc <<"\n";
        return ERROR_IN_COMMAND_LINE;
    }

    // instantiation & initialization of wim2d
    WimDiscr<double> wim(vm);
    //WimDiscr<float> wim(vm);

    wim.initStandAlone();

#if 0
    //test fsd:
    double dave,dmax;
    std::vector<double> dmaxvec={15,20,30,40,80,160,200,220};
    int mom=2;

    for (int j=0; j<dmaxvec.size(); ++j)
    {
       dmax = dmaxvec[j];
       wim.floeScalingSmooth(dmax,mom,dave);
       //wim.floeScaling(dmax,2,dave);
       std::cout<<std::endl
                <<"<D> (m): "<<std::setprecision(10)<<dave
                <<std::endl;
    }
    std::abort();
#endif

    // run the simulation
    wim.run();
      /* NB calls "assign" subroutine
       * - some quantities need reassignment at each call to wim
       *   eg conc, thickness, nfloes/dmax
       */

#if 0
    //std::vector<float> vec1 = {1.,2.,3.,4.,5.};
    std::vector<std::complex<float> > vec1 = {(1.,2.),2.,3.,4.,5.};
    std::vector<std::complex<float> > vec2 = {(1.5,1.),4.,7.,1.,2.};
    // float sum_of_elems = std::accumulate(vec1.begin(),vec1.end(),0.);
    //std::transform(vec1.begin(), vec1.end(), vec2.begin(), vec1.begin(), std::multiplies<std::complex<float> >());
    //std::complex<float> ans = (1.,1.);
    //ans = std::accumulate(vec1.begin(), vec1.end(), std::complex<float>{1.0}, std::multiplies<std::complex<float> >{});

    std::transform(vec1.begin(), vec1.end(), vec1.begin(), vec1.begin(), [](std::complex<float> a, std::complex<float> b) -> std::complex<float> { return a+b; });

    //std::cout<<"ANS= "<< ans <<"\n";

    for (auto const& it : vec1)
    {
        std::cout<<"Real= "<< std::real(it) <<"\n";
        std::cout<<"Imag= "<< it.imag() <<"\n";
    }
    //float sum_of_elems = std::accumulate(vec1.begin(),vec1.end(),0.);

    //std::cout<<"sum= "<< sum_of_elems <<"\n";

    std::complex<float> a; //(1,2);
    std::complex<float> b; //(3,5);

    a.real(1);
    a.imag(2);
    b = a;

    std::cout<<"A*B= "<< b <<"\n";

    //using namespace std::literals;
    //std::complex<double> z1 = 1i * 1i;
#endif

}
