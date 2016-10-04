/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#include <laplacian.hpp>
#include <netcdf>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

using namespace netCDF;
using namespace netCDF::exceptions;

// We are reading 2D data, a 6 x 12 grid.
static const int NX = 6;
static const int NY = 12;

// Return this in event of a problem.
static const int NC_ERR = 2;

namespace Nextsim
{
    po::options_description descrOptions();
}

int main(int argc, char** argv )
{
    using namespace Nextsim;

    Environment env(argc, argv, descrOptions());

    //namespace po = boost::program_options;
    //po::options_description desc = descrOptions();
    //po::variables_map vm;

    if ( Environment::vm().count("help")  )
    {
        std::cout<< "BOOST VERSION= "<< BOOST_LIB_VERSION <<"\n";
        std::cout << "Basic Command Line Parameter Application" <<"\n"
                  << descrOptions() << "\n";
        //return SUCCESS;
        return 0;
    }

    //std::cout<<"nx= "<< vm["nx"]. as<int>() <<"\n";

    try
    {
        // This is the array we will read.
        int dataIn[NX][NY];
        std::string filenc = (boost::format( "%1%/data/simple_xy.nc" ) % Environment::nextsimDir().string()).str();
        //std::cout<<"filenc= "<< filenc <<"\n";

        // Open the file for read access
        //NcFile dataFile("../data/simple_xy.nc", NcFile::read);
        NcFile dataFile(filenc, NcFile::read);

        // Retrieve the variable named "data"
        NcVar data=dataFile.getVar("data");
        if(data.isNull()) return NC_ERR;
        data.getVar(dataIn);

        // Check the values.
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
            {
                //std::cout<<"simple_xy.nc= "<< dataIn[i][j] <<"\n";
                if (dataIn[i][j] != i * NY + j)
                    return NC_ERR;
            }

        // The netCDF file is automatically closed by the NcFile destructor
        //cout << "*** SUCCESS reading example file simple_xy.nc!" << endl;
    }
    catch(NcException& e)
    {
        e.what();
        cout<<"FAILURE*************************************"<<endl;
        return NC_ERR;
    }

    Laplacian LP;
    LP.run();
}
