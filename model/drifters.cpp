/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   drifters.cpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Sun Feb 19 09:31:24 CET 2017
 */

#include <drifters.hpp>

/**
 * @class Drifters
 * @brief Manage drifters (base class)
 *
 * @see
 *
*/

namespace Nextsim
{
    /*
     * Constructors and destructors
     * We have one constructor for regularly spaced drifters and one for
     * drifter positions read in from file. In addition we have a private
     * function (maskXY), called only by the constructors, to mask out drifters
     * placed outside the ice cover.
     * There is also a default, constructor which initialises things to zero and false.
     */

	Drifters::Drifters()
    {
        M_is_initialised = false;
        M_no_drifters    = 0;
        M_X.resize(0);
        M_Y.resize(0);
    }

	Drifters::Drifters(double spacing, GmshMesh const &mesh, std::vector<double> &conc, double climit)
    {
        // Calculate the grid spacing assuming a regular grid
        std::vector<double> RX = mesh.coordX();
        std::vector<double> RY = mesh.coordY();
        auto xcoords = std::minmax_element( RX.begin(), RX.end() );
        auto ycoords = std::minmax_element( RY.begin(), RY.end() );

        int ncols = (int) ( 0.5 + ( *xcoords.second - *xcoords.first )/spacing );
        int nrows = (int) ( 0.5 + ( *ycoords.second - *ycoords.first )/spacing );

        std::vector<double> X(ncols*nrows);
        std::vector<double> Y(ncols*nrows);

        int i=0;
        int y = *ycoords.first;
        for (int row=0; row<nrows; ++row)
        {
            int x = *xcoords.first;
            for (int col=0; col<ncols; ++col)
            {
                X[i] = x;
                Y[i] = y;
                ++i;
                x += spacing;
            }
            y += spacing;
        }

        // Apply mask using conc and climit, and save to M_X and M_Y
        maskXY(mesh, X, Y, conc, climit);

        M_is_initialised = true;
    }

    Drifters::Drifters(std::string dirname, std::string gridFile, std::string dimNameX, std::string dimNameY, std::string latName, std::string lonName, GmshMesh const &mesh, std::vector<double> &conc, double climit)
    {
        // Load the grid from file
        // Check file
        std::string filename = (boost::format( "%1%/%2%/%3%" )
                                % Environment::simdataDir().string()
                                % dirname
                                % gridFile
                                ).str();
        if ( ! boost::filesystem::exists(filename) )
            throw std::runtime_error("Drifters::Drifters: File not found: " + filename);

        // Open file
        netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

        // Read the dimension of the grid
        netCDF::NcDim dim;
        dim = dataFile.getDim(dimNameX);
        int ncols = dim.getSize();

        dim = dataFile.getDim(dimNameY);
        int nrows = dim.getSize();

        int gridSize = ncols*nrows;

        // Read the lat and lon
        netCDF::NcVar VLAT = dataFile.getVar(latName);
        netCDF::NcVar VLON = dataFile.getVar(lonName);

        // Read data into LON & LAT
        std::vector<double> LAT(gridSize);
        std::vector<double> LON(gridSize);
        VLAT.getVar(&LAT[0]);
        VLON.getVar(&LON[0]);

        // Calculate x and y
        std::vector<double> X(gridSize);
        std::vector<double> Y(gridSize);

        mapx_class *map;
        filename = Environment::nextsimDir().string() + "/data/" + Environment::vm()["simul.proj_filename"].as<std::string>();
        std::vector<char> str(filename.begin(), filename.end());
        str.push_back('\0');

        map = init_mapx(&str[0]);

        for (int i=0; i<gridSize; ++i)
            forward_mapx(map, LAT[i], LON[i], &X[i], &Y[i]);

        close_mapx(map);

        // Apply mask using conc and climit, and save to M_X and M_Y
        maskXY(mesh, X, Y, conc, climit);

        M_is_initialised = true;
    }

    // Mask out the initial X and Y values so we only have drifters where there's ice
    void Drifters::maskXY(GmshMesh const &mesh, std::vector<double> &X, std::vector<double> &Y, std::vector<double> &conc, double clim)
    {
        // Interpolate the concentration onto the drifter positions
        int gridSize = X.size();
        double* interp_drifter_out;
        InterpFromMeshToMesh2dx(&interp_drifter_out,
                                &mesh.indexTr()[0],&mesh.coordX()[0],&mesh.coordY()[0],
                                mesh.numNodes(),mesh.numTriangles(),
                                &conc[0],
                                mesh.numTriangles(),1,
                                &X[0],&Y[0],gridSize,
                                true, 0.);

        // Add drifter positions where conc > conc_lim
        M_X.resize(0); // this shouldn't be necessary!
        M_Y.resize(0);
        for ( int i=0; i<gridSize; ++i )
        {
            if ( interp_drifter_out[i] > clim )
            {
                M_X.push_back(X[i]);
                M_Y.push_back(Y[i]);
            }
        }
        M_no_drifters = M_X.size();

        xDelete<double>(interp_drifter_out);
    }

    // Check to see if we're properly initialised
    bool Drifters::isInitialised()
    {
        return M_is_initialised;
    }

    /*
     * There are two use-cases for moving drifters. The "normal" case where you
     * overwrite the old coordinates with the new ones and the case where you
     * want to store the new coordinates in a seperate variable. This latter
     * case is useful when writing to file
     */

    // Move drifters puting output into X and Y
    void Drifters::move(GmshMesh const &mesh, std::vector<double> const &UT, std::vector<double> &X, std::vector<double> &Y)
    {
        // Do nothing if we don't have to
        if ( M_no_drifters == 0 )
            return;

        // Interpolate the total displacement onto the drifter positions
        int nb_var=2;
        int numNodes = mesh.numNodes();
        std::vector<double> interp_drifter_in(nb_var*numNodes);

        // Interpolate the velocity
        for (int i=0; i<numNodes; ++i)
        {
            interp_drifter_in[nb_var*i]   = UT[i];
            interp_drifter_in[nb_var*i+1] = UT[i+numNodes];
        }

        double* interp_drifter_out;
        InterpFromMeshToMesh2dx(&interp_drifter_out,
                                &mesh.indexTr()[0],&mesh.coordX()[0],&mesh.coordY()[0],
                                numNodes,mesh.numTriangles(),
                                &interp_drifter_in[0],
                                numNodes,nb_var,
                                &M_X[0],&M_Y[0],M_no_drifters,
                                true, 0.);

        for ( int i=0; i<M_no_drifters; ++i )
        {
            X[i] = M_X[i] + interp_drifter_out[i];
            Y[i] = M_Y[i] + interp_drifter_out[i+M_no_drifters];
        }

        xDelete<double>(interp_drifter_out);
    }

    // Move drifters and replace the old coordinates with the new ones
    void Drifters::move(GmshMesh const &mesh, std::vector<double> const &UT)
    {
        move(mesh, UT, M_X, M_Y);
    }

    /*
     * File operations
     *
     * We use netCDF files to save a fixed number of drifters. The file
     * contains two vectors of latitude and longitude positions and can be
     * appended to through time. File name contains time stamp for the creation
     * date (format: yyyymmdd, we assume we won't be creating multiple files
     * with the same prefix per day.
     */

    // Initialise the netCDF file
    void Drifters::initNetCDF(std::string file_prefix, double current_time)
    {
        // Construct the filename
        boost::gregorian::date now = Nextsim::parse_date(current_time);
        std::stringstream filename;

        filename << file_prefix;
        filename << now.year() << setw(2) << setfill('0') << now.month().as_number() << now.day();
        filename << ".nc";
        M_filename = filename.str();

        /* Not sure if I want to check if the file exists
        // Throw an error if the file exists
        if ( ! boost::filesystem::exists(filename) )
            throw std::runtime_error("Drifters::initNetCDF: Trying to initialise a file that already exists: " + filename);
        */

        // Create the netCDF file.
        netCDF::NcFile dataFile(filename.str(), netCDF::NcFile::replace);

        // Create the time dimension
        netCDF::NcDim tDim = dataFile.addDim("time"); // unlimited

        // Create the time variable
        netCDF::NcVar time = dataFile.addVar("time", netCDF::ncFloat, tDim);
        time.putAtt("standard_name","time");
        time.putAtt("long_name","simulation time");
        time.putAtt("units","days since 1900-01-01 00:00:00");
        time.putAtt("calendar","standard");

        M_nc_step=0;

        // Create the vector dimension
        netCDF::NcDim vecDim = dataFile.addDim("x", M_no_drifters);
        netCDF::NcDim emptyDim = dataFile.addDim("y", 1);
        std::vector<netCDF::NcDim> dims2(3);
        dims2[0] = tDim;
        dims2[1] = vecDim;
        dims2[2] = emptyDim;

        // Create the output variables
        // Longitude
        netCDF::NcVar lon = dataFile.addVar("longitude", netCDF::ncFloat, dims2);
        lon.putAtt("standard_name","longitude");
        lon.putAtt("long_name","longitude");
        lon.putAtt("units","degrees_north");

        // Latitude
        netCDF::NcVar lat = dataFile.addVar("latitude", netCDF::ncFloat, dims2);
        lat.putAtt("standard_name","latitude");
        lat.putAtt("long_name","latitude");
        lat.putAtt("units","degrees_north");

        dataFile.putAtt("Conventions", "CF-1.6");
        dataFile.putAtt("institution", "NERSC, Thormoehlens gate 47, N-5006 Bergen, Norway");
        dataFile.putAtt("source", "neXtSIM model fields");
    }

    // Write data to the netCDF file
    void Drifters::appendNetCDF(double current_time, GmshMesh const &mesh, std::vector<double> const &UT)
    {
        // Move the drifters before export, but save the result in a different variable
        std::vector<double> X(M_no_drifters);
        std::vector<double> Y(M_no_drifters);
        move(mesh, UT, X, Y);

        // Calculate lat and lon
        mapx_class *map;
        std::string mpp_file = Environment::nextsimDir().string() + "/data/" + Environment::vm()["simul.proj_filename"].as<std::string>();
        std::vector<char> str(mpp_file.begin(), mpp_file.end());
        str.push_back('\0');

        map = init_mapx(&str[0]);

        std::vector<double> lat(M_no_drifters);
        std::vector<double> lon(M_no_drifters);
        for (int i=0; i<M_no_drifters; ++i)
            inverse_mapx(map, X[i], Y[i], &lat[i], &lon[i]);

        close_mapx(map);

        // Open the netCDF file
        netCDF::NcFile dataFile(M_filename, netCDF::NcFile::write);

        // Append to time
        std::vector<size_t> start;
        start.push_back(M_nc_step);

        std::vector<size_t> count;
        count.push_back(1);

        M_nc_step++;

        netCDF::NcVar time = dataFile.getVar("time");
        time.putVar(start, count, &current_time);

        // Append to the output variables
        start.push_back(0);
        count.push_back(M_no_drifters);

        start.push_back(0);
        count.push_back(1);

        netCDF::NcVar data;

        data = dataFile.getVar("longitude");
        data.putVar(start, count, &lon[0]);

        data = dataFile.getVar("latitude");
        data.putVar(start, count, &lat[0]);
    }
}
