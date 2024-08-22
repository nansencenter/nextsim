#!/bin/bash

########################################################################
# Example of bash script to set the compilation environment
# now that all Makefiles are rid of the hard-coded `CXX = mpicxx`
# Also use of following variables:
#  MPI_DIR
#  MPI_INC_DIR
#  MPI_LIB_DIR
#
# => adapt this file to your system and source it prior to using
#    the `make` command to clean or compile neXtSIM
#
# L. Brodeau, October 2021
########################################################################

# DO NOT EXPORT "INTEL_ROOT" !

GMSH_VERSION="3.0.6"
BOOST_VERSION="1.67"

export NEXTSIMDIR=`pwd`

# AeroBulk turbulent air-sea flux computation
#--------------------------------------------
l_aerobulk=true ; # Call AeroBulk to compute air-sea fluxes ? If true, give the appropriate value for "AEROBULK_DIR" in the arch CASE block...
if ${l_aerobulk}; then
    export USE_AEROBULK=true
    # => then, later before launching neXtSIM, pick an algorithm by setting `ocean_bulk_formula=<algo>` of the `thermo` block in the config file
    #   ==> algos are: 'ecmwf', 'coare3p0', 'coare3p0', 'ncar' or 'andreas' ('ecmwf recomended if using an ECMWF-based atmo forcing)
    #   ==> use 'nextsim' if you want to use the old neXtSIM bulk formulae
else
    unset USE_AEROBULK ; # that's the whole point of this if/else/fi block, their could be "remnant" values in USE_AEROBULK...
fi

# Coupling with another GCM component via OASIS (WW3, NEMO, etc)
# --------------------------------------------------------------
l_cpl_oasis=false ; # neXtSIM is going to be coupled to something via OASIS ? If true, give the appropriate value for "OASIS_DIR" in the arch CASE block...
if ${l_cpl_oasis}; then
    export USE_OASIS=true
else
    unset USE_OASIS ; # that's the whole point of this if/else/fi block, their could be "remnant" values in USE_OASIS...
fi

# Remeshing with MMG
l_mmg=true ; # The remeshing part can be achieved by bamg but anisotropic remeshing can only be achieved with mmg2d
if ${l_mmg}; then
    export USE_MMG=true
else
    unset USE_MMG ; # that's the whole point of this if/else/fi block, their could be "remnant" values in USE_MMG...
fi

############################################
# Defaults before host-specific adjustment #
############################################
### Linux / GNU / OpenMPI specific part:
export CC=mpicc.openmpi
export CXX=mpicxx.openmpi
export FC=mpifort.openmpi
#
export CFLAGS="-O3 -fPIC";                      # NOTE: many other flags are still hard-coded in the Makefiles!
export CCFLAGS="${CFLAGS}";                     # only for `mapx` ???
export CXXFLAGS="-O3 -pthread -fPIC -fopenmp "  # NOTE: many other flags are still hard-coded in the Makefiles!
export LDFLAGS=""
export FFLAGS="-O2 -fopenmp -lstdc++ -fPIC"
#
export MPI_DIR="/usr" ; # default OpenMPI
#
export NETCDF_DIR="/usr"
export NETCDF_CXX_DIR="/usr"
#
export AEROBULK_DIR="${HOME}/DEV/aerobulk"
export OASIS_DIR="${HOME}/src/oasis3-mct"
#
NXTSM_DEP_DIR="/opt/nextsim_gnu" ; # path to directory containing compiled BOOST and GMSH (with the relevant compiler!)
#
#############################################

######################################################################
# Host-specific adjustment of previously defined variables if needed #
######################################################################
case `hostname | cut -d. -f2` in
    "ige-meom-cal1" ) NXTSM_DEP_DIR="/mnt/meom/workdir/brodeau/opt/nextsim_gnu"
                      export AEROBULK_DIR="${NXTSM_DEP_DIR}/aerobulk"
                      ;;
    "occigen" )       NXTSM_DEP_DIR="/store/CT1/ige2071/brodeau/opt/nextsim_gnu"
                      ;;
    "fram" )          echo "Someone fixes me!!!! Env. variables for 'fram' and Gnu compilers...." ; exit
                      ;;
    *               ) echo;
                      echo "WARNING: Unknow machine with HOSTNAME = `hostname` !"
                      echo
                      sleep 2
                      ;;
esac

if ${l_aerobulk} || ${l_cpl_oasis}; then
    export LDFLAGS+="-lgfortran"
fi

# Normally the following 2 are pretty standard:
export MPI_LIB_DIR=${MPI_DIR}/lib
export MPI_INC_DIR=${MPI_DIR}/include

# Third-party software dependencies, compiled with relevant compiler!
export GMSH_DIR=${NXTSM_DEP_DIR}/gmsh-${GMSH_VERSION}

export PARMMG_DIR=${NEXTSIMDIR}/deps/parmmg2d/build

export BOOST_DIR=${NXTSM_DEP_DIR}/boost-${BOOST_VERSION}
export BOOST_INCDIR=${BOOST_DIR}/include
export BOOST_LIBDIR=${BOOST_DIR}/lib

if [ ! -f ./model/version.hpp ]; then
    cd model
    ./version.sh
    cd ../
fi
