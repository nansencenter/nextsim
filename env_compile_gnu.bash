#!/bin/bash

########################################################################
# Example of bash script to set the compilation environment
# now that all Makefiles are rid of the hard-coded `CXX = mpicxx`
# Also use of following variables:
#  MPI_DIR
#  MPI_INC_DIR
#  MPI_LIB_DIR
#
# L. Brodeau, August 2021
########################################################################

GMSH_VERSION="3.0.6"
BOOST_VERSION="1.67"

export NEXTSIMDIR=`pwd`

#export USE_AEROBULK=true
#export AEROBULK_DIR=${HOME}/DEV/aerobulk
unset USE_AEROBULK
unset AEROBULK_DIR

# Defaults before arch-scanning:
export GNU_COMP_DIR="/usr"
export GNU_MPI_DIR="/usr"
export NETCDF_DIR="/usr"
export GNU_COMP_VERSION="9"

case `hostname | cut -d. -f2` in
    "merlat"        ) NXTSM_DEP_DIR="/opt/nextsim_gnu"
                      ;;
    "ige-meom-cal1" ) NXTSM_DEP_DIR="/mnt/meom/workdir/brodeau/opt/nextsim_gnu"
                      ;;
    "occigen" )       NXTSM_DEP_DIR="/store/CT1/ige2071/brodeau/opt/nextsim_gnu"
                      ;;
    *               ) echo "Unknow HOSTNAME `hostname` !"
                      exit
                      ;;
esac

# Get appropriate Gnu compiler include files:
export GNU_COMP_INC_DIR=${GNU_COMP_DIR}/include/c++/${GNU_COMP_VERSION}
export CPATH=${GNU_COMP_INC_DIR}:${CPATH}

#######################################################################
export CC=mpicc
export CXX=mpicxx
export FC=mpifort
#
export CFLAGS="-O2"
export CXXFLAGS="-O2 -pthread -I${GNU_COMP_INC_DIR}"
export LDFLAGS=""

export MPI_DIR=${GNU_MPI_DIR}

export LD_EXTRA_AEROBULK="-lgfortran" ; # if relevant...
#######################################################################

export GMSH_DIR=${NXTSM_DEP_DIR}/gmsh-${GMSH_VERSION}

export BOOST_DIR=${NXTSM_DEP_DIR}/boost-${BOOST_VERSION}
export BOOST_INCDIR=${BOOST_DIR}/include
export BOOST_LIBDIR=${BOOST_DIR}/lib

export NETCDF_CXX_DIR=${NETCDF_DIR}

export MPI_LIB_DIR=${MPI_DIR}/lib
export MPI_INC_DIR=${MPI_DIR}/include


#if [ ! -f ./model/version.hpp ]; then
cd model
./version.sh
cd ../
#fi





