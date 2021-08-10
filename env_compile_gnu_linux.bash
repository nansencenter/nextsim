#!/bin/bash

########################################################################
# Example of bash script to set the compilation environment
# now that all Makefiles are rid of the hard-coded `CXX = mpicxx`
# Also use of following variables:
#  MPI_DIR
#  MPI_INC_DIR
#  MPI_LIB_DIR
#
# => adapt this file tou your system and source it prior to using
#    the `make` command to clean or compile neXtSIM
#
# L. Brodeau, August 2021
########################################################################

GMSH_VERSION="3.0.6"
BOOST_VERSION="1.67"

export NEXTSIMDIR=`pwd`

############################################
# Defaults before host-specific adjustment #
############################################
### Linux / GNU / OpenMPI specific part:
export CC=mpicc.openmpi
export CXX=mpicxx.openmpi
export FC=mpifort.openmpi
#
export CFLAGS="-O0 -fPIC";                      # NOTE: many other flags are still hard-coded in the Makefiles!
export CCFLAGS="${CFLAGS}";                     # only for `mapx` ???
export CXXFLAGS="-O0 -pthread -fPIC -fopenmp "  # NOTE: many other flags are still hard-coded in the Makefiles!
export LDFLAGS=""
#
export MPI_DIR="/usr" ; # default OpenMPI
export MPI_LIB_DIR=${MPI_DIR}/lib
export MPI_INC_DIR=${MPI_DIR}/include
#
export NETCDF_DIR="/usr"
export NETCDF_CXX_DIR=${NETCDF_DIR}
#
export NXTSM_DEP_DIR="/opt/nextsim_gnu" ; # path to directory containing compiled BOOST and GMSH (with the relevant compiler!)
#############################################

######################################################################
# Host-specific adjustment of previously defined variables if needed #
######################################################################
case `hostname | cut -d. -f2` in
    "ige-meom-cal1" ) NXTSM_DEP_DIR="/mnt/meom/workdir/brodeau/opt/nextsim_gnu"
                      ;;
    "occigen" )       NXTSM_DEP_DIR="/store/CT1/ige2071/brodeau/opt/nextsim_gnu"
                      ;;
    *               ) echo;
                      echo "WARNING: Unknow machine with HOSTNAME = `hostname` !"
                      echo
                      sleep 2
                      ;;
esac

export GMSH_DIR=${NXTSM_DEP_DIR}/gmsh-${GMSH_VERSION}

export BOOST_DIR=${NXTSM_DEP_DIR}/boost-${BOOST_VERSION}
export BOOST_INCDIR=${BOOST_DIR}/include/boost
export BOOST_LIBDIR=${BOOST_DIR}/lib

#if [ ! -f ./model/version.hpp ]; then
cd model
./version.sh
cd ../
#fi





