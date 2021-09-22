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

#export USE_AEROBULK=true
#export AEROBULK_DIR="${HOME}/DEV/aerobulk"

############################################
# Defaults before host-specific adjustment #
############################################
### Linux / Intel / IntelMPI specific part:
export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort
#
export CFLAGS="-xHost -O3 -fPIC";                      # NOTE: many other flags are still hard-coded in the Makefiles!
export CCFLAGS="${CFLAGS}";                     # only for `mapx` ???
export CXXFLAGS="-xHost -O3 -pthread -fPIC -qopenmp";  # NOTE: many other flags are still hard-coded in the Makefiles!
export LDFLAGS=""
#
export INTEL_ROOT="/opt/intel/oneapi";                        # root directory of the Intel suite (compiler and MPI)
export INTEL_COMP_DIR="${INTEL_ROOT}/compiler/latest/linux";  #     "              " compiler
#
export MPI_DIR="${INTEL_ROOT}/mpi/latest"
export MPI_LIB_DIR=${MPI_DIR}/lib
export MPI_INC_DIR=${MPI_DIR}/include
#
export NETCDF_DIR="/opt/hdf5_netcdf4_intel"
export NETCDF_CXX_DIR=${NETCDF_DIR}
#
export NXTSM_DEP_DIR="/opt/nextsim_intel" ; # path to directory containing compiled BOOST and GMSH (with the relevant compiler!)
#############################################

######################################################################
# Host-specific adjustment of previously defined variables if needed #
######################################################################
case `hostname | cut -d. -f2` in
    "merlat" )        export LDFLAGS="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -limf"
                      ;;
    "ige-meom-cal1" ) export NXTSM_DEP_DIR="/mnt/meom/workdir/brodeau/opt/nextsim_intel"
                      export INTEL_ROOT="/mnt/meom/workdir/brodeau/opt/intel/oneapi"
                      export INTEL_COMP_DIR="${INTEL_ROOT}/compiler/latest/linux"
                      export NETCDF_DIR=/mnt/meom/workdir/brodeau/opt/hdf5_netcdf4_intel_par
                      ;;
    "occigen" )       export NXTSM_DEP_DIR="/store/CT1/ige2071/brodeau/opt/nextsim_intel"
                      export INTEL_ROOT="/opt/software/common/intel"
                      export INTEL_COMP_DIR="${INTEL_ROOT}/compilers_and_libraries_2019.4.243/linux"
                      export MPI_DIR="${INTEL_ROOT}/impi/2019.4.243/intel64"
                      export MPI_LIB_DIR=${MPI_DIR}/lib
                      export MPI_INC_DIR=${MPI_DIR}/include
                      export NETCDF_DIR="/store/CT1/hmg2840/lbrodeau/opt/hdf5_netcdf4_intel_mpi"
                      ;;
    "fram" )          export NXTSM_DEP_DIR="/cluster/projects/nn9878k/brodeau/opt/nextsim_intel"
                      export INTEL_ROOT="/cluster/software/ifort/2018.1.163-GCC-6.4.0-2.28"
                      export INTEL_COMP_DIR="${INTEL_ROOT}/compilers_and_libraries_2018.1.163/linux"
                      export MPI_DIR="/cluster/software/impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28/intel64"
                      export MPI_LIB_DIR=${MPI_DIR}/lib
                      export MPI_INC_DIR=${MPI_DIR}/include
                      export NETCDF_DIR="/cluster/software/netCDF/4.4.1.1-intel-2018a-HDF5-1.8.19"
                      export NETCDF_CXX_DIR="/cluster/software/netCDF-C++4/4.3.0-intel-2018a-HDF5-1.8.19"
                      ;;

    *               ) echo;
                      echo "WARNING: Unknow machine with HOSTNAME = `hostname` !"
                      echo
                      sleep 2
                      ;;
esac

export GMSH_DIR=${NXTSM_DEP_DIR}/gmsh-${GMSH_VERSION}

export BOOST_DIR=${NXTSM_DEP_DIR}/boost-${BOOST_VERSION}
#export BOOST_INCDIR=${BOOST_DIR}/include/boost
export BOOST_INCDIR=${BOOST_DIR}/include
export BOOST_LIBDIR=${BOOST_DIR}/lib

#export LD_EXTRA_AEROBULK="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -lifcore"

#if [ ! -f ./model/version.hpp ]; then
cd model
./version.sh
cd ../
#fi
