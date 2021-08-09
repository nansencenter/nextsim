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

case `hostname | cut -d. -f2` in
    "ige-mcpc-36"   ) NXTSM_DEP_DIR="/opt/nextsim_intel"
                      INTEL_ROOT="/opt/intel/oneapi"
                      export INTEL_DIR="${INTEL_ROOT}/compiler/latest/linux"
                      export INTEL_MPI_DIR="${INTEL_ROOT}/mpi/latest"
                      export NETCDF_DIR=/opt/hdf5_netcdf4_intel
                      ;;
    "merlat"        ) NXTSM_DEP_DIR="/opt/nextsim_intel"
                      INTEL_ROOT="/opt/intel/oneapi"
                      export INTEL_DIR="${INTEL_ROOT}/compiler/latest/linux"
                      export INTEL_MPI_DIR="${INTEL_ROOT}/mpi/latest"
                      export NETCDF_DIR=/opt/hdf5_netcdf4_intel
                      #
                      ;;
    "ige-meom-cal1" ) NXTSM_DEP_DIR="/mnt/meom/workdir/brodeau/opt/nextsim_intel"
                      INTEL_ROOT="/mnt/meom/workdir/brodeau/opt/intel/oneapi"
                      export INTEL_DIR="${INTEL_ROOT}/compiler/latest/linux"
                      export INTEL_MPI_DIR="${INTEL_ROOT}/mpi/latest"
                      export NETCDF_DIR=/mnt/meom/workdir/brodeau/opt/hdf5_netcdf4_intel_par
                      #export NETCDF_DIR=/mnt/meom/workdir/brodeau/opt/hdf5_netcdf4_intel_seq
                      ;;
    "occigen" )       NXTSM_DEP_DIR="/store/CT1/ige2071/brodeau/opt/nextsim_intel"
                      INTEL_ROOT="/opt/software/common/intel"
                      export INTEL_DIR="${INTEL_ROOT}/compilers_and_libraries_2019.4.243/linux"
                      export INTEL_MPI_DIR="${INTEL_ROOT}/impi/2019.4.243/intel64"
                      export NETCDF_DIR="/store/CT1/hmg2840/lbrodeau/opt/hdf5_netcdf4_intel_mpi"
                      ;;
    #
    *               ) echo "Unknow HOSTNAME `hostname` !"
                      exit
                      ;;
esac



#################### Loading INTEL system #############################
export PATH=${INTEL_DIR}/bin/intel64:${PATH}
export LD_LIBRARY_PATH=${INTEL_DIR}/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
export CPATH=${INTEL_DIR}/include:${CPATH}
export LDFLAGS="-L${INTEL_DIR}/compiler/lib/intel64_lin -limf"
#
# Intel MPI:
export PATH=${INTEL_MPI_DIR}/bin:${PATH}
export LD_LIBRARY_PATH=${INTEL_MPI_DIR}/lib:${LD_LIBRARY_PATH}
export CPATH=${INTEL_MPI_DIR}/include:${CPATH}
#######################################################################

export MPI_DIR=${INTEL_MPI_DIR}

export LD_EXTRA_AEROBULK="-L${INTEL_DIR}/compiler/lib/intel64_lin -lifcore"


#######################################################################
export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort
#
export CFLAGS="-xHost -O3"
export CXXFLAGS="-xHost -O3"
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





