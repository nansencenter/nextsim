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

# DO NOT EXPORT "INTEL_ROOT" !

GMSH_VERSION="3.0.6"
BOOST_VERSION="1.67"

export NEXTSIMDIR=`pwd`

# Uncomment these 2 if you want AeroBulk to compute air-sea fluxes:
export USE_AEROBULK=true
export AEROBULK_DIR="${HOME}/DEV/aerobulk"

# This neXtSIM is going to be coupled to NEMO via OASIS so:
export USE_OASIS=true
export OASIS_DIR="${HOME}/src/oasis3-mct"

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
export FFLAGS="-O2 -qopenmp -lstdc++ -fPIC"
#
INTEL_ROOT="/opt/intel/oneapi";                        # root directory of the Intel suite (compiler and MPI)
export INTEL_COMP_DIR="${INTEL_ROOT}/compiler/latest/linux";  #     "              " compiler
#
export MPI_DIR="${INTEL_ROOT}/mpi/latest"
#
#
NXTSM_DEP_DIR="/opt/nextsim_intel" ; # path to directory containing compiled BOOST and GMSH (with the relevant compiler!)
#############################################

######################################################################
# Host-specific adjustment of previously defined variables if needed #
######################################################################
case `hostname | cut -d. -f2` in
    "merlat" )        export LDFLAGS="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -limf"
                      ;;
    "ige-mcpc-36"   ) NXTSM_DEP_DIR="/opt/nextsim_intel"
                      INTEL_ROOT="/opt/intel/oneapi"
                      export INTEL_COMP_DIR="${INTEL_ROOT}/compiler/latest/linux"
                      export MPI_DIR="${INTEL_ROOT}/mpi/latest"
                      export NETCDF_DIR=/opt/hdf5_netcdf4_intel
                      export NETCDF_CXX_DIR=${NETCDF_DIR}
                      ;;
    "ige-meom-cal1" ) NXTSM_DEP_DIR="/mnt/meom/workdir/brodeau/opt/nextsim_intel"
                      INTEL_ROOT="/mnt/meom/workdir/brodeau/opt/intel/oneapi"
                      export INTEL_COMP_DIR="${INTEL_ROOT}/compiler/latest/linux"
                      export MPI_DIR="${INTEL_ROOT}/mpi/latest"
                      export NETCDF_DIR=/mnt/meom/workdir/brodeau/opt/hdf5_netcdf4_intel_par
                      export NETCDF_CXX_DIR=${NETCDF_DIR}
                      #
                      export LD_EXTRA_OASIS="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -lifcore -lifport"
                      #export FFLAGS="-O2 -qopenmp -lstdc++ -fPIC -I`which ifort | sed -e "s|bin/intel64/ifort|compiler/include/intel64|g"`" ; # only for C++ to find module "iso_c_binding":
                      ;;
    "occigen" )       NXTSM_DEP_DIR="/store/CT1/hmg2840/lbrodeau/opt/nextsim_intel"
                      INTEL_ROOT="/opt/software/common/intel"
                      export INTEL_COMP_DIR="${INTEL_ROOT}/compilers_and_libraries_2019.4.243/linux"
                      export MPI_DIR="${INTEL_ROOT}/impi/2019.4.243/intel64"
                      export NETCDF_DIR="/store/CT1/hmg2840/lbrodeau/opt/hdf5_netcdf4_intel_mpi"
                      export NETCDF_CXX_DIR=${NETCDF_DIR}
                      #
                      export AEROBULK_DIR="/store/CT1/hmg2840/lbrodeau/DEV/aerobulk"
                      export OASIS_DIR="/store/CT1/hmg2840/lbrodeau/src/oasis3-mct"
                      #
                      export LD_EXTRA_OASIS="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -lifcore -lifport"
                      export LDFLAGS="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -lifcore -lifport"
                      #export FFLAGS="-O2 -qopenmp -lstdc++ -fPIC -I`which ifort | sed -e "s|bin/intel64/ifort|compiler/include/intel64|g"`" ; # only for C++ to find module "iso_c_binding"
                      ;;
    "fram" )          NXTSM_DEP_DIR="/cluster/projects/nn9878k/brodeau/opt/nextsim_intel"
                      INTEL_VERSION="2018.1.163"
                      INTEL_ROOT="/cluster/software/ifort/${INTEL_VERSION}-GCC-6.4.0-2.28"
                      export INTEL_COMP_DIR="${INTEL_ROOT}/compilers_and_libraries_${INTEL_VERSION}/linux"
                      export MPI_DIR="/cluster/software/impi/${INTEL_VERSION}-iccifort-${INTEL_VERSION}-GCC-6.4.0-2.28/intel64"
                      export NETCDF_DIR="/cluster/software/netCDF/4.4.1.1-intel-2018a-HDF5-1.8.19"
                      export NETCDF_CXX_DIR="/cluster/software/netCDF-C++4/4.3.0-intel-2018a-HDF5-1.8.19"
                      #
                      export AEROBULK_DIR="/cluster/projects/nn9878k/brodeau/src/aerobulk"
                      export OASIS_DIR="/cluster/projects/nn9878k/brodeau/src/oasis3-mct"
                      #
                      #export LD_EXTRA_OASIS="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -lifcore -lifport"
                      export LDFLAGS="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -lifcore -lifport"
                      #export FFLAGS="-O2 -qopenmp -lstdc++ -fPIC -I`which ifort | sed -e "s|bin/intel64/ifort|compiler/include/intel64|g"`" ; # only for C++ to find module "iso_c_binding":
                      ;;

    *               ) echo;
                      echo "WARNING: Unknow machine with HOSTNAME = `hostname` !"
                      echo
                      sleep 2
                      ;;
esac

# Normally the following 2 are pretty standard:
export MPI_LIB_DIR=${MPI_DIR}/lib
export MPI_INC_DIR=${MPI_DIR}/include


#################### Loading INTEL system #############################
export PATH=${INTEL_COMP_DIR}/bin/intel64:${PATH}
export LD_LIBRARY_PATH=${INTEL_COMP_DIR}/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
export CPATH=${INTEL_COMP_DIR}/include:${CPATH}
#
# M P I
export PATH=${MPI_DIR}/bin:${PATH}
export LD_LIBRARY_PATH=${MPI_LIB_DIR}:${LD_LIBRARY_PATH}
export CPATH=${MPI_INC_DIR}:${CPATH}
#######################################################################



# Third-party software dependencies, compiled with relevant compiler!
export GMSH_DIR=${NXTSM_DEP_DIR}/gmsh-${GMSH_VERSION}

export BOOST_DIR=${NXTSM_DEP_DIR}/boost-${BOOST_VERSION}
export BOOST_INCDIR=${BOOST_DIR}/include
export BOOST_LIBDIR=${BOOST_DIR}/lib

#export LD_EXTRA_AEROBULK="-L${INTEL_DIR}/compiler/lib/intel64_lin -lifcore -lifport"
#export LD_EXTRA_AEROBULK="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -lifcore"

if [ ! -f ./model/version.hpp ]; then
    cd model
    ./version.sh
    cd ../
fi
