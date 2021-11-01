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
export NETCDF_DIR="/opt/netcdf_intel"
export NETCDF_CXX_DIR="/opt/netcdf_intel"
#
export AEROBULK_DIR="${HOME}/DEV/aerobulk"
export OASIS_DIR="${HOME}/src/oasis3-mct"
#
NXTSM_DEP_DIR="/opt/nextsim_intel" ; # path to directory containing compiled BOOST and GMSH (with the relevant compiler!)
#
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
                      export LDFLAGS="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -lifcore -lifport"
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
                      export LDFLAGS="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -lifcore -lifport"
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
                      export LDFLAGS="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -lifcore -lifport"
                      ;;
    *               ) echo;
                      echo "WARNING: Unknow machine with HOSTNAME = `hostname` !"
                      echo
                      sleep 2
                      ;;
esac

if ${l_aerobulk} || ${l_cpl_oasis}; then
    export LDFLAGS="-L${INTEL_COMP_DIR}/compiler/lib/intel64_lin -lifcore -lifport"
fi

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

if [ ! -f ./model/version.hpp ]; then
    cd model
    ./version.sh
    cd ../
fi
