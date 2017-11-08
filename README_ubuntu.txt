export BOOST_DIR=/usr/lib/x86_64-linux-gnu
# =============================================================


# =============================================================
CLONE GIT REPO
git clone https://tdcwilliams@bitbucket.org/asamake/nextsim.git
# =============================================================


# =============================================================
BOOST
sudo apt-get install libboost-dev
sudo apt-get install libboost-chrono-dev
sudo apt-get install libboost-date-time-dev
sudo apt-get install libboost-filesystem-dev
sudo apt-get install libboost-iostreams-dev
sudo apt-get install libboost-log-dev
sudo apt-get install libboost-locale-dev
sudo apt-get install libboost-math-dev
sudo apt-get install libboost-mpi-dev
sudo apt-get install libboost-program-options-dev
sudo apt-get install libboost-regex-dev
sudo apt-get install libboost-serialization-dev
sudo apt-get install libboost-system-dev
sudo apt-get install libboost-timer-dev
# =============================================================


# =============================================================
COMPILERS
*gcc 4.8
sudo apt-get install gcc-4.8
sudo apt-get install libgcc-4.8-dev
cd /usr/bin
sudo ln -sf gcc-4.8 gcc
sudo ln -sf gcc-ar-4.8 gcc-ar
which gcc
which mpicc
gcc --version
mpicc --version

*gfortran 4.8
sudo apt-get install gfortran-4.8
sudo ln -sf gfortran-4.8 gfortran
# =============================================================


# =============================================================
OTHER TOOLS/LIBRARIES
sudo apt-get install cmake valgrind
sudo apt-get install liblapack-dev libblas-dev
# =============================================================


# =============================================================
PETSC
1) download v3.6.4 from:
http://www.mcs.anl.gov/petsc/download/
- After extraction of files, make two copies
- one for normal running and one for debugging
2) copy scripts/petsc/pconfigure.sh to this dir
3) unset PETSC_DIR
   unset PETSC_ARCH
4) sh pconfigure.sh
5) follow instructions for compilation/installation
eg:
a) make PETSC_DIR=/home/timill/Packages/cpp/petsc-3.6.4 PETSC_ARCH=arch-linux2-c-opt all
b) sudo make PETSC_DIR=/home/timill/Packages/cpp/petsc-3.6.4 PETSC_ARCH=arch-linux2-c-opt install
5) In .bash_profile, add:
export PETSC_ARCH=arch-linux2-c-opt
export PETSC_DIR=/opt/local/petsc

PETSC-DEBUG
1) copy scripts/petsc/pconfigure_debug.sh to this dir
2) unset PETSC_DIR
   unset PETSC_ARCH
3) edit this file if necessary, then do
   sh pconfigure_debug.sh
4) follow instructions for compilation/installation
eg
a) make PETSC_DIR=/home/timill/Packages/cpp/petsc-3.6.4-debug PETSC_ARCH=arch-linux2-c-debug all
b) sudo make PETSC_DIR=/home/timill/Packages/cpp/petsc-3.6.4-debug PETSC_ARCH=arch-linux2-c-debug install
5) in .bash_profile add
export PETSC_ARCH=arch-linux2-c-debug
export PETSC_DIR=/opt/local/petsc-debug
# =============================================================


# =============================================================
GMSH
1) svn --username=gmsh co https://onelab.info/svn/gmsh/trunk gmsh
      password: gmsh
2) mkdir Build
   cd Build
   cp $NEXTSIMDIR/scripts/gmsh/gconfigure.sh .
   ./gconfigure.sh
3) make
   sudo make install
# =============================================================


# =============================================================
NETCDF
1) sudo apt-get install libhdf5-dev libhdf5-openmpi-dev
2) download latest stable c version of netcdf from
   http://www.unidata.ucar.edu/downloads/netcdf/index.jsp and unzip it
3) copy configure_ubuntu.sh from /nextsim/scripts/netcdf
   ./configure_ubuntu.sh
   make
   make check
   sudo make install
# =============================================================


# =============================================================
NETCDF-CXX
1) download latest stable c++ version of netcdf from
   http://www.unidata.ucar.edu/downloads/netcdf/index.jsp and unzip it
2) copy configure_ubuntu.sh from /nextsim/scripts/netcdf-cxx
   ./configure_ubuntu.sh
   make
   make check
   sudo make install
3) sudo cp /opt/local/netcdf/include/* /opt/local/netcdf-cxx/include
# =============================================================


# =============================================================
ENVIRONMENT VARIABLES
# Add these lines (with correct paths) to your .bash_profile:

# for nextSIM in C++
export NEXTSIMDIR=$HOME/Developer/nextsim/
export GMSH_DIR=/opt/local/gmsh

export NETCDF_DIR=/opt/local/netcdf-cxx

export PETSC_DIR=/opt/local/petsc
export PETSC_ARCH=arch-linux2-c-opt
# export PETSC_DIR=/opt/local/petsc-debug
# export PETSC_ARCH=arch-linux2-c-debug

export BOOST_INCDIR=/usr/include/boost
export BOOST_LIBDIR=/usr/lib/x86_64-linux-gnu

export OPENMPI_LIB_DIR=/usr/lib/x86_64-linux-gnu
export OPENMPI_INCLUDE_DIR=/usr/lib/openmpi/include
export OPENMPI_INCLUDE_DIR=/usr/lib/openmpi/lib

export LD_LIBRARY_PATH="$BOOST_LIBDIR"
export LD_LIBRARY_PATH=$NEXTSIMDIR/lib:$LD_LIBRARY_PATH
# =============================================================


# =============================================================
COMPILE NEXTSIM
1) compile libraries
cd $NEXTSIMDIR
make
2) compile model
cd $NEXTSIMDIR/model
make
# =============================================================


# =============================================================
NEXTWIM
1) in .bash_profile, add
export USE_NEXTWIM=1
(compiler just checks if it's defined)
2) now compiled automatically with first "make"
3) To run with it set
simul.use_wim=True (as in coupling_wim.cfg)
# =============================================================


# =============================================================
MATLAB

Installation with Bjerknes license
1) get uib email, username and password
2) get instructions and download activation code from
https://it.uib.no/en/Matlab

R2017b:
# complains it needs g++-4.9 for mex to compile c++,
# but it still compiles with g++-5
# c++ code itself won't compile with g++-4.8/4.9 (boost)

3) Test mex compilation with:
   cd matlab
   make
# =============================================================
