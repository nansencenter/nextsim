export BOOST_DIR=/usr/lib/x86_64-linux-gnu
# =============================================================


# =============================================================
CLONE GIT REPO
git clone https://tdcwilliams@bitbucket.org/asamake/nextsim.git
# =============================================================


# =============================================================
BOOST

Compile from source:
1) Download version 1.60 of boost on http://www.boost.org
   (It is better to restart from here if you had an upgrade of your os)
2) copy bconfigure_ubuntu.sh and binstall_ubuntu.sh from /nextsim/scripts/boost/
   to your boost directory
3) type the command "unset BOOST_DIR"
4) check if bconfigure.sh corresponds to your architecture
5) From boost directory, type: "./bconfigure.sh"
6) Add the line "using mpi ;" at the end of the file project-config.jam
7) From boost directory, type: "./binstall.sh" (sudo password is required during the process)

Compile with apt:
NB only boost v1.58 is available so won't run parallel code
sudo apt-get install libboost-dev \
                     libboost-chrono-dev \
                     libboost-date-time-dev \
                     libboost-filesystem-dev \
                     libboost-iostreams-dev \
                     libboost-log-dev \
                     libboost-locale-dev \
                     libboost-math-dev \
                     libboost-mpi-dev \
                     libboost-program-options-dev \
                     libboost-regex-dev \
                     libboost-serialization-dev \
                     libboost-system-dev \
                     libboost-timer-dev
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

# Comment: It can also be compiled with the default compilers from APT:
# sudo apt install build-essential
# GNU Fortran (Ubuntu 5.4.0-6ubuntu1~16.04.5) 5.4.0 20160609
# gcc (Ubuntu 5.4.0-6ubuntu1~16.04.5) 5.4.0 20160609
# =============================================================


# =============================================================
OTHER TOOLS/LIBRARIES
sudo apt-get install cmake \
                     valgrind \
                     tau \
                     liblapack-dev \
                     libblas-dev \
                     nco \
                     ncview \
                     libnetcdf-dev \
                     libnetcdf-c++4-dev \
                     subversion
# =============================================================


# =============================================================
GMSH
1) svn --username=gmsh co https://onelab.info/svn/gmsh/trunk gmsh
      password: gmsh
2) cd gmsh
   mkdir Build
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
3) copy configure_ubuntu.sh from nextsim/scripts/netcdf
   ./configure_ubuntu.sh
   make
   make check
   sudo make install

# Comment: it can also be compiled with the default version from APT 
# sudo apt install libnetcdf-dev
# libnetcdf-dev Version: 1:4.4.0-2
# =============================================================


# =============================================================
NETCDF-CXX
1) download latest stable c++ version of netcdf from
   http://www.unidata.ucar.edu/downloads/netcdf/index.jsp and unzip it
2) copy configure_ubuntu.sh from nextsim/scripts/netcdf-cxx
   ./configure_ubuntu.sh
   make
   make check
   sudo make install
3) sudo cp /opt/local/netcdf/include/* /opt/local/netcdf-cxx/include

# Comment: it can also be compiled with the default version from APT 
# sudo apt install libnetcdf-c++4-dev
# libnetcdf-c++4-dev Version: 4.2.1-3
# =============================================================


# =============================================================
ENVIRONMENT VARIABLES
# Add these lines (with correct paths) to your .bash_profile:

# for nextSIM in C++
export NEXTSIMDIR=$HOME/Developer/nextsim/
export GMSH_DIR=/opt/local/gmsh

export NETCDF_DIR=/opt/local/netcdf-cxx

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
