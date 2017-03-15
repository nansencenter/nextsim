# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact

### Installation (that has worked with osX 10.10 and 10.11) ###

##------- Installation of the needed compilers ----------

####### Install macport #######
	by following the instructions on https://www.macports.org/install.php
	Check that "/opt/local/bin” and “/opt/local/sbin” are defined in PATH (“echo $PATH” in a new command window)

####### Install gcc48 (or more recent versions) with macport #######
	by doing "sudo port search gcc”, "sudo port install gcc48”, "sudo port select --list gcc" and “sudo port select --set gcc mp-gcc48"

####### Install openmpi-gcc48 (or a version compatible with your gcc) with macport #######
	by doing "sudo port search openmp”,"sudo port install openmpi-gcc48”, "sudo port select --list mpi”, “sudo port select --set mpi openmpi-gcc48-fortran” and “sudo ln -s /opt/local/bin/mpicxx /opt/local/bin/mpic++”


##------- Installation of the needed libraries ----------

####### Install boost #######
	1) Download version 1.55 of boost on http://www.boost.org (It is better to restart from here if you had an upgrade of your os)
	2) copy bconfigure.sh and binstall.sh from /nextsim/scripts/boost/ to your boost directory
	3) Add the line “using mpi ;” at the end of the file tools/build/v2/user-config.jam
	NOte: if you install the version 1.56 or more recent of Boost, the file user-config.jam is locate in tools/build/example and you NEED to copy this file in your home directory before the compilation.
	4) type the command “unset BOOST_DIR”
	5) check if bconfigure.sh corresponds to your architecture
	6) From boost directory, type: “./bconfigure.sh”
	7) From boost directory, type: “./binstall.sh” (sudo password is required during the process)

####### Install petsc #######
    First, install cmake with macports if not already installed: "sudo port install cmake" (note: make sure this is used and not ISSM version)
    0) Download from http://www.mcs.anl.gov/petsc/download/
	   or from http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/
            - choose 3.6.XX version
    1) Check mpicc, mpicxx, mpiexec use correct version of gcc (gcc --version)
	2) Uninstall other version of petsc by doing (It is better to restart from here if you had an upgrade of your os)
			sudo rm -rf /opt/local/petsc
		or by changing the prefix in pconfigure.sh, for example:
			--prefix=/opt/local/petsc/3.6
	3) copy the contents of /nextsim/scripts/petsc/ in petsc directory
	4) type the command “unset PETSC_DIR”
	5) check if pconfigure.sh corresponds to your architecture and change the permission to make it executable.
	6) From petsc directory, type: “./pconfigure.sh” (sudo password is required during the process)
	7) From petsc directory, follow the instruction for the rest of the compilation (something like: “make PETSC_DIR=/Users/syloui/Developer/petsc-3.6.1 PETSC_ARCH=arch-darwin-c-opt all”, “sudo make PETSC_DIR=/Users/syloui/Developer/petsc-3.6.1 PETSC_ARCH=arch-darwin-c-opt install” and “make PETSC_DIR=/opt/local/petsc PETSC_ARCH="" test”)

Note: for debugging, use another copy of petsc compiled with --download-mpich and --with-debugging and using another path, for example using export PETSC_PREFIX=/opt/local/petsc-debug. Set this path in your .bashprofile by using export PETSC_DIR=/opt/local/petsc-debug before recompiling the whole neXtSIM code (core + model).

You can run the code with the options:
-malloc_debug for memory error

####### Install netcdf #######
1) install hdf5 via macport: "sudo port install hdf5"
2) download latest stable c version of netcdf on http://www.unidata.ucar.edu/downloads/netcdf/index.jsp and unzip it
3) copy configure_c.sh from /nextsim/scripts/netcdf
4) From netcdf directory, type: “./configure_c.sh” then "make", "make check" and finally "sudo make install"
5) download latest stable c-xx version of netcdf on http://www.unidata.ucar.edu/downloads/netcdf/index.jsp and unzip it
6) copy configure_cxx.sh from /nextsim/scripts/netcdf
7) From netcdf-cxx directory, type: “./configure_cxx.sh”then "make" "make check" and finally "sudo make install"

####### install gmsh from the source code #######

1) Download the source code from svn
svn co https://onelab.info/svn/gmsh/trunk gmsh
username: gmsh
password: gmsh

2) In the gmsh directory do

i) full install (OSX):
* mkdir Build
edit CMakeLists.txt:
* paste set(CMAKE_MACOSX_RPATH 1)
  under set(CMAKE_LEGACY_CYGWIN_WIN32 0)
* cd Build
* cmake ..
* ccmake ..
   change  CMAKE_BUILD_TYPE      to Release
   change  CMAKE_INSTALL_PREFIX  to /opt/local/gmsh
   change  ENABLE_BUILD_DYNAMIC  to ON
   change  ENABLE_BUILD_LIB      to ON
   change  ENABLE_BUILD_SHARED   to ON
   "c" to configure
   "g" to generate and exit
* make -j8
* sudo make install
* shouldn't need to do this, but...
   sudo mkdir /opt/local/gmsh/lib
   sudo cp libGmsh.a *.dylib /opt/local/gmsh/lib

ii) quick install:
mkdir lib
cd lib
cmake -DDEFAULT=0 -DENABLE_BUILD_LIB=1 ..
make -j 32 lib
sudo make install/fast

####### Set the PATH correctly #######

# Add these lines to your .bash_profile:

# for nextSIM in C++
export NEXTSIMDIR=$HOME/Developer/nextsim/
export GMSH_DIR=/opt/local/gmsh or /usr/local

export NETCDF_DIR=/opt/local/netcdf-cxx

export PETSC_DIR=/opt/local/petsc
export PETSC_ARCH=arch-darwin-c-opt

export BOOST_DIR=/opt/local/boost

export DYLD_LIBRARY_PATH="/opt/local/boost/lib"
export DYLD_LIBRARY_PATH=NEXTSIMDIR/lib:$DYLD_LIBRARY_PATH

# Check your environment variables
Check if you have the environment variable WIM2D_PATH (echo $WIM2D_PATH)?
if so, you will need to unset it (unset WIM2D_PATH)

##-------  Compile neXtSIM itself ---------

# Open a new command window

# Either modify the Makefile to have the right link to openmpi or do the following:

sudo rm -rf /opt/local/include/openmpi-mp
sudo ln -sf /opt/local/include/openmpi-gcc48 /opt/local/include/openmpi-mp

sudo rm -rf /opt/local/lib/openmpi-mp
sudo ln -sf /opt/local/lib/openmpi-gcc48 /opt/local/lib/openmpi-mp

# Go to $NEXTSIMDIR and type "make"

##------- For the model application -------

# go to $NEXTSIM_DIR/model

# Type "make"

# type “bin/nextsim.exec --configfile=nextsim.cfg”

##-------  Compile neXtWIM itself ---------

First go to the WIM2D github repo ( = $WIM2D_PATH - probably need to set
this)

Compile the WIM library with:

cd CXX/Build
make vclean
make lib

back in nextsim/model

compile with

make clean
make wim

simul.use_wim=True (as in coupling_wim.cfg)

##------- Compile with OASIS3-MCT support -------

You must have OASIS3-MCT installed. Instructions to come, but the neccesary libraries are installed on johansen in /Data/sim/packages/oasis3-mct
You must have Fortran netCDF libraries installed. Instructions to come, but the neccesary libraries are installed on johansen in /Data/sim/packages/netcdf-fortran

In order to compile with OASIS3-MCT support you must set the following environment variables: USE_OASIS, OASIS_DIR, and NETCDF_FOR_DIR. On johansen do:
export USE_OASIS="true" # or any non-empty value
export OASIS_DIR=/Data/sim/packages/oasis3-mct/
export NETCDF_FOR_DIR=/Data/sim/packages/netcdf-fortran

Compile the OASIS C++ module:
cd $NEXTSIMDIR
make

Compile the model with OASIS3-MCT support:
cd $NEXTSIMDIR/model
make

##-------Debugging-------

# install gdb on mac following those links:
http://www.patosai.com/blog/post/installing-gdb-on-mac-os-x-yosemite

# install valgrind on mac following those links:
http://ranf.tl/2014/11/28/valgrind-on-mac-os-x-10-10-yosemite/
http://superuser.com/questions/630674/valgrind-installation-errors-on-osx-10-8

##-------Profiling-------
Comparison of available profilers:
http://gernotklingler.com/blog/gprof-valgrind-gperftools-evaluation-tools-application-level-cpu-profiling-linux/

# install google perf tools on mac following those links:
http://brewformulas.org/Google-perftool

You need a recent version of XQuartz (X11)

You also need graphviz
install graphviz

You also need gv:
brew install -v ghostscript gv

Then follow read this how to:
http://goog-perftools.sourceforge.net/doc/cpu_profiler.html

For CPU analysis, I think we need to compile with at least -g option and link with -lprofiler:

I tried adding
-pg -g -DNDEBUG to CXXFLAGS
-lprofiler to LDFLAGS

Option 1)
Then run it by:
CPUPROFILE=profile.log ./bin/nextsim.exec --config-files=nextsim_BK.cfg

Option 2)
You can apply the profiler on a small part of the code by compiling with -DWITHGPERFTOOLS and adding those lines where needed:
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#ifdef WITHGPERFTOOLS
    ProfilerStart("profile.log");
#endif

#ifdef WITHGPERFTOOLS
    ProfilerStop();
#endif

You then simply run:
./bin/nextsim.exec --config-files=nextsim_BK.cfg
at it will save the log info into profile.log file.

To print the graph, you can use :
pprof --gv ./bin/nextsim.exec profile.log

You can also have a text file :
pprof --text ./bin/nextsim.exec profile.log


Unfortunately, I haven't manage to output correctly the symbolic information, even when adding -g option to the core and contrib code, or even when using petsc in debud mode. I only get info with the adresses of the functions which is not very usefull. For example:
618   0.8%  26.1%      618   0.8% 0x000000010213ed6e
    616   0.8%  26.9%      616   0.8% 0x000000010218586a
    581   0.7%  27.6%      581   0.7% 0x00007fff93faa348
    568   0.7%  28.3%      568   0.7% 0x0000000102185807
    524   0.7%  29.0%      524   0.7% 0x00007fff946e448a
    397   0.5%  29.5%      397   0.5% 0x00000001021d7b5b
    393   0.5%  30.0%      393   0.5% 0x00000001021ed9b4
    387   0.5%  30.5%      387   0.5% 0x00007fff93fa7e0f
    376   0.5%  30.9%      376   0.5% 0x00007fff93fabcd1



###### Use gperftools ######

# before compiling the libraries and model, we need first
export NEXTSIM_BUILD_TYPE=Debug

# Then, we can compile and run nextsim

# Finally, run one of the following lines

You can also have a pdf file
pprof --pdf --functions --focus=run --cum --drop_negative --nodecount=50 bin/nextsim.exec profile.log > profile.pdf
or
pprof --pdf --functions --focus=run bin/nextsim.exec profile.log > profile.pdf
or
pprof --pdf --functions --focus=run --edgefraction=1e-02 --nodefraction=1e-02 bin/nextsim.exec profile.log > profile.pdf
or
pprof --pdf --functions --focus=run --cum --drop_negative --nodecount=20 bin/nextsim.exec profile.log > profile.pdf
