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

### Installation (that has worked with osX 10.10) ###

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
	3) Add the line “using mpi ;” at the end of the file tools/build/v2/user-config.jam (This file could also be copied in your home directory for boost from the version 1.56)
	4) type the command “unset BOOST_DIR”
	5) check if bconfigure.sh corresponds to your architecture
	6) From boost directory, type: “./bconfigure.sh”
	7) From boost directory, type: “./binstall.sh” (sudo password is required during the process)

####### Install petsc #######
	0) Desinstall other version of petsc by doing (It is better to restart from here if you had an upgrade of your os)
			sudo rm -rf /opt/local/petsc
		or by changing the prefix in pconfigure.sh, for example:
			--prefix=/opt/local/petsc/3.6	
	0bis) Desinstall older installation of suitesparse/cholmod/umfpack in /usr/local (It depends on how you installed it. For example you will have to type: make uninstall in SuiteSparse/, or to use: sudo brew remove suite-sparse)
	1) Download the version 3.6.1 of petsc 
	from http://www.mcs.anl.gov/petsc/download/
	or from http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/
	2) copy the content of /nextsim/scripts/petsc/ in petsc directory
	3) type the command “unset PETSC_DIR”
	4) check if pconfigure.sh corresponds to your architecture and change the permission to make it executable.
	5) From petsc directory, type: “./pconfigure.sh” (sudo password is required during the process)
	6) From petsc directory, follow the instruction for the rest of the compilation (something like: “make PETSC_DIR=/Users/syloui/Developer/petsc-3.6.1 PETSC_ARCH=arch-darwin-c-opt all”, “sudo make PETSC_DIR=/Users/syloui/Developer/petsc-3.6.1 PETSC_ARCH=arch-darwin-c-opt install” and “make PETSC_DIR=/opt/local/petsc PETSC_ARCH="" test”)


####### Install netcdf #######
1) instal hdf5 via macport: "sudo port install hdf5"
2) download latest stable c version of netcdf on http://www.unidata.ucar.edu/downloads/netcdf/index.jsp and unzip it
3) copy configure_c.sh from /nextsim/scripts/netcdf
4) From netcdf directory, type: “./configure_c.sh” then "make" "make check" and finally "sudo make install"
5) download latest stable c-xx version of netcdf on http://www.unidata.ucar.edu/downloads/netcdf/index.jsp and unzip it
6) copy configure_cxx.sh from /nextsim/scripts/netcdf
7) From netcdf-cxx directory, type: “./configure_cxx.sh”then "make" "make check" and finally "sudo make install"

####### install gmsh from the source code #######

1) Download the source code from http://geuz.org/gmsh/

2) In the gmsh directory do

mkdir lib
cd lib
cmake -DDEFAULT=0 -DENABLE_BUILD_LIB=1 ..
make -j 32 lib
sudo make install/fast

####### Set the PATH correctly #######

# Add these lines to your .bash_profile:

# for nextSIM in C++
export NEXTSIMDIR=$HOME/Developer/nextsim/
export GMSH_DIR=/usr/local/

export NETCDF_DIR=/opt/local/netcdf-cxx

export PETSC_DIR=/opt/local/petsc
export PETSC_ARCH=arch-darwin-c-opt

export BOOST_DIR=/opt/local/boost

export DYLD_LIBRARY_PATH="/opt/local/boost/lib"
export DYLD_LIBRARY_PATH=NEXTSIMDIR/lib:$DYLD_LIBRARY_PATH

##-------  Compile neXtSIM itself --------- 

# Open a new command window

# Either modify the Makefile to have the right link to openmpi or do the following:

sudo rm -rf /opt/local/include/openmpi-mp
sudo ln -sf /opt/local/include/openmpi-gcc48 /opt/local/include/openmpi-mp

sudo rm -rf /opt/local/lib/openmpi-mp
sudo ln -sf /opt/local/lib/openmpi-gcc48 /opt/local/lib/openmpi-mp

# Type make in nextsim

## For the mode application (run nextsim)

# go to nextsim/model

# Type make

# Go to bin, type ./nextsim.exec

## run

# type “bin/nextsim.exec” from the nextsim/research directory
