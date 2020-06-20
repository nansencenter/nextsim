# makefile includes for hexagon, portland compilers
# Standard compilers and linkers
COMP=/opt/local/bin
CF90=$(COMP)/gfortran-mp-6
CF77=$(COMP)/gfortran-mp-6
CC  =$(COMP)/g++-mp-6
LD  =$(COMP)/gfortran-mp-6

# CPP, used internally by compilers - note that some compilers
# need special options for this, consult compiler manual if in trouble
CPP =/usr/bin/cpp -P -traditional

# Flags for fixed/free format
#F90FLG= -Mfree
F90FLG= -free
#F90FLG= -ffree-form
#F77FLG= -Mfixed
#F90FLG= -ffixed-form
#F77FLG= -ffixed-form

# Compiler flags, real8 version
#FLAGS= -fast -fastsse -byteswapio -g -r8 -Minform=inform
#FFLAGS= -O2 -fconvert=big-endian -fdefault-real-8 -fcheck=all -g -fbacktrace
# use -C for debugging
FFLAGS= -O2 -fconvert=big-endian -fdefault-real-8
CFLAGS= -O2 -fast -fastsse -Wall 

#SFLAGS= -dynamiclib
SFLAGS= -shared

#Netcdf, FFTW and lapack Libraries
#NETCDF_DIR = /y1/aydogdu/Anaconda/anaconda/pkgs/netcdf-fortran-4.4.4-3
NETCDF_DIR = /opt/local
NETCDF_INC = $(NETCDF_DIR)/include
NETCDF_LIB = $(NETCDF_DIR)/lib

FFTW_DIR   = /y1/aydogdu/Anaconda/anaconda/pkgs/fftw-3.3.7-0
FFTW_DIR   = /usr/local
FFTW_INC   = $(FFTW_DIR)/include
FFTW_LIB   = $(FFTW_DIR)/lib

BOOST_LIBDIR=/opt/local/boost/lib
BOOST_INCDIR=/opt/local/boost/include

INCLUDE    = -I$(FFTW_INC) -I$(NETCDF_INC)
LIBS       = -L$(FFTW_LIB) -lfftw3 -L$(NETCDF_LIB) -lnetcdf -lnetcdff 

CINCLUDE   = -I$(BOOST_INCDIR) -I../include
CLIBS      = -L../lib -lpseudo2D -L$(BOOST_LIBDIR) -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi
LINKFLAGS  = -L/opt/local/lib/gcc6 -lstdc++ -lgfortran

# Some fortran compilers have iargc as built in, 
# others as library routine. NB: CPP flags and libs mustr be consistent
CPPFLAGS=-DIARGC -DFFTW -DLAPACK
