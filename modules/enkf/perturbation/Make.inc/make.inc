# makefile includes for hexagon, portland compilers
# Standard compilers and linkers
CF90= ftn
CF77= ftn
CC= cc
LD= ftn

# CPP, used internally by compilers - note that some compilers
# need special options for this, consult compiler manual if in trouble
CPP= /usr/bin/cpp -P -traditional

#NCAR graphics compiler wrappers
NCARGCF90=ncargf90
NCARGCF77=ncargf77
NCARGCC=ncargcc
NCARGLD=ncargf90

# Flags for fixed/free format
F90FLG= -Mfree
F77FLG= -Mfixed

## uncomment below for debugging and set MPI = NO in makefile
#  CF90 = pgfortran
#  CF77 = pgfortran
#  LD = pgfortran
#  DEBUG_FLAGS = -g -Minform=inform
#  PAR0 =
#  OPTO = 
#  INC_NETCDF = -I/local/netcdf-3.6.2-pgi/include
#  LIB_NETCDF = /local/netcdf-3.6.2-pgi/lib/libnetcdf.a
#  INC_FFTW = -I/opt/fftw/3.2.2/include
#  LIB_FFTW = /opt/fftw/3.2.2/lib/libfftw3.a
#
# Compiler flags, real8 version
#FLAGS= -fast -fastsse -byteswapio -g -r8 -Minform=inform
# use -C for debugging
#FFLAGS= -fast -fastsse -byteswapio -g -r8 -C
FFLAGS= -fast -fastsse -byteswapio -g -r8 
CFLAGS= -fast -fastsse
#FFLAGS= -byteswapio -g -r8 -Minform=inform
#CFLAGS= 

# Compiler flags, real4 version (needed for NCARG routines)
FFLAGSR4= -fast -fastsse -byteswapio -g -r4
CFLAGSR4= -fast -fastsse

# Linker flags
LINKFLAGS= $(FFLAGS)

#Netcdf, FFTW and lapack Libraries
NETCDF_INC =


INCLUDE= $(NETCDF_INC) $(FFTW_INCLUDE)
LIBS =  $(NETCDF_POST_LINK_OPTS) $(FFTW_POST_LINK_OPTS) $(LIB_LAPACK)

# Some fortran compilers have iargc as built in, 
# others as library routine. NB: CPP flags and libs mustr be consistent
CPPFLAGS=-DIARGC -DFFTW -DLAPACK
