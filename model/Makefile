KERNEL := $(shell uname -s)

PROGNAME=nextsim.exec

CXX = g++

# setting the C++ standard according to the gcc compiler version (from gcc-5.2, the default is C++14)
ifeq ($(shell echo `$(CXX) -dumpversion | cut -f1-2 -d.` \>= 5.2 | sed -e 's/\.//g' | bc),1)
	CXXFLAGS += -std=c++14
else
	CXXFLAGS += -std=c++11
endif

# add g++ option flags
CXXFLAGS += -ftemplate-depth-256 -Wno-inline \
		-fPIC -fopenmp \
		-DHAVE_CONFIG_H -D_MULTITHREADING_
ifneq ($(XSIMCOMP_VERBOSE),no)
	CXXFLAGS += -v
endif

ifdef USE_OASIS
	CXXFLAGS += -DOASIS
	CXXFLAGS += -I $(NEXTSIMDIR)/modules/oasis/include
	LDFLAGS += -lgfortran
	LDFLAGS += -L $(NEXTSIMDIR)/lib -loasis
	CHAN = MPI1
	#LIBPSMILE = $(OASIS_DIR)/lib/libpsmile.${CHAN}.a $(OASIS_DIR)/lib/libmct.a $(OASIS_DIR)/lib/libmpeu.a $(OASIS_DIR)/lib/libscrip.a
endif

ifneq (,$(strip $(filter DEBUG Debug debug PROFILE Profile profile,$(NEXTSIM_BUILD_TYPE))))
	#ifeq ($(NEXTSIM_BUILD_TYPE),$(filter $(NEXTSIM_BUILD_TYPE),Debug debug))
	CXXFLAGS := $(filter-out -O3 -pthread,$(CXXFLAGS))
	CXXFLAGS += -g -O0 -DNDEBUG
ifneq (,$(strip $(filter PROFILE Profile profile,$(NEXTSIM_BUILD_TYPE))))
	CXXFLAGS += -DWITHGPERFTOOLS
endif
ifneq ($(KERNEL),Linux)
	CXXFLAGS += -Wl,-no_pie
endif
else
	CXXFLAGS += -O3 -pthread
endif

# add include paths
CXXFLAGS += -I $(NEXTSIMDIR)/core/include
CXXFLAGS += -isystem $(NEXTSIMDIR)/contrib/bamg/include # suppress annoying compilation warnings from -I
CXXFLAGS += -isystem $(NEXTSIMDIR)/contrib/mapx/include # suppress annoying compilation warnings from -I
# CXXFLAGS += -I $(NEXTSIMDIR)/contrib/interp/include

ifdef USE_ENSEMBLE
	CXXFLAGS += -DENSEMBLE
	CXXFLAGS += -I $(NEXTSIMDIR)/modules/enkf/perturbation/include
endif

ifdef USE_AEROBULK
		CXXFLAGS += -I $(NEXTSIMDIR)/modules/aerobulk/include
		CXXFLAGS += -DAEROBULK
endif

# openmpi
CXXFLAGS += -I $(OPENMPI_INCLUDE_DIR)/

# petsc
CXXFLAGS +=	-I $(PETSC_DIR)/include

# boost
CXXFLAGS += -I $(BOOST_INCDIR)/ -I .

# netcdf
CXXFLAGS +=	-I $(NETCDF_DIR)/include

# gmsh
CXXFLAGS += -I $(GMSH_DIR)/include/gmsh

CXXFLAGS +=	-I /opt/local/include

ifeq ($(KERNEL),Linux)
	#CXXFLAGS += -std=c++0x -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -g -lm -pthread -v #-MMD -MP -lm -pthread -v

else

	ifeq ($(CXX),clang)
	CXXFLAGS += -stdlib=libc++
	endif

	CXXFLAGS += -I /usr/local/include #-I /opt/local/include/openmpi-mp

	#LDFLAGS += -Wl,-rpath,/usr/local/lib #-Wl,-rpath,/opt/local/lib/openmpi-mp
	#LDFLAGS += -L /usr/local/lib #-L /opt/local/lib/openmpi-mp -lmpi_cxx -lmpi -ldl -lstdc++ -lpthread

	ifeq ($(CXX),clang)
	LDFLAGS += -stdlib=libc++
	endif

endif

LDFLAGS += -L /usr/local/lib

LDFLAGS += -Wl,-rpath,$(OPENMPI_LIB_DIR)/
ifndef MACHINE_HEXAGON
	LDFLAGS += -L $(OPENMPI_LIB_DIR)/ -lmpi_cxx -lmpi -ldl -lstdc++ #-lpthread
else
	LDFLAGS += -L $(OPENMPI_LIB_DIR)/ -lmpichcxx -lmpich -ldl -lstdc++ #-lpthread #-lssl -luuid -lpthread -lrt
        LDFLAGS += -Wl,-rpath,$(BLAS_LAPACK_DIR)/lib
        LDFLAGS += -L $(BLAS_LAPACK_DIR)/lib -lsci_gnu_mp
endif

LDFLAGS += -Wl,-rpath,$(NETCDF_DIR)/lib -L $(NETCDF_DIR)/lib -lnetcdf_c++4

LDFLAGS += -Wl,-rpath,$(BOOST_LIBDIR)
LDFLAGS += -L $(BOOST_LIBDIR) -lboost_program_options -lboost_filesystem -lboost_system -lboost_serialization -lboost_mpi -lboost_date_time

LDFLAGS += -Wl,-rpath,$(PETSC_DIR)/lib
LDFLAGS += -L $(PETSC_DIR)/lib -lpetsc

LDFLAGS += -Wl,-rpath,$(NEXTSIMDIR)/lib
LDFLAGS += -L $(NEXTSIMDIR)/lib -lbamg
#LDFLAGS += -L $(NEXTSIMDIR)/lib -linterp
LDFLAGS += -L $(NEXTSIMDIR)/lib -lmapx
#LDFLAGS += -L $(NEXTSIMDIR)/lib -loasis

ifdef USE_ENSEMBLE 
	LDFLAGS += -L $(NEXTSIMDIR)/lib -lpseudo2D
	LDFLAGS += -lgfortran
	LDFLAGS += -L${FFTW_DIR}/lib -lfftw3
endif

ifdef USE_AEROBULK
		LDFLAGS += -L $(NEXTSIMDIR)/lib -laerobulk_nextsim
		LDFLAGS += -L /opt/local/aerobulk/lib -laerobulk
		LDFLAGS += -lgfortran
endif

ifneq (,$(strip $(filter DEBUG Debug debug PROFILE Profile profile,$(NEXTSIM_BUILD_TYPE))))
#ifeq ($(NEXTSIM_BUILD_TYPE),$(filter $(NEXTSIM_BUILD_TYPE),Debug debug))
	LDFLAGS += -Wl,-rpath,/opt/local/lib
ifneq (,$(strip $(filter PROFILE Profile profile,$(NEXTSIM_BUILD_TYPE))))
	LDFLAGS += -L /opt/local/lib -lprofiler
endif
endif

LDFLAGS += -L $(NEXTSIMDIR)/lib -lnextsim

OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
BINARYDIR=bin/

# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)

OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)

# Rules to always execute.
.PHONY: exec all clean mrproper libs cleanlibs mrproperlibs

# Default action.
exec: $(PROGNAME)

# Delete the object files.
clean:
	@echo " "
	$(RM) $(OBJS) $(DEPS)
	@echo " "

mrproper: clean
	$(RM) $(BINARYDIR)$(PROGNAME)
	@echo " "

# Rule for making the actual target
lines="=========="
Lines=$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)
$(PROGNAME): $(OBJS) #$(CCOBJS)
	@mkdir -p $(BINARYDIR)
	@echo " "
	@echo $(Lines)$(Lines)
	@echo "Creating executable: $(BINARYDIR)$(PROGNAME)"
	@echo $(Lines)$(Lines)
	@echo " "
	$(CXX) $(CXXFLAGS) -o $(BINARYDIR)$@ $^ $(LDFLAGS)
	@echo " "
	@echo $(Lines)$(Lines)
	@echo "Created executable: $(BINARYDIR)$(PROGNAME)"
	@echo $(Lines)$(Lines)
	@echo " "

# Rules for object files from cpp files
$(OBJECTDIR)%.o: %.cpp
	@mkdir -p $(OBJECTDIR) # $(dir $@)
	@echo " "
	$(CXX) -o $@ -c $< $(CXXFLAGS)

# Make dependancy rules
$(DEPENDDIR)%.d: %.cpp
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'

# The compilation depends on this Makefile.
$(OBJS): Makefile

# Make the required libraries
libs:
	cd ..; make

# Compile the libraries as well as the executable
all: libs exec

# Clean the required libraries
cleanlibs:
	cd ..; make clean

# Clean libs and current dir
cleanall: cleanlibs clean

# Clean recompile
fresh: cleanall all

# Fully clean the required libraries
mrproperlibs:
	cd ..; make mrproper

-include $(DEPS)
