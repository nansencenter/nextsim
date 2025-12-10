#KERNEL := $(shell uname -s)

PROGNAME=nextsim.exec

CXXFLAGS += -std=c++14

# add g++ option flags
CXXFLAGS += -ftemplate-depth-256 -Wno-inline \
		-DHAVE_CONFIG_H -D_MULTITHREADING_
ifdef NEXTSIM_COMPILE_VERBOSE
	CXXFLAGS += -v
endif

# our libraries (include paths only, linking comes later)
CXXFLAGS += -I$(NEXTSIMDIR)/core/include
CXXFLAGS += -I$(NEXTSIMDIR)/contrib/bamg/include
CXXFLAGS += -I$(NEXTSIMDIR)/contrib/mapx/include

# our optional libraries
# -oasis
ifdef USE_OASIS
	CXXFLAGS += -DOASIS
	CXXFLAGS += -I$(NEXTSIMDIR)/modules/oasis/include
	CHAN = MPI1
endif

# -ensemble
ifdef USE_ENSEMBLE
	CXXFLAGS += -DENSEMBLE
	CXXFLAGS += -I$(NEXTSIMDIR)/modules/enkf/perturbation/include
endif

# external libraries (include paths)
# -boost
ifeq ($(BOOST_DIR),/usr)
	CXXFLAGS += -I$(BOOST_INCDIR)/ -I.
else
	CXXFLAGS += -isystem $(BOOST_INCDIR)/ -I.
endif

# -gmsh
CXXFLAGS += -isystem $(GMSH_DIR)/include/gmsh

# -netcdf
ifeq ($(NETCDF_CXX_DIR),)
	NETCDF_CXX_DIR=$(NETCDF_DIR)
endif
ifeq ($(NETCDF_FOR_DIR),)
	NETCDF_FOR_DIR=$(NETCDF_DIR)
endif

ifeq ($(NETCDF_DIR),/usr)
	CXXFLAGS += -I$(NETCDF_DIR)/include
	CXXFLAGS += -I$(NETCDF_CXX_DIR)/include
else
	CXXFLAGS += -isystem $(NETCDF_DIR)/include
	CXXFLAGS += -isystem $(NETCDF_CXX_DIR)/include
endif

# -mmg
ifdef USE_MMG
    CXXFLAGS += -DMMG
    CXXFLAGS += -I$(PARMMG_DIR)/include
endif

# -aerobulk
ifdef USE_AEROBULK
	CXXFLAGS += -isystem $(AEROBULK_DIR)/include
	CXXFLAGS += -DAEROBULK
endif

# debug build
ifneq (,$(strip $(filter DEBUG Debug debug PROFILE Profile profile,$(NEXTSIM_BUILD_TYPE))))
	CXXFLAGS := $(filter-out -O3 -pthread,$(CXXFLAGS))
	CXXFLAGS += -g -O0 -DNDEBUG
	LDFLAGS += -Wl,-rpath,/opt/local/lib

	# profiling build
	ifneq (,$(strip $(filter PROFILE Profile profile,$(NEXTSIM_BUILD_TYPE))))
		CXXFLAGS += -DWITHGPERFTOOLS
		LDFLAGS += -L/opt/local/lib -lprofiler
	endif
endif

# ========================================
# LINKING ORDER - CRITICAL FOR STATIC LIBS
# Most dependent libraries FIRST
# ========================================

# Our static libraries first (they depend on everything below)
LDFLAGS += -L$(NEXTSIMDIR)/lib -lnextsim -lbamg -lmapx

# Optional libraries
ifdef USE_OASIS
	LDFLAGS += -loasis
	LDFLAGS += -L$(OASIS_DIR)/lib -lpsmile.$(CHAN) -lmct -lmpeu -lscrip $(LD_EXTRA_OASIS)
endif

ifdef USE_ENSEMBLE
	LDFLAGS += -lpseudo2D
endif

# MMG and METIS
ifdef USE_MMG
    LDFLAGS += -L$(PARMMG_DIR)/lib -lmmg -lparmmg2d -lmetis
endif

# External libraries in dependency order
# Gmsh
LDFLAGS += -Wl,-rpath,$(GMSH_DIR)/lib
LDFLAGS += -L$(GMSH_DIR)/lib -lGmsh

# Boost
LDFLAGS += -Wl,-rpath,$(BOOST_LIBDIR)
LDFLAGS += -L$(BOOST_LIBDIR) -lboost_program_options -lboost_filesystem -lboost_system -lboost_serialization -lboost_mpi -lboost_date_time

# NetCDF
LDFLAGS += -Wl,-rpath,$(NETCDF_DIR)/lib -L$(NETCDF_DIR)/lib -lnetcdf
LDFLAGS += -Wl,-rpath,$(NETCDF_CXX_DIR)/lib -L$(NETCDF_CXX_DIR)/lib -lnetcdf_c++4
ifdef USE_OASIS
	LDFLAGS += -L$(NETCDF_FOR_DIR)/lib -lnetcdff
endif

# Aerobulk
ifdef USE_AEROBULK
	LDFLAGS += -L$(AEROBULK_DIR)/lib -laerobulk_cxx -laerobulk
	LDFLAGS += $(LD_EXTRA_AEROBULK)
endif

# System libraries last
LDFLAGS += -ldl -lstdc++

OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
BINARYDIR=bin/

# C++ files
CXXSRCDIR=.
CXXHDRDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
# We must exclude the version.hpp file from the list of header files because otherwise we get a circular dependency
CXXHDR=$(filter-out $(CXXHDRDIR)/version.hpp, $(wildcard $(CXXHDRDIR)/*.hpp))

OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)

# Rules to always execute.
.PHONY: exec clean mrproper all cleanall mrproperall

# Default action.
exec: $(PROGNAME)

# Create a header file with the git version
version.hpp: version.sh $(CXXSRC) $(CXXHDR)
	$(SHELL) -x $<

# Delete the object files.
clean:
	@echo " "
	$(RM) $(OBJS) $(DEPS) version.hpp
	@echo " "

mrproper: clean
	$(RM) $(BINARYDIR)$(PROGNAME)
	@echo " "

# Rule for making the actual target
lines="=========="
Lines=$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)
$(PROGNAME): $(OBJS)
	@mkdir -p $(BINARYDIR)
	@echo " "
	@echo $(Lines)$(Lines)
	@echo "Creating executable: $(BINARYDIR)$(PROGNAME)"
	@echo $(Lines)$(Lines)
	@echo " "
	$(CXX) $(CXXFLAGS) -o $(BINARYDIR)$@ $^ $(LDFLAGS)
	chmod u+x $(BINARYDIR)$(PROGNAME)
	@echo " "
	@echo $(Lines)$(Lines)
	@echo "Created executable: $(BINARYDIR)$(PROGNAME)"
	@echo $(Lines)$(Lines)
	@echo " "

# Rules for object files from cpp files
$(OBJECTDIR)%.o: %.cpp
	@mkdir -p $(OBJECTDIR)
	@echo " "
	$(CXX) -o $@ -c $< $(CXXFLAGS)

# Make dependancy rules
$(DEPENDDIR)%.d: %.cpp
	@mkdir -p $(DEPENDDIR)
	$(CXX) -MM -MG -MT $(OBJECTDIR)$*.o $(CXXFLAGS) $< -MF $@

# The compilation depends on this Makefile.
$(OBJS): Makefile

# Make everything
all:
	cd ..; $(MAKE) all

# Clean everything
cleanall:
	cd ..; $(MAKE) clean

# Properly clean everything
mrproperall:
	cd ..; $(MAKE) mrproper

# Properly clean & recompile
fresh:
	cd ..; $(MAKE) fresh

# Include dependencies only if the goal needs it
MAKEFILE_TARGETS_WITHOUT_INCLUDE := clean mrproper cleanall mrproperall fresh
ifeq ($(filter $(MAKECMDGOALS),$(MAKEFILE_TARGETS_WITHOUT_INCLUDE)),)
	-include $(DEPS)
endif
