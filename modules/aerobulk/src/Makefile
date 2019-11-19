# ========================================================
# kernel type
KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	LIBRARYNAME=libaerobulk_nextsim.so.1.0
else
	LIBRARYNAME=libaerobulk_nextsim.dylib
endif

# ========================================================

# ========================================================
# output names/locations
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
LIBRARYDIR=$(NEXTSIMDIR)/lib
# ========================================================

# ========================================================
# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)

CXX = mpicxx

# setting the C++ standard according to the gcc compiler version (from gcc-5.2, the default is C++14)
ifeq ($(shell echo `$(CXX) -dumpversion | cut -f1-2 -d.` \>= 5.2 | sed -e 's/\.//g' | bc),1)
	CXXFLAGS += -std=c++14
else
	CXXFLAGS += -std=c++11
endif

# add g++ option flags
CXXFLAGS += -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread -fopenmp -v #-pedantic

# add include paths
CXXFLAGS += -I $(NEXTSIMDIR)/modules/aerobulk/include

# ========================================================


# ========================================================
# Fortran files
FCSRCDIR=.
FCSRC=$(wildcard $(FCSRCDIR)/*.f90)
FCOBJS=$(FCSRC:%.f90=$(OBJECTDIR)%.o)
FCDEPS=$(FCSRC:%.f90=$(DEPENDDIR)%.d)

FC = mpif90
FFLAGS = -std=f2008 -fopenmp -lstdc++ -fPIC
FFLAGS += -I $(AEROBULK_DIR)
FFLAGS += -I $(OPENMPI_INCLUDE_DIR)/
# ========================================================

# ========================================================
# Libraries
LDFLAGS += -L $(AEROBULK_DIR)/lib

# ========================================================


# ========================================================
# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: lib

# Delete the object and dependancy files.
clean:
	$(RM) $(OBJS) $(FCOBJS) $(DEPS) $(FCDEPS)

# Delete the object and dependancy files, and the executable
vclean: clean
	$(RM) $(LIBRARYDIR)/libaerobulk_nextsim.*

# ================================================================
# Rules for .F90 files
$(FCOBJS): $(OBJECTDIR)%.o : $(FCSRCDIR)/%.F90
	@mkdir -p $(OBJECTDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Make dependancy rules for .F90 files
$(FCDEPS): $(DEPENDDIR)%.d : $(FCSRCDIR)/%.F90
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(FC) -M $(FFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'
# ================================================================


# ================================================================
# Rules for .cpp files
$(OBJECTDIR)%.o: $(CXXSRCDIR)/%.cpp
	@mkdir -p $(OBJECTDIR)
	$(CXX) -o $@ -c $< $(CXXFLAGS)

# Make dependancy rules for .cpp files
$(DEPENDDIR)%.d: $(CXXSRCDIR)/%.cpp
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'
# ================================================================

# ========================================================
# make aerobulk library

LDFLAGSLIB = $(LDFLAGS)
ifeq ($(KERNEL),Linux)
	LDFLAGSLIB += -shared -Wl,-soname,libaerobulk_nextsim.so.1
else
	LDFLAGSLIB += -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-current_version,1.0,-install_name,$(LIBRARYDIR)/$(LIBRARYNAME)
endif

ifeq ($(KERNEL),Linux)
lib: $(OBJS) $(FCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)/$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)/$(LIBRARYNAME) $^ $(LDFLAGSLIB)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)/libaerobulk_nextsim.so
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)/libaerobulk_nextsim.so.1
else
lib: $(OBJS) $(FCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)/$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)/$(LIBRARYNAME) $^ $(LDFLAGSLIB)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)/libaerobulk_nextsim.dylib.0
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)/libaerobulk_nextsim.dylib.1
endif
# ========================================================

# The compilation depends on this Makefile.
$(OBJS): Makefile

-include $(DEPS)
-include $(FCDEPS)
