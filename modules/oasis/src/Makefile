# ========================================================
# kernel type
KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	LIBRARYNAME=liboasis.so.1.0
else
	LIBRARYNAME=liboasis.dylib
endif

# ========================================================

# ========================================================
# output names/locations
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
BINARYDIR=bin/
LIBRARYDIR=$(NEXTSIMDIR)/lib/
# ========================================================

# ========================================================
# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)


CXX = g++
CXXFLAGS += -std=c++11 -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread -fopenmp -v #-pedantic

CXXFLAGS += -I $(OPENMPI_INCLUDE_DIR)/
CXXFLAGS += -I $(NEXTSIMDIR)/modules/oasis/include

# ========================================================


# ========================================================
# Fortran files
FCSRCDIR=.
FCSRC=$(wildcard $(FCSRCDIR)/*.F90)
FCOBJS=$(FCSRC:%.F90=$(OBJECTDIR)%.o)
FCDEPS=$(FCSRC:%.F90=$(DEPENDDIR)%.d)

FC = gfortran
FFLAGS = -std=f2008 -fopenmp -lstdc++ -fPIC
FFLAGS += -I $(NEXTSIMDIR)/modules/oasis/include
FFLAGS += -I $(OASIS_DIR)/build/lib/mct -I $(OASIS_DIR)/build/lib/psmile.MPI1
FFLAGS += -I $(NETCDF_FOR_DIR)/include
# ========================================================

# ========================================================
# Libraries
LDFLAGS += -fopenmp
LDFLAGS += -Wl,-rpath,$(OPENMPI_LIB_DIR)/
LDFLAGS += -L $(OPENMPI_LIB_DIR)/ -lmpi_cxx -lmpi -ldl -lstdc++ -lpthread
#LDFLAGS += -Wl,-rpath,$(OASIS_DIR)/lib
#LDFLAGS += -L /Users/abdama/COUPLE/oasis3-mct/build_darwin/lib -lmct -lmpeu -lpsmile.MPI1 -lscrip

LDFLAGS += -Wl,-rpath,$(NETCDF_FOR_DIR)/lib
LDFLAGS += -L $(NETCDF_FOR_DIR)/lib -lnetcdff
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
	$(RM) $(BINARYDIR)$(PROG) $(LIBRARYDIR)/*

# Delete the executable
execlean:
	$(RM) $(PROGOBJS) $(BINARYDIR)$(PROG)

# ========================================================


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
# make oasis library

LDFLAGSLIB = $(LDFLAGS)
ifeq ($(KERNEL),Linux)
	LDFLAGSLIB += -shared -Wl,-soname,liboasis.so.1
	#LDFLAGSLIB += -lm -ldl #-lhwloc
else
	LDFLAGSLIB += -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-current_version,1.0,-install_name,$(LIBRARYDIR)$(LIBRARYNAME)
	#LDFLAGSLIB += -ldl -lstdc++
endif

ifeq ($(KERNEL),Linux)
lib: $(OBJS) $(FCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGSLIB)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)liboasis.so
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)liboasis.so.1
else
lib: $(OBJS) $(FCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGSLIB)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)liboasis.dylib.0
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)liboasis.dylib.1
endif
# ========================================================

# The compilation depends on this Makefile.
$(OBJS): Makefile

-include $(DEPS)
-include $(FCDEPS)
