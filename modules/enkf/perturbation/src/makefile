# ========================================================
# kernel type
KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	LIBRARYNAME=libpseudo2D.so.1.0
else
	LIBRARYNAME=libpseudo2D.dylib
endif

# ========================================================

# ========================================================
# output names/locations
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
BINARYDIR=bin/
LIBRARYDIR=$(NEXTSIMDIR)/lib/


CPPFLAGS=-DIARGC -DFFTW -DLAPACK
# ========================================================

# ========================================================
# C++ files
CXXSRCDIR=./
CXXSRC=$(wildcard $(CXXSRCDIR)*.cpp)
OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)


# setting the C++ standard according to the gcc compiler version (from gcc-5.2, the default is C++14)
ifeq ($(shell echo `$(CXX) -dumpversion | cut -f1-2 -d.` \>= 5.2 | sed -e 's/\.//g' | bc),1)
CXXFLAGS += -std=c++14
else
CXXFLAGS += -std=c++11
endif

# add g++ option flags
CXXFLAGS += -ftemplate-depth-256 -Wno-inline -v #-pedantic
# CXXFLAGS += -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread

# add include paths
CXXFLAGS += -I$(BOOST_INCDIR)/
CXXFLAGS += -I$(NEXTSIMDIR)/modules/enkf/perturbation/include
LNKFLAGS  = -lstdc++ -lgfortran

# ========================================================
# Fortran files
FCSRCDIR=./
FCSRC=$(wildcard $(FCSRCDIR)*.F90)
#FCOBJS=$(FCSRC:%.F90=$(OBJECTDIR)%.o)
FCOBJS= \
	$(OBJECTDIR)memory_stack.o \
	$(OBJECTDIR)m_set_random_seed.o \
	$(OBJECTDIR)sangoma_base.o \
	$(OBJECTDIR)mod_pseudo.o \
	$(OBJECTDIR)mod_random_forcing.o \
	$(OBJECTDIR)p_pseudo2D_fld.o
FCDEPS=$(FCSRC:%.F90=$(DEPENDDIR)%.d)
MODS=$(FCSRC:%.F90=$(FCSRCDIR)%.mod)

FFTW_INC   = $(FFTW_DIR)/include
FFTW_LIB   = $(FFTW_DIR)/lib

#FC = gfortran
#FC = ifort
#FFLAGS = -fPIC -O3 -r8
#FFLAGS += -I $(OPENMPI_INCLUDE_DIR)
FC = mpif90
F90FLG= -free
FFLAGS = -fPIC -O3 -fconvert=big-endian -fdefault-real-8
FFLAGS += -I$(FFTW_INC)/
# ========================================================

# ========================================================
# Libraries
CXXLIBS  += -L$(BOOST_LIBDIR) -L$(OPENMPI_LIB_DIR)/
CXXLIBS  += -L$(NEXTSIMDIR)/lib
CXXLIBS  += -lboost_program_options -lboost_filesystem -lboost_system
CXXLIBS  += -lboost_serialization -lboost_mpi -lboost_date_time
CXXLIBS  += -lmpi_cxx -lmpi -ldl -lstdc++ -lpthread
CXXLIBS  += -Wl,-rpath,$(OPENMPI_LIB_DIR)/

LDFLAGS  += -L$(FFTW_LIB) -lfftw3

# ========================================================
# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: lib

# Delete the object and dependancy files.
clean:
	$(RM) $(OBJS) $(FCOBJS) $(DEPS) $(FCDEPS) $(MODS)

# Delete the object and dependancy files, and the executable
mrproper: clean
	$(RM) $(BINARYDIR)$(PROG) $(LIBRARYDIR)*

# Delete the executable
execlean:
	$(RM) $(PROGOBJS) $(BINARYDIR)$(PROG)

## ================================================================
## Rules for .F90 files
$(FCOBJS): $(OBJECTDIR)%.o: $(FCSRCDIR)%.F90
	@mkdir -p $(OBJECTDIR)
	$(FC) $(CPPFLAGS) $(F90FLG) $(FFLAGS) $(LDFLAGS) -c $< -o $@

## Make dependency rules for .F90 files
#$(FCDEPS): $(DEPENDDIR)%.d : $(FCSRCDIR)%.F90
#	@mkdir -p $(DEPENDDIR)
#	@$(SHELL) -ec '$(FC) -M $(FFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'
# ================================================================


# ================================================================
# Rules for .cpp files
$(OBJS): $(OBJECTDIR)%.o: $(CXXSRCDIR)%.cpp
	@mkdir -p $(OBJECTDIR)
	$(CXX) -o $@ -c $< $(CXXFLAGS)
#	$(CXX) $(CXXFLAGS) $(CXXLIBS) $(LNKFLAGS) -c $< -o $@

# Make dependency rules for .cpp files
#$(DEPENDDIR)%.d: $(CXXSRCDIR)%.cpp
#	@mkdir -p $(DEPENDDIR)
#	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'
## ================================================================

# ========================================================
# make enkf perturbation library

LDFLAGSLIB = $(LDFLAGS)
ifeq ($(KERNEL),Linux)
	LDFLAGSLIB += -shared -Wl,-soname,libpseudo2D.so.1
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
	$(CXX) $(LNKFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGSLIB)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)/libpseudo2D.so
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)/libpseudo2D.so.1
else
lib: $(FCOBJS) $(OBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(LNKFLAGS) $(CXXFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGSLIB)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)/libpseudo2D.dylib.0
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)/libpseudo2D.dylib.1
endif
# ========================================================

# The compilation depends on this Makefile.
#$(OBJS): Makefile

#-include $(DEPS)
#-include $(FCDEPS)
