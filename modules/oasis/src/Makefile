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


CXX = mpicxx

CXXFLAGS += -std=c++14

# add g++ option flags
CXXFLAGS += -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread -fopenmp -v #-pedantic

# add include paths
CXXFLAGS += -I $(NEXTSIMDIR)/modules/oasis/include

# ========================================================


# ========================================================
# Fortran files
FCSRCDIR=.
FCSRC=$(wildcard $(FCSRCDIR)/*.F90)
FCOBJS=$(FCSRC:%.F90=$(OBJECTDIR)%.o)
FCDEPS=$(FCSRC:%.F90=$(DEPENDDIR)%.d)

FC = mpif90
FFLAGS = -std=f2008 -fopenmp -lstdc++ -fPIC
FFLAGS += -I $(NEXTSIMDIR)/modules/oasis/include
FFLAGS += -I $(OASIS_DIR)/build/lib/mct -I $(OASIS_DIR)/build/lib/psmile.MPI1
FFLAGS += -I $(NETCDF_FOR_DIR)/include
# ========================================================

# ========================================================
# Libraries
LDFLAGS += -fopenmp
LDFLAGS += -ldl -lstdc++ -lpthread
LDFLAGS += -L $(OASIS_DIR)/lib
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
CHAN = MPI1
LIBPSMILE = $(OASIS_DIR)/lib/libpsmile.${CHAN}.a $(OASIS_DIR)/lib/libmct.a $(OASIS_DIR)/lib/libmpeu.a $(OASIS_DIR)/lib/libscrip.a

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
	# I'm aving problems with the shared library code, so using a static library in the meantime
	# Compiler asks to recompile with -fPIC, but doing so doesn't help
	# $(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LIBPSMILE) $(LDFLAGSLIB)
	# @ln -sf $(LIBRARYNAME) $(LIBRARYDIR)liboasis.so
	# @ln -sf $(LIBRARYNAME) $(LIBRARYDIR)liboasis.so.1
	ar rcs $(LIBRARYDIR)/liboasis.a $^
else
lib: $(OBJS) $(FCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LIBPSMILE) $(LDFLAGSLIB)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)liboasis.dylib.0
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)liboasis.dylib.1
endif
# ========================================================

# The compilation depends on this Makefile.
$(OBJS): Makefile

-include $(DEPS)
-include $(FCDEPS)
