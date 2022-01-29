KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	LIBRARYNAME=libnextsim.so.1.0
else
	LIBRARYNAME=libnextsim.dylib
endif

LIB=$(LIBRARYNAME)
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
LIBRARYDIR=$(NEXTSIMDIR)/lib/



CXXFLAGS += -std=c++14

# add g++ option flags
CXXFLAGS += -pedantic -ftemplate-depth-256 -Wno-inline -fpermissive

ifdef USE_OASIS
	CXXFLAGS += -DOASIS
	CXXFLAGS += -I$(NEXTSIMDIR)/modules/oasis/include
	LDFLAGS += -lgfortran
	LDFLAGS += -L$(NEXTSIMDIR)/lib -loasis
	CHAN = MPI1
	#LIBPSMILE = $(OASIS_DIR)/lib/libpsmile.${CHAN}.a $(OASIS_DIR)/lib/libmct.a $(OASIS_DIR)/lib/libmpeu.a $(OASIS_DIR)/lib/libscrip.a
endif

ifneq (,$(strip $(filter DEBUG Debug debug,$(NEXTSIM_BUILD_TYPE))))
	#ifeq ($(NEXTSIM_BUILD_TYPE),$(filter $(NEXTSIM_BUILD_TYPE),Debug debug))
	#CXXFLAGS := $(filter-out -O3 -pthread,$(CXXFLAGS))
	CXXFLAGS += -g -O0 -DNDEBUG
ifneq ($(KERNEL),Linux)
	CXXFLAGS += -Wl,-no_pie
endif
else
	CXXFLAGS += -O3 -pthread
endif

ifdef NEXTSIM_COMPILE_VERBOSE
	CXXFLAGS += -v
endif

# add include paths
CXXFLAGS += -I$(NEXTSIMDIR)/core/include -I/usr/local/include -I$(MPI_INC_DIR)/
CXXFLAGS += -I$(BOOST_INCDIR)/ -I$(GMSH_DIR)/include/gmsh
CXXFLAGS += -I$(NEXTSIMDIR)/contrib/mapx/include

ifeq ($(KERNEL),Linux)
	LDFLAGS += -shared -Wl,-soname,libnextsim.so.1
else
	ifeq ($(CXX),clang)
		CXXFLAGS += -stdlib=libc++
	endif

	LDFLAGS += -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-current_version,1.0,-install_name,$(LIBRARYDIR)$(LIBRARYNAME)

	ifeq ($(CXX),clang)
		LDFLAGS += -stdlib=libc++
	endif
endif

LDFLAGS += -ldl -lstdc++ #-lpthread #-lm -lrt

#boost
LDFLAGS += -Wl,-rpath,$(BOOST_LIBDIR)
LDFLAGS += -L$(BOOST_LIBDIR) -lboost_program_options -lboost_filesystem -lboost_system -lboost_serialization -lboost_mpi

# gmsh
LDFLAGS += -Wl,-rpath,$(GMSH_DIR)/lib
LDFLAGS += -L$(GMSH_DIR)/lib -lGmsh

# mapx
LDFLAGS += -L$(NEXTSIMDIR)/lib -lmapx


# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)

OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)

# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: $(LIB)

# Delete the object files.
clean:
	$(RM) $(OBJS) $(DEPS) # $(CCOBJS) $(CCDEPS)

mrproper: clean
	$(RM) $(LIBRARYDIR)$(LIBRARYNAME)

# Rule for making the actual target
ifeq ($(KERNEL),Linux)
$(LIBRARYNAME): $(OBJS) #$(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$@ $^ $(LDFLAGS)
	@ln -sf $@ $(LIBRARYDIR)libnextsim.so
	@ln -sf $@ $(LIBRARYDIR)libnextsim.so.1
else
$(LIBRARYNAME): $(OBJS) #$(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$@ $^ $(LDFLAGS)
	@ln -sf $@ $(LIBRARYDIR)libnextsim.dylib.0
	@ln -sf $@ $(LIBRARYDIR)libnextsim.dylib.1
endif

# Rules for object files from cpp files
$(OBJECTDIR)%.o: %.cpp
	@mkdir -p $(OBJECTDIR) # $(dir $@)
	$(CXX) -o $@ -c $< $(CXXFLAGS)

# Make dependancy rules
$(DEPENDDIR)%.d: %.cpp
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'


# The compilation depends on this Makefile.
$(OBJS): Makefile

# Include the dependency files
-include $(DEPS)
