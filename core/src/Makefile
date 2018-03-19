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


CXX = g++

# setting the C++ standard according to the gcc compiler version (from gcc-5.2, the default is C++14)
ifeq ($(shell echo `$(CXX) -dumpversion | cut -f1-2 -d.` \>= 5.2 | sed -e 's/\.//g' | bc),1)
	CXXFLAGS += -std=c++14
else
	CXXFLAGS += -std=c++11
endif

# add g++ option flags
CXXFLAGS += -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -fopenmp

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

ifneq ($(XSIMCOMP_VERBOSE),no)
	CXXFLAGS += -v
endif

# add include paths
CXXFLAGS += -I $(NEXTSIMDIR)/core/include -I /usr/local/include -I $(OPENMPI_INCLUDE_DIR)/
CXXFLAGS += -I $(PETSC_DIR)/include -I $(BOOST_INCDIR) -I $(GMSH_DIR)/include/gmsh
CXXFLAGS += -I $(NEXTSIMDIR)/contrib/mapx/include

ifeq ($(KERNEL),Linux)
	#CXXFLAGS += -std=c++0x -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -g -lm -pthread #-MMD -MP
	LDFLAGS += -shared -Wl,-soname,libnextsim.so.1
else
	ifeq ($(CXX),clang)
		CXXFLAGS += -stdlib=libc++
	endif

	LDFLAGS += -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-current_version,1.0,-install_name,$(LIBRARYDIR)$(LIBRARYNAME)

	#LDFLAGS += -Wl,-rpath,/usr/local/lib #-Wl,-rpath,$(OPENMPI_LIB_DIR)/
	#LDFLAGS += -L /usr/local/lib #-L $(OPENMPI_LIB_DIR)/ -lmpi_cxx -lmpi -ldl -lstdc++ -lpthread

	ifeq ($(CXX),clang)
		LDFLAGS += -stdlib=libc++
	endif
endif

LDFLAGS += -Wl,-rpath,$(OPENMPI_LIB_DIR)/
ifndef MACHINE_HEXAGON
	LDFLAGS += -L $(OPENMPI_LIB_DIR)/ -lmpi_cxx -lmpi -ldl -lstdc++ #-lpthread
else
	LDFLAGS += -L $(OPENMPI_LIB_DIR)/ -lmpichcxx -lmpich -ldl -lstdc++ #-lpthread #-lm -lrt
endif

#boost
LDFLAGS += -Wl,-rpath,$(BOOST_LIBDIR)
LDFLAGS += -L $(BOOST_LIBDIR) -lboost_program_options -lboost_filesystem -lboost_system -lboost_serialization -lboost_mpi

# petsc
LDFLAGS += -Wl,-rpath,$(PETSC_DIR)/lib
LDFLAGS += -L $(PETSC_DIR)/lib -lpetsc
ifeq ($(KERNEL),Linux)
	LDFLAGS += -L $(PETSC_DIR)/lib -lparmetis
	LDFLAGS += -L $(PETSC_DIR)/lib -lmetis
endif

# gmsh
LDFLAGS += -Wl,-rpath,$(GMSH_DIR)/lib
LDFLAGS += -L $(GMSH_DIR)/lib -lGmsh

# mapx
LDFLAGS += -L $(NEXTSIMDIR)/lib -lmapx


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
