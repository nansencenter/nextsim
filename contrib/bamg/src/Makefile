KERNEL := $(shell uname -s)
ifeq ($(KERNEL),Linux)
   LIBEXT=so.1.0
else
   LIBEXT=dylib
endif

# default stuff
LIBRARYDIR=$(NEXTSIMDIR)/lib/
LIB0=libbamg
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/

# compile tools using omp
tools: LIBRARYDIR=$(NEXTSIMTOOLS_ROOT_DIR)/lib/nextsim/

# compile tools without using omp
tools-no-omp: LIBRARYDIR=$(NEXTSIMTOOLS_ROOT_DIR)/lib/nextsim/
tools-no-omp: LIB0=libbamg_no_omp
tools-no-omp: CXXFLAGS+= -DBAMG_NO_OMP

LIBRARYNAME=$(LIB0).$(LIBEXT)
CXX = g++

# setting the C++ standard according to the gcc compiler version (from gcc-5.2, the default is C++14)
ifeq ($(shell echo `$(CXX) -dumpversion | cut -f1-2 -d.` \>= 5.2 | sed -e 's/\.//g' | bc),1)
   CXXFLAGS += -std=c++14
else
   CXXFLAGS += -std=c++11
endif

CXXFLAGS += -pedantic -ftemplate-depth-256 \
		-Wno-inline -fPIC \
		-D_MULTITHREADING_

ifneq (,$(strip $(filter DEBUG Debug debug,$(NEXTSIM_BUILD_TYPE))))
   # debugging
   CXXFLAGS += -g -O0 -DDEBUGGING
   ifneq ($(KERNEL),Linux)
      CXXFLAGS += -Wl,-no_pie
   endif
else
   # optimised
   CXXFLAGS += -O3
   CXXFLAGS += -pthread  
endif

ifneq ($(XSIMCOMP_VERBOSE),no)
   CXXFLAGS += -v
endif


CXXFLAGS += -I $(NEXTSIMDIR)/contrib/bamg/include
ifeq ($(KERNEL),Linux)
   CXXFLAGS += -I $(OPENMPI_INCLUDE_DIR)
   LDFLAGS += -shared -Wl,-soname,$(LIB0).so.1
   LDFLAGS += -lm -ldl #-lhwloc
else
   ifeq ($(CXX),clang)
      CXXFLAGS += -stdlib=libc++
   endif
   
   CXXFLAGS += -I /usr/local/include
   LDFLAGS += -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-current_version,1.0,-install_name,$(LIBRARYDIR)$(LIBRARYNAME)
   LDFLAGS += -ldl -lstdc++
   
   ifeq ($(CXX),clang)
      LDFLAGS += -stdlib=libc++
   endif
endif

# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)

OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)

# Rules to always execute.
.PHONY: default lib clean mrproper tools tools-no-omp

# Default action.
default: lib

tools: clean lib
tools-no-omp: clean lib

# Delete the object files.
clean:
	$(RM) $(OBJS) $(DEPS)

mrproper: clean
	$(RM) $(LIBRARYDIR)$(LIBRARYNAME)

# Rule for making the actual target
ifeq ($(KERNEL),Linux)
lib: $(OBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGS)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)$(LIB0).so
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)$(LIB0).so.1
else
lib: $(OBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGS)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)$(LIB0).dylib.0
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)$(LIB0).dylib.1
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
