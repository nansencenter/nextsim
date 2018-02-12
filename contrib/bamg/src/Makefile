KERNEL := $(shell uname -s)
ifeq ($(KERNEL),Linux)
   LIBEXT=so.1.0
else
   LIBEXT=dylib
endif

ifndef NEXTSIM_TOOLS_COMP_LIBS
   LIBRARYDIR=$(NEXTSIMDIR)/lib/
else
   LIBRARYDIR=$(NEXTSIMTOOLS_ROOT_DIR)/lib/nextsim/
endif

ifndef BAMG_NO_OMP
   # using omp
   LIB0=libbamg
   OBJECTDIR=$(NEXTSIMDIR)/objs/
   DEPENDDIR=$(NEXTSIMDIR)/.deps/
else
   # not using omp
   LIB0=libbamg_no_omp
   OBJECTDIR=./objs/
   DEPENDDIR=./.deps/
endif

LIB=$(LIB0).$(LIBEXT)
CXX = g++

# setting the C++ standard according to the gcc compiler version (from gcc-5.2, the default is C++14)
ifeq ($(shell echo `$(CXX) -dumpversion | cut -f1-2 -d.` \>= 5.2 | sed -e 's/\.//g' | bc),1)
   CXXFLAGS += -std=c++14
else
   CXXFLAGS += -std=c++11
endif

CXXFLAGS += -pedantic -ftemplate-depth-256 \
		-Wno-inline -fPIC
ifndef BAMG_NO_PTHREAD
   CXXFLAGS += -D_MULTITHREADING_
endif
ifndef BAMG_NO_OMP
   CXXFLAGS += -DBAMG_NO_OMP
endif

ifneq (,$(strip $(filter DEBUG Debug debug,$(NEXTSIM_BUILD_TYPE))))
   # debugging
   CXXFLAGS += -g -O0 -DDEBUGGING
   ifneq ($(KERNEL),Linux)
      CXXFLAGS += -Wl,-no_pie
   endif
else
   # optimised
   CXXFLAGS += -O3
   ifndef BAMG_NO_PTHREAD
      CXXFLAGS += -pthread  
   endif
endif

ifneq ($(XSIMCOMP_VERBOSE),no)
   CXXFLAGS += -v
endif


ifeq ($(KERNEL),Linux)
   CXXFLAGS += -I $(NEXTSIMDIR)/contrib/bamg/include
   ifndef BAMG_NO_OMP
      # CXXFLAGS += -I /usr/lib/openmpi/include
      CXXFLAGS += -I $(OPENMPI_INCLUDE_DIR)
   endif
   LDFLAGS += -shared -Wl,-soname,$(LIB0).so.1
   LDFLAGS += -lm -ldl #-lhwloc
else
   ifeq ($(CXX),clang)
      CXXFLAGS += -stdlib=libc++
   endif
   
   CXXFLAGS += -I $(NEXTSIMDIR)/contrib/bamg/include -I /usr/local/include
   LDFLAGS += -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-current_version,1.0,-install_name,$(LIBRARYDIR)$(LIB)
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
$(LIB): $(OBJS) #$(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$@ $^ $(LDFLAGS)
	@ln -sf $@ $(LIBRARYDIR)$(LIB0).so
	@ln -sf $@ $(LIBRARYDIR)$(LIB0).so.1
else
$(LIB): $(OBJS) #$(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$@ $^ $(LDFLAGS)
	@ln -sf $@ $(LIBRARYDIR)$(LIB0).dylib.0
	@ln -sf $@ $(LIBRARYDIR)$(LIB0).dylib.1
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
