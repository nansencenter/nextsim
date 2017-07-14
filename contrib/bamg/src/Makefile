KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	LIBRARYNAME=libbamg.so.1.0
else
	LIBRARYNAME=libbamg.dylib
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

CXXFLAGS += -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread -v -D_MULTITHREADING_

ifneq (,$(strip $(filter DEBUG Debug debug,$(NEXTSIM_BUILD_TYPE))))
	#ifeq ($(NEXTSIM_BUILD_TYPE),$(filter $(NEXTSIM_BUILD_TYPE),Debug debug))
	CXXFLAGS := $(filter-out -O3 -pthread,$(CXXFLAGS))
	CXXFLAGS += -g -O0 -Wl,-no_pie -DNDEBUG
endif


ifeq ($(KERNEL),Linux)
	CXXFLAGS += -I $(NEXTSIMDIR)/contrib/bamg/include -I /usr/lib/openmpi/include
	LDFLAGS += -shared -Wl,-soname,libbamg.so.1
	LDFLAGS += -lm -ldl #-lhwloc
else

	ifeq ($(CXX),clang)
	CXXFLAGS += -stdlib=libc++
	endif

	CXXFLAGS += -I $(NEXTSIMDIR)/contrib/bamg/include -I /usr/local/include
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
	@ln -sf $@ $(LIBRARYDIR)libbamg.so
	@ln -sf $@ $(LIBRARYDIR)libbamg.so.1
else
$(LIBRARYNAME): $(OBJS) #$(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$@ $^ $(LDFLAGS)
	@ln -sf $@ $(LIBRARYDIR)libbamg.dylib.0
	@ln -sf $@ $(LIBRARYDIR)libbamg.dylib.1
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
