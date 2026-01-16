KERNEL := $(shell uname -s)

LIBRARYNAME=libnextsim.a
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
LIBRARYDIR=$(NEXTSIMDIR)/lib/


CXXFLAGS += -std=c++14

# add g++ option flags
CXXFLAGS += -pedantic -ftemplate-depth-256 -Wno-inline -fpermissive

ifdef USE_OASIS
	CXXFLAGS += -DOASIS
	CXXFLAGS += -I$(NEXTSIMDIR)/modules/oasis/include
	CHAN = MPI1
endif

ifneq (,$(strip $(filter DEBUG Debug debug,$(NEXTSIM_BUILD_TYPE))))
	CXXFLAGS += -g -O0 -DNDEBUG
ifneq ($(KERNEL),Linux)
	CXXFLAGS += -Wl,-no_pie
endif
endif

ifdef NEXTSIM_COMPILE_VERBOSE
	CXXFLAGS += -v
endif

# add include paths
CXXFLAGS += -I$(NEXTSIMDIR)/core/include -I/usr/local/include -I$(MPI_INC_DIR)/
CXXFLAGS += -I$(BOOST_INCDIR)/ -I$(GMSH_DIR)/include/gmsh
CXXFLAGS += -I$(NEXTSIMDIR)/contrib/mapx/include

ifneq ($(KERNEL),Linux)
	ifeq ($(CXX),clang)
		CXXFLAGS += -stdlib=libc++
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
all: $(LIBRARYNAME)

# Delete the object files.
clean:
	$(RM) $(OBJS) $(DEPS) # $(CCOBJS) $(CCDEPS)

mrproper: clean
	$(RM) $(LIBRARYDIR)$(LIBRARYNAME)

# Rule for making the actual target
$(LIBRARYNAME): $(OBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating static library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(AR) rcs -o $(LIBRARYDIR)$@ $^

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
