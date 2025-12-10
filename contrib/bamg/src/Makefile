KERNEL := $(shell uname -s)

LIBRARYNAME=libbamg.a

OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
LIBRARYDIR=$(NEXTSIMDIR)/lib/

CXXFLAGS := $(filter-out -fopenmp -pthread,$(CXXFLAGS))

CXXFLAGS += -std=c++14

CXXFLAGS += -pedantic -ftemplate-depth-256 -Wno-inline -D_MULTITHREADING_

ifneq (,$(strip $(filter DEBUG Debug debug,$(NEXTSIM_BUILD_TYPE))))
	CXXFLAGS += -g -O0 -DNDEBUG
ifneq ($(KERNEL),Linux)
	CXXFLAGS += -Wl,-no_pie
endif
endif

ifdef NEXTSIM_COMPILE_VERBOSE
	CXXFLAGS += -v
endif

ifeq ($(KERNEL),Linux)
	CXXFLAGS += -I$(NEXTSIMDIR)/contrib/bamg/include -I$(MPI_INC_DIR)
else
	ifeq ($(CXX),clang)
		CXXFLAGS += -stdlib=libc++
	endif
	CXXFLAGS += -I$(NEXTSIMDIR)/contrib/bamg/include -I/usr/local/include
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
$(LIBRARYNAME): $(OBJS) #$(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "=================================================="
	@echo Creating static library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "=================================================="
	$(AR) rcs $(LIBRARYDIR)$@ $^

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
