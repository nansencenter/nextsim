# ========================================================
# kernel type
KERNEL := $(shell uname -s)

# main output
LIBRARYNAME=bamg.so
# ========================================================

# ========================================================
# output names/locations
NEWSRCDIR=./pysrc/
OBJECTDIR=./objs/
DEPENDDIR=./.deps/
LIBRARYDIR=$(NEXTSIMDIR)/python/lib/
# ========================================================

# ========================================================
# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)


CXX = gcc
#CXXFLAGS += -std=c++0x -std=c++11

# setting the C++ standard according to the gcc compiler version (from gcc-5.2, the default is C++14)
ifeq ($(shell echo `$(CXX) -dumpversion | cut -f1-2 -d.` \>= 5.2 | sed -e 's/\.//g' | bc),1)
	CXXFLAGS += -std=c++14
else
	CXXFLAGS += -std=c++11
endif

# add g++ option flags
CXXFLAGS += -Wfatal-errors
CXXFLAGS += -pthread -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -fPIC
# CXXFLAGS += -pedantic
ifneq ($(XSIMCOMP_VERBOSE),no)
	CXXFLAGS += -v
	CFLAGS   += -v
endif


# add include paths

CXXFLAGS += -I .
CXXFLAGS += -I ../include
# CXXFLAGS += -I $(BOOST_DIR)/include
# CXXFLAGS += -I $(PETSC_DIR)/include
# CXXFLAGS += -I $(GMSH_DIR)/include/gmsh
# CXXFLAGS += -I $(OPENMPI_INCLUDE_DIR)
# CXXFLAGS += -I $(NEXTSIMDIR)/contrib/mapx/include
CXXFLAGS += -isystem $(NEXTSIMDIR)/contrib/bamg/include # turn off annoying bamg warnings during compilation
# CXXFLAGS += -I $(NEXTSIMDIR)/core/include
# CXXFLAGS += -I $(NEXTSIMDIR)/modules/wim/include
CXXFLAGS += -I/Home/timill/python/miniconda2/include/python2.7
# ========================================================

# ========================================================
# pyx files (cython)
PYXSRCDIR=.
PYXSRC=$(wildcard $(PYXSRCDIR)/*.pyx)
PYXSRC2=$(PYXSRC:%.pyx=$(NEWSRCDIR)%.cpp)
PYXOBJS=$(PYXSRC:%.pyx=$(OBJECTDIR)%.o)
PYXDEPS=$(PYXSRC:%.pyx=$(DEPENDDIR)%.d)


PYXC = cython
PYXFLAGS = --cplus
#PYXFLAGS += -Wno-inline -fPIC -O3 # -MMD -MP -lm #-Wall -pedantic
#PYXFLAGS += -I $(NEXTSIMDIR)/modules/wim/include
PYXFLAGS += -I .
PYXFLAGS += -I ../include
# CFLAGS += -I $(CCSRCDIR)
# ========================================================

# ========================================================
# C files
CCSRCDIR=.
CCSRC=$(wildcard $(CCSRCDIR)/*.c)
CCOBJS=$(CCSRC:%.c=$(OBJECTDIR)%.o)
CCDEPS=$(CCSRC:%.c=$(DEPENDDIR)%.d)


CC = gcc
CFLAGS = $(CXXFLAGS)
CFLAGS += -I .
CFLAGS += -I ../include
CFLAGS += -I $(CCSRCDIR)
CFLAGS += -Wstrict-prototypes
# ========================================================

# ========================================================
# Libraries
LC = g++
LDFLAGS = -pthread -shared -L/Home/timill/python/miniconda2/lib -Wl,-rpath=/Home/timill/python/miniconda2/lib,--no-as-needed
LDFLAGS += -L/Home/timill/python/miniconda2/lib -lpython2.7
LDFLAGS += -L$(NEXTSIMDIR)/lib -lbamg
# ========================================================


# ========================================================
# Rules to always execute.
.PHONY: default pyx lib clean

# Default action.
default: lib

# Delete the object and dependancy files.
clean:
	@echo " "
	$(RM) $(OBJS) $(CCOBJS) $(DEPS) $(CCDEPS) $(PYXOBJS) 
	@echo " "

fresh: clean lib

# Delete the object and dependancy files, and the library
vclean: clean
	$(RM) $(PYXSRC2)
# ========================================================


# ================================================================
# Rules for .c files
$(CCOBJS): $(OBJECTDIR)%.o : $(CCSRCDIR)/%.c
	@echo " "
	@mkdir -p $(OBJECTDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Make dependancy rules for .c files
$(CCDEPS): $(DEPENDDIR)%.d : $(CCSRCDIR)/%.c
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CC) -M $(CFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'
# ================================================================


# ================================================================
# Rules for .cpp files
$(OBJS): $(OBJECTDIR)%.o: $(CXXSRCDIR)/%.cpp
	@echo " "
	@mkdir -p $(OBJECTDIR)
	$(CXX) -o $@ -c $< $(CXXFLAGS)

# Make dependancy rules for .cpp files
$(DEPENDDIR)%.d: $(CXXSRCDIR)/%.cpp
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'
# ================================================================


# ================================================================
# Rules for .pyx files

# 1. make .cpp from .pyx
$(PYXSRC2): $(NEWSRCDIR)%.cpp: $(PYXSRCDIR)/%.pyx
	@echo " "
	@mkdir -p $(NEWSRCDIR)
	$(PYXC) -o $@ $< $(PYXFLAGS)

# 2. make .o from .cpp
# - depends on the new .cpp files made by cython
$(PYXOBJS): $(PYXSRC2)

# - compile the new .cpp files made by cython
$(PYXOBJS): $(OBJECTDIR)%.o: $(NEWSRCDIR)%.cpp
	@echo " "
	@mkdir -p $(OBJECTDIR)
	$(CXX) -o $@ -c $< $(CXXFLAGS)

# # Make dependancy rules for .pyx files
# $(DEPENDDIR)%.d: $(PYXSRCDIR)/%.pyx
# @mkdir -p $(DEPENDDIR)
# @$(SHELL) -ec '$(PYXC) -M $(PYXFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'
# ================================================================

default: lib

# compile the cython files
pyx: $(PYXSRC2) $(PYXOBJS)

# compile the wrapper files
wrappers: $(OBJS) $(CCOBJS)

# ========================================================
# make library
lines="=========="
Lines=$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)

lib: $(OBJS) $(CCOBJS) $(PYXOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo " "
	@echo $(Lines)$(Lines)
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo $(Lines)$(Lines)
	@echo " "
	$(LC) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGS)
	@echo " "
	@echo $(Lines)$(Lines)
	@echo Created shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo $(Lines)$(Lines)
	@echo " "
# ========================================================

# The compilation depends on this Makefile.
$(OBJS): Makefile
$(CCOBJS): Makefile
$(PYXOBJS): Makefile

-include $(DEPS)
-include $(CCDEPS)
