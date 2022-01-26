# ========================================================
# kernel type
KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	LIBRARYNAME=libwim.so.1.0
else
	LIBRARYNAME=libwim.dylib
endif
# ========================================================

# ========================================================
# output names/locations
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
BINARYDIR=../run/bin/
LIBRARYDIR=$(NEXTSIMDIR)/lib/
# ========================================================

# ========================================================
# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)


CXXFLAGS += -std=c++14

# add g++ option flags
CXXFLAGS += -ftemplate-depth-256 -Wno-inline
ifneq ($(XSIMCOMP_VERBOSE),no)
	CXXFLAGS += -v
	CFLAGS   += -v
endif
CXXFLAGS += -Wfatal-errors
# CXXFLAGS += -pedantic

ifneq (,$(strip $(filter DEBUG Debug debug,$(NEXTSIM_BUILD_TYPE))))
	CXXFLAGS += -g -O0 -DDEBUGGING
ifneq ($(KERNEL),Linux)
	CXXFLAGS += -Wl,-no_pie
endif
endif

# add include paths
# CXXFLAGS += -I.
CXXFLAGS += -I$(NEXTSIMDIR)/modules/wim/include

CXXFLAGS += -I$(BOOST_INCDIR)
CXXFLAGS += -I$(MPI_INC_DIR)
CXXFLAGS += -I$(NEXTSIMDIR)/contrib/mapx/include
CXXFLAGS += -I$(NEXTSIMDIR)/contrib/bamg/include
CXXFLAGS += -I$(GMSH_DIR)/include/gmsh
CXXFLAGS += -I$(NEXTSIMDIR)/core/include
# ========================================================


# ========================================================
# C files
CCSRCDIR=.
CCSRC=$(wildcard $(CCSRCDIR)/*.c)
CCOBJS=$(CCSRC:%.c=$(OBJECTDIR)%.o)
CCDEPS=$(CCSRC:%.c=$(DEPENDDIR)%.d)


CFLAGS += -Wno-inline # -MMD -MP -lm #-Wall -pedantic
CFLAGS += -I$(NEXTSIMDIR)/modules/wim/include
# ========================================================

# ========================================================
# Libraries
LDFLAGS += -fopenmp
LDFLAGS += -ldl -lstdc++ -lpthread
LDFLAGS += -Wl,-rpath,$(BOOST_LIBDIR)
LDFLAGS += -L$(BOOST_LIBDIR) -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi -lboost_date_time
LDFLAGS += -L$(NEXTSIMDIR)/lib -lmapx -lbamg -lnextsim
# ========================================================


# ========================================================
# Rules to always execute.
.PHONY: default all lib clean mrproper single

# Default action.
default: lib

# make lib and stand-alone executable
all: lib single

# Delete the object and dependancy files.
clean:
	@echo " "
	$(RM) $(OBJS) $(CCOBJS) $(DEPS) $(CCDEPS) $(PROGOBJS)
	@echo " "

fresh: clean all

# Delete the object and dependancy files, and the executable
mrproper: clean
	$(RM) $(BINARYDIR)$(PROG) $(LIBRARYDIR)/*
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
# Programs
PROGDIR=$(CXXSRCDIR)/programs
PROG=WIM2d_single_call_cpp.exec

# Rules to compile all programs
ALLPROGS=$(wildcard $(PROGDIR)/*.cpp)
PROGOBJS=$(ALLPROGS:$(PROGDIR)/%.cpp=$(OBJECTDIR)/%.o)
$(PROGOBJS): $(OBJS) $(CCOBJS)
$(PROGOBJS): $(OBJECTDIR)%.o : $(PROGDIR)/%.cpp
	@mkdir -p $(OBJECTDIR)
	$(CXX) -o $@ -c $< $(CXXFLAGS)
# ================================================================


# ========================================================
# make executable 1
lines="=========="
Lines=$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)$(lines)
single: $(PROGOBJS) $(OBJS) $(CCOBJS)
	@mkdir -p $(BINARYDIR)
	@echo " "
	@echo $(Lines)$(Lines)
	@echo "Creating executable: $(BINARYDIR)$(PROG)"
	@echo $(Lines)$(Lines)
	@echo " "
	$(CXX) $(CXXFLAGS) -o $(BINARYDIR)$(PROG) $^ $(LDFLAGS)
	@echo " "
	@echo $(Lines)$(Lines)
	@echo "Created executable: $(BINARYDIR)$(PROG)"
	@echo $(Lines)$(Lines)
	@echo " "
# ========================================================


# ========================================================
# make wim library

LDFLAGSLIB = $(LDFLAGS)
ifeq ($(KERNEL),Linux)
	LDFLAGSLIB += -shared -Wl,-soname,libwim.so.1
	#LDFLAGSLIB += -lm -ldl #-lhwloc
else
	LDFLAGSLIB += -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-current_version,1.0,-install_name,$(LIBRARYDIR)$(LIBRARYNAME)
	#LDFLAGSLIB += -ldl -lstdc++
endif

LOBJS := $(filter-out $(OBJECTDIR)/main.o, $(OBJS))

ifeq ($(KERNEL),Linux)
lib: $(OBJS) $(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo " "
	@echo $(Lines)$(Lines)
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo $(Lines)$(Lines)
	@echo " "
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGSLIB)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libwim.so
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libwim.so.1
	@echo " "
	@echo $(Lines)$(Lines)
	@echo Created shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo $(Lines)$(Lines)
	@echo " "
else
lib: $(OBJS) $(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo " "
	@echo $(Lines)$(Lines)
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo $(Lines)$(Lines)
	@echo " "
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGSLIB)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libwim.dylib.0
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libwim.dylib.1
	@echo " "
	@echo $(Lines)$(Lines)
	@echo Created shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo $(Lines)$(Lines)
	@echo " "
endif
# ========================================================

# The compilation depends on this Makefile.
$(OBJS): Makefile

-include $(DEPS)
-include $(CCDEPS)
