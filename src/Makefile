KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	LIBRARYNAME=libnextsim.so.1.0
else
	LIBRARYNAME=libnextsim.dylib
endif

NEXTSIMDIR=..

CXX = g++
CXXFLAGS += -std=c++0x -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -g -pthread -v

ifeq ($(KERNEL),Linux)
	#CXX = g++-4.8
	#CXXFLAGS += -std=c++0x -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -g -lm -pthread -v #-MMD -MP
	CXXFLAGS += -I $(NEXTSIMDIR)/include -I /usr/lib/openmpi/include -I $(PETSC_DIR)/include
	LDFLAGS += -shared -Wl,-soname,libnextsim.so.1 -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi -lmpi_cxx -lmpi -lm -ldl -lhwloc
	LDFLAGS += -Wl,-rpath,$(PETSC_DIR)/lib
	LDFLAGS += -L $(PETSC_DIR)/lib -lpetsc
else
	CXXFLAGS += -I $(NEXTSIMDIR)/include -I /usr/lib/openmpi/include -I /opt/local/include/openmpi-mp -I $(PETSC_DIR)/include -I $(BOOST_DIR)/include
	LDFLAGS += -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-current_version,1.0,-install_name,libnextsim.dylib
	LDFLAGS += -Wl,-rpath,/opt/local/lib/openmpi-mp
	LDFLAGS += -L /opt/local/lib/openmpi-mp -lmpi_cxx -lmpi -ldl -lstdc++
	LDFLAGS += -Wl,-rpath,$(BOOST_DIR)/lib
	LDFLAGS += -L $(BOOST_DIR)/lib -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi
	LDFLAGS += -Wl,-rpath,$(PETSC_DIR)/lib
	LDFLAGS += -L $(PETSC_DIR)/lib -lpetsc
endif


LIB=$(LIBRARYNAME)
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
LIBRARYDIR=$(NEXTSIMDIR)/lib/

# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
# OBJS=$(CXXSRC:.cpp=.o)
# DEPS=$(CXXSRC:.cpp=.d)

# OBJS=$(join $(addsuffix $(OBJECTDIR), $(dir $(CXXSRC))), $(notdir $(CXXSRC:.cpp=.o)))
# DEPS=$(join $(addsuffix $(DEPENDDIR), $(dir $(CXXSRC))), $(notdir $(CXXSRC:.cpp=.d)))

OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)


# C files
# CCSRCDIR=src
# CCSRC=$(wildcard $(CCSRCDIR)/*.c)
# CCOBJS=$(CCSRC:.c=.o)
# CCDEPS=$(CCSRC:.c=.d)

# CCOBJS:=$(addprefix $(OBJECTDIR), $(CCOBJS))
# CCDEPS:=$(addprefix $(OBJECTDIR), $(CCDEPS))



# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: $(LIB)

# Delete the object files.
clean:
	$(RM) $(OBJS) $(DEPS) # $(CCOBJS) $(CCDEPS)

mrproper: clean
	$(RM) $(LIBRARYDIR)$(LIBRARYNAME)


# default: $(CCOBJS)
# 	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# $(OBJECTDIR)%.o: %.c
# 	$(CC) -o $@ -c $< $(CFLAGS)

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
