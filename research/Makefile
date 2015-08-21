KERNEL := $(shell uname -s)

PROGNAME=nextsim.exec

NEXTSIMDIR=..

CXX = g++
CXXFLAGS += -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -g -pthread -v -DHAVE_CONFIG_H

ifeq ($(KERNEL),Linux)
	#CXX = g++-4.8
	#CXXFLAGS += -std=c++0x -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -g -lm -pthread -v #-MMD -MP -lm -pthread -v
	CXXFLAGS += -I. -I $(NEXTSIMDIR)/include -I /usr/lib/openmpi/include -I $(PETSC_DIR)/include #-I headers
	LDFLAGS += -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi -lmpi_cxx -lmpi -lm -ldl -lhwloc
	LDFLAGS += -Wl,-rpath,$(PETSC_DIR)/lib -Wl,-rpath,$(NEXTSIMDIR)/lib
	LDFLAGS += -L $(PETSC_DIR)/lib -lpetsc -L$(NEXTSIMDIR)/lib -lnextsim
else

	ifeq ($(CXX),clang)
	CXXFLAGS += -stdlib=libc++
	endif

	CXXFLAGS += -I $(NEXTSIMDIR)/include -I $(NEXTSIMDIR)/bamg/include -I $(NEXTSIMDIR)/mapx/include
	CXXFLAGS += -I /usr/local/include -I /opt/local/include/openmpi-mp
	CXXFLAGS +=	-I $(PETSC_DIR)/include -I $(BOOST_DIR)/include #-I /Users/abdama/Downloads/ISSM
	#CXXFLAGS +=	-I $(PETSC_DIR)/include -I $(BOOST_DIR)/include -I /Users/abdama/ISSM
	CXXFLAGS +=	-I $(NETCDF_DIR)/include -I $(GMSH_DIR)/include/gmsh

	LDFLAGS += -Wl,-rpath,/usr/local/lib -Wl,-rpath,/opt/local/lib/openmpi-mp
	LDFLAGS += -L /usr/local/lib -L /opt/local/lib/openmpi-mp -lmpi_cxx -lmpi -ldl -lstdc++
	#LDFLAGS += -Wl,-rpath,/opt/local/issm/lib -L /opt/local/issm/lib -lISSMOverload -lISSMCore # -lISSMModules -lISSMOverload
	LDFLAGS += -Wl,-rpath,/opt/local/netcdf-cxx/lib -L /opt/local/netcdf-cxx/lib -lnetcdf_c++4

	ifeq ($(CXX),clang)
	LDFLAGS += -stdlib=libc++
	endif

	LDFLAGS += -Wl,-rpath,$(BOOST_DIR)/lib
	LDFLAGS += -L $(BOOST_DIR)/lib -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi
	LDFLAGS += -Wl,-rpath,$(PETSC_DIR)/lib -Wl,-rpath,$(NEXTSIMDIR)/lib
	LDFLAGS += -L $(PETSC_DIR)/lib -lpetsc -L $(NEXTSIMDIR)/lib -lbamg -L $(NEXTSIMDIR)/lib -lmapx
	#LDFLAGS += -lm -pthread -L /opt/local/mapx/lib -lmapx
	LDFLAGS += -L $(NEXTSIMDIR)/lib -lnextsim

	#LDFLAGS += -Wl,-rpath,$(GMSH_DIR)/lib
	#LDFLAGS += -L $(GMSH_DIR)/lib -lGmsh

endif

#-shared -Wl,-soname,libnextsim.so.1 -o libnextsim.so.1.0.0

EXEC=$(PROGNAME)
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
BINARYDIR=bin/

# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
#OBJS=$(CXXSRC:.cpp=.o)
#DEPS=$(CXXSRC:.cpp=.d)

# OBJS:=$(addprefix $(OBJECTDIR), $(OBJS))
# DEPS:=$(addprefix $(OBJECTDIR), $(DEPS))

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
all: $(EXEC)

# Delete the object files.
clean:
	$(RM) $(OBJS) $(DEPS) # $(CCOBJS) $(CCDEPS)

mrproper: clean
	$(RM) $(BINARYDIR)$(PROGNAME)


# default: $(CCOBJS)
# 	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# $(OBJECTDIR)%.o: %.c
# 	$(CC) -o $@ -c $< $(CFLAGS)

# Rule for making the actual target
$(PROGNAME): $(OBJS) #$(CCOBJS)
	@mkdir -p $(BINARYDIR)
	$(CXX) $(CXXFLAGS) -o $(BINARYDIR)$@ $^ $(LDFLAGS)


# $(OBJECTDIR)%.o: %.cpp
# 	@mkdir -p $(OBJECTDIR)
# 	$(CXX) -o $@ -c $< $(CXXFLAGS)

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

-include $(DEPS)
