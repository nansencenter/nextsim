KERNEL := $(shell uname -s)

PROGNAME=nextsim.exec

#NEXTSIMDIR=..

CXX = g++
#CXXFLAGS += -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -g -pthread -fopenmp -v -DHAVE_CONFIG_H
CXXFLAGS += -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread -fopenmp -v -DHAVE_CONFIG_H

ifeq ($(KERNEL),Linux)
	#CXX = g++-4.8
	#CXXFLAGS += -std=c++0x -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -g -lm -pthread -v #-MMD -MP -lm -pthread -v
	CXXFLAGS += -I. -I $(NEXTSIMDIR)/core/include -I /usr/lib/openmpi/include -I $(PETSC_DIR)/include #-I headers
	LDFLAGS += -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi -lmpi_cxx -lmpi -lm -ldl -lhwloc
	LDFLAGS += -Wl,-rpath,$(PETSC_DIR)/lib -Wl,-rpath,$(NEXTSIMDIR)/lib
	LDFLAGS += -L $(PETSC_DIR)/lib -lpetsc -L$(NEXTSIMDIR)/lib -lnextsim
else

	ifeq ($(CXX),clang)
	CXXFLAGS += -stdlib=libc++
	endif

	CXXFLAGS += -I $(NEXTSIMDIR)/core/include -I $(NEXTSIMDIR)/contrib/bamg/include
	CXXFLAGS += -I $(NEXTSIMDIR)/contrib/mapx/include -I $(NEXTSIMDIR)/contrib/interp/include
	CXXFLAGS += -I /usr/local/include -I /opt/local/include/openmpi-mp
	CXXFLAGS +=	-I $(PETSC_DIR)/include -I $(BOOST_DIR)/include -I .
	CXXFLAGS +=	-I $(NETCDF_DIR)/include -I $(GMSH_DIR)/include/gmsh
	CXXFLAGS +=	-I /opt/local/netcdf/include/

	LDFLAGS += -Wl,-rpath,/usr/local/lib -Wl,-rpath,/opt/local/lib/openmpi-mp
	LDFLAGS += -L /usr/local/lib -L /opt/local/lib/openmpi-mp -lmpi_cxx -lmpi -ldl -lstdc++ -lpthread
	LDFLAGS += -Wl,-rpath,$(NETCDF_DIR)/lib -L $(NETCDF_DIR)/lib -lnetcdf_c++4

	ifeq ($(CXX),clang)
	LDFLAGS += -stdlib=libc++
	endif

	LDFLAGS += -Wl,-rpath,$(BOOST_DIR)/lib
	LDFLAGS += -L $(BOOST_DIR)/lib -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi -lboost_date_time

	LDFLAGS += -Wl,-rpath,$(PETSC_DIR)/lib
	LDFLAGS += -L $(PETSC_DIR)/lib -lpetsc

	LDFLAGS += -Wl,-rpath,$(NEXTSIMDIR)/lib
	LDFLAGS += -L $(NEXTSIMDIR)/lib -lbamg
	LDFLAGS += -L $(NEXTSIMDIR)/lib -linterp
	LDFLAGS += -L $(NEXTSIMDIR)/lib -lmapx
	LDFLAGS += -L $(NEXTSIMDIR)/lib -lnextsim

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
