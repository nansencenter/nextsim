KERNEL := $(shell uname -s)

PROGNAME=nextsim.exec
WIM_PROG=nextwim.exec
FPROGNAME=FWrap.x

CXX = g++
CXXFLAGS += -std=c++11 -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread -fopenmp -v -DHAVE_CONFIG_H -D_MULTITHREADING_

CXXFLAGS += -pg -g -DNDEBUG

CXXFLAGS += -I $(NEXTSIMDIR)/core/include -I $(NEXTSIMDIR)/contrib/bamg/include
CXXFLAGS += -I $(NEXTSIMDIR)/contrib/mapx/include -I $(NEXTSIMDIR)/contrib/interp/include

# ifdef WIM2D_PATH
# 	CXXFLAGS += -I $(WIM2D_PATH)/CXX/include
# else
# 	CXXFLAGS += -I $(NEXTSIMDIR)/modules/wim/include
# endif

CXXFLAGS += -I $(OPENMPI_INCLUDE_DIR)/
CXXFLAGS +=	-I $(PETSC_DIR)/include -I $(BOOST_DIR)/include -I .
CXXFLAGS +=	-I $(NETCDF_DIR)/include -I $(GMSH_DIR)/include/gmsh

FC = gfortran
FFLAGS = -std=f2008 -fopenmp -lstdc++

ifeq ($(KERNEL),Linux)
	#CXXFLAGS += -std=c++0x -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -g -lm -pthread -v #-MMD -MP -lm -pthread -v

else

	ifeq ($(CXX),clang)
	CXXFLAGS += -stdlib=libc++
	endif

	CXXFLAGS += -I /usr/local/include #-I /opt/local/include/openmpi-mp
	#CXXFLAGS +=	-I /opt/local/netcdf/include/

	#LDFLAGS += -Wl,-rpath,/usr/local/lib #-Wl,-rpath,/opt/local/lib/openmpi-mp
	#LDFLAGS += -L /usr/local/lib #-L /opt/local/lib/openmpi-mp -lmpi_cxx -lmpi -ldl -lstdc++ -lpthread

	ifeq ($(CXX),clang)
	LDFLAGS += -stdlib=libc++
	endif

endif

# ifdef WIM2D_PATH
# 	LDFLAGS += -Wl,-rpath,$(WIM2D_PATH)/CXX/lib
# 	LDFLAGS += -L $(WIM2D_PATH)/CXX/lib -lwim
# endif

LDFLAGS += -L /usr/local/lib

LDFLAGS += -lprofiler

LDFLAGS += -Wl,-rpath,$(OPENMPI_LIB_DIR)/
ifndef MACHINE_HEXAGON
	LDFLAGS += -L $(OPENMPI_LIB_DIR)/ -lmpi_cxx -lmpi -ldl -lstdc++ -lpthread
else
	LDFLAGS += -L $(OPENMPI_LIB_DIR)/ -lmpichcxx -lmpich -ldl -lstdc++ -lpthread #-lssl -luuid -lpthread -lrt
endif

LDFLAGS += -Wl,-rpath,$(NETCDF_DIR)/lib -L $(NETCDF_DIR)/lib -lnetcdf_c++4

LDFLAGS += -Wl,-rpath,$(BOOST_DIR)/lib
LDFLAGS += -L $(BOOST_DIR)/lib -lboost_program_options -lboost_filesystem -lboost_system -lboost_serialization -lboost_mpi -lboost_date_time

# LDFLAGS += -Wl,-rpath,$(PETSC_DIR)/lib
# LDFLAGS += -L $(PETSC_DIR)/lib -lpetsc

LDFLAGS += -Wl,-rpath,$(NEXTSIMDIR)/lib
LDFLAGS += -L $(NEXTSIMDIR)/lib -lbamg
LDFLAGS += -L $(NEXTSIMDIR)/lib -linterp
LDFLAGS += -L $(NEXTSIMDIR)/lib -lmapx

# ifndef WIM2D_PATH
# 	LDFLAGS += -L $(NEXTSIMDIR)/lib -lwim
# endif

LDFLAGS += -L $(NEXTSIMDIR)/lib -lnextsim


EXEC=$(PROGNAME)
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
BINARYDIR=bin/

# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)

OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)

# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: $(EXEC)

# Delete the object files.
clean:
	$(RM) $(OBJS) $(DEPS)

mrproper: clean
	$(RM) $(BINARYDIR)$(PROGNAME)

# Rule for making the actual target
wim:	CXXFLAGS += -I $(WIM2D_PATH)/CXX/include -DWAVES
wim:    LDFLAGS += -Wl,-rpath,$(WIM2D_PATH)/CXX/lib -L $(WIM2D_PATH)/CXX/lib -lwim

# ifeq ($(USE_WIM),0)
# CXXFLAGS += -I $(NEXTSIMDIR)/modules/wim/include
# LDFLAGS += -L $(NEXTSIMDIR)/lib -lwim
# USE_WIM2 := 0
# else
# CXXFLAGS += -I $(WIM2D_PATH)/CXX/include -DWAVES
# LDFLAGS += -Wl,-rpath,$(WIM2D_PATH)/CXX/lib
# LDFLAGS += -L $(WIM2D_PATH)/CXX/lib -lwim
# USE_WIM2 := 1
# endif


$(PROGNAME): $(OBJS) #$(CCOBJS)
	@mkdir -p $(BINARYDIR)
	$(CXX) $(CXXFLAGS) -o $(BINARYDIR)$@ $^ $(LDFLAGS)


wim: $(OBJS) #$(CCOBJS)
	@mkdir -p $(BINARYDIR)
	$(CXX) $(CXXFLAGS) -o $(BINARYDIR)$(WIM_PROG) $^ $(LDFLAGS)

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

# The fortran wrapper
fortran: $(FPROGNAME)

$(FPROGNAME): $(filter-out $(OBJECTDIR)$(CXXSRCDIR)/main.o, $(OBJS)) finiteelement_mod.f90
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(BINARYDIR)$@ $^ FE.f90

-include $(DEPS)
