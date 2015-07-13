PROGNAME=nextsim.exec

CXX = g++-4.8
CXXFLAGS += -std=c++0x -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -MMD -MP -lm -pthread -v #-Wall
CXXFLAGS += -I. -I include -I /usr/lib/openmpi/include -I /opt/local/petsc/include #-I headers
LDFLAGS += -lboost_program_options -lboost_filesystem -lboost_system -lmpi_cxx -lmpi -ldl -lhwloc -lboost_mpi -L/opt/local/petsc/lib -lpetsc

EXEC=$(PROGNAME)
OBJECTDIR=objs/
BINARYDIR=bin/

# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
OBJS=$(CXXSRC:.cpp=.o)
DEPS=$(CXXSRC:.cpp=.d)

OBJS:=$(addprefix $(OBJECTDIR), $(OBJS))
DEPS:=$(addprefix $(OBJECTDIR), $(DEPS))

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

$(PROGNAME): $(OBJS) #$(CCOBJS)
	@mkdir -p $(BINARYDIR)
	$(CXX) $(CXXLAGS) -o $(BINARYDIR)$@ $^ $(LDFLAGS)


$(OBJECTDIR)%.o: %.cpp
	@mkdir -p $(OBJECTDIR)
	$(CXX) -o $@ -c $< $(CXXFLAGS)

# The compilation depends on this Makefile.
$(OBJS): Makefile

-include $(DEPS)
