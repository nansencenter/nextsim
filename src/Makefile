LIBRARYNAME=libnextsim.so.1.0

NEXTSIMDIR=..

CXX = g++-4.8
CXXFLAGS += -std=c++0x -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -g -lm -pthread -v #-MMD -MP -lm -pthread -v
CXXFLAGS += -I $(NEXTSIMDIR)/include -I /usr/lib/openmpi/include -I $(PETSC_DIR)/include #-I headers
LDFLAGS += -shared -Wl,-soname,libnextsim.so.1 -lboost_program_options -lboost_filesystem -lboost_system -lmpi_cxx -lmpi -ldl -lhwloc -lboost_mpi
LDFLAGS += -Wl,-rpath,$(PETSC_DIR)/lib
LDFLAGS += -L $(PETSC_DIR)/lib -lpetsc

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
$(LIBRARYNAME): $(OBJS) #$(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	$(CXX) $(CXXLAGS) -o $(LIBRARYDIR)$@ $^ $(LDFLAGS)
	@ln -sf $@ $(LIBRARYDIR)libnextsim.so
	@ln -sf $@ $(LIBRARYDIR)libnextsim.so.1


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
