KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	LIBRARYNAME=libwim.so.1.0
else
	LIBRARYNAME=libwim.dylib
endif

LIB=$(LIBRARYNAME)
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
LIBRARYDIR=$(NEXTSIMDIR)/lib/

CXX = g++
CXXFLAGS += -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread -v

CC = gcc
CFLAGS += -pedantic -Wno-inline -fPIC -O3 #-v # -MMD -MP -lm #-Wall

CXXFLAGS += -I $(OPENMPI_INCLUDE_DIR)/
CXXFLAGS += -I $(BOOST_DIR)/include
CFLAGS += -I $(NEXTSIMDIR)/modules/wim/include

ifeq ($(KERNEL),Linux)
	CXXFLAGS += -I $(NEXTSIMDIR)/modules/wim/include -I /usr/lib/openmpi/include
	LDFLAGS += -shared -Wl,-soname,libwim.so.1
	LDFLAGS += -lm -ldl #-lhwloc
	LDFLAGS += -Wl,-rpath,$(BOOST_DIR)/lib
	LDFLAGS += -L $(BOOST_DIR)/lib -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi -lboost_date_time
else

	ifeq ($(CXX),clang)
	CXXFLAGS += -stdlib=libc++
	endif

	CXXFLAGS += -I $(NEXTSIMDIR)/modules/wim/include -I /usr/local/include
	LDFLAGS += -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-current_version,1.0,-install_name,$(LIBRARYDIR)$(LIBRARYNAME)
	LDFLAGS += -ldl -lstdc++
	LDFLAGS += -Wl,-rpath,$(BOOST_DIR)/lib
	LDFLAGS += -L $(BOOST_DIR)/lib -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi -lboost_date_time

	ifeq ($(CXX),clang)
	LDFLAGS += -stdlib=libc++
	endif
endif

# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)

OBJS=$(CXXSRC:%.cpp=$(OBJECTDIR)%.o)
DEPS=$(CXXSRC:%.cpp=$(DEPENDDIR)%.d)


# C files
CCSRCDIR=.
CCSRC=$(wildcard $(CCSRCDIR)/*.c)

COBJS=$(CCSRC:%.c=$(OBJECTDIR)%.o)
CDEPS=$(CCSRC:%.c=$(DEPENDDIR)%.d)


# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: $(LIB)

# Delete the object files.
clean:
	$(RM) $(OBJS) $(DEPS) $(COBJS) $(CDEPS)

mrproper: clean
	$(RM) $(LIBRARYDIR)$(LIBRARYNAME)

$(COBJS): $(OBJECTDIR)%.o : $(CCSRCDIR)/%.c
	@mkdir -p $(OBJECTDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(CDEPS): $(DEPENDDIR)%.d : $(CCSRCDIR)/%.c
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CC) -M $(CFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'

# Rule for making the actual target
ifeq ($(KERNEL),Linux)
$(LIBRARYNAME): $(OBJS) $(COBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$@ $^ $(LDFLAGS)
	@ln -sf $@ $(LIBRARYDIR)libwim.so
	@ln -sf $@ $(LIBRARYDIR)libwim.so.1
else
$(LIBRARYNAME): $(OBJS) $(COBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$@ $^ $(LDFLAGS)
	@ln -sf $@ $(LIBRARYDIR)libwim.dylib.0
	@ln -sf $@ $(LIBRARYDIR)libwim.dylib.1
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
-include $(CDEPS)
