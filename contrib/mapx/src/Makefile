KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	LIBRARYNAME=libmapx.so.1.0
else
	LIBRARYNAME=libmapx.dylib
endif

OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
LIBRARYDIR=$(NEXTSIMDIR)/lib/
tools: LIBRARYDIR=$(NEXTSIMTOOLS_ROOT_DIR)/lib/nextsim/

CC = gcc

# add gcc option flags
CCFLAGS += -pedantic -Wno-inline -fPIC
ifneq ($(XSIMCOMP_VERBOSE),no)
   	CXXFLAGS += -v
endif

ifneq (,$(strip $(filter DEBUG Debug debug,$(NEXTSIM_BUILD_TYPE))))
	#ifeq ($(NEXTSIM_BUILD_TYPE),$(filter $(NEXTSIM_BUILD_TYPE),Debug debug))
	#CCFLAGS := $(filter-out -O3 -pthread,$(CXXFLAGS))
	CCFLAGS += -g -O0 -DDEBUGGING
ifneq ($(KERNEL),Linux)
	CCFLAGS += -Wl,-no_pie
endif
else
	CCFLAGS += -O3 -pthread  
endif


ifeq ($(KERNEL),Linux)
	CCFLAGS += -I $(NEXTSIMDIR)/contrib/mapx/include -I /usr/lib/openmpi/include
	LDFLAGS += -shared -Wl,-soname,libmapx.so.1
	LDFLAGS += -lm -ldl #-lhwloc
else

	ifeq ($(CC),clang)
	CCFLAGS += -stdlib=libc++
	endif

	CCFLAGS += -I $(NEXTSIMDIR)/contrib/mapx/include -I /usr/local/include
	LDFLAGS += -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-current_version,1.0,-install_name,$(LIBRARYDIR)$(LIBRARYNAME)
	LDFLAGS += -ldl -lstdc++

	ifeq ($(CC),clang)
	LDFLAGS += -stdlib=libc++
	endif

endif

# C files
CCSRCDIR=.
CCSRC=$(wildcard $(CCSRCDIR)/*.c)

OBJS=$(CCSRC:%.c=$(OBJECTDIR)%.o)
DEPS=$(CCSRC:%.c=$(DEPENDDIR)%.d)

# Rules to always execute.
.PHONY: default lib tools clean mrproper

# Default action.
default: lib

tools: clean lib

# Delete the object files.
clean:
	$(RM) $(OBJS) $(DEPS) # $(CCOBJS) $(CCDEPS)

mrproper: clean
	$(RM) $(LIBRARYDIR)$(LIBRARYNAME)

# Rule for making the actual target
ifeq ($(KERNEL),Linux)
lib: $(OBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CC) $(CCFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGS)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libmapx.so
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libmapx.so.1
else
lib: $(OBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CC) $(CCFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGS)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libmapx.dylib.0
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libmapx.dylib.1
endif

# Rules for object files from cpp files
$(OBJECTDIR)%.o: %.c
	@mkdir -p $(OBJECTDIR) # $(dir $@)
	$(CC) -o $@ -c $< $(CCFLAGS)

# Make dependancy rules
$(DEPENDDIR)%.d: %.c
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CC) -M $(CCFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'


# The compilation depends on this Makefile.
$(OBJS): Makefile

# Include the dependency files
-include $(DEPS)
