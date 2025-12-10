KERNEL := $(shell uname -s)

LIBRARYNAME=libmapx.a

OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
LIBRARYDIR=$(NEXTSIMDIR)/lib/


# add gcc option flags
CCFLAGS += -pedantic -Wno-inline
ifdef NEXTSIM_COMPILE_VERBOSE
	CXXFLAGS += -v
endif

ifneq (,$(strip $(filter DEBUG Debug debug,$(NEXTSIM_BUILD_TYPE))))
	CCFLAGS += -g -O0 -DNDEBUG
ifneq ($(KERNEL),Linux)
	CCFLAGS += -Wl,-no_pie
endif
endif


ifeq ($(KERNEL),Linux)
	CCFLAGS += -I$(NEXTSIMDIR)/contrib/mapx/include -I$(MPI_INC_DIR)
else

	ifeq ($(CC),clang)
		CCFLAGS += -stdlib=libc++
	endif
	CCFLAGS += -I$(NEXTSIMDIR)/contrib/mapx/include -I/usr/local/include
endif

# C files
CCSRCDIR=.
CCSRC=$(wildcard $(CCSRCDIR)/*.c)

OBJS=$(CCSRC:%.c=$(OBJECTDIR)%.o)
DEPS=$(CCSRC:%.c=$(DEPENDDIR)%.d)

# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: $(LIBRARYNAME)

# Delete the object files.
clean:
	$(RM) $(OBJS) $(DEPS) # $(CCOBJS) $(CCDEPS)

mrproper: clean
	$(RM) $(LIBRARYDIR)$(LIBRARYNAME)

# Rule for making the actual target
$(LIBRARYNAME): $(OBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "=================================================="
	@echo Creating static library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "=================================================="
	$(AR) rcs $(LIBRARYDIR)$@ $^

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
