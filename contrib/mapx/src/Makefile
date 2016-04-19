KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	LIBRARYNAME=libmapx.so.1.0
else
	LIBRARYNAME=libmapx.dylib
endif

LIB=$(LIBRARYNAME)
OBJECTDIR=$(NEXTSIMDIR)/objs/
DEPENDDIR=$(NEXTSIMDIR)/.deps/
LIBRARYDIR=$(NEXTSIMDIR)/lib/

CC = gcc
CCFLAGS += -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread -v

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
.PHONY: all clean mrproper

# Default action.
all: $(LIB)

# Delete the object files.
clean:
	$(RM) $(OBJS) $(DEPS) # $(CCOBJS) $(CCDEPS)

mrproper: clean
	$(RM) $(LIBRARYDIR)$(LIBRARYNAME)

# Rule for making the actual target
ifeq ($(KERNEL),Linux)
$(LIBRARYNAME): $(OBJS) #$(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CC) $(CCFLAGS) -o $(LIBRARYDIR)$@ $^ $(LDFLAGS)
	@ln -sf $@ $(LIBRARYDIR)libmapx.so
	@ln -sf $@ $(LIBRARYDIR)libmapx.so.1
else
$(LIBRARYNAME): $(OBJS) #$(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CC) $(CCFLAGS) -o $(LIBRARYDIR)$@ $^ $(LDFLAGS)
	@ln -sf $@ $(LIBRARYDIR)libmapx.dylib.0
	@ln -sf $@ $(LIBRARYDIR)libmapx.dylib.1
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
