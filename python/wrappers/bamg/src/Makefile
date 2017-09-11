# ========================================================
# kernel type
KERNEL := $(shell uname -s)

# main output
LIBRARYNAME=bamg.so
LIBRARYDIR=$(NEXTSIMDIR)/python/lib

# ========================================================


# ========================================================
# Rules to always execute.
.PHONY: default lib install clean vclean fresh

# Default action.
default: lib

# Delete the object and dependancy files.
clean:
	@echo " "
	rm -rf build
	@echo " "

fresh: clean lib

# Delete the object and dependancy files, and the library
vclean: clean
	rm -rf lib
	@echo " "

lib:
	@echo " "
	@echo "Running setup.py... "
	@echo " "
	python setup.py build_ext --inplace
	@echo " "

install: lib
	@mkdir -p $(LIBRARYDIR)
	cp lib/$(LIBRARYNAME) $(LIBRARYDIR)
# ========================================================
