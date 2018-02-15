# @file   Makefile
# @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
# @date   Tue May 10 10:52:04 2016

.PHONY: all clean mrproper All Clean fresh core bamg mapx model wim oasis tools clean-tools

# compile all the libraries; these are saved to $(NEXTSIMDIR)/lib as usual
all: bamg mapx core
	#@cd $(NEXTSIMDIR)/contrib/interp/src; make
ifdef USE_NEXTWIM
	@cd $(NEXTSIMDIR)/modules/wim/src; make
endif
ifdef USE_OASIS
	@cd $(NEXTSIMDIR)/modules/oasis/src; make
endif

clean:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make clean
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make clean
	#@cd $(NEXTSIMDIR)/contrib/interp/src; make clean
	@cd $(NEXTSIMDIR)/modules/wim/src; make clean
	@cd $(NEXTSIMDIR)/core/src; make clean
ifdef USE_OASIS
	@cd $(NEXTSIMDIR)/modules/oasis/src; make clean
endif

mrproper:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make mrproper
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make mrproper
	#@cd $(NEXTSIMDIR)/contrib/interp/src; make mrproper
	@cd $(NEXTSIMDIR)/modules/wim/src; make mrproper
	@cd $(NEXTSIMDIR)/core/src; make mrproper
ifdef USE_OASIS
	@cd $(NEXTSIMDIR)/modules/oasis/src; make mrproper
endif
	@cd $(NEXTSIMDIR)/model; make mrproper;


# rules to compile model code as well as lib's
# - NB doesn't work on osx
# - still need to do "cd model;make clean;make" afterwards
All: all model

Clean: clean
	@cd $(NEXTSIMDIR)/model; make clean;

# fresh compile (clean first)
# - NB doesn't work on osx
# - still need to do "cd model;make clean;make" afterwards
fresh: mrproper All

bamg:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make
mapx:
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make
core:
	@cd $(NEXTSIMDIR)/core/src; make;
model:
	@cd $(NEXTSIMDIR)/model; make;
wim:
	@cd $(NEXTSIMDIR)/modules/wim/src; make
oasis:
	@cd $(NEXTSIMDIR)/modules/oasis/src; make

tools:
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make tools
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make tools; make tools-no-omp

clean-tools:
	rm -f $(NEXTSIMTOOLS_ROOT_DIR)/lib/nextsim/*
