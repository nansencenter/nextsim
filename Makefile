# @file   Makefile
# @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
# @date   Tue May 10 10:52:04 2016

.PHONY: all contrib modules model core clean cleanmodel mrproper fresh


all: contrib modules model core

contrib:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make

modules:
ifdef USE_OASIS
	@cd $(NEXTSIMDIR)/modules/oasis/src; make
endif
ifdef USE_ENSEMBLE
	@cd $(NEXTSIMDIR)/modules/enkf/perturbation/src; make
	@cd $(NEXTSIMDIR)/modules/enkf/gridutils-c; make
	@cd $(NEXTSIMDIR)/modules/enkf/enkf-c; make
endif

model: core
	@cd $(NEXTSIMDIR)/model; make;

core:
	@cd $(NEXTSIMDIR)/core/src; make

clean: cleanmodel
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make clean
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make clean
	@cd $(NEXTSIMDIR)/modules/oasis/src; make clean
	@cd $(NEXTSIMDIR)/core/src; make clean
	@cd $(NEXTSIMDIR)/modules/enkf/perturbation/src; make clean
	@cd $(NEXTSIMDIR)/modules/enkf/gridutils-c; make clean
	@cd $(NEXTSIMDIR)/modules/enkf/enkf-c; make clean

cleanmodel:
	@cd $(NEXTSIMDIR)/model; make clean;

mrproper: clean
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make mrproper
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make mrproper
	@cd $(NEXTSIMDIR)/modules/oasis/src; make mrproper
	@cd $(NEXTSIMDIR)/modules/enkf/perturbation/src; make mrproper
	@cd $(NEXTSIMDIR)/core/src; make mrproper
	@cd $(NEXTSIMDIR)/model; make mrproper
	rm -r objs  || true
	rm -r lib   || true
	rm -r .deps || true

fresh: mrproper all
