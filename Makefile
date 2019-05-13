# @file   Makefile
# @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
# @date   Tue May 10 10:52:04 2016

.PHONY: all contrib modules model core clean cleanmodel mrproper fresh


all: contrib modules model core

contrib:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; $(MAKE)
	@cd $(NEXTSIMDIR)/contrib/mapx/src; $(MAKE)

modules:
ifdef USE_OASIS
	@cd $(NEXTSIMDIR)/modules/oasis/src; $(MAKE)
endif
ifdef USE_ENSEMBLE
	@cd $(NEXTSIMDIR)/modules/enkf/perturbation/src; $(MAKE)
endif

model: core
	@cd $(NEXTSIMDIR)/model; $(MAKE);

core:
	@cd $(NEXTSIMDIR)/core/src; $(MAKE)

clean: cleanmodel
	@cd $(NEXTSIMDIR)/contrib/bamg/src; $(MAKE) clean
	@cd $(NEXTSIMDIR)/contrib/mapx/src; $(MAKE) clean
	@cd $(NEXTSIMDIR)/modules/oasis/src; $(MAKE) clean
	@cd $(NEXTSIMDIR)/modules/enkf/perturbation/src; $(MAKE) clean
	@cd $(NEXTSIMDIR)/core/src; $(MAKE) clean

cleanmodel:
	@cd $(NEXTSIMDIR)/model; $(MAKE) clean;

mrproper: clean
	@cd $(NEXTSIMDIR)/contrib/bamg/src; $(MAKE) mrproper
	@cd $(NEXTSIMDIR)/contrib/mapx/src; $(MAKE) mrproper
	@cd $(NEXTSIMDIR)/modules/oasis/src; $(MAKE) mrproper
	@cd $(NEXTSIMDIR)/modules/enkf/perturbation/src; $(MAKE) mrproper
	@cd $(NEXTSIMDIR)/core/src; $(MAKE) mrproper
	@cd $(NEXTSIMDIR)/model; $(MAKE) mrproper
	rm -r objs  || true
	rm -r lib   || true
	rm -r .deps || true

fresh: mrproper all
