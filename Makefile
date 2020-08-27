# @file   Makefile
# @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
# @date   Tue May 10 10:52:04 2016

.PHONY: all contrib modules model core clean cleanmodel mrproper fresh


all:
	$(MAKE) contrib
	$(MAKE) modules
	$(MAKE) model

contrib:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; $(MAKE)
	@cd $(NEXTSIMDIR)/contrib/mapx/src; $(MAKE)

modules:
ifdef USE_OASIS
	@cd $(NEXTSIMDIR)/modules/oasis/src; $(MAKE)
endif
ifdef USE_ENSEMBLE
	@cd $(NEXTSIMDIR)/modules/enkf/perturbation/src; make
	@cd $(NEXTSIMDIR)/modules/enkf/gridutils-c; $(MAKE)
	@cd $(NEXTSIMDIR)/modules/enkf/enkf-c; $(MAKE)
endif

docker: core modules
	@cd $(NEXTSIMDIR)/model; $(MAKE);

model: core modules contrib
	@cd $(NEXTSIMDIR)/model; $(MAKE);

core:
	@cd $(NEXTSIMDIR)/core/src; $(MAKE)

clean: cleanmodel
	@cd $(NEXTSIMDIR)/contrib/bamg/src; $(MAKE) clean
	@cd $(NEXTSIMDIR)/contrib/mapx/src; $(MAKE) clean
ifdef USE_OASIS
	@cd $(NEXTSIMDIR)/modules/oasis/src; $(MAKE) clean
endif
ifdef USE_ENSEMBLE
	@cd $(NEXTSIMDIR)/modules/enkf/perturbation/src; $(MAKE) clean
	@cd $(NEXTSIMDIR)/modules/enkf/gridutils-c; $(MAKE) clean
	@cd $(NEXTSIMDIR)/modules/enkf/enkf-c; $(MAKE) clean
endif
	@cd $(NEXTSIMDIR)/core/src; $(MAKE) clean

cleanmodel:
	@cd $(NEXTSIMDIR)/model; $(MAKE) clean;

mrproper: clean
	@cd $(NEXTSIMDIR)/contrib/bamg/src; $(MAKE) mrproper
	@cd $(NEXTSIMDIR)/contrib/mapx/src; $(MAKE) mrproper
ifdef USE_OASIS
	@cd $(NEXTSIMDIR)/modules/oasis/src; $(MAKE) mrproper
endif
ifdef USE_ENSEMBLE
	@cd $(NEXTSIMDIR)/modules/enkf/perturbation/src; $(MAKE) mrproper
endif
	@cd $(NEXTSIMDIR)/core/src; $(MAKE) mrproper
	@cd $(NEXTSIMDIR)/model; $(MAKE) mrproper
	rm -rf objs
	rm -rf lib
	rm -rf .deps

fresh:
	$(MAKE) mrproper
	$(MAKE) all
