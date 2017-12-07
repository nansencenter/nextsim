# @file   Makefile
# @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
# @date   Tue May 10 10:52:04 2016

all:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make
	#@cd $(NEXTSIMDIR)/contrib/interp/src; make
	@cd $(NEXTSIMDIR)/core/src; make;
ifdef USE_NEXTWIM
	@cd $(NEXTSIMDIR)/modules/wim/src; make
endif
ifdef USE_OASIS
	@cd $(NEXTSIMDIR)/modules/oasis/src; make
endif
	@cd $(NEXTSIMDIR)/model; make;

clean:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make clean
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make clean
	#@cd $(NEXTSIMDIR)/contrib/interp/src; make clean
	@cd $(NEXTSIMDIR)/modules/wim/src; make clean
	@cd $(NEXTSIMDIR)/core/src; make clean
ifdef USE_OASIS
	@cd $(NEXTSIMDIR)/modules/oasis/src; make clean
endif
	@cd $(NEXTSIMDIR)/model; make clean;

mrproper: clean
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make clean mrproper
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make clean mrproper
	#@cd $(NEXTSIMDIR)/contrib/interp/src; make clean mrproper
	@cd $(NEXTSIMDIR)/modules/wim/src; make clean mrproper
	@cd $(NEXTSIMDIR)/core/src; make clean mrproper
ifdef USE_OASIS
	@cd $(NEXTSIMDIR)/modules/oasis/src; make clean mrproper
endif
	@cd $(NEXTSIMDIR)/model; make clean mrproper;
