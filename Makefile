# @file   Makefile
# @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
# @date   Tue May 10 10:52:04 2016

.PHONY: all clean mrproper All Clean Mrproper fresh

all:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make
	#@cd $(NEXTSIMDIR)/contrib/interp/src; make
	@cd $(NEXTSIMDIR)/core/src; make;
ifeq ($(USE_NEXTWIM),1)
	@cd $(NEXTSIMDIR)/modules/wim/src; make;
endif

clean:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make clean
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make clean
	#@cd $(NEXTSIMDIR)/contrib/interp/src; make clean
	@cd $(NEXTSIMDIR)/core/src; make clean

mrproper: clean
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make clean mrproper
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make clean mrproper
	#@cd $(NEXTSIMDIR)/contrib/interp/src; make clean mrproper
	@cd $(NEXTSIMDIR)/core/src; make clean mrproper

wim:
	@cd $(NEXTSIMDIR)/modules/wim/src; make;

# rules to compile model code as well as lib's
All: all
	@cd $(NEXTSIMDIR)/model; make;

Clean: clean
	@cd $(NEXTSIMDIR)/model; make clean;

Mrproper: mrproper
	@cd $(NEXTSIMDIR)/model; make mrproper;

# fresh compile (clean first)
fresh: Clean All
