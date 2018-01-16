# VPATH=.

#NEXTSIMDIR=.

all:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make
	@cd $(NEXTSIMDIR)/contrib/interp/src; make
	@cd $(NEXTSIMDIR)/core/src; make;

clean:
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make clean
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make clean
	@cd $(NEXTSIMDIR)/contrib/interp/src; make clean
	@cd $(NEXTSIMDIR)/core/src; make clean;

mrproper: clean
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make clean mrproper
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make clean mrproper
	@cd $(NEXTSIMDIR)/contrib/interp/src; make clean mrproper
	@cd $(NEXTSIMDIR)/core/src; make clean mrproper
# rules to compile model code as well as lib's
All: all
	@cd $(NEXTSIMDIR)/model; make;

Clean: clean
	@cd $(NEXTSIMDIR)/model; make clean;

Mrproper: mrproper
	@cd $(NEXTSIMDIR)/model; make mrproper;

# fresh compile (clean first)
fresh: Clean All
