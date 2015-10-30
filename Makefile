# VPATH=.

#NEXTSIMDIR=.

all:
	@cd $(NEXTSIMDIR)/core/src; make;
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make
	@cd $(NEXTSIMDIR)/contrib/interp/src; make

clean:
	@cd $(NEXTSIMDIR)/core/src; make clean;
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make clean
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make clean
	@cd $(NEXTSIMDIR)/contrib/interp/src; make clean

mrproper: clean
	@cd $(NEXTSIMDIR)/core/src; make clean mrproper
	@cd $(NEXTSIMDIR)/contrib/bamg/src; make clean mrproper
	@cd $(NEXTSIMDIR)/contrib/mapx/src; make clean mrproper
	@cd $(NEXTSIMDIR)/contrib/interp/src; make clean mrproper
