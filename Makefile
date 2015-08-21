# VPATH=.

NEXTSIMDIR=.

all:
	@cd $(NEXTSIMDIR)/src; make;
	@cd $(NEXTSIMDIR)/bamg/src; make
	@cd $(NEXTSIMDIR)/mapx/src; make

clean:
	@cd $(NEXTSIMDIR)/src; make clean;
	@cd $(NEXTSIMDIR)/bamg/src; make clean
	@cd $(NEXTSIMDIR)/mapx/src; make clean

mrproper: clean
	@cd $(NEXTSIMDIR)/src; make clean mrproper
	@cd $(NEXTSIMDIR)/bamg/src; make clean mrproper
	@cd $(NEXTSIMDIR)/mapx/src; make clean mrproper
