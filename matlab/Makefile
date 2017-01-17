KERNEL := $(shell uname -s)
NODENAME := $(shell uname -n)

ifeq ($(KERNEL),Linux)
	ifeq ($(NODENAME),johansen.ad.nersc.no)
		MEX = /opt/MATLAB/R2014b/bin/mex
	else
		MEX = mex
	endif
else
	MEX = /Applications/MATLAB_R2014b.app/bin/mex
endif

MEX += -largeArrayDims

LDFLAGS = -L../lib -I../contrib/mapx/include -lmapx

all: *.c
	$(foreach file, $^, $(MEX) $(file) $(LDFLAGS);)
	#for file in $^; do $(MEX) $$file $(LDFLAGS); done

	 
