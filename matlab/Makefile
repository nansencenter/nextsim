KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	MEX = mex
else
	MEX = /Applications/MATLAB_R2016b.app/bin/mex
endif

MEX += -largeArrayDims

LDFLAGS = -L../lib -I../contrib/mapx/include -lmapx

all: *.c
	$(foreach file, $^, $(MEX) $(file) $(LDFLAGS);)
	#for file in $^; do $(MEX) $$file $(LDFLAGS); done

	 
