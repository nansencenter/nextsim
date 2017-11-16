KERNEL := $(shell uname -s)
NODENAME := $(shell uname -n)

ifeq ($(KERNEL),Linux)
	ifeq ($(NODENAME),johansen.ad.nersc.no)
		MEX = /opt/MATLAB/R2014b/bin/mex
	else
		MEX = mex
	endif
else
	MEX = /Applications/MATLAB_R2016b.app/bin/mex
endif

MEX += -largeArrayDims

# may need to specify g++ version for the cpp files
# - to do this set $MEXCPP to path to correct g++
MEX2 = $(MEX)
ifdef MEXCPP
	MEX2 += GCC=$(MEXCPP)
endif

LDFLAGS = -L../lib -I../contrib/mapx/include -lmapx -I../contrib/bamg/include -lbamg

all: c cpp

c: *.c
	$(foreach file, $^, $(MEX) $(file) $(LDFLAGS);)

cpp: *.cpp
	$(foreach file, $^, $(MEX2) $(file) $(LDFLAGS);)

clean:
	rm -f *.mexa64 *.meximac64
