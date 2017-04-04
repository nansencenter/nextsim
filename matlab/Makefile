KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	MEX = mex
endif

LDFLAGS = -L../lib -I../contrib/mapx/include -lmapx

all: check-MEX compile

compile: *.c
	$(foreach file, $^, $(MEX) -largeArrayDims $(file) $(LDFLAGS);)

check-MEX:
	@if test "$(MEX)" = "" ; then \
		echo "environment variable MEX is empty"; \
		echo "for OSX, need to export MEX=[path to mex executable]"; \
		exit 1;\
	fi
