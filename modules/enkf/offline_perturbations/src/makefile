INCPATH = .:../include

OBJS=../objs


.SUFFIXES:
.SUFFIXES: .o .F90 .


CF90 = mpiifort #compiler -- F90
CF77 = $(CF90)  #Compiler -- F77
LD =  $(CF90)   #Linker 
#
#
# Paralellization opts
# PARO = -lmpi -parallel
#
# Size defaults - Change to real_size 64 for default double...
# SIZEO = -real_size 64 -double_size 64
# SIZEO = -r8 
#
# #rch opts
# ARCHO= -xCORE-AVX2 -O2 -no-prec-div -qopt-prefetch=2 -auto-p32 -no-ansi-alias -qopt-mem-layout-trans=2 
#
# #Optimalization opts
OPTO= -g -fpic -nofor-main -convert big_endian -assume byterecl -cm -vec_reprot0 -w -O2 -r8 -g -traceback
#
#
F77FLG = -nofree
F90FLG = -free
CPPFLAGS += $(CPPARCH) -DFFTW -DIARGC -DLAPACK # -DQMPI
#
#
# This uses the OpenMPI implementation of mpi-2
FLIBS = -mkl=parallel -lfftw3 -lgsl -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -lsz -lz -ldl -lm -liomp5 -lpthread -lcurl
#
# Include dir for header and module files
INCLUDEDIR= -I/local/netcdf/include -I/local/fftw/include
#
#
# # Put together flags
FFLAGS    = $(SIZEO) $(OPTO) $(ARCHO) $(PARO)  $(INLO) $(DIVO) $(DEBUG_FLAGS) $(INCLUDEDIR)
#FLINKFLAGS = -shared#$(SIZEO) $(OPTO) $(PARO) $(INLO) $(DIVO) $(DEBUG_FLAGS) -nofor-main

ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif


# Rules for running cpp and updating files in include directory
.F90.o:
	cd $(OBJS) ; $(CF90) -c $(CPPFLAGS) $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o ../src/$<

.F.o:
	cd $(OBJS) ; $(CF77) -c $(CPPFLAGS) $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o ../src/$<


P_E = ../bin/p_pseudo2D


all  : $(P_E)

############################################################################# 
#	memory_stack.o 
FOBJECTS= \
	m_set_random_seed.o \
	mod_pseudo.o \
	mod_random_forcing.o \
	main_pseudo2D.o

$(P_E): $(FOBJECTS)
	cd $(OBJS) ; $(CF90) $(FLINKFLAGS) $(FOBJECTS) $(FLIBS) -o $(P_E)
#############################################################################

bin  : $(P_E)

clean:
	cd $(OBJS); rm -f *.o *.mod
	rm -f $(P_E)
