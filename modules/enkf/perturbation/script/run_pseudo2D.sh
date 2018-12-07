#!/bin/bash

BOOST_DIR=/opt/local/boost
BOOST_LIBDIR=${BOOST_DIR}/lib

WD=/y1/aydogdu/GIT/nextsim/modules/enkf/perturbation
BIN=./bin/p_pseudo2D 
LIB=./lib/libpseudo2D.dylib
SRC=${WD}/src

export DYLD_LIBRARY_PATH=${BOOST_LIBDIR}:${WD}/lib

echo "number of args: $#"

# Compiles the whole code if there are more than 2 arguments
[[ $# -ge 3 ]] &&\
       	cd ${SRC} && make clean && make all && cd ../

# Compiles the if fortran library is missing
[[ !  -f ${LIB} ]] &&\
       	cd ${SRC} && make clean && make lib && make bin && cd ../

# Compiles the binary if the executable is missing
[[ !  -f ${BIN} ]] &&\
       	cd ${SRC} && make bin && cd ../

[[ $# -eq 0 ]] && \
	echo "Usage: args -> run and/or make and/or plot " &&\
	echo "RUN: ./script/run_pseudo2D.sh run plot make" &&\
	echo "                  OR                       " &&\
	echo "RUN: ./script/run_pseudo2D.sh run make make"

# Runs if there is only one argument
[[ $# -ge 1 ]] &&\
       	cd ${WD} && ${BIN}

# Plots if one of the arguments is 'plot'
[[ $# -ge 2 ]] && ( [[ $1 = 'plot' ]] || [[ $2 = 'plot' ]] ) &&\
       source activate FERRET &&\
       pyferret -png -script POST/frt.pseudo.jnl &&\
       cd data/ASR  && \
       pyferret -png -script frt.asr_forcing.jnl stereo noprint &&\
       source deactivate FERRET
