FROM nansencenter/boost_petsc_gmsh:0.0.6

RUN apt-get clean \
&&  apt-get update && apt-get install -yq --no-install-recommends \
    libnetcdf-dev \
    libnetcdf-c++4-dev \
    ssh \
    valgrind \
    libhdf5-dev\
    libfftw3-dev \
    bc \
&&  apt-get clean \
&&  rm -rf /var/lib/apt/lists/*

ENV OPENMPI_INCLUDE_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include \
    OPENMPI_LIB_DIR=/usr/lib/x86_64-linux-gnu/openmpi/lib \
    BOOST_INCDIR=/opt/local/boost/include \
    BOOST_LIBDIR=/opt/local/boost/lib \
    PETSC_DIR=/opt/local/petsc \
    GMSH_DIR=/opt/local/gmsh \
    PATH=$PATH:/nextsim/model/bin \
    NEXTSIMDIR=/nextsim \
    NEXTSIM_MESH_DIR=/mesh \
    NEXTSIM_DATA_DIR=/data \
    LIBRARY_PATH=/opt/local/mapx/lib:/opt/local/bamg/lib 

# copy mapx source, compile and copy libs and include into
# /opt/local/mapx/lib and /opt/local/mapx/include
COPY contrib/mapx $NEXTSIMDIR/contrib/mapx
WORKDIR $NEXTSIMDIR/contrib/mapx/src

RUN make -j8 \
&&  mkdir -p /opt/local/mapx/lib  \
&&  cp $NEXTSIMDIR/lib/libmapx* /opt/local/mapx/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/mapx/include /opt/local/mapx

# copy bamg source, compile and copy libs and include into
# /opt/local/bamg/lib and /opt/local/bamg/include
COPY contrib/bamg $NEXTSIMDIR/contrib/bamg
WORKDIR $NEXTSIMDIR/contrib/bamg/src
RUN make -j8 \
&&  mkdir -p /opt/local/bamg/lib  \
&&  cp $NEXTSIMDIR/lib/libbamg* /opt/local/bamg/lib/ \
&&  cp -r $NEXTSIMDIR/contrib/bamg/include /opt/local/bamg

# --- scripts below are added for enkf ---------
# copy core source and compile
COPY core $NEXTSIMDIR/core
WORKDIR $NEXTSIMDIR/core/src
RUN make -j8

ENV USE_ENSEMBLE=1 \
FFTW_DIR=/usr
# copy modules and compile
COPY modules $NEXTSIMDIR/modules
# COPY scripts/ensemble/pseudo2D.nml ${NEXTSIMDIR}/.
WORKDIR $NEXTSIMDIR/modules/enkf/perturbation/src
RUN make all 
WORKDIR $NEXTSIMDIR/modules/enkf/gridutils-c
RUN make all -j8 \
&& cp $NEXTSIMDIR/modules/enkf/gridutils-c/libgu.so $NEXTSIMDIR/lib/libgu.so
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NEXTSIMDIR/lib

WORKDIR $NEXTSIMDIR/modules/enkf/enkf-c
RUN make all -j8

# copy model source and compile
COPY model $NEXTSIMDIR/model
WORKDIR $NEXTSIMDIR/model
RUN make -j8

WORKDIR $NEXTSIMDIR
